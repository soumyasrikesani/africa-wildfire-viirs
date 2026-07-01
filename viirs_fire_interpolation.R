library(sf)
library(tidyverse)
library(lubridate)
library(spdep)
library(igraph)
library(concaveman)
library(terra)
library(gstat)
library(ggnetwork)

# Read in data (only geo_area is currently used)
africa <- st_read("ESRI_Africa_countries.shp")
africa <- st_transform(africa, 4326)
ghana <- filter(africa, COUNTRY == "Ghana")
geo_area <- st_read("EW_Guinean_Lowland_Forests.shp")
# Reproject to WGS84
geo_area <- st_transform(geo_area, 4326)
# Reproject to Albers Equal Area Africa projection
geo_area_aea <- st_transform(geo_area, st_crs("ESRI:102022"))

# Geographic area name to include in output files
geo_name <- "ugf"
# Maximum distance to consider two points part of the same fire "event"
maxdist <- 1000
# Maximum time difference to consider two points part of the same fire "event"
maxtime <- 2
# Fire confidence levels to include
# NOTE: try re-running with only n and h
fireconf <- c("n", "h", "l")

#####################################################################
# Loop #1 Process VIIRS active fire data for each year and identify
# clusters of active fires occuring at similar locations and time.
#####################################################################

for(curyear in 2016:2016) {

  ###################################################
  # Step 1.1 Read in and process the active fire data 
  ###################################################
  print(paste0("year=", curyear))
  infile1 <- paste0("africa_fires_", curyear, ".csv")
  firepts1 <- read_csv(infile1) %>%
    # For duplicate geometries, keep only one of each
    distinct(longitude, latitude, .keep_all = TRUE) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  # Clip to area of interest and add day of year and year fields
  firepts1 <- firepts1 %>%
    filter(lengths(st_covered_by(., geo_area)) > 0) %>%
    mutate(fdoy = yday(acq_date),
           fyear = year(acq_date)) %>%
    filter(confidence %in% fireconf)
  # Convert to sf and reproject to AEA
  fire_clip <- st_as_sf(firepts1, coords = c("longitude", "latitude"))
  fire_clip <- st_transform(fire_clip, st_crs("ESRI:102022"))

  #############################################################################
  # Step 1.2 Create a graph object based on geographic distances and time 
  # differences
  #############################################################################
  
  # matrix of active fire coordinates
  fcoords <- st_coordinates(fire_clip)
  # ID numbers for active fire points
  IDs <- row.names(as.data.frame(fcoords))
  # Use spdep to create a neighbor object based on the Delauny triangulation
  
  duplicate_rows <- as.data.frame(fcoords) %>%
    group_by(across(everything())) %>% # Groups by all columns
    filter(n() > 1) %>%                # Keeps only groups with more than 1 row
    ungroup()     
  
  fires_nb <- tri2nb(fcoords, row.names=IDs)
  # Extract day of year
  fdate <- fire_clip$fdoy
  # Set up a template for the output data frame
  outdf <- data.frame("from" = integer(),
                      "to" = integer(),
                      "sdist" = double(),
                      "tdist" = double())
  
  # Loop through all active fires
  for(x in 1:length(IDs)) {
    #print(x)
    # Extract the neighbors of the current active fire point
    nlist <- fires_nb[x]
    nneigh <- 0
    # Loop through the neighbors
    for(y in 1:length(nlist[[1]])) {
      # Compute geographic distance and time difference
      sdist <- sqrt((fcoords[x, 1] - fcoords[nlist[[1]][y], 1])^2 + 
                      (fcoords[x, 2] - fcoords[nlist[[1]][y], 2])^2)
      tdist <- abs(fdate[x] - fdate[nlist[[1]][y]])
      # Only include points less than the threshold distance and time
      if((sdist < maxdist && tdist <= maxtime)) {
        nneigh <- nneigh + 1
        outdf <- outdf %>%
          add_row("from" = x,
                  "to" = nlist[[1]][y],
                  "sdist" = sdist,
                  "tdist" = tdist)
      }
    }
    # Do this if no points meet the distance and time criteria
    if(nneigh == 0) {
      outdf <- outdf %>%
        add_row("from" = x,
                "to" = x,
                "sdist" = 0,
                "tdist" = 0)
      
    }
  }

  # Create an igraph object
  firegraph <- graph_from_data_frame(outdf, directed=FALSE, vertices=NULL)
  # Take a look at the igraph object
  #min(components(firegraph)$csize)
  #length(components(firegraph)$membership)
  #components(firegraph)$membership[1:20]
  
  # Number of the cluster to which each active fire belongs
  clustnum <- components(firegraph)$membership
  # Size of the cluster to which each active fire belongs
  clustsize <- components(firegraph)$csize[clustnum]
  # Add these variables to the main active fire data frame
  fire_clip$clustnum <- clustnum
  fire_clip$clustsize <- clustsize
  
  #########################################
  # Step 1.3 Save files with processed data
  #########################################
  
  # Output the main active fire data as a shapefile
  outname <- paste0(geo_name, "_fires_", curyear, ".shp")
  st_write(fire_clip, outname, append = TRUE)
  # Output the graph information as an edgelist
  outname2 <- paste0(geo_name, "_edges_", curyear)
  write_graph(firegraph, outname2, format = "edgelist")

}

#####################################################################
# Loop #2 Convert the active fire points into gridded burned area
# estimates
#####################################################################

for(curyear in 2013:2024) {

  ###################################
  # Step 2.1 Read in and process data
  ###################################
  
  # Read in the shapefiles created in the previous loop
  inname <- paste0(geo_name, "_fires_", curyear, ".shp")
  fire_clip <- st_read(inname)
  # Separate into large (3 or more points) and small clusters
  fires_lg <- filter(fire_clip, clustsize >= 3)
  fires_sm <- filter(fire_clip, clustsize < 3)
  
  #########################################################################
  # Step 2.2 Convert large point cluster into polygons using the concaveman
  # concave hull algorithm
  #########################################################################
  
  # Create an empty list
  polylist <- list()
  # Loop through all the large fires clusters
  for(clust in unique(fires_lg$clustnum)) {
    print(clust)
    curpts <- filter(fires_lg, clustnum == clust)
    # Extract concave hull and add it as a list element
    curpoly <- concaveman(curpts, concavity = 2)
    polylist[[length(polylist) + 1]] <- curpoly
  }
  
  # Convert the list of polygons to a multipolygon sf object
  fire_poly <- bind_rows(polylist)
  # Union the polygons into a single feature
  fire_poly <- st_union(fire_poly)
  # Convert back to sf object
  fire_poly <- st_as_sf(fire_poly)
  
  #########################################################################
  # Step 2.3 Interpolate missing data (date or burning and fire radiative
  # power) inside the convex hull polygons
  #########################################################################
  
  # Create blank 375 m raster grid
  blankraster <- rast(vect(geo_area_aea), resolution = 375)
  
  # Rasterize date of burning (all active fires, large cluster, small clusters)
  fire_doy_all <- terra::rasterize(vect(fire_clip), blankraster, field = "fdoy")
  fire_doy_lg <- terra::rasterize(vect(fires_lg), blankraster, field = "fdoy")
  fire_doy_sm <- terra::rasterize(vect(fires_sm), blankraster, field = "fdoy")
  # Rasterize FRP (all active fires, large cluster, small clusters)
  fire_frp_all <- terra::rasterize(vect(fire_clip), blankraster, field = "frp")
  fire_frp_lg <- terra::rasterize(vect(fires_lg), blankraster, field = "frp")
  fire_frp_sm <- terra::rasterize(vect(fires_sm), blankraster, field = "frp")
  
  # Dataset for intepolation - need to do this for gstat() to work with the
  # terra::interpolate() function
  fires_lg_xy <- fires_lg %>%
    mutate(x = sf::st_coordinates(.)[,1],
           y = sf::st_coordinates(.)[,2]) %>%
    data.frame()
  
  # Nearest neighbor interpolation for burn date
  m1doy <- gstat(formula=fdoy~1, locations=~x+y, data=fires_lg_xy, 
                 nmax=1, set=list(idp = 0))
  
  # Nearest neighbor interpolation for FRP
  m1frp <- gstat(formula=frp~1, locations=~x+y, data=fires_lg_xy, 
                 nmax=1, set=list(idp = 0))
  
  # Generate interpolated grid for burn date and FRP
  doy_interp <- interpolate(blankraster, m1doy, debug.level=0)
  frp_interp <- interpolate(blankraster, m1frp, debug.level=0)

  # NOTES: add code for cross-validation of gstat() interpolation models
  # Also try IDW (or multiple nearest neighbors?) for FRP
  
  # Rasterize the burned area polygons
  poly_rast <- terra::rasterize(vect(fire_poly), blankraster, touches = FALSE)
  # Inside the polygons, use interpolated values where there are missing data,
  # otherwise use the active fire observations
  doy_comb <- ifel(is.na(fire_doy_all) & poly_rast == 1 & not.na(poly_rast), 
                   doy_interp[[1]], 
                   fire_doy_lg)
  frp_comb <- ifel(is.na(fire_frp_all) & poly_rast == 1 & not.na(poly_rast), 
                   frp_interp[[1]], 
                   fire_frp_lg)
  # Put the small fires back into the dataset
  doy_comb <- merge(doy_comb, fire_doy_sm)
  frp_comp <- merge(frp_comb, fire_frp_sm)
  
  #########################################################################
  # Step 2.4 Output the burned area polygons and raster data
  #########################################################################  
  
  outdoy <- paste0(geo_name, "_activedoy_", curyear, ".tif")
  outdoy2 <- paste0(geo_name, "_interpdoy_", curyear, ".tif")
  outfrp <- paste0(geo_name, "_activefrp_", curyear, ".tif")
  outfrp2 <- paste0(geo_name, "_interpfrp_", curyear, ".tif")
  outpoly <- paste0(geo_name, "_polyfires_", curyear, ".shp")
  
  # Interpolated burn dates
  writeRaster(doy_comb, outdoy2, overwrite = T)
  # Observed burn dates
  writeRaster(fire_doy_all, outdoy, overwrite = T)
  # Interpolated FRP
  writeRaster(frp_comb, outfrp2, overwrite = T)
  # Observed FRP
  writeRaster(fire_frp_all, outfrp, overwrite = T)
  # Burned area polygons
  st_write(fire_poly, outpoly, append = FALSE)
  
}

