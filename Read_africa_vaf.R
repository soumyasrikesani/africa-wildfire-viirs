library(raster)
library(tidyverse)
library(sf)

curdir <- "e:/work/projects/VIIRS_fire/downloads"
tempdir <- "e:/work/projects/VIIRS_fire/temp"

africa_data <- read.csv("Africa_country_table.csv")

zipfiles <- list.files(path = curdir,
                       pattern = "*.zip",
                       full.names = TRUE)

years <- 2012:2021
for(i in 1:length(zipfiles)) {
  curzip <- zipfiles[i]
  curyear <- years[i]
  africa_files <- paste0("viirs-snpp/", curyear, "/viirs-snpp_", 
                         curyear, "_", africa_data[,1], ".csv")
  unzip(curzip, 
        files = africa_files, 
        exdir = tempdir,
        junkpaths = TRUE)
  csvfile <- list.files(path = tempdir,
                        pattern = "*.csv",
                        full.names = TRUE)
  for(j in 1:length(csvfile)) {
    print(csvfile[j])
    incsv <- read.csv(csvfile[j])
    if(j == 1) {
      outcsv <- incsv
    } else {
      outcsv <- bind_rows(outcsv, incsv)
    }
  }
  outname <- paste0("africa_fires_", curyear, ".csv")
  write.csv(outcsv, file = outname, row.names=F)
  unlink(csvfile)
}

africa_bnd <- st_read("ESRI_Africa_countries.shp")
africa_proj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
africa_aea <- st_transform(africa_bnd, crs=africa_proj)
africa_grid <- raster(ext=extent(africa_aea), resolution=10000, crs=africa_proj)
africa_mask <- rasterize(africa_aea, africa_grid)

years <- 2012:2021
for(i in 1:length(years)) {
  firefile <- paste0("africa_fires_", years[i], ".csv")
  print(firefile)
  annfires <- read.csv(file.path(".", "output", firefile))
  annfire_pt <- st_as_sf(annfires, coords = c("longitude", "latitude"), crs = 4326)
  annfire_aea <- st_transform(annfire_pt, crs = africa_proj)

  annfire_aea <- annfire_aea %>%
    filter(satellite == "N",
           confidence == "h" | confidence == "n")
  
  africa_fd <- rasterize(annfire_aea, africa_grid, field='satellite',  fun='count')
  africa_frp <- rasterize(annfire_aea, africa_grid, field='frp',  fun=mean)
  africa_fd[is.na(africa_fd)] <- 0
  #africa_frp[is.na(africa_fd)] <- 0
  africa_fd <- mask(africa_fd, africa_mask)
  africa_frp <- mask(africa_frp, africa_mask)
  if(i == 1) {
    annual_fd <- africa_fd
    annual_frp <- africa_frp
  } else {
    annual_fd <- addLayer(annual_fd, africa_fd)
    annual_frp <- addLayer(annual_frp, africa_frp)
  }
}
writeRaster(annual_fd, "africa_viirs_fire_annual.tif", overwrite = TRUE)

writeRaster(annual_frp, "africa_viirs_frp_annual.tif", overwrite = TRUE)

annual_frp_msk <- annual_frp
annual_frp_msk[annual_fd < 50] <- NA
annual_frp_mean <- mean(annual_frp_msk[[2:10]])
annual_frp_anom <- annual_frp_msk - annual_frp_mean

writeRaster(annual_frp_anom, "africa_viirs_frp_anomalies.tif", overwrite = TRUE)

annual_fm <- mean(log(annual_fd[[2:10]] + 1))
annual_fa <- log(annual_fd[[2:10]] + 1) - annual_fm
writeRaster(annual_fa, "africa_viirs_fire_anomalies.tif")
annual_fa2 <- annual_fa
annual_fa2[annual_fa2 > 3] <- 3
annual_fa2[annual_fa2 < -3] <- -3

#ugcountries <- c("Ghana", "Côte d'Ivoire", "Liberia", "Sierra Leone", "Guinea",
#                 "Benin", "Togo", "Nigeria")
ugcountries <- c("Congo DRC", "Congo", "Gabon", "Cameroon")
ugregion <- africa_aea %>%
  filter(COUNTRY %in% ugcountries)

annual_fa3 <- crop(annual_fa2, ugregion)

protected <- st_read("africa_protected_areas.shp")
protected_aea <- st_transform(protected, crs(africa_aea))
protected_ug <- st_crop(protected_aea, ugregion)

annfireanom <- rasterdf(annual_fa3[[8]])
ggplot() +
  geom_raster(data = annfireanom, aes(x = x, y = y, fill = value)) +
  geom_sf(data = protected_ug, color = "black", fill = NA, size = 0.5) +
  scale_fill_gradient2(name = "Degrees C", low = "blue", mid = "lightyellow", 
                       high = "red") +  coord_sf(expand = F) +
  #facet_wrap(~ variable, ncol = 3) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())




rasterdf <- function(x, aggregate = 1) {
  resampleFactor <- aggregate        
  inputRaster <- x    
  inCols <- ncol(inputRaster)
  inRows <- nrow(inputRaster)
  # Compute numbers of columns and rows in the new raster for mapping
  resampledRaster <- raster(ncol=(inCols / resampleFactor), 
                            nrow=(inRows / resampleFactor))
  # Match to the extent of the original raster
  extent(resampledRaster) <- extent(inputRaster)
  # Resample data on the new raster
  y <- resample(inputRaster,resampledRaster,method='ngb')
  
  # Extract cell coordinates into a data frame
  coords <- xyFromCell(y, seq_len(ncell(y)))
  # Extract layer names
  dat <- stack(as.data.frame(getValues(y)))
  # Add names - 'value' for data, 'variable' to indicate different raster layers
  # in a stack
  names(dat) <- c('value', 'variable')
  dat <- cbind(coords, dat)
  dat
}
