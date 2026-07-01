# Wildfire Encroachment into African Moist Tropical Forests

**DSA 5900 Professional Practicum | University of Oklahoma | Spring 2026**

**Team:** Soumya Sri Kesani, Poorna Chandan Reddy Pandem  
**Sponsor:** Dr. Michael Wimberly and Dr. Gopichandh Danala, Data Institute for Societal Challenges (DISC)  
**Supervisor:** Dr. Matt J. Beattie

---

## Overview

Wildfires burning beneath the closed canopy of West Africa's moist tropical forests are largely invisible to standard global monitoring tools. This project builds a spatiotemporal detection and mapping pipeline to close that discovery gap across the **East-West Guinean Lowland Forests** ecoregion.

Using NASA VIIRS SNPP 375m active fire detections (2012-2024), we engineer a graph-based clustering engine that groups raw thermal hotspots into contiguous fire events. These events are converted into burned area polygons and interpolated onto 375m raster surfaces of burn date and Fire Radiative Power (FRP). The pipeline is benchmarked against MODIS MCD64A1 and validated against Landsat-derived dNBR masks from Google Earth Engine.

---

## Key Results

| Metric | Value |
|--------|-------|
| Study Period | 2012-2024 |
| Total Detections Processed | ~1.85 million |
| Optimal Config | D = 2000m, T = 2 days, NH confidence |
| Mean Kappa (6 validation sites) | 0.588 |
| Mean Precision | 0.856 |
| Discovery Gain over MODIS | 14.1% (128,067 Ha) |
| Active Encroachment Frontiers | Tinte Bepo (p=0.039), Opro River (p=0.047) |
| Peak Fire Year | 2015 |
| Lowest Fire Year | 2021 |

---

## Repository Structure

```
africa-wildfire-viirs/
│
├── process_africa_vaf.R          # Data ingestion from NASA FIRMS
├── Read_africa_vaf.R             # Data reading and rasterization
├── viirs_fire_interpolation.R    # Core clustering and interpolation pipeline
├── MK.R                          # Mann-Kendall trend analysis
├── seasonality.R                 # STL decomposition and seasonal analysis
├── burnvalidation.R              # dNBR-based validation
├── summay_analysis2k2.R          # Summary statistics
│
├── Africa_country_table.csv      # African country reference table
├── EW_Guinean_Lowland_Forests.*  # Ecoregion boundary shapefile
│
├── *.png                         # Output maps and figures
└── DSA_5900_Group_6_Final_Report_with_Changes_Compressed.pdf
```

---

## Methodology

1. **Data Ingestion** - Annual VIIRS SNPP 375m active fire CSVs downloaded from NASA FIRMS and clipped to the ecoregion boundary.

2. **Spatiotemporal Clustering** - Delaunay triangulation via `tri2nb` creates a sparse neighbor graph. Edges are retained only if spatial separation < 2000m and temporal separation <= 2 days. Connected components define discrete fire events.

3. **Polygon Construction** - Clusters with 3 or more detections are converted to concave hull polygons using `concaveman` (concavity = 2) to approximate irregular fire perimeters.

4. **Interpolation** - Nearest-neighbor interpolation via `gstat` fills undetected pixels inside event polygons, producing continuous 375m rasters of burn date and FRP.

5. **Validation** - Six Ghana forest reserve sites validated against Landsat dNBR masks extracted in Google Earth Engine. Metrics: Kappa, Precision, Recall, F1, Omission, and Commission.

6. **Trend Analysis** - Non-parametric Mann-Kendall tests applied to 12 Ghana reserves to identify monotonic encroachment trends.

---

## Dependencies

All analysis was conducted in **R**. Key packages:

- `terra`, `sf` - spatial data handling
- `spdep`, `igraph` - graph-based clustering
- `concaveman` - concave hull generation
- `gstat` - spatial interpolation
- `Kendall`, `trend` - Mann-Kendall testing
- `tidyverse`, `lubridate` - data wrangling
- `ggplot2`, `patchwork` - visualization

---

## Data Sources

- [NASA FIRMS VIIRS SNPP 375m Active Fire Archive](https://firms.modaps.eosdis.nasa.gov/)
- [MODIS MCD64A1 Burned Area Product](https://lpdaac.usgs.gov/)
- Landsat 8/9 Collection 2 L2 via Google Earth Engine
- WWF Terrestrial Ecoregions (Olson et al., 2001)
- Ghana Forest Reserve Boundaries

> Raw fire CSV files (~1.85M rows) are not included in this repository due to size. Download directly from NASA FIRMS.

---

## References

- Schroeder, W., et al. (2014). The new VIIRS 375m active fire detection data product. *Remote Sensing of Environment*, 143, 85-96.
- Wimberly, M.C., et al. (2024). Increasing fire activity in African tropical forests is associated with deforestation and climate change. *Geophysical Research Letters*, 50, e2023GL106240.
- Olson, D.M., et al. (2001). Terrestrial ecoregions of the world. *BioScience*, 51(11), 933-938.
