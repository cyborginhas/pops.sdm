# Updating pops.sdm version 1 to include the following in version 2:

1. Facilitate debugging (chunking existing code into R functions throughout) X
2. Dynamic base raster creation from elevation (conus 30m SRTM & worldclim elevation 30s), storage to research drive and metadata creation. X
3. Helper functions to confirm if base rasters already exist at requested resolution; resolutions below 30m (conus) and 1000m (elevation), but trigger warnings to indicate bilinear/near resampling.
4. Dynamic retrieval of existing predictor variables, storage to research drive, and metadata creation.
5. Updated to include dynamic retrieval of several other conus fine-extent predictor variables (bioclimatic, ndvi, evi, dist2rails, dist2roads, population, PAR, soils, gdd, precip_timing). 
7. Increased flexibility in study analysis (domain) cropping; removed vector domain scripts at specific geopolitical boundaries, and replaced with user-provided extent as terra SpatVect object; cropping occurs
occurs after predictors have been processed (downloaded & backed up - one time; re-processed if need variable needed/subsetting - one time; reprojected and resampled - one time for specific common extents e.g. 30m, 
100m, 250m, 500m and dynamically to accomodate other resolutions). 
8. Cluster analysis added; and can be called to dynamically updated the clusters considered (optimal clusters depend on predictors considered, their resolution, etc.); Cluster
analysis also exports a list of predictor combinations for get_BestVars2.R script. X
9. Dynamic retrieval, storage, and metadata creation of host points by species name, additional databases considered (neon), and scripts to overcome api limits (gbif, rinat) X

Low priority:

10. Variogram analysis script added (optional); can be used to determine an appropriate resolution given the data considered; Helper function will throw a warning if resolution selected is not between
nugget and range for the predictor(s) selected.
11. Species-environment plots added (optional); can be used to help determine which algorithms are appropriated.
