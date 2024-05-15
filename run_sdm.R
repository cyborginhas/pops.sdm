#' 1. Load the required packages
library(flexsdm)
library(terra)
library(data.table)
library(dplyr)
library(ks)
library(raster)
library(spatialEco)

#' 2. Load the data and crop it to the appropriate extent
og_path <- "Z:/pops_pesthostuse/pops.sdm/Data/"

## rasters <- get_rasters(path, domain = "USA")
## files <- unlist(lapply(rasters, function(x) {
##   get_filename(x, "USA", path)
## }))

## # Pull out all the rasters that are 30m
## files_30m <- files[grep("30m", files)]
## files_30m <- files_30m[-1]

# Create a directory for predictors
pred_path <- "D:/blaginh/new_sdm/predictors/"
pred_path_fullextent <- "D:/blaginh/new_sdm/predictors/full_extent/"
pred_path_trunc <- "D:/blaginh/new_sdm/predictors/truncated/"
studyext_path <- "D:/blaginh/new_sdm/study_extent/"

## dir.create(pred_path_fullextent, showWarnings = FALSE, recursive = TRUE)
## dir.create(pred_path_trunc, showWarnings = FALSE, recursive = TRUE)
## dir.create(studyext_path, showWarnings = FALSE, recursive = TRUE)

# Copy the study extent to the study extent directory
ext25 <- vect("D:/blaginh/new_sdm/study_extent/ext25.gpkg")
pa <- vect("D:/blaginh/new_sdm/study_extent/pa.gpkg")

# Copy the 30m rasters to the predictors directory; test on one raster
## for (file in files_30m) {
##   s <- Sys.time()
##   # Copy the file to the predictors directory
##   file.copy(file, pred_path)
##   copy <- paste0(pred_path, basename(file))
##   r <- terra::rast(copy)
##   # Replace the path with the new path
##   new_filename <- paste0(pred_path_fullextent, basename(copy))
##   # Replace reproj_30m.tif with 30m_ext25.tif
##   new_filename <- gsub("reproj_30m.tif", "30m_ext25.tif", new_filename)
##   dtype <- terra::datatype(r)
##   # Crop the raster to the study extent
##   r2 <- terra::crop(r, ext25)
##   writeRaster(r2, new_filename, overwrite = TRUE, datatype = dtype, gdal = "COMPRESS=ZSTD")
##   t <- Sys.time() - s
##   print(paste0("Time taken to crop ", basename(new_filename), ": ", t))
##   tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
## }

# Function to get species occurrences
format_species_name <- function(species) {
  # Replace " " with "_" in species name
  species <- gsub(" ", "_", species)
  # ensure the first letter is in uppercase and the rest are in lowercase
  species <- tolower(species)
  # make first letter uppercase
  species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
  return(species)
}

get_species_occurrences <- function(species, path) {
  files <- list.files(paste0(path, "/Table/"), pattern = species, full.names = TRUE, recursive = TRUE)
  occs <- lapply(files, function(x) {
    fread(x, colClasses = "character")
  })
  occs <- rbindlist(occs, fill = TRUE)
  occs[, `:=`(lat = as.numeric(lat), lon = as.numeric(lon), date = as.numeric(date), p_a = as.numeric(p_a))]
  occs <- occs[date > 1990]
  #occs <- occs[p_a > 0]
  return(occs)
}

species <- "Ailanthus altissima"
species <- format_species_name(species)

## occs <- get_species_occurrences(species, og_path)
## # convert to vector with WGS84 CRS using terra vect()
## occs_pts <- vect(occs, crs = crs(ext25), geom = c("lon", "lat"))
## # crop to study extent
## occs_ext25 <- crop(occs_pts, ext25)
## # convert back to data.table
## occs_ext25$lon <- crds(occs_ext25)[, 1]
## occs_ext25$lat <- crds(occs_ext25)[, 2]
## occs_ext25 <- as.data.table(occs_ext25)
## # crop to pa extent
## occs_pa <- crop(occs_pts, pa)
## occs_pa$lon <- crds(occs_pa)[, 1]
## occs_pa$lat <- crds(occs_pa)[, 2]
## occs_pa <- as.data.table(occs_pa)

## # Export the occurrence data to a csv file
response_path <- "D:/blaginh/new_sdm/response/"
## dir.create(paste0(response_path, "full_extent/"),
##   showWarnings = FALSE,
##   recursive = TRUE
## )
## dir.create(paste0(response_path, "truncated/"),
##   showWarnings = FALSE,
##   recursive = TRUE
## )
## # write to file
## write.csv(occs_ext25, paste0(
##   response_path, "full_extent/", species,
##   "_occ_ext25.csv"
## ), row.names = FALSE)
## write.csv(occs_pa, paste0(response_path, "truncated/", species, "_occ_pa.csv"),
##   row.names = FALSE
## )

#' 3. Load the cropped predictors and response data
somevar <- list.files(pred_path_trunc, full.names = TRUE)
# Continuous predictors
somevar_cont <- somevar[!grepl("nlcd|BioClimComposite|wc2.1|dist2|human_pop_density|theta_r|theta_s|hb_0_5|_n_0_5", somevar)]
somevar_cont <- rast(somevar_cont)
# Categorical predictors
somevar_cat <- somevar[grepl("nlcd", somevar)]
somevar_cat <- rast(somevar_cat)
names(somevar_cat) <- "nlcd"

# Bias corrections
somevar_bias <- somevar[grepl("dist2|human_pop_density", somevar)]
somevar_bias <- rast(somevar_bias)

# Response data
sp_region <- pa
species <- "Ailanthus_altissima"
spp1 <- fread(paste0(response_path, "truncated/", species, "_occ_pa.csv"))

# Re-label lon and lat as x and y
setnames(spp1, c("lon", "lat"), c("x", "y"))

#' 4. Preprocess the data: geographic data fitering to offset sampling bias
base_raster <- somevar_cat[[1]]
# Ifel >0 is 1, else 0
base_raster <- ifel(base_raster > 0, 1, 0)
# Turn all na into 0
base_raster <- ifel(is.na(base_raster), 0, base_raster)
factors <- c(1, seq(3, 33, by = 9))
## # Empty list to store the filtered data
## filt_geo <- list()

## s <- Sys.time()
## for (i in seq_along(factors)) {
##   filt_geo[[i]] <- occfilt_geo(
##     data = spp1,
##     x = "x",
##     y = "y",
##     env_layer = base_raster,
##     method = c("cellsize", factor = factors[i]), # coarser resolution than the provided raster
##     prj = crs(base_raster)
##   )
## }
## t <- Sys.time() - s
## print(paste0("Time taken to filter geographic data: ", t))

## # Write out the filtered data
## dir.create(paste0(response_path, "truncated/geofilter/"), showWarnings = FALSE, recursive = TRUE)

## for (i in seq_along(factors)) {
##   write.csv(filt_geo[[i]], paste0(response_path, "truncated/geofilter/", species, "_occ_pa_geofilter_factor", factors[i], ".csv"), row.names = FALSE)
## }

## # Read in filtered data
## filt_geo <- list()
## for (i in factors) {
##   filt_geo[[i]] <- fread(paste0(response_path, "truncated/", species, "_occ_pa_geofilter_factor", i, ".csv"))
## }

#' 6. Create a spatial-block partition for the data based on continuous predictors
no_parts <- 5
## sp_part3 <- part_sblock(
##   env_layer = base_raster,
##   data = filt_geo[[1]],
##   x = "x",
##   y = "y",
##   pr_ab = "p_a",
##   min_res_mult = 10, # Minimum value used for multiplying raster resolution and define the finest resolution to be tested
##   max_res_mult = 1000, # Maximum value used for multiplying raster resolution and define the coarsest resolution to be tested
##   num_grids = 50, # Number of grid to be tested between min_res_mult X (raster resolution) and max_res_mult X (raster resolution)
##   n_part = no_parts, # Number of partitions
##   prop = 0.5, # Proportion of points used for testing autocorrelation between groups (0-1)
##   min_occ = floor(length(filt_geo$fkey)/(no_parts*2)) # Minimum number of occurrences to be used in each partition
## )

## grid_env <- get_block(env_layer = base_raster, best_grid = sp_part3$grid)

## # Write out the spatial block partition
blockcv_path <- "D:/blaginh/new_sdm/blockCV/"
## dir.create(paste0(blockcv_path,"/truncated"), showWarnings = FALSE, recursive = TRUE)
## writeRaster(grid_env, paste0(blockcv_path, "truncated/", species, "_spatialblock5.tif"), overwrite = TRUE)
# Read in the spatial block partition
grid_env <- rast(paste0(blockcv_path, "truncated/", species, "_spatialblock5.tif"))

#' 7. Sample background points (using tidysdm because it is faster than flexsdm)
## base_raster_part <- list()
## pts_part <- list()
## p_a_part <- list()

## for(i in 1:no_parts){
##   base_raster_part[[i]] <- ifel(grid_env == i, 1, NA)
##   pts_part[[i]] <- sp_part3$part[sp_part3$part$.part == i, c("x", "y")]
##   p_a_part[[i]] <- tidysdm::sample_pseudoabs(
##     data = pts_part[[i]],
##     raster = base_raster_part[[i]],
##     method = c('dist_max',10000),
##     n = length(pts_part[[i]]$x),
##     return_pres = TRUE
##   )
##   p_a_part[[i]]$.part <- i
## }

## # Combine the background points
## p_a <- rbindlist(p_a_part)
## names(p_a) <- c("x", "y", "pr_ab", ".part")

## # If pr_ab == "presence", then pr_ab = 1, else 0
## p_a$pr_ab <- ifelse(p_a$pr_ab == "presence", 1, 0)

## # Write out full dataset
## write.csv(p_a, paste0(response_path, "truncated/", species, "_p_a.csv"), row.names = FALSE)

# Read in the full dataset
response_path <- "D:/blaginh/new_sdm/response/"
occs_pa <- fread(paste0(response_path, "truncated/", species, "_p_a.csv"))

#' 8. Create target background layer
#' Read in csv file using fread
target_group_path <- "D:/blaginh/new_sdm/target_group_lyr/"
## target_group_dt <-fread(paste0(target_group_path, "tree_observations_uscanada.csv"))

## ext_vals <- round(ext(pa), 3)


## # Drop data outside the truncated extent
## target_group_dt_trunc <- target_group_dt[
##   target_group_dt$longitude >= ext_vals[1] &
##     target_group_dt$longitude <= ext_vals[2] &
##     target_group_dt$latitude >= ext_vals[3] &
##     target_group_dt$latitude <= ext_vals[4]
## ]

## # Write out the target group data
## dir.create(paste0(target_group_path, "truncated"), showWarnings = FALSE, recursive = TRUE)
## write.csv(target_group_dt_trunc, paste0(target_group_path, "truncated/", species, "_target_group_trunc.csv"), row.names = FALSE)


## # Drop data outside the full extent
## ext_vals <- round(ext(ext25), 3)


## target_group_dt_full <- target_group_dt[
##   target_group_dt$longitude >= ext_vals[1] &
##     target_group_dt$longitude <= ext_vals[2] &
##     target_group_dt$latitude >= ext_vals[3] &
##     target_group_dt$latitude <= ext_vals[4]
## ]

## # Write out the target group data
## dir.create(paste0(target_group_path, "full_extent"), showWarnings = FALSE, recursive = TRUE)
## write.csv(target_group_dt_full, paste0(target_group_path, "full_extent/", species, "_target_group_full.csv"), row.names = FALSE)

## # Remove dups
## tgt_grp_dt_trunc_clean <- unique(target_group_dt_trunc[ ,.(genus, species, latitude, longitude)])

## # Aggregate by coodinates
## sum_records <- tgt_grp_dt_trunc_clean[, .(count = .N), by = .(latitude, longitude)]

## # Create a scale.
## scale <- length(sum_records[, count]) / sum(sum_records[, count])

## sum_records$scale <- sum_records[,count] * scale

## sum_records <- as.matrix(sum_records)

## # Do a 2d kernel density estimation.
## target_density <- kde(sum_records[,c(2,1)], w=sum_records[, 4])


## # Create raster.
## target_raster <- raster(target_density, crs=crs(base_raster))

## # Clip data to the same resolution/extent.
## # If 1, is 1, else 0
## base_raster_na <- raster(ifel(base_raster > 0, 1, NA))
## target_raster <- raster::resample(target_raster, raster(base_raster), method='bilinear')
## target_raster <- mask(target_raster, base_raster_na)

## # Normalize bias file between 0 and 1.
## target_raster <- target_raster - minValue(target_raster)
## target_raster <- raster.transformation(rast(target_raster), trans="norm")

## # Create a file pathway and export the raster.
## writeRaster(target_raster, paste0(target_group_path, "truncated/TargetGroup_biasfile_trunc.tif"),
##             overwrite = TRUE, gdal = "COMPRESS=ZSTD")

## # Read in the target group bias file
target_group_bias <- rast(paste0(target_group_path, "truncated/TargetGroup_biasfile_trunc.tif"))
plot(target_group_bias)

#' 9. Conduct cluster analysis to remove collinearity
clusters <- cluster_analysis(somevar_cont, 20000, 0.6, somevar_cat)
combos <- generate_combinations(clusters$var, clusters$cluster)
predictors <- extract_predictors(combos)
predictors

#' 10. Fit maxent models across all predictor sets



#' 11. Pull out the best model based on TSS (independent testing set)

#' 12. Predict the best model

predictors <- split(predictors, ceiling(seq_along(predictors) / 1000))

# Chunks of predictors
for (j in 6:length(predictors)) {
  predictors2 <- predictors[[j]]
  n_cores <- detectCores()
  cluster <- makeCluster(n_cores - 1)
  registerDoParallel(cluster)

  mglm <- foreach(i = 1:length(predictors2), .packages = c("glmnet", "caret")) %dopar% {
    flexsdm::fit_glm(
      data = occs_pa,
      response = "pr_ab",
      predictors = predictors2[[i]],
      partition = ".part",
      thr = "max_sens_spec" # same as max TSS
    )
  }
  stopCluster(cl = cluster)

  # Pull out models with highest AIC and TSS
  aic <- lapply(mglm, function(x) x$model$aic)
  model_bestaic <- mglm[[which.min(aic)]]
  saveRDS(model_bestaic, paste0(path, "hostmap/", species, "_mglm_aic_", j, ".RData"))

  tss <- lapply(mglm, function(x) x$performance$TSS_mean)
  model_besttss <- mglm[[which.max(tss)]]
  saveRDS(model_besttss, paste0(path, "hostmap/", species, "_mglm_tss_", j, ".RData"))
  beepr::beep(1)
  tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  gc()
  rm(mglm)
}

# Re-do with lc
occs_pa$landcoverrc <- as.factor(occs_pa$landcoverrc)

# Chunks of predictors
for (j in 1:length(predictors)) {
  predictors2 <- predictors[[j]]
  n_cores <- detectCores()
  cluster <- makeCluster(n_cores - 1)
  registerDoParallel(cluster)

  mglm <- foreach(i = 1:length(predictors2), .packages = c("glmnet", "caret")) %dopar% {
    flexsdm::fit_glm(
      data = occs_pa,
      response = "pr_ab",
      predictors = predictors2[[i]],
      predictors_f = c("landcoverrc"),
      partition = ".part",
      thr = "max_sens_spec" # same as max TSS
    )
  }
  stopCluster(cl = cluster)

  # Pull out models with highest AIC and TSS
  aic <- lapply(mglm, function(x) x$model$aic)
  model_bestaic <- mglm[[which.min(aic)]]
  saveRDS(model_bestaic, paste0(path, "hostmap/", species, "lc_mglm_aic_", j, ".RData"))

  tss <- lapply(mglm, function(x) x$performance$TSS_mean)
  model_besttss <- mglm[[which.max(tss)]]
  saveRDS(model_besttss, paste0(path, "hostmap/", species, "lc_mglm_tss_", j, ".RData"))
  beepr::beep(1)
  tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  gc()
  rm(mglm)
}


mraf <- fit_raf(
  occs_pa,
  "pr_ab",
  predictors2[[1]],
  partition = ".part",
  thr = "max_sens_spec"
  # metric = "TSS",
  # grid = tune_grid
)

mraf$performance$TSS_mean


tune_grid <-
  expand.grid(mtry = seq(1, 5, 1))


mglm_sum <- as.data.table(mglm_sum)
mglm_sum[order(-TSS_mean)]


# 8. Ensemble models
eglm <-
  esm_glm(
    data = occs_pa3,
    response = "pr_ab",
    predictors = preds,
    partition = ".part",
    thr = "max_sens_spec"
  )

egbm <- esm_gbm(
  data = occs_pa3,
  response = "pr_ab",
  predictors = preds,
  partition = ".part",
  thr = "max_sens_spec"
)

esvm <- esm_svm(
  data = occs_pa3,
  response = "pr_ab",
  predictors = preds,
  partition = ".part",
  thr = "max_sens_spec"
)

# 9. Predict models
eglm <-
  esm_glm(
    data = hespero_pa3,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part",
    thr = "max_sens_spec"
  )

mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

# 9. Summarize models
merge_df <- sdm_summarize(models = list(mglm, mgbm, msvm))

knitr::kable(
  merge_df %>% dplyr::select(
    model,
    AUC = AUC_mean,
    TSS = TSS_mean,
    JACCARD = JACCARD_mean,
    BOYCE = BOYCE_mean,
    IMAE = IMAE_mean
  )
)

#' Predict models

mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = vars,
  con_thr = TRUE,
  predict_area = ca
)
