#' Run the species distribution model: currently runs for Maxent only. Will add in other models later.

#' 1. Load the required packages & set path

library(flexsdm)
library(terra)
library(data.table)
library(parallel)

path <- "Z:/pops_pesthostuse/pops.sdm/Data/"
setwd("D:/blaginh/sdm_truncated_extent/")
domain <- "USA"
extent <- "D:/blaginh/sdm_truncated_extent/flexsdm_results/1_Inputs/3_Calibration_area/studyext.gpkg"
res <- 30

#' 2. Set up SDM directories
flexsdm::sdm_directory(
  algorithm = TRUE
)


#' 3. Create base raster for the study extent
source("C:/Users/blaginh/Documents/GitHub/pops.sdm/pops.sdm/sdm_helpers.R")
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/get_envi_chunked.R") # get_topo_global # nolint
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/raster_base.R") # base_raster at correct resolution and extent # nolint

res <- fix_resolution(res, domain)
base <- cropped_base_raster(domain, res, path, extent)

#' Write out the base raster
cropped_base_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/3_Calibration_area/base_study_extent.tif") # nolint
writeRaster(base, cropped_base_path, overwrite = TRUE,
            gdal = "COMPRESS=ZSTD", datatype = "INT1U")

#' 3. Copy predictors and response data to the appropriate directories
## predictors <- subset_rasters(path, domain = domain, all = FALSE, pop = FALSE, roads = FALSE, rails = FALSE, gdd = FALSE)

## files <- as.vector(unlist(lapply(predictors, function(x) {
##   get_filename(x, "USA", path)
## })))

#' Retain only the rasters at the correct resolution
## subset_files <- files[grep(paste0(res, "m"), files)]

# Create a directory for predictors
## pred_copies_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/copies/") # nolint
pred_cropped_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/") # nolint
## dir.create(pred_copies_path, showWarnings = TRUE, recursive = TRUE)
## dir.create(pred_cropped_path, showWarnings = TRUE, recursive = TRUE)

# Copy and crop the predictors
## cropped_predictors <- list()

## for (i in seq_along(subset_files)) {
##   cropped_predictors[[i]] <- crop_predictor(
##     file = subset_files[i],
##     pred_copies_path = pred_copies_path,
##     pred_cropped_path = pred_cropped_path,
##     extent = base
##   )
## }

cropped_predictors <- list.files(pred_cropped_path, full.names = TRUE, pattern = "tif", recursive = TRUE)
cropped_predictors <- terra::rast(cropped_predictors)

#' 4. Load the species occurrence data
species <- "Ailanthus altissima"
species <- format_species_name(species)
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/get_pts_v2.r")
occs_pa <- batch_get_pts(species, "species", path, conus = TRUE)
occs_pa <- prep_occurrences(occs_pa, base, year = 1990)
occs_pa <- export_occurrences(occs_pa)

#' 5. Geographic filtering of occurrence data to address spatial autocorrelation
cellsizes <- get_geo_cellsizes(domain = "USA", res = 30)
filt_geo <- list()
for (i in seq_along(cellsizes)) {
  filt_geo[[i]] <- geo_filter_occs(occs_pa, base, cellsizes[i],
                                   species)
}



#' 6. Create a spatial-block partition for the data based on continuous predictors

#' pull out the continuous predictors
#' use datatype to determine which predictors are continuous
i <- 1

env_layer <- list()
for (i in 1:nlyr(cropped_predictors)) {
  dt <- terra::datatype(cropped_predictors[[i]])
  # If datatype is FLT4S, then place in env_layers, otherwise skip
  if (dt == "FLT4S" | dt == "INT2U" | dt == "INT2S") {
    env_layer[[i]] <- cropped_predictors[[i]]
  } else {
  }
}

env_layer <- terra::rast(env_layer)

no_parts <- 5
s <- Sys.time()
sp_part3 <- part_sblock(
  env_layer = env_layer,
  data = filt_geo[[11]],
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  min_res_mult = 2, # Minimum value used for multiplying raster resolution and define the finest resolution to be tested
  max_res_mult = 300, # Maximum value used for multiplying raster resolution and define the coarsest resolution to be tested
  num_grids = 80, # Number of grid to be tested between min_res_mult X (raster resolution) and max_res_mult X (raster resolution)
  n_part = no_parts, # Number of partitions
  prop = 1, # Proportion of points used for testing autocorrelation between groups (0-1)
  min_occ = 10 # Minimum number of occurrences to be used in each partition
)

sp_part4 <- part_senv(
  env_layer = cropped_predictors,
  data = filt_geo[[1]],
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  min_n_groups = 2,
  max_n_groups = nlyr(cropped_predictors),
  prop = 1,
  min_occ = floor(length(filt_geo$fkey)/(no_parts*1.1))
)

plot(regions, col = gray.colors(9))
points(sp_part4$part[c("x", "y")],
  col = hcl.colors(length(unique(filt_geo[[1]]$x)))[sp_part4$part$.part],
  cex = 1,
  pch = 19
)

## grid_env <- get_block(env_layer = base_raster, best_grid = sp_part3$grid)

## # Write out the spatial block partition
blockcv_path <- "D:/blaginh/new_sdm/blockCV/"
## dir.create(paste0(blockcv_path,"/truncated"), showWarnings = FALSE, recursive = TRUE)
## writeRaster(grid_env, paste0(blockcv_path, "truncated/", species, "_spatialblock5.tif"), overwrite = TRUE)
# Read in the spatial block partition
grid_env <- rast(paste0(blockcv_path, "truncated/", species, "_spatialblock5.tif"))

#' 7. Sample background points (using tidysdm because it is faster than flexsdm)
occs_p <- fread(input = paste0(response_path, "truncated/final/", species, "_p.csv"))

#' Random points
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = "random",
    rlayer = grid_env,
    maskval = i)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_random.csv"),
                row.names = FALSE)

#' Buffer points
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = c("thickening", width = 10000),
    rlayer = grid_env,
    maskval = i)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_thickening.csv"),
                row.names = FALSE)

## # Biased: Population density
bias <- rast(bias_files[3])
base_raster_part <- list()
bias_raster_part <- list()
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  base_raster_part[[i]] <- ifel(grid_env == i, 1, NA)
  bias_raster_part[[i]] <- ifel(grid_env == i, bias, NA)
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = "biased",
    rlayer = base_raster_part[[i]],
    rbias = bias_raster_part[[i]],
    maskval = 1)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_thickening.csv"),
                row.names = FALSE)


## # Biased: Target Group Density
bias <- rast(bias_files[grep("target", bias_files)])
base_raster_part <- list()
bias_raster_part <- list()
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  base_raster_part[[i]] <- ifel(grid_env == i, 1, NA)
  bias_raster_part[[i]] <- ifel(grid_env == i, bias, NA)
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = "biased",
    rlayer = base_raster_part[[i]],
    rbias = bias_raster_part[[i]],
    maskval = 1)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_target_group.csv"),
                row.names = FALSE)


## # Biased: Roads
bias <- rast(bias_files[grep("roads", bias_files)])
base_raster_part <- list()
bias_raster_part <- list()
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  base_raster_part[[i]] <- ifel(grid_env == i, 1, NA)
  bias_raster_part[[i]] <- ifel(grid_env == i, bias, NA)
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = "biased",
    rlayer = base_raster_part[[i]],
    rbias = bias_raster_part[[i]],
    maskval = 1)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_roads.csv"),
                row.names = FALSE)

## # Biased: Rails
bias <- rast(bias_files[grep("rails", bias_files)])
base_raster_part <- list()
bias_raster_part <- list()
pts_part <- list()
p_a_part <- list()

for (i in 1:no_parts) {
  base_raster_part[[i]] <- ifel(grid_env == i, 1, NA)
  bias_raster_part[[i]] <- ifel(grid_env == i, bias, NA)
  pts_part[[i]] <- occs_p[occs_p$.part == i, c("x", "y")]
  p_a_part[[i]] <- sample_background(
    data = pts_part[[i]],
    x = "x",
    y = "y",
    n  = length(pts_part[[i]]$x),
    method = "biased",
    rlayer = base_raster_part[[i]],
    rbias = bias_raster_part[[i]],
    maskval = 1)
  p_a_part[[i]]$.part <- i
}

## # Combine the background points
p_a <- rbindlist(p_a_part)
p_a <- rbindlist(list(p_a, occs_p))

## # Write out full dataset
fwrite(p_a, paste0(response_path, "truncated/final/", species, "_p_a_rails.csv"),
                row.names = FALSE)

# Read in the full dataset

#' Read in the bias files
bias_path <- "D:/blaginh/new_sdm/bias_files/"
bias_files <- list.files(bias_path, full.names = TRUE, pattern = "tif", recursive = TRUE)
# Remove "raw"
bias_files <- bias_files[!grepl("raw", bias_files)]

#' 9. Conduct cluster analysis to remove collinearity
source("create_clusters.R")
clusters <- cluster_analysis(somevar_cont, 20000, 0.6, somevar_cat)
combos <- generate_combinations(clusters$var, clusters$cluster)
predictors <- extract_predictors(combos)
length(predictors)

#' 10. Fit maxent models across all predictor sets
ncores <- (detectCores()/2) - 1
makeCluster(ncores)
model_name <- paste0("uniform_", species, "_")
random_pts <-fread(paste0(response_path, "truncated/final/", species, "_p_a_random.csv")
)




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
