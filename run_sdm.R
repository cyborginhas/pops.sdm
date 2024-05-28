#' Run the species distribution model: currently runs for Maxent only. Will add in other models later.

#' 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

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
base <- crop_base_raster(domain, res, path, extent)

#' 4. Copy predictors, crop to the study extent, and transform to z-scores
predictors <- subset_rasters(path, domain = domain, all = FALSE, gdd = FALSE, roads = FALSE, rails = FALSE, pop = FALSE)

files <- as.vector(unlist(lapply(predictors, function(x) {
  get_filename(x, "USA", path)
})))

#' Retain only the rasters at the correct resolution
subset_files <- files[grep(paste0(res, "m"), files)]

# Create a directory for predictors
pred_copies_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/copies/") # nolint
pred_cropped_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/") # nolint
dir.create(pred_copies_path, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_cropped_path, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(pred_cropped_path, "/transformed"), showWarnings = FALSE,
           recursive = TRUE)

# Copy and crop the predictors
cropped_predictors <- list()

for (i in seq_along(subset_files)) {
  cropped_predictors[[i]] <- crop_predictor(
    file = subset_files[i],
    pred_copies_path = pred_copies_path,
    pred_cropped_path = pred_cropped_path,
    extent = base
  )
}

#' 4. Load the species occurrence data
species <- "Ailanthus altissima"
species <- format_species_name(species)
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/get_pts_v2.r")
occs_file <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/original/",
                   species, "_train_occs.csv")

if(!file.exists(occs_file)) {
  occs_pa <- batch_get_pts(species, "species", path, conus = TRUE)
  occs_pa <- prep_occurrences(occs_pa, extent = base, year = 1990)
  occs_pa <- export_occurrences(occs_pa)
} else {
  occs_pa <- fread(occs_file)
}

#' 5. Geographic filtering of occurrence data to address spatial autocorrelation
cellsizes <- get_geo_cellsizes(domain = "USA", res = 30)
filtgeo_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/filtered/",
                        species, "_occ_pa_geofilter_", cellsizes, "m2.csv")

if(!all(file.exists(filtgeo_files))) {
  filt_geo <- list()
  for (i in seq_along(cellsizes)) {
    filt_geo[[i]] <- geo_filter_occs(
    occs_pa, extent = base, d = cellsizes[i],
    species
  )
}
} else {
  filt_geo <- lapply(filtgeo_files, fread)
}

#' 6. Create a spatial-block partition for the data based on continuous predictors

# Isolate numerical predictors for spatial block partitioning
env_layer <- list()

for (i in seq_along(cropped_predictors)) {
  env_layer[[i]] <- get_numerical_rasters(cropped_predictors[[i]])
}

env_layer <- rast(env_layer)

# Partition the filtered data into spatial blocks & write out
part_dt_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_data_filtgeo", cellsizes, "m2.csv")
part_raster_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_grid_filtgeo", cellsizes, "m2.tif")

partitioned_dt <- list()
partitioned_raster <- list()
partitioned_data <- list()

if(!all(file.exists(part_dt_files))) {
  for (i in seq_along(filt_geo)) {
    partitioned_data[[i]] <- spatial_block_partition(
      filt_geo[[i]], extent = base, env_layer, species
    )
  }
} else {
  for (i in seq_along(part_dt_files)) {
    partitioned_dt[[i]] <- fread(part_dt_files[i])
    partitioned_raster[[i]] <- rast(part_raster_files[i])
    partitioned_data[[i]] <- list(partitioned_dt[[i]], partitioned_raster[[i]])
  }
}

#' 7. Make background points using various methods
# Random background points
random_files <- paste0(getwd(),
                       "/flexsdm_results/1_Inputs/1_Occurrences/background/random/", #nolint
                       species, "_bg_pts_filtgeo_", cellsizes, "m2.csv")
random_pts <- list()

if (!all(file.exists(random_files))) {
  for (i in seq_along(partitioned_data)) {
    random_pts[[i]] <- create_random_bg_pts(
      partitioned_data[[i]], species
    )
  }
} else {
  random_pts <- lapply(random_files, fread)
}

#' Thicken background points around occurrence data
thickened_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/thicken/", #nolint
                          species, "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(thickened_files))) {
thickened_pts <- list()
for (i in seq_along(partitioned_data)) {
  thickened_pts[[i]] <- create_thickened_bg_pts(
    partitioned_data[[i]], base, res, species
  )
}
} else {
  thickened_pts <- lapply(thickened_files, fread)
}

#' Background points with target groups as a bias
rbias_type <- "target_group"
rbias_lyr <- target_group_rbias_lyr(path, extent = base, domain, res)
tg_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", rbias_type, "/", #nolint
                          species,"_", rbias_type, "_bg_pts_filtgeo_", cellsizes, "m2.csv")


if (!all(file.exists(tg_files))) {
  tg_pts <- list()
  for (i in seq_along(partitioned_data)) {
    tg_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species
    )
  }
} else {
  tg_pts <- lapply(tg_files, fread)
}

#' Background points with population density as a bias layer
rbias_type <- "pop_density"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type)

pop_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(pop_files))) {
  pop_pts <- list()
  for (i in seq_along(partitioned_data)) {
    pop_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species
    )
  }
} else {
  pop_pts <- lapply(pop_files, fread)
}

#' Background points with dist2roads as a bias layer
rbias_type <- "distance_roads"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type)

roads_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(roads_files))) {
  roads_pts <- list()
  for (i in seq_along(partitioned_data)) {
    roads_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species
    )
  }
} else {
  roads_pts <- lapply(roads_files, fread)
}

#' Background points with dist2rails as a bias layer

rbias_type <- "distance_rails"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type)

rails_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(rails_files))) {
  rails_pts <- list()
  for (i in seq_along(partitioned_data)) {
    rails_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species
    )
  }
} else {
  rails_pts <- lapply(rails_files, fread)
}

#' 9. Conduct cluster analysis to remove collinearity
source("create_clusters.R")
clusters <- cluster_analysis(somevar_cont, 20000, 0.6, somevar_cat)
combos <- generate_combinations(clusters$var, clusters$cluster)
predictors <- extract_predictors(combos)
length(predictors)

#' 10. Fit maxent models across all predictor sets
ncores <- (detectCores() / 2) - 1
makeCluster(ncores)
model_name <- paste0("uniform_", species, "_")
random_pts <- fread(paste0(response_path, "truncated/final/", species, "_p_a_random.csv"))




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
