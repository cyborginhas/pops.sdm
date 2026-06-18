#' Run the species distribution model: currently runs for Maxent only. Will add in other models later.

#' 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)

path <- "Z:/pops_pesthostuse/pops.sdm/Data/"
setwd("D:/blaginh/sdm_full_extent/")
domain <- "USA"
extent <- "D:/blaginh/sdm_full_extent/misc/ext25.gpkg"
res <- 30

#' 2. Set up SDM directories
flexsdm::sdm_directory(
  algorithm = TRUE
)

#' 3. Create base raster for the study extent
source("C:/Users/blaginh/Documents/GitHub/pops.sdm/pops.sdm/sdm_helpers.R")
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/get_envi_chunked.R") # get_topo_global # nolint
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/raster_base.R") # base_raster at correct resolution and extent # nolint
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/part_sblock.R") # nolint
source("C:/Users/blaginh/Documents/Github/pops.sdm/pops.sdm/sample_background.R") # nolint

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
source("C:/Users/blaginh/Documents/Github/pops.sdm/get_pts_v2.r")
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

# Pull in aggregate data
fn <- subset_files[grep("nlcd", subset_files, invert = TRUE)]
fn <- basename(fn)
fn <- paste0(pred_cropped_path, "/transformed/aggregated/", fn)
fn <- gsub(".tif", "_cropped_zscore_agg.tif", fn)

if (!all(file.exists(fn))) {
  for (i in seq_along(cropped_predictors)) {
    cropped_predictors[[i]] <- aggregate_predictor(cropped_predictors[[i]], factor = 10)
  }
} else {
  agg_predictors <- lapply(fn, rast)
  }

agg_predictors <- rast(agg_predictors)

## #' Reduce number of agg_predictors to consider using cluster
## #' analysis to remove collinearity
## categorical_lyrs <- list()
## for (i in seq_along(cropped_predictors)) {
##   categorical_lyrs[[i]] <- get_categorical_rasters(cropped_predictors[[i]])
## }
## categorical_lyrs <- rast(categorical_lyrs)
## clusters <- cluster_analysis(env_layer, 200000, 0.85, categorical_lyrs)
## clusters <- clusters[var != "NLCD Land Cover Class",]


## #' Keep one predictor from each cluster
## selected_vars <- list()
## clusterids <- unique(clusters$cluster)

## for (i in seq_along(clusterids)) {
##   selected_vars[[i]] <- clusters[cluster == clusterids[i], ]$var[1]
## }
## selected_vars <- unlist(selected_vars)
## # Pull out the selected predictors from agg_predictors
## selected_predictors <- list()
## for (i in seq_along(selected_vars)) {
##   selected_predictors[[i]] <- agg_predictors[[selected_vars[i]]]
## }
## selected_predictors <- rast(selected_predictors)

# Partition the filtered data into spatial blocks & write out
part_dt_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_data_filtgeo", cellsizes, "m2.csv")
part_raster_files <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/",
                            species, "_part_grid_filtgeo", cellsizes, "m2.tif")

partitioned_dt <- list()
partitioned_raster <- list()
partitioned_data <- list()

if(!all(file.exists(part_dt_files))) {
  for (i in 1:length(filt_geo)) {
    partitioned_data[[i]] <- spatial_block_partition(
      filt_geo[[i]], selected_predictors, species, d = cellsizes[i]
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
      partitioned_data[[i]], species, res
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

partitioned_dt <- list()
partitioned_raster <- list()
partitioned_data <- list()

for (i in seq_along(part_dt_files)) {
    partitioned_dt[[i]] <- fread(part_dt_files[i])
    partitioned_raster[[i]] <- rast(part_raster_files[i])
    partitioned_data[[i]] <- list(partitioned_dt[[i]], partitioned_raster[[i]])
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
      partitioned_data[[i]], rbias, rbias_type, species, res,
      extent = base
    )
  }
} else {
  tg_pts <- lapply(tg_files, fread)
}

#' Background points with population density as a bias layer
rbias_type <- "pop_density"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type, res)

pop_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")


if (!all(file.exists(pop_files))) {
  pop_pts <- list()
  for (i in seq_along(partitioned_data)) {
    pop_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species,
      res = res, extent = base
    )
  }
} else {
  pop_pts <- lapply(pop_files, fread)
}

#' Background points with dist2roads as a bias layer
rbias_type <- "distance_roads"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type, res)

roads_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(roads_files))) {
  roads_pts <- list()
  for (i in seq_along(partitioned_data)) {
    roads_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species, 
      res, extent = base
    )
  }
} else {
  roads_pts <- lapply(roads_files, fread)
}

#' Background points with dist2rails as a bias layer
rbias_type <- "distance_rails"
rbias_lyr <- human_factors_rbias_lyr(path, extent = base, domain, rbias_type, res)

rails_files <- paste0(getwd(),
                          "/flexsdm_results/1_Inputs/1_Occurrences/background/", #nolint
                    rbias_type, "/", species, "_", rbias_type,
                    "_bg_pts_filtgeo_", cellsizes, "m2.csv")

if (!all(file.exists(rails_files))) {
  rails_pts <- list()
  for (i in seq_along(partitioned_data)) {
    rails_pts[[i]] <- create_biased_bg_pts(
      partitioned_data[[i]], rbias_lyr, rbias_type, species,
      res, extent = base
    )
  }
} else {
  rails_pts <- lapply(rails_files, fread)
}

#' 9. Conduct cluster analysis to remove collinearity
## categorical_lyrs <- list()
## for (i in seq_along(cropped_predictors)) {
##   categorical_lyrs[[i]] <- get_categorical_rasters(cropped_predictors[[i]])
## }
## categorical_lyrs <- rast(categorical_lyrs)
## clusters <- cluster_analysis(env_layer, 200000, 0.6, categorical_lyrs)
## clusters[var == "NLCD Land Cover Class",]$var <- "landcoverrc"
## combos <- generate_combinations(clusters)
## dir.create(paste0(getwd(), "/flexsdm_results/1_Inputs/3_Clusters/"), showWarnings = FALSE, recursive = TRUE)
## fwrite(clusters, paste0(getwd(), "/flexsdm_results/1_Inputs/3_Clusters/", species, "_clusters.csv"))
## fwrite(combos, paste0(getwd(), "/flexsdm_results/1_Inputs/3_Clusters/", species, "_combos.csv"))

#' 10. Extract predictors for model fitting
bg_method <- c("random", "thicken", "target_group", "pop_density", "distance_roads", "distance_rails")

fn <- vector("list", length(bg_method))
for (i in seq_along(bg_method)) {
  fn[[i]] <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/background/", 
                    bg_method[i], "/", species, "_", bg_method[i], "_bg_pts_filtgeo_",
                    cellsizes, "m2_wpreds.csv")
}
# Collapse the list
fn <- unlist(fn)

random_pts_wdata <- list()
thickened_pts_wdata <- list()
tg_pts_wdata <- list()
pop_pts_wdata <- list()
roads_pts_wdata <- list()
rails_pts_wdata <- list()

cropped_predictors <- rast(cropped_predictors)

#if(!all(file.exists(fn))) {
#for (i in seq_along(cellsizes)) {
  #random_pts_wdata[[i]] <- append_env_data(random_pts[[i]], cropped_predictors, bg_method = "random", species)
  #thickened_pts_wdata[[i]] <- append_env_data(thickened_pts[[i]], cropped_predictors, bg_method = "thicken", species)
  #tg_pts_wdata[[i]] <- append_env_data(tg_pts[[i]], cropped_predictors, bg_method = "target_group", species)
  #pop_pts_wdata[[i]] <- append_env_data(pop_pts[[i]], cropped_predictors, bg_method = "pop_density", species)
  #roads_pts_wdata[[i]] <- append_env_data(roads_pts[[i]], cropped_predictors, bg_method = "distance_roads", species)
  #rails_pts_wdata[[i]] <- append_env_data(rails_pts[[i]], cropped_predictors, bg_method = "distance_rails", species)
#}
#} else {
#random_pts_wdata <- lapply(fn[grep("random", fn)], fread)
#thickened_pts_wdata <- lapply(fn[grep("thicken", fn)], fread)
#tg_pts_wdata <- lapply(fn[grep("target_group", fn)], fread)
#pop_pts_wdata <- lapply(fn[grep("pop_density", fn)], fread)
#roads_pts_wdata <- lapply(fn[grep("distance_roads", fn)], fread)
#rails_pts_wdata <- lapply(fn[grep("distance_rails", fn)], fread)
#}

#' 11. Create testing set

testingset <- create_testing_set(path, extent = base, species)

#' plot the testing set
testingset <- testingset[setid == 1]

#' 12. Fit models
#' Tune hyperparameters for maxent
gridtest <-
  expand.grid(
    regmult = seq(0.1, 4.1, 0.5),
    classes = c("l", "lh", "lq", "lqh", "lqp", "lqph", "lqt", "lqth", "lqthp")
  )

#' Prep data for maxent
tg_pts_wdata <- "C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_target_group_bg_pts_filtgeo_210m2_wpreds.csv"
tg_pts_wdata <- fread(tg_pts_wdata)
occs <- splitpts4sdm(tg_pts_wdata)

#Combos with categorical predictors
combos <- fread("C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_combos.csv")
combos <- data.table::as.data.table(combos)
colnames(combos) <- paste0("cluster_", 1:12)
combos_f <- combos[!is.na(cluster_12),]
combos_num <- combos[is.na(cluster_12),]

# Create outpath
model_outpath <- paste0(getwd(), "/flexsdm_results/2_Outputs/0_Model_performance/")
dir.create(paste0(model_outpath, "/failed/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(model_outpath, "/completed/"), showWarnings = FALSE, recursive = TRUE)


## max_t1 <- flexsdm::tune_max(
##   data = occs$response,
##   response = "pr_ab",
##   predictors = c(
##     "elevation", "par", "silt_0_5", "bio19",
##     "ph_0_5", "theta_s_0_5"
##   ),
##   predictors_f = c("landcoverrc"),
##   background = occs$bg,
##   partition = ".part",
##   grid = gridtest,
##   thr = c("max_sens_spec"),
##   metric = "TSS",
##   clamp = TRUE,
##   pred_type = "cloglog",
##   n_cores = detectCores()/2
## )
## saveRDS(max_t1, "C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_maxent.rds")
## print(paste0("Model complete"))

best_predictors <- list.files("C:/Users/blaginh/Desktop/maximize/", pattern = ".tif", full.names = TRUE)
best_predictors <- terra::rast(best_predictors)
names(best_predictors)[names(best_predictors) == "NLCD Land Cover Class"] <- "landcoverrc"
max_t1 <- readRDS("C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_maxent.rds")
#summary_predictors(best_predictors, max_t1)
str(max_t1$model)

agg_base <- terra::aggregate(base, fact = 4000)
# Set everything to 1 

dir.create("C:/Users/blaginh/Desktop/maximize/tiles/tiles/", showWarnings = FALSE, recursive = TRUE)
tiles <- terra::makeTiles(best_predictors, agg_base,
                   "C:/Users/blaginh/Desktop/maximize/tiles/tiles/bestpredictors_.tif", na.rm = TRUE, 
                   wopt = list(gdal = "COMPRESS=ZSTD"))

# Read in one tile
tile <- terra::rast("C:/Users/blaginh/Desktop/maximize/tiles/tiles/bestpredictors_1.tif")

dir.create("C:/Users/blaginh/Desktop/maximize/output/", showWarnings = FALSE, recursive = TRUE)


tiles <- list.files("C:/Users/blaginh/Desktop/maximize/tiles/tiles/", pattern = ".tif", full.names = TRUE, na.rm = TRUE)
tiles <- lapply(tiles, terra::rast)
basename(sources(tiles[[1]]))

# a list of models
for (i in seq_along(tiles)) {
    s <- Sys.time()
    ind_p <- sdm_predict(
    nchunk = 1,
    models = max_t1,
    pred = tiles[[i]],
    thr = c("max_sens_spec"),
    con_thr = FALSE,
    predict_area = NULL
  )
  fn <- basename(sources(tiles[[i]]))
  writeRaster(ind_p$max, paste0("C:/Users/Desktop/maximize/output/Ailanthus_altissima_prediction_",fn))
  t <- Sys.time() - s
  print(t)
  tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
}
