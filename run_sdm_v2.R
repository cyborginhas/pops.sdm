#' Load the required packages
library(flexsdm)
library(terra)
library(spatialEco)
library(dplyr)
library(sf)
library(data.table)
library(foreach)
library(doParallel)

path <- "~/Desktop/pops_pesthostuse/"
species <- "Acer rubrum"
sp_region <- vect(paste0(path, "/studyext/ext25_states_mask.gpkg"))

format_species_name <- function(species) {
  # Replace " " with "_" in species name
  species <- gsub(" ", "_", species)
  # ensure the first letter is in uppercase and the rest are in lowercase
  species <- tolower(species)
  # make first letter uppercase
  species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
  return(species)
}

species <- format_species_name(species)

#' 1a. import environmental data
var_paths <- list.files(paste0(path, "/predictors"),
  recursive = TRUE, full.names = TRUE, pattern = "1000m"
)

#' remove bio1|bio2|bio3|bio4|bio6|bio7|bio12|bio14|bio15|bio19
var_paths <- var_paths[!grepl("bio1_|bio2_|bio3_|bio4_|bio6_|bio7_|bio12|bio14|bio15|bio19", var_paths)]
var_paths <- var_paths[!grepl("built|cropland|decid|everg|grassl|pasture|shrubs|wetland", var_paths)]
vars <- rast(var_paths)
names(vars)[1:3] <- c("soils_pH", "soils_1500kPa", "soils_33kPa")

#' 2a. import species occurence data
occs_paths <- list.files(paste0(path, "Table/Global"),
  pattern = species, recursive = TRUE,
  full.names = TRUE
)

occs <- lapply(occs_paths, function(x) {
  dt <- fread(x)
  # convert all columns to character
  dt <- dt[, lapply(.SD, as.character)]
  dt <- unique(dt)
})

occs <- rbindlist(occs)

# Convert date to date format, p_a, lat, lon to numeric
occs[, c("date", "p_a", "lat", "lon") := lapply(.SD, as.numeric),
  .SDcols = c("date", "p_a", "lat", "lon")
]

# Drop rows with NA in lat and lon
occs <- occs[!is.na(lat) & !is.na(lon)]
occs <- occs[p_a == 1]

# Convert to terra object; retain lat, lon, p_a
occs <- unique(occs[, .(x = lon, y = lat, pr_ab = p_a)])

# 3. delimit calibration area using buffer method
ca <- calib_area(
  data = occs,
  x = "x",
  y = "y",
  method = c("buffer", width = 100000),
  crs = crs(vars[[1]])
)

#' Crop the calibration area to fall within the sp_region outline
ca <- terra::intersect(ca, sp_region)
writeVector(ca, paste0(path, "/hostmap/", species, "_calib_area.gpkg"))
# Convert to terra object; retain lat, lon, p_a
occs.pts <- vect(occs, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")
occs.pts <- cbind(occs.pts, occs[, .(x, y)])

# Keep all records within the study extent polygon
occs.pts <- terra::intersect(occs.pts, ca)

# Create grid polygons using sf
bbox <- st_bbox(occs.pts, crs = st_crs(occs.pts))
bbox <- st_as_sfc(bbox)



occs <- as.data.table(occs.pts)

#' 4. Filter occurrence data by cell size
occs_f <- occfilt_geo(occs, "x", "y",
  env_layer = vars[[1]],
  method = c("cellsize", factor = "1")
)

#' 5. Spatial block partitioning
occ_part <- occs_f %>%
  part_sblock(
    data = .,
    env_layer = vars[[1]],
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    n_part = 4,
    min_res_mult = 3,
    max_res_mult = 200,
    num_grids = 10,
    prop = 1
  )

occs_f <- occ_part$part
occ_part$best_part_info
block_layer <- get_block(env_layer = vars, best_grid = occ_part$grid)

# 6. Pseudoabsence sampling
psa <- lapply(1:4, function(x) {
  sample_pseudoabs(
    data = occs_f,
    x = "x",
    y = "y",
    n = sum(occs_f$.part == x),
    method = "random",
    rlayer = block_layer,
    maskval = x,
    calibarea = ca
  )
}) %>%
  bind_rows()

psa <- sdm_extract(data = psa, x = "x", y = "y", env_layer = block_layer)
occs_pa <- bind_rows(occs_f, psa)

# 7. Extract environmental data
occs_pa2 <-
  sdm_extract(
    data = occs_pa,
    x = "x",
    y = "y",
    env_layer = vars,
  )

saveRDS(occs_pa2, paste0(path, "/hostmap/", species, "_occs_pa2.RData"))

# 8. Run create_clusters.R
# exclude landcoverrc
vars2 <- vars[[!grepl("landcoverrc", names(vars))]]
cluster_dt <- cluster_analysis(vars2, sample_size = ncell(vars[[1]]), mincor = 0.5)
combinations <- generate_combinations(cluster_dt$var, cluster_dt$cluster)
predictors <- extract_predictors(combinations)
saveRDS(predictors, paste0(path, "/hostmap/", species, "_predictors.RData"))

# 8. Fit models - GLM, GBM, SVM with max sens/spec threshold
# How many cores does your CPU have
n_cores <- detectCores()
n_cores
# Register cluster
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

s <- Sys.time()
mglm <- foreach(i = 1:length(predictors), .packages = c("glmnet", "caret")) %dopar% {
  flexsdm::fit_glm(
    data = occs_pa2,
    response = "pr_ab",
    predictors = predictors[[i]],
    partition = ".part",
    thr = "max_sens_spec" # same as max TSS
    )
}
t <- Sys.time() - s
stopCluster(cl = cluster)
print(t)

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
