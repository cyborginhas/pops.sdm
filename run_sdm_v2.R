# Load the required packages
#devtools::install_github('sjevelazco/flexsdm')
library(flexsdm)
library(terra)
library(spatialEco)
library(dplyr)
library(bigmemory)
library(sf)
path <- "/Volumes/cmjone25/pops_pesthostuse/pops.sdm/Data/"
# # 1a. import environmental data
# somevar <- system.file("external/somevar.tif", package = "flexsdm")
# somevar <- terra::rast(somevar) 
# names(somevar) <- c("aet", "cwd", "tmx", "tmn")

# 1b. import environmental data
somevar <- list.files(paste0(path,"Raster/USA"), recursive = TRUE, full.names = TRUE)
# Retain 30m resolution only
somevar <- somevar[grep("30m", somevar)]
# Keep "bioclimatic", "SRTMGL1", "landsat", "modis", "nlcd", "polaris", "tiger", "meta_pop"
numvar <- somevar[grep("bioclimatic|SRTMGL1|landsat|modis|polaris|tiger|meta_pop", somevar)]
numvar <- terra::rast(numvar)
r <- numvar[[1]]

# Define the bounding box
bbox <- ext(r)

# Create the grid
grid <- st_make_grid(bbox, cellsize = c(10, 10))

# Convert to sf object
grid <- st_sf(grid)

# Pull out first polygon
poly <- grid[1, ]

r1 <- extract(r, grid[1,], raw = TRUE, na.rm = TRUE)

# Create a big.matrix
large_matrix <- big.matrix(nrow = 92124, ncol = 214317, init = 0)



# 2. import species occurence data (presence-only)
data(hespero)
hespero <- hespero %>% dplyr::select(-id)

# 3. Import study area - California ecoregions
regions <- system.file("external/regions.tif", package = "flexsdm")
regions <- terra::rast(regions)
regions <- as.polygons(regions)
sp_region <- terra::subset(regions, regions$category == "SCR") # ecoregion where *Hesperocyparis stephensonii* is found

# visualize the species occurrences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Hesperocyparis stephensonii occurrences"
)
points(hespero[, c("x", "y")], col = "black", pch = 16)
cols <- rep("gray80", 8)
cols[regions$category == "SCR"] <- "yellow"
terra::inset(
  regions,
  loc = "bottomleft",
  scale = .3,
  col = cols
)

# 4. delimit calibration area using buffer method

ca <- calib_area(
  data = hespero,
  x = "x",
  y = "y",
  method = c('buffer', width=25000),
  crs = crs(somevar)
)

# visualize the species occurrences & calibration area
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Calibration area and occurrences")
plot(ca, add=TRUE)
points(hespero[,c("x", "y")], col = "black", pch = 16)


# 5. Sample the same number of species presences w/in the calibration area
set.seed(10)
psa <- sample_pseudoabs(
  data = hespero,
  x = "x",
  y = "y",
  n = sum(hespero$pr_ab),
  method = "random",
  rlayer = somevar,
  calibarea = ca
)

# Visualize species presences and pseudo-absences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  xlim = c(289347, 353284),
  ylim = c(-598052,  -520709),
  main = "Presence = yellow, Pseudo-absence = black")
plot(ca, add=TRUE)
points(psa[,c("x", "y")], cex=0.8, pch=16, col = "black") # Pseudo-absences
points(hespero[,c("x", "y")], col = "yellow", pch = 16, cex = 1.5) # Presences

# Bind a presences and pseudo-absences
hespero_pa <- bind_rows(hespero, psa)
hespero_pa # Presence-Pseudo-absence database

set.seed(10)

# 6. Setup data partitioning - repeated K-fold method
hespero_pa2 <- part_random(
  data = hespero_pa,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 10)
)

# 7. Extract environmental data

hespero_pa3 <-
  sdm_extract(
    data = hespero_pa2,
    x = 'x',
    y = 'y',
    env_layer = somevar,
    variables = c('aet', 'cwd', 'tmx', 'tmn')
  )

# 8. Fit models - GLM, GBM, SVM with max sens/spec threshold
mglm <-
  fit_glm(
    data = hespero_pa3,
    response = 'pr_ab',
    predictors = c('aet', 'cwd', 'tmx', 'tmn'),
    partition = '.part',
    thr = 'max_sens_spec' # same as max TSS
  )

mgbm <- fit_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

msvm <-  fit_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# 8. Ensemble models

eglm  <-
  esm_glm(
    data = hespero_pa3,
    response = 'pr_ab',
    predictors = c('aet', 'cwd', 'tmx', 'tmn'),
    partition = '.part',
    thr = 'max_sens_spec'
  )
  
egbm <- esm_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

esvm <-  esm_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)
  
# 9. Predictions

eglm_pred <- sdm_predict(
  models = eglm ,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

egbm_pred <- sdm_predict(
  models = egbm ,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

esvm_pred <- sdm_predict(
  models = esvm,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

# 10. Summarize models
merge_df <- sdm_summarize(models = list(mglm, mgbm, msvm, eglm, egbm, esvm))

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



