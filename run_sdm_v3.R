#' 1. Load the required packages
library(flexsdm)
library(terra)
library(spatialEco)
library(dplyr)
library(sf)
library(data.table)
library(foreach)
library(doParallel)

#' 2. Load the data
path <- "~/Desktop/pops_pesthostuse/"
sp_region <- vect(paste0(path, "/studyext/ext25_states_mask.gpkg"))
predictors <- readRDS(paste0(path, "hostmap/all_predictors.Rdata"))

format_species_name <- function(species) {
  # Replace " " with "_" in species name
  species <- gsub(" ", "_", species)
  # ensure the first letter is in uppercase and the rest are in lowercase
  species <- tolower(species)
  # make first letter uppercase
  species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
  return(species)
}

species <- "Ailanthus altissima"
species <- format_species_name(species)

ca <- vect(paste0(path, "hostmap/", species,  "_calib_area.gpkg"))
occs_pa <- readRDS(paste0(path, "hostmap/", species, "_occs_pa2.RData"))

#' 2. Fit models - GLM, GBM, SVM with max sens/spec threshold
# Split predictors into chunks
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
  tmpFiles(orphan=TRUE,old=TRUE,remove = TRUE)
  gc()
  rm(mglm)
}

#Re-do with lc
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
  tmpFiles(orphan=TRUE,old=TRUE,remove = TRUE)
  gc()
  rm(mglm)
}


mraf <- fit_raf(
  occs_pa,
  "pr_ab",
  predictors2[[1]],
  partition = ".part",
  thr = "max_sens_spec"
  #metric = "TSS",
  #grid = tune_grid
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
