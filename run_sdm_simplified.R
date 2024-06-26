#' Run the species distribution model: currently runs for Maxent only. Will add in other models later.

#' 1. Load the required packages & set path
library(flexsdm)
library(terra)
library(data.table)
library(parallel)
library(spatialEco)


#' Prep data for maxent
tg_pts_wdata <- "C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_target_group_bg_pts_filtgeo_210m2_wpreds.csv"
tg_pts_wdata <- fread(tg_pts_wdata)
occs <- splitpts4sdm(tg_pts_wdata)

max_t1 <- readRDS("C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_maxent.rds")
tiles <- list.files("C:/Users/blaginh/Desktop/maximize/tiles/tiles/", pattern = ".tif", full.names = TRUE)
tiles <- lapply(tiles, terra::rast)

# a list of models
for (i in seq_along(tiles)) {
  tryCatch({
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
    writeRaster(ind_p$max, paste0("C:/Users/blaginh/Desktop/maximize/output/Ailanthus_altissima_prediction_",fn))
    t <- Sys.time() - s
    print(t)
    tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  }, error = function(e) {
    print(paste("Error occurred for tile", i))
    print(e)
  })
}

