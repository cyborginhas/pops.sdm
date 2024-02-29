# Load packages
```{r load-packages, message = FALSE, warning = FALSE}
require(geodata)
require(terra)
require(sbtools)
require(rgeoboundaries)
require(zen4R)
```
```{r get-envi, message = FALSE, warning = FALSE}

#' @title Append to metadata file
#' @description This function checks for a pops.sdm.metadata file in the
#' "Original" data folder. If the file exists, the function appends the
#' metadata to the file. If the file does not exist, the function creates the
#' file and adds the metadata to the file.
#' @param file_path The path of the predictor variable file exported to the
#' "Original" data folder.
#' @param file The terra raster or vector object exported to the file_path.
#' Several metadata fields are extracted from the file and added to the
#' metadata file. Including the number of layers, resolution, extent,
#' projection.
#' @param predictor A character string specifying the broad category of the
#' predictor variable. For example, "biovars" or "elevation".
#' @param extent A character string specifying the extent of the predictor
#' variable as "global" or "conus".
#' @param path A character string specifying the path to write the metadata file.

create_metadata <- function(file_path, file, predictor, extent, path) {
  file_path <- gsub(path,"", file_path)
  if (file.exists(paste0(path, "Original/pops_sdm_orig_data_meta.csv"))) {
    metadata <- read.csv(paste0(path, "Original/pops_sdm_orig_data_meta.csv"), header = TRUE)
    if (class(file) == "SpatRaster") {
      metadata <- rbind(metadata, data.frame(file_path = file_path,
                                             nlayers = terra::nlyr(file),
                                             res = paste(terra::res(file), collapse = "x"),
                                             extent = extent,
                                             projection = terra::crs(file),
                                             names = paste(names(file), collapse = ","),
                                             predictor = predictor,
                                             download_date = as.character(Sys.Date())))
    } else if (class(file) == "SpatVector") {
      metadata <- rbind(metadata, data.frame(file_path = file_path,
                                             nlayers = NA,
                                             res = NA,
                                             extent = extent,
                                             projection = terra::crs(file),
                                             names = NA,
                                             predictor = predictor,
                                             download_date = as.character(Sys.Date())))
    }
    write.csv(metadata, paste0(path, "Original/pops_sdm_orig_data_meta.csv"), row.names = FALSE)
  } else {
    metadata <- data.frame(file_path = file_path,
                           nlayers = terra::nlyr(file),
                           res = if (class(file)[1]=="SpatRaster") paste0(terra::res(file), 
                           collapse = "x") else NA,
                           extent = extent,
                           projection = terra::crs(file),
                           names = if (class(file)[1]=="SpatRaster") paste(names(file), collapse = ",") else NA,
                           predictor = predictor,
                           download_date = as.character(Sys.Date()))
    write.csv(metadata, paste0(path, "Original/pops_sdm_orig_data_meta.csv"), row.names = FALSE)
  }
}
