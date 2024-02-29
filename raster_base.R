#' Function to create a 1 arc second (30m) base raster for CONUS;
#' and aggregate to common resolutions: 30m, 100m, 250m, 500m, 1000m
#' @param path Path to write the data to
#' @return A print message indicating that the base rasters have been created
#' @export
base_conus <- function(path) {
  # Setup base raster filenames for check
  proj <- c(rep("epsg4326", 5))
  common_res_names <- c(30, 100, 250, 500, 1000)
  fn <- paste0(path, "Raster/USA/pops_sdm_base/base_", common_res_names, "m_conus_",
               proj, ".tif")
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 30m base raster for CONUS from 1s L48 boundary and 1s elevation
    l48 <- terra::vect(get_l48_boundary())
    elev <- get_topo_conus(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    l48 <- terra::project(l48, y = terra::crs(elev))
    dir.create(paste0(path, "Raster/USA/pops_sdm_base/"), showWarnings = FALSE,
               recursive = TRUE)
    elev <- terra::crop(elev, l48, mask = TRUE, filename =
                        paste0(path, "Raster/USA/pops_sdm_base/base_30m_conus_epsg4326.tif"), #nolint line length
                        overwrite = TRUE, wopt = list(gdal = "COMPRESS=ZSTD",
                                                      datatype = "INT1U"))

    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_names <- common_res_names[-1]
    common_res_vals <- c(base_res * 3, base_res * 7.5, base_res * 15,
                    base_res * 30)
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(path, "Raster/USA/pops_sdm_base/base_",
                              common_res_names[x], "m_conus_epsg4326.tif")
      terra::project(elev, y = "EPSG:4326", threads = TRUE, method = "near",
                     filename = base_filename, res = common_res_vals[x],
                     overwrite = TRUE, wopt = list(gdal = "COMPRESS=ZSTD",
                                                  datatype = "INT1U"))
    })
    msg <- paste0("Base rasters created for CONUS at resolutions: ",
                  paste0(common_res, collapse = ", "), "m")
  } else {
    msg <- paste0("Base rasters already exist for CONUS at resolutions: ",
                  paste0(common_res, collapse = ", "), "m")
  }
  return(msg)
}

#' Function to create a 30 arc second (1km) base raster for the world;
#' and aggregate to common resolutions: 2.5km, 5km, 10km.
#' @param path Path to write the data to
#' @return A print message indicating that the base rasters have been created.
#' @export

base_global <- function(path) {
  # Setup base raster filenames for check
  proj <- c(rep("epsg4326", 5))
  common_res_names <- c(1000, 2500, 5000)
  fn <- paste0(path, "Raster/Global/pops_sdm_base/base_", common_res_names, "m_conus_",
               proj, ".tif")
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 1000m base raster for globe from wc2.1_30s elevation
    elev <- get_topo_global(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    dir.create(paste0(path, "Raster/Global/pops_sdm_base/"), showWarnings = FALSE,
               recursive = TRUE)
    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_vals <- c(base_res, base_res * 2.5, base_res * 5)
    common_res_names <- common_res_names
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(path, "Raster/Global/pops_sdm_base/base_",
                              common_res_names[x], "m_global_epsg4326.tif")
      terra::project(elev, y = "EPSG:4326", threads = TRUE, method = "near",
                     filename = base_filename, res = common_res_vals[x],
                     overwrite = TRUE, wopt = list(gdal = "COMPRESS=ZSTD",
                                                  datatype = "INT1U"))
    })
    msg <- paste0("Base rasters created for the Globe at resolutions: ",
                  paste0(common_res, collapse = ", "), "m")
  } else {
    msg <- paste0("Base rasters already exist for the Globe at resolutions: ",
                  paste0(common_res, collapse = ", "), "m")
  }
  return(msg)
}

#' @description Function to check if the resolution provided by the user
#' has already been created for the domain. If the resolution has not been
#' created, the function will create it, and return it. If the resolution
#' has been created, the function will return the existing raster. Resolutions
#' requested below 30m for CONUS and below 1000m for the world will be
#' returned as the finest resolution base raster for the domain.
#' @param domain A character string specifying the domain for the analysis.
#' Either "global" or "conus".
#' @param res A numeric value specifying the resolution of the base raster to
#' create.
#' @param path A character string specifying the path to write the base raster
#' to.
#' @return A base raster for the given domain and resolution.
#' @export

get_baseraster <- function(country, conus = FALSE, res = 1000, path = NULL) {
  domain <- tolower(domain)
  if (domain == "global") {
    basefiles0 <- list.files(paste0(path, "basefiles/rasters"), pattern =
                               (domain), full.names = TRUE)
  } else if (domain == "conus") {
    domain <- domain
    basefiles0 <- list.files(paste0(path, "basefiles/rasters"), pattern =
                               (domain), full.names = TRUE)
  } else {
    stop("Domain must be 'global' or 'conus'")
  }

  # Create base raster if it does not exist
  if (length(basefiles0) == 0) {
    # If base raster does not exist, create it
    if (country == "world") {
      world_1km_baseraster(path)
      basefiles0 <- list.files(paste0(path, "basefiles/rasters"),
                               pattern = (country), full.names = TRUE)
    } else {
      ctry_30m_baseraster(path, country, conus, proj)
      if (conus == TRUE) {
        basefiles0 <- list.files(paste0(path, "basefiles/rasters"),
                                 pattern = (country), full.names = TRUE)
      } else {
        basefiles0 <- list.files(paste0(path, "basefiles/rasters"),
                                 pattern = (country), full.names = TRUE)
      }
    }
  }
  # Create raster at target resolution
  basefiles1 <- basefiles0[grepl(paste0(res, "m"), basefiles0)]
  if (length(basefiles1) == 0) {
    # If resolutions do not match, create base raster at target resolution
    r <- terra::rast(basefiles0[grep("base_", basefiles0)])
    country <- ifelse(conus == TRUE, "conus", country)
    base_raster <- terra::project(r, y = terra::crs(r), threads = TRUE,
                                  method = "near", filename =
                                    paste0(path, "basefiles/rasters/coarsened_",
                                           res, "m_", country, ".tif"),
                                  res = res, overwrite = TRUE, wopt =
                                    list(gdal = "COMPRESS=ZSTD",
                                    datatype = "INT1U"))
  } else {
    base_raster <- terra::rast(basefiles1)
  }
  return(base_raster)
}
