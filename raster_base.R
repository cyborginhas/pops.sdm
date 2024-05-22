#' Function to create a 1 arc second (30m) base raster for USA;
#' and aggregate to common resolutions: 30m, 100m, 250m, 500m, 1000m
#' @param path Path to write the data to
#' @return A print message indicating that the base rasters have been created
#' @export

base_conus <- function(path) {
  # Setup base raster filenames for check
  proj <- c(rep("epsg4326", 5))
  common_res_names <- c(30, 100, 250, 500, 1000)
  fn <- paste0(
    path, "Raster/USA/pops_sdm_base/base_", common_res_names, "m_conus_",
    proj, ".tif"
  )
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 30m base raster for USA from 1s L48 boundary and 1s elevation
    l48 <- terra::vect(get_l48_boundary())
    elev <- get_topo_conus(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    l48 <- terra::project(l48, y = terra::crs(elev))
    dir.create(paste0(path, "Raster/USA/pops_sdm_base/"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    elev <- terra::crop(elev, l48,
      mask = TRUE, filename =
        paste0(path, "Raster/USA/pops_sdm_base/base_30m_conus_epsg4326.tif"), # nolint line length
      overwrite = TRUE, wopt = list(
        gdal = "COMPRESS=ZSTD",
        datatype = "INT1U"
      )
    )

    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_names <- common_res_names[-1]
    common_res_vals <- c(
      base_res * 3, base_res * 7.5, base_res * 15,
      base_res * 30
    )
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(
        path, "Raster/USA/pops_sdm_base/base_",
        common_res_names[x], "m_conus_epsg4326.tif"
      )
      terra::project(elev,
        y = "EPSG:4326", threads = TRUE, method = "near",
        filename = base_filename, res = common_res_vals[x],
        overwrite = TRUE, wopt = list(
          gdal = "COMPRESS=ZSTD",
          datatype = "INT1U"
        )
      )
    })
    msg <- paste0(
      "Base rasters created for USA at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  } else {
    msg <- paste0(
      "Base rasters already exist for USA at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  }
  base <- lapply(fn, terra::rast)
  return(base)
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
  fn <- paste0(
    path, "Raster/Global/pops_sdm_base/base_", common_res_names, "m_conus_",
    proj, ".tif"
  )
  # Check if all fn do not exist, create base rasters
  if (!all(file.exists(fn))) {
    # Create 1000m base raster for globe from wc2.1_30s elevation
    elev <- get_topo_global(path)[[1]]
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")
    dir.create(paste0(path, "Raster/Global/pops_sdm_base/"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    # Base rasters at common resolutions: 30m, 100m, 250m, 500m, 1000m EPSG:4326
    base_res <- terra::res(elev)[1]
    common_res_vals <- c(base_res, base_res * 2.5, base_res * 5)
    common_res_names <- common_res_names
    lapply(seq_along(common_res_vals), function(x) {
      base_filename <- paste0(
        path, "Raster/Global/pops_sdm_base/base_",
        common_res_names[x], "m_global_epsg4326.tif"
      )
      terra::project(elev,
        y = "EPSG:4326", threads = TRUE, method = "near",
        filename = base_filename, res = common_res_vals[x],
        overwrite = TRUE, wopt = list(
          gdal = "COMPRESS=ZSTD",
          datatype = "INT1U"
        )
      )
    })
    msg <- paste0(
      "Base rasters created for the Globe at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  } else {
    msg <- paste0(
      "Base rasters already exist for the Globe at resolutions: ",
      paste0(common_res_names, collapse = ", "), "m"
    )
  }
  base <- lapply(fn, terra::rast)
  return(base)
}