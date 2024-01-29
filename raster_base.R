#' Function to create a 1 arc second (30m) base raster for CONUS
#' or country using the OpenTopography API. This function will download
#' the data in 450000 km^2 tiles. For each tile, the function will download
#' the data, read it into R, clamp the data to 1 and convert to INT1U, and
#' write the data to disk. The function will then mosaic all tiles into
#' a single raster, crop to country or CONUS, write to disk the base raster with
#' its original CRS (epsg:4326), and reproject to 30m resolution. Finally, the 
#' function will delete the downloaded tiles.

#' @param path Path to write the data to
#' @param country The ISO-3 code for the country to create the base raster for
#' @param conus A logical value specifying whether to create the base raster
#' for CONUS when United States is specified as the country.
#' @param proj A character string specifying the projection to reproject the
#' base raster to.
#' @return A 30 meter base raster for the given country or CONUS.
#' @export

ctry_30m_baseraster <- function(path, country, conus = FALSE) {
  country <- tolower(country)
  # Create country boundary sf object
  if (country == "usa" & conus == TRUE) {
    ctry_sf <- rgeoboundaries::geoboundaries(country,adm_lvl = 1)
    ctry_sf <- ctry_sf[ctry_sf$shapeISO != "US-PR" & ctry_sf$shapeISO
                       != "US-VI" & ctry_sf$shapeISO != "US-GU" &
                         ctry_sf$shapeISO != "US-MP" & ctry_sf$shapeISO
                       != "US-AS" & ctry_sf$shapeISO != "US-AK" &
                         ctry_sf$shapeISO != "US-HI",]
    ctry_sf <- sf::st_union(ctry_sf)
    country <- "conus"
  } else {
    ctry_sf <- rgeoboundaries::geoboundaries(country, adm_lvl = 0)
  }
  # Make grid of 450,000 sq km cells
  ctry_sf_bb <- sf::st_transform(ctry_sf, crs = 3857)
  ctry_sf_bb <- sf::st_make_grid(ctry_sf_bb, cellsize = 670820)
  ctry_sf_bb <- sf::st_transform(ctry_sf_bb, crs = 4326)
  ctry_sf_bb <- lapply(seq_along(ctry_sf_bb), function(x){
    ctry_sf_bb[[x]] <- sf::st_bbox(ctry_sf_bb[[x]])
  })
  # Download tiles
  options(timeout = 60 * 5)
  dir.create(paste0(path, "SRTMGL1"), showWarnings = FALSE)
  lapply(seq_along(ctry_sf_bb), function(x) {
    url <- paste0("https://portal.opentopography.org/API/",
                  "globaldem?demtype=","SRTMGL1",
                  "&south=", ctry_sf_bb[[x]][2],
                  "&north=", ctry_sf_bb[[x]][4],
                  "&west=", ctry_sf_bb[[x]][1],
                  "&east=", ctry_sf_bb[[x]][3],
                  "&outputFormat=GTiff&API_Key=",
                  "00e4abc8a10c75653b02fa0c46aae93a")
    # Use curl to request download url
    dl_url <- curl::curl_fetch_memory(url)$url
    download.file(dl_url, destfile = paste0(path, "SRTMGL1/srtm_", x,".tif"))
    # Clamp to 1 & convert to INT1U, and write to disk when download is complete
    if (file.exists(paste0(path, "SRTMGL1/srtm_30m_", x, ".tif"))) {
      srtm <- terra::rast(paste0(path, "SRTMGL1/srtm_30m_", x, ".tif"))
      srtm <- terra::clamp(srtm, 1, 1, datatype = "INT1U")
      terra::writeRaster(srtm, paste0(path, "SRTMGL1/srtm_", x, ".tif"),
                         overwrite = TRUE, datatype = "INT1U")
    } else {
      print(paste0("Skipping tile ", x, " because its outside srtm bounds."))
    }
  })
  # Read in tiles
  srtm_files <- list.files(paste0(path, "SRTMGL1"), pattern = ".tif",
                           full.names = TRUE)
  # If there are more than 1 file, mosaic them
  if (length(srtm_files) > 1) {
    srtm <- terra::sprc(srtm_files)
    srtm <- terra::mosaic(srtm, fun = "first", wopt = list(datatype = "INT1U"))
  } else {
    srtm <- terra::rast(srtm_files)
  }
  # Mask to country
  srtm <- terra::crop(srtm, terra::vect(ctry_sf), mask = TRUE, filename =
                        paste0(path, "basefiles/rasters/base_1s_", country,
                               ".tif"), overwrite = TRUE, wopt =
                        list(gdal = "COMPRESS=ZSTD", datatype = "INT1U"))
  srtm <- terra::project(srtm, y = proj, threads = TRUE, method = "near",
                         filename = paste0(path, "basefiles/rasters/base_30m_",
                                           country, ".tif"), overwrite = TRUE,
                         wopt = list(gdal = "COMPRESS=ZSTD",
                                     datatype = "INT1U"), res = 30)
  unlink(paste0(path, "SRTMGL1"), recursive = TRUE)
  return(srtm)
}

#' This function uses elevation_global() from the geodata package
#' to download a 30s global elevation raster. Then, the function clamps the
#' raster to 1 and converts to INT1U. The function then writes the raster
#' to disk. Next, the function reprojects the raster to EPSG:4087 and writes
#' the raster to disk.
#' @param path A character string specifying the path to download the
#' base raster to.

world_1km_baseraster <- function(path) {
  # Download 30s global elevation raster
  srtm <- geodata::elevation_global(res = 0.5, mask=TRUE, path = tempdir())
  # Clamp to 1 & convert to INT1U
  srtm <- terra::clamp(srtm, 1, 1, datatype = "INT1U")
  dir.create(paste0(path, "basefiles/rasters"), showWarnings = FALSE)
  # Write to disk
  srtm <- terra::writeRaster(srtm, paste0(path, "basefiles/rasters/",
                                          "base_30s_world.tif"),
                             overwrite = TRUE, datatype = "INT1U",
                             gdal = "COMPRESS=ZSTD")
  # Reproject to EPSG:4087
  srtm <- terra::project(srtm, y = "EPSG:4087", threads = TRUE, method = "near",
                         filename = paste0(path, "basefiles/rasters/",
                                           "base_1000m_world.tif"), res = 1000,
                         overwrite = TRUE, wopt = list(gdal = "COMPRESS=ZSTD",
                                                       datatype = "INT1U"))
  return(srtm)
}

#' @title Check for and create base raster for a given domain and resolution
#' @description This function checks if the finest resolution base raster for
#' a given domain exists. If it does not, the function creates the base raster
#' using the appropriate function. If the base raster does exist, the function
#' checks if the resolution of the base raster matches the resolution provided
#' by the user. If the resolutions do not match, the function creates the
#' target resolution base raster. If the resolutions do match, the function
#' reads the existing base raster.
#' @param country The ISO-3 code for the country to create the base raster for;
#' if left NULL, it is assumed that the analysis is global.
#' @param conus A logical value specifying whether to create the base raster
#' for CONUS when United States is specified as the country.
#' @param res A numeric value specifying the resolution of the base raster to
#' create.
#' @param proj A character string specifying the projection to reproject the
#' base raster to for country analyses.
#' @param path A character string specifying the path to write the base raster
#' to.
#' @return A base raster for the given domain and resolution.
#' @export

get_baseraster <- function(country, conus = FALSE, res = 1000, path = NULL) {
  if (is.null(country)) {
    country <- "world"
    basefiles0 <- list.files(paste0(path, "basefiles/rasters"), pattern =
                               (country), full.names = TRUE)
  } else if (country == "USA" & conus == TRUE) {
    country <- country
    basefiles0 <- list.files(paste0(path, "basefiles/rasters"), pattern =
                               ("conus"), full.names = TRUE)
  } else {
    country <- country
    basefiles0 <- list.files(paste0(path, "basefiles/rasters"), pattern =
                               (country), full.names = TRUE)
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