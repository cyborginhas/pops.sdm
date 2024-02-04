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

#'@title Functions to obtain predictors for the world and CONUS
#' Function to obtain the biovariable raster layers available for the world
#' (worldclim 30s resolution) or CONUS (USGS biocomposite ~15s resolution)
#' @param conus A logical value specifying whether to obtain the biovars for
#' CONUS when TRUE or the world when FALSE.
#' @param path A character string specifying the path to write the biovars

get_biovars <- function(conus = FALSE, path = NULL) {
  # If conus is F get biovars for the world, else get biovars for the US
  if (conus == FALSE) {
    biovars <- geodata::worldclim_global(var = "bio", res = 0.5,
                                         path = tempdir())
    dir.create(paste0(path, "Original/wc2.1_30s_bio"), showWarnings = FALSE)
    fn <- paste0(path, "Original/wc2.1_30s_bio/wc2.1_30s_bio.tif")
    terra::writeRaster(biovars, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                       filename = fn, datatype = "FLT4S")
    create_metadata(file_path = fn, file = biovars, predictor = "biovars",
                    extent = "global", path = path)
  } else {
    sbtools::item_file_download(sb_id = "4ff32906e4b0e183ef5a2f16", dest_dir =
                                  tempdir(), names =
                                  "BioClimComposite_1971_2000_400m.tif",
                                destinations =  paste0(tempdir(), "BioClimComposite_1971_2000_400m.tif"), # nolint: line_length_linter.
                                overwrite_file = TRUE)
    biovars <- terra::rast(paste0(tempdir(),
                                  "BioClimComposite_1971_2000_400m.tif"))
    names(biovars) <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
    dir.create(paste0(path, "Original/BioClimComposite"), showWarnings = FALSE)
    fn <- paste0(path, "Original/BioClimComposite/BioClimComposite_1971_2000_400m.tif")
    terra::writeRaster(biovars, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                       filename = fn, datatype = "FLT4S")
    create_metadata(file_path = fn, file = biovars, predictor = "biovars",
                    extent = "conus", path = path)
  }
  return(biovars)
}

#' Function to obtain elevation data for CONUS using the OpenTopography API
#' @param key A character string specifying the OpenTopography API key
#' @param path A character string specifying the path to write the biovars
#' @param loc An sf object specifying the extent of the domain to obtain
#' elevation data for

get_opentopo <- function(key, path = NULL, loc = NULL) {
  # Make grid of 450,000 sq km cells
  loc_bb <- sf::st_transform(loc, crs = 3857)
  loc_bb <- sf::st_make_grid(loc_bb, cellsize = 670820)
  loc_bb <- sf::st_transform(loc_bb, crs = 4326)
  loc_bb <- lapply(seq_along(loc_bb), function(x) {
    loc_bb[[x]] <- sf::st_bbox(loc_bb[[x]])
  })
  # Download tiles
  dir.create(paste0(path, "Original/SRTMGL1/tiles"), showWarnings = FALSE,
             recursive = TRUE)
  lapply(seq_along(loc_bb), function(x) {
    url <- paste0("https://portal.opentopography.org/API/",
                  "globaldem?demtype=", "SRTMGL1", "&south=", loc_bb[[x]][2],
                  "&north=", loc_bb[[x]][4], "&west=", loc_bb[[x]][1], "&east=",
                  loc_bb[[x]][3],"&outputFormat=GTiff&API_Key=", key)
    # Use curl to request download url
    dl_url <- curl::curl_fetch_memory(url)$url
    if (.Platform$OS.type == "windows") {
      download.file(dl_url, destfile = paste0(path, "Original/SRTMGL1/tiles/srtm_",
                                              x, ".tif"), mode = "wb")
    } else {
      download.file(dl_url, destfile = paste0(path, "Original/SRTMGL1/tiles/srtm_", x,
                                              ".tif"))
    }
  })
  # Read in tiles
  srtm_files <- list.files(paste0(path, "Original/SRTMGL1/tiles"), pattern = ".tif",
                           full.names = TRUE)
  srtm_files <- terra::sprc(srtm_files)
  srtm <- terra::mosaic(srtm_files, wopt = list(datatype = "INT2S"))
  names(srtm) <- "elevation"
  srtm <- terra::crop(srtm, terra::vect(loc), mask = TRUE, filename =
                        paste0(path, "Original/SRTMGL1/SRTMGL1_1s_l48.tif"),
                      overwrite = TRUE, wopt = list(datatype = "INT2S",
                                                    gdal = "COMPRESS=ZSTD"))
  # Delete tiles folder
  unlink(paste0(path, "Original/SRTMGL1/tiles"), recursive = TRUE)
  return(srtm)
}

#' Function to obtain the boundary of the contiguous United States (USA)
#' @return An sf object representing the boundary of CONUS

get_l48_boundary <- function() {
  locs <- rgeoboundaries::geoboundaries(country = "USA", adm_lvl = 1, quiet = TRUE,
                                        overwrite = TRUE)
  locs <- locs[locs$shapeISO != "US-PR" & locs$shapeISO != "US-VI" &
                 locs$shapeISO != "US-GU" & locs$shapeISO != "US-MP" &
                 locs$shapeISO != "US-AS" & locs$shapeISO != "US-AK" &
                 locs$shapeISO != "US-HI", ]
  locs <- sf::st_union(locs)
  return(locs)
}

#' Function to obtain an elevation raster layer for the world (worldclim 30s)
#' or CONUS (NASA SRTM 1 arc-second)
#' @param conus A logical value specifying whether to obtain the elevation for
#' CONUS when TRUE or the world when FALSE.
#' @param key when conus is TRUE, a character string specifying the
#' OpenTopography API key is needed.
#' @param path A character string specifying the path to write the biovars

get_elevation <- function(conus = FALSE, path, key) {
  if (conus == FALSE) {
    elevation <- geodata::elevation_global(res = 0.5, path = tempdir())
    dir.create(paste0(path, "Original/wc2.1_30s_elev"), showWarnings = FALSE)
    fn <- paste0(path, "Original/wc2.1_30s_elev/wc2.1_30s_elev.tif")
    terra::writeRaster(elevation, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                       filename = fn, datatype = "INT2S")
    create_metadata(file_path = fn, file = elevation, predictor = "elevation",
                    extent = "global", path = path)
  } else {
    conus <- get_l48_boundary
    elevation <- get_opentopo(key = key, path = path, loc = conus)
    fn <- paste0(path, "Original/SRTMGL1/SRTMGL1_1s_l48.tif")
    create_metadata(file_path = fn, file = elevation, predictor = "elevation",
                    extent = "conus", path = path)
  }
  return(elevation)
}

#' Function to obtain monthly tavg normals for the world or CONUS
#' (worldclim 30s resolution). Monthly tavg normals are used to calculate
#' growing degree days rasters.
#' @param path A character string specifying the path to write the tavg normals

get_tavg <- function(path) {
  tavg <- geodata::worldclim_global(var = "tavg", res = 0.5, path = tempdir())
  dir.create(paste0(path, "Original/wc2.1_30s_tavg"), showWarnings = FALSE)
  fn <- paste0(path, "Original/wc2.1_30s_tavg/wc2.1_30s_tavg.tif")
  terra::writeRaster(tavg, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                     filename = fn, datatype = "FLT4S")
  create_metadata(file_path = fn, file = tavg, predictor = "gdd",
                  extent = "global", path = path)
  return(tavg)
}

#' Function to obtain monthly precip normals for the world or CONUS
#' (worldclim 30s resolution). Monthly precip normals are used to calculate
#' precipitation timing rasters.
#' @param path A character string specifying the path to write the precip normals

get_precip <- function(path) {
  precip <- geodata::worldclim_global(var = "prec", res = 0.5, path = tempdir())
  dir.create(paste0(path, "Original/wc2.1_30s_prec"), showWarnings = FALSE)
  fn <- paste0(path, "Original/wc2.1_30s_prec/wc2.1_30s_prec.tif")
  terra::writeRaster(precip, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                     filename = fn, datatype = "INT2S")
  create_metadata(file_path = fn, file = precip, predictor = "precip",
                  extent = "global", path = path)
return(precip)
}

#' Function to obtain landcover data for the world (30s ESA World Cover) and
#' CONUS (30m NLCD)
#' @param conus A logical value specifying whether to obtain the landcover for
#' CONUS when TRUE or the world when FALSE.
#' @param path A character string specifying the path to write the landcover

get_landcover <- function(conus = FALSE, path) {
    if (conus == FALSE) {
        types <- c("built", "cropland", "grassland", "shrubs", "trees", "wetland")
        landcover <- lapply(seq_along(types), function(x) {
            geodata::landcover(var = types[[x]], path = tempdir())
        })
        dir.create(paste0(path, "Original/WorldCover"), showWarnings = FALSE)
        fn <- paste0(path, "Original/WorldCover/", types, "_30s.tif")
        lapply(seq_along(fn), function(x) {
            terra::writeRaster(landcover[[x]], overwrite = TRUE, 
            gdal = "COMPRESS=ZSTD", filename = fn[x], datatype = "INT1U")
        })
        lapply(seq_along(landcover), function(x) {
            create_metadata(file_path = fn[x], file = landcover[[x]],
            predictor = "landcover", extent = "global", path = path)
        })
    } else {
        options(timeout = 60 * 5)
        url <- "https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2021_land_cover_l48_20230630.zip" # nolint: line_length_linter.
        # Download zip file
        download.file(url, destfile = paste0(tempdir(), "/", basename(url)), mode = "wb")
        unzip(paste0(tempdir(), "/", basename(url)), exdir = paste0(tempdir(), "/nlcd"))
        # Read in raster
        landcover <- terra::rast(paste0(tempdir(),"/nlcd/nlcd_2021_land_cover_l48_20230630.img"))
        dir.create(paste0(path, "Original/nlcd"), showWarnings = FALSE, recursive = TRUE)
        fn <- paste0(path, "Original/nlcd/nlcd_2021_land_cover_l48_20230630.tif")
        terra::writeRaster(landcover, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = fn, datatype = "INT1U")
        create_metadata(file_path = fn, file = landcover, predictor = "landcover",
        extent = "conus", path = path)
    }
    return(landcover)
}

#' Function to download human population density data for the world and CONUS
#' (30s GHS_POP_E2020_GLOBE)
#' @param path A character string specifying the path to write the population data

get_pop <- function(path) {
    pop <- geodata::population(year = 2020, res = 0.5, path = tempdir())
    dir.create(paste0(path, "Original/gpw_v4_population_density"), showWarnings = FALSE)
    fn <- paste0(path, "Original/gpw_v4_population_density/gpw_v4_population_density_rev11_2020_30s.tif") # nolint: line_length_linter.
    terra::writeRaster(pop, overwrite = TRUE, gdal = "COMPRESS=ZSTD", filename = fn, datatype = "FLT4S")
    create_metadata(file_path = fn, file = pop, predictor = "pop", extent = "global", path = path)
    return(pop)
}

#' Function to download roads data for the world and CONUS
#' (10m rnaturalearth)
#' @param path A character string specifying the path to write the roads

get_roads <- function(path) {
  roads <- rnaturalearth::ne_download(type = "roads", scale = 10)
  dir.create(paste0(path, "Original/ne_roads"), showWarnings = FALSE)
  fn <- paste0(path, "Original/ne_roads/ne_10m_roads.gpkg")
  sf::st_write(roads, dsn = fn)
  create_metadata(file_path = fn, file = terra::vect(roads), predictor = "roads", 
  extent = "global", path = path)
  return(roads)
}

#' Function to download railroads data for the world and CONUS
#' (10m rnaturalearth)
#' @param path A character string specifying the path to write the railroads

get_railroads <- function(path) {
  railroads <- rnaturalearth::ne_download(type = "railroads", scale = 10)
  dir.create(paste0(path, "Original/ne_railroads"), showWarnings = FALSE)
  fn <- paste0(path, "Original/ne_railroads/ne_10m_railroads.gpkg")
  sf::st_write(railroads, dsn = fn)
  create_metadata(file_path = fn, file = terra::vect(railroads), predictor = "railroads",
  extent = "global", path = path)
  return(railroads)
}

#' Function to download soil properties (pH, H2O-33 and H2O-1500) for the world
#' and CONUS (250m OpenLandMap)
#' @param path A character string specifying the path to write the soil properties

get_soil_vars <- function(path){
  dl_urls <- c("https://doi.org/10.5281/zenodo.2525664", "https://doi.org/10.5281/zenodo.2784001")
  names(dl_urls) <- c("ph", "watercontent")
  
  lapply(seq_along(dl_urls), function(x) {
    dir.create(paste0(tempdir(), "/OpenLandMap/soils/", names(dl_urls[x])), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(path, "Original/OpenLandMap/soils/", names(dl_urls[x])), showWarnings = FALSE, recursive = TRUE)
    zen4R::download_zenodo(dl_urls[x], path = paste0(tempdir(), "/OpenLandMap/soils/", names(dl_urls)[x]),
                           timeout = 60 * 30)
                           })
  # List files in tempdir
  files <- list.files(paste0(tempdir(), "/OpenLandMap/soils/"), full.names = TRUE, recursive = TRUE)
  
  # Keep only .tif files
  files <- files[grepl(".tif$", files)]
  files <- files[!grepl("_md_", files)]
  files <- sort(files)
  
  # Create export file names
  fn <- c(paste0(path, "Original/OpenLandMap/soils/ph/", basename(files[grepl("sol_ph", files)])),
  paste0(path, "Original/OpenLandMap/soils/watercontent/", basename(files[grepl("sol_watercontent", files)])))
  fn <- sort(fn)
  
  # For each file, check if it is pH or water content, and write to appropriate folder
  for (x in seq_along(files)) {
    if (grepl("sol_ph", files[x])) {
      terra::writeRaster(terra::rast(files[x]), overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                         filename = fn[x], datatype = "INT1U")
    } else if (grepl("sol_watercontent", files[x])) {
      terra::writeRaster(terra::rast(files[x]), overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                         filename = fn[x], datatype = "INT1U")
    } else {
      print(paste0("Skipped file ", files[x], " because it is not pH or water content"))
    }
    print(paste0("Finished writing ", files[x]))
  }
  lapply(seq_along(fn), function(x) {
    create_metadata(file_path = fn[x], file = terra::rast(fn[x]), predictor = "soil",
    extent = "global", path = path)
  })
}
```
