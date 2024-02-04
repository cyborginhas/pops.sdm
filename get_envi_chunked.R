# Load packages
require(geodata)
require(terra)
require(sbtools)
require(rnaturalearth)
require(zen4R)
require(rgee)
require(foreach)

#' Function to get biovariables for the world (worldclim 30s)
#' @param path A character string specifying the path to write the biovars
#' @return A raster with 19 biovariable layers for the world.
#' @export

get_biovars_global <- function(path) {
  fn <- paste0(path, "Original/wc2.1_30s_bio/wc2.1_30s_bio_", 1:19, ".tif")
  if (!all(file.exists(fn))) {
    biovars <- geodata::worldclim_global(
      var = "bio", res = 0.5,
      path = tempdir()
    )
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    # Write each layer to file
    lapply(1:terra::nlyr(biovars), function(x) {
      terra::writeRaster(biovars[[x]],
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = paste0(dirname(fn[1]),
          names(biovars)[x], ".tif"
        ),
        datatype = "FLT4S"
      )
    })
  } else {
    biovars <- lapply(fn, terra::rast)
  }
  return(biovars)
}

#' @description Function to get biovariables for CONUS (USGS biocomposite 15s)
#' @param path A character string specifying the path to write the biovars.
#' @return A raster with 19 biovariable layers for CONUS.
#' @export

get_biovars_conus <- function(path) {
  bio_names <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
  fn <- paste0(path, "Original/BioClimComposite/BioClimComposite_1971_2000_400m_", bio_names, ".tif") # nolint: line_length_linter.
  if (!all(file.exists(fn))) {
    sbtools::item_file_download(
      sb_id = "4ff32906e4b0e183ef5a2f16", dest_dir =
      tempdir(), names = "BioClimComposite_1971_2000_400m.tif", # nolint
      destinations = paste0(tempdir(), "BioClimComposite_1971_2000_400m.tif"),
      overwrite_file = TRUE
    )
    biovars <- terra::rast(paste0(
      tempdir(),
      "BioClimComposite_1971_2000_400m.tif" # nolint
    ))
    names(biovars) <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    # Write each layer to file
    lapply(1:terra::nlyr(biovars), function(x) {
      terra::writeRaster(biovars[[x]],
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = paste0(dirname(fn[1]), "/BioClimComposite_1971_2000_400m_",
                          names(biovars)[x], ".tif"),
        datatype = "FLT4S"
      )
    })
  } else {
    biovars <- lapply(fn, terra::rast)
  }
  return(biovars)
}

#' @description Function to obtain a vector of CONUS.
#' @return A vector of CONUS.
#' @export
get_l48_boundary <- function() {
  locs <- rnaturalearth::ne_states(country = "United States of America",
                                   returnclass = "sf")
  locs <- unique(locs["iso_3166_2"])
  locs <- locs[locs$iso_3166_2 != "US-PR" & locs$iso_3166_2 != "US-VI" &
                 locs$iso_3166_2 != "US-GU" & locs$iso_3166_2 != "US-MP" &
                 locs$iso_3166_2 != "US-AS" & locs$iso_3166_2 != "US-AK" &
                 locs$iso_3166_2 != "US-HI", ]
  locs <- sf::st_union(locs)
  plot(locs)
  return(locs)
}

#' @description Function to create terrain variables (slope, aspect,
#' & hillshade) from an elevation layer.
#' @param path A character string with the path to write the terrain variables.
#' @return A hillshade raster.
#' @export

get_hillshade <- function(path, elevation) {
  fn <- paste0(c("slope", "aspect", "hillshade"), ".tif")
  # Calculate slope & export
  slp <- terra::terrain(elevation, v = "slope", unit = "radians", neighbors = 8)
  names(slp) <- "slope"
  terra::writeRaster(slp,
    filename = paste0(path, fn[1]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  # Calculate aspect & export
  asp <- terra::terrain(elevation,
    v = "aspect", unit = "radians",
    neighbors = 4
  )
  names(asp) <- "aspect"
  terra::writeRaster(asp,
    filename = paste0(path, fn[2]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  # Calculate hillshade & export
  shd <- terra::shade(slp, asp, direction = 180)
  names(shd) <- "hillshade"
  terra::writeRaster(shd,
    filename = paste0(path, fn[3]), overwrite = TRUE, datatype = "FLT4S",
    gdal = "COMPRESS=ZSTD"
  )
  return(shd)
}

#' @description Function to obtain a global elevation raster (30s WorldClim).
#' and to calculate slope, aspect, and hillshade from the elevation layer.
#' @param path A character string specifying the path to write the elevation.
#' @return A list with elevation & hillshade rasters.
#' @export

get_topo_global <- function(path) {
  fn <- c(
    paste0(path, "Original/wc2.1_30s_elev/wc2.1_30s_elev.tif"),
    paste0(path, "Raster/Global/wc2.1_30s_elev/wc2.1_30s_hillshade.tif")
  )
  if (!all(file.exists(fn))) {
    elevation <- geodata::elevation_global(res = 0.5, path = tempdir())
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(elevation,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn, datatype = "INT2S"
    )
    shd <- get_hillshade(path = dirname(fn[2]),
                         elevation)
    topo <- list(elevation, shd)
  } else {
    topo <- list(terra::rast(fn[1]), terra::rast(fn[2]))
  }
  return(topo)
}

#' @description Function to get a elevation layer for CONUS (NASA SRTM 1s),
#' and to calculate slope, aspect, and hillshade from the elevation layer.
#' @param key A character string specifying the OpenTopography API key.
#' @param path A character string specifying the path to write the elevation.
#' @export

get_topo_conus <- function(path, key) {
  fn <- c(
    paste0(path, "Original/SRTMGL1/SRTMGL1_1s_l48.tif"),
    paste0(path, "Raster/USA/SRTMGL1/SRTMGL1_1s_l48_hillshade.tif")
  )
  if (!all(file.exists(fn))) {
    # Get bounding box of CONUS
    loc <- get_l48_boundary()
    # Make grid of 450,000 sq km cells
    loc_bb <- sf::st_transform(loc, crs = 3857)
    loc_bb <- sf::st_make_grid(loc_bb, cellsize = 670820)
    loc_bb <- sf::st_transform(loc_bb, crs = 4326)
    loc_bb <- lapply(seq_along(loc_bb), function(x) {
      loc_bb[[x]] <- sf::st_bbox(loc_bb[[x]])
    })
    # Download tiles
    dir.create(paste0(dirname(fn[1]), "/tiles"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    lapply(seq_along(loc_bb), function(x) {
      url <- paste0(
        "https://portal.opentopography.org/API/",
        "globaldem?demtype=", "SRTMGL1", "&south=", loc_bb[[x]][2],
        "&north=", loc_bb[[x]][4], "&west=", loc_bb[[x]][1], "&east=",
        loc_bb[[x]][3], "&outputFormat=GTiff&API_Key=", key
      )
      # Use curl to request download url
      dl_url <- curl::curl_fetch_memory(url)$url
      if (.Platform$OS.type == "windows") {
        download.file(dl_url, destfile = paste0(
          path, "Original/SRTMGL1/tiles/srtm_",
          x, ".tif"
        ), mode = "wb")
      } else {
        download.file(dl_url, destfile = paste0(
          path, "Original/SRTMGL1/tiles/srtm_", x,
          ".tif"
        ))
      }
    })
    # Read in tiles
    srtm_files <- list.files(paste0(path, "Original/SRTMGL1/tiles"),
      pattern = ".tif",
      full.names = TRUE
    )
    srtm_files <- terra::sprc(srtm_files)
    srtm <- terra::mosaic(srtm_files, wopt = list(datatype = "INT2S"))
    names(srtm) <- "elevation"
    srtm <- terra::crop(srtm, terra::vect(loc),
      mask = TRUE, filename = fn[1],
      overwrite = TRUE, wopt = list(
        datatype = "INT2S",
        gdal = "COMPRESS=ZSTD"
      )
    )
    # Delete tiles folder
    unlink(paste0(path, "Original/SRTMGL1/tiles"), recursive = TRUE)
    # Calculate hillshade
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    shd <- get_hillshade(path = paste0(dirname(fn[2]), "/SRTMGL1_1s_l48_"),
                         srtm)
    topo <- list(srtm, shd)
  } else {
    topo <- list(terra::rast(fn[1]), terra::rast(fn[2]))
  }
  return(topo)
}

#' @description Function to get landcover for the world (30s ESA World Cover)
#' @param path A character string specifying the path to write the landcover
#' @return A list with landcover rasters: built, cropland, grassland, shrubs,
#' trees, wetland.
#' @export

get_landcover_global <- function(path) {
  types <- c("built", "cropland", "grassland", "shrubs", "trees", "wetland")
  fn <- paste0(path, "Original/WorldCover/", types, "_30s.tif")
  if (!all(file.exists(fn))) {
    landcover <- lapply(seq_along(types), function(x) {
      geodata::landcover(var = types[[x]], path = tempdir())
    })
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    lapply(seq_along(fn), function(x) {
      terra::writeRaster(landcover[[x]],
        overwrite = TRUE,
        gdal = "COMPRESS=ZSTD", filename = fn[x], datatype = "INT1U"
      )
    })
  } else {
    landcover <- lapply(seq_along(fn), function(x) {
      terra::rast(fn[x])
    })
  }
  return(landcover)
}

#' @description Function to obtain landcover data for CONUS (30m NLCD)
#' @param path A character string specifying the path to write the landcover
#' @return A list with landcover rasters: built, deciduous, evergreen, trees,
#' shrubs, grass, pasture, cropland, cultivated, wetland.
#' @export

get_landcover_conus <- function(path) {
  landcover_types <- c("built", "deciduous", "evergreen", "trees", "shrubs",
                       "grassland", "pasture", "cropland", "cultivated",
                       "wetland")
  fn <- c(paste0(path, "Original/nlcd/nlcd_2021_land_cover_l48_20230630.tif"),
          paste0(path, "Raster/USA/nlcd/nlcd_2021_land_cover_l48_",
                 landcover_types,
                 ".tif"))
  if (!all(file.exists(fn))) {
    url <- "https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2021_land_cover_l48_20230630.zip" # nolint: line_length_linter.
    # Download zip file
    options(timeout = 60 * 5)
    download.file(url, destfile = paste0(tempdir(), "/", basename(url)),
                  mode = "wb")
    unzip(paste0(tempdir(), "/", basename(url)), exdir = paste0(tempdir(),
                                                                "/nlcd"))
    # Read in raster
    landcover <- terra::rast(paste0(tempdir(), "/nlcd/nlcd_2021_land_cover_l48_20230630.img")) # nolint: line_length_linter.
    dir.create(dirname(fn[1]), showWarnings = FALSE,
               recursive = TRUE)
    terra::writeRaster(landcover,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn, datatype = "INT1U"
    )
    # Define landcover types and corresponding reclassification rules
    reclass_rules <- list(
      built = rbind(c(1, 20, 0), c(21, 24, 1), c(25, 255, 0)), #21, 22, 23, 24
      decid = rbind(c(1, 40, 0), c(41, 41, 1), c(42, 42, 0),
                    c(43, 43, 1), c(44, 255, 0)), #41, 43
      everg = rbind(c(1, 41, 0), c(42, 43, 1), c(44, 255, 0)), #42, 43
      trees = rbind(c(1, 40, 0), c(41, 43, 1), c(44, 255, 0)), #41, 42, 43
      shrub = rbind(c(1, 51, 0), c(52, 52, 1), c(53, 255, 0)), #52
      grass = rbind(c(1, 70, 0), c(71, 71, 1), c(72, 255, 0)), #71
      pastr = rbind(c(1, 80, 0), c(81, 81, 1), c(82, 255, 0)), #81
      cropl = rbind(c(1, 81, 0), c(82, 82, 1), c(83, 255, 0)), #82
      culti = rbind(c(1, 80, 0), c(81, 82, 1), c(83, 255, 0)), #81, 82
      wetld = rbind(c(1, 89, 0), c(90, 90, 1), c(91, 94, 0), c(95, 95, 1),
                    c(96, 255, 0)) #90, 95
    )
    for (i in seq_along(landcover_types)) {
      reclass_rule <- reclass_rules[[i]]
      dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
      output_file <- paste0(dirname(fn[2]), "/nlcd_2021_land_cover_l48_",
                            landcover_types[i], ".tif")
      start_time <- Sys.time() # Start time
      r <- terra::classify(landcover, reclass_rule, right = NA, others = NA)
      l48 <- get_l48_boundary()
      l48 <- sf::st_transform(l48, crs = terra::crs(r))
      r <- terra::crop(r, terra::vect(l48), mask = TRUE)
      names(r) <- landcover_types[i]
      terra::writeRaster(r,
        overwrite = TRUE, gdal = "COMPRESS=ZSTD",
        filename = output_file, datatype = "INT1U"
      )
      end_time <- Sys.time() # End time
      time_taken <- end_time - start_time
      print(paste("Time taken for", landcover_types[i], ":", time_taken))
    }
  } else {
    landcover <- lapply(fn[-1], terra::rast)
  }
  return(landcover)
}

#' @description Function to obtain monthly tavg  for the world (worldclim 30s)
#' and then calculate growing degree days.
#' @param path A character string specifying the path to write the tavg normals
#' @export

get_gdd_global <- function(path) {
  fn <- c(
    paste0(path, "Original/wc2.1_30s_tavg/wc2.1_30s_tavg.tif"), # original
    paste0(path, "Raster/Global/wc2.1_30s_tavg/wc2.1_30s_gdd.tif")
  ) # gdd
  if (!all(file.exists(fn))) {
    tavg <- geodata::worldclim_global(var = "tavg", res = 0.5, path = tempdir())
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(tavg,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn[1], datatype = "FLT4S"
    )
    tbase <- 5
    gdd <- terra::tapp(x = tavg, tbase = tbase, fun = function(tavg, tbase) {
      days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      xbase <- tavg - tbase
      xbase[(xbase < 0)] <- 0
      xbase <- sum(xbase * days)
      return(xbase)
    })
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    names(gdd) <- "gdd_tbase5"
    terra::writeRaster(gdd,
      filename = fn[2], overwrite = TRUE, datatype = "FLT4S",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    gdd <- terra::rast(fn[2])
  }
  return(gdd)
}

#' @description Function to obtain monthly precip normals for the world
#' (Worldclim 30s)
#' Monthly precip normals are used to calculate precipitation timing.
#' @param path A character string specifying the path to write precip normals
#' @export

get_prectiming_global <- function(path) {
  fn <- c(paste0(path, "Original/wc2.1_30s_prec/wc2.1_30s_prec.tif"),
          paste0(path, "Raster/Global/wc2.1_30s_prec/wc2.1_30s_prec_timing.tif")) # nolint: line_length_linter.
  if (!all(file.exists(fn))) {
    precip <- geodata::worldclim_global(var = "prec", res = 0.5,
                                        path = tempdir())
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(precip,
      overwrite = TRUE, gdal = "COMPRESS=ZSTD",
      filename = fn[1], datatype = "INT2S"
    )
    precip <- terra::app(x = precip, fun = function(x) {
      return(sum(x[[12]], x[[1]], x[[2]]) - sum(x[[6]], x[[7]], x[[8]]))
    })
    names(precip) <- "PrecipTiming"
    precip <- terra::writeRaster(precip,
      filename = paste0(dirname(fn[2]), "/wc2.1_30s_prec_timing.tif"),
      overwrite = TRUE, datatype = "INT2S", gdal = "COMPRESS=ZSTD"
    )
  } else {
    precip <- terra::rast(fn[2])
  }
  return(precip)
}

#' @description Function to download human pop. density for the world
#' (30s GHS_POP_E2020_GLOBE)
#' @param path A character string specifying the path to write the pop data
#' @return A raster with human population density for the world
#' @export

get_pop_global <- function(path) {
  fn <- paste0(path, "Original/gpw_v4_population_density/gpw_v4_population_density_rev11_2020_30s.tif") # nolint: line_length_linter.
  if (!file.exists(fn)) {
    pop <- geodata::population(year = 2020, res = 0.5, path = tempdir())
    dir.create(dirname(fn), showWarnings = FALSE)
    terra::writeRaster(pop, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                       filename = fn, datatype = "FLT4S")
  } else {
    pop <- terra::rast(fn)
  }
  return(pop)
}

#' @description Function to calc distance to roads/railroads for the world (30s)
#' @param lines A sf object of roads/railroads.
#' @param domain A character string specifying the domain (world or conus)
#' @param type A character string specifying the type of network (roads/rails)
#' @return A raster with distance to roads or railroads for the world
#' @import rgee
#' @export

dist_networks <- function(user, lines, domain, type) {
  ee_Initialize(user, drive = TRUE)
  ee_Authenticate()
  # Remove all columns except geometry
  lines <- lines[, 1]
  # Split into 1200 chunks if type == rails, else split into 3000 chunks
  if (type == "rails") {
    lines <- split(lines, ceiling(1:nrow(lines) / 1200)) #nolint
  } else {
    lines <- split(lines, ceiling(1:nrow(lines) / 3000)) #nolint
  }
  # Create rgee folder in assets
  asset <- ee_get_assethome()
  ee_manage_create(paste0(asset, "/rgee/"))
  # Merge the assets
  assetids <- paste0(asset, "/rgee/", type, "_", 1:length(lines)) #nolint
  # Convert sf to ee$FeatureCollection
  lapply(seq_along(lines), function(x) {
    sf_as_ee(
      x = lines[[x]],
      via = "getInfo_to_asset",
      assetId = assetids[x],
      overwrite = TRUE,
      monitoring = TRUE,
      proj = "EPSG:4326"
    )
  })

  assetids <- ee_manage_assetlist(paste0(ee_get_assethome(), "/rgee/"))$ID

  mergefeaturecollections <- function(assetids) {
    # Start with an empty FeatureCollection
    mergedcollection <- ee$FeatureCollection(list())

    # Merge each asset into the mergedCollection
    for (assetid in assetids) {
      collection <- ee$FeatureCollection(assetid)
      mergedcollection <- mergedcollection$merge(collection)
    }
    return(mergedcollection)
  }

  combinedassets <- mergefeaturecollections(assetids)

  # Load the domain boundaries
  if (domain == "world") {
    domainboundaries <- ee$FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  } else if (domain == "conus") {
    # Load CONUS boundaries
    conus <- get_l48_boundary()
    sf_as_ee(
      x = conus,
      via = "getInfo_to_asset",
      assetId = paste0(ee, "/rgee/conus"),
      overwrite = TRUE,
      monitoring = TRUE,
      proj = "EPSG:4326"
    )
    # Replace with appropriate FeatureCollection for CONUS
    domainboundaries <- ee$FeatureCollection(paste0(ee, "/rgee/conus"))
  }

  # Calculate the distance to the nearest feature
  # If statement to determine distance value
  if (domain == "conus") {
    distancetoassets <- combinedassets$distance(ee$Number(2000252))
  } else if (type == "world") {
    distancetoassets <- combinedassets$distance(ee$Number(20002520))
  } else {
    stop("Invalid domain")
  }

  # Crop distanceToAssets with domainboundaries
  croppeddistance <- distancetoassets$clip(domainboundaries)

  # Convert the cropped image to UInt32
  croppeddistanceuint32 <- croppeddistance$toUint32()

  # Export to Google Drive
  exporttask <- ee$batch$Export$image$toDrive(
    image = croppeddistanceuint32,
    description = paste0("dist2", type, "_", domain),
    folder = "eu_distancee",
    fileNamePrefix = paste0("dist2", type, "_", domain),
    scale = 1000,
    fileFormat = "GeoTIFF",
    maxPixels = 1e13
  )
  exporttask$start()
  while (exporttask$status()$state != "COMPLETED") {
    Sys.sleep(300)
  }
}

#' @description Function to get roads data for the world (10m rnaturalearth)
#' @param path A character string specifying the path to write the roads
#' @param gee_path A character string specifying the path to write the roads
#' @param user A character string specifying the Google Earth Engine user name.
#' @return A euclidean distance raster with roads for the world.
#' @export

get_roads_global <- function(path, gee_path, user) {
  fn <- c(
    paste0(path, "Original/ne_roads/ne_10m_roads.gpkg"),
    paste0(path, "Raster/Global/ne_roads/ne_roads_distance_world_30s.tif")
  )
  if (!all(file.exists(fn))) {
    roads <- rnaturalearth::ne_download(type = "roads", scale = 10)
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    sf::st_write(roads, dsn = fn[1])
    dist_networks(user, lines = roads, domain = "world", type = "roads")
    roads <- list.files(gee_path, pattern = "dist2roads_world",
                        full.names = TRUE)
    roads <- terra::sprc(roads)
    roads <- terra::mosaic(roads)
    names(roads) <- "dist2roads"
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(roads,
      filename = fn[2], overwrite = TRUE, datatype = "INT4U",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    roads <- terra::rast(fn[2])
  }
  return(roads)
}

#' @description Function to download railroads for the world (10m rnaturalearth)
#' @param path A character string specifying the path to write the railroads
#' @param gee_path A character string specifying the path to write the railroads
#' @param user A character string specifying the Google Earth Engine user name.
#' @return A euclidean distance raster with railroads for the world.
#' @export

get_rails_global <- function(path, gee_path, user) {
  fn <- c(
    paste0(path, "Original/ne_railroads/ne_10m_railroads.gpkg"),
    paste0(path, "Raster/Global/ne_railroads/ne_railroads_dist_world_30s.tif")
  )
  if (!all(file.exists(fn))) {
    railroads <- rnaturalearth::ne_download(type = "railroads", scale = 10)
    dir.create(dirname(fn[1]), showWarnings = FALSE)
    sf::st_write(railroads, dsn = fn[1])
    dist_networks(user, lines = railroads, domain = "world", type = "rails")
    railroads <- list.files(gee_path, pattern = "dist2rails_world",
                            full.names = TRUE)
    railroads <- terra::sprc(railroads)
    railroads <- terra::mosaic(railroads)
    names(railroads) <- "dist2rails"
    dir.create(dirname(fn[2]), showWarnings = FALSE, recursive = TRUE)
    terra::writeRaster(railroads,
      filename = fn[2], overwrite = TRUE, datatype = "INT4U",
      gdal = "COMPRESS=ZSTD"
    )
  } else {
    railroads <- terra::rast(fn[2])
  }
  return(railroads)
}

#' @description Function to get soil pH & wc 33kpa & wc 1500kpa for the world
#' (250m OpenLandMap)
#' @param path A character string specifying the path to write soil vars
#' @param return A list pH & wc 33kpa & wc 1500kpa rasters for the world
#' @export

get_soilvars_global <- function(path) {
  # Create file names
  depths <- c("b0..0cm", "b10..10cm", "b30..30cm", "b60..60cm", "b100..100cm",
              "b200..200cm")
  orig_root <- "Original/OpenLandMap/soils/"
  fn <- c(paste0(path, orig_root, "sol_ph.h2o_usda.4c1a2a_m_250m_", depths,
                 "_1950..2017_v0.2.tif"),
  paste0(path, orig_root,"sol_watercontent.33kPa_usda.4b1c_m_250m_", depths, "_1950..2017_v0.1.tif"), # nolint: line_length_linter
  paste0(path, orig_root,"sol_watercontent.1500kPa_usda.3c2a1a_m_250m_", depths, "_1950..2017_v0.1.tif"), # nolint: line_length_linter
  paste0(path, "Raster/Global/OpenLandMap/soils/sol_ph.h2o_usda.4c1a2a_m_250m_mean_depth_1950..2017_v0.2.tif"), # nolint: line_length_linter
  paste0(path, "Raster/Global/OpenLandMap/soils/sol_watercontent.33kPa_usda.4b1c_m_250m_mean_depth_1950..2017_v0.1.tif"), # nolint: line_length_linter
  paste0(path, "Raster/Global/OpenLandMap/soils/sol_watercontent.1500kPa_usda.3c2a1a_m_250m_mean_depth_1950..2017_v0.1.tif")) # nolint: line_length_linter
  # check if file exists
  if (!all(file.exists(fn))) {
    # Download soil pH and water content
    dl_urls <- c("https://doi.org/10.5281/zenodo.2525664",
                 "https://doi.org/10.5281/zenodo.2784001")
    # Get list of files
    doi <- gsub("https://doi.org/", "", dl_urls)
    names(dl_urls) <- c("ph", "watercontent")
    lapply(seq_along(dl_urls), function(x) {
      files <- zen4R::get_zenodo(doi[x])
      files <- files$files
      files <- unlist(lapply(seq_along(files), function(y) {
        files[[y]]$filename[1]
      }))
      files <- files[grepl(".tif$", files)]
      files <- files[!grepl("_md_", files)]
      dir.create(paste0(tempdir(), "/OpenLandMap/soils/", names(dl_urls[x])),
                 showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(path, "Original/OpenLandMap/soils/", names(dl_urls[x])),
                 showWarnings = FALSE, recursive = TRUE)
      zen4R::download_zenodo(dl_urls[x],
        path = paste0(tempdir(), "/OpenLandMap/soils/", names(dl_urls)[x]),
        timeout = 60 * 30, files = list(files)
      )
    })
    # List files of .tifs in tempdir
    files <- list.files(paste0(tempdir(), "/OpenLandMap/soils/"),
                        full.names = TRUE, recursive = TRUE)
    files <- files[grepl(".tif$", files)]
    files <- sort(files)
    dir.create(dirname(fn[1]), showWarnings = FALSE, recursive = TRUE)
    # For each file, check if pH or water content, & write to correct folder
    for (x in seq_along(files)) {
      filename <- paste0(dirname(fn[1]), "/", basename(files[x]))
      terra::writeRaster(terra::rast(files[x]), overwrite = TRUE,
                         gdal = "COMPRESS=ZSTD", filename = filename,
                         datatype = "INT1U")
    }
    soilvars <- list.files(dirname(fn[1]), full.names = TRUE, recursive = TRUE)
    dir.create(dirname(fn[grep("Raster", fn)][1]), showWarnings = FALSE,
               recursive = TRUE)
    # Compute mean for ph and water content
    variables <- c("_ph", "33kPa", "1500kPa")
    output_filenames <- c("sol_ph.h2o_usda.4c1a2a_m_250m_mean_depth_1950..2017_v0.2.tif", # nolint: line_length_linter
                          "sol_watercontent.33kPa_usda.4b1c_m_250m_mean_depth_1950..2017_v0.1.tif", # nolint: line_length_linter
                          "sol_watercontent.1500kPa_usda.3c2a1a_m_250m_mean_depth_1950..2017_v0.1.tif") # nolint: line_length_linter
    soilvars <- lapply(seq_along(variables), function(i) {
      var <- variables[i]
      raster <- terra::rast(soilvars[grepl(var, soilvars)])
      raster <- terra::mean(raster, na.rm = TRUE)
      names(raster) <- var
      terra::writeRaster(raster, overwrite = TRUE, gdal = "COMPRESS=ZSTD",
                         filename = paste0(path, "Raster/Global/OpenLandMap/soils/", output_filenames[i]), # nolint: line_length_linter
                         datatype = "INT1U")
      return(raster)
    })
  } else {
    soilvars <- lapply(fn, terra::rast)
  }
  return(soilvars)
}

#' @description Function to reproject, resample, and extend a raster to a base
#' raster; and in the case of the USA, coarsen the raster to common
#' resolutions: 100m, 250m, 500m, and 1000m.
#' @param pred The predictor raster to reproject, resample, and extend.
#' @param domain A character string of Global or USA
#' @param path A character string specifying the path to write the reprojected
#' and resampled rasters.
#' @param base The base raster to reproject and resample the predictor to.
#' @param crop A logical specifying whether to crop the raster to the L48
#' boundary BEFORE reprojecting and resampling. Only in cases where the
#' predictor is only available at the global resolution.
#' @return A print statement of the time taken to reproject and resample the
#' predictor raster.
#' @export

match_to_base <- function(pred, base, domain, path, crop) {
  filename <- names(pred)
  dir <- gsub(".*\\/", "", dirname(terra::sources(pred)))
  filename <- paste0(dir, "/", dir, "_",filename)
  if (domain == "Global" && crop == FALSE) {
    filename <- paste0(path, "Raster/Global/", filename, "_reproj_1000m.tif")
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  } else if (domain == "USA" && crop == FALSE) {
    filename <- paste0(path, "Raster/USA/", filename, "_reproj_",c(30,100,250,500,1000),"m.tif")
    dir.create(dirname(filename[1]), showWarnings = FALSE, recursive = TRUE)
  } else if (domain == "Global" && crop == TRUE) {
    filename <- c(paste0(path, "Raster/Global/", filename, "_reproj_1000m.tif"), 
    paste0(path, "Raster/USA/", filename, "_reproj_",c(30,100,250,500,1000),"m.tif"))
    dir.create(dirname(filename[1]), showWarnings = FALSE, recursive = TRUE)
    dir.create(dirname(filename[6]), showWarnings = FALSE, recursive = TRUE)
  }

  if (!all(file.exists(filename))) {
    # If datatype is INT1U method is nearest neighbor, else bilinear
    d_type <- terra::datatype(pred)
    method <- ifelse(d_type == "INT1U", "near", "bilinear")

    # If crop is TRUE, crop to l48 boundary
    if (crop == TRUE) {
      l48 <- get_l48_boundary()
      l48 <- terra::vect(l48)
      l48 <- terra::as.polygons(terra::ext(l48) * 1.1)
      l48 <- terra::project(l48, pred)
      pred <- terra::crop(pred, l48)
    } else {
      pred <- pred
    }

    # Reproject and resample
    start_time <- Sys.time()
    r <- terra::project(pred, terra::crs(base), method = method)
    terra::resample(r, base, method = method, filename = filename, overwrite = TRUE,
                    wopt = list(datatype = d_type, gdal = "COMPRESS=ZSTD"))
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    msg <- print(paste("Time taken for", basename(filename), ":", time_taken))
  } else {
    msg <- print(paste0(basename(filename), " has already been created"))
  }
  return(msg)
}

#' @description Function to get all the get functions
#' @param path A character string specifying the path to write the rasters
#' @export

reproject_and_extend <- function(path) {
  get_functions <- c("get_biovars_global", "get_topo_global",
                     "get_landcover_global", "get_gdd_global",
                     "get_prectiming_global", "get_pop_global",
                     "get_roads_global", "get_rails_global",
                     "get_soilvars_global", "get_biovars_conus",
                     "get_topo_conus", "get_landcover_conus")

  get_functions <- data.frame(do.call(rbind, strsplit(get_functions, "_")))
  colnames(get_functions) <- c("get", "type", "domain")
  types <- unique(get_functions$type)
  base_global <- terra::rast(paste0(path, "Raster/Global/baserasters_pops_sdm/base_1000m_world.tif")) # nolint: line_length_linter
  base_conus <- terra::rast(paste0(path, "Raster/USA/baserasters_pops_sdm/base_30m_conus.tif")) # nolint: line_length_linter
  for (i in seq_along(types)) {
    dt <- get_functions[get_functions$type == types[i], ]
    if (nrow(dt) > 1) {
      global <- dt[dt$domain == "global", ]
      global <- paste0(global$get, "_", global$type, "_", global$domain)
      global <- get(global)
      pred <- global(path)
      for (j in seq_along(pred)) {
        match_to_base(pred[[j]], base_global, "Global", path, crop = FALSE)
      }
      conus <- dt[dt$domain == "conus", ]
      conus <- paste0(conus$get, "_", conus$type, "_", conus$domain)
      conus <- get(conus)
      pred <- conus(path)
      for (j in seq_along(pred)) {
        match_to_base(pred[[j]], base_conus, "USA", path, crop = FALSE)
      }
    } else if (nrow(dt) == 1 && dt$domain == "global") {
      global <- paste0(dt$get, "_", dt$type, "_", dt$domain)
      global <- get(global)
      pred <- global(path)
      for (j in seq_along(pred)) {
        match_to_base(pred[[j]], base_global, "Global", path, crop = FALSE)
        match_to_base(pred[[j]], base_conus, "USA", path, crop = TRUE)
      }
    } else if (nrow(dt) == 1 && dt$domain == "conus") {
      conus <- paste0(dt$get, "_", dt$type, "_", dt$domain)
      conus <- get(conus)
      pred <- conus(path)
      for (j in seq_along(pred)) {
        match_to_base(pred[[j]], base_conus, "USA", path, crop = FALSE)
      }
    }
  }
}

reproject_and_extend(path = "D:/blaginh/pops.sdm/Data/")
