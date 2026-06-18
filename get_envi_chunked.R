invisible(lapply(c(
  "geodata", "terra", "sbtools", "rnaturalearth", "zen4R", "rgee",
  "foreach", "doParallel", "reticulate"
), require, character.only = TRUE))

py_require("earthengine-api")

.zstd <- "COMPRESS=ZSTD"
.gdal_big <- c("TILED=YES", "BIGTIFF=YES")
.gdal_big_zstd <- c("COMPRESS=ZSTD", "ZSTD_LEVEL=1", "TILED=YES", "BIGTIFF=YES")
.month_days <- c(31, 28.2425, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

.makedir <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)
.rasts <- function(x) lapply(x, terra::rast)

.write <- function(x, filename, datatype = "FLT4S", gdal = .zstd, ...) {
  .makedir(dirname(filename[1]))
  terra::writeRaster(
    x,
    filename = filename, overwrite = TRUE,
    datatype = datatype, gdal = gdal, ...
  )
}

.cluster <- function(n) {
  cl <- parallel::makeCluster(n)
  doParallel::registerDoParallel(cl)
  cl
}

.stopcl <- function(cl) try(parallel::stopCluster(cl), silent = TRUE)

.need <- function(x, nm) {
  if (is.null(x)) stop(nm, " is needed when cached files don't already exist.")
  x
}

#' Determine the smallest safe terra datatype for a raster layer
#'
#' @param raster A single-layer `terra::SpatRaster`.
#' @return A character string giving the terra datatype.
#' @export
get_datatype <- function(raster) {
  if (terra::nlyr(raster) > 1) stop("Input raster must have only one layer.")

  r <- terra::global(raster, fun = "range", na.rm = TRUE, par = TRUE)
  lo <- r[1, 1]
  hi <- r[1, 2]
  samples <- terra::spatSample(raster, size = 100)
  samples <- samples[!is.na(samples)]
  is_int <- all(samples == floor(samples))

  if (lo >= 0 && hi <= 255 && is_int) {
    dt <- "INT1U"
  } else if (lo >= -128 && hi <= 127 && is_int) {
    dt <- "INT1S"
  } else if (lo >= -32768 && hi <= 32767 && is_int) {
    dt <- "INT2S"
  } else if (lo >= -2147483648 && hi <= 2147483647 && is_int) {
    dt <- "INT4S"
  } else if (lo >= 0 && hi <= 65535 && is_int) {
    dt <- "INT2U"
  } else if (lo >= 0 && hi <= 4294967295 && is_int) {
    dt <- "INT4U"
  } else {
    dt <- "FLT4S"
  }
  assigned_dt <- terra::datatype(raster)

  cat(
    "Value range of raster layer:", lo, "to", hi, "\n",
    "Raster-assigned data type:", assigned_dt, "\n",
    "Value-based data type:", dt, "\n"
  )

  dt <- data.frame(
    dt = dt, min = lo, max = hi,
    assigned_dt = assigned_dt, source = terra::sources(raster)
  )
}

#' Get global WorldClim bioclimatic variables
#'
#' @param path Path to the data directory.
#' @return A `SpatRaster` or list of 19 `SpatRaster` layers.
#' @export
get_biovars_global <- function(path) {
  fn <- paste0(path, "Original/wc2.1_30s_bio/wc2.1_30s_bio_", 1:19, ".tif")
  if (all(file.exists(fn))) {
    return(.rasts(fn))
  }

  x <- geodata::worldclim_global(var = "bio", res = 0.5, path = tempdir())

  lapply(seq_len(terra::nlyr(x)), function(i) {
    .write(x[[i]], file.path(dirname(fn[1]), paste0(names(x)[i], ".tif")))
  })

  x
}

#' Get CONUS USGS biocomposite bioclimatic variables
#'
#' @param path Path to the data directory.
#' @return A list of 19 `SpatRaster` layers.
#' @export
get_biovars_conus <- function(path) {
  nm <- c(paste0("bio", 1:4), "bio4a", paste0("bio", 5:19))
  fn <- paste0(
    path,
    "Original/BioClimComposite/BioClimComposite_1971_2000_400m_",
    nm,
    ".tif"
  )

  if (all(file.exists(fn))) {
    return(.rasts(fn)[-5])
  }

  sbtools::item_file_download(
    sb_id = "4ff32906e4b0e183ef5a2f16",
    dest_dir = tempdir(),
    names = "BioClimComposite_1971_2000_400m.tif",
    destinations = file.path(tempdir(), "BioClimComposite_1971_2000_400m.tif"),
    overwrite_file = TRUE
  )

  x <- terra::rast(file.path(tempdir(), "BioClimComposite_1971_2000_400m.tif"))
  names(x) <- nm

  lapply(seq_len(terra::nlyr(x)), function(i) {
    .write(
      x[[i]],
      file.path(
        dirname(fn[1]),
        paste0("BioClimComposite_1971_2000_400m_", names(x)[i], ".tif")
      )
    )
  })

  x[-5]
}

#' Get downscaled CONUS bioclimatic variables
#'
#' @param path Path to the data directory.
#' @return A list of 19 `SpatRaster` layers.
#' @export
get_biovars_conus_ds <- function(path) {
  .rasts(paste0(path, "Raster/USA/bioclimatic/bio", 1:19, ".tif"))
}

#' Get lower-48 United States boundary
#'
#' @return An `sf` geometry for CONUS.
#' @export
get_l48_boundary <- function() {
  locs <- rnaturalearth::ne_states(
    country = "United States of America",
    returnclass = "sf"
  )

  locs <- unique(locs["iso_3166_2"])
  locs <- locs[!locs$iso_3166_2 %in% c(
    "US-PR", "US-VI", "US-GU", "US-MP", "US-AS", "US-AK", "US-HI"
  ), ]

  sf::st_union(locs)
}

#' Get global land boundary
#'
#' @return An `sf` geometry for global land.
#' @export
get_global_boundary <- function() {
  locs <- rnaturalearth::ne_countries(
    scale = 50,
    returnclass = "sv"
  )
  # Get rid of antarctica
  locs <- locs[locs$continent != "Antarctica", ]
  locs <- terra::aggregate(locs, dissolve = TRUE)
  locs <- sf::st_as_sf(locs)
}

#' Create slope, aspect, and hillshade rasters
#'
#' @param path Path prefix for output terrain files.
#' @param elevation A `SpatRaster` elevation layer.
#' @return A list containing slope and aspect rasters.
#' @export
get_terrain <- function(path, elevation) {
  slp <- terra::terrain(elevation, v = "slope", unit = "radians", neighbors = 8)
  asp <- terra::terrain(elevation,
    v = "aspect", unit = "radians",
    neighbors = 4
  )
  shd <- terra::shade(slp, asp, direction = 180)

  names(slp) <- "slope"
  names(asp) <- "aspect"
  names(shd) <- "hillshade"

  .write(slp, paste0(path, "slope.tif"))
  .write(asp, paste0(path, "aspect.tif"))
  .write(shd, paste0(path, "hillshade.tif"))

  list(slp, asp)
}

#' Get global topographic rasters
#'
#' @param path Path to the data directory.
#' @return A list containing elevation, slope, and aspect rasters.
#' @export
get_topo_global <- function(path) {
  fn <- c(
    paste0(path, "Original/wc2.1_30s_elev/wc2.1_30s_elev.tif"),
    paste0(path, "Raster/Global/wc2.1_30s_elev/wc2.1_30s_slope.tif"),
    paste0(path, "Raster/Global/wc2.1_30s_elev/wc2.1_30s_aspect.tif")
  )

  if (all(file.exists(fn))) {
    return(.rasts(fn))
  }

  e <- geodata::elevation_global(res = 0.5, path = tempdir())
  .write(e, fn[1], datatype = "INT2S")

  sa <- get_terrain(dirname(fn[2]), e)
  list(e, sa[[1]], sa[[2]])
}


#' Get CONUS SRTM topographic rasters
#'
#' @param path Path to the data directory.
#' @param key OpenTopography API key.
#' @return A list containing elevation, slope, and aspect rasters.
#' @export
get_topo_conus <- function(path, key = NULL) {
  fn <- c(
    paste0(path, "Original/SRTMGL1/SRTMGL1_1s_l48.tif"),
    paste0(path, "Raster/USA/SRTMGL1/SRTMGL1_1s_l48_slope.tif"),
    paste0(path, "Raster/USA/SRTMGL1/SRTMGL1_1s_l48_aspect.tif")
  )

  if (all(file.exists(fn))) {
    return(.rasts(fn))
  }

  key <- .need(key, "key")
  loc <- get_l48_boundary()
  bbs <- sf::st_make_grid(sf::st_transform(loc, 3857), cellsize = 670820)
  bbs <- lapply(sf::st_transform(bbs, 4326), sf::st_bbox)

  tile_dir <- file.path(dirname(fn[1]), "tiles")
  .makedir(tile_dir)

  lapply(seq_along(bbs), function(i) {
    bb <- bbs[[i]]

    url <- paste0(
      "https://portal.opentopography.org/API/globaldem?demtype=SRTMGL1",
      "&south=", bb[2],
      "&north=", bb[4],
      "&west=", bb[1],
      "&east=", bb[3],
      "&outputFormat=GTiff&API_Key=", key
    )

    dest <- file.path(tile_dir, paste0("srtm_", i, ".tif"))
    dl <- curl::curl_fetch_memory(url)$url

    if (.Platform$OS.type == "windows") {
      download.file(dl, dest, mode = "wb")
    } else {
      download.file(dl, dest)
    }
  })

  srtm <- terra::mosaic(
    terra::sprc(list.files(tile_dir, pattern = "\\.tif$", full.names = TRUE)),
    wopt = list(datatype = "INT2S")
  )

  names(srtm) <- "elevation"

  srtm <- terra::crop(
    srtm,
    terra::vect(loc),
    mask = TRUE,
    filename = fn[1],
    overwrite = TRUE,
    wopt = list(datatype = "INT2S", gdal = .zstd)
  )

  unlink(tile_dir, recursive = TRUE)

  sa <- get_terrain(paste0(dirname(fn[2]), "/SRTMGL1_1s_l48_"), srtm)
  list(srtm, sa[[1]], sa[[2]])
}

.base_rasters <- function(
    path, domain, res_names, res_mult,
    topo_fun, stem, key = NULL) {
  root <- paste0(path, "Raster/", domain, "/pops_sdm_base/")
  fn <- paste0(root, "base_", res_names, "m_", stem, "_epsg4326.tif")

  if (!all(file.exists(fn))) {
    elev <- if (identical(stem, "conus")) {
      topo_fun(path, key = key)[[1]]
    } else {
      topo_fun(path)[[1]]
    }
    elev <- terra::clamp(elev, 1, 1, datatype = "INT1U")

    .makedir(root)

    base_res <- terra::res(elev)[1]

    lapply(seq_along(res_names), function(i) {
      terra::project(
        elev,
        y = "EPSG:4326",
        threads = TRUE,
        method = "near",
        filename = fn[i],
        res = base_res * res_mult[i],
        overwrite = TRUE,
        wopt = list(gdal = .zstd, datatype = "INT1U")
      )
    })

    message(
      "Base rasters created for ", stem, " at resolutions: ",
      paste0(res_names, collapse = ", "), "m"
    )
  } else {
    message(
      "Base rasters already exist for ", stem, " at resolutions: ",
      paste0(res_names, collapse = ", "), "m"
    )
  }

  .rasts(fn)
}

#' Create CONUS base rasters
#'
#' Creates a 1 arc-second base raster for CONUS and common resolutions:
#' 30m, 100m, 250m, 500m, and 1000m.
#'
#' @param path Path to the data directory.
#' @param key OpenTopography API key, needed only if CONUS elevation is absent.
#' @return A list of CONUS base rasters.
#' @export
base_conus <- function(path, key = NULL) {
  .base_rasters(
    path = path,
    domain = "USA",
    res_names = c(30, 100, 250, 500, 1000),
    res_mult = c(1, 3, 7.5, 15, 30),
    topo_fun = get_topo_conus,
    stem = "conus",
    key = key
  )
}

#' Create global base rasters
#'
#' Creates a 30 arc-second base raster for the world and common resolutions:
#' 1000m, 2500m, and 5000m.
#'
#' @param path Path to the data directory.
#' @return A list of global base rasters.
#' @export
base_global <- function(path) {
  .base_rasters(
    path = path,
    domain = "Global",
    res_names = c(1000, 2500, 5000),
    res_mult = c(1, 2.5, 5),
    topo_fun = get_topo_global,
    stem = "global"
  )
}

#' Get global landcover rasters
#'
#' @param path Path to the data directory.
#' @return A list of landcover rasters.
#' @export
get_landcover_global <- function(path) {
  types <- c("built", "cropland", "grassland", "shrubs", "trees", "wetland")
  fn <- paste0(path, "Original/WorldCover/", types, "_30s.tif")

  if (all(file.exists(fn))) {
    return(.rasts(fn))
  }

  x <- lapply(types, geodata::landcover, path = tempdir())

  lapply(seq_along(fn), function(i) {
    .write(x[[i]], fn[i], datatype = "INT1U")
  })
}

#' Get CONUS NLCD landcover raster
#'
#' @param path Path to the data directory.
#' @return A landcover `SpatRaster` or list containing the raster.
#' @export
get_landcover_conus <- function(path) {
  fn <- c(
    paste0(path, "Original/nlcd/nlcd_2021_land_cover_l48_20230630.tif"),
    paste0(path, "Raster/USA/nlcd/nlcd_2021_land_cover_l48.tif")
  )

  if (all(file.exists(fn))) {
    return(.rasts(fn[-1]))
  }

  url <- paste0(
    "https://s3-us-west-2.amazonaws.com/mrlc/",
    "nlcd_2021_land_cover_l48_20230630.zip"
  )
  zip <- file.path(tempdir(), basename(url))

  options(timeout = 60 * 5)
  download.file(url, zip, mode = "wb")
  unzip(zip, exdir = file.path(tempdir(), "nlcd"))

  x <- terra::rast(file.path(
    tempdir(),
    "nlcd/nlcd_2021_land_cover_l48_20230630.img"
  ))
  .write(x, fn[1], datatype = "INT1U")

  rcl <- matrix(c(
    0, 20, 0, 21, 24, 1, 25, 40, 0, 41, 41, 2, 42, 42, 3, 43, 43, 4,
    44, 51, 0, 52, 52, 5, 53, 70, 0, 71, 71, 6, 72, 80, 0, 81, 81, 7,
    82, 82, 8, 83, 89, 0, 90, 90, 9, 91, 94, 0, 95, 95, 9, 96, 255, 0
  ), ncol = 3, byrow = TRUE)

  x <- terra::classify(x, rcl, right = NA)
  .write(x, fn[2], datatype = "INT1U")
}

#' Calculate growing degree days
#'
#' @param tavg_files Monthly temperature raster paths or 12-layer raster path.
#' @param fn Output filename.
#' @param ncore Number of cores.
#' @param tbase Base temperature.
#' @param layered Logical; whether `tavg_files` is one multi-layer raster.
#' @param memfrac terra memory fraction per worker.
#' @param out_gdal GDAL options for the output raster.
#' @return A growing degree day `SpatRaster`.
#' @keywords internal
calc_gdd <- function(tavg_files, fn, ncore = 6, tbase = 5, layered = FALSE,
                     memfrac = NULL, out_gdal = .gdal_big) {
  .makedir(dirname(fn))

  tmp <- file.path(dirname(fn), "tmp_gdd_monthly")
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE)
  .makedir(tmp)

  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  fn_m <- file.path(tmp, paste0("gdd_", sprintf("%02d", 1:12), ".tif"))
  nw <- min(ncore, 12, max(1, parallel::detectCores() - 1))
  mf <- if (is.null(memfrac)) min(0.25, 0.9 / nw) else memfrac

  cl <- .cluster(nw)
  on.exit(.stopcl(cl), add = TRUE)

  foreach::foreach(i = seq_len(12), .packages = "terra") %dopar% {
    terra::terraOptions(
      memfrac = mf, todisk = TRUE,
      progress = 0, tempdir = tmp
    )

    x <- if (layered) {
      terra::rast(tavg_files)[[i]]
    } else {
      terra::rast(tavg_files[i])
    }

    terra::writeRaster(
      terra::classify(x - tbase, rbind(c(-Inf, 0, 0))) * .month_days[i],
      filename = fn_m[i],
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = .gdal_big
    )

    rm(x)
    terra::tmpFiles(remove = TRUE)
    gc()
    NULL
  }

  gdd <- sum(
    terra::rast(fn_m),
    na.rm = FALSE,
    filename = fn,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = out_gdal,
      names = "gdd_tbase5",
      progress = 1
    )
  )

  names(gdd) <- "gdd_tbase5"
  gdd
}

#' Get global growing degree days
#'
#' @param path Path to the data directory.
#' @param ncore Number of cores.
#' @return A growing degree day `SpatRaster`.
#' @export
get_gdd_global <- function(path, ncore = 6) {
  fn_tavg <- paste0(path, "Original/wc2.1_30s_tavg/wc2.1_30s_tavg.tif")
  fn_gdd <- paste0(path, "Raster/Global/wc2.1_30s_gdd/wc2.1_30s_gdd.tif")

  if (file.exists(fn_gdd)) {
    return(terra::rast(fn_gdd))
  }

  if (!file.exists(fn_tavg)) {
    x <- geodata::worldclim_global(var = "tavg", res = 0.5, path = tempdir())
    .write(x, fn_tavg, datatype = "FLT4S", gdal = .gdal_big_zstd)
  }

  if (terra::nlyr(terra::rast(fn_tavg)) != 12) {
    stop("Expected 12 monthly tavg layers.")
  }

  calc_gdd(
    tavg_files = fn_tavg,
    fn = fn_gdd,
    ncore = ncore,
    layered = TRUE,
    memfrac = 0.15,
    out_gdal = .gdal_big_zstd
  )
}

#' Get downscaled CONUS growing degree days
#'
#' @param path Path to the data directory.
#' @param ncore Number of cores.
#' @return A growing degree day `SpatRaster`.
#' @export
get_gdd_conus_ds <- function(path, ncore = 6) {
  fn <- paste0(path, "Raster/USA/gdd/gdd.tif")

  if (file.exists(fn)) {
    return(terra::rast(fn))
  }

  files <- list.files(
    paste0(path, "Raster/USA/daymet/tmean/downscaled"),
    pattern = "\\.tif$",
    full.names = TRUE
  )

  month_id <- as.integer(gsub(".*?([0-9]{1,2})\\.tif$", "\\1", basename(files)))

  if (any(is.na(month_id))) {
    stop("Could not infer month number from one or more tmean filenames.")
  }

  o <- order(month_id)
  files <- files[o]
  month_id <- month_id[o]

  if (length(files) != 12 || !identical(month_id, 1:12)) {
    stop("Expected exactly one tmean raster for each month 1:12.")
  }

  calc_gdd(files, fn, ncore = ncore)
}

#' Get global population density
#'
#' @param path Path to the data directory.
#' @return A population density `SpatRaster`.
#' @export
get_pop_global <- function(path) {
  fn <- file.path(
    path,
    "Original",
    "gpw_v4_population_density",
    "gpw_v4_population_density_rev11_2020_30s.tif"
  )

  if (file.exists(fn)) {
    return(terra::rast(fn))
  }

  x <- geodata::population(year = 2020, res = 0.5, path = tempdir())
  .write(x, fn)
  x
}

#' Get CONUS population density
#'
#' @param path Path to the data directory.
#' @return A list containing a population density raster.
#' @export
get_pop_conus <- function(path) {
  .rasts(paste0(path, "Original/meta_population/meta_population.tif"))
}
#' Get network vector data
#'
#' Downloads or reads cached vector network data for either roads or rails.
#' For `domain = "USA"`, this uses TIGER/Line data. For `domain = "Global"`,
#' this uses Natural Earth data.
#'
#' @param path Path to the data directory.
#' @param network_type Network type. One of `"roads"` or `"rails"`.
#' @param domain Spatial domain. One of `"USA"` or `"Global"`.
#' @param year TIGER/Line vintage year. Used only when `domain = "USA"`.
#' @return A network `SpatVector`.
#' @export
get_network <- function(
    path,
    network_type = c("rails", "roads"),
    domain = c("USA", "Global"),
    year = 2024) {
  network_type <- match.arg(network_type)
  domain <- match.arg(domain)

  if (domain == "USA") {
    dirs <- list(
      rails = file.path(path, "Original/tiger_rails"),
      roads = file.path(path, "Original/tiger_roads")
    )

    dir.create(dirs[[network_type]], recursive = TRUE, showWarnings = FALSE)

    out <- file.path(
      dirs[[network_type]],
      sprintf(
        if (network_type == "rails") {
          "tiger_railroads_l48_%s.gpkg"
        } else {
          "tiger_roads_l48_%s.gpkg"
        },
        year
      )
    )

    if (file.exists(out)) {
      return(terra::vect(out))
    }

    if (network_type == "rails") {
      b <- get_l48_boundary() |>
        sf::st_transform(5070) |>
        sf::st_buffer(200000)

      x <- tigris::rails(
        filter_by = sf::st_bbox(sf::st_transform(b, 4326)),
        year = year
      ) |>
        sf::st_transform(5070) |>
        sf::st_intersection(b)
    } else {
      states <- c(
        "AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID",
        "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI",
        "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
        "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD",
        "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY", "DC"
      )

      x <- do.call(rbind, lapply(states, \(s) {
        f <- file.path(dirs$roads, sprintf("tiger_roads_%s_%s.gpkg", s, year))

        if (!file.exists(f)) {
          tigris::primary_secondary_roads(s, year = year, progress = FALSE) |>
            sf::st_transform(5070) |>
            sf::st_write(f, delete_dsn = TRUE, quiet = TRUE)
        }

        sf::st_read(f, quiet = TRUE)
      }))
      file.remove(file.path(
        dirs$roads,
        sprintf("tiger_roads_%s_%s.gpkg", states, year)
      ))
    }

    sf::st_write(x, out, delete_dsn = TRUE, quiet = TRUE)
    return(terra::vect(out))
  }

  if (domain == "Global") {
    dirs <- list(
      rails = file.path(path, "Original/naturalearth_rails"),
      roads = file.path(path, "Original/naturalearth_roads")
    )

    dir.create(dirs[[network_type]], recursive = TRUE, showWarnings = FALSE)

    ne_type <- c(rails = "railroads", roads = "roads")[[network_type]]
    out <- file.path(dirs[[network_type]], sprintf("ne_%s.gpkg", ne_type))

    if (!file.exists(out)) {
      x <- rnaturalearth::ne_download(
        scale = 10,
        type = ne_type,
        category = "cultural",
        returnclass = "sf"
      )

      sf::st_write(x, out, delete_dsn = TRUE, quiet = TRUE)
    }

    return(terra::vect(out))
  }
}

#' Get network-distance raster
#'
#' Creates or reads a cached distance-to-network raster for roads or rails.
#' For `domain = "USA"`, the raster is based on TIGER/Line data and a 30 m
#' CONUS base raster. For `domain = "Global"`, the raster is based on Natural
#' Earth data and a 1000 m global base raster.
#'
#' @param path Path to the data directory.
#' @param network_type Network type. One of `"roads"` or `"rails"`.
#' @param domain Spatial domain. One of `"USA"` or `"Global"`.
#' @param year TIGER/Line vintage year. Used only when `domain = "USA"`.
#' @return A network-distance `SpatRaster`.
#' @export
get_network_distance <- function(
    path,
    network_type = c("rails", "roads"),
    domain = c("USA", "Global"),
    year = 2024) {
  network_type <- match.arg(network_type)
  domain <- match.arg(domain)
  cat("Processing distance to", network_type, "for", domain, "...\n")

  cfg <- list(
    USA = list(
      folder = "USA",
      res = 30,
      stem = "conus",
      crs = "EPSG:5070",
      topo_fun = get_topo_conus,
      boundary_fun = get_l48_boundary,
      base_domain = "USA",
      prefix = "tiger"
    ),
    Global = list(
      folder = "Global",
      res = 1000,
      stem = "global",
      crs = "EPSG:8857",
      topo_fun = get_topo_global,
      boundary_fun = get_global_boundary,
      base_domain = "Global",
      prefix = "ne"
    )
  )[[domain]]

  base <- file.path(
    path,
    sprintf(
      "Raster/%s/pops_sdm_base/base_%sm_%s_%s.tif",
      cfg$folder,
      cfg$res,
      cfg$stem,
      sub("EPSG:", "", cfg$crs)
    )
  )

  if (!file.exists(base)) {
    dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)

    r <- .base_rasters(
      path = path,
      domain = cfg$base_domain,
      res_names = cfg$res,
      res_mult = 1,
      topo_fun = cfg$topo_fun,
      stem = cfg$stem,
      key = NULL
    )[[1]]

    terra::project(
      r,
      cfg$crs,
      method = "near",
      mask = TRUE,
      filename = base,
      overwrite = TRUE,
      wopt = list(gdal = .zstd, datatype = "INT1U"),
      threads = TRUE
    )
  }

  out <- file.path(
    path,
    sprintf(
      "Raster/%s/%s_%s/%s_dist_%sm_%s.tif",
      cfg$folder,
      cfg$prefix,
      network_type,
      network_type,
      cfg$res,
      sub("EPSG:", "", cfg$crs)
    )
  )

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(out)) {
    return(terra::rast(out))
  }

  net <- get_network(
    path = path,
    network_type = network_type,
    domain = domain,
    year = year
  )

  if (terra::crs(net) != cfg$crs) {
    net <- terra::project(net, cfg$crs)
  }

  tmp_dir <- file.path(path, "Raster", cfg$folder, paste0(
    cfg$prefix,
    "_", network_type
  ), "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  tmp_file <- file.path(
    tmp_dir,
    sprintf("%s_rasterized_%s.tif", network_type, sub("EPSG:", "", cfg$crs))
  )

  if (file.exists(tmp_file)) {
    r <- terra::rast(tmp_file)
  } else {
    r <- terra::rasterize(
      net,
      terra::rast(base),
      background = NA,
      touches = TRUE,
      filename = tmp_file,
      overwrite = TRUE,
      wopt = list(datatype = "INT1U")
    )
  }

  r <- terra::classify(r, matrix(c(0, NA), ncol = 2))

  d <- terra::distance(
    r,
    maxdist = if (domain == "USA") {
      c(roads = 200000, rails = 300000)[[network_type]]
    } else {
      c(roads = 5000000, rails = 5000000)[[network_type]]
    }
  )

  rm(r, net)
  gc()

  b <- cfg$boundary_fun() |>
    sf::st_transform(cfg$crs) |>
    terra::vect()

  terra::crop(
    d,
    b,
    mask = TRUE,
    filename = out,
    overwrite = TRUE,
    wopt = list(gdal = .zstd)
  )
  rm(d)
  gc()
}

#' Get global OpenLandMap soil variables
#'
#' @param path Path to the data directory.
#' @param ncores Number of cores.
#' @return A list of soil variable rasters.
#' @export
get_soilvars_global <- function(path, ncores = 2) {
  depths <- c(
    "b0..0cm", "b10..10cm", "b30..30cm",
    "b60..60cm", "b100..100cm", "b200..200cm"
  )
  og <- "Original/OpenLandMap_soils/"
  new <- "Raster/Global/OpenLandMap_soils/"

  outs <- c(
    "sol_ph_median_250m_1950_2017.tif",
    "sol_wc33kPa_median_250m_1950_2017.tif",
    "sol_wc1500kPa_median_250m_1950_2017.tif"
  )

  fn <- c(
    paste0(
      path, og, "sol_ph.h2o_usda.4c1a2a_m_250m_", depths,
      "_1950..2017_v0.2.tif"
    ),
    paste0(
      path, og, "sol_watercontent.33kPa_usda.4b1c_m_250m_",
      depths, "_1950..2017_v0.1.tif"
    ),
    paste0(
      path, og, "sol_watercontent.1500kPa_usda.3c2a1a_m_250m_",
      depths, "_1950..2017_v0.1.tif"
    ),
    paste0(path, new, outs)
  )

  if (all(file.exists(fn))) {
    return(.rasts(fn[grep("Raster", fn)]))
  }

  dois <- c(
    ph = "10.5281/zenodo.2525664",
    watercontent = "10.5281/zenodo.2784001"
  )

  for (doi in dois) {
    files <- suppressMessages(zen4R::get_zenodo(doi))$files
    files <- unlist(lapply(files, function(z) z$filename[1]))
    files <- files[grepl("tif$", files) & !grepl("_md_", files)]

    tmp <- file.path(path, og, "tmp")
    .makedir(tmp)
    .makedir(file.path(path, og))

    if (all(file.exists(file.path(path, og, basename(files))))) next

    cl <- parallel::makeCluster(pmin(length(files), ncores))
    on.exit(.stopcl(cl), add = TRUE)

    parallel::clusterExport(cl, c("doi", "tmp", "files"), envir = environment())

    zen4R::download_zenodo(
      doi,
      path = tmp,
      timeout = 60 * 30,
      parallel = TRUE,
      parallel_handler = parallel::parLapply,
      cl = cl,
      files = files,
      quiet = TRUE
    )

    .stopcl(cl)
  }

  files <- list.files(file.path(path, og), full.names = TRUE, recursive = TRUE)
  .makedir(dirname(fn[1]))

  for (f in files) {
    out <- file.path(dirname(fn[1]), basename(f))
    if (file.exists(out)) next

    r <- terra::rast(f)
    dt <- get_datatype(r)

    cat("Processing file:", basename(f), "with datatype", dt, "\n")

    .write(r, out, datatype = dt)
    file.remove(f)

    rm(r)
    gc()
  }

  soilvars <- list.files(dirname(fn[1]), full.names = TRUE, recursive = TRUE)
  vars <- c("_ph", "33kPa", "1500kPa")

  .makedir(file.path(path, new))

  cl <- .cluster(pmin(ncores, length(vars)))
  on.exit(.stopcl(cl), add = TRUE)

  foreach::foreach(
    i = seq_along(vars),
    .packages = "terra",
    .export = c("soilvars", "path", "new", "outs", "vars")
  ) %dopar% {
    output_file <- paste0(path, new, outs[i])

    if (file.exists(output_file)) next

    r <- terra::mean(
      terra::rast(soilvars[grepl(vars[i], soilvars)]),
      na.rm = TRUE
    )

    names(r) <- vars[i]

    terra::writeRaster(
      r,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD",
      filename = output_file,
      datatype = "FLT4S"
    )

    rm(r)
    gc()
  }

  .stopcl(cl)

  .rasts(paste0(path, new, outs))
}

#' Get CONUS POLARIS soil variables
#'
#' @param path Path to the data directory.
#' @return A list of soil variable rasters.
#' @export
get_soilvars_conus <- function(path) {
  .rasts(paste0(
    path,
    "Original/polaris_soils/",
    c("clay", "om", "ph", "sand", "silt", "theta_r", "theta_s"),
    "_mean_0_5.tif"
  ))
}

#' Get CONUS PAR raster
#'
#' @param path Path to the data directory.
#' @return A list containing the PAR raster.
#' @export
get_par_conus <- function(path) {
  .rasts(paste0(path, "Original/modis_par/modis_par.tif"))
}

#' Get CONUS Landsat NDVI raster
#'
#' @param path Path to the data directory.
#' @return A list containing the NDVI raster.
#' @export
get_landsat_ndvi_conus <- function(path) {
  .rasts(paste0(path, "Original/landsat_ndvi/landsat_ndvi.tif"))
}

#' Get CONUS Landsat EVI raster
#'
#' @param path Path to the data directory.
#' @return A list containing the EVI raster.
#' @export
get_landsat_evi_conus <- function(path) {
  .rasts(paste0(path, "Original/landsat_evi/landsat_evi.tif"))
}

#' Create output filenames for resampled predictors
#'
#' @param pred A predictor `SpatRaster`.
#' @param domain Either `"Global"` or `"USA"`.
#' @param path Path to the data directory.
#' @return A character vector of output filenames.
#' @export
get_filename <- function(pred, domain, path) {
  dir <- gsub(".*/", "", dirname(terra::sources(pred)))

  res <- switch(domain,
    Global = c(1000, 2500, 5000),
    USA = c(30, 100, 250, 500, 1000),
    stop("Invalid domain")
  )

  root <- switch(domain,
    Global = "Raster/Global/",
    USA = "Raster/USA/"
  )

  fn <- paste0(
    path,
    root,
    dir,
    "/",
    dir,
    "_",
    names(pred),
    "_reproj_",
    res,
    "m.tif"
  )

  .makedir(dirname(fn[1]))
  fn
}

#' Match a predictor raster to base rasters
#'
#' @param pred A predictor `SpatRaster`.
#' @param base A list of base rasters.
#' @param domain Either `"Global"` or `"USA"`.
#' @param path Path to the data directory.
#' @param crop Logical; TRUE == crop to the L48 boundary before matching.
#' @return A printed status message.
#' @export
match_to_base <- function(pred, base, domain, path, crop) {
  fn <- get_filename(pred, domain, path)

  if (all(file.exists(fn))) {
    return(print(paste0(basename(fn), " has already been created")))
  }

  d_type <- terra::datatype(pred)
  method <- ifelse(d_type == "INT1U", "near", "bilinear")

  if (isTRUE(crop)) {
    l48 <- terra::project(terra::vect(get_l48_boundary()), terra::crs(pred))
    pred <- terra::crop(pred, l48, mask = TRUE)
  }

  for (i in seq_along(base)) {
    start <- Sys.time()

    r <- terra::project(
      pred,
      base[[i]],
      method = method,
      threads = TRUE,
      wopt = list(datatype = d_type)
    )

    terra::crop(
      r,
      base[[i]],
      mask = TRUE,
      extend = TRUE,
      filename = fn[i],
      overwrite = TRUE,
      wopt = list(datatype = d_type, gdal = .zstd)
    )

    terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)

    elapsed <- Sys.time() - start

    print(paste(
      "Time taken for",
      basename(fn[i]),
      ":",
      elapsed,
      attr(elapsed, "units")
    ))
  }

  print(paste0(basename(fn), " has been created"))
}

#' Get predictor rasters by domain
#'
#' @param path Path to the data directory.
#' @param domain Either `"Global"` or `"USA"`.
#' @param key OpenTopo API key needed if CONUS elevation is absent.
#' @return A named list of predictor rasters.
#' @export
get_rasters <- function(path, domain = "Global", key = NULL) {
  if (domain == "Global") {
    x <- c(
      topo = get_topo_global(path),
      landcover = get_landcover_global(path),
      gdd = get_gdd_global(path),
      pop = get_pop_global(path),
      roads = get_network_distance(path,
        network_type = "roads",
        domain = domain
      ),
      rails = get_network_distance(path,
        network_type = "rails",
        domain = domain
      ),
      soilvars = get_soilvars_global(path),
      biovars_ = get_biovars_global(path),
      gdd = get_gdd_global(path)
    )

    names(x)[grep("biovars_1$", names(x))] <- "biovars_01"
  } else if (domain == "USA") {
    x <- c(
      topo = get_topo_conus(path, key = key),
      landcover = get_landcover_conus(path),
      biovars_400m_ = get_biovars_conus(path),
      biovars_ds_ = get_biovars_conus_ds(path),
      pop = get_pop_conus(path),
      biotic_evi = get_landsat_evi_conus(path),
      biotic_ndvi = get_landsat_ndvi_conus(path),
      par = get_par_conus(path),
      soilvars = get_soilvars_conus(path),
      roads = get_network_distance(path,
        network_type = "roads",
        domain = domain
      ),
      rails = get_network_distance(path,
        network_type = "rails",
        domain = domain
      ),
      gdd_ds = get_gdd_conus_ds(path),
      gdd_1km = get_gdd_global(path)
    )
    names(x)[grep("biovars_ds_1$", names(x))] <- "biovars_ds_01"
    names(x)[grep("biovars_400m_1$", names(x))] <- "biovars_400m_01"
  } else {
    stop("Invalid domain")
  }

  x
}

#' Batch match all predictors to base rasters
#'
#' @param path Path to the data directory.
#' @param key OpenTopography API key needed only if CONUS elevation is absent.
#' @return Runs raster matching for global and CONUS predictors.
#' @export
batch_match_to_base <- function(path, key = NULL) {
  bases <- list(
    Global = base_global(path),
    USA = base_conus(path, key = key)
  )

  lapply(
    get_rasters(path, "Global"),
    match_to_base,
    base = bases$Global,
    domain = "Global",
    path = path,
    crop = FALSE
  )

  lapply(
    get_rasters(path, "USA", key = key),
    match_to_base,
    base = bases$USA,
    domain = "USA",
    path = path,
    crop = FALSE
  )
}

#' Get resampled rasters
#'
#' @param path Path to the data directory.
#' @param domain Either `"Global"` or `"USA"`.
#' @param res Numeric vector of resolutions.
#' @param key OpenTopography API key needed only if CONUS elevation is absent.
#' @return A nested list of resampled rasters.
#' @export
get_resampled_rasters <- function(path, domain, res, key = NULL) {
  if (!domain %in% c("Global", "USA")) stop("Invalid domain")

  files <- lapply(
    get_rasters(path, domain, key = key),
    get_filename,
    domain = domain,
    path = path
  )

  lapply(files, function(x) {
    lapply(res, function(y) terra::rast(x[grep(paste0(y, "m"), x)]))
  })
}

#' Expand grouped raster variable names
#'
#' @param downscaled Logical; whether to use downscaled climate rasters.
#' @param domain Either `"Global"` or `"USA"`.
#' @param var_names Character vector of variable names.
#' @return Expanded character vector of raster names.
#' @export
fix_raster_names <- function(downscaled, domain, var_names) {
  add_bio <- function(v, idx) {
    var_names <<- var_names[var_names != v]

    prefix <- if (domain != "USA") {
      "biovars_"
    } else if (isTRUE(downscaled)) {
      "biovars_ds_"
    } else {
      "biovars_400m_"
    }

    if (domain == "USA") {
      print(paste(
        "For",
        v,
        "selecting",
        ifelse(isTRUE(downscaled), "downscaled", "resampled"),
        "rasters"
      ))
    }

    var_names <<- c(var_names, paste0(prefix, idx))
  }

  if ("precip" %in% var_names) add_bio("precip", 13:19)
  if ("temp" %in% var_names) add_bio("temp", c("01", 2:12))

  var_names
}

#' Subset predictor rasters
#'
#' @param path Path to the data directory.
#' @param domain Either `"Global"` or `"USA"`.
#' @param all Logical; whether to return all predictors.
#' @param temp Logical; include temperature predictors.
#' @param precip Logical; include precipitation predictors.
#' @param topo Logical; include topographic predictors.
#' @param land Logical; include landcover predictors.
#' @param soils Logical; include soil predictors.
#' @param pop Logical; include population predictors.
#' @param gdd Logical; include growing degree day predictors.
#' @param biotic Logical; include NDVI and EVI predictors. USA only.
#' @param light Logical; include PAR predictors. USA only.
#' @param rails Logical; include rail-distance predictors.
#' @param roads Logical; include road-distance predictors.
#' @param downscaled Logical; whether to use downscaled climate rasters for USA.
#' @param key OpenTopography API key needed only if CONUS elevation is absent.
#' @return A named list of selected predictor rasters.
#' @export
subset_rasters <- function(path, domain, all = TRUE, temp = TRUE, precip = TRUE,
                           topo = TRUE, land = TRUE, soils = TRUE, pop = TRUE,
                           gdd = TRUE, biotic = TRUE, light = TRUE,
                           rails = TRUE, roads = TRUE, downscaled = TRUE,
                           key = NULL) {
  rasters <- get_rasters(
    path,
    domain,
    key = key
  )

  if (!isTRUE(all)) {
    vars <- if (domain == "Global") {
      c(
        temp = temp,
        precip = precip,
        topo = topo,
        landcover = land,
        soilvars = soils,
        pop = pop,
        gdd = gdd,
        rails = rails,
        roads = roads
      )
    } else if (domain == "USA") {
      c(
        temp = temp,
        precip = precip,
        topo = topo,
        landcover = land,
        soilvars = soils,
        pop = pop,
        gdd = gdd,
        biotic = biotic,
        par = light,
        rails = rails,
        roads = roads
      )
    } else {
      stop("Invalid domain")
    }

    var_names <- names(vars)[vars]

    print(paste("Selecting predictors:", paste(var_names, collapse = ", ")))

    var_names <- fix_raster_names(downscaled, domain, var_names)
    rasters <- rasters[grep(paste(var_names, collapse = "|"), names(rasters))]
  } else {
    print("Selecting all predictors")

    rasters <- rasters[
      grep(ifelse(isTRUE(downscaled), "_400m_", "_ds_"),
        names(rasters),
        invert = TRUE
      )
    ]

    print(paste(
      "For temp & precip, selecting",
      ifelse(isTRUE(downscaled), "resampled", "downscaled"),
      "rasters"
    ))
  }

  rasters
}
