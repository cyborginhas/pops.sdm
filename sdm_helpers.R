#' Description: Checks and helpers for run_sdm.R

#' @description Function to check if the domain provided by the user is in the
#' correct format. If the domain is not in the correct format, the function
#' returns an error message.
#' @param domain The domain provided by the user. Must be "USA" or "Global".

check_domain <- function(domain) {
  if (domain != "USA" && domain != "Global") {
    stop("Domain must be 'USA' or 'Global'")
  }
}

#' @description Function to update the resolution based on the user input,
#' and the domain provided by the user. If the domain is "USA", the
#' resolution is set to can only be 30, 100, 250, 500, 1000; if the
#' domain is "Global", the resolution can only be 1000, 2500, 5000.
#' @param res The resolution provided by the user.
#' @param domain The domain provided by the user.

fix_resolution <- function(res, domain) {
  if (domain == "USA") {
    common_res_names <- c(30, 100, 250, 500, 1000)
    res <- common_res_names[which.min(abs(common_res_names - res))]
    print(paste0("Resolution set to: ", res, "m"))
  } else if (domain == "Global") {
    common_res_names <- c(1000, 2500, 5000)
    res <- common_res_names[which.min(abs(common_res_names - res))]
    print(paste0("Resolution set to: ", res, "m"))
  }
  return(res)
}


#' @description This function checks if the species name provided by the user is
#' in the correct format. If the species name is not in the correct format, the
#' function returns the species name in the correct format.
#' @param species The species name provided by the user.

format_species_name <- function(species) {
  # Replace " " with "_" in species name
  species <- gsub(" ", "_", species)
  # ensure the first letter is in uppercase and the rest are in lowercase
  species <- tolower(species)
  # make first letter uppercase
  species <- paste0(
    toupper(substr(species, 1, 1)),
    substr(species, 2, nchar(species))
  )
  return(species)
}

#' @description This function reads in predictors and crops them to the domain
#' of interest. The function returns a list of cropped predictors.
#' @param file The predictor file to be cropped.
#' @param pred_copies_path The directory path to copy the predictor files
#' will be copied to.
#' @param pred_cropped_path The directory path to the cropped predictor
#' files will be saved.
#' @param extent The base raster cropped to the study extent.

crop_predictor <- function(
    file,
    pred_copies_path,
    pred_cropped_path,
    extent) {
  # Check if the file exists, if it does not, copy, crop & zscore
  cropped_file <- paste0(pred_cropped_path, basename(file))
  cropped_file <- gsub(".tif", "_cropped.tif", cropped_file)
  zscore_file <- gsub("/cropped/", "/cropped/transformed/", cropped_file)
  zscore_file <- gsub(".tif", "_zscore.tif", zscore_file)

  if (!file.exists(cropped_file)) {
    print(paste0("Cropping ", basename(file)))
    s <- Sys.time()
    # Copy the file to the predictors directory
    file.copy(file, pred_copies_path)
    copy <- paste0(pred_copies_path, basename(file))
    r <- terra::rast(copy)
    # Replace the path with the new path
    new_filename <- paste0(pred_cropped_path, basename(copy))
    # Add cropped to the filename
    new_filename <- gsub(".tif", "_cropped.tif", new_filename)
    dtype <- terra::datatype(r)
    # Crop the raster to the study extent
    terra::crop(r, extent,
      filename = new_filename,
      wopt = list(
        gdal = "COMPRESS=ZSTD",
        datatype = dtype, overwrite = TRUE
      )
    )
    cropped_raster <- terra::rast(new_filename)
    # Throw an error if copy does not contain the correct path
    if (!grepl("1_Inputs/2_Predictors/1_Current/copies", copy)) {
      stop("Attemping to delete a file that is not in the correct directory.")
    }
    # Change file access permissions
    Sys.chmod(copy, mode = "0777")
    # Remove the file
    file.remove(copy)
    t <- Sys.time() - s
    print(paste0(
      "Time taken to crop ", basename(new_filename), ": ",
      t, attr(t, "units")
    ))
    if (dtype != "INT1U") {
      zscore_raster <- spatialEco::raster.transformation(
        cropped_raster,
        trans = "std"
      )
      terra::writeRaster(zscore_raster, zscore_file,
        overwrite = TRUE, gdal = "COMPRESS=ZSTD", datatype = "FLT4S"
      )
    } else {
      zscore_raster <- cropped_raster
    }
    terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  } else {
    cropped_raster <- terra::rast(cropped_file)
    dtype <- terra::datatype(cropped_raster)
    if (dtype != "INT1U") {
      zscore_raster <- terra::rast(zscore_file)
    } else {
      zscore_raster <- cropped_raster
    }
  }
  return(zscore_raster)
}

#' @description Function to crop species occurrence data to a specified extent
#' and year.
#' @param occs A data table with the species occurrence data retrieved from
#' batch_get_pts.
#' @param path The path to the species occurrence data. Should point to the main Data folder.
#' @param extent The cropped base map raster.
#' @param year The year of interest. If the year is not provided, the function
#' will crop the species occurrence data to the extent of interest.
#' @return A cropped data table with the species occurrence data.
#' @export crop_species_occurrences

prep_occurrences <- function(occs, extent, year = NULL) {
  occs$lat <- as.numeric(occs$lat)
  occs$lon <- as.numeric(occs$lon)
  occs_pts <- terra::vect(occs,
    crs = terra::crs(extent),
    geom = c("lon", "lat")
  )
  occs_extent <- terra::crop(occs_pts, extent)
  occs_extent$lon <- terra::crds(occs_extent)[, 1]
  occs_extent$lat <- terra::crds(occs_extent)[, 2]
  occs_extent <- data.table::as.data.table(occs_extent)

  # Filter by year if provided
  if (!is.null(year)) {
    occs_pa$date <- as.numeric(format(as.Date(occs_pa$date,
      format = "%Y"
    ), "%Y"))
    occs_extent <- occs_extent[date >= year]
  }
  return(occs_extent)
}

#' @description Function to retrieve base raster for domain and resolution
#' and extent of interest.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster to create.
#' @param path Path for the Data folder
#' @param extent Path to vector file of extent of interest to crop the base raster.
#' The cropped based raster will be saved locally in specified output path.
#' @note Will add option to specify resolutions beyond 30m, 100m, 250m, 500m, 1000m
#' in future versions.

crop_base_raster <- function(domain, res, path, extent) {
  cropped_base_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/3_Calibration_area/base_study_extent.tif") # no lint
  if (!file.exists(cropped_base_path)) {
    # Retrieve base raster
    if (domain == "USA") {
      common_res_names <- c(30, 100, 250, 500, 1000)
      base <- base_conus(path)
      # Pull the base raster at the specified resolution
      base <- base[[which(common_res_names == res)]]
    } else if (domain == "Global") {
      common_res_names <- c(1000, 2500, 5000)
      base <- base_global(path)
      # Pull the base raster at the specified resolution
      base <- base[[which(common_res_names == res)]]
    } else {
      stop("Domain must be 'USA' or 'Global'")
    }
    # Crop the base raster to the extent of interest
    extent <- terra::vect(extent)
    extent <- terra::project(extent, y = terra::crs(base))
    base <- terra::crop(base, extent)
    terra::writeRaster(base, cropped_base_path,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD", datatype = "INT1U"
    )
  } else {
    base <- terra::rast(cropped_base_path)
  }
  return(base)
}

#' @description Function to export the cropped species occurrence
#' to the local output path. With option to filter NEON to use
#' as independent validation data. Data filtered so that no
#' more than one occurrence per lat/lon pair is included.
#' @param occs The cropped species occurrence data.

export_occurrences <- function(occs) {
  sciname <- unique(occs$sciname)
  sciname <- format_species_name(sciname)
  neon <- occs[db == "neon", ]
  occs <- occs[db != "neon", ]
  neon <- unique(neon[, .(x = lon, y = lat, pr_ab = as.integer(p_a))])
  neon$id <- 1:nrow(neon)
  occs <- unique(occs[, .(x = lon, y = lat, pr_ab = as.integer(p_a))])
  occs$id <- 1:nrow(occs)
  # Create outpath
  occs_dirpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/original/"
  )
  dir.create(occs_dirpath, showWarnings = FALSE, recursive = TRUE)
  write.csv(occs,
    file = paste0(occs_dirpath, sciname, "_train_occs.csv"),
    row.names = FALSE
  )
  write.csv(neon,
    file = paste0(occs_dirpath, sciname, "_test_occs.csv"),
    row.names = FALSE
  )
  return(occs)
}

#' @description Function to help determine the range of values to apply for
#' geographic filtering. The function returns the range of values to apply
#' for geographic filtering. with the flexsdm::occfilt_geo function.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' If the domain is "USA", the range of values is set to 90m intervals;
#' if the domain is "Global", the range of values is set to 150m intervals.
#' @param res The resolution of the analysis.

get_geo_cellsizes <- function(domain, res) {
  if (domain == "USA") {
    range <- seq(res, 1000, by = 3 * res)
  } else {
    range <- seq(res, 5000, by = 150)
  }
  return(range)
}

#' @description Function to apply geographic filtering to the species occurrence
#' data. The function returns the filtered species occurrence data.
#' @param occs The species occurrence data.
#' @param extent The base raster cropped to the extent of interest.
#' @param d The distance to apply for geographic filtering.
#' @param species The species name.

geo_filter_occs <- function(occs, extent, d, species) {
  filt_geo <- flexsdm::occfilt_geo(
    data = occs,
    x = "x",
    y = "y",
    env_layer = extent,
    method = c("defined", d = d / 1000)
  )
  filt_geo$d <- d
  # Write out the filtered data
  occsfilt_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/filtered/"
  )
  dir.create(occsfilt_path, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(filt_geo, paste0(
    occsfilt_path,
    species, "_occ_pa_geofilter_", d, "m2.csv"
  ),
  row.names = FALSE
  )
  return(filt_geo)
}

#' @description Function to remove the predictors that are categorical
#' and retain a list of the integer and continuous predictors.
#' @param predictor A predictor to be used in the analysis.

get_numerical_rasters <- function(predictor) {
  dt <- terra::datatype(predictor)
  if (dt != "INT1U") {
    predictor <- predictor
  } else {
    predictor <- NULL
  }
  return(predictor)
}


#' @description Function to keep the predictors that are categorical
#' @param predictor A predictor to be used in the analysis.

get_categorical_rasters <- function(predictor) {
  dt <- terra::datatype(predictor)
  if (dt == "INT1U") {
    predictor <- predictor
  } else {
    predictor <- NULL
  }
  return(predictor)
}

#' @description Function to partition the data using  environmenal
#' and spatial cross-validation, and the inputs from geo_filter_occs.
#' The function returns the partitioned data.
#' @param occs The filtered species occurrence data.
#' @param env_layer The list of predictors.

spatial_block_partition <- function(occs, extent, env_layer, species) {
  s <- Sys.time()
  # Partition the data
  part_data <- flexsdm::part_sblock(
    env_layer = env_layer,
    data = occs,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    n_part = 4,
    min_res_mult = 3,
    max_res_mult = 200,
    num_grids = 30,
    prop = 0.5,
    min_occ = 10
  )

  # Write out the partitioned data
  # Create outpath
  part_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/partitioned/"
  )
  dir.create(part_path, showWarnings = FALSE, recursive = TRUE)

  # Write out the partitioned data
  names(part_data) <- paste0(names(part_data), "_filtgeo_", d, "m2")
  d <- unique(occs$d)
  part_data[[2]]$d <- d
  # Best partition
  data.table::fwrite(part_data[[2]], paste0(
    part_path, species, "_best_part_filtgeo", d, "m2.csv"
  ))
  # Partition data
  part_data[[1]]$d <- d
  data.table::fwrite(part_data[[1]], paste0(
    part_path, species, "_part_data_filtgeo", d, "m2.csv"
  ))
  names(part_data[[3]]) <- paste0(names(part_data[[3]]), "_filtgeo_", d, "m2")
  # Partition grid
  terra::writeRaster(part_data[[3]], paste0(
    part_path, species, "_part_grid_filtgeo", d, "m2.tif"
  ), overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=ZSTD")

  t <- Sys.time() - s
  print(paste0(
    "Time taken to partition ", length(occs$x),
    " points: ", t, attr(t, "units")
  ))
  return(part_data)
}

#' @description A function to create randomly sampled background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param species The species name.
#' associated grid.

create_random_bg_pts <- function(partitioned_data, species) {
  occs_p <- partitioned_data[[1]]
  grid_env <- partitioned_data[[3]]
  d <- unique(occs_p$d)

  bg_pts <- flexsdm::sample_background(
    data = occs_p,
    x = "x",
    y = "y",
    n = length(occs_p$x),
    method = "random",
    rlayer = grid_env
  )

  # Extract grid values
  part <- terra::extract(grid_env, bg_pts[, c("x", "y")], ID = FALSE)
  part <- part[[1]]
  bg_pts$.part <- part

  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/random/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  return(occs_p)
}


#' @description A function to create thickened background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param base The base raster cropped to the extent of interest.
#' @param res The resolution of the base raster.
#' @param species The species name.
#' associated grid.

create_thickened_bg_pts <- function(
    partitioned_data,
    base, res, species) {
  s <- Sys.time()
  occs_p <- partitioned_data[[1]]
  grid_env <- partitioned_data[[2]]

  # Disaggregate the spatial block grid

  if (res == 30) {
    factor <- floor(terra::res(grid_env)[1] / terra::res(base)[1] / 5)
    grid_env <- terra::disagg(grid_env, factor)
  } else {
    grid_env <- terra::resample(grid_env, base, method = "near")
  }

  d <- unique(occs_p$d)
  no_parts <- length(unique(occs_p$.part))

  part_pts <- list()
  bg_pts <- list()
  grid_parts <- list()

  for (i in 1:no_parts) {
    s <- Sys.time()
    part_pts[[i]] <- occs_p[occs_p$.part == i, ]
    grid_parts[[i]] <- terra::ifel(grid_env == i, 1, NA)
    bg_pts[[i]] <- flexsdm::sample_background(
      data = part_pts[[i]],
      x = "x",
      y = "y",
      n = length(part_pts[[i]]$x),
      method = c("thickening", width = 10000),
      rlayer = grid_parts[[i]],
      maskval = 1
    )
    bg_pts[[i]]$.part <- i
    t <- Sys.time() - s
    print(paste0(
      "Time taken to create background points for part ", i, ": ",
      t, " ", attr(t, "units")
    ))
  }

  bg_pts <- data.table::rbindlist(bg_pts, use.names = TRUE, fill = TRUE)
  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/thicken/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  t <- Sys.time() - s
  print(paste0(
    "Time taken to create background points: ", t, " ", attr(t, "units"),
    " for ", length(occs_p$x), " points for filtgeo", d, "m2"
  ))
  return(occs_p)
}

#' @description Function to create a target group bias file from
#' iNaturalist tree occurrence data. The function returns the target
#' group bias file to be used with the sample_background function.
#' @param path A path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster.

target_group_rbias_lyr <- function(path, extent, domain, res) {
  # Copy the tree data to the target group directory
  path2 <- gsub(
    "/.*",
    "/auto_arborist_cvpr2022_v0.15/data/tree_classification/inat/metadata/", # nolint
    path
  )

  tg_outpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/target_group/"
  ) # nolint

  dir.create(tg_outpath, showWarnings = FALSE, recursive = TRUE)

  fn <- paste0(tg_outpath, "TargetGroup_biasfile.tif")
  if (!file.exists(fn)) {
    if (domain == "USA") {
      file.copy(
        paste0(path2, "tree_observations_uscanada.csv"),
        paste0(tg_outpath, "tree_observations_uscanada.csv")
      )
    } else {
      file.copy(
        paste0(path2, "tree_observations_global.csv"),
        paste0(tg_outpath, "tree_observations_global.csv")
      )
    }

    # Read in the tree data
    trees <- data.table::fread(paste0(
      tg_outpath,
      "tree_observations_uscanada.csv"
    ))

    base_extent <- terra::ext(extent)

    # Crop the tree data to the extent of interest
    trees <- trees[longitude >= base_extent[1] & longitude <= base_extent[2] &
      latitude >= base_extent[3] & latitude <= base_extent[4], ]

    trees <- unique(trees[, .(longitude, latitude,
      sciname = paste0(genus, " ", species)
    )])

    trees <- unique(trees[, .(count = .N), by = .(longitude, latitude)])

    # Write out the tree data
    data.table::fwrite(trees, paste0(
      tg_outpath,
      "tree_observations_cropped.csv"
    ),
    row.names = FALSE
    )

    # Delete the original tree data
    file.remove(paste0(tg_outpath, "tree_observations_uscanada.csv"))

    # Calculate scale for the tree data
    scale <- 1 / sum(trees$count)
    trees$scale <- trees$count * scale
    trees_crds <- as.matrix(trees[, .(longitude, latitude)])

    if (res == 30) {
      base_extent <- terra::aggregate(extent, factor = 2)
    } else {
      base_extent <- extent
    }

    # Do a 2d kernel density estimation
    target_density <- ks::kde(trees_crds,
      w = trees$scale,
      gridsize = c(nrow(base_extent), ncol(base_extent))
    )

    target_raster <- base_extent
    terra::values(target_raster) <- target_density$estimate

    # Normalize the target bias file between 0 and 1
    min <- terra::global(target_raster, "min")
    min <- min$min

    f <- function(i) i - min
    target_raster <- terra::app(target_raster, f)
    target_raster <- spatialEco::raster.transformation(target_raster,
      trans = "norm"
    )

    if (res == 30) {
      target_raster <- terra::resample(target_raster, extent)
    } else {
      target_raster <- target_raster
    }

    terra::writeRaster(target_raster, paste0(
      tg_outpath,
      "TargetGroup_biasfile.tif"
    ),
    overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD"
    )
  } else {
    target_raster <- terra::rast(fn)
  }
  return(target_raster)
}

#' @description Function to create a population density bias file from
#' the Population Density raster. The function returns the population
#' density bias file to be used with the sample_background function.
#' @param path A path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param rbias_type The type of bias file to use. Options are
#' "pop_density", "distance_roads", or "distance_rails".
#' @param domain The domain of interest. Options are "USA" or "Global".
#' @param res The resolution of the base raster.

human_factors_rbias_lyr <- function(path, extent, domain, rbias_type, res) {
  # Create the output path
  hf_outpath <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/", rbias_type, "/"
  )
  dir.create(hf_outpath, showWarnings = FALSE, recursive = TRUE)
  fn <- paste0(hf_outpath, paste0(rbias_type, "_biasfile.tif"))

  if (!file.exists(fn)) {
    # Pull in the population density data from the Data folder
    if (rbias_type == "pop_density") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = TRUE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = FALSE, roads = FALSE,
        downscaled = FALSE
      )
    } else if (rbias_type == "distance_roads") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = FALSE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = FALSE, roads = TRUE,
        downscaled = FALSE
      )
    } else if (rbias_type == "distance_rails") {
      hf <- subset_rasters(path, domain,
        all = FALSE, temp = FALSE,
        precip = FALSE, topo = FALSE, land = FALSE, soils = FALSE, pop = FALSE,
        gdd = FALSE, biotic = FALSE, light = FALSE, rails = TRUE, roads = FALSE,
        downscaled = FALSE
      )
    }
    hf_files <- as.vector(unlist(lapply(hf, function(x) {
      get_filename(x, "USA", path)
    })))
    subset_files <- hf_files[grep(res, hf_files)]
    hf <- terra::rast(subset_files)
    # Crop the population density to the extent of interest and write it out
    hf <- terra::crop(hf[[1]], extent, filename = paste0(
      hf_outpath, rbias_type, "_original_cropped.tif"
    ), wopt = list(
      gdal = "COMPRESS=ZSTD", datatype = "FLT4S",
      overwrite = TRUE
    ))
    # plot all values <1000
    hf2 <- terra::ifel(hf[[1]] < 1000, NA, hf)

    hf <- log(hf)

    if (rbias_type == "distance_roads" | rbias_type == "distance_rails") {
      remove_infs <- function(x) {
        x[x == -Inf] <- NA
        return(x)
      }
      hf <- terra::app(hf, remove_infs)
      hf <- terra::app(hf, function(x) -x)
    } else {
      hf <- hf
    }
    # Normalize bias file to between 0 and 1.
    min <- terra::global(hf, "min", na.rm = TRUE)$min
    hf <- terra::app(hf, function(x) x - min)

    # Normalize bias file to between 0 and 1.
    hf <- spatialEco::raster.transformation(hf, trans = "norm")

    # Add 0.1 to the bias file to avoid 0 values
    hf <- hf + 0.1
    # Convert NA values to 0
    hf <- terra::app(hf, function(x) ifelse(is.na(x), 0, x))
    # Normalize
    hf <- spatialEco::raster.transformation(hf, trans = "norm")
    terra::writeRaster(hf, paste0(
      hf_outpath, rbias_type, "_biasfile.tif"
    ), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD")
    unlink(paste0(hf_outpath, rbias_type, "_original_cropped.tif"))
  } else {
    hf <- terra::rast(fn)
  }

  return(hf)
}


#' @description Function to create sample background points using the
#' @description A function to create thickened background points data
#' and write the data to the local output path.
#' @param partitioned_data A list of partitioned data occurrences, and the
#' associated spatial block grid.
#' @param rbias_lyr The bias file: target group raster, population raster, or
#' distance rasters.
#' @param rbias_type The type of bias file to use. Options are "target_group",
#' "pop_density", "distance_roads", or "distance_rails".
#' @param species The species name.
#' associated grid.


create_biased_bg_pts <- function(
    partitioned_data, rbias_lyr,
    rbias_type, species) {
  s <- Sys.time()
  occs_p <- partitioned_data[[1]]
  grid_env <- partitioned_data[[2]]

  # Disaggregate the spatial block grid
  d <- unique(occs_p$d)
  no_parts <- length(unique(occs_p$.part))

  part_pts <- list()
  bg_pts <- list()
  grid_parts <- list()

  for (i in 1:no_parts) {
    part_pts[[i]] <- occs_p[occs_p$.part == i, ]
    grid_parts[[i]] <- terra::ifel(grid_env == i, 1, NA)
    bg_pts[[i]] <- flexsdm::sample_background(
      data = part_pts[[i]],
      x = "x",
      y = "y",
      n = length(part_pts[[i]]$x),
      method = "biased",
      rlayer = grid_parts[[i]],
      maskval = 1,
      rbias = rbias_lyr
    )
    bg_pts[[i]]$.part <- i
  }

  bg_pts <- data.table::rbindlist(bg_pts, use.names = TRUE, fill = TRUE)
  # Combine the background points with the occurrences
  occs_p <- data.table::rbindlist(list(occs_p, bg_pts),
    use.names = TRUE,
    fill = TRUE
  )
  occs_p$d <- d

  # Write out the background points
  bg_path <- paste0(
    getwd(),
    "/flexsdm_results/1_Inputs/1_Occurrences/background/", rbias_type, "/"
  )
  dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(occs_p, paste0(
    bg_path, species, "_", rbias_type, "_bg_pts_filtgeo_", d, "m2.csv"
  ), row.names = FALSE)
  t <- Sys.time() - s
  print(paste0(
    "Time taken to create background points: ", t, " ", attr(t, "units"),
    " for ", length(occs_p$x), " points for filtgeo", d, "m2"
  ))
  return(occs_p)
}

#' @description Function to conduct cluster analysis on the predictors
#' based on the correlation threshold. A sample_size is used to reduce
#' the time taken to perform the cluster analysis.
#' @param data A list of  predictors that have been z-score
#' normalized.
#' @param sample_size The number of samples to use in the cluster analysis.
#' @param mincor The minimum correlation threshold to use in the cluster analysis.
#' @param categorical_vars A list of categorical variables to include in the

cluster_analysis <- function(data, sample_size, mincor, categorical_vars) {
  # Perform cluster analysis
  tryCatch(
    {
      sample <- terra::spatSample(data,
        size = sample_size, method = "regular",
        as.df = TRUE, na.rm = TRUE, values = TRUE,
        xy = TRUE
      )

      # Keep all columns except 1:2 and cat_cols
      ccres <- klaR::corclust(sample[, -c(1:2)])
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      message("Please increase the sample_size.")
    }
  )
  # Create hierarchical cluster tree based correlation threshold
  cluster_dt <- klaR::cvtree(ccres, mincor = mincor)
  vars <- rownames(cluster_dt$correlations)
  cluster_dt <- as.data.table(cluster_dt$correlations)
  cluster_dt$var <- vars

  # Reassign clusters
  cluster_dt[, clustercount := .N, by = .(cluster)]
  singles <- cluster_dt[cluster_dt$clustercount == 1 &
    cluster_dt$av.cor2closest < mincor]
  multis <- cluster_dt[cluster_dt$clustercount > 1]
  multisfix <- multis[av.cor2closest > mincor, ]
  multisfix[, cluster := pmin(cluster, closest)]
  multis <- multis[av.cor2closest < mincor, ]
  cluster_dt <- rbind(singles, multis, multisfix)
  cluster_dt <- cluster_dt[, .(var, cluster)]
  cluster_dt <- cluster_dt[order(cluster_dt$cluster, cluster_dt$var)]
  # Add categorical variables to the cluster assignments
  if (!is.null(categorical_vars)) {
    cat_vars <- names(categorical_vars)
    cat_clusters <- max(cluster_dt$cluster) + 1:length(cat_vars)
    cat_vars <- data.table(var = cat_vars, cluster = cat_clusters)
    cluster_dt <- rbind(cluster_dt, cat_vars)
  }
  return(cluster_dt)
}


#' @description Generate combinations of variables based on clusters
#' @param vars A character vector containing the variable names
#' @param clusters A character vector containing the cluster assignments
#' @return A list containing data frames of variable combinations for each
#' cluster count

library(dplyr)
library(purrr)

# Sample data frame
df <- data.frame(
  var = c("A", "B", "C", "D", "E", "F"),
  cluster = c(1, 1, 2, 2, 3, 3)
)

# Function to generate combinations
generate_combinations <- function(df) {
  # Ensure the data frame is sorted by cluster
  df <- dplyr::arrange(df, cluster)

  # Split the data frame by cluster
  split_df <- split(df$var, df$cluster)

  # Get unique clusters
  unique_clusters <- unique(df$cluster)

  # Helper function to generate all valid combinations
  generate_valid_combinations <- function(vars_by_cluster) {
    combs <- list()
    for (i in 1:length(vars_by_cluster)) {
      cluster_combinations <- combn(seq_along(vars_by_cluster), i, simplify = FALSE)
      for (cluster_comb in cluster_combinations) {
        valid_comb <- do.call(expand.grid, vars_by_cluster[cluster_comb])
        combs <- append(combs, list(valid_comb))
      }
    }
    combs
  }

  # Generate all valid combinations of variables from different clusters
  all_combinations <- generate_valid_combinations(split_df)

  # Combine all combinations into a single data frame
  final_combinations <- purrr::map_dfr(all_combinations, dplyr::as_tibble)

  return(final_combinations)
}

#' @description Function to append enviromental data to the partitioned
#' presence-background data
#' @param part_data A list of partitioned data presence-background data.
#' @param env_layer The list of predictors to extract.
#' @param bg_method The method used to sample the background points: random,
#' thicken, or biased: target_group, pop_density, distance_roads, or
#' distance_rails.
#' @param species The species name.

append_env_data <- function(part_data, env_layer, bg_method, species) {
  s <- Sys.time()
  d <- unique(part_data$d)
  # Convert to spatvect
  part_data <- terra::vect(part_data, crs = terra::crs(env_layer), geom = c("x", "y"))
  # Extract the environmental data
  part_data <- terra::extract(env_layer, part_data, ID = FALSE, xy = TRUE, bind = TRUE, method = "simple")
  part_data <- data.table::as.data.table(part_data)
  fwrite(part_data, paste0(
    getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/background/", bg_method, "/",
    species, "_", bg_method, "_bg_pts_filtgeo_", d, "m2_wpreds.csv"
  ), row.names = FALSE)
  t <- Sys.time() - s
  print(paste0(
    "Time taken to append predictor data: ", t, " ", attr(t, "units")
  ))
  return(part_data)
}

#' @description Function to pull in systematically sampled occurrence data
#' from NEON and VMI NPS Inventory data and create testing sets for
#' independent validation.
#' @param path The path to the Data folder.
#' @param extent The base raster cropped to the extent of interest.
#' @param species The species name.


create_testing_set <- function(path, extent, species) {
  test_data_path <- paste0(getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/original/")
  neon <- data.table::fread(paste0(test_data_path, species, "_test_occs.csv"))
  vmi <- data.table::fread(paste0(path, "/Original/vmi_nps_inventory/data/cleaned/csv/aial.csv"))
  vmi <- unique(vmi[, .(x = lon, y = lat, pr_ab = p_a)])
  vmi <- vmi[, .(pr_ab = max(pr_ab)), by = .(x, y)]
  test_data <- data.table::rbindlist(list(neon, vmi), fill = TRUE, use.names = TRUE)
  test_data <- test_data[, .(x, y, pr_ab)]
  test_data <- terra::vect(test_data, crs = terra::crs(extent), geom = c("x", "y"))
  test_data <- terra::crop(test_data, extent)
  test_data$x <- terra::crds(test_data)[, 1]
  test_data$y <- terra::crds(test_data)[, 2]
  pr_ab_1 <- data.table::as.data.table(test_data[test_data$pr_ab == 1])
  pr_ab_0 <- data.table::as.data.table(test_data[test_data$pr_ab == 0])
  # Randomly select length(pr_ab_1) points from pr_ab_0 ten times
  testing_sets <- list()

  for (i in 1:5) {
    pr_ab_0_sample <- pr_ab_0[sample(1:nrow(pr_ab_0), nrow(pr_ab_1)), ]
    testing_set <- data.table::rbindlist(list(pr_ab_1, pr_ab_0_sample), fill = TRUE)
    testing_sets[[i]] <- testing_set
    testing_sets[[i]]$setid <- i
  }
  testing_sets <- data.table::rbindlist(testing_sets, fill = TRUE)
  fwrite(testing_sets, paste0(
    getwd(), "/flexsdm_results/1_Inputs/1_Occurrences/original/",
    species, "_testing_sets.csv"
  ), row.names = FALSE)
  return(testing_sets)
}

#' Function to setup sdm prepped datasets for
#' analysis.
#' @param data The target group background points data
splitpts4sdm <- function(data) {
  #' Remove unnecessary columns: d, x, y using data.table
  data <- data[, c("d", "x", "y") := NULL]
  #' Rename NLCD Land Cover Class to landcoverrc
  setnames(data, "NLCD Land Cover Class", "landcoverrc")
  #' Convert to tibble
  data <- tibble::as_tibble(data)
  #' Train on .part == 1:4; test on .part == 5
  no_parts <- length(unique(data$.part))
  data <- data[data$.part %in% 1:(no_parts - 1), ]
  #' Pull out presence-absence data
  response <- data[data$pr_ab == 1, ]
  #' Add id column
  response$id <- 1:nrow(response)
  #' Pull out background data; pr_ab == 0
  bg <- data[data$pr_ab == 0, ]
  return(list(response = response, bg = bg))
}


#' @description Function to aggregate transformed predictors by a factor of 3
#' to speed up the analysis.
#' @param predictor A cropped and transformed predictor.
#' @param factor The factor to aggregate the predictor by.

aggregate_predictor <- function(predictor, factor) {
  s <- Sys.time()
  fn <- basename(terra::sources(predictor))
  fn <- gsub(".tif", "_agg.tif", fn)
  dt <- terra::datatype(predictor)
  agg_outpath <- paste0(
      getwd(),
      "/flexsdm_results/1_Inputs/2_Predictors/1_Current/cropped/transformed/aggregated/" # nolint
    )
    dir.create(agg_outpath, showWarnings = FALSE, recursive = TRUE)
  if (dt != "INT1U") {
    predictor <- terra::aggregate(predictor, fact = factor, fun = "mean")
    terra::writeRaster(predictor, paste0(
      agg_outpath, fn
    ), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD")
  } else {
    predictor <- NULL
  }
  t <- Sys.time() - s
  print(paste0(
    "Time taken to aggregate predictor: ", fn, " ", t, " ", attr(t, "units")
  ))
  return(predictor)
}
