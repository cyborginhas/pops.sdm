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

#' @description Function to pull out the predictors that are categorical
#' and create a list of the integer and continuous predictors.
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

#' @description Function to partition the data using  environmenal
#' and spatial cross-validation, and the inputs from geo_filter_occs.
#' The function returns the partitionned data.
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
    n_part = 5,
    min_res_mult = 3,
    max_res_mult = 300,
    num_grids = 45,
    prop = 0.75,
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
