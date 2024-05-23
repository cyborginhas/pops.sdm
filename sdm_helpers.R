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
  species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
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
  terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)
  return(cropped_raster)
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

cropped_base_raster <- function(domain, res, path, extent) {
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
    method = c("defined", d = d/1000)
  )

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

#' @description Function to partition the data using spatial
#' block cross-validation, and the inputs from geo_filter_occs.
#' The function returns the partitioned data, and a raster of
#' the spatial blocks.
#' @param occs The filtered species occurrence data.
#' @param extent The base raster cropped to the extent of interest.
#' @param npart The number of partitions to create.
#' @param res The resolution of the analysis.

## npart <- 5
## extent <- base
## occs <- filt_geo[[1]]
## occs$id <- 1:nrow(occs)


## partition_data <- function(occs, npart, res) {
##   sp_part3 <- part_sblock(
##     env_layer = extent,
##     data = occs,
##     x = "x",
##     y = "y",
##     pr_ab = "pr_ab",
##     min_res_mult = 10, # Minimum value used for multiplying raster resolution and define the finest resolution to be tested
##     max_res_mult = 100, # Maximum value used for multiplying raster resolution and define the coarsest resolution to be tested
##     num_grids = 5, # Number of grid to be tested between min_res_mult X (raster resolution) and max_res_mult X (raster resolution)
##     n_part = npart, # Number of partitions
##     prop = 0.5, # Proportion of points used for testing autocorrelation between groups (0-1)
##     min_occ = floor(length(occs) / (npart * 2)) # Minimum number of occurrences to be used in each partition
##   )

##   grid_env <- get_block(env_layer = extent, best_grid = sp_part3$grid)
##   return(list(sp_part3, grid_env))
## }

## sp_part3 <- part_sblock(
##   env_layer = base_raster,
##   data = filt_geo[[1]],
##   x = "x",
##   y = "y",
##   pr_ab = "p_a",
##   min_res_mult = 10, # Minimum value used for multiplying raster resolution and define the finest resolution to be tested
##   max_res_mult = 1000, # Maximum value used for multiplying raster resolution and define the coarsest resolution to be tested
##   num_grids = 50, # Number of grid to be tested between min_res_mult X (raster resolution) and max_res_mult X (raster resolution)
##   n_part = n_parts, # Number of partitions
##   prop = 0.5, # Proportion of points used for testing autocorrelation between groups (0-1)
##   min_occ = floor(length(filt_geo$fkey)/(no_parts*2)) # Minimum number of occurrences to be used in each partition
## )

## grid_env <- get_block(env_layer = base_raster, best_grid = sp_part3$grid)