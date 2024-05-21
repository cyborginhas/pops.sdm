#' Description: Functions for checking inputs to sdm_prep functions

#' @description This function checks if the resolution provided by the user is
#' appropriate for the domain. If the resolution is not appropriate, the
#' function returns the finest resolution base raster for the domain.
#' @param country The ISO-3 code for the country to create the base raster for;
#' if left NULL, it is assumed that the analysis is global.
#' @param res A numeric value specifying the resolution of the base raster to
#' create.
#' @export

check_resolution <- function(country, res) {
  if (domain == "global" && res < 1000) {
    print("Resolution must be >= 1000m for global analysis. Returning 1000m base
    raster.")
    res <- 1000
  } else if (domain != "conus" && res < 30) {
    print("Resolution must be >= 30m for L48 analysis. Returning 30m base 
    raster.")
    res <- 30
  } else {
    res <- res
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

#' @description This function uses the species name provided by format_species_name
#' and the path to the species occurrence data to read in the species occurrence
#' data. The function reads in the species occurrence data and returns a data
#' table with the species occurrence data.
#' @param species The species name returned by format_species_name.
#' @param path The path to the species occurrence data.


get_species_occurrences <- function(species, path) {
  species <- format_species_name(species)
  files <- list.files(paste0(path, "/Table/"), pattern = species, full.names = TRUE, recursive = TRUE)
  occs <- lapply(files, function(x) {
    fread(x, colClasses = "character")
  })
  occs <- rbindlist(occs, fill = TRUE)
  occs[, `:=`(lat = as.numeric(lat), lon = as.numeric(lon), date = as.numeric(date), p_a = as.numeric(p_a))]
  occs <- occs[date > 1990]
  #occs <- occs[p_a > 0]
  return(occs)
}

#' @description After reading in the species occurrence data, this function
#' crops the species occurrence data to the domain of interest. The function
#' returns the cropped species occurrence data.
#' @param species The species name returned by format_species_name.
#' @param path The path to the species occurrence data. Should point to the
#' main Data folder.
#' @param extent The SpatVector object representing the domain of
#' interest.

crop_species_occurrences <- function(species, path, extent) {
  occs <- get_species_occurrences(species, path)
  occs_pts <- vect(occs, crs = crs(extent), geom = c("lon", "lat"))
  extent <- vect(extent)
  extent <- project(extent, crs(occs_pts))
  occs_extent <- crop(occs_pts, extent)
  occs_extent$lon <- crds(occs_extent)[, 1]
  occs_extent$lat <- crds(occs_extent)[, 2]
  occs_extent <- as.data.table(occs_extent)
}

#' @description This function
