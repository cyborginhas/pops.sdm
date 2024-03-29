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
