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
  if (country == "world" && res < 1000) {
    print("Resolution must be >= 1000m for global analysis. Returning 1000m base
    raster.")
    res <- 1000
  } else if (country != "world" && res < 30) {
    print("Resolution must be >= 30m for country analysis. Returning 30m base 
    raster.")
    res <- 30
  } else {
    res <- res
  }
  return(res)
}

#' @description This function checks if the country provided by the user is a
#' valid ISO-3 code. If the country is not a valid ISO-3 code, the function
#' returns an error.
#' @param country The ISO-3 code for the country to create the base raster for;
#' if left NULL, it is assumed that the analysis is global.
#' @export

check_country <- function(country) {
  if (nchar(country) > 3) {
    stop("Country is not a valid ISO-3 code.")
  } else if (is.null(country)) {
    print("No country ISO-3 provided. Assuming global analysis.")
  } else {
    country <- country
  }
  return(country)
}

