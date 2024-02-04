#require(rnaturalearth)
#' @export
#' @export names a character vector of country names

country <- function(names){
  country <- rnaturalearth::ne_countries(scale=10, country=names)
  return(terra::vect(country))
}
