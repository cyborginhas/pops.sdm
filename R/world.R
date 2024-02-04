#require(rnaturalearth)
#' @export
#' @param drive_path Letter of drive that points to existing data

world <- function(){
  countries <- terra::vect(paste0(drive_path, '/Data/Original/ne_10m_admin_0_countries_lakes/ne_10m_admin_0_countries_lakes.shp'))
  return(countries)
}
