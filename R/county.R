#require(rnaturalearth)
#' @export
#' @param drive_path Letter of drive that points to existing data
#' @param names the state of interest
#' @param cty_names the county or counties of interest

county <- function(names, cty_names, drive_path){
  usa <- terra::vect(paste0(drive_path,'/Data/Vector/USA/us_lower_48_counties.gpkg'))
  ste <- usa[usa$STATE_NAME==names,]
  cty <- ste[ste$NAME%in%cty_names,]
  return(cty)
}
