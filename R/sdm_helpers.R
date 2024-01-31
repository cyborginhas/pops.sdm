#' @export

load_raster <- function(drive_path, fact, res) {
  base.rast <- terra::rast(drive_path)
  base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
  base.rast <- terra::aggregate(x = base.rast, fact = fact)
  base.rast <- base.rast / base.rast
  terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
}
