#' @export
#' @param drive_path Letter of drive that points to existing data
#' @param fact Factor to aggregate raster by
#' @param res Resolution of raster

scratch <- "C:/Users/blaginh/Desktop/test/"

s <- Sys.time()
geodata::soil_world(var = "bdod", depth = 5, stat = "mean", path = scratch)
e <- Sys.time()

s2 <- Sys.time()
geodata::worldclim_global(var='bio', res=.5, path=scratch)
e2 <- Sys.time()

rasterbase <- function(res){
  geodir <- paste0(drive_path,'/Data/Raster/')
  if(res==1000){
    if(file.exists(paste(geodir, 'Global/base_1000m.tif', sep=''))==F){
      base.rast <- geodata::worldclim_global(var='bio', res=.5, path=paste(geodir, 'Global/',sep=''))
      base.rast <- base.rast[[1]]
      base.rast <- base.rast/base.rast
      terra::writeRaster(base.rast, paste(geodir, 'Global/base_1000m.tif', sep=''))
    }
    else{
      base.rast <- terra::rast(paste(geodir, 'Global/base_1000m.tif', sep=''))
    }
  }

  if(res == 500){
    if(file.exists(paste(geodir, 'USA/base_', res, 'm.tif', sep=''))==F){

    }
  }


  if(res<1000 & res>33){
    if(file.exists(paste(geodir, 'USA/base_', res, 'm.tif', sep=''))==F){
      if(res!=33){
        if(res==500){
          base.rast <- terra::rast("Q:\\Shared drives\\Data\\Raster\\USA\\soils\\250m\\ph.mean.tif")
          base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
          base.rast <- terra::aggregate(x=base.rast, fact=2)
          base.rast <- base.rast/base.rast
          terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
        }
        if(res==250){
          base.rast <- terra::rast("Q:\\Shared drives\\Data\\Raster\\USA\\soils\\250m\\ph.mean.tif")
          base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
          base.rast <- base.rast/base.rast
          terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
        }
        if(res==100){
          if(file.exists('Q:\\Shared drives\\Data\\Raster\\USA\\base_33m.tif')){
            base.rast <- terra::rast('Q:\\Shared drives\\Data\\Raster\\USA\\base_33m.tif')
            base.rast <- aggregate(base.rast, 3)
            terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''), overwrite=T)
          }
          if(!file.exists('Q:\\Shared drives\\Data\\Raster\\USA\\base_33m.tif')){
          base.rast <- terra::rast(paste(geodir, 'Global\\pop\\GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif', sep=''))
          base.rast <- terra::crop(base.rast, terra::ext(pops.sdm::l48()))
          terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
          }
        }

      }
      if(res==33){
        if(file.exists(paste('Q:\\Shared drives\\Data\\Raster\\USA\\landcover\\nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))){
         base.rast <- terra::rast(paste('Q:\\Shared drives\\Data\\Raster\\USA\\landcover\\nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))
         base.rast <- terra::app(base.rast!=c(0,11), min)
         terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''), overwrite=T)
        }
        if(!file.exists(paste('Q:\\Shared drives\\Data\\Raster\\USA\\landcover\\nlcd_2019_land_cover_l48_20210604_1s.tif', sep=''))){
        base.rast <- pops.sdm::rasterbase(100)
        base.rast <- terra::disagg(x=base.rast, fact=3)
        terra::writeRaster(base.rast, paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
        }
      }
    }
    if(file.exists(paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))==T){
      base.rast <- terra::rast(paste(geodir, 'USA\\base_', res, 'm.tif', sep=''))
    }
  }
  return(base.rast)
}
