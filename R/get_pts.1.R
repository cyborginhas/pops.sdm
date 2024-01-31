# require(BIEN)
# require(rinat)
# require(terra)
# require(geodata)
# require(sp)
#' @export
#' @param drive_path Letter of drive that points to existing data

get_pts.1 <- function(spname, domain=NULL, drive_path){

  if(is.null(domain)){x0 <- terra::ext(c(-180, 180, -90, 90)); bb <- c(x0[3], x0[1], x0[4], x0[2])}
  if(!is.null(domain)){x1 <- terra::ext(domain); bb <- c(x1[3], x1[1], x1[4], x1[2])}

  sp.inat <- rinat::get_inat_obs(taxon_name=spname, maxresults=9999, quality='research', geo=T, bounds=bb)

  if(grepl(' ', spname)){
    sp.bien <- BIEN::BIEN_occurrence_species(species=spname, natives.only = F)
    # sp.gbif <- geodata::sp_occurrence(genus=strsplit(spname, split=' ')[[1]][1],
    #                                   species=strsplit(spname, split=' ')[[1]][2],
    #                                   ext=bb, geo=T, removeZeros=T, download=T, fixnames=T)
  }

  if(!grepl(' ', spname)){
    sp.bien <- BIEN::BIEN_occurrence_genus(genus=spname, natives.only = F)
    #sp.gbif <- geodata::sp_occurrence(genus=spname, species='*', ext=bb, geo=T, removeZeros=T, download=T, fixnames=T)
  }

  sp.pts <- data.frame('longitude'=as.numeric(c(sp.inat$longitude, sp.bien$longitude)),
                       'latitude'=as.numeric(c(sp.inat$latitude, sp.bien$latitude)))
  sp.pts <- unique(sp.pts[which(unique(!is.na(sp.pts$longitude), !is.na(sp.pts$latitude))),])
  sp.pts <- terra::vect(sp::SpatialPoints(sp.pts, sp::CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')))

  if(is.null(domain)){sp.pts <- terra::crop(x=sp.pts, y=pops.sdm::world(drive_path))}
  if(!is.null(domain)){sp.pts <- terra::crop(x=sp.pts, y=domain)}

  return(sp.pts)
}
