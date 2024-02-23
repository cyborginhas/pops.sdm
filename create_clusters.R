# Biovars clusters - global
path <- "Z:/pops_pesthostuse/pops.sdm/Data/"
biovars <- get_biovars_global(path)
biovars_names <- unlist(lapply(biovars, names))
biovars_files <- unlist(lapply(biovars, function(x){
    terra::sources(x)
}))
biocl_global <- data.frame(var = biovars_names, cluster = NA, domain = "Global", paths = biovars_files)
biocl_global$cluster[which(biocl_global$var %in% c(
    "wc2.1_30s_bio_1",
    "wc2.1_30s_bio_5",
    "wc2.1_30s_bio_8",
    "wc2.1_30s_bio_10"
))] <- "Temp 1"
biocl_global$cluster[which(biocl_global$var %in% c(
    "wc2.1_30s_bio_2",
    "wc2.1_30s_bio_3",
    "wc2.1_30s_bio_4",
    "wc2.1_30s_bio_6",
    "wc2.1_30s_bio_7",
    "wc2.1_30s_bio_9",
    "wc2.1_30s_bio_11"
))] <- "Temp 2"
biocl_global$cluster[which(biocl_global$var %in% c(
    "wc2.1_30s_bio_12",
    "wc2.1_30s_bio_13",
    "wc2.1_30s_bio_16",
    "wc2.1_30s_bio_18",
    "wc2.1_30s_bio_19"
))] <- "Precip 1"
biocl_global$cluster[which(biocl_global$var %in% c(
    "wc2.1_30s_bio_15",
    "wc2.1_30s_bio_14",
    "wc2.1_30s_bio_17"
))] <- "Precip 2"
biocl_global

# Biovars clusters - USA
biovars <- get_biovars_conus(path)
biovars_names <- unlist(lapply(biovars, names))
biovars_files <- unlist(lapply(biovars, function(x){
    terra::sources(x)
}))
biocl_conus <- data.frame(var = biovars_names, cluster = NA, domain = "USA", paths = biovars_files)
biocl_conus$cluster[which(biocl_conus$var %in% c(
    "bio1",
    "bio5",
    "bio8",
    "bio10"
))] <- "Temp 1"
biocl_conus$cluster[which(biocl_conus$var %in% c(
    "bio2",
    "bio3",
    "bio4",
    "bio4a",
    "bio6",
    "bio7",
    "bio9",
    "bio11"
))] <- "Temp 2"
biocl_conus$cluster[which(biocl_conus$var %in% c(
    "bio12",
    "bio13",
    "bio16",
    "bio18",
    "bio19"
))] <- "Precip 1"
biocl_conus$cluster[which(biocl_conus$var %in% c(
    "bio15",
    "bio14",
    "bio17"
))] <- "Precip 2"
biocl_conus

# Elevation clusters - global
elev <- get_topo_global(path)
elev_names <- unlist(lapply(elev, names))
elev_files <- unlist(lapply(elev, function(x){
    terra::sources(x)
}))
elevcl_global <- data.frame(var = elev_names, cluster = NA, domain = "Global", paths = elev_files)
elevcl_global[elevcl_global$var == "wc2.1_30s_elev",]$cluster <- "Elevation 1"
elevcl_global[elevcl_global$var == "hillshade",]$cluster <- "Elevation 2"

# Elevation clusters - USA
elev <- get_topo_conus(path)
elev_names <- unlist(lapply(elev, names))
elev_files <- unlist(lapply(elev, function(x){
    terra::sources(x)
}))
elevcl_conus <- data.frame(var = elev_names, cluster = NA, domain = "USA", paths = elev_files)
elevcl_conus[elevcl_conus$var == "elevation",]$cluster <- "Elevation 1"
elevcl_conus[elevcl_conus$var == "hillshade",]$cluster <- "Elevation 2"

# GDD - global
gdd <- get_gdd_global(path)
gdd_names <- unlist(lapply(gdd, names))
gdd_files <- unlist(lapply(gdd, function(x){
    terra::sources(x)
}))
gddcl <- data.frame(var = gdd_names, cluster = "Temp 1", domain = "Global", paths = gdd_files)

# Land cover - global
lc <- get_landcover_global(path)
lc_names <- unlist(lapply(lc, names))
lc_files <- unlist(lapply(lc, function(x){
    terra::sources(x)
}))
lccl_global <- data.frame(var = lc_names, cluster = "Landcover 1", domain = "Global", paths = lc_files)

# Land cover - USA
lc <- get_landcover_conus(path)
lc_names <- unlist(lapply(lc, names))
lc_files <- unlist(lapply(lc, function(x){
    terra::sources(x)
}))
lccl_conus <- data.frame(var = lc_names, cluster = "Landcover 1", domain = "USA", paths = lc_files)

# Population density - global
pop <- get_pop_global(path)
pop_names <- unlist(lapply(pop, names))
pop_files <- unlist(lapply(pop, function(x){
    terra::sources(x)
}))
popcl <- data.frame(var = pop_names, cluster = "Population", domain = "Global", paths = pop_files)

# Precipitation Timing - global
ppt <- get_prectiming_global(path)
ppt_names <- unlist(lapply(ppt, names))
ppt_files <- unlist(lapply(ppt, function(x){
    terra::sources(x)
}))
pptcl <- data.frame(var = ppt_names, cluster = "Precip 2", domain = "Global", paths = ppt_files)

# Railroads - global
rr <- get_rails_global(path)
rr_names <- unlist(lapply(rr, names))
rr_files <- unlist(lapply(rr, function(x){
    terra::sources(x)
}))
rrcl <- data.frame(var = rr_names, cluster = "Roads/Rails", domain = "Global", paths = rr_files)

# Roads - global
roads <- get_roads_global(path)
roads_names <- unlist(lapply(roads, names))
roads_files <- unlist(lapply(roads, function(x){
    terra::sources(x)
}))
rcl <- data.frame(var = roads_names, cluster = "Roads/Rails", domain = "Global", paths = roads_files)

# Soil - global
soil <- get_soilvars_global(path)
