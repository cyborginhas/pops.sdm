# This script creates the bias files that are used for Ailanthus altissima SDM.

# Load packages.
# library(raster)

library(spatialEco)
library(data.table)
library(terra)
library(spatstat)

# Set working directory.
setwd("D:/blaginh/new_sdm/")

# Predictors: predictors/truncated/
# Occurrence data: response/truncated/Ailanthus_altissima_occ_pa.csv
# Tree data: target_group_lyr/tree_observations_uscanada.csv

# Read in environmental rasters.
bioclim_filenames <- list.files("predictors/truncated/", full.names = TRUE,
                                pattern = "bioclimatic")

bio_stack <- rast(bioclim_filenames)

# Read in over 1000 occurrence data.
host <- read.csv("response/truncated/Ailanthus_altissima_occ_pa.csv")
host_species <- unique(host$sciname)

# Read in all tree occurrence data.
trees <- fread("bias_files/target_group_lyr/tree_observations_uscanada.csv")

# Create a uniform bias file for modeling biased/unbiased
# Copy a raster from the stack.
uniform_biasfile <- bio_stack[[1]]

# Set all values to 1 (equal probability).
values(uniform_biasfile) <- 1

# Mask the bias file.
uniform_biasfile <- mask(x = uniform_biasfile, mask = bio_stack[[1]])

# Save the bias file for later use.
dir.create("bias_files/uniform/", recursive = TRUE)
## writeRaster(uniform_biasfile, "bias_files/uniform/Uniform_biasfile.tif",
##             overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=ZSTD")


# Create a bias file for the target group.
# Retain records within study extent.
xmin <- ext(bio_stack[[1]])[1]
xmax <- ext(bio_stack[[1]])[2]
ymin <- ext(bio_stack[[1]])[3]
ymax <- ext(bio_stack[[1]])[4]

trees <- trees[trees$longitude >= xmin & trees$longitude <= xmax &
                 trees$latitude >= ymin & trees$latitude <= ymax, ]

# Retain "genus", "species", "latitude", and "longitude" columns.
trees <- unique(trees[, c("genus", "species", "latitude", "longitude")])

# Count the number of unique species within each cell.
sum_records <- trees[, .(n = .N), by = .(longitude, latitude)]
scale <- length(sum_records[,3]) / sum(sum_records[,3])
sum_records$scaled <- sum_records[,3] * scale


# Create sciname column.
trees$sciname <- paste(trees$genus, trees$species, sep = " ")
# Aggregate by coordinates & compute scale.
trees_pts <- vect(trees, geom = c("longitude", "latitude"), crs = "+init=EPSG:4326")

# Define a function that counts the number of unique categories within each cell.
count_unique <- function(x) {
  sum(!is.na(unique(x)))
}

# Rasterize the points file while applying the count_unique function
rasterized <- rasterize(trees_pts, bio_stack[[1]], field="sciname", fun=count_unique, background = NA)

# Compute the scale.
scale <- 1 / global(rasterized, sum, na.rm = TRUE)[1]

# Multiply the rasterized file by the scale.
rasterized <- rasterized * scale$sum
rasterized_df <- cbind(crds(rasterized, na.rm=TRUE), values(rasterized,na.rm=TRUE))
rasterized_df <- as.data.table(rasterized_df)
rasterized_pts <- st_as_sf(rasterized_df, coords = c("x", "y"), crs = 4326)
target_density <- sf.kde(rasterized_pts, rasterized_pts$sciname, ref = bio_stack[[1]], res = res(bio_stack[[1]])[1])


# Do a 2d kernel density estimation.
# Convert sum_records to a sf point object.
sum_records_pts <- st_as_sf(sum_records, coords = c("longitude", "latitude"), crs = 4326)
target_density <- sf.kde(sum_records_pts, sum_records_pts$scale, res = 0.0002694946)


# Create raster.
target_raster <- rast(target_density)

# Define in OSGB.
crs(target_raster) <- '+init=EPSG:27700'

# Clip data to the same resolution/extent.
target_raster <- resample(target_raster, bio_stack, method='bilinear')

# Mask bias file.
for (x in 1:nlayers(bio_stack)){
  target_raster <- mask(target_raster, bio_stack[[x]])
}

# Normalize bias file between 0 and 1.
target_raster <- target_raster - minValue(target_raster)
target_raster <- raster.transformation(target_raster, trans="norm")

# Create a file pathway and export the raster.
writeRaster(target_raster, "Data/Bias Files/TargetGroup_biasfile.tif", overwrite=TRUE)




################################################################################
                   #### Create buffer bias files ####


# Register parallel cores.
library(foreach)
library(doParallel)
registerDoParallel(cores = 8)

# Create a directory to store bias files.
if (!dir.exists("Data/Bias Files/Buffer 10k")){
  dir.create("Data/Bias Files/Buffer 10k", recursive = TRUE)
}

# Create a blank raster to use as a template for each bias file.
blank_raster <- raster(ncol=700, nrow=1300, crs=crs(bio_stack))
extent(blank_raster) <- extent(bio_stack)

# Loop through each hoverfly species to create a unique bias file.
foreach(hoverfly=hoverfly_species,
        .inorder = FALSE,
        .errorhandling = "stop",
        .packages = c("raster", "dplyr", "spatialEco", "ks")) %dopar%
  {
  # Create a species specific dataframe.
  dataframe <- filter(hoverflies, Species == hoverfly)
  # Create spatial points of all the occurrences.
  points <- SpatialPoints(dataframe[,2:3], proj4string = CRS('+init=EPSG:27700'))
  # Rasterize the points using the blank raster.
  raster <- rasterize(points, blank_raster, field=1)
  # Add a 10km buffer around each point.
  raster_buffer_10k <- buffer(raster, 10000)
  # Mask the rasters.
  for (x in 1:nlayers(bio_stack)){
    raster_buffer_10k <- mask(raster_buffer_10k, bio_stack[[x]])
  }
  # Create a name for the pathway and export the bias file.
  pathway <- paste("Data/Bias Files/Buffer 10k/", hoverfly, "_Buffer10k_biasfile.tif", sep="")
  writeRaster(raster_buffer_10k, pathway, overwrite=TRUE)
}

# Create a directory to store bias files.
if (!dir.exists("Data/Bias Files/Virtual Buffer 10k")){
  dir.create("Data/Bias Files/Virtual Buffer 10k", recursive = TRUE)
}

# Loop through each virtual species to create a unique bias file.
foreach(hoverfly=virtual_species,
        .inorder = FALSE,
        .errorhandling = "stop",
        .packages = c("raster", "dplyr", "spatialEco", "ks")) %dopar%
  {
    # Create a species specific dataframe.
    dataframe <- filter(virtual_occurences, random_name == hoverfly) %>% select(random_name, PO_x_biased, PO_y_biased)
    # Create spatial points of all the occurrences.
    points <- SpatialPoints(dataframe[,2:3], proj4string = CRS('+init=EPSG:27700'))
    # Rasterize the points using the blank raster.
    raster <- rasterize(points, blank_raster, field=1)
    # Add a 10km buffer around each point.
    raster_buffer_10k <- buffer(raster, 10000)
    # Mask the raster.
    for (x in 1:nlayers(bio_stack)){
      raster_buffer_10k <- mask(raster_buffer_10k, bio_stack[[x]])
    }
    # Create a name for the pathway and export the bias file.
    pathway <- paste("Data/Bias Files/Virtual Buffer 10k/", hoverfly, "_Buffer10k_biasfile.tif", sep="")
    writeRaster(raster_buffer_10k, pathway, overwrite=TRUE)
  }




################################################################################
             #### Create Population Density Bias File ####


# Load in Gridded Human Population Density of the World.
Pop_Dens_East_Europe <- raster("Data/Population Density/gpw_v4_population_density_rev11_2000_30_sec_3.asc")
Pop_Dens_West_Europe <- raster("Data/Population Density/gpw_v4_population_density_rev11_2000_30_sec_2.asc")

# Merge east and west europe into one raster. # Takes a long time.
Europe <- raster::merge(Pop_Dens_West_Europe, Pop_Dens_East_Europe)

# Crop to roughly uk extent in WGS.
WGS_extent <- c(-8, 2, 49, 62)
Pop_Dens_UK_WGS <- crop(Europe, WGS_extent)
plot(Pop_Dens_UK_WGS)

# Project into OSGB.
Pop_Dens_UK <- projectRaster(Pop_Dens_UK_WGS, crs ='+init=EPSG:27700', method = 'bilinear')

# Clip map to UK resolution.
Pop_Dens_UK <- resample(Pop_Dens_UK, bio_stack, method="bilinear")

# Mask map and plot.
for (x in 1:nlayers(bio_stack)){
  Pop_Dens_UK <- mask(Pop_Dens_UK, bio_stack[[x]])
}

# Plot UK population density.
plot(Pop_Dens_UK)

# Values are logged to reduce overcompensation.
UK_Log <- log(Pop_Dens_UK)

# '-inf' values from logging 0 values are set to NAs.
remove_infs <- function(x) { x[x==-Inf] <- NA; return(x)}
UK_Log <- calc(UK_Log, remove_infs)
summary(UK_Log)

#Values are shifted to a positive number (min value is approx -6.3).
UK_Log_Pos <- UK_Log - minValue(UK_Log)

# Normalize bias file to between 0 and 1.
UK_Log_Pos <- raster.transformation(UK_Log_Pos, trans="norm")
plot(UK_Log_Pos)

# Export bias file.
writeRaster(UK_Log_Pos, "Data/Bias Files/PopDens_biasfile.tif", overwrite=TRUE)




################################################################################
            #### Create a bias file from travel time ####


# Load in Gridded Travel Time of the World.
Travel_Time <- raster("Data/Travel Time/acc_50k.tif")

# Crop to roughly the UK.
Travel_Time_UK <<- crop(Travel_Time, WGS_extent)

# Project into OSGB.
Travel_Time_UK <- projectRaster(Travel_Time_UK, crs ='+init=EPSG:27700', method = 'bilinear')

# Clip map to UK resolution.
Travel_Time_UK <- resample(Travel_Time_UK, bio_stack, method="bilinear")

# Mask map.
for (x in 1:nlayers(bio_stack)){
  Travel_Time_UK <- mask(Travel_Time_UK, bio_stack[[x]])
}

# Plot raw travel time.
plot(Travel_Time_UK)

# Values are logged to reduce overcompensation.
Travel_Time_Log <- log(Travel_Time_UK)

# '-inf' values from logging 0 values are set to NAs.
Travel_Time_Log <- calc(Travel_Time_Log, remove_infs)
summary(Travel_Time_Log)
plot(Travel_Time_Log)

# Values are negated to reverse the relationship (so that the most remote 
# locations have the lowest probability of being selected).
Travel_Time_Log <- Travel_Time_Log * -1

# Values are shifted to a positive number (min value is approx -6.6).
Travel_Time_Log_Pos <- Travel_Time_Log - minValue(Travel_Time_Log)

# Normalize bias file to between 0 and 1.
Travel_Time_Log_Pos <- raster.transformation(Travel_Time_Log_Pos, trans="norm")
plot(Travel_Time_Log_Pos)

# Export bias file.
writeRaster(Travel_Time_Log_Pos, "Data/Bias Files/TravelTime_biasfile.tif", overwrite=TRUE)

