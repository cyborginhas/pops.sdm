# This script creates the bias files that are used for Ailanthus altissima SDM.
# Load packages.
# library(raster)

library(spatialEco)
library(data.table)
library(terra)
library(spatstat)
library(sf)

# Set working directory.
setwd("D:/blaginh/new_sdm/")

# Predictors: predictors/truncated/
# Occurrence data: response/truncated/Ailanthus_altissima_occ_pa.csv
# Tree data: target_group_lyr/tree_observations_uscanada.csv

# Read in environmental rasters.
bioclim_filenames <- list.files("predictors/truncated/",
  full.names = TRUE,
  pattern = "bioclimatic"
)

bio_stack <- rast(bioclim_filenames)

# Read in over 1000 occurrence data.
host <- read.csv("response/truncated/Ailanthus_altissima_occ_pa.csv")
host_species <- unique(host$sciname)

# Read in all tree occurrence data.
# trees <- fread("bias_files/target_group_lyr/tree_observations_uscanada.csv")

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
## xmin <- ext(bio_stack[[1]])[1]
## xmax <- ext(bio_stack[[1]])[2]
## ymin <- ext(bio_stack[[1]])[3]
## ymax <- ext(bio_stack[[1]])[4]

## trees <- trees[trees$longitude >= xmin & trees$longitude <= xmax &
##   trees$latitude >= ymin & trees$latitude <= ymax, ]

## # Retain "genus", "species", "latitude", and "longitude" columns.
## trees <- unique(trees[, c("genus", "species", "latitude", "longitude")])

## # Create sciname column.
## trees$sciname <- paste(trees$genus, trees$species, sep = " ")


# Write out the tree data.
## fwrite(trees, "bias_files/target_group_lyr/tree_observations_pa_ext.csv",
##           row.names = FALSE)

# Read in the tree data.
## trees <- fread("bias_files/target_group_lyr/tree_observations_pa_ext.csv")
## trees <- trees[, .(count = .N), by = .(longitude, latitude)]

## # Convert trees to points
## trees_pts <- vect(trees, geom = c("longitude", "latitude"))

# function to count unique values when rasterizing
## count_unique <- function(x) {
##   length(unique(x))
## }


# Create a raster from the tree data.
## scale <- 1 / sum(trees_pts$count)
## trees_pts$scale <- trees_pts$count * scale

# Do a 2d kernel density estimation.
# target_density <- ks::kde(crds(trees_pts), w=trees_pts$scale, gridsize = c(9455,21618))
# target_raster <- rast(bio_stack[[1]])
# values(target_raster) <- target_density$estimate

## # Write out the target bias file.
## dir.create("bias_files/target_group_lyr/", recursive = TRUE)
## # writeRaster(target_raster, "bias_files/target_group_lyr/TargetGroup_biasfile_raw.tif", datatype = "FLT4S", gdal = "COMPRESS=ZSTD", overwrite = TRUE)
## target_raster <- rast("bias_files/target_group_lyr/TargetGroup_biasfile_raw.tif")

## # Normalize the target bias file between 0 and 1.
## min <- global(target_raster, "min")
## min <- min$min

## f <- function(i) i - min
## target_raster <- app(target_raster, f)
## target_raster <- raster.transformation(target_raster, trans = "norm")

## writeRaster(target_raster, "bias_files/target_group_lyr/TargetGroup_biasfile.tif",
##   overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD"
## )


################################################################################
#### Create buffer bias files ####
## # Create a directory to store bias files.
## dir.create("bias_files/buffers", recursive = TRUE)

## # Create spatial points of all the occurrences.
## host_pts <- vect(host, geom = c("lon", "lat"), crs = "+init=EPSG:4326")

## # Add 5km, 10km, 15km buffer around each point.
## buffers <- c(2000, 8000, 14000, 20000)

## for(i in 1:length(buffers)) {
##   buffer <- buffer(host_pts, buffers[i])
##   buffer$group <- 1
##   buffer <- aggregate(buffer, "group")
##   raster_buffer <- rasterize(buffer, bio_stack[[1]], field = 1)
##   raster_buffer <- mask(raster_buffer, bio_stack[[1]])
##   writeRaster(raster_buffer, paste0("bias_files/buffers/buffer_", buffers[i],"m_biasfile.tif"),
##               overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=ZSTD")
## }

################################################################################
#### Create Population Density Bias File ####
## pop_dens <- rast("predictors/truncated/meta_population_human_pop_density_30m_ext25.tif")

## # Values are logged to reduce overcompensation.
## pop_dens <- log(pop_dens)
## min <- min(values(pop_dens, na.rm = TRUE))
## f <- function(i) i - min
## pop_dens <- app(pop_dens, f)

## # Normalize bias file to between 0 and 1.
## pop_dens <- raster.transformation(pop_dens, trans = "norm")

## dir.create("bias_files/population", recursive = TRUE)
## writeRaster(pop_dens, "bias_files/population/PopDens_biasfile.tif", overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=ZSTD")

################################################################################
#### Create a bias file from dist2roads ####

# Load in Gridded Travel Time of the World.
dist2roads <- rast("predictors/truncated/tiger_roads_dist2roads_30m_ext25.tif")

# Values are logged to reduce overcompensation.
dist2roads <- log(dist2roads)

# Set all -Inf values to NA.
remove_infs <- function(x) { x[x == -Inf] <- NA; return(x) }
dist2roads <- app(dist2roads, remove_infs)

# Values are negated to reverse the relationship (so that the most remote


# Values are shifted to a positive number (min value is approx -6.6).
min <- global(dist2roads, "min", na.rm = TRUE)$min
f <- function(i) i - min
dist2roads <- app(dist2roads, f)

# Normalize bias file to between 0 and 1.
dist2roads <- raster.transformation(dist2roads, trans = "norm")

# Export bias file.
setwd("D:/blaginh/new_sdm/")
dir.create("bias_files/dist2roads", recursive = TRUE)
writeRaster(dist2roads, "bias_files/dist2roads/dist2roads_bias_file.tif", overwrite = TRUE,
            datatype = "FLT4S", gdal = "COMPRESS=ZSTD")

# Repeat for dist2rails
dist2rails <- rast("predictors/truncated/tiger_railroads_dist2rails_30m_ext25.tif")

# Values are logged to reduce overcompensation.
dist2rails <- log(dist2rails)

# Set all -Inf values to NA.
dist2rails <- app(dist2rails, remove_infs)

# Values are negated to reverse the relationship (so that the most remote
# areas have the highest values).
dist2rails <- app(dist2rails, function(x) -x)

# Values are shifted to a positive number (min value is approx -6.6).
min <- global(dist2rails, "min", na.rm = TRUE)$min
dist2rails <- app(dist2rails, function(x) x - min)

# Normalize bias file to between 0 and 1.
dist2rails <- raster.transformation(dist2rails, trans = "norm")

# Export bias file.
dir.create("bias_files/dist2rails", recursive = TRUE)
writeRaster(dist2rails, "bias_files/dist2rails/dist2rails_bias_file.tif", overwrite = TRUE,
            datatype = "FLT4S", gdal = "COMPRESS=ZSTD")