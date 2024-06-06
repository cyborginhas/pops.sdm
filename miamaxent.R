# Load library
library(MIAmaxent)
library(terra)
library(data.table)

# Load data
occs_file <- "D:/blaginh/sdm_full_extent/flexsdm_results/1_Inputs/1_Occurrences/background/target_group/Ailanthus_altissima_target_group_bg_pts_filtgeo_120m2_wpreds.csv"

# Copy to "C:/Users/blaginh/Desktop/maximize/"
#file.copy(occs_file, "C:/Users/blaginh/Desktop/maximize/")
occs_file <- "C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_target_group_bg_pts_filtgeo_120m2_wpreds.csv"
occs_120_tg_data <- fread(occs_file)

# Set names .part to "part"
occs_120_tg_data <- occs_120_tg_data[, part := .part]

# Use setnames to change NLCD Land Cover Class to landcoverrc
setnames(occs_120_tg_data, "NLCD Land Cover Class", "landcoverrc")

# Convert "NLCD Land Cover Class" to "landcoverrc"
occs_120_tg_data <- occs_120_tg_data[, .(
    pr_ab, x, y, part, d, elevation, slope, aspect,
    bio19, bio15, evi, par, clay_0_5, silt_0_5, ph_0_5, theta_s_0_5, landcoverrc
)]

# Only keep complete cases
occs_120_tg_data <- occs_120_tg_data[complete.cases(occs_120_tg_data), ]

# Convert pr_ab 0 to NA
occs_120_tg_data$pr_ab[occs_120_tg_data$pr_ab == 0] <- NA

# Convert pr_ab to RV
setnames(occs_120_tg_data, "pr_ab", "RV")

# Convert landcoverrc to factor
occs_120_tg_data$landcoverrc <- as.factor(occs_120_tg_data$landcoverrc)

# Filter by part column
occs_part <- list()
for (i in 1:4) {
    occs_part[[i]] <- occs_120_tg_data[part == i]
}

# Drop d, x, y, and part in for loop
for (i in 1:4) {
    occs_part[[i]] <- occs_part[[i]][, .(RV, elevation, slope, aspect, bio19, bio15, evi, par, clay_0_5, silt_0_5, ph_0_5, theta_s_0_5, landcoverrc)]
}

# Partition
# Define a function to perform variable selection and model fitting
perform_variable_selection <- function(occurrence) {
    # Derive variables
    occs_dvs <- deriveVars(occurrence, transformtype = c("L", "M", "D", "HF", "HR", "T", "B"))

    # Select dependent variables for environmental variables
    occs_dvs_select <- selectDVforEV(occs_dvs$dvdata, alpha = 0.001, quiet = TRUE)

    # Select environmental variables
    occs_EVselect <- selectEV(occs_dvs_select$dvdata, alpha = 0.001, interaction = TRUE)

    # Return the selected variables and selected dependent variables
    return(list(occs_EVselect = occs_EVselect, occs_dvs_select = occs_dvs_select, occs_dvs = occs_dvs))
}


# Call the function with the occurrence data
part123 <- rbindlist(occs_part[1:3])
part134 <- rbindlist(occs_part[c(1, 3, 4)])
part124 <- rbindlist(occs_part[c(1, 2, 4)])
part234 <- rbindlist(occs_part[2:4])
full <- rbindlist(occs_part)

# Perform variable selection per validation part
part123_select <- perform_variable_selection(part123)
part134_select <- perform_variable_selection(part134)
part124_select <- perform_variable_selection(part124)
part234_select <- perform_variable_selection(part234)
full_select  <- perform_variable_selection(full)

# Per part plot the Dsq
plot(part123_select$occs_EVselect$selection$round, part123_select$occs_EVselect$selection$Dsq,
    xlab = "round", ylab = "Dsq", main = "Part 1, 2, 3"
)
part123_select$occs_EVselect$selection$Dsq
index123 <- which(part123_select$occs_EVselect$selection$Dsq > 0.017)[1]
model123_formula <- formula(paste("~", paste(part123_select$occs_EVselect$selection$variables[index123], collapse = " + ")))
model123 <- chooseModel(part123_select$occs_EVselect$dvdata, model123_formula)
varimp123 <- calculateFTVA(part123_select$occs_EVselect, model123_formula)

plot(part134_select$occs_EVselect$selection$round, part134_select$occs_EVselect$selection$Dsq,
    xlab = "round", ylab = "Dsq", main = "Part 1, 3, 4"
)
part134_select$occs_EVselect$selection$Dsq
index134 <- which(part134_select$occs_EVselect$selection$Dsq > 0.017)[1]
model134_formula <- formula(paste("~", paste(part134_select$occs_EVselect$selection$variables[index134], collapse = " + ")))
model134 <- chooseModel(part134_select$occs_EVselect$dvdata, model134_formula)
varimp134 <- calculateFTVA(part134_select$occs_EVselect, model134_formula)

plot(part124_select$occs_EVselect$selection$round, part124_select$occs_EVselect$selection$Dsq,
    xlab = "round", ylab = "Dsq", main = "Part 1, 2, 4"
)
part124_select$occs_EVselect$selection$Dsq
index124 <- which(part124_select$occs_EVselect$selection$Dsq > 0.017)[1]
model124_formula <- formula(paste("~", paste(part124_select$occs_EVselect$selection$variables[index124], collapse = " + ")))
model124 <- chooseModel(part124_select$occs_EVselect$dvdata, model124_formula)
varimp124 <- calculateFTVA(part124_select$occs_EVselect, model124_formula)

plot(part234_select$occs_EVselect$selection$round, part234_select$occs_EVselect$selection$Dsq,
    xlab = "round", ylab = "Dsq", main = "Part 2, 3, 4"
)
part234_select$occs_EVselect$selection$Dsq
index234 <- which(part234_select$occs_EVselect$selection$Dsq > 0.017)[1]
model234_formula <- formula(paste("~", paste(part234_select$occs_EVselect$selection$variables[index234], collapse = " + ")))
model234 <- chooseModel(part234_select$occs_EVselect$dvdata, model234_formula)
varimp234 <- calculateFTVA(part234_select$occs_EVselect, model234_formula)

plot(full_select$occs_EVselect$selection$round, full_select$occs_EVselect$selection$Dsq,
    xlab = "round", ylab = "Dsq", main = "Full"
)
full_select$occs_EVselect$selection$Dsq
index_full <- which(full_select$occs_EVselect$selection$Dsq > 0.021)[1]
model_full_formula <- formula(paste("~", paste(full_select$occs_EVselect$selection$variables[index_full], collapse = " + ")))
model_full <- chooseModel(full_select$occs_EVselect$dvdata, model_full_formula)
varimp_full <- calculateFTVA(full_select$occs_EVselect, model_full_formula)

# Rbind the varimp
varimp <- rbind(varimp123, varimp134, varimp124, varimp234, varimp_full)
varimp <- data.table(varimp)
# calculate the mean of the varimp by pred
varimp_median <- varimp[, .(mean_varimp = median(FTVA)), by = variable]
# order by mean_varimp
varimp_median <- varimp_median[order(-mean_varimp)]

# Get preds in varimp_median
preds <- varimp_median$variable
# grep out any preds with ":"
preds <- preds[!grepl(":", preds)]
# gsub . with _
preds <- gsub("\\.", "_", preds)
## # rename raster names
## names(cropped_predictors)[names(cropped_predictors) == "NLCD Land Cover Class"] <- "landcoverrc"

## best_predictors <- list()
## i <- 
## for (i in seq_along(preds)) {
##     best_predictors[[i]] <- terra::sources(cropped_predictors[[names(cropped_predictors) == preds[i]]])
##     file.copy(best_predictors[[i]], paste0("C:/Users/blaginh/Desktop/maximize/", basename(best_predictors[[i]])))
##     best_predictors[[i]] <- terra::rast(paste0("C:/Users/blaginh/Desktop/maximize/", basename(best_predictors[[i]])))
## }

best_predictors <- list.files("C:/Users/blaginh/Desktop/maximize/", full.names = TRUE, pattern = ".tif")
best_predictors <- terra::rast(best_predictors)
names(best_predictors)[names(best_predictors) == "NLCD Land Cover Class"] <- "landcoverrc"

# Test the model
# Replace "." in column names with "_"
occs_part <- lapply(seq_along(occs_part), function(x) {
    part <- occs_part[[x]]
    names(part) <- gsub("_", ".", names(part))
    part
})

# Test AUC
testauc123 <- testAUC(model123, part123_select$occs_dvs$transformations, occs_part[[4]], plot = FALSE)
testauc134 <- testAUC(model134, part134_select$occs_dvs$transformations, occs_part[[2]], plot = FALSE)
testauc124 <- testAUC(model124, part124_select$occs_dvs$transformations, occs_part[[3]], plot = FALSE)
testauc234 <- testAUC(model234, part234_select$occs_dvs$transformations, occs_part[[1]], plot = FALSE)

# Get preds in varimp_median
preds <- varimp_full$variable
# grep out any preds with ":"
preds <- preds[!grepl(":", preds)]
# gsub . with _
preds <- gsub("\\.", "_", preds)

best_predictors_full <- list()


for (i in seq_along(preds)) {
    best_predictors_full[[i]] <- best_predictors[[names(best_predictors) == preds[i]]]
    #file.copy(best_predictors[[i]], paste0("C:/Users/blaginh/Desktop/maximize/", basename(best_predictors[[i]])))
    #best_predictors[[i]] <- terra::rast(paste0("C:/Users/blaginh/Desktop/maximize/", basename(best_predictors[[i]])))
}
best_predictors_full <- terra::rast(best_predictors_full)
#' Read in testing set
testing <- fread("C:/Users/blaginh/Desktop/maximize/Ailanthus_altissima_testing_sets.csv")
testing_pts <- terra::vect(testing, crs = terra::crs(best_predictors_full), geom = c("x", "y"))

# Extract the environmental data
testing_data <- terra::extract(best_predictors_full, testing_pts, ID = FALSE, xy = TRUE, bind = TRUE, method = "simple")
testing_data <- data.table::as.data.table(testing_data)
names(testing_data) <- gsub("_", ".", names(testing_data))
testAUC(model_full, full_select$occs_dvs$transformations, testing_data[setid == 5,], plot = FALSE)

varimp <- varimp_full
predictors <- best_predictors
model <- model_full
cycle <- "full"
occs_dvs_select <- full_select$occs_dvs_select

# Create stack of predictors based on the best predictors per part
create_stack <- function(varimp, model, predictors, cycle, occs_dvs_select) {
    s<- Sys.time()
    varimp_pred <- varimp$variable
    varimp_pred <- varimp_pred[!grepl(":", varimp_pred)]
    varimp_pred <- gsub("\\.", "_", varimp_pred)

    stack <- list()
    for (i in seq_along(varimp_pred)) {
        stack[[i]] <- raster::raster(predictors[[names(predictors) == varimp_pred[i]]])
    }

    stack <- raster::stack(stack)

    preds <- projectModel(
        model = model,
        transformations = occs_dvs_select$transformations,
        data = stack
    )
    raster::writeRaster(preds, paste0("C:/Users/blaginh/Desktop/maximize/",
                        "projectModel_", cycle, ".tif"), overwrite = TRUE)
    t <- Sys.time() - s
    print(t)
    return(preds)
}

stack123 <- create_stack(varimp123, model123, best_predictors, "part123", part123_select$occs_dvs_select)
stack134 <- create_stack(varimp134, model134, best_predictors, "part134", part134_select$occs_dvs_select)
stack124 <- create_stack(varimp124, model124, best_predictors, "part124", part124_select$occs_dvs_select)
stack234 <- create_stack(varimp234, model234, best_predictors, "part234", part234_select$occs_dvs_select)
