# Load required libraries
library(data.table)
library(klaR)
library(terra)

#' @description Perform z-score normalization on a dataset
#' @param var A numeric vector containing the variable to be normalized
#' @return A data table containing the normalized variables
#' @export

z_score_normalization <- function(var) {
  mean_var <- mean(var, na.rm = TRUE)
  sd_var <- sd(var, na.rm = TRUE)
  normalized_var <- (var - mean_var) / sd_var
  return(normalized_var)
}

#' @description Perform cluster analysis on a dataset
#' @param data A SpatRaster object containing the binary/integer/numeric
#' variables to be used for clustering
#' @param sample_size The number of samples to use for the cluster analysis;
#' if memory is an issue, reduce this number.
#' @param mincor The minimum correlation threshold for clustering
#' @param categorical_vars A SpatRaster object containing categorical variables;
#' these variables will be treated as classes in the cluster analysis.
#' @return A data table containing the variables and their assigned clusters

cluster_analysis <- function(data, sample_size, mincor, categorical_vars) {
  # Perform cluster analysis
  tryCatch(
      {
        sample <- terra::spatSample(data,
          size = sample_size, method = "regular",
          as.df = TRUE, na.rm = TRUE, values = TRUE,
          xy = TRUE
        )
        # Z-score normalization: (x - mean(x)) / sd(x)) apply to each column
        col_names <- names(sample)
        sample <- as.data.table(lapply(sample, z_score_normalization))
        names(sample) <- paste0(col_names, "_zscore")
        # Keep all columns except 1:2 and cat_cols
        ccres <- klaR::corclust(sample[,-c(1:2)])
      },
      error = function(e) {
        message("An error occurred: ", e$message)
        message("Please increase the sample_size.")
      }
    )
  # Create hierarchical cluster tree based correlation threshold
  cluster_dt <- klaR::cvtree(ccres, mincor = mincor)
  vars <- rownames(cluster_dt$correlations)
  cluster_dt <- as.data.table(cluster_dt$correlations)
  cluster_dt$var <- vars

  # Reassign clusters
  cluster_dt[, clustercount := .N, by = .(cluster)]
  singles <- cluster_dt[cluster_dt$clustercount == 1 &
    cluster_dt$av.cor2closest < mincor]
  multis <- cluster_dt[cluster_dt$clustercount > 1]
  multisfix <- multis[av.cor2closest > mincor, ]
  multisfix[, cluster := pmin(cluster, closest)]
  multis <- multis[av.cor2closest < mincor, ]
  cluster_dt <- rbind(singles, multis, multisfix)
  cluster_dt <- cluster_dt[, .(var, cluster)]
  cluster_dt <- cluster_dt[order(cluster_dt$cluster, cluster_dt$var)]
  # Add categorical variables to the cluster assignments
  if (!is.null(categorical_vars)) {
    cat_vars <- names(categorical_vars)
    cat_clusters <- max(cluster_dt$cluster) + 1:length(cat_vars)
    cat_vars <- data.table(var = cat_vars, cluster = cat_clusters)
    cluster_dt <- rbind(cluster_dt, cat_vars)
  }
  return(cluster_dt)
}


#' @description Generate combinations of variables based on clusters
#' @param vars A character vector containing the variable names
#' @param clusters A character vector containing the cluster assignments
#' @return A list containing data frames of variable combinations for each
#' cluster count

generate_combinations <- function(vars, clusters) {
  # Ensure clusters are treated as characters for consistent indexing
  clusters <- as.character(clusters)
  # Split the variables by their cluster
  split_vars <- split(vars, clusters)
  # Prepare a list to store the combinations for each cluster count
  all_combinations <- list()
  # Get all unique cluster numbers
  unique_clusters <- unique(clusters)
  # Generate combinations for selecting clusters
  for (num_clusters in seq_along(unique_clusters)) {
    cluster_combos <- combn(unique_clusters, num_clusters, simplify = FALSE)
    # Generate variable combinations for each cluster combination
    for (cluster_combo in cluster_combos) {
      # Get the variables for the current combination of clusters
      vars_in_combo <-
        lapply(
          cluster_combo,
          function(cluster) {
            split_vars[[as.character(cluster)]]
          }
        )
      if (length(vars_in_combo) > 0) {
        var_combos <- expand.grid(vars_in_combo, stringsAsFactors = FALSE)
        # Create a unique key for the combination
        key <- paste(cluster_combo, collapse = "-")
        all_combinations[[key]] <- var_combos
      }
    }
  }
  return(all_combinations)
}

#' @description Extract predictors from a list of combinations to be used in
#' in best variables selection
#' @param combinations A list containing data frames of variable combinations
#' @return A list containing vectors of variable names for each combination

extract_predictors <- function(combinations) {
  predictors <- lapply(combinations, function(combo) {
    combo$Var1
  })
  return(predictors)
}
