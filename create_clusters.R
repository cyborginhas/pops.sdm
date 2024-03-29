# Load required libraries
library(data.table)
library(klaR)
library(terra)

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
  # Sample data
  sample <- terra::spatSample(data, size = sample_size, method = "regular",
                              as.df = TRUE, na.rm = TRUE, values = TRUE,
                              xy = TRUE)
  # Perform cluster analysis
  if (!is.null(categorical_vars)) {
    ccres <- klaR::corclust(sample[, -c(1:2)], categorical_vars)
  } else {
    ccres <- klar::corclust(sample[, -c(1:2)])
  }
  # Create hierarchical cluster tree based correlation threshold
  cluster_dt <- klar::cvtree(ccres, mincor = mincor)
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
        lapply(cluster_combo,
               function(cluster) {
                 split_vars[[as.character(cluster)]]
               })
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
  predictors_list <- list() # Initialize an empty list to store the vectors

  # Iterate over each element in the combinations list
  combination_names <- names(combinations)
  for (name in combination_names) {
    combination <- combinations[[name]] # Access the combination data frame

    # Iterate over each row in the combination data frame
    for (i in seq_len(nrow(combination))) {
      row_as_vector <- as.character(unlist(combination[i, ]))
      predictors_list[[paste(name, i, sep = "_")]] <- row_as_vector
    }
  }
  return(predictors_list)
}
