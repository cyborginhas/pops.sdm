library(klaR)

cluster_analysis <- function(data, sample_size = 1000000, mincor = 0.5) {
    # Sample the data
    somevars_sample <- spatSample(data, size = sample_size, method = "regular", as.df = TRUE, na.rm = TRUE, values = TRUE, xy = TRUE)
    
    # Convert landcover to factor
    #somevars_sample$landcoverrc <- as.factor(somevars_sample$landcoverrc)
    # Cluster analysis
    #ccres <- corclust(somevars_sample[, -c(1:2)], cl = somevars_sample$landcoverrc)
    ccres <- corclust(somevars_sample[, -c(1:2)]) 
    cors <- ccres$cor
    cluster_dt <- cvtree(ccres, mincor = mincor)
    vars <- rownames(cluster_dt$correlations)
    cluster_dt <- as.data.table(cluster_dt$correlations)
    cluster_dt$var <- vars
    
    # Reassign clusters
    cluster_dt[, clustercount := .N, by = .(cluster)]
    singles <- cluster_dt[cluster_dt$clustercount == 1 & cluster_dt$av.cor2closest < mincor]
    multis <- cluster_dt[cluster_dt$clustercount > 1]
    multisfix <- multis[av.cor2closest > mincor, ]
    multisfix[, cluster := pmin(cluster, closest)]
    multis <- multis[av.cor2closest < mincor, ]
    cluster_dt <- rbind(singles, multis, multisfix)
    cluster_dt <- cluster_dt[, .(var, cluster)]
    cluster_dt <- cluster_dt[order(cluster_dt$cluster, cluster_dt$var)]
    #lc_row <- as.data.table(cbind(var = "landcoverrc", cluster = max(cluster_dt$cluster) + 1))
    #cluster_dt <- rbind(cluster_dt, lc_row)
    
    return(cluster_dt)
}
vars <- cluster_dt$var
clusters <- cluster_dt$cluster

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
  for (num_clusters in 1:length(unique_clusters)) {
    cluster_combos <- combn(unique_clusters, num_clusters, simplify = FALSE)
    
    # Generate variable combinations for each cluster combination
    for (cluster_combo in cluster_combos) {
      # Get the variables for the current combination of clusters
      vars_in_combo <- lapply(cluster_combo, function(cluster) split_vars[[as.character(cluster)]])
      
      # Generate all combinations of these variables across the selected clusters
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

# Function to extract each row of the combinations as a vector of predictors
extract_predictors <- function(combinations) {
  predictors_list <- list() # Initialize an empty list to store the vectors
  
  # Iterate over each element in the combinations list
  combination_names <- names(combinations)
  for (name in combination_names) {
    combination <- combinations[[name]] # Access the combination data frame
    
    # Iterate over each row in the combination data frame
    for (i in 1:nrow(combination)) {
      row_as_vector <- as.character(unlist(combination[i, ])) # Convert row to vector
      predictors_list[[paste(name, i, sep = "_")]] <- row_as_vector # Store the vector
    }
  }
  
  return(predictors_list)
}
