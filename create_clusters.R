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
    lc_row <- as.data.table(cbind(var = "landcoverrc", cluster = max(cluster_dt$cluster) + 1))
    cluster_dt <- rbind(cluster_dt, lc_row)
    
    return(cluster_dt)
}

generate_combinations <- function(vars, clusters) {
  # Combine variables and clusters into a data frame
  var_clusters <- data.frame(var = vars, cluster = clusters)
  
  # Split variables by clusters
  split_vars <- split(var_clusters$var, var_clusters$cluster)
  
  # Generate all combinations of selecting one variable from each unique cluster
  # If the actual clusters are not consecutive numbers or there are gaps, this still works
  combo_list <- lapply(split_vars, function(x) x)
  
  # Use do.call with expand.grid to dynamically generate combinations
  combos <- do.call(expand.grid, combo_list)
  
  # Rename columns for clarity based on cluster numbers
  cluster_names <- paste0('Cluster', seq_along(combo_list))
  colnames(combos) <- cluster_names
  
  return(combos)
}
