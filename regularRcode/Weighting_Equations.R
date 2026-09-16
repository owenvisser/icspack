make_weights <- function(clusterID,
                         clusterSize,
                         K,
                         L,
                         weight_type = "no_weight") {
  
  if (weight_type == "no_weight") {
    
    w <- rep(1, length(clusterID))
    
  } else if (weight_type == "CW") {
    
    w <- 1 / clusterSize
    
  } else if (weight_type == "PPW") {
    
    pair <- paste(K, L, sep = "_")
    
    n_ikl <- ave(
      pair,
      clusterID,
      pair,
      FUN = length
    )
    
    w <- 1 / as.numeric(n_ikl)
    
  } else if (weight_type == "OPW") {
    
    pair <- paste(K, L, sep = "_")
    
    n_ikl <- ave(
      pair,
      clusterID,
      pair,
      FUN = length
    )
    
    N_iKL <- ave(
      pair,
      clusterID,
      FUN = function(z) length(unique(z))
    )
    
    w <- 1 / (as.numeric(N_iKL) * as.numeric(n_ikl))
    
  } else if (weight_type == "MOPW") {
    
    pair <- paste(K, L, sep = "_")
    
    n_ikl <- ave(
      pair,
      clusterID,
      pair,
      FUN = length
    )
    
    N_iK <- ave(
      K,
      clusterID,
      FUN = function(z) length(unique(z))
    )
    
    N_iL <- ave(
      L,
      clusterID,
      FUN = function(z) length(unique(z))
    )
    
    w <- 1 / (as.numeric(N_iK) * as.numeric(N_iL) * as.numeric(n_ikl))
    
  } else {
    
    stop("weight_type must be one of: no_weight, CW, PPW, OPW, MOPW")
    
  }
  
  as.numeric(w)
}