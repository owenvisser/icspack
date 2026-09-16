pearson_association <- function(clusterID,
                                Y,
                                X,
                                omega = NULL) {
  
  if (is.null(omega)) {
    omega <- rep(1, length(X))
  }
  
  clusterID <- as.factor(clusterID)
  M <- length(unique(clusterID))
  
  H <- cbind(
    m10 = X,
    m01 = Y,
    m11 = X * Y,
    m20 = X^2,
    m02 = Y^2
  )
  
  mhat <- colSums(omega * H) / sum(omega)
  
  m10 <- mhat["m10"]
  m01 <- mhat["m01"]
  m11 <- mhat["m11"]
  m20 <- mhat["m20"]
  m02 <- mhat["m02"]
  
  cov_xy <- m11 - m10 * m01
  var_x  <- m20 - m10^2
  var_y  <- m02 - m01^2
  
  rho_hat <- cov_xy / sqrt(var_x * var_y)
  
  U <- sweep(H, 2, mhat, "-")
  omega_U <- omega * U
  
  S_i <- rowsum(omega_U, clusterID)
  
  A_hat <- -(sum(omega) / M) * diag(5)
  B_hat <- crossprod(S_i) / M
  
  Sigma_hat <- solve(A_hat) %*% B_hat %*% solve(A_hat) / M
  
  I_hat <- matrix(c(
    -m01 / sqrt(var_x * var_y) + cov_xy * m10 * var_y / (var_x * var_y)^(3 / 2),
    -m10 / sqrt(var_x * var_y) + cov_xy * m01 * var_x / (var_x * var_y)^(3 / 2),
    1 / sqrt(var_x * var_y),
    -0.5 * cov_xy * var_y / (var_x * var_y)^(3 / 2),
    -0.5 * cov_xy * var_x / (var_x * var_y)^(3 / 2)
  ), nrow = 1)
  
  var_rho <- as.numeric(I_hat %*% Sigma_hat %*% t(I_hat))
  
  se_rho <- sqrt(var_rho)
  
  ci_rho <- c(
    lower = rho_hat - qnorm(0.975) * se_rho,
    upper = rho_hat + qnorm(0.975) * se_rho
  )
  
  list(
    rho = as.numeric(rho_hat),
    variance = var_rho,
    se = se_rho,
    ci = ci_rho,
    mhat = mhat,
    Sigma = Sigma_hat,
    gradient = I_hat
  )
}





spearman_association <- function(clusterID,
                                 Y,
                                 X,
                                 omega = NULL) {
  
  if (is.null(omega)) {
    omega <- rep(1, length(X))
  }
  
  clusterID <- as.factor(clusterID)
  M <- length(unique(clusterID))
  
  weighted_midranks <- function(z, omega) {
    
    ord <- order(z)
    z_ord <- z[ord]
    w_ord <- omega[ord]
    
    total_w <- sum(w_ord)
    
    unique_z <- unique(z_ord)
    
    R_ord <- numeric(length(z_ord))
    
    for (val in unique_z) {
      
      idx <- which(z_ord == val)
      
      F_left <- sum(w_ord[z_ord < val]) / total_w
      F_right <- sum(w_ord[z_ord <= val]) / total_w
      
      R_ord[idx] <- 0.5 * (F_left + F_right)
    }
    
    R <- numeric(length(z))
    R[ord] <- R_ord
    
    R
  }
  
  RX <- weighted_midranks(X, omega)
  RY <- weighted_midranks(Y, omega)
  
  H <- cbind(
    r10 = RX,
    r01 = RY,
    r11 = RX * RY,
    r20 = RX^2,
    r02 = RY^2
  )
  
  rhat <- colSums(omega * H) / sum(omega)
  
  r10 <- rhat["r10"]
  r01 <- rhat["r01"]
  r11 <- rhat["r11"]
  r20 <- rhat["r20"]
  r02 <- rhat["r02"]
  
  cov_r <- r11 - r10 * r01
  var_rx <- r20 - r10^2
  var_ry <- r02 - r01^2
  
  rho_hat <- cov_r / sqrt(var_rx * var_ry)
  
  U <- sweep(H, 2, rhat, "-")
  omega_U <- omega * U
  
  S_i <- rowsum(omega_U, clusterID)
  
  A_hat <- -(sum(omega) / M) * diag(5)
  B_hat <- crossprod(S_i) / M
  
  Sigma_hat <- solve(A_hat) %*% B_hat %*% solve(A_hat) / M
  
  I_hat <- matrix(c(
    -r01 / sqrt(var_rx * var_ry) + cov_r * r10 * var_ry / (var_rx * var_ry)^(3 / 2),
    -r10 / sqrt(var_rx * var_ry) + cov_r * r01 * var_rx / (var_rx * var_ry)^(3 / 2),
    1 / sqrt(var_rx * var_ry),
    -0.5 * cov_r * var_ry / (var_rx * var_ry)^(3 / 2),
    -0.5 * cov_r * var_rx / (var_rx * var_ry)^(3 / 2)
  ), nrow = 1)
  
  var_rho <- as.numeric(I_hat %*% Sigma_hat %*% t(I_hat))
  
  se_rho <- sqrt(var_rho)
  
  ci_rho <- c(
    lower = rho_hat - qnorm(0.975) * se_rho,
    upper = rho_hat + qnorm(0.975) * se_rho
  )
  
  list(
    rho = as.numeric(rho_hat),
    variance = var_rho,
    se = se_rho,
    ci = ci_rho,
    rhat = rhat,
    RX = RX,
    RY = RY,
    Sigma = Sigma_hat,
    gradient = I_hat
  )
}




phi_association <- function(clusterID,
                            Y,
                            X,
                            omega = NULL) {
  
  if (is.null(omega)) {
    omega <- rep(1, length(X))
  }
  
  clusterID <- as.factor(clusterID)
  M <- length(unique(clusterID))
  
  if (!all(X %in% c(0, 1))) {
    stop("X must be binary and coded as 0/1.")
  }
  
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary and coded as 0/1.")
  }
  
  H <- cbind(
    p11 = as.numeric(X == 1 & Y == 1),
    p10 = as.numeric(X == 1 & Y == 0),
    p01 = as.numeric(X == 0 & Y == 1),
    p00 = as.numeric(X == 0 & Y == 0)
  )
  
  phat <- colSums(omega * H) / sum(omega)
  
  p11 <- phat["p11"]
  p10 <- phat["p10"]
  p01 <- phat["p01"]
  p00 <- phat["p00"]
  
  p1_dot <- p11 + p10
  p0_dot <- p01 + p00
  p_dot1 <- p11 + p01
  p_dot0 <- p10 + p00
  
  numerator <- p11 * p00 - p10 * p01
  denominator <- sqrt(p1_dot * p0_dot * p_dot1 * p_dot0)
  
  phi_hat <- numerator / denominator
  
  U <- sweep(H, 2, phat, "-")
  omega_U <- omega * U
  
  S_i <- rowsum(omega_U, clusterID)
  
  A_hat <- -(sum(omega) / M) * diag(4)
  B_hat <- crossprod(S_i) / M
  
  Sigma_hat <- solve(A_hat) %*% B_hat %*% solve(A_hat) / M
  
  phi_function <- function(p) {
    
    p11 <- p[1]
    p10 <- p[2]
    p01 <- p[3]
    p00 <- p[4]
    
    p1_dot <- p11 + p10
    p0_dot <- p01 + p00
    p_dot1 <- p11 + p01
    p_dot0 <- p10 + p00
    
    numerator <- p11 * p00 - p10 * p01
    denominator <- sqrt(p1_dot * p0_dot * p_dot1 * p_dot0)
    
    numerator / denominator
  }
  
  eps <- 1e-6
  
  I_hat <- numeric(4)
  
  for (q in 1:4) {
    
    p_plus <- phat
    p_minus <- phat
    
    p_plus[q] <- p_plus[q] + eps
    p_minus[q] <- p_minus[q] - eps
    
    I_hat[q] <- (phi_function(p_plus) - phi_function(p_minus)) / (2 * eps)
  }
  
  I_hat <- matrix(I_hat, nrow = 1)
  
  var_phi <- as.numeric(I_hat %*% Sigma_hat %*% t(I_hat))
  
  se_phi <- sqrt(var_phi)
  
  ci_phi <- c(
    lower = phi_hat - qnorm(0.975) * se_phi,
    upper = phi_hat + qnorm(0.975) * se_phi
  )
  
  list(
    phi = as.numeric(phi_hat),
    variance = var_phi,
    se = se_phi,
    ci = ci_phi,
    phat = phat,
    Sigma = Sigma_hat,
    gradient = I_hat
  )
}




# Running a grid for the pieces I want.
# 
run_association_grid <- function(dat,
                                 Y_vars,
                                 X_vars,
                                 clusterID,
                                 clusterSize,
                                 K_vars,
                                 L_vars,
                                 weight_type = "no_weight",
                                 association_type = "pearson") {
  
  cid <- dat[[clusterID]]
  cis <- dat[[clusterSize]]
  
  results <- list()
  counter <- 1
  
  if(length(Y_vars) != length(L_vars)) stop("Y_vars and L_vars must match length")
  if(length(X_vars) != length(K_vars)) stop("X_vars and K_vars must match length")
  
  for (l in seq_along(Y_vars)) {
    
    for (k in seq_along(X_vars)) {
      
      if (association_type == "pearson") {
        
        omega <- make_weights(
          clusterID = cid,
          clusterSize = cis,
          K = dat[[K_vars[k]]],
          L = dat[[L_vars[l]]],
          weight_type = weight_type
        )
        
        fit <- pearson_association(
          clusterID = cid,
          Y = dat[[Y_vars[l]]],
          X = dat[[X_vars[k]]],
          omega = omega
        )
        
        estimate <- fit$rho
        
      } else if (association_type == "spearman") {
        
        omega <- make_weights(
          clusterID = cid,
          clusterSize = cis,
          K = dat[[K_vars[k]]],
          L = dat[[L_vars[l]]],
          weight_type = weight_type
        )
        
        fit <- spearman_association(
          clusterID = cid,
          Y = dat[[Y_vars[l]]],
          X = dat[[X_vars[k]]],
          omega = omega
        )
        
        estimate <- fit$rho
        
      } else if (association_type == "phi") {
        
        omega <- make_weights(
          clusterID = cid,
          clusterSize = cis,
          K = dat[[K_vars[k]]],
          L = dat[[L_vars[l]]],
          weight_type = weight_type
        )
        
        fit <- phi_association(
          clusterID = cid,
          Y = dat[[Y_vars[l]]],
          X = dat[[X_vars[k]]],
          omega = omega
        )
        
        estimate <- fit$phi
        
      } else {
        
        stop("association_type must be one of: pearson, spearman, phi")
        
      }
      
      results[[counter]] <- data.frame(
        Y = Y_vars[l],
        L = L_vars[l],
        X = X_vars[k],
        K = K_vars[k],
        weight_type = weight_type,
        association_type = association_type,
        estimate = estimate,
        variance = fit$variance,
        se = fit$se,
        lower_95 = unname(fit$ci[1]),
        upper_95 = unname(fit$ci[2])
      )
      
      counter <- counter + 1
    }
  }
  
  do.call(rbind, results)
}
