HiperGator = T

if(HiperGator){
  # On HiperGator
  library(MASS)
  library(parallel)
} else {
  library(MASS)
}


simulate_one_cluster <- function(cluster_id, p, disc = F, cuts_x, cuts_y) {
  
  Sigma_uv <- matrix(c(p$sigma_u^2,
                       p$sigma_u * p$sigma_v * p$rho_uv,
                       p$sigma_u * p$sigma_v * p$rho_uv,
                       p$sigma_v^2),
                     nrow = 2)
  uv <- matrix(MASS::mvrnorm(n = 1, mu = c(p$mu_u, p$mu_v), Sigma = Sigma_uv))
  
  Sigma_xy <- matrix(c(p$sigma_x^2,
                       p$sigma_x * p$sigma_y * p$rho_xy,
                       p$sigma_x * p$sigma_y * p$rho_xy,
                       p$sigma_y^2),
                     ncol = 2)
  xy <- matrix(MASS::mvrnorm(n = p$nmax,
                             mu = c(p$alpha_x + p$beta_x * uv[1],
                                    p$alpha_y + p$beta_y * uv[2]),
                             Sigma = Sigma_xy),
               ncol = 2)
  
  x_std <- (xy[,1] - (p$alpha_x + p$beta_x * p$mu_u)) / sqrt(p$beta_x^2 * p$sigma_u^2 + p$sigma_x^2)
  y_std <- (xy[,2] - (p$alpha_y + p$beta_y * p$mu_v)) / sqrt(p$beta_y^2 * p$sigma_v^2 + p$sigma_y^2)
  
  if(p$tau_x == 1) { 
    z_x = 0 
  } else {
    z_x <- p$tau_x - cut(x_std,breaks = cuts_x, labels = FALSE, right = FALSE)
  }
  
  if(p$tau_y == 1) { 
    z_y = 0 
  } else {
    z_y <- p$tau_y - cut(y_std,breaks = cuts_y, labels = FALSE, right = FALSE)
  }
  
  Rij <- pmin(plogis(p$eta_0 + p$eta_x * xy[,1]), plogis(p$eta_0 + p$eta_y * xy[,2]))
  
  n_idx <- sapply(Rij, function(j) rbinom(n = 1, size = 1, prob = j))
  
  data.frame(
    cluster = cluster_id,
    size = sum(n_idx==1),
    j = seq_len(p$nmax),
    X = xy[, 1],
    Y = xy[, 2],
    U = uv[1], 
    V = uv[2],
    R = n_idx,
    X_sev = z_x,
    Y_sev = z_y,
    Prob_R = Rij
  )
  
}

simulate_dataset <- function(M, p, d, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cuts_x <- c(-Inf, qnorm(seq(1, p$tau_x - 1) / p$tau_x), Inf)
  cuts_y <- c(-Inf, qnorm(seq(1, p$tau_y - 1) / p$tau_y), Inf)
  
  do.call(rbind, lapply(seq_len(M), function(i) simulate_one_cluster(i, p, d, cuts_x, cuts_y)))
}

make_label <- function(x) {
  x <- format(x, scientific = FALSE, trim = TRUE)
  x <- gsub("-", "m", x)
  x <- gsub("\\.", "p", x)
  x
}

simulate_replicates <- function(M,
                                p,
                                d,
                                n_reps = 1000,
                                base_seed = 1,
                                save_dir = ".",
                                chunk_size = 500) {
  
  sim_label <- paste0(
    "rxy_", make_label(p$rho_xy),
    "_ruv_", make_label(p$rho_uv),
    "_etax_", make_label(p$eta_x),
    "_etay_", make_label(p$eta_y),
    "_M_", M
  )
  
  n_chunks <- ceiling(n_reps / chunk_size)
  out_files <- vector("character", n_chunks)
  
  for (chunk_id in seq_len(n_chunks)) {
    
    rep_start <- (chunk_id - 1) * chunk_size + 1
    rep_end   <- min(chunk_id * chunk_size, n_reps)
    reps_this_chunk <- rep_start:rep_end
    
    reps <- vector("list", length(reps_this_chunk))
    
    for (k in seq_along(reps_this_chunk)) {
      
      b <- reps_this_chunk[k]
      
      reps[[k]] <- simulate_dataset(
        M = M,
        p = p,
        d = d,
        seed = base_seed + b - 1
      )
      
      reps[[k]]$replicate <- b
    }
    
    out_file <- file.path(
      save_dir,
      paste0(
        "sim_reps_",
        sim_label,
        "_",
        chunk_id,
        ".rds"
      )
    )
    
    out <- list(
      label = sim_label,
      M = M,
      n_reps_total = n_reps,
      chunk_id = chunk_id,
      chunk_size = chunk_size,
      rep_start = rep_start,
      rep_end = rep_end,
      base_seed = base_seed,
      seeds = base_seed + reps_this_chunk - 1,
      params = p,
      disc = d,
      data = reps
    )
    
    saveRDS(out, out_file)
    out_files[chunk_id] <- out_file
    
    cat("Saved chunk", chunk_id, "of", n_chunks, "\n")
    cat("Saved:", out_file, "\n\n")
    
    rm(reps, out)
    gc()
  }
  
  invisible(out_files)
}


## ------------------------------------------------------------
## Usage
## ------------------------------------------------------------



## Global parameter vector (everything shared across functions)
params <- list(
  #-----maximum cluster size-----
  nmax = 100, 
  
  #-----cluster level latent parameters-----
  # Controls the latent mean of "Tooth Health", higher values indicate better tooth health of the average person.
  
  mu_u = 0, # mean of U 
  mu_v = 0, # mean of V
  sigma_u = 1, #marginal var of U
  sigma_v = 1, #marginal var of V
  rho_uv = 0.5, ###################### correlation control of U and V at latent level 
  
  #-----tooth level latent parameters-----
  
  alpha_x = 0, #mean shift for U
  alpha_y = 0, #mean shift for V
  
  beta_x = 1, #another knob, controls importance of U mean
  beta_y = 1, #another knob, controls importance of V mean
  
  sigma_x = 0.5, #Tooth level variability of U
  sigma_y = 0.5, #Tooth level variability of V
  rho_xy = 0.5, ######################## correlation control of tooth level outcomes
  
  #-----discretization-----
  
  tau_x = 5, #The number of categories of X
  tau_y = 5, # num. cat. of Y
  
  #-----retention-----
  
  eta_0 = 3, # baseline rate of teeth retaining. log(2.7) ~ 0.993 prob of retaining (~0.7% trauma injury assuming) 
  eta_x = 4, # the change in log-odds per one increase in X and Y
  eta_y = 4  # the change in log-odds per one increase in X and Y
)



grid <- expand.grid(
  eta_x  = c(0, 4),
  eta_y  = c(0, 4),
  rho_xy = c(0, 0.5),
  rho_uv = c(0, 0.5),
  M      = c(20, 100),
  KEEP.OUT.ATTRS = FALSE
)

if(HiperGator){
  save_dir <- "/orange/somnath.datta/NHANES/Simulation_RDS"
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  } else {
  save_dir <- "C:/Users/owvis/OneDrive - University of Florida/nHANES/Simulation_RDS"
}

reps <- 10000
chunk_size <- 250
cores <- 32  # change this yourself as needed

run_one_grid_scenario <- function(s) {
  
  p_s <- params
  
  p_s$eta_y  <- grid$eta_y[s]
  p_s$eta_x  <- grid$eta_x[s]
  p_s$rho_xy <- grid$rho_xy[s]
  p_s$rho_uv <- grid$rho_uv[s]
  
  cat("Starting scenario", s, "of", nrow(grid), "\n")
  cat("eta_y:", p_s$eta_y, "\n")
  cat("eta_x:", p_s$eta_x, "\n")
  cat("rho_xy:", p_s$rho_xy, "\n")
  cat("rho_uv:", p_s$rho_uv, "\n")
  cat("M:", grid$M[s], "\n\n")
  
  out_files_s <- simulate_replicates(
    M = grid$M[s],
    p = p_s,
    d = FALSE,
    n_reps = reps,
    base_seed = 100000 * s + 1,
    save_dir = save_dir,
    chunk_size = chunk_size
  )
  
  cat("Finished scenario", s, "of", nrow(grid), "\n\n")
  
  data.frame(
    scenario = s,
    eta_y = p_s$eta_y,
    eta_x = p_s$eta_x,
    rho_xy = p_s$rho_xy,
    rho_uv = p_s$rho_uv,
    M = grid$M[s],
    chunk_id = seq_along(out_files_s),
    out_file = out_files_s
  )
}

out_list <- parallel::mclapply(
  seq_len(nrow(grid)),
  run_one_grid_scenario,
  mc.cores = cores
)

grid_index <- do.call(rbind, out_list)

saveRDS(
  grid_index,
  file.path(save_dir, "simulation_grid_index.rds")
)

# 
# 
# # ======================================
# # 
# # ## Below is some testing code; simulation run is above.
# # 
# # ======================================
# # Simple Test
# sampleset <- simulate_dataset(1000, params, d=F)
# 
# 
# plot_size_hist <- function(p, dat) {
#   cl_size <- unique(dat[c("cluster","size")])
# 
#   breaks_vec <- seq(-0.5, 100.5, by = 1)
# 
#   hist(cl_size$size,
#        breaks = breaks_vec,
#        right = FALSE,
#        col = "gray",
#        main = "",
#        xlab = "Cluster Size",
#        ylab = "",
#        xlim = c(0, p$nmax),
#        xaxt = "n")
#   axis(1, at = seq(0, p$nmax, by = 4))
# 
#   invisible(cl_size)
# }
# 
# plot_size_hist(params, sampleset)

# simulate_replicates(
#   M = 100,
#   p = params,
#   d = FALSE, #d = discretization. going to say no to discretization initially, for now it's just a clinically meaningful marker.
#   n_reps = 1000,
#   save_dir = "C:/Users/owvis/OneDrive - University of Florida/nHANES/Simulation_RDS"
# )
# 
# 

# # Quick theoretical calculation:
# # 
# ## Global parameter vector (everything shared across functions)
# params <- list(
#   #-----maximum cluster size-----
#   nmax = 28, 
#   
#   #-----cluster level latent parameters-----
#   # Controls the latent mean of "Tooth Health", higher values indicate better tooth health of the average person.
#   
#   mu_u = 0, # mean of U 
#   mu_v = 0, # mean of V
#   sigma_u = 1, #marginal var of U
#   sigma_v = 1, #marginal var of V
#   rho_uv = 0.5, ###################### correlation control of U and V at latent level 
#   
#   #-----tooth level latent parameters-----
#   
#   alpha_x = 0, #mean shift for U
#   alpha_y = 0, #mean shift for V
#   
#   beta_x = 1, #another knob, controls importance of U mean
#   beta_y = 1, #another knob, controls importance of V mean
#   
#   sigma_x = 0.5, #Tooth level variability of U
#   sigma_y = 0.5, #Tooth level variability of V
#   rho_xy = 0.5, ######################## correlation control of tooth level outcomes
#   
#   #-----discretization-----
#   
#   tau_x = 5, #The number of categories of X
#   tau_y = 5, # num. cat. of Y
#   
#   #-----retention-----
#   
#   eta_0 = 3, # baseline rate of teeth retaining. log(2.7) ~ 0.993 prob of retaining (~0.7% trauma injury assuming) 
#   eta_x = 4, # the change in log-odds per one increase in X and Y
#   eta_y = 4  # the change in log-odds per one increase in X and Y
# )
# 
# 
# 
# 
# # simulate_replicates(
# #   M = 100,
# #   p = params,
# #   d = FALSE, #d = discretization. going to say no to discretization initially, for now it's just a clinically meaningful marker.
# #   n_reps = 1000,
# #   save_dir = "C:/Users/owvis/OneDrive - University of Florida/nHANES/Simulation_RDS"
# # )
# 
# grid <- expand.grid(
#   eta_x  = c(0, 4),
#   rho_xy = c(0, 0.5),
#   rho_uv = c(0, 0.5),
#   M      = c(20, 50),
#   KEEP.OUT.ATTRS = FALSE
# )
# 
# 
# 
# theoretical_correlation <- function(p) {
#   cov_xy <- p$beta_x * p$beta_y * p$rho_uv * p$sigma_u * p$sigma_v +
#     p$rho_xy * p$sigma_x * p$sigma_y
#   
#   var_x <- p$beta_x^2 * p$sigma_u^2 + p$sigma_x^2
#   var_y <- p$beta_y^2 * p$sigma_v^2 + p$sigma_y^2
#   
#   cov_xy / sqrt(var_x * var_y)
# }
# 
# 
# theory_grid <- grid
# 
# theory_grid$rho_true <- apply(theory_grid, 1, function(row) {
#   p <- params
#   
#   p$eta_x  <- as.numeric(row["eta_x"])
#   p$rho_xy <- as.numeric(row["rho_xy"])
#   p$rho_uv <- as.numeric(row["rho_uv"])
#   
#   theoretical_correlation(p)
# })
# 
# theory_grid <- theory_grid[order(
#   theory_grid$M,
#   theory_grid$eta_x,
#   theory_grid$rho_uv,
#   theory_grid$rho_xy
# ), ]
# 
# theory_grid
