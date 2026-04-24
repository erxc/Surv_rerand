library(survival)

#' @param design Character string indicating the randomization design ("SRS", "ReM", or "SReM")

#' IPCW Kaplan-Meier Estimator for Rerandomization
run_IPCW <- function(time, event, treat, W_cens, Z_rr, t_eval = c(1, 2, 3, 4), t_grid = seq(1, 4, length.out = 50), c_threshold = 1.83) {
  
  n <- length(treat)
  t_grid <- sort(unique(c(t_grid, t_eval)))
  pi_1 <- 0.5
  
  idx1 <- which(treat == 1)
  time1 <- time[idx1]
  event1 <- event[idx1]
  W1 <- W_cens[idx1, , drop=FALSE]
  
  # 1. Fit Cox Model for Censoring on A = 1
  # Surv(time, 1 - event) to model the censoring hazard
  cox_C <- coxph(Surv(time1, 1 - event1) ~ ., data = as.data.frame(W1))
  basehaz_C <- basehaz(cox_C, centered = FALSE)
  
  # Helper to get S_C(t | W_i)
  get_S_C <- function(t_target) {
    bh <- basehaz_C$hazard[max(which(basehaz_C$time <= t_target + 1e-8))]
    if (length(bh) == 0) bh <- 0
    exp(-bh * exp(as.matrix(W1) %*% coef(cox_C)))
  }
  
  # 2. Compute the IPCW-KM Point Estimates
  # Since the weights depend on time, we must compute the weighted Nelson-Aalen jumps
  jumps_T <- sort(unique(time1[event1 == 1]))
  dLambda_ipcw <- numeric(length(jumps_T))
  
  for (j in 1:length(jumps_T)) {
    s_jump <- jumps_T[j]
    # Evaluate S_C at s-
    S_C_minus <- get_S_C(s_jump - 1e-8)
    weights <- 1 / pmax(S_C_minus, 1e-4) # Stabilize extreme weights
    
    # Weighted at-risk and event counts
    Y_weighted <- sum(weights[time1 >= s_jump])
    dN_weighted <- sum(weights[time1 == s_jump & event1 == 1])
    
    dLambda_ipcw[j] <- dN_weighted / pmax(Y_weighted, 1e-8)
  }
  
  # Product integral for survival
  Lambda_ipcw_step <- stepfun(jumps_T, c(0, cumsum(dLambda_ipcw)))
  S_T_ipcw_func <- function(t) exp(-Lambda_ipcw_step(t))
  
  S_T_eval <- sapply(t_eval, S_T_ipcw_func)
  S_T_grid <- sapply(t_grid, S_T_ipcw_func)
  
  # 3. EIF Extraction Function (Weighted Martingale)
  compute_ipcw_eif <- function(t_target, S_target) {
    eif_full <- rep(0, n)
    valid_jumps <- jumps_T[jumps_T <= t_target]
    if (length(valid_jumps) == 0) return(eif_full)
    
    eif_1 <- numeric(length(idx1))
    
    for (j in 1:length(valid_jumps)) {
      s_jump <- valid_jumps[j]
      S_C_minus <- get_S_C(s_jump - 1e-8)
      weights <- 1 / pmax(S_C_minus, 1e-4)
      
      Y_weighted <- sum(weights[time1 >= s_jump])
      S_0_weighted <- Y_weighted / length(idx1) 
      
      # Patient specific martingale increments at s_jump
      dN_i <- as.numeric(time1 == s_jump & event1 == 1)
      Y_i <- as.numeric(time1 >= s_jump)
      dLambda_s <- dLambda_ipcw[j]
      
      # The weighted contribution to the EIF integral
      increment <- (weights * dN_i / S_0_weighted) - (weights * Y_i * dLambda_s / S_0_weighted)
      eif_1 <- eif_1 + increment
    }
    
    eif_full[idx1] <- (1 / pi_1) * (-S_target * eif_1)
    return(eif_full)
  }
  
  # Compute EIF matrices
  eif_matrix_eval <- sapply(1:length(t_eval), function(k) compute_ipcw_eif(t_eval[k], S_T_eval[k]))
  eif_matrix_grid <- sapply(1:length(t_grid), function(k) compute_ipcw_eif(t_grid[k], S_T_grid[k]))
  
  # 4. Geometric Projection Mechanics
  if (design == "SRS") {

    # Under Simple Random Sampling, the deflation factor is mathematically zero.
    # Pass NULL to the projection arguments to return standard unadjusted variance.
    inf_eval <- compute_inference(S_hat = S_T_eval, eif_matrix = eif_matrix_eval, n = n)
    inf_grid <- compute_inference(S_hat = S_T_grid, eif_matrix = eif_matrix_grid, n = n)

  } else if (design %in% c("ReM", "SReM")) {

    # Under Rerandomization architectures, compute the explicit geometric projection
    Z_centered <- scale(Z_rr, center = TRUE, scale = FALSE)
    Sigma_B <- cov(Z_rr) / (pi_1 * (1 - pi_1))
    Sigma_B_inv <- solve(Sigma_B)

    # kappa(c) limits the deflation based on the Mahalanobis acceptance threshold
    kappa_c <- pchisq(c_threshold, df = ncol(Z_rr) + 2) / pchisq(c_threshold, df = ncol(Z_rr))

    # cross_cov_eval <- t(eif_matrix_eval) %*% Z_centered / n
    # cross_cov_grid <- t(eif_matrix_grid) %*% Z_centered / n
    
    cross_cov_eval <- (t(eif_matrix_eval) %*% Z_centered / n) / pi_1
    cross_cov_grid <- (t(eif_matrix_grid) %*% Z_centered / n) / pi_1

    # Feed projection components to the inference engine
    inf_eval <- compute_inference(S_hat = S_T_eval, eif_matrix = eif_matrix_eval,
                                  cross_cov_matrix = cross_cov_eval, Sigma_B_inv = Sigma_B_inv, kappa_c = kappa_c, n = n)
    inf_grid <- compute_inference(S_hat = S_T_grid, eif_matrix = eif_matrix_grid,
                                  cross_cov_matrix = cross_cov_grid, Sigma_B_inv = Sigma_B_inv, kappa_c = kappa_c, n = n)
  } else {
    stop("Invalid design specified.")
  }
  
  list(
    Estimate = S_T_eval,
    Var_Deflated = inf_eval$Pointwise_Var,
    CI_Lower = inf_eval$CI_Lower,
    CI_Upper = inf_eval$CI_Upper,
    Uniform_CB_Lower = inf_grid$CB_Lower,
    Uniform_CB_Upper = inf_grid$CB_Upper,
    Uniform_Critical_Value = inf_grid$Supremum_Quantile
  )
}
