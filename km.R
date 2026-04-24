library(survival)

#' @param design Character string indicating the randomization design ("SRS", "ReM", or "SReM")

#' Unadjusted Kaplan-Meier Estimator for Rerandomization
run_KM <- function(time, event, treat, Z_rr, t_eval = c(1, 2, 3, 4), t_grid = seq(1, 4, length.out = 50), c_threshold = 1.83) {
  
  n <- length(treat)
  t_grid <- sort(unique(c(t_grid, t_eval)))
  pi_1 <- 0.5
  
  # Isolate A = 1
  idx1 <- which(treat == 1)
  time1 <- time[idx1]
  event1 <- event[idx1]
  
  # 1. Fit Kaplan-Meier for Outcome T and Censoring C
  fit_T <- survfit(Surv(time1, event1) ~ 1)
  fit_C <- survfit(Surv(time1, 1 - event1) ~ 1)
  
  S_T_eval <- summary(fit_T, times = t_eval, extend = TRUE)$surv
  S_T_grid <- summary(fit_T, times = t_grid, extend = TRUE)$surv
  
  # 2. EIF Extraction Function (Martingale Integral)
  compute_km_eif <- function(t_target, S_target) {
    eif_full <- rep(0, n)
    
    jumps <- fit_T$time[fit_T$n.event > 0 & fit_T$time <= t_target]
    if (length(jumps) == 0) return(eif_full)
    
    # Evaluate at s-
    S_T_minus <- summary(fit_T, times = jumps - 1e-8, extend = TRUE)$surv
    S_C_minus <- summary(fit_C, times = jumps - 1e-8, extend = TRUE)$surv
    
    d_Lambda <- summary(fit_T, times = jumps)$n.event / summary(fit_T, times = jumps)$n.risk
    
    denom <- S_T_minus * S_C_minus
    denom[denom == 0] <- 1e-8
    
    eif_1 <- sapply(1:length(idx1), function(i) {
      x_i <- time1[i]; d_i <- event1[i]
      part_N <- sum((jumps == x_i & d_i == 1) / denom[jumps == x_i])
      part_Y <- sum((d_Lambda / denom)[jumps <= x_i])
      return(-S_target * (part_N - part_Y))
    })
    
    eif_full[idx1] <- (1 / pi_1) * eif_1
    return(eif_full)
  }
  
  # Compute EIF matrices
  eif_matrix_eval <- sapply(1:length(t_eval), function(k) compute_km_eif(t_eval[k], S_T_eval[k]))
  eif_matrix_grid <- sapply(1:length(t_grid), function(k) compute_km_eif(t_grid[k], S_T_grid[k]))
  
  # 3. Geometric Projection Mechanics & Route to Universal Inference Engine
  if (design == "SRS") {
    
    # Under Simple Random Sampling, the deflation factor is mathematically zero.
    inf_eval <- compute_inference(S_hat = S_T_eval, eif_matrix = eif_matrix_eval, n = n)
    inf_grid <- compute_inference(S_hat = S_T_grid, eif_matrix = eif_matrix_grid, n = n)
    
  } else if (design %in% c("ReM", "SReM")) {
    
    # Under Rerandomization architectures, compute the explicit geometric projection
    Z_centered <- scale(Z_rr, center = TRUE, scale = FALSE)
    Sigma_B <- cov(Z_rr) / (pi_1 * (1 - pi_1))
    Sigma_B_inv <- solve(Sigma_B)
    
    # kappa(c) limits the deflation based on the Mahalanobis acceptance threshold
    kappa_c <- pchisq(c_threshold, df = ncol(Z_rr) + 2) / pchisq(c_threshold, df = ncol(Z_rr))
    
    # STRICT CORRECTION: Divide by pi_1 to recover the potential outcome expectation
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