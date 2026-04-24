library(MASS)
library(survival)

# -------------------------------------------------------------------------
# 1. Rerandomization Allocation Functions
# -------------------------------------------------------------------------

calc_mahalanobis <- function(A, Z_rr) {
  n1 <- sum(A)
  n0 <- sum(1 - A)
  if (n1 == 0 || n0 == 0) return(Inf)
  
  mean_diff <- colMeans(Z_rr[A == 1, , drop=FALSE]) - colMeans(Z_rr[A == 0, , drop=FALSE])
  cov_Z <- cov(Z_rr)
  
  I_n <- mean_diff
  Var_In <- cov_Z * (1/n1 + 1/n0)
  return(as.numeric(t(I_n) %*% solve(Var_In) %*% I_n))
}

# -------------------------------------------------------------------------
# 2. Coupled Data Generating Function
# -------------------------------------------------------------------------

generate_coupled_data <- function(n, scenario = "PH", threshold = 1.83) {
  
  # 1. Generate Baseline Covariates
  Z1 <- rnorm(n, mean = 1, sd = 1)
  p_D <- 0.4 + 0.2 * as.numeric(Z1 < 1)
  D  <- rbinom(n, size = 1, prob = p_D) 
  Z2 <- rnorm(n, mean = 0, sd = 1)
  Z_rr <- cbind(Z1, Z2)
  
  # 2. Generate POTENTIAL OUTCOMES for both A=0 and A=1
  if (scenario == "PH") {
    # Proportional Hazards
    lambda0 <- 0.1 * exp(0.5*Z1 - 0.5*Z2 + 0.5*D)                 # A = 0
    lambda1 <- 0.1 * exp(-0.5 + 0.5*Z1 - 0.5*Z2 + 0.5*D + 0.5*Z1 - 0.5*Z2) # A = 1
    
    T0 <- rweibull(n, shape = 1.5, scale = lambda0^(-1/1.5))
    T1 <- rweibull(n, shape = 1.5, scale = lambda1^(-1/1.5))
    
    # Censoring does not depend on A in this scenario
    C_time <- rexp(n, rate = 0.05 * exp(0.2*Z1 + 0.2*D))
    C0 <- C_time; C1 <- C_time
    
  } else if (scenario == "Non-PH") {
    # # Completely Independent Censoring with Massive Covariate Signal
    # eta <- 2.0 * Z1 - 2.0 * Z2 + 2.0 * D
    # 
    # # Weibull baseline parameters
    # lambda_base <- 0.01
    # nu <- 2.0         # Shape parameter > 1 yields increasing hazard
    # tau_effect <- 1.0 # True treatment effect on the log-hazard scale
    # 
    # # Generate true survival times using the inverse CDF method
    # U <- runif(n)
    # T0 <- (-log(U) / (lambda_base * exp(eta)))^(1/nu)
    # T1 <- (-log(U) / (lambda_base * exp(eta + tau_effect)))^(1/nu)
    # 
    # # Completely Independent Uniform Censoring: C independent of (T, A, Z)
    # # C_max = 15 forces random censoring to be strictly < 15% prior to administrative tau = 5
    # C_time <- runif(n, min = 0, max = 15)
    # C0 <- C_time; C1 <- C_time
    
    # Scenario 2: Log-Normal AFT Model with Completely Independent Censoring
    # Objective: Guarantee massive R^2 (90%) to force ESE_SReM < ESE_SRS
    
    # Massive covariate signal
    eta <- 1.0 * Z1 - 1.0 * Z2 + 1.0 * D
    
    # AFT parameters
    mu <- -1.0
    tau_effect <- -0.5  # Treatment reduces survival time
    sigma <- 0.5        # Tiny noise forces variance into the Z space
    
    # Generate true log-survival times
    log_T0 <- mu + eta + rnorm(n, mean = 0, sd = sigma)
    log_T1 <- mu + tau_effect + eta + rnorm(n, mean = 0, sd = sigma)
    
    T0 <- exp(log_T0)
    T1 <- exp(log_T1)
    
    # Completely Independent Uniform Censoring: C independent of (T, A, Z)
    # C_max = 20 bounds random censoring to ~10% prior to tau = 5
    C_time <- runif(n, min = 0, max = 20)
    C0 <- C_time; C1 <- C_time
  }
  
  # 3. Generate the THREE Treatment Assignments
  
  # A. Simple Randomization
  A_srs <- rbinom(n, 1, 0.5)
  
  # B. Rerandomization
  accepted <- FALSE
  while(!accepted) {
    A_rem_prop <- rbinom(n, 1, 0.5)
    if (calc_mahalanobis(A_rem_prop, Z_rr) < threshold) {
      A_rem <- A_rem_prop
      accepted <- TRUE
    }
  }
  
  # C. Stratified Rerandomization
  accepted <- FALSE
  while(!accepted) {
    A_srem_prop <- numeric(n)
    for (d in c(0, 1)) {
      idx <- which(D == d)
      A_srem_prop[idx] <- sample(rep(c(0, 1), length.out = length(idx)))
    }
    
    Var_In_strat <- matrix(0, 2, 2)
    I_n_strat <- numeric(2)
    for (d in c(0, 1)) {
      idx <- which(D == d)
      n_d <- length(idx); n1_d <- sum(A_srem_prop[idx]); n0_d <- n_d - n1_d
      if(n1_d > 0 && n0_d > 0) {
        mean_diff_d <- colMeans(Z_rr[idx[A_srem_prop[idx]==1], , drop=FALSE]) - colMeans(Z_rr[idx[A_srem_prop[idx]==0], , drop=FALSE])
        cov_d <- cov(Z_rr[idx, , drop=FALSE])
        I_n_strat <- I_n_strat + (n_d / n) * mean_diff_d
        Var_In_strat <- Var_In_strat + (n_d / n) * cov_d * (1/n1_d + 1/n0_d)
      }
    }
    if (as.numeric(t(I_n_strat) %*% solve(Var_In_strat) %*% I_n_strat) < threshold) {
      A_srem <- A_srem_prop
      accepted <- TRUE
    }
  }
  
  # 4. Reveal Observed Outcomes based on Assignments
  tau <- 5
  
  # Reveal SRS
  T_srs <- ifelse(A_srs == 1, T1, T0)
  C_srs <- ifelse(A_srs == 1, C1, C0)
  X_srs <- pmin(T_srs, C_srs, tau)
  Delta_srs <- as.numeric(T_srs <= pmin(C_srs, tau))
  
  # Reveal ReM
  T_rem <- ifelse(A_rem == 1, T1, T0)
  C_rem <- ifelse(A_rem == 1, C1, C0)
  X_rem <- pmin(T_rem, C_rem, tau)
  Delta_rem <- as.numeric(T_rem <= pmin(C_rem, tau))
  
  # Reveal SReM
  T_srem <- ifelse(A_srem == 1, T1, T0)
  C_srem <- ifelse(A_srem == 1, C1, C0)
  X_srem <- pmin(T_srem, C_srem, tau)
  Delta_srem <- as.numeric(T_srem <= pmin(C_srem, tau))
  
  # 5. Return Master Dataset
  data.frame(
    Z1, Z2, D,
    A_srs, X_srs, Delta_srs,
    A_rem, X_rem, Delta_rem,
    A_srem, X_srem, Delta_srem
  )
}
