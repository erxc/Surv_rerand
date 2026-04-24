# Required for isotonic regression
# install.packages("Iso")
library(Iso)
library(MASS)

expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(u) log(u) - log(1 - u)

#' Compute Pointwise CIs and Fixed-Width Uniform CBs
#' @param S_hat Numeric vector: Point estimates of the survival function over the grid.
#' @param eif_matrix Matrix (n x m): Empirical influence functions at the m grid points.
#' @param cross_cov_matrix Matrix (m x p): Covariance between EIFs and rerandomization covariates.
#' @param Sigma_B_inv Matrix (p x p): Inverse covariance of the balancing covariates.
#' @param kappa_c Numeric: The variance reduction multiplier.
#' @param n Integer: Sample size.
#' @param alpha Numeric: Significance level (default 0.05).
compute_inference <- function(S_hat, eif_matrix, cross_cov_matrix = NULL, Sigma_B_inv = NULL, kappa_c = NULL, n, alpha = 0.05) {
  
  m <- length(S_hat)
  z_crit <- qnorm(1 - alpha / 2)
  
  # 1. Compute the Covariance Operator Sigma(s,t)
  # -------------------------------------------------------------------
  # Unadjusted empirical covariance matrix (m x m)
  Sigma_unadj <- (t(eif_matrix) %*% eif_matrix) / n
  
  # Apply geometric deflation if ReM/SReM parameters are provided
  if (!is.null(cross_cov_matrix) && !is.null(Sigma_B_inv) && !is.null(kappa_c)) {
    deflation_factor <- 1 - kappa_c
    projection_matrix <- cross_cov_matrix %*% Sigma_B_inv %*% t(cross_cov_matrix)
    Sigma_final <- Sigma_unadj - deflation_factor * projection_matrix
  } else {
    Sigma_final <- Sigma_unadj # DML or SRS case where covariance is unchanged
  }
  
  # Ensure positive semi-definiteness for simulation due to numerical floating point errors
  Sigma_final <- (Sigma_final + t(Sigma_final)) / 2
  eigen_decomp <- eigen(Sigma_final)
  eigen_decomp$values[eigen_decomp$values < 0] <- 0
  Sigma_final <- eigen_decomp$vectors %*% diag(eigen_decomp$values) %*% t(eigen_decomp$vectors)
  
  # Extract pointwise variance
  sigma2_t <- diag(Sigma_final)
  sigma_t <- sqrt(pmax(sigma2_t, 1e-12))
  
  # 2. Logit-Transformed Pointwise Confidence Intervals
  # -------------------------------------------------------------------
  S_bound <- pmax(pmin(S_hat, 1 - 1e-8), 1e-8) # Prevent logit(0) or logit(1)
  sigma_tilde <- sigma_t / (S_bound * (1 - S_bound))
  
  margin <- z_crit * sigma_tilde / sqrt(n)
  
  # Apply expit to transformed bounds
  ci_lower <- expit(logit(S_bound) - margin)
  ci_upper <- expit(logit(S_bound) + margin)
  
  # Boundary condition enforcement per manuscript
  ci_lower[S_hat == 0] <- 0
  ci_upper[S_hat == 1] <- 1
  
  # 3. Fixed-Width Uniform Confidence Bands (Simulating Gaussian Process)
  # -------------------------------------------------------------------
  B_sim <- 2000
  # Simulate mean-zero Gaussian paths with covariance Sigma_final
  gp_paths <- mvrnorm(n = B_sim, mu = rep(0, m), Sigma = Sigma_final)
  
  # Compute the supremum of the absolute value of the paths over [0, tau]
  supremum_norms <- apply(abs(gp_paths), 1, max)
  c_hat_alpha <- quantile(supremum_norms, 1 - alpha)
  
  # Construct fixed-width bounds
  cb_lower_raw <- S_hat - c_hat_alpha / sqrt(n)
  cb_upper_raw <- S_hat + c_hat_alpha / sqrt(n)
  
  # Clamp to [0, 1]
  cb_lower_raw <- pmax(cb_lower_raw, 0)
  cb_upper_raw <- pmin(cb_upper_raw, 1)
  
  # 4. Isotonic Regression for Monotonicity
  # -------------------------------------------------------------------
  # Survival curves are monotonically decreasing, so bounds must be decreasing.
  # isoreg operates on non-decreasing sequences, so we negate, fit, and negate back.
  cb_lower_iso <- -isoreg(1:m, -cb_lower_raw)$yf
  cb_upper_iso <- -isoreg(1:m, -cb_upper_raw)$yf
  
  # -------------------------------------------------------------------
  # Return structured inference output
  # -------------------------------------------------------------------
  list(
    Pointwise_Var = sigma2_t,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    CB_Lower = cb_lower_iso,
    CB_Upper = cb_upper_iso,
    Supremum_Quantile = c_hat_alpha
  )
}
