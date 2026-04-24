library(survival)
library(survSuperLearner)
# library(Iso)  # Assumes your compute_inference engine is loaded
# library(MASS)

#' Cross-Fitted DML Estimator
#' Implements Westling Theorem 1 with strictly V=5 and known pi=0.5
run_DML <- function(time, event, treat, W, 
                    t_eval = c(1, 2, 3, 4), 
                    t_grid = seq(1, 4, length.out = 50),
                    event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
                    # event.SL.library = c("survSL.km", "survSL.coxph"),
                    cens.SL.library = c("survSL.km", "survSL.coxph")) {
  
  n <- length(time)
  t_grid <- sort(unique(c(t_grid, t_eval)))
  W_df <- as.data.frame(W)
  
  V <- 5           # Strictly 5-fold cross-fitting
  pi_1 <- 0.5      # True randomized propensity 
  
  # The integration grid must contain all observed event times AND evaluation targets
  jumps <- sort(unique(time[event == 1]))
  if (jumps[1] > 0) jumps <- c(0, jumps)
  eval_times <- sort(unique(c(jumps, t_grid, t_eval)))
  
  phi_matrix_grid <- matrix(0, nrow = n, ncol = length(t_grid))
  folds <- sample(rep(1:V, length.out = n))
  
  # ---------------------------------------------------------
  # 1. The 5-Fold Cross-Fitting Loop
  # ---------------------------------------------------------
  for (v in 1:V) {
    idx_train <- which(folds != v)
    idx_test  <- which(folds == v)
    
    # Isolate training data strictly to A = 1
    train_A1 <- idx_train[treat[idx_train] == 1]
    
    time_train <- time[train_A1]
    event_train <- event[train_A1]
    W_train <- W_df[train_A1, , drop = FALSE]
    W_test <- W_df[idx_test, , drop = FALSE]
    
    # ---------------------------------------------------------
    # A. Execute survSuperLearner (The proper way)
    # ---------------------------------------------------------
    # Fits both S_0 and G_0 simultaneously via the iterative algorithm
    fit_SL <- survSuperLearner(
      time = time_train, 
      event = event_train, 
      X = W_train, 
      newX = W_test, 
      new.times = eval_times, 
      event.SL.library = event.SL.library,
      cens.SL.library = cens.SL.library
    )
    
    # Extract prediction matrices (rows = patients, cols = eval_times)
    pred_S_matrix <- fit_SL$event.SL.predict
    pred_G_matrix <- fit_SL$cens.SL.predict
    
    # ---------------------------------------------------------
    # B. Extract EIF on Validation Fold (Martingale Integral)
    # ---------------------------------------------------------
    for (i_local in 1:length(idx_test)) {
      i <- idx_test[i_local]
      y_i <- time[i]
      d_i <- event[i]
      a_i <- treat[i]
      
      S_T <- pred_S_matrix[i_local, ]
      G_C <- pred_G_matrix[i_local, ]
      
      for (k in 1:length(t_grid)) {
        t_target <- t_grid[k]
        
        idx_t <- match(t_target, eval_times)
        S_T_t <- S_T[idx_t]
        
        # Bypass augmentation if A_i == 0
        if (a_i == 0) {
          phi_matrix_grid[i, k] <- S_T_t
          next
        }
        
        limit <- min(t_target, y_i)
        valid_jumps_idx <- which(eval_times <= limit & eval_times %in% jumps)
        
        # 1. Event Term
        event_term <- 0
        if (y_i <= t_target && d_i == 1) {
          idx_Y <- match(y_i, eval_times)
          # Strictly enforce left-limits S(X-|Z)
          S_T_Y_minus <- ifelse(idx_Y == 1, 1, S_T[idx_Y - 1])
          G_C_Y_minus <- ifelse(idx_Y == 1, 1, G_C[idx_Y - 1])
          denom_Y <- max(S_T_Y_minus * G_C_Y_minus, 1e-8)
          event_term <- 1 / denom_Y
        }
        
        # 2. Integral Term
        integral_term <- 0
        if (length(valid_jumps_idx) > 0) {
          for (j in valid_jumps_idx) {
            S_C_u_minus <- ifelse(j == 1, 1, G_C[j - 1])
            S_T_u_minus <- ifelse(j == 1, 1, S_T[j - 1])
            
            dLambda_u <- 1 - (S_T[j] / max(S_T_u_minus, 1e-8))
            
            # Strictly enforce left-limits S(s-|Z) for the compensator denominator
            denom_u <- max(S_T_u_minus * S_C_u_minus, 1e-8)
            
            integral_term <- integral_term + (dLambda_u / denom_u)
          }
        }
        
        # 3. Theorem 1 Combination
        bracket <- 1 - (1 / pi_1) * (event_term - integral_term)
        phi_matrix_grid[i, k] <- S_T_t * bracket
      }
    }
  }
  
  # ---------------------------------------------------------
  # 2. Aggregation & Formatting
  # ---------------------------------------------------------
  S_hat_grid <- colMeans(phi_matrix_grid)
  eif_matrix_grid <- sweep(phi_matrix_grid, 2, S_hat_grid, FUN = "-")
  
  eval_idx <- match(t_eval, t_grid)
  S_hat_eval <- S_hat_grid[eval_idx]
  eif_matrix_eval <- eif_matrix_grid[, eval_idx, drop = FALSE]
  
  # Route to Universal Inference Engine (Deflation is 0 for DML)
  inf_eval <- compute_inference(S_hat = S_hat_eval, eif_matrix = eif_matrix_eval, n = n)
  inf_grid <- compute_inference(S_hat = S_hat_grid, eif_matrix = eif_matrix_grid, n = n)
  
  list(
    Estimate = S_hat_eval,
    Var_True = inf_eval$Pointwise_Var,
    CI_Lower = inf_eval$CI_Lower,
    CI_Upper = inf_eval$CI_Upper,
    Uniform_CB_Lower = inf_grid$CB_Lower,
    Uniform_CB_Upper = inf_grid$CB_Upper,
    Uniform_Critical_Value = inf_grid$Supremum_Quantile
  )
}