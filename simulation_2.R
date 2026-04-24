# =====================================================================
# MASTER SCRIPT: SCENARIO 2 (Non-Proportional Hazards)
# Estimators: KM and DML
# =====================================================================

rm(list = ls())

library(MASS)
library(survival)
library(survSuperLearner)
library(data.table) # For fast list-to-dataframe binding

# SOURCE YOUR HELPER SCRIPTS HERE
source("data-generation.R")
source("inference.R")
source("dml.R")
source("km.R")

# ---------------------------------------------------------
# 1. Simulation Parameters
# ---------------------------------------------------------
set.seed(2026)
n <- 100            # Sample size (Run separately for n=100)
n_reps <- 1000      # Monte Carlo Replicates
t_eval <- c(1, 2, 3, 4)
t_grid <- seq(1, 4, length.out = 50)
t_grid <- sort(unique(c(t_grid, t_eval)))

# Pre-allocate a list to store results
results_list <- list()
counter <- 1

# ---------------------------------------------------------
# 2. The Monte Carlo Loop
# ---------------------------------------------------------
cat("Starting Simulation Scenario 2 (n =", n, ")...\n")

for (m in 1:n_reps) {
  if (m %% 50 == 0) cat("Processing Replicate:", m, "/", n_reps, "\n")
  
  # A. Generate Coupled Data (Outputs SRS, ReM, SReM together)
  dat <- generate_coupled_data(n = n, scenario = "Non-PH", threshold = 1.83)
  
  Z_rr <- cbind(dat$Z1, dat$Z2)
  W_all <- cbind(dat$Z1, dat$Z2, dat$D) # Covariates for censoring and DML
  
  # B. Define the designs to loop over
  designs <- c("SRS", "ReM", "SReM")
  
  # Safe execution wrapper to prevent crashing
  tryCatch({
    
    for (design in designs) {
      
      # Extract design-specific observed data
      A_obs <- dat[[paste0("A_", tolower(design))]]
      X_obs <- dat[[paste0("X_", tolower(design))]]
      Delta_obs <- dat[[paste0("Delta_", tolower(design))]]
      
      # ---------------------------------------------------
      # RUN ESTIMATOR 1: KM
      # ---------------------------------------------------
      res_km <- run_KM(
        time = X_obs, event = Delta_obs, treat = A_obs, 
        Z_rr = Z_rr
      )
      
      # ---------------------------------------------------
      # RUN ESTIMATOR 2: DML
      # ---------------------------------------------------
      res_dml <- run_DML(
        time = X_obs, event = Delta_obs, treat = A_obs, 
        W = W_all
      )
      
      # ---------------------------------------------------
      # Package Results for this Replicate/Design
      # ---------------------------------------------------
      for (k in 1:length(t_eval)) {
        
        # KM Row
        results_list[[counter]] <- data.table(
          Replicate = m,
          n = n,
          Design = design,
          Estimator = "KM",
          Time = t_eval[k],
          Estimate = res_km$Estimate[k],
          Var_Estimated = res_km$Var_Deflated[k],
          CI_Lower = res_km$CI_Lower[k],
          CI_Upper = res_km$CI_Upper[k],
          Uniform_CB_Lower = res_km$Uniform_CB_Lower[which(t_grid == t_eval[k])],
          Uniform_CB_Upper = res_km$Uniform_CB_Upper[which(t_grid == t_eval[k])]
        )
        counter <- counter + 1
        
        # DML Row
        results_list[[counter]] <- data.table(
          Replicate = m,
          n = n,
          Design = design,
          Estimator = "DML",
          Time = t_eval[k],
          Estimate = res_dml$Estimate[k],
          Var_Estimated = res_dml$Var_True[k],
          CI_Lower = res_dml$CI_Lower[k],
          CI_Upper = res_dml$CI_Upper[k],
          Uniform_CB_Lower = res_dml$Uniform_CB_Lower[which(t_grid == t_eval[k])], # If DML supports grid
          Uniform_CB_Upper = res_dml$Uniform_CB_Upper[which(t_grid == t_eval[k])]
        )
        counter <- counter + 1
      }
    }
    
  }, error = function(e) {
    cat("Error in replicate", m, ":", conditionMessage(e), "\n")
    # Skip to next replicate on failure, leaving out this m from the list
  })
}

# ---------------------------------------------------------
# 3. Aggregate and Save to CSV
# ---------------------------------------------------------
cat("Compiling results...\n")
final_results <- rbindlist(results_list)

file_name <- paste0("Simulation_Results_Scenario2_AFT_20260422_n", n, ".csv")
write.csv(final_results, file_name, row.names = FALSE)
cat("Successfully saved to", file_name, "\n")
