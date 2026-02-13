
library(mcstatsim)

# Define a simple simulation function
sim_function <- function(a, b) {
  Sys.sleep(0.2)  # Simulate a time-consuming process
  return(data.frame(result = a + b))
}

# Generate a grid of parameters
params <- expand.grid(a = 1:3, b = 4:6)

# Run simulations
results <- runsim(n = 5, grid_params = params, sim_func = sim_function)


## Working example (more detailed)

# (1) Generate fully observed data -> `data_complete`

# (2) estimate the beta coefficients values from `data_complete`

# (3) Introduce missigness under MCAR to complete dataset generated in (1) -> `data_missing`

# (4) impute the dataset generated at (3) -> `data_imputed`

# (5) Use the following simulation targets to compute the distortion between beta from `data_complete` and beta `data_imputed`: 
  
# ============================================================
# FULL LAB CODE — DO NOT MODIFY
# ============================================================

pacman::p_load(
  mcstatsim, survival, dplyr, ggplot2, tidyr,
  VIM, simputation, missForest, missRanger, missCforest
)

set.seed(123)

# ------------------------------------------------------------
# Data generator
# ------------------------------------------------------------

gencox <- function(n = 300, maxTime = 7, logHR = 0.5) {
  
  lambda <- 0.1
  rho <- 1.6
  rateC <- 0.09
  
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5)
  
  U <- runif(n)
  
  Tlat <- (-log(U) /
             (lambda * exp(logHR * (x1 + x2 + x3))))^(1 / rho)
  
  Ctimes <- rexp(n, rate = rateC)
  
  time <- pmin(Tlat, Ctimes)
  status <- as.numeric(Tlat <= Ctimes)
  
  time <- ifelse(time > maxTime, maxTime, time)
  status <- ifelse(time >= maxTime, 1, status)
  
  data <- data.frame(time, status, x1, x2, x3)
  
  # IMPORTANT: stable factor levels (prevents imputation + model issues)
  data$x3 <- factor(data$x3, levels = c(0, 1))
  
  data
}

# ------------------------------------------------------------
# Reference Cox estimate
# ------------------------------------------------------------

estimate_coxest <- function(data) {
  
  mod <- survival::Surv(time, status) ~ x1 + x2 + x3
  
  coefs <- summary(
    survival::coxph(mod, data = data)
  )$coef
  
  coefs[, 1]
}

# ------------------------------------------------------------
# Introduce MCAR missingness
# ------------------------------------------------------------

introduce_MCAR <- function(x, covariates, p = 0.3) {
  
  stopifnot(is.data.frame(x), p >= 0, p <= 1)
  stopifnot(all(covariates %in% names(x)))
  
  x[covariates] <- lapply(
    x[covariates],
    function(z) {
      idx <- sample.int(length(z), floor(p * length(z)))
      z[idx] <- NA
      z
    }
  )
  
  x
}

# ------------------------------------------------------------
# Imputation methods (patched to avoid missCforest warnings)
# ------------------------------------------------------------

imputer <- function(data, method) {
  
  stopifnot(is.data.frame(data))
  method <- as.character(method)
  
  supported_methods <- c(
    "knn", "cart", "missforest", "missranger", "misscforest", "complete"
  )
  stopifnot(method %in% supported_methods)
  
  # keep a copy of original factor levels (important!)
  x3_levels <- NULL
  if ("x3" %in% names(data) && is.factor(data$x3)) {
    x3_levels <- levels(data$x3)
  }
  
  out <- switch(
    
    method,
    
    knn = {
      VIM::kNN(data)[names(data)]
    },
    
    cart = {
      simputation::impute_cart(data, . ~ .)
    },
    
    missforest = {
      missForest::missForest(
        data,
        verbose = FALSE
      )$ximp
    },
    
    missranger = {
      missRanger::missRanger(
        data,
        pmm.k = 5,
        num.trees = 100,
        verbose = 0
      )
    },
    
    misscforest = {
      
      # ---- KEY PATCH ----
      # missCforest / cforest can return multi-column predictions for factors,
      # producing: "number of items to replace..." warnings.
      # Solution: temporarily convert factors to character (or numeric),
      # impute, then restore factors.
      
      data2 <- data
      
      # convert all factors to character for imputation stability
      fac_cols <- names(which(sapply(data2, is.factor)))
      if (length(fac_cols) > 0) {
        data2[fac_cols] <- lapply(data2[fac_cols], as.character)
      }
      
      imp <- suppressWarnings(
        missCforest::missCforest(data2)
      )
      
      # restore factor columns
      if (length(fac_cols) > 0) {
        imp[fac_cols] <- lapply(imp[fac_cols], as.factor)
      }
      
      imp
    },
    
    complete = {
      data[complete.cases(data), ]
    }
  )
  
  # Re-stabilize x3 factor levels after imputation / filtering
  if ("x3" %in% names(out)) {
    if (!is.factor(out$x3)) out$x3 <- as.factor(out$x3)
    if (!is.null(x3_levels)) out$x3 <- factor(out$x3, levels = x3_levels)
  }
  
  out
}

# ------------------------------------------------------------
# Evaluation
# ------------------------------------------------------------

evaluate_coxest <- function(data, truelogHR) {
  
  mod <- survival::Surv(time, status) ~ x1 + x2 + x3
  
  fit <- survival::coxph(mod, data = data)
  
  coefs <- summary(fit)$coef
  
  estimates <- coefs[, 1]
  se <- coefs[, 3]
  
  ci_lower <- estimates - 1.96 * se
  ci_upper <- estimates + 1.96 * se
  
  data.frame(
    estimates = estimates,
    bias = mcstatsim::calc_bias(estimates, truelogHR)$bias,
    coverage = mcstatsim::calc_coverage(ci_lower, ci_upper, truelogHR)$coverage,
    rmse = mcstatsim::calc_rmse(estimates, truelogHR)$rmse
  )
}

# ------------------------------------------------------------
# One simulation replicate
# ------------------------------------------------------------

simcox <- function(n, logHR, pmiss, method) {
  
  data_complete <- gencox(
    n = n,
    logHR = logHR
  )
  
  # reference coefficients from complete data
  truelogHR <- estimate_coxest(data_complete)
  
  # induce missingness (MCAR) in x2
  data_missing <- introduce_MCAR(
    data_complete,
    covariates = "x2",
    p = pmiss
  )
  
  # impute / complete cases
  data_imputed <- imputer(
    data_missing,
    method
  )
  
  # evaluate Cox on the imputed/complete dataset
  res <- evaluate_coxest(
    data_imputed,
    truelogHR = truelogHR
  )
  
  cbind(
    n = n,
    pmiss = pmiss,
    method = method,
    res,
    row.names = NULL
  )
}

# ------------------------------------------------------------
# Simulation grid
# ------------------------------------------------------------

params <- expand.grid(
  n = c(200, 500),
  logHR = 0.5,
  pmiss = c(0.2, 0.5),
  method = c(
    "knn",
    "cart",
    "missforest",
    "missranger",
    "misscforest",
    "complete"
  )
)

# ------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------

sim_res <- mcstatsim::runsim(
  n = 10,
  grid_params = params,
  sim_func = simcox,
  show_progress = TRUE,
  num_cores = 4
)

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------

sim_res$bias <- abs(sim_res$bias)

sim_res2 <- tidyr::gather(
  sim_res,
  metric,
  value,
  c(bias, rmse)
)

ggplot(sim_res2, aes(value, method)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_grid(pmiss ~ metric, scales = "free") +
  theme_bw()
