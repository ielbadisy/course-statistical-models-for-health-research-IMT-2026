############################################################
# PROJECT SOLUTION
############################################################

pacman::p_load(survival, torch, survdnn, ggplot2, dplyr, table1, unsurv)

torch::torch_set_num_threads(4); torch::torch_manual_seed(1)


melanoma2 <- read.csv("datasets/melanoma2.csv") 


melanoma2 <- melanoma2 |> mutate(OS = as.numeric(OS), Event_of_OS = as.integer(Event_of_OS)) |>
  mutate(across(c(Sex, Histology, Comorbidity, PS,
                  Lung_Metastasis, GG_metastasis, Bone_metastasis, Other,
                  Therapeutic_line, Total_doses, Response_type), as.factor))



# 1) SurvDNN (AFT)
form <- Surv(OS, Event_of_OS) ~ Sex + Age + Histology + Comorbidity + PS +
  Lung_Metastasis + GG_metastasis + Bone_metastasis + Other +
  Therapeutic_line + Total_doses

mod <- survdnn(
  form, melanoma2,
  hidden = c(32, 64, 16), epochs = 500, loss = "aft", lr = 1e-3,
  activation = "gelu", dropout = 0.10, verbose = TRUE, .seed = 1
  )

print(summary(mod))

# 2) Predict survival curves on common grid
tmax  <- floor(max(melanoma2$OS, na.rm = TRUE))
times <- seq(0, tmax, by = 1)
S <- as.matrix(predict(mod, melanoma2, type = "survival", times = times))

# 3) Clustering (unsurv API)
set.seed(123)
fit <- unsurv(S, times, K = 3, K_max = 6)
print(fit)

cluster_vec <- fit$clusters

melanoma2$cluster <- factor(cluster_vec)

table(melanoma2$cluster)

# 4) Table 1
table1(
  ~ Sex + Age + Histology + Comorbidity + PS +
    Lung_Metastasis + GG_metastasis + Bone_metastasis + Other +
    Therapeutic_line + Total_doses + Response_type | cluster,
  data = melanoma2
  )

# 5) Dynamic RMST + tau*
# RMST_i(tau) = ∫_0^tau S_i(t) dt  (area under the predicted survival curve)
rmst_dynamic <- function(S_mat, time) {
  stopifnot(is.matrix(S_mat), length(time) == ncol(S_mat))
  stopifnot(all(diff(time) > 0))
  
  n <- nrow(S_mat)
  m <- ncol(S_mat)
  
  # Step 1) widths of each interval [t_j, t_{j+1}]
  dt <- diff(time)                         # length m-1
  
  # Step 2) trapezoid height per interval: (S(t_j) + S(t_{j+1})) / 2
  S_left  <- S_mat[, 1:(m - 1), drop = FALSE]
  S_right <- S_mat[, 2:m,       drop = FALSE]
  S_mid   <- (S_left + S_right) / 2        # n x (m-1)
  
  # Step 3) area increment per interval: height * width
  area_incr <- sweep(S_mid, 2, dt, `*`)    # multiply each column j by dt[j]
  
  # Step 4) cumulative sum of increments gives RMST at each grid point
  RMST_no0 <- t(apply(area_incr, 1, cumsum))  # n x (m-1)
  
  # Step 5) RMST at tau = 0 is 0 → prepend a zero column to align with `time`
  RMST <- cbind(0, RMST_no0)               # n x m
  colnames(RMST) <- time
  
  RMST
}

RMST_mat <- rmst_dynamic(S, times)

# Choose clinical horizon tau* (e.g., 24 months), but cap at observed max grid
tau_star <- min(24, max(times))

# Find exact grid match; if not present, use nearest grid point (explicit)
k_tau <- which(times == tau_star)
if (length(k_tau) != 1) {
  k_tau <- which.min(abs(times - tau_star))
  tau_star <- times[k_tau]
}

# Individual RMST at tau*
melanoma2$RMST_tau <- RMST_mat[, k_tau]

aggregate(RMST_tau ~ cluster, data = melanoma2, FUN = mean)

# 7) Plots (mean + individual)
# Build a long data.frame: one row per (patient i, horizon tau_j)
n <- nrow(melanoma2)
m <- length(times)

rmst_long <- data.frame(
  id      = rep(seq_len(n), times = m),
  tau     = rep(times, each = n),
  rmst    = as.vector(RMST_mat),
  cluster = rep(melanoma2$cluster, times = m)
)

# Mean dynamic RMST per cluster (one mean curve per facet)
print(
  ggplot(rmst_long, aes(x = tau, y = rmst, group = cluster)) +
    stat_summary(fun = mean, geom = "line") +
    facet_wrap(~ cluster) +
    labs(x = "Tau", y = "RMST(tau)", title = "Mean dynamic RMST by cluster") +
    theme_minimal()
  )

# Individual RMST trajectories (many faint lines per facet)
print(
  ggplot(rmst_long, aes(x = tau, y = rmst, group = id)) +
    geom_line(alpha = 0.15) +
    facet_wrap(~ cluster) +
    labs(x = "Tau", y = "RMST(tau)", title = "Individual dynamic RMST curves by cluster") +
    theme_minimal()
)
############################################################
# END
############################################################
