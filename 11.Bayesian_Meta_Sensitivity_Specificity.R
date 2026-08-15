# ============================================================
# Bayesian bivariate diagnostic meta-analysis
# Outcomes: Sensitivity and Specificity only
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------
required_packages <- c("rjags", "coda", "ggplot2")

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

library(rjags)
library(coda)
library(ggplot2)

# ------------------------------------------------------------
# 2. Paths and colour palette
# ------------------------------------------------------------
input_file <- "C:/Users/30399/Desktop/LLM/main.txt"

output_dir <- "C:/Users/30399/Desktop/LLM/11.Bayesian_Meta_Sensitivity_Specificity"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------
dat_raw <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  fileEncoding = "GB18030",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Fixed columns:
# Column 1  = First author
# Column 2  = Year
# Column 12 = TP
# Column 13 = FP
# Column 14 = FN
# Column 15 = TN

dat <- data.frame(
  Study = paste0(
    trimws(as.character(dat_raw[[1]])),
    " (",
    trimws(as.character(dat_raw[[2]])),
    ")"
  ),
  TP = suppressWarnings(as.numeric(dat_raw[[12]])),
  FP = suppressWarnings(as.numeric(dat_raw[[13]])),
  FN = suppressWarnings(as.numeric(dat_raw[[14]])),
  TN = suppressWarnings(as.numeric(dat_raw[[15]]))
)

dat <- dat[
  complete.cases(dat[, c("TP", "FP", "FN", "TN")]) &
    dat$TP >= 0 &
    dat$FP >= 0 &
    dat$FN >= 0 &
    dat$TN >= 0,
]

dat$Disease_N <- dat$TP + dat$FN
dat$NonDisease_N <- dat$FP + dat$TN
dat$Total_N <- dat$TP + dat$FP + dat$FN + dat$TN

dat <- dat[
  dat$Disease_N > 0 &
    dat$NonDisease_N > 0,
]

rownames(dat) <- NULL

if (nrow(dat) < 2) {
  stop("Fewer than two eligible studies were identified.")
}

write.csv(
  dat,
  file.path(output_dir, "Bayesian_Analysis_Input_Data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 4. Observed study-level diagnostic estimates
# ------------------------------------------------------------
dat$Sensitivity <- dat$TP / dat$Disease_N
dat$Specificity <- dat$TN / dat$NonDisease_N
dat$False_Positive_Rate <- 1 - dat$Specificity

# ------------------------------------------------------------
# 5. JAGS Bayesian bivariate random-effects model
# ------------------------------------------------------------
# theta[i,1] = logit(Sensitivity)
# theta[i,2] = logit(Specificity)

jags_model <- "
model {

  for (i in 1:N) {

    TP[i] ~ dbin(se[i], diseased_n[i])
    TN[i] ~ dbin(sp[i], nondiseased_n[i])

    logit(se[i]) <- theta[i, 1]
    logit(sp[i]) <- theta[i, 2]

    theta[i, 1:2] ~ dmnorm(mu[], Omega[,])
  }

  # Vague priors for pooled logit sensitivity and specificity
  mu[1] ~ dnorm(0, 0.001)
  mu[2] ~ dnorm(0, 0.001)

  # Between-study heterogeneity
  sd_sensitivity ~ dunif(0, 5)
  sd_specificity ~ dunif(0, 5)

  # Correlation between logit sensitivity and logit specificity
  rho ~ dunif(-0.95, 0.95)

  # Precision matrix
  Omega[1, 1] <- 1 / (pow(sd_sensitivity, 2) * (1 - pow(rho, 2)))
  Omega[2, 2] <- 1 / (pow(sd_specificity, 2) * (1 - pow(rho, 2)))
  Omega[1, 2] <- -rho / (sd_sensitivity * sd_specificity * (1 - pow(rho, 2)))
  Omega[2, 1] <- Omega[1, 2]

  # Pooled diagnostic performance
  pooled_sensitivity <- ilogit(mu[1])
  pooled_specificity <- ilogit(mu[2])

}
"

jags_data <- list(
  N = nrow(dat),
  TP = dat$TP,
  TN = dat$TN,
  diseased_n = dat$Disease_N,
  nondiseased_n = dat$NonDisease_N
)

parameters_to_monitor <- c(
  "pooled_sensitivity",
  "pooled_specificity",
  "mu",
  "sd_sensitivity",
  "sd_specificity",
  "rho"
)

# ------------------------------------------------------------
# 6. Run MCMC
# ------------------------------------------------------------
set.seed(20260806)

jags_fit <- jags.model(
  file = textConnection(jags_model),
  data = jags_data,
  n.chains = 4,
  n.adapt = 5000
)

update(
  object = jags_fit,
  n.iter = 20000
)

mcmc_samples <- coda.samples(
  model = jags_fit,
  variable.names = parameters_to_monitor,
  n.iter = 40000,
  thin = 10
)

saveRDS(
  mcmc_samples,
  file.path(output_dir, "Bayesian_Bivariate_MCMC_Samples.rds")
)

posterior <- as.matrix(mcmc_samples)

# ------------------------------------------------------------
# 7. Posterior summaries
# ------------------------------------------------------------
posterior_summary <- function(x) {
  c(
    Posterior_Mean = mean(x),
    Posterior_Median = median(x),
    CrI_2.5 = unname(quantile(x, 0.025)),
    CrI_97.5 = unname(quantile(x, 0.975))
  )
}

result_table <- rbind(
  Sensitivity = posterior_summary(posterior[, "pooled_sensitivity"]),
  Specificity = posterior_summary(posterior[, "pooled_specificity"]),
  SD_Logit_Sensitivity = posterior_summary(posterior[, "sd_sensitivity"]),
  SD_Logit_Specificity = posterior_summary(posterior[, "sd_specificity"]),
  Correlation_Logit_Sensitivity_Specificity = posterior_summary(posterior[, "rho"])
)

result_table <- data.frame(
  Outcome = rownames(result_table),
  result_table,
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  result_table,
  file.path(output_dir, "Bayesian_Posterior_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 8. MCMC convergence diagnostics
# ------------------------------------------------------------
rhat_result <- gelman.diag(
  mcmc_samples,
  multivariate = FALSE
)$psrf

effective_sample_size <- effectiveSize(mcmc_samples)

diagnostic_table <- data.frame(
  Parameter = rownames(rhat_result),
  Rhat = rhat_result[, 1],
  Rhat_Upper_95CI = rhat_result[, 2],
  Effective_Sample_Size = effective_sample_size[rownames(rhat_result)],
  row.names = NULL
)

write.csv(
  diagnostic_table,
  file.path(output_dir, "Bayesian_MCMC_Diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 9. Posterior density plots
# ------------------------------------------------------------
plot_posterior_density <- function(values, outcome_name, colour_value, output_file) {

  summary_values <- posterior_summary(values)

  density_data <- density(values)

  plot_data <- data.frame(
    Value = density_data$x,
    Density = density_data$y
  )

  subtitle_text <- sprintf(
    "Posterior median = %.3f (95%% CrI %.3f to %.3f)",
    summary_values["Posterior_Median"],
    summary_values["CrI_2.5"],
    summary_values["CrI_97.5"]
  )

  p <- ggplot(plot_data, aes(x = Value, y = Density)) +
    geom_area(
      fill = colour_value,
      alpha = 0.35
    ) +
    geom_line(
      colour = colour_value,
      linewidth = 1
    ) +
    geom_vline(
      xintercept = summary_values["Posterior_Median"],
      colour = pal[1],
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    labs(
      title = paste0("Bayesian Posterior Distribution: ", outcome_name),
      subtitle = subtitle_text,
      x = outcome_name,
      y = "Posterior density"
    ) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    theme_classic(base_size = 14, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", colour = pal[1]),
      plot.subtitle = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(16, 35, 16, 16)
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 7,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )
}

plot_posterior_density(
  values = posterior[, "pooled_sensitivity"],
  outcome_name = "Sensitivity",
  colour_value = pal[2],
  output_file = file.path(
    output_dir,
    "Figure_Bayesian_Posterior_Sensitivity.pdf"
  )
)

plot_posterior_density(
  values = posterior[, "pooled_specificity"],
  outcome_name = "Specificity",
  colour_value = pal[3],
  output_file = file.path(
    output_dir,
    "Figure_Bayesian_Posterior_Specificity.pdf"
  )
)

# ------------------------------------------------------------
# 10. Bayesian SROC and posterior AUC
# ------------------------------------------------------------
# SROC is calculated over False Positive Rate (FPR), from 0 to 1.
# The integration direction is explicitly increasing, ensuring AUC > 0.

set.seed(20260806)

n_draws_sroc <- min(3000, nrow(posterior))
selected_rows <- sample(
  seq_len(nrow(posterior)),
  size = n_draws_sroc,
  replace = FALSE
)

fpr_grid <- seq(0.001, 0.999, length.out = 300)

sroc_matrix <- matrix(
  NA,
  nrow = n_draws_sroc,
  ncol = length(fpr_grid)
)

auc_values <- numeric(n_draws_sroc)

for (j in seq_len(n_draws_sroc)) {

  draw_id <- selected_rows[j]

  mu_sens <- posterior[draw_id, "mu[1]"]
  mu_spec <- posterior[draw_id, "mu[2]"]

  sd_sens <- posterior[draw_id, "sd_sensitivity"]
  sd_spec <- posterior[draw_id, "sd_specificity"]

  rho_value <- posterior[draw_id, "rho"]

  # FPR = 1 - Specificity
  # logit(FPR) = -logit(Specificity)
  mu_fpr <- -mu_spec

  # Conditional slope for logit sensitivity given logit FPR
  beta_sroc <- -rho_value * sd_sens / sd_spec

  logit_fpr <- qlogis(fpr_grid)

  predicted_logit_sensitivity <- mu_sens +
    beta_sroc * (logit_fpr - mu_fpr)

  predicted_sensitivity <- plogis(predicted_logit_sensitivity)

  sroc_matrix[j, ] <- predicted_sensitivity

  # Trapezoidal numerical integration over increasing FPR:
  # AUC is therefore always non-negative.
  auc_values[j] <- sum(
    diff(fpr_grid) *
      (
        head(predicted_sensitivity, -1) +
          tail(predicted_sensitivity, -1)
      ) / 2
  )
}

sroc_summary <- data.frame(
  False_Positive_Rate = fpr_grid,
  Median_Sensitivity = apply(sroc_matrix, 2, median),
  Lower_95CrI = apply(sroc_matrix, 2, quantile, probs = 0.025),
  Upper_95CrI = apply(sroc_matrix, 2, quantile, probs = 0.975)
)

auc_summary <- posterior_summary(auc_values)

auc_table <- data.frame(
  Measure = "Posterior_AUC",
  Posterior_Mean = auc_summary["Posterior_Mean"],
  Posterior_Median = auc_summary["Posterior_Median"],
  CrI_2.5 = auc_summary["CrI_2.5"],
  CrI_97.5 = auc_summary["CrI_97.5"]
)

write.csv(
  auc_table,
  file.path(output_dir, "Bayesian_SROC_AUC_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sroc_summary,
  file.path(output_dir, "Bayesian_SROC_Curve_Data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Short subtitle prevents text truncation in the PDF
auc_subtitle <- sprintf(
  "Posterior AUC = %.3f (95%% CrI %.3f to %.3f)",
  auc_summary["Posterior_Median"],
  auc_summary["CrI_2.5"],
  auc_summary["CrI_97.5"]
)

p_sroc <- ggplot() +
  geom_ribbon(
    data = sroc_summary,
    aes(
      x = False_Positive_Rate,
      ymin = Lower_95CrI,
      ymax = Upper_95CrI
    ),
    fill = pal[5],
    alpha = 0.75
  ) +
  geom_line(
    data = sroc_summary,
    aes(
      x = False_Positive_Rate,
      y = Median_Sensitivity
    ),
    colour = pal[1],
    linewidth = 1.2
  ) +
  geom_point(
    data = dat,
    aes(
      x = False_Positive_Rate,
      y = Sensitivity
    ),
    shape = 21,
    size = 3,
    stroke = 0.7,
    fill = pal[3],
    colour = pal[1],
    alpha = 0.85
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.7
  ) +
  labs(
    title = "Bayesian Summary Receiver Operating Characteristic Curve",
    subtitle = auc_subtitle,
    x = "False Positive Rate (1 - Specificity)",
    y = "Sensitivity"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1],
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      colour = "black",
      margin = margin(b = 12)
    ),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(18, 42, 18, 18)
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_Bayesian_SROC_and_AUC.pdf"
  ),
  plot = p_sroc,
  width = 11,
  height = 9,
  units = "in",
  device = "pdf",
  limitsize = FALSE
)

# ------------------------------------------------------------
# 11. MCMC trace plots
# ------------------------------------------------------------
pdf(
  file = file.path(
    output_dir,
    "Figure_MCMC_Trace_Plots.pdf"
  ),
  width = 12,
  height = 9,
  family = "sans"
)

plot(
  mcmc_samples,
  trace = TRUE,
  density = FALSE
)

dev.off()

# ------------------------------------------------------------
# 12. Export concise report
# ------------------------------------------------------------
sink(
  file.path(
    output_dir,
    "Bayesian_Analysis_Report.txt"
  )
)

cat("Bayesian Bivariate Diagnostic Meta-analysis\n")
cat("Outcomes: Sensitivity and Specificity\n")
cat("Number of studies:", nrow(dat), "\n\n")

cat("Posterior summary:\n")
print(result_table, row.names = FALSE)

cat("\nPosterior SROC AUC:\n")
print(auc_table, row.names = FALSE)

cat("\nMCMC diagnostics:\n")
print(diagnostic_table, row.names = FALSE)

sink()

message("Analysis completed successfully.")
message("Output directory: ", output_dir)