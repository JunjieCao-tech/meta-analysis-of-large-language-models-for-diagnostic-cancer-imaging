# ============================================================
# Model-switch sensitivity analysis for diagnostic meta-analysis
# Models:
# (1) Bivariate random-effects model
#     - Sensitivity and Specificity only
# (2) Univariate random-effects model (GLMM)
#     - Sensitivity, Specificity, Accuracy, PPV, NPV
# (3) Fixed-effect model
#     - Sensitivity, Specificity, Accuracy, PPV, NPV
#
# Output files:
# 1. Bivariate_Random_Effects_Results.csv
# 2. Fixed_Effect_Results.csv
# 3. Model_Switch_All_Model_Results.csv
# 4. Model_Switch_Comparison.csv
# 5. Model_Switch_Data_Used.csv
# 6. Model_Switch_Sensitivity_Plots.pdf
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Install and load packages
# ------------------------------------------------------------
required_packages <- c(
  "meta",
  "metafor",
  "ggplot2"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

library(meta)
library(metafor)
library(ggplot2)

# ------------------------------------------------------------
# 2. Paths and colour palette
# ------------------------------------------------------------
input_file <- "C:/Users/30399/Desktop/LLM/main.txt"

output_dir <- "C:/Users/30399/Desktop/LLM/Model_switch_sensitivity"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

# Required palette
pal <- c(
  "#3C3B8B",
  "#6656A6",
  "#917BBD",
  "#B19CCB",
  "#E5CBE1"
)

# ------------------------------------------------------------
# 3. Read original data
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
# Column 11 = Total sample size
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
  Total_N_Original = suppressWarnings(
    as.numeric(dat_raw[[11]])
  ),
  TP = suppressWarnings(
    as.numeric(dat_raw[[12]])
  ),
  FP = suppressWarnings(
    as.numeric(dat_raw[[13]])
  ),
  FN = suppressWarnings(
    as.numeric(dat_raw[[14]])
  ),
  TN = suppressWarnings(
    as.numeric(dat_raw[[15]])
  )
)

# ------------------------------------------------------------
# 4. Clean data and calculate diagnostic indicators
# ------------------------------------------------------------
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
    dat$NonDisease_N > 0 &
    dat$Total_N > 0,
]

if (nrow(dat) < 2) {
  stop("Fewer than two eligible studies were identified.")
}

rownames(dat) <- NULL

# Study-level diagnostic indicators
dat$Sensitivity <- dat$TP / dat$Disease_N
dat$Specificity <- dat$TN / dat$NonDisease_N
dat$Accuracy <- (dat$TP + dat$TN) / dat$Total_N

dat$PPV_Denominator <- dat$TP + dat$FP
dat$NPV_Denominator <- dat$TN + dat$FN

dat$PPV <- ifelse(
  dat$PPV_Denominator > 0,
  dat$TP / dat$PPV_Denominator,
  NA
)

dat$NPV <- ifelse(
  dat$NPV_Denominator > 0,
  dat$TN / dat$NPV_Denominator,
  NA
)

# Apply continuity correction only for the bivariate model.
# This avoids undefined logits when TP/FN/FP/TN contains zero.
dat$TP_CC <- dat$TP
dat$FP_CC <- dat$FP
dat$FN_CC <- dat$FN
dat$TN_CC <- dat$TN

zero_cell_studies <- (
  dat$TP == 0 |
    dat$FP == 0 |
    dat$FN == 0 |
    dat$TN == 0
)

dat$TP_CC[zero_cell_studies] <- dat$TP_CC[zero_cell_studies] + 0.5
dat$FP_CC[zero_cell_studies] <- dat$FP_CC[zero_cell_studies] + 0.5
dat$FN_CC[zero_cell_studies] <- dat$FN_CC[zero_cell_studies] + 0.5
dat$TN_CC[zero_cell_studies] <- dat$TN_CC[zero_cell_studies] + 0.5

# Export data used
write.csv(
  dat,
  file = file.path(
    output_dir,
    "Model_Switch_Data_Used.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

# Extract pooled result from meta::metaprop object
extract_metaprop_result <- function(
    meta_object,
    outcome_name,
    model_name,
    model_type
) {
  
  if (model_type == "Random") {
    pooled_logit <- meta_object$TE.random
    lower_logit <- meta_object$lower.random
    upper_logit <- meta_object$upper.random
  } else {
    pooled_logit <- meta_object$TE.common
    lower_logit <- meta_object$lower.common
    upper_logit <- meta_object$upper.common
  }
  
  data.frame(
    Outcome = outcome_name,
    Model = model_name,
    Number_of_Studies = meta_object$k,
    
    Pooled_Estimate = plogis(pooled_logit),
    Lower_95CI = plogis(lower_logit),
    Upper_95CI = plogis(upper_logit),
    
    Tau = ifelse(
      model_type == "Random",
      meta_object$tau,
      0
    ),
    
    Tau_Squared = ifelse(
      model_type == "Random",
      meta_object$tau^2,
      0
    ),
    
    I_Squared_Percent = meta_object$I2 * 100,
    
    Q_Statistic = meta_object$Q,
    Q_P_Value = meta_object$pval.Q,
    
    Correlation = NA_real_,
    
    stringsAsFactors = FALSE
  )
}

# Convert result values to a formatted result string
create_result_string <- function(
    estimate,
    lower_ci,
    upper_ci
) {
  sprintf(
    "%.3f (95%% CI %.3f to %.3f)",
    estimate,
    lower_ci,
    upper_ci
  )
}

# ------------------------------------------------------------
# 6. Bivariate random-effects model
# ------------------------------------------------------------
# A bivariate random-effects model is fitted using metafor::rma.mv.
#
# Logit sensitivity:
# logit[TP / (TP + FN)]
#
# Logit false-positive rate:
# logit[FP / (FP + TN)]
#
# Specificity is transformed back as:
# Specificity = 1 - False Positive Rate

dat$Logit_Sensitivity <- log(
  (dat$TP_CC + 0.5) / (dat$FN_CC + 0.5)
)

dat$Var_Logit_Sensitivity <- 1 / (dat$TP_CC + 0.5) +
  1 / (dat$FN_CC + 0.5)

dat$Logit_FPR <- log(
  (dat$FP_CC + 0.5) / (dat$TN_CC + 0.5)
)

dat$Var_Logit_FPR <- 1 / (dat$FP_CC + 0.5) +
  1 / (dat$TN_CC + 0.5)

bivariate_long <- rbind(
  data.frame(
    Study = dat$Study,
    Outcome = "Sensitivity",
    yi = dat$Logit_Sensitivity,
    vi = dat$Var_Logit_Sensitivity
  ),
  data.frame(
    Study = dat$Study,
    Outcome = "False Positive Rate",
    yi = dat$Logit_FPR,
    vi = dat$Var_Logit_FPR
  )
)

bivariate_long$Study <- factor(
  bivariate_long$Study
)

bivariate_long$Outcome <- factor(
  bivariate_long$Outcome,
  levels = c(
    "Sensitivity",
    "False Positive Rate"
  )
)

# Multivariate random-effects model
bivariate_model <- rma.mv(
  yi = yi,
  V = vi,
  mods = ~ Outcome - 1,
  random = ~ Outcome | Study,
  struct = "UN",
  method = "REML",
  data = bivariate_long,
  sparse = TRUE
)

# Fixed effects on logit scale
beta_biv <- as.numeric(bivariate_model$beta)
se_biv <- sqrt(diag(bivariate_model$vb))

# The first coefficient is pooled logit sensitivity.
# The second coefficient is pooled logit false-positive rate.
pooled_logit_sens <- beta_biv[1]
pooled_logit_fpr <- beta_biv[2]

se_logit_sens <- se_biv[1]
se_logit_fpr <- se_biv[2]

# Pooled sensitivity
biv_sens_estimate <- plogis(pooled_logit_sens)
biv_sens_lower <- plogis(
  pooled_logit_sens - 1.96 * se_logit_sens
)
biv_sens_upper <- plogis(
  pooled_logit_sens + 1.96 * se_logit_sens
)

# Pooled specificity
# Specificity = 1 - FPR
biv_spec_estimate <- 1 - plogis(pooled_logit_fpr)

# Note the reversed limits after transformation:
# Lower specificity corresponds to upper FPR.
biv_spec_lower <- 1 - plogis(
  pooled_logit_fpr + 1.96 * se_logit_fpr
)

biv_spec_upper <- 1 - plogis(
  pooled_logit_fpr - 1.96 * se_logit_fpr
)

# Between-study covariance matrix
tau2_sens <- bivariate_model$tau2[1]
tau2_fpr <- bivariate_model$tau2[2]

rho_bivariate <- bivariate_model$rho[1]

# Approximate I2, based on between-study variance and
# mean within-study variance on the logit scale.
mean_vi_sens <- mean(dat$Var_Logit_Sensitivity)
mean_vi_fpr <- mean(dat$Var_Logit_FPR)

i2_sens <- 100 * tau2_sens / (
  tau2_sens + mean_vi_sens
)

i2_fpr <- 100 * tau2_fpr / (
  tau2_fpr + mean_vi_fpr
)

# Q statistic from the multivariate model
biv_q <- bivariate_model$QE
biv_q_p <- bivariate_model$QEp

bivariate_results <- data.frame(
  Outcome = c(
    "Sensitivity",
    "Specificity"
  ),
  
  Model = c(
    "Bivariate random-effects",
    "Bivariate random-effects"
  ),
  
  Number_of_Studies = c(
    nrow(dat),
    nrow(dat)
  ),
  
  Pooled_Estimate = c(
    biv_sens_estimate,
    biv_spec_estimate
  ),
  
  Lower_95CI = c(
    biv_sens_lower,
    biv_spec_lower
  ),
  
  Upper_95CI = c(
    biv_sens_upper,
    biv_spec_upper
  ),
  
  Tau = c(
    sqrt(tau2_sens),
    sqrt(tau2_fpr)
  ),
  
  Tau_Squared = c(
    tau2_sens,
    tau2_fpr
  ),
  
  I_Squared_Percent = c(
    i2_sens,
    i2_fpr
  ),
  
  Q_Statistic = c(
    biv_q,
    biv_q
  ),
  
  Q_P_Value = c(
    biv_q_p,
    biv_q_p
  ),
  
  Correlation = c(
    rho_bivariate,
    rho_bivariate
  ),
  
  stringsAsFactors = FALSE
)

bivariate_results$Pooled_Result <- mapply(
  create_result_string,
  bivariate_results$Pooled_Estimate,
  bivariate_results$Lower_95CI,
  bivariate_results$Upper_95CI
)

# Export bivariate results
write.csv(
  bivariate_results,
  file = file.path(
    output_dir,
    "Bivariate_Random_Effects_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 7. Univariate random-effects models: GLMM
# ------------------------------------------------------------

# Sensitivity
uni_re_sensitivity <- metaprop(
  event = TP,
  n = Disease_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

# Specificity
uni_re_specificity <- metaprop(
  event = TN,
  n = NonDisease_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

# Accuracy
uni_re_accuracy <- metaprop(
  event = TP + TN,
  n = Total_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

# PPV: exclude studies with TP + FP = 0
dat_ppv <- dat[
  !is.na(dat$PPV) &
    dat$PPV_Denominator > 0,
]

uni_re_ppv <- metaprop(
  event = TP,
  n = PPV_Denominator,
  studlab = Study,
  data = dat_ppv,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

# NPV: exclude studies with TN + FN = 0
dat_npv <- dat[
  !is.na(dat$NPV) &
    dat$NPV_Denominator > 0,
]

uni_re_npv <- metaprop(
  event = TN,
  n = NPV_Denominator,
  studlab = Study,
  data = dat_npv,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

univariate_re_results <- rbind(
  extract_metaprop_result(
    uni_re_sensitivity,
    "Sensitivity",
    "Univariate random-effects (GLMM)",
    "Random"
  ),
  extract_metaprop_result(
    uni_re_specificity,
    "Specificity",
    "Univariate random-effects (GLMM)",
    "Random"
  ),
  extract_metaprop_result(
    uni_re_accuracy,
    "Accuracy",
    "Univariate random-effects (GLMM)",
    "Random"
  ),
  extract_metaprop_result(
    uni_re_ppv,
    "PPV",
    "Univariate random-effects (GLMM)",
    "Random"
  ),
  extract_metaprop_result(
    uni_re_npv,
    "NPV",
    "Univariate random-effects (GLMM)",
    "Random"
  )
)

univariate_re_results$Pooled_Result <- mapply(
  create_result_string,
  univariate_re_results$Pooled_Estimate,
  univariate_re_results$Lower_95CI,
  univariate_re_results$Upper_95CI
)

# ------------------------------------------------------------
# 8. Fixed-effect models
# ------------------------------------------------------------

# Sensitivity
fixed_sensitivity <- metaprop(
  event = TP,
  n = Disease_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "Inverse",
  method.tau = "DL",
  common = TRUE,
  random = FALSE,
  incr = 0.5,
  allstudies = TRUE
)

# Specificity
fixed_specificity <- metaprop(
  event = TN,
  n = NonDisease_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "Inverse",
  method.tau = "DL",
  common = TRUE,
  random = FALSE,
  incr = 0.5,
  allstudies = TRUE
)

# Accuracy
fixed_accuracy <- metaprop(
  event = TP + TN,
  n = Total_N,
  studlab = Study,
  data = dat,
  sm = "PLOGIT",
  method = "Inverse",
  method.tau = "DL",
  common = TRUE,
  random = FALSE,
  incr = 0.5,
  allstudies = TRUE
)

# PPV
fixed_ppv <- metaprop(
  event = TP,
  n = PPV_Denominator,
  studlab = Study,
  data = dat_ppv,
  sm = "PLOGIT",
  method = "Inverse",
  method.tau = "DL",
  common = TRUE,
  random = FALSE,
  incr = 0.5,
  allstudies = TRUE
)

# NPV
fixed_npv <- metaprop(
  event = TN,
  n = NPV_Denominator,
  studlab = Study,
  data = dat_npv,
  sm = "PLOGIT",
  method = "Inverse",
  method.tau = "DL",
  common = TRUE,
  random = FALSE,
  incr = 0.5,
  allstudies = TRUE
)

fixed_effect_results <- rbind(
  extract_metaprop_result(
    fixed_sensitivity,
    "Sensitivity",
    "Fixed-effect",
    "Fixed"
  ),
  extract_metaprop_result(
    fixed_specificity,
    "Specificity",
    "Fixed-effect",
    "Fixed"
  ),
  extract_metaprop_result(
    fixed_accuracy,
    "Accuracy",
    "Fixed-effect",
    "Fixed"
  ),
  extract_metaprop_result(
    fixed_ppv,
    "PPV",
    "Fixed-effect",
    "Fixed"
  ),
  extract_metaprop_result(
    fixed_npv,
    "NPV",
    "Fixed-effect",
    "Fixed"
  )
)

fixed_effect_results$Pooled_Result <- mapply(
  create_result_string,
  fixed_effect_results$Pooled_Estimate,
  fixed_effect_results$Lower_95CI,
  fixed_effect_results$Upper_95CI
)

# Export fixed-effect results
write.csv(
  fixed_effect_results,
  file = file.path(
    output_dir,
    "Fixed_Effect_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 9. Combine all model results
# ------------------------------------------------------------
all_model_results <- rbind(
  bivariate_results,
  univariate_re_results,
  fixed_effect_results
)

all_model_results$Pooled_Percent <- round(
  all_model_results$Pooled_Estimate * 100,
  2
)

all_model_results$Lower_95CI_Percent <- round(
  all_model_results$Lower_95CI * 100,
  2
)

all_model_results$Upper_95CI_Percent <- round(
  all_model_results$Upper_95CI * 100,
  2
)

all_model_results$Pooled_Result_Percent <- sprintf(
  "%.1f%% (95%% CI %.1f%% to %.1f%%)",
  all_model_results$Pooled_Percent,
  all_model_results$Lower_95CI_Percent,
  all_model_results$Upper_95CI_Percent
)

all_model_results$Pooled_Estimate <- round(
  all_model_results$Pooled_Estimate,
  4
)

all_model_results$Lower_95CI <- round(
  all_model_results$Lower_95CI,
  4
)

all_model_results$Upper_95CI <- round(
  all_model_results$Upper_95CI,
  4
)

all_model_results$Tau <- round(
  all_model_results$Tau,
  4
)

all_model_results$Tau_Squared <- round(
  all_model_results$Tau_Squared,
  4
)

all_model_results$I_Squared_Percent <- round(
  all_model_results$I_Squared_Percent,
  2
)

all_model_results$Q_Statistic <- round(
  all_model_results$Q_Statistic,
  4
)

all_model_results$Q_P_Value <- format.pval(
  all_model_results$Q_P_Value,
  digits = 4,
  eps = 0.0001
)

all_model_results$Correlation <- round(
  all_model_results$Correlation,
  4
)

write.csv(
  all_model_results,
  file = file.path(
    output_dir,
    "Model_Switch_All_Model_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 10. Create model-comparison table
# ------------------------------------------------------------
comparison_outcomes <- c(
  "Sensitivity",
  "Specificity",
  "Accuracy",
  "PPV",
  "NPV"
)

comparison_table <- data.frame(
  Outcome = comparison_outcomes,
  Bivariate_RE = NA_real_,
  Bivariate_RE_Lower_95CI = NA_real_,
  Bivariate_RE_Upper_95CI = NA_real_,
  
  Univariate_RE = NA_real_,
  Univariate_RE_Lower_95CI = NA_real_,
  Univariate_RE_Upper_95CI = NA_real_,
  
  Fixed_Effect = NA_real_,
  Fixed_Effect_Lower_95CI = NA_real_,
  Fixed_Effect_Upper_95CI = NA_real_,
  
  stringsAsFactors = FALSE
)

for (i in seq_along(comparison_outcomes)) {
  
  current_outcome <- comparison_outcomes[i]
  
  biv_row <- bivariate_results[
    bivariate_results$Outcome == current_outcome,
  ]
  
  uni_row <- univariate_re_results[
    univariate_re_results$Outcome == current_outcome,
  ]
  
  fixed_row <- fixed_effect_results[
    fixed_effect_results$Outcome == current_outcome,
  ]
  
  if (nrow(biv_row) == 1) {
    comparison_table$Bivariate_RE[i] <- biv_row$Pooled_Estimate
    comparison_table$Bivariate_RE_Lower_95CI[i] <- biv_row$Lower_95CI
    comparison_table$Bivariate_RE_Upper_95CI[i] <- biv_row$Upper_95CI
  }
  
  if (nrow(uni_row) == 1) {
    comparison_table$Univariate_RE[i] <- uni_row$Pooled_Estimate
    comparison_table$Univariate_RE_Lower_95CI[i] <- uni_row$Lower_95CI
    comparison_table$Univariate_RE_Upper_95CI[i] <- uni_row$Upper_95CI
  }
  
  if (nrow(fixed_row) == 1) {
    comparison_table$Fixed_Effect[i] <- fixed_row$Pooled_Estimate
    comparison_table$Fixed_Effect_Lower_95CI[i] <- fixed_row$Lower_95CI
    comparison_table$Fixed_Effect_Upper_95CI[i] <- fixed_row$Upper_95CI
  }
}

comparison_table$Difference_Univariate_RE_minus_Bivariate_RE <-
  comparison_table$Univariate_RE -
  comparison_table$Bivariate_RE

comparison_table$Difference_Fixed_Effect_minus_Bivariate_RE <-
  comparison_table$Fixed_Effect -
  comparison_table$Bivariate_RE

comparison_table$Difference_Fixed_Effect_minus_Univariate_RE <-
  comparison_table$Fixed_Effect -
  comparison_table$Univariate_RE

comparison_table[, -1] <- lapply(
  comparison_table[, -1],
  function(x) round(x, 4)
)

write.csv(
  comparison_table,
  file = file.path(
    output_dir,
    "Model_Switch_Comparison.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 11. Prepare plotting data
# ------------------------------------------------------------
plot_model_results <- all_model_results[
  ,
  c(
    "Outcome",
    "Model",
    "Pooled_Estimate",
    "Lower_95CI",
    "Upper_95CI"
  )
]

plot_model_results$Model <- factor(
  plot_model_results$Model,
  levels = c(
    "Bivariate random-effects",
    "Univariate random-effects (GLMM)",
    "Fixed-effect"
  )
)

plot_model_results$Outcome <- factor(
  plot_model_results$Outcome,
  levels = c(
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "PPV",
    "NPV"
  )
)

difference_plot_data <- comparison_table[
  ,
  c(
    "Outcome",
    "Difference_Univariate_RE_minus_Bivariate_RE",
    "Difference_Fixed_Effect_minus_Bivariate_RE",
    "Difference_Fixed_Effect_minus_Univariate_RE"
  )
]

difference_plot_data <- rbind(
  data.frame(
    Outcome = difference_plot_data$Outcome,
    Comparison = "Univariate RE - Bivariate RE",
    Difference = difference_plot_data$Difference_Univariate_RE_minus_Bivariate_RE
  ),
  data.frame(
    Outcome = difference_plot_data$Outcome,
    Comparison = "Fixed Effect - Bivariate RE",
    Difference = difference_plot_data$Difference_Fixed_Effect_minus_Bivariate_RE
  ),
  data.frame(
    Outcome = difference_plot_data$Outcome,
    Comparison = "Fixed Effect - Univariate RE",
    Difference = difference_plot_data$Difference_Fixed_Effect_minus_Univariate_RE
  )
)

difference_plot_data <- difference_plot_data[
  !is.na(difference_plot_data$Difference),
]

difference_plot_data$Outcome <- factor(
  difference_plot_data$Outcome,
  levels = c(
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "PPV",
    "NPV"
  )
)

difference_plot_data$Comparison <- factor(
  difference_plot_data$Comparison,
  levels = c(
    "Univariate RE - Bivariate RE",
    "Fixed Effect - Bivariate RE",
    "Fixed Effect - Univariate RE"
  )
)

# ------------------------------------------------------------
# 12. Create model-switch sensitivity plots with readable labels
# ------------------------------------------------------------

# Model order
plot_model_results$Model <- factor(
  plot_model_results$Model,
  levels = c(
    "Bivariate random-effects",
    "Univariate random-effects (GLMM)",
    "Fixed-effect"
  ),
  labels = c(
    "Bivariate RE",
    "Univariate RE (GLMM)",
    "Fixed effect"
  )
)

plot_model_results$Outcome <- factor(
  plot_model_results$Outcome,
  levels = c(
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "PPV",
    "NPV"
  )
)

# One-line label: estimate (95% CI lower to upper)
plot_model_results$Value_Label <- sprintf(
  "%.3f (%.3f–%.3f)",
  plot_model_results$Pooled_Estimate,
  plot_model_results$Lower_95CI,
  plot_model_results$Upper_95CI
)

# Put labels slightly to the right of confidence intervals
plot_model_results$Label_X <- pmin(
  plot_model_results$Upper_95CI + 0.018,
  1.03
)

# ------------------------------------------------------------
# Figure 1: Faceted horizontal forest plot
# Each diagnostic outcome has its own panel.
# ------------------------------------------------------------
p_model_comparison <- ggplot(
  plot_model_results,
  aes(
    x = Pooled_Estimate,
    y = Model,
    colour = Model,
    shape = Model
  )
) +
  geom_errorbar(
    aes(
      xmin = Lower_95CI,
      xmax = Upper_95CI
    ),
    orientation = "y",
    width = 0.16,
    linewidth = 0.9
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_text(
    aes(
      x = Label_X,
      label = Value_Label
    ),
    colour = "black",
    hjust = 0,
    size = 3.3,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Outcome,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "Bivariate RE" = pal[1],
      "Univariate RE (GLMM)" = pal[2],
      "Fixed effect" = pal[3]
    )
  ) +
  scale_shape_manual(
    values = c(
      "Bivariate RE" = 16,
      "Univariate RE (GLMM)" = 17,
      "Fixed effect" = 15
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1.20),
    breaks = seq(0, 1, 0.1),
    labels = function(x) paste0(x * 100, "%")
  ) +
  labs(
    title = "Pooled Diagnostic Performance Across Models",
    subtitle = "Points indicate pooled estimates; horizontal bars indicate 95% confidence intervals",
    x = "Pooled Estimate",
    y = NULL,
    colour = "Model",
    shape = "Model"
  ) +
  theme_classic(
    base_size = 13,
    base_family = "sans"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1],
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      colour = "black",
      margin = margin(b = 12)
    ),
    strip.background = element_rect(
      fill = pal[5],
      colour = pal[1],
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      colour = pal[1],
      size = 12
    ),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(
      colour = "black",
      size = 11
    ),
    legend.position = "bottom",
    plot.margin = margin(15, 100, 15, 15)
  )

# ------------------------------------------------------------
# Figure 2: Difference between models
# ------------------------------------------------------------
difference_plot_data$Difference_Label <- sprintf(
  "%.3f",
  difference_plot_data$Difference
)

difference_plot_data$Label_Position <- ifelse(
  difference_plot_data$Difference >= 0,
  -0.3,
  1.3
)

p_difference_comparison <- ggplot(
  difference_plot_data,
  aes(
    x = Outcome,
    y = Difference,
    fill = Comparison
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(
      label = Difference_Label,
      vjust = Label_Position
    ),
    position = position_dodge(width = 0.75),
    colour = "black",
    size = 3.3,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Univariate RE - Bivariate RE" = pal[2],
      "Fixed Effect - Bivariate RE" = pal[3],
      "Fixed Effect - Univariate RE" = pal[4]
    )
  ) +
  scale_y_continuous(
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Differences in Pooled Estimates Between Models",
    subtitle = "Positive values indicate a higher estimate for the model listed first",
    x = "Diagnostic Outcome",
    y = "Difference in Pooled Estimate",
    fill = "Comparison"
  ) +
  theme_classic(
    base_size = 13,
    base_family = "sans"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1]
    ),
    plot.subtitle = element_text(colour = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "bottom",
    plot.margin = margin(20, 35, 20, 20)
  )

# ------------------------------------------------------------
# Export both figures into one PDF
# ------------------------------------------------------------
pdf(
  file = file.path(
    output_dir,
    "Model_Switch_Sensitivity_Plots.pdf"
  ),
  width = 14,
  height = 13,
  family = "sans"
)

print(p_model_comparison)
print(p_difference_comparison)

dev.off()

# ------------------------------------------------------------
# 13. Export a concise text report
# ------------------------------------------------------------
sink(
  file.path(
    output_dir,
    "Model_Switch_Sensitivity_Report.txt"
  )
)

cat("Model-switch sensitivity analysis\n")
cat("=================================\n\n")

cat("Number of eligible studies:", nrow(dat), "\n\n")

cat("Models evaluated:\n")
cat("1. Bivariate random-effects model\n")
cat("2. Univariate random-effects GLMM\n")
cat("3. Fixed-effect inverse-variance model\n\n")

cat("Bivariate random-effects results:\n")
print(bivariate_results, row.names = FALSE)

cat("\nUnivariate random-effects results:\n")
print(univariate_re_results, row.names = FALSE)

cat("\nFixed-effect results:\n")
print(fixed_effect_results, row.names = FALSE)

cat("\nModel-switch comparison:\n")
print(comparison_table, row.names = FALSE)

sink()

# ------------------------------------------------------------
# 14. Print key results in console
# ------------------------------------------------------------
print(
  all_model_results[
    ,
    c(
      "Outcome",
      "Model",
      "Number_of_Studies",
      "Pooled_Result",
      "I_Squared_Percent"
    )
  ],
  row.names = FALSE
)

message("Model-switch sensitivity analysis completed successfully.")
message("Output directory: ", output_dir)