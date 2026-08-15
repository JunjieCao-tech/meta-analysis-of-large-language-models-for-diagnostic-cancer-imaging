# ============================================================
# Meta-regression Sensitivity Analysis
# Seven subgroup factors, univariate and bivariate models
#
# Outcomes:
# Sensitivity, Specificity, Accuracy, PPV, NPV
#
# Models:
# (1) Univariate random-effects meta-regression:
#     Sensitivity, Specificity, Accuracy, PPV, NPV
# (2) Bivariate random-effects meta-regression:
#     Sensitivity and Specificity
# (3) Subgroup-specific bivariate SROC AUC
#
# Seven subgroup factors:
# Publication year
# Sample size
# Tumour organ
# Study design
# Input category
# Imaging modality
# Country
#
# All plots are in English and saved as separate PDF files.
# Font: Arial, 14 pt
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

# -------------------- 1. Packages --------------------

required_packages <- c(
  "metafor",
  "ggplot2",
  "Matrix",
  "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(metafor)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

# -------------------- 2. Paths --------------------

work_dir <- "C:/Users/30399/Desktop/LLM"

if (!dir.exists(work_dir)) {
  work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
}

if (!dir.exists(work_dir)) {
  stop("Working directory was not found.")
}

input_file <- file.path(work_dir, "main.txt")

if (!file.exists(input_file)) {
  stop("main.txt was not found.")
}

output_dir <- file.path(work_dir, "6.Meta_regression")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Minimum number of studies required in each subgroup
min_studies_per_group <- 3

# -------------------- 3. Colours and font --------------------

pal <- c(
  "#3C3B8B",
  "#6656A6",
  "#917BBD",
  "#B19CCB",
  "#E5CBE1"
)

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

font_family <- "Arial"

text_size_14 <- 14 / .pt
text_size_12 <- 12 / .pt

open_pdf <- function(filename, width, height) {
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = "Arial",
    pointsize = 14,
    onefile = TRUE
  )
}

# -------------------- 4. Helper functions --------------------

safe_num <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("%", "", x)
  suppressWarnings(as.numeric(x))
}

logit <- function(x) {
  log(x / (1 - x))
}

inv_logit <- function(x) {
  plogis(x)
}

safe_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  ifelse(is.finite(x), x, NA_real_)
}

bind_nonempty <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  
  if (length(x) == 0) {
    return(data.frame())
  }
  
  do.call(rbind, x)
}

# -------------------- 5. Read data --------------------

dat_raw <- tryCatch(
  read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    fileEncoding = "UTF-8-BOM",
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  error = function(e) NULL
)

if (is.null(dat_raw) || nrow(dat_raw) == 0) {
  dat_raw <- tryCatch(
    read.delim(
      input_file,
      header = TRUE,
      sep = "\t",
      quote = "\"",
      fileEncoding = "GB18030",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
}

if (is.null(dat_raw) || nrow(dat_raw) == 0) {
  stop("main.txt could not be read.")
}

if (ncol(dat_raw) < 15) {
  stop("main.txt has fewer than 15 columns.")
}

# -------------------- 6. Prepare diagnostic data --------------------
# Fixed column positions based on main.txt:
# 1  = First author
# 2  = Year
# 3  = Country
# 5  = Tumour organ
# 6  = Imaging modality
# 7  = Input category
# 10 = Study design
# 11 = Total sample size
# 12 = TP
# 13 = FP
# 14 = FN
# 15 = TN

dat <- data.frame(
  Author = as.character(dat_raw[[1]]),
  Year = safe_num(dat_raw[[2]]),
  Country = as.character(dat_raw[[3]]),
  Organ = as.character(dat_raw[[5]]),
  Modality = as.character(dat_raw[[6]]),
  InputType = as.character(dat_raw[[7]]),
  Design = as.character(dat_raw[[10]]),
  N = safe_num(dat_raw[[11]]),
  TP = safe_num(dat_raw[[12]]),
  FP = safe_num(dat_raw[[13]]),
  FN = safe_num(dat_raw[[14]]),
  TN = safe_num(dat_raw[[15]]),
  stringsAsFactors = FALSE
)

dat$Author[is.na(dat$Author) | dat$Author == ""] <- "Unknown"

dat$Study <- make.unique(
  paste0(
    trimws(dat$Author),
    " (",
    dat$Year,
    ")"
  ),
  sep = "_"
)

dat <- dat[
  is.finite(dat$TP) &
    is.finite(dat$FP) &
    is.finite(dat$FN) &
    is.finite(dat$TN) &
    dat$TP >= 0 &
    dat$FP >= 0 &
    dat$FN >= 0 &
    dat$TN >= 0 &
    (dat$TP + dat$FN) > 0 &
    (dat$FP + dat$TN) > 0 &
    (dat$TP + dat$FP) > 0 &
    (dat$TN + dat$FN) > 0,
  ,
  drop = FALSE
]

if (nrow(dat) < 4) {
  stop("Fewer than four eligible studies were found.")
}

dat$Total_from_2x2 <- dat$TP + dat$FP + dat$FN + dat$TN

dat$N[is.na(dat$N) | dat$N <= 0] <- dat$Total_from_2x2[
  is.na(dat$N) | dat$N <= 0
]

# -------------------- 7. Continuity correction --------------------
# Add 0.5 to all cells only when at least one cell equals zero.

dat$Has_zero <- with(
  dat,
  TP == 0 | FP == 0 | FN == 0 | TN == 0
)

dat$TP_cc <- ifelse(dat$Has_zero, dat$TP + 0.5, dat$TP)
dat$FP_cc <- ifelse(dat$Has_zero, dat$FP + 0.5, dat$FP)
dat$FN_cc <- ifelse(dat$Has_zero, dat$FN + 0.5, dat$FN)
dat$TN_cc <- ifelse(dat$Has_zero, dat$TN + 0.5, dat$TN)

# -------------------- 8. Calculate five outcomes --------------------

dat$Sensitivity <- with(
  dat,
  TP_cc / (TP_cc + FN_cc)
)

dat$Specificity <- with(
  dat,
  TN_cc / (TN_cc + FP_cc)
)

dat$Accuracy <- with(
  dat,
  (TP_cc + TN_cc) /
    (TP_cc + FP_cc + FN_cc + TN_cc)
)

dat$PPV <- with(
  dat,
  TP_cc / (TP_cc + FP_cc)
)

dat$NPV <- with(
  dat,
  TN_cc / (TN_cc + FN_cc)
)

dat$FPR <- with(
  dat,
  FP_cc / (FP_cc + TN_cc)
)

# -------------------- 9. Create seven subgroup variables --------------------

# (1) Publication year
dat$Year_group <- ifelse(
  dat$Year >= 2023 & dat$Year <= 2025,
  "2023-2025",
  ifelse(dat$Year == 2026, "2026", NA_character_)
)

dat$Year_group <- factor(
  dat$Year_group,
  levels = c("2023-2025", "2026")
)

# (2) Sample size
median_N <- median(dat$N, na.rm = TRUE)

dat$Sample_size_group <- ifelse(
  dat$N <= median_N,
  paste0("≤ Median (", round(median_N), ")"),
  paste0("> Median (", round(median_N), ")")
)

dat$Sample_size_group <- factor(
  dat$Sample_size_group,
  levels = c(
    paste0("≤ Median (", round(median_N), ")"),
    paste0("> Median (", round(median_N), ")")
  )
)

# (3) Tumour organ
organ_text <- tolower(trimws(dat$Organ))

dat$Tumour_organ_group <- ifelse(
  grepl("thyroid|甲状腺", organ_text),
  "Thyroid",
  ifelse(
    grepl("breast|乳腺|乳房", organ_text),
    "Breast",
    "Other"
  )
)

dat$Tumour_organ_group <- factor(
  dat$Tumour_organ_group,
  levels = c("Thyroid", "Breast", "Other")
)

# (4) Study design
design_text <- tolower(trimws(dat$Design))

dat$Study_design_group <- ifelse(
  grepl("prospective|前瞻", design_text),
  "Prospective",
  ifelse(
    grepl("retrospective|回顾", design_text),
    "Retrospective",
    NA_character_
  )
)

dat$Study_design_group <- factor(
  dat$Study_design_group,
  levels = c("Prospective", "Retrospective")
)

# (5) Input category
input_text <- tolower(trimws(dat$InputType))

dat$Input_type_group <- ifelse(
  grepl(
    "multimodal|multi-input|multi input|多模|多输入|image.*text|text.*image",
    input_text
  ),
  "Multimodal / multi-input",
  ifelse(
    grepl("text|文本", input_text),
    "Text",
    ifelse(
      grepl("image|imaging|影像", input_text),
      "Imaging",
      NA_character_
    )
  )
)

dat$Input_type_group <- factor(
  dat$Input_type_group,
  levels = c(
    "Text",
    "Imaging",
    "Multimodal / multi-input"
  )
)

# (6) Imaging modality
modality_text <- tolower(trimws(dat$Modality))

dat$Imaging_modality_group <- ifelse(
  grepl("ultrasound|超声", modality_text),
  "Ultrasound",
  "Other"
)

dat$Imaging_modality_group <- factor(
  dat$Imaging_modality_group,
  levels = c("Ultrasound", "Other")
)

# (7) Country
country_text <- tolower(trimws(dat$Country))

dat$Country_group <- ifelse(
  grepl("china|中国|mainland china", country_text),
  "China",
  "Other"
)

dat$Country_group <- factor(
  dat$Country_group,
  levels = c("China", "Other")
)

# -------------------- 10. Define subgroup factors --------------------

subgroups <- list(
  Publication_year = list(
    variable = "Year_group",
    title = "Publication year"
  ),
  Sample_size = list(
    variable = "Sample_size_group",
    title = "Sample size"
  ),
  Tumour_organ = list(
    variable = "Tumour_organ_group",
    title = "Tumour organ"
  ),
  Study_design = list(
    variable = "Study_design_group",
    title = "Study design"
  ),
  Input_category = list(
    variable = "Input_type_group",
    title = "Input category"
  ),
  Imaging_modality = list(
    variable = "Imaging_modality_group",
    title = "Imaging modality"
  ),
  Country = list(
    variable = "Country_group",
    title = "Country"
  )
)

write.csv(
  dat,
  file.path(
    output_dir,
    "Meta_Regression_Data_Used.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 11. Get numerator and denominator --------------------

get_events_totals <- function(data, outcome) {
  
  if (outcome == "Sensitivity") {
    return(list(
      events = data$TP_cc,
      totals = data$TP_cc + data$FN_cc
    ))
  }
  
  if (outcome == "Specificity") {
    return(list(
      events = data$TN_cc,
      totals = data$TN_cc + data$FP_cc
    ))
  }
  
  if (outcome == "Accuracy") {
    return(list(
      events = data$TP_cc + data$TN_cc,
      totals = data$TP_cc + data$FP_cc +
        data$FN_cc + data$TN_cc
    ))
  }
  
  if (outcome == "PPV") {
    return(list(
      events = data$TP_cc,
      totals = data$TP_cc + data$FP_cc
    ))
  }
  
  if (outcome == "NPV") {
    return(list(
      events = data$TN_cc,
      totals = data$TN_cc + data$FN_cc
    ))
  }
  
  stop("Unknown outcome.")
}

# -------------------- 12. Prepare data for one moderator --------------------

prepare_moderator_data <- function(data, moderator) {
  
  data_sub <- data[
    !is.na(data[[moderator]]) &
      trimws(as.character(data[[moderator]])) != "",
    ,
    drop = FALSE
  ]
  
  data_sub[[moderator]] <- factor(
    as.character(data_sub[[moderator]])
  )
  
  level_counts <- table(data_sub[[moderator]])
  
  valid_levels <- names(level_counts)[
    level_counts >= min_studies_per_group
  ]
  
  data_sub <- data_sub[
    data_sub[[moderator]] %in% valid_levels,
    ,
    drop = FALSE
  ]
  
  data_sub[[moderator]] <- droplevels(
    factor(data_sub[[moderator]])
  )
  
  data_sub
}

# -------------------- 13. Univariate random-effects model --------------------

fit_univariate_re <- function(data, outcome) {
  
  events_totals <- get_events_totals(data, outcome)
  
  yi <- logit(events_totals$events / events_totals$totals)
  
  vi <- 1 / events_totals$events +
    1 / (events_totals$totals - events_totals$events)
  
  rma(
    yi = yi,
    vi = vi,
    method = "REML"
  )
}

fit_univariate_meta_regression <- function(
    data,
    outcome,
    moderator
) {
  
  events_totals <- get_events_totals(data, outcome)
  
  yi <- logit(events_totals$events / events_totals$totals)
  
  vi <- 1 / events_totals$events +
    1 / (events_totals$totals - events_totals$events)
  
  Moderator <- droplevels(
    factor(data[[moderator]])
  )
  
  rma(
    yi = yi,
    vi = vi,
    mods = ~ Moderator,
    method = "REML"
  )
}

extract_univariate_pool <- function(fit, outcome) {
  
  estimate_logit <- as.numeric(fit$b[1, 1])
  se_logit <- as.numeric(fit$se[1])
  
  data.frame(
    Outcome = outcome,
    Estimate = inv_logit(estimate_logit),
    Lower_95CI = inv_logit(
      estimate_logit - 1.96 * se_logit
    ),
    Upper_95CI = inv_logit(
      estimate_logit + 1.96 * se_logit
    ),
    Tau2 = safe_value(fit$tau2),
    I2_percent = safe_value(fit$I2),
    QE = safe_value(fit$QE),
    QE_p = safe_value(fit$QEp),
    stringsAsFactors = FALSE
  )
}

run_univariate_meta_regression <- function(
    data,
    outcome,
    moderator,
    moderator_title
) {
  
  data_sub <- prepare_moderator_data(
    data,
    moderator
  )
  
  if (
    nrow(data_sub) < 4 ||
    nlevels(data_sub[[moderator]]) < 2
  ) {
    return(NULL)
  }
  
  fit_null <- tryCatch(
    fit_univariate_re(data_sub, outcome),
    error = function(e) NULL
  )
  
  fit_moderator <- tryCatch(
    fit_univariate_meta_regression(
      data_sub,
      outcome,
      moderator
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_null) || is.null(fit_moderator)) {
    return(NULL)
  }
  
  tau2_null <- safe_value(fit_null$tau2)
  tau2_residual <- safe_value(fit_moderator$tau2)
  
  pseudo_r2 <- NA_real_
  
  if (
    is.finite(tau2_null) &&
    tau2_null > 0 &&
    is.finite(tau2_residual)
  ) {
    pseudo_r2 <- 100 *
      (tau2_null - tau2_residual) /
      tau2_null
    
    pseudo_r2 <- max(0, pseudo_r2)
  }
  
  contribution <- data.frame(
    Moderator = moderator_title,
    Outcome = outcome,
    Studies = nrow(data_sub),
    Number_of_groups = nlevels(data_sub[[moderator]]),
    Groups = paste(
      levels(data_sub[[moderator]]),
      collapse = "; "
    ),
    Tau2_null = tau2_null,
    Tau2_residual = tau2_residual,
    Pseudo_R2_percent = pseudo_r2,
    QM = safe_value(fit_moderator$QM),
    QM_p = safe_value(fit_moderator$QMp),
    stringsAsFactors = FALSE
  )
  
  coefficient_result <- data.frame(
    Moderator = moderator_title,
    Outcome = outcome,
    Term = rownames(fit_moderator$b),
    Beta_logit = as.numeric(fit_moderator$b[, 1]),
    SE = as.numeric(fit_moderator$se),
    Z = as.numeric(fit_moderator$zval),
    P = as.numeric(fit_moderator$pval),
    stringsAsFactors = FALSE
  )
  
  subgroup_list <- list()
  
  for (level_now in levels(data_sub[[moderator]])) {
    
    data_group <- data_sub[
      data_sub[[moderator]] == level_now,
      ,
      drop = FALSE
    ]
    
    fit_group <- tryCatch(
      fit_univariate_re(data_group, outcome),
      error = function(e) NULL
    )
    
    if (!is.null(fit_group)) {
      
      group_result <- extract_univariate_pool(
        fit_group,
        outcome
      )
      
      group_result$Moderator <- moderator_title
      group_result$Subgroup <- level_now
      group_result$Studies <- nrow(data_group)
      
      subgroup_list[[level_now]] <- group_result
    }
  }
  
  list(
    contribution = contribution,
    coefficients = coefficient_result,
    subgroup = bind_nonempty(subgroup_list)
  )
}

# -------------------- 14. Bivariate random-effects model --------------------
# Bivariate model:
# logit(Sensitivity) and logit(False-positive rate)

create_bivariate_long_data <- function(
    data,
    moderator = NULL
) {
  
  sensitivity_logit <- logit(
    data$TP_cc / (data$TP_cc + data$FN_cc)
  )
  
  fpr_logit <- logit(
    data$FP_cc / (data$FP_cc + data$TN_cc)
  )
  
  sensitivity_variance <- 1 / data$TP_cc +
    1 / data$FN_cc
  
  fpr_variance <- 1 / data$FP_cc +
    1 / data$TN_cc
  
  long_data <- rbind(
    data.frame(
      Study = data$Study,
      Outcome = "Sensitivity",
      yi = sensitivity_logit,
      vi = sensitivity_variance,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Study = data$Study,
      Outcome = "FPR",
      yi = fpr_logit,
      vi = fpr_variance,
      stringsAsFactors = FALSE
    )
  )
  
  long_data$Study <- factor(long_data$Study)
  
  long_data$Outcome <- factor(
    long_data$Outcome,
    levels = c("Sensitivity", "FPR")
  )
  
  if (!is.null(moderator)) {
    long_data$Moderator <- factor(
      c(
        as.character(data[[moderator]]),
        as.character(data[[moderator]])
      )
    )
  }
  
  long_data
}

fit_bivariate_re <- function(data) {
  
  long_data <- create_bivariate_long_data(data)
  
  rma.mv(
    yi = yi,
    V = vi,
    mods = ~ Outcome - 1,
    random = ~ Outcome | Study,
    struct = "UN",
    method = "REML",
    data = long_data
  )
}

fit_bivariate_meta_regression <- function(
    data,
    moderator
) {
  
  long_data <- create_bivariate_long_data(
    data,
    moderator = moderator
  )
  
  rma.mv(
    yi = yi,
    V = vi,
    mods = ~ Outcome * Moderator - 1,
    random = ~ Outcome | Study,
    struct = "UN",
    method = "REML",
    data = long_data
  )
}

extract_bivariate_parameters <- function(fit) {
  
  coefficient_names <- rownames(fit$b)
  
  sensitivity_id <- grep(
    "^OutcomeSensitivity$",
    coefficient_names
  )[1]
  
  fpr_id <- grep(
    "^OutcomeFPR$",
    coefficient_names
  )[1]
  
  if (is.na(sensitivity_id) || is.na(fpr_id)) {
    return(NULL)
  }
  
  tau_sensitivity2 <- safe_value(fit$tau2[1])
  tau_fpr2 <- safe_value(fit$tau2[2])
  rho <- safe_value(fit$rho)
  
  covariance <- rho *
    sqrt(tau_sensitivity2 * tau_fpr2)
  
  list(
    sensitivity_logit = as.numeric(
      fit$b[sensitivity_id, 1]
    ),
    fpr_logit = as.numeric(
      fit$b[fpr_id, 1]
    ),
    tau_sensitivity2 = tau_sensitivity2,
    tau_fpr2 = tau_fpr2,
    rho = rho,
    covariance = covariance
  )
}

calculate_sroc_auc <- function(fit) {
  
  parameters <- extract_bivariate_parameters(fit)
  
  if (is.null(parameters)) {
    return(NA_real_)
  }
  
  fpr_grid <- seq(
    0.0001,
    0.9999,
    length.out = 5000
  )
  
  if (
    !is.finite(parameters$tau_fpr2) ||
    parameters$tau_fpr2 <= 0 ||
    !is.finite(parameters$covariance)
  ) {
    
    sensitivity_curve <- rep(
      inv_logit(parameters$sensitivity_logit),
      length(fpr_grid)
    )
    
  } else {
    
    sroc_slope <- parameters$covariance /
      parameters$tau_fpr2
    
    predicted_sensitivity_logit <-
      parameters$sensitivity_logit +
      sroc_slope * (
        logit(fpr_grid) -
          parameters$fpr_logit
      )
    
    sensitivity_curve <- inv_logit(
      predicted_sensitivity_logit
    )
  }
  
  auc <- sum(
    diff(fpr_grid) *
      (
        sensitivity_curve[-1] +
          sensitivity_curve[-length(sensitivity_curve)]
      ) / 2
  )
  
  max(0, min(1, auc))
}

extract_bivariate_pool <- function(
    fit,
    subgroup_name
) {
  
  parameters <- extract_bivariate_parameters(fit)
  
  if (is.null(parameters)) {
    return(NULL)
  }
  
  coefficient_names <- rownames(fit$b)
  
  sensitivity_id <- grep(
    "^OutcomeSensitivity$",
    coefficient_names
  )[1]
  
  fpr_id <- grep(
    "^OutcomeFPR$",
    coefficient_names
  )[1]
  
  sensitivity_logit <- as.numeric(
    fit$b[sensitivity_id, 1]
  )
  
  fpr_logit <- as.numeric(
    fit$b[fpr_id, 1]
  )
  
  sensitivity_se <- as.numeric(fit$se[sensitivity_id])
  fpr_se <- as.numeric(fit$se[fpr_id])
  
  pooled_specificity <- 1 - inv_logit(fpr_logit)
  
  specificity_lower <- 1 - inv_logit(
    fpr_logit + 1.96 * fpr_se
  )
  
  specificity_upper <- 1 - inv_logit(
    fpr_logit - 1.96 * fpr_se
  )
  
  joint_heterogeneity <-
    parameters$tau_sensitivity2 *
    parameters$tau_fpr2 *
    (1 - parameters$rho^2)
  
  data.frame(
    Subgroup = subgroup_name,
    
    Sensitivity = inv_logit(sensitivity_logit),
    Sensitivity_lower_95CI = inv_logit(
      sensitivity_logit - 1.96 * sensitivity_se
    ),
    Sensitivity_upper_95CI = inv_logit(
      sensitivity_logit + 1.96 * sensitivity_se
    ),
    
    Specificity = pooled_specificity,
    Specificity_lower_95CI = specificity_lower,
    Specificity_upper_95CI = specificity_upper,
    
    Tau2_sensitivity = parameters$tau_sensitivity2,
    Tau2_FPR = parameters$tau_fpr2,
    Rho_sensitivity_FPR = parameters$rho,
    Joint_heterogeneity = joint_heterogeneity,
    AUC = calculate_sroc_auc(fit),
    
    stringsAsFactors = FALSE
  )
}

# -------------------- 15. Bivariate meta-regression --------------------

run_bivariate_meta_regression <- function(
    data,
    moderator,
    moderator_title
) {
  
  data_sub <- prepare_moderator_data(
    data,
    moderator
  )
  
  if (
    nrow(data_sub) < 4 ||
    nlevels(data_sub[[moderator]]) < 2
  ) {
    return(NULL)
  }
  
  fit_null <- tryCatch(
    fit_bivariate_re(data_sub),
    error = function(e) NULL
  )
  
  fit_moderator <- tryCatch(
    fit_bivariate_meta_regression(
      data_sub,
      moderator
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_null) || is.null(fit_moderator)) {
    return(NULL)
  }
  
  parameters_null <- extract_bivariate_parameters(
    fit_null
  )
  
  parameters_residual <- extract_bivariate_parameters(
    fit_moderator
  )
  
  if (
    is.null(parameters_null) ||
    is.null(parameters_residual)
  ) {
    return(NULL)
  }
  
  joint_heterogeneity_null <-
    parameters_null$tau_sensitivity2 *
    parameters_null$tau_fpr2 *
    (1 - parameters_null$rho^2)
  
  joint_heterogeneity_residual <-
    parameters_residual$tau_sensitivity2 *
    parameters_residual$tau_fpr2 *
    (1 - parameters_residual$rho^2)
  
  joint_pseudo_r2 <- NA_real_
  
  if (
    is.finite(joint_heterogeneity_null) &&
    joint_heterogeneity_null > 0 &&
    is.finite(joint_heterogeneity_residual)
  ) {
    
    joint_pseudo_r2 <- 100 *
      (
        joint_heterogeneity_null -
          joint_heterogeneity_residual
      ) /
      joint_heterogeneity_null
    
    joint_pseudo_r2 <- max(0, joint_pseudo_r2)
  }
  
  moderator_coefficients <- grep(
    "Moderator",
    rownames(fit_moderator$b)
  )
  
  omnibus_test <- tryCatch(
    anova(
      fit_moderator,
      btt = moderator_coefficients
    ),
    error = function(e) NULL
  )
  
  subgroup_list <- list()
  
  for (level_now in levels(data_sub[[moderator]])) {
    
    data_group <- data_sub[
      data_sub[[moderator]] == level_now,
      ,
      drop = FALSE
    ]
    
    fit_group <- tryCatch(
      fit_bivariate_re(data_group),
      error = function(e) NULL
    )
    
    if (!is.null(fit_group)) {
      
      group_result <- extract_bivariate_pool(
        fit_group,
        level_now
      )
      
      if (!is.null(group_result)) {
        group_result$Moderator <- moderator_title
        group_result$Studies <- nrow(data_group)
        
        subgroup_list[[level_now]] <- group_result
      }
    }
  }
  
  subgroup_result <- bind_nonempty(subgroup_list)
  
  max_auc_difference <- NA_real_
  
  if (
    nrow(subgroup_result) >= 2 &&
    sum(is.finite(subgroup_result$AUC)) >= 2
  ) {
    
    max_auc_difference <- max(
      subgroup_result$AUC,
      na.rm = TRUE
    ) - min(
      subgroup_result$AUC,
      na.rm = TRUE
    )
  }
  
  contribution <- data.frame(
    Moderator = moderator_title,
    Studies = nrow(data_sub),
    Number_of_groups = nlevels(data_sub[[moderator]]),
    Groups = paste(
      levels(data_sub[[moderator]]),
      collapse = "; "
    ),
    
    Joint_heterogeneity_null =
      joint_heterogeneity_null,
    
    Joint_heterogeneity_residual =
      joint_heterogeneity_residual,
    
    Joint_Pseudo_R2_percent =
      joint_pseudo_r2,
    
    QM = ifelse(
      is.null(omnibus_test),
      NA_real_,
      safe_value(omnibus_test$QM)
    ),
    
    QM_p = ifelse(
      is.null(omnibus_test),
      NA_real_,
      safe_value(omnibus_test$QMp)
    ),
    
    Overall_AUC = calculate_sroc_auc(fit_null),
    
    Maximum_subgroup_AUC_difference =
      max_auc_difference,
    
    stringsAsFactors = FALSE
  )
  
  coefficient_result <- data.frame(
    Moderator = moderator_title,
    Term = rownames(fit_moderator$b),
    Beta_logit = as.numeric(fit_moderator$b[, 1]),
    SE = as.numeric(fit_moderator$se),
    Z = as.numeric(fit_moderator$zval),
    P = as.numeric(fit_moderator$pval),
    stringsAsFactors = FALSE
  )
  
  list(
    contribution = contribution,
    coefficients = coefficient_result,
    subgroup = subgroup_result
  )
}

# -------------------- 16. Overall univariate random-effects models --------------------

outcomes <- c(
  "Sensitivity",
  "Specificity",
  "Accuracy",
  "PPV",
  "NPV"
)

overall_univariate_list <- list()

for (outcome_now in outcomes) {
  
  fit_now <- tryCatch(
    fit_univariate_re(dat, outcome_now),
    error = function(e) NULL
  )
  
  if (!is.null(fit_now)) {
    overall_univariate_list[[outcome_now]] <-
      extract_univariate_pool(
        fit_now,
        outcome_now
      )
  }
}

overall_univariate_results <- bind_nonempty(
  overall_univariate_list
)

# -------------------- 17. Overall bivariate model --------------------

overall_bivariate_fit <- tryCatch(
  fit_bivariate_re(dat),
  error = function(e) NULL
)

if (is.null(overall_bivariate_fit)) {
  stop("The overall bivariate model could not be fitted.")
}

overall_bivariate_results <- extract_bivariate_pool(
  overall_bivariate_fit,
  "Overall"
)

overall_bivariate_results$Studies <- nrow(dat)

# -------------------- 18. Run univariate meta-regressions --------------------

cat("\nRunning univariate meta-regressions...\n")

univariate_contribution_list <- list()
univariate_coefficient_list <- list()
univariate_subgroup_list <- list()

counter <- 1

for (subgroup_name in names(subgroups)) {
  
  moderator_variable <- subgroups[[subgroup_name]]$variable
  moderator_title <- subgroups[[subgroup_name]]$title
  
  cat("Moderator: ", moderator_title, "\n")
  
  for (outcome_now in outcomes) {
    
    result_now <- run_univariate_meta_regression(
      data = dat,
      outcome = outcome_now,
      moderator = moderator_variable,
      moderator_title = moderator_title
    )
    
    if (!is.null(result_now)) {
      
      univariate_contribution_list[[counter]] <-
        result_now$contribution
      
      univariate_coefficient_list[[counter]] <-
        result_now$coefficients
      
      univariate_subgroup_list[[counter]] <-
        result_now$subgroup
      
      counter <- counter + 1
    }
  }
}

univariate_contribution_results <- bind_nonempty(
  univariate_contribution_list
)

univariate_coefficient_results <- bind_nonempty(
  univariate_coefficient_list
)

univariate_subgroup_results <- bind_nonempty(
  univariate_subgroup_list
)

# -------------------- 19. Run bivariate meta-regressions --------------------

cat("\nRunning bivariate meta-regressions and SROC AUC analyses...\n")

bivariate_contribution_list <- list()
bivariate_coefficient_list <- list()
bivariate_subgroup_list <- list()

counter <- 1

for (subgroup_name in names(subgroups)) {
  
  moderator_variable <- subgroups[[subgroup_name]]$variable
  moderator_title <- subgroups[[subgroup_name]]$title
  
  cat("Moderator: ", moderator_title, "\n")
  
  result_now <- run_bivariate_meta_regression(
    data = dat,
    moderator = moderator_variable,
    moderator_title = moderator_title
  )
  
  if (!is.null(result_now)) {
    
    bivariate_contribution_list[[counter]] <-
      result_now$contribution
    
    bivariate_coefficient_list[[counter]] <-
      result_now$coefficients
    
    bivariate_subgroup_list[[counter]] <-
      result_now$subgroup
    
    counter <- counter + 1
  }
}

bivariate_contribution_results <- bind_nonempty(
  bivariate_contribution_list
)

bivariate_coefficient_results <- bind_nonempty(
  bivariate_coefficient_list
)

bivariate_subgroup_results <- bind_nonempty(
  bivariate_subgroup_list
)

# -------------------- 20. Save CSV files --------------------

write.csv(
  overall_univariate_results,
  file.path(
    output_dir,
    "Overall_Univariate_Random_Effects_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  overall_bivariate_results,
  file.path(
    output_dir,
    "Overall_Bivariate_Random_Effects_and_AUC.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  univariate_contribution_results,
  file.path(
    output_dir,
    "Univariate_Meta_Regression_Contribution.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  univariate_coefficient_results,
  file.path(
    output_dir,
    "Univariate_Meta_Regression_Coefficients.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  univariate_subgroup_results,
  file.path(
    output_dir,
    "Univariate_Subgroup_Pooled_Results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  bivariate_contribution_results,
  file.path(
    output_dir,
    "Bivariate_Meta_Regression_Contribution_and_AUC.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  bivariate_coefficient_results,
  file.path(
    output_dir,
    "Bivariate_Meta_Regression_Coefficients.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  bivariate_subgroup_results,
  file.path(
    output_dir,
    "Bivariate_Subgroup_Sensitivity_Specificity_AUC.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 21. Write text report --------------------

report_file <- file.path(
  output_dir,
  "Meta_Regression_Contribution_Report.txt"
)

sink(report_file)

cat("============================================================\n")
cat("Meta-regression Contribution Analysis\n")
cat("============================================================\n\n")

cat("Number of eligible studies: ", nrow(dat), "\n\n")

cat("Seven subgroup factors:\n")
cat("- Publication year\n")
cat("- Sample size\n")
cat("- Tumour organ\n")
cat("- Study design\n")
cat("- Input category\n")
cat("- Imaging modality\n")
cat("- Country\n\n")

cat("------------------------------------------------------------\n")
cat("Overall univariate random-effects results\n")
cat("------------------------------------------------------------\n")
print(overall_univariate_results, row.names = FALSE)

cat("\n\n")

cat("------------------------------------------------------------\n")
cat("Overall bivariate random-effects results\n")
cat("------------------------------------------------------------\n")
print(overall_bivariate_results, row.names = FALSE)

cat("\n\n")

cat("------------------------------------------------------------\n")
cat("Univariate meta-regression contribution results\n")
cat("------------------------------------------------------------\n")

if (nrow(univariate_contribution_results) > 0) {
  print(univariate_contribution_results, row.names = FALSE)
} else {
  cat("No eligible univariate meta-regression result was available.\n")
}

cat("\n\n")

cat("------------------------------------------------------------\n")
cat("Bivariate meta-regression contribution and AUC results\n")
cat("------------------------------------------------------------\n")

if (nrow(bivariate_contribution_results) > 0) {
  print(bivariate_contribution_results, row.names = FALSE)
} else {
  cat("No eligible bivariate meta-regression result was available.\n")
}

cat("\n\n")

cat("------------------------------------------------------------\n")
cat("Interpretation\n")
cat("------------------------------------------------------------\n")

cat(
  "Pseudo_R2_percent represents the reduction in tau-squared after adding one moderator.\n"
)

cat(
  "A larger Pseudo_R2_percent indicates a greater contribution to between-study heterogeneity.\n"
)

cat(
  "QM_p is the omnibus P value for the moderator.\n"
)

cat(
  "QM_p < 0.05 suggests that the moderator is statistically associated with variation in diagnostic performance.\n"
)

cat(
  "Joint_Pseudo_R2_percent represents the reduction in joint Sensitivity/Specificity heterogeneity.\n"
)

cat(
  "AUC is calculated only from bivariate Sensitivity/Specificity SROC models.\n"
)

cat(
  "Accuracy, PPV, and NPV do not have independent standard SROC AUC values.\n"
)

cat(
  "All meta-regressions are exploratory single-moderator analyses and should not be interpreted as causal effects.\n"
)

sink()

# -------------------- 22. Figure 1: Univariate Pseudo R2 --------------------

if (nrow(univariate_contribution_results) > 0) {
  
  p_univariate_r2 <- ggplot(
    univariate_contribution_results,
    aes(
      x = Moderator,
      y = Pseudo_R2_percent,
      fill = Outcome
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.80),
      width = 0.70,
      colour = "white"
    ) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 0.55
    ) +
    scale_fill_manual(
      values = c(
        "Sensitivity" = pal[1],
        "Specificity" = pal[2],
        "Accuracy" = pal[3],
        "PPV" = pal[4],
        "NPV" = pal[5]
      )
    ) +
    labs(
      title = "Moderator Contribution to Between-study Heterogeneity",
      subtitle = "Univariate random-effects meta-regression",
      x = "Subgroup factor",
      y = "Pseudo R² (%)",
      fill = "Outcome"
    ) +
    theme_classic(
      base_family = font_family,
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        colour = pal[1]
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        colour = pal[2]
      ),
      axis.title = element_text(
        face = "bold",
        colour = pal[1]
      ),
      axis.text.x = element_text(
        angle = 32,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      legend.position = "bottom"
    )
  
  open_pdf(
    file.path(
      output_dir,
      "Figure_1_Univariate_Meta_Regression_Pseudo_R2.pdf"
    ),
    width = 17,
    height = 10
  )
  
  print(p_univariate_r2)
  dev.off()
}

# -------------------- 23. Figure 2: Bivariate Joint Pseudo R2 --------------------

if (nrow(bivariate_contribution_results) > 0) {
  
  p_bivariate_r2 <- ggplot(
    bivariate_contribution_results,
    aes(
      x = Moderator,
      y = Joint_Pseudo_R2_percent
    )
  ) +
    geom_col(
      fill = pal[2],
      width = 0.65
    ) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 0.55
    ) +
    labs(
      title = "Moderator Contribution in the Bivariate Model",
      subtitle = "Joint heterogeneity reduction for Sensitivity and Specificity",
      x = "Subgroup factor",
      y = "Joint Pseudo R² (%)"
    ) +
    theme_classic(
      base_family = font_family,
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        colour = pal[1]
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        colour = pal[2]
      ),
      axis.title = element_text(
        face = "bold",
        colour = pal[1]
      ),
      axis.text.x = element_text(
        angle = 32,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      )
    )
  
  open_pdf(
    file.path(
      output_dir,
      "Figure_2_Bivariate_Meta_Regression_Joint_Pseudo_R2.pdf"
    ),
    width = 17,
    height = 10
  )
  
  print(p_bivariate_r2)
  dev.off()
}

# -------------------- 24. Figure 3: Maximum subgroup AUC difference --------------------

if (nrow(bivariate_contribution_results) > 0) {
  
  p_auc_difference <- ggplot(
    bivariate_contribution_results,
    aes(
      x = Moderator,
      y = Maximum_subgroup_AUC_difference
    )
  ) +
    geom_col(
      fill = pal[3],
      width = 0.65
    ) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 0.55
    ) +
    labs(
      title = "Maximum AUC Difference Between Subgroups",
      subtitle = "AUC from bivariate Sensitivity/Specificity SROC models",
      x = "Subgroup factor",
      y = "Maximum subgroup AUC difference"
    ) +
    theme_classic(
      base_family = font_family,
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        colour = pal[1]
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        colour = pal[2]
      ),
      axis.title = element_text(
        face = "bold",
        colour = pal[1]
      ),
      axis.text.x = element_text(
        angle = 32,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      )
    )
  
  open_pdf(
    file.path(
      output_dir,
      "Figure_3_Maximum_Subgroup_AUC_Difference.pdf"
    ),
    width = 17,
    height = 10
  )
  
  print(p_auc_difference)
  dev.off()
}

# -------------------- 25. Figure 4: Subgroup-specific AUC --------------------

if (nrow(bivariate_subgroup_results) > 0) {
  
  bivariate_subgroup_results$Subgroup_label <- paste0(
    bivariate_subgroup_results$Moderator,
    ": ",
    bivariate_subgroup_results$Subgroup
  )
  
  plot_height <- max(
    10,
    min(
      28,
      0.45 * nrow(bivariate_subgroup_results) + 5
    )
  )
  
  p_subgroup_auc <- ggplot(
    bivariate_subgroup_results,
    aes(
      x = reorder(Subgroup_label, AUC),
      y = AUC,
      fill = Moderator
    )
  ) +
    geom_col(
      width = 0.70,
      colour = "white"
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1)
    ) +
    labs(
      title = "Bivariate SROC AUC by Subgroup",
      subtitle = "Sensitivity and Specificity bivariate random-effects models",
      x = "Subgroup",
      y = "SROC AUC",
      fill = "Subgroup factor"
    ) +
    theme_classic(
      base_family = font_family,
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        colour = pal[1]
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        colour = pal[2]
      ),
      axis.title = element_text(
        face = "bold",
        colour = pal[1]
      ),
      axis.text.x = element_text(
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 12
      ),
      legend.position = "bottom"
    )
  
  open_pdf(
    file.path(
      output_dir,
      "Figure_4_Bivariate_SROC_AUC_by_Subgroup.pdf"
    ),
    width = 18,
    height = plot_height
  )
  
  print(p_subgroup_auc)
  dev.off()
}

# -------------------- 26. Completion message --------------------

cat("\n============================================================\n")
cat("Meta-regression analysis completed.\n")
cat("Output directory:\n", output_dir, "\n")
cat("Number of eligible studies:", nrow(dat), "\n")
cat("============================================================\n")

cat("\nMain CSV files:\n")
cat("- Univariate_Meta_Regression_Contribution.csv\n")
cat("- Univariate_Subgroup_Pooled_Results.csv\n")
cat("- Bivariate_Meta_Regression_Contribution_and_AUC.csv\n")
cat("- Bivariate_Subgroup_Sensitivity_Specificity_AUC.csv\n")

cat("\nSeparate PDF figures:\n")
cat("- Figure_1_Univariate_Meta_Regression_Pseudo_R2.pdf\n")
cat("- Figure_2_Bivariate_Meta_Regression_Joint_Pseudo_R2.pdf\n")
cat("- Figure_3_Maximum_Subgroup_AUC_Difference.pdf\n")
cat("- Figure_4_Bivariate_SROC_AUC_by_Subgroup.pdf\n")