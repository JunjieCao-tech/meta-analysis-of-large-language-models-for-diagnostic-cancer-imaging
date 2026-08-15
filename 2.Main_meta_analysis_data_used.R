# =========================================================
# Main Meta-analysis of Diagnostic Performance
# Outcomes: Sensitivity / Specificity / Accuracy / PPV / NPV
# Data source: main.txt only (reported values and reported 95% CI)
# Output: Five separate vector PDF forest plots
# Font: Arial, 14 pt
# =========================================================

# -------------------- 0. Install and load packages --------------------
pkgs <- c("meta", "dplyr")

new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

library(meta)
library(dplyr)

# -------------------- 1. File path --------------------
work_dir <- "C:/Users/30399/Desktop/LLM"
data_file <- file.path(work_dir, "main.txt")

if (!file.exists(data_file)) {
  stop("Cannot find main.txt: ", data_file)
}

# -------------------- 2. Colour palette --------------------
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

col_fixed  <- pal[2]
col_random <- pal[1]
col_study  <- pal[4]

# -------------------- 3. Read main.txt --------------------
dat <- read.delim(
  file = data_file,
  header = TRUE,
  sep = "\t",
  fileEncoding = "GB18030",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "第一作者", "年份",
  "Sensitivity", "Sens CI下限", "Sens CI上限",
  "Specificity", "Spec CI下限", "Spec CI上限",
  "Accuracy", "Acc CI下限", "Acc CI上限",
  "PPV", "PPV CI下限", "PPV CI上限",
  "NPV", "NPV CI下限", "NPV CI上限"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in main.txt:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# -------------------- 4. Helper functions --------------------

# Convert percentage/text values to proportions
to_prop <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "-", ".", "NULL")] <- NA
  x <- gsub("%", "", x)
  x <- gsub(",", "", x)
  
  out <- suppressWarnings(as.numeric(x))
  
  # Values > 1 are regarded as percentages
  out[!is.na(out) & out > 1] <- out[!is.na(out) & out > 1] / 100
  
  return(out)
}

# Correct 0 and 1 values only for logit transformation
bound_prop <- function(p) {
  p <- as.numeric(p)
  p[p <= 0] <- 0.0001
  p[p >= 1] <- 0.9999
  return(p)
}

# Derive logit-scale standard error directly from reported 95% CI
se_from_ci <- function(lower, upper) {
  se <- rep(NA_real_, length(lower))
  
  valid <- !is.na(lower) &
    !is.na(upper) &
    lower > 0 & lower < 1 &
    upper > 0 & upper < 1 &
    upper > lower
  
  se[valid] <- (
    qlogis(upper[valid]) - qlogis(lower[valid])
  ) / (2 * 1.96)
  
  return(se)
}

# Prepare an outcome using only reported values in main.txt
prepare_metric <- function(data, metric_name, value_col, lower_col, upper_col) {
  
  study <- paste0(
    gsub('"', "", as.character(data$第一作者)),
    " (",
    as.character(data$年份),
    ")"
  )
  
  study <- make.unique(study)
  
  effect <- to_prop(data[[value_col]])
  lower_ci <- to_prop(data[[lower_col]])
  upper_ci <- to_prop(data[[upper_col]])
  
  effect_analysis <- bound_prop(effect)
  se_logit <- se_from_ci(lower_ci, upper_ci)
  
  result <- data.frame(
    Study = study,
    Metric = metric_name,
    Effect = effect,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    Effect_analysis = effect_analysis,
    SE_logit = se_logit,
    stringsAsFactors = FALSE
  )
  
  # Only retain studies with reported effect and reported valid 95% CI
  result <- result %>%
    filter(
      !is.na(Effect),
      !is.na(Lower_CI),
      !is.na(Upper_CI),
      !is.na(Effect_analysis),
      !is.na(SE_logit)
    )
  
  return(result)
}

# -------------------- 5. Prepare five outcomes --------------------
dat_sens <- prepare_metric(
  data = dat,
  metric_name = "Sensitivity",
  value_col = "Sensitivity",
  lower_col = "Sens CI下限",
  upper_col = "Sens CI上限"
)

dat_spec <- prepare_metric(
  data = dat,
  metric_name = "Specificity",
  value_col = "Specificity",
  lower_col = "Spec CI下限",
  upper_col = "Spec CI上限"
)

dat_acc <- prepare_metric(
  data = dat,
  metric_name = "Accuracy",
  value_col = "Accuracy",
  lower_col = "Acc CI下限",
  upper_col = "Acc CI上限"
)

dat_ppv <- prepare_metric(
  data = dat,
  metric_name = "PPV",
  value_col = "PPV",
  lower_col = "PPV CI下限",
  upper_col = "PPV CI上限"
)

dat_npv <- prepare_metric(
  data = dat,
  metric_name = "NPV",
  value_col = "NPV",
  lower_col = "NPV CI下限",
  upper_col = "NPV CI上限"
)

analysis_list <- list(
  Sensitivity = dat_sens,
  Specificity = dat_spec,
  Accuracy = dat_acc,
  PPV = dat_ppv,
  NPV = dat_npv
)

# -------------------- 6. Save included data --------------------
used_data <- bind_rows(analysis_list)

write.csv(
  used_data,
  file = file.path(work_dir, "Main_meta_analysis_data_used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nIncluded studies:\n")
for (metric in names(analysis_list)) {
  cat(metric, ": ", nrow(analysis_list[[metric]]), "\n", sep = "")
}

# -------------------- 7. Meta-analysis function --------------------
run_meta <- function(metric_data, metric_name) {
  
  if (nrow(metric_data) < 2) {
    warning(metric_name, ": fewer than 2 eligible studies; skipped.")
    return(NULL)
  }
  
  meta::metagen(
    TE = qlogis(metric_data$Effect_analysis),
    seTE = metric_data$SE_logit,
    studlab = metric_data$Study,
    sm = "PLOGIT",
    common = TRUE,
    random = TRUE,
    method.tau = "REML",
    title = metric_name
  )
}

# -------------------- 8. Register Arial --------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

# -------------------- 9. Generate forest plots --------------------
meta_results <- list()

for (metric in names(analysis_list)) {
  
  metric_data <- analysis_list[[metric]]
  meta_obj <- run_meta(metric_data, metric)
  
  if (is.null(meta_obj)) {
    next
  }
  
  meta_results[[metric]] <- meta_obj
  
  pdf_file <- file.path(
    work_dir,
    paste0("Meta_", metric, "_Fixed_and_Random.pdf")
  )
  
  grDevices::cairo_pdf(
    filename = pdf_file,
    width = 16,
    height = max(10, 5 + nrow(metric_data) * 0.42),
    family = "Arial",
    onefile = TRUE
  )
  
  forest(
    meta_obj,
    common = TRUE,
    random = TRUE,
    backtransf = TRUE,
    
    xlab = paste0(metric, " (Proportion)"),
    smlab = metric,
    
    leftcols = c("studlab", "effect", "ci"),
    leftlabs = c("Study", "Effect", "95% CI"),
    rightcols = FALSE,
    
    digits = 2,
    
    fontsize = 14,
    fs.study = 14,
    fs.axis = 14,
    fs.heading = 14,
    fs.summary = 14,
    fs.xlab = 14,
    
    ff.study = "plain",
    ff.axis = "plain",
    ff.heading = "bold",
    ff.smlab = "bold",
    ff.xlab = "plain",
    ff.common = "bold",
    ff.random = "bold",
    
    col.square = col_study,
    col.square.lines = pal[1],
    
    col.common = col_fixed,
    col.random = col_random,
    
    col.diamond.common = col_fixed,
    col.diamond.random = col_random,
    
    col.diamond.lines.common = pal[1],
    col.diamond.lines.random = pal[1],
    
    text.common = "Fixed-effect model",
    text.random = "Random-effects model",
    
    print.tau2 = TRUE,
    print.I2 = TRUE,
    print.pval.Q = TRUE
  )
  
  dev.off()
  
  cat("Saved: ", pdf_file, "\n", sep = "")
}

# -------------------- 10. Save results text --------------------
summary_file <- file.path(
  work_dir,
  "Main_Meta_Analysis_Results.txt"
)

sink(summary_file)

cat("Main Meta-analysis Results\n")
cat("Fixed-effect and Random-effects Models\n")
cat("Data source: Reported values and reported 95% CI in main.txt only\n")
cat("====================================================\n")

for (metric in names(meta_results)) {
  cat("\n\nOutcome: ", metric, "\n", sep = "")
  cat("----------------------------------------------------\n")
  print(summary(meta_results[[metric]]))
}

sink()

cat("\nCompleted.\n")
cat("Data used: ", file.path(work_dir, "Main_meta_analysis_data_used.csv"), "\n", sep = "")
cat("Results: ", summary_file, "\n", sep = "")