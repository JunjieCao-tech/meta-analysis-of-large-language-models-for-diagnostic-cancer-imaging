# ============================================================
# Clinical Application Analysis
# Pooled LR+, LR-, Forest Plots and Fagan Nomograms
# Separate PDF output
# English labels
# Arial 14 where available
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

# ============================================================
# 1. Install and load packages
# ============================================================

required_packages <- c("meta", "showtext", "sysfonts")

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(meta)
  library(showtext)
  library(sysfonts)
})

# ============================================================
# 2. Working directory
# ============================================================

work_dir <- "C:/Users/30399/Desktop/LLM"

if (!dir.exists(work_dir)) {
  work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
}

if (!dir.exists(work_dir)) {
  stop("Working directory was not found.")
}

input_file <- file.path(work_dir, "main.txt")
output_dir <- file.path(work_dir, "10.Clinical_application")

if (!file.exists(input_file)) {
  stop("main.txt was not found.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================
# 3. Font settings
# ============================================================

arial_regular <- "C:/Windows/Fonts/arial.ttf"
arial_bold <- "C:/Windows/Fonts/arialbd.ttf"

if (file.exists(arial_regular)) {
  
  sysfonts::font_add(
    family = "Arial",
    regular = arial_regular,
    bold = if (file.exists(arial_bold)) {
      arial_bold
    } else {
      arial_regular
    }
  )
  
  font_family <- "Arial"
  
} else {
  
  font_family <- "sans"
  warning("Arial was not found. The sans font will be used.")
}

showtext::showtext_auto(enable = TRUE)

# ============================================================
# 4. General settings
# ============================================================

pal <- c(
  "#3C3B8B",
  "#6656A6",
  "#917BBD",
  "#B19CCB",
  "#E5CBE1"
)

open_pdf <- function(filename, width = 12, height = 10) {
  
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = "sans",
    pointsize = 14,
    onefile = TRUE
  )
}

safe_numeric <- function(x) {
  
  x <- trimws(as.character(x))
  x <- gsub("%", "", x)
  x <- gsub(",", "", x)
  
  suppressWarnings(as.numeric(x))
}

format_number <- function(x, digits = 3) {
  
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

format_p <- function(p) {
  
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 0.001) {
    return("<0.001")
  }
  
  formatC(p, format = "f", digits = 3)
}

# ============================================================
# 5. Read data
# ============================================================

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
  
  dat_raw <- read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    fileEncoding = "GB18030",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

if (ncol(dat_raw) < 15) {
  stop("main.txt has fewer than 15 columns.")
}

# ============================================================
# 6. Extract diagnostic data
# ============================================================

dat <- data.frame(
  Author = as.character(dat_raw[[1]]),
  Year = safe_numeric(dat_raw[[2]]),
  N = safe_numeric(dat_raw[[11]]),
  TP = safe_numeric(dat_raw[[12]]),
  FP = safe_numeric(dat_raw[[13]]),
  FN = safe_numeric(dat_raw[[14]]),
  TN = safe_numeric(dat_raw[[15]]),
  stringsAsFactors = FALSE
)

dat$Author <- trimws(dat$Author)

dat$Author[
  is.na(dat$Author) | dat$Author == ""
] <- "Unknown"

dat$Study <- ifelse(
  is.na(dat$Year),
  dat$Author,
  paste0(dat$Author, " (", dat$Year, ")")
)

dat <- dat[
  is.finite(dat$TP) &
    is.finite(dat$FP) &
    is.finite(dat$FN) &
    is.finite(dat$TN) &
    dat$TP >= 0 &
    dat$FP >= 0 &
    dat$FN >= 0 &
    dat$TN >= 0,
]

dat$Total_2x2 <- dat$TP + dat$FP + dat$FN + dat$TN

dat <- dat[
  dat$Total_2x2 > 0 &
    (dat$TP + dat$FN) > 0 &
    (dat$FP + dat$TN) > 0,
]

if (nrow(dat) < 2) {
  stop("Fewer than two eligible studies were identified.")
}

# ============================================================
# 7. Continuity correction
# ============================================================

dat$Any_zero <- dat$TP == 0 |
  dat$FP == 0 |
  dat$FN == 0 |
  dat$TN == 0

dat$TP_cc <- ifelse(dat$Any_zero, dat$TP + 0.5, dat$TP)
dat$FP_cc <- ifelse(dat$Any_zero, dat$FP + 0.5, dat$FP)
dat$FN_cc <- ifelse(dat$Any_zero, dat$FN + 0.5, dat$FN)
dat$TN_cc <- ifelse(dat$Any_zero, dat$TN + 0.5, dat$TN)

dat$Disease_total <- dat$TP_cc + dat$FN_cc
dat$Non_disease_total <- dat$FP_cc + dat$TN_cc

dat$Sensitivity <- dat$TP_cc / dat$Disease_total
dat$Specificity <- dat$TN_cc / dat$Non_disease_total

# ============================================================
# 8. Study-specific likelihood ratios
# ============================================================

dat$LR_positive <- dat$Sensitivity / (1 - dat$Specificity)
dat$LR_negative <- (1 - dat$Sensitivity) / dat$Specificity

dat$log_LR_positive <- log(dat$LR_positive)
dat$log_LR_negative <- log(dat$LR_negative)

dat$var_log_LR_positive <- (
  1 / dat$TP_cc -
    1 / dat$Disease_total +
    1 / dat$FP_cc -
    1 / dat$Non_disease_total
)

dat$var_log_LR_negative <- (
  1 / dat$FN_cc -
    1 / dat$Disease_total +
    1 / dat$TN_cc -
    1 / dat$Non_disease_total
)

dat$SE_log_LR_positive <- sqrt(dat$var_log_LR_positive)
dat$SE_log_LR_negative <- sqrt(dat$var_log_LR_negative)

dat$LR_positive_lower <- exp(
  dat$log_LR_positive -
    1.96 * dat$SE_log_LR_positive
)

dat$LR_positive_upper <- exp(
  dat$log_LR_positive +
    1.96 * dat$SE_log_LR_positive
)

dat$LR_negative_lower <- exp(
  dat$log_LR_negative -
    1.96 * dat$SE_log_LR_negative
)

dat$LR_negative_upper <- exp(
  dat$log_LR_negative +
    1.96 * dat$SE_log_LR_negative
)

write.csv(
  dat,
  file.path(output_dir, "Study_specific_likelihood_ratios.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 9. Meta-analysis of LR+ and LR-
# ============================================================

meta_LR_positive <- meta::metagen(
  TE = dat$log_LR_positive,
  seTE = dat$SE_log_LR_positive,
  studlab = dat$Study,
  data = dat,
  sm = "RR",
  common = FALSE,
  random = TRUE,
  method.tau = "REML",
  method.random.ci = "HK",
  hakn = TRUE
)

meta_LR_negative <- meta::metagen(
  TE = dat$log_LR_negative,
  seTE = dat$SE_log_LR_negative,
  studlab = dat$Study,
  data = dat,
  sm = "RR",
  common = FALSE,
  random = TRUE,
  method.tau = "REML",
  method.random.ci = "HK",
  hakn = TRUE
)

# ============================================================
# 10. Forest plot function
# ============================================================

draw_forest_plot <- function(
    meta_object,
    output_file,
    xlab,
    effect_label
) {
  
  plot_height <- max(
    11,
    6 + meta_object$k * 0.50
  )
  
  open_pdf(
    output_file,
    width = 17,
    height = plot_height
  )
  
  par(
    family = font_family,
    mar = c(5, 4, 4, 2)
  )
  
  meta::forest(
    meta_object,
    
    common = FALSE,
    random = TRUE,
    prediction = TRUE,
    backtransf = TRUE,
    
    xlab = xlab,
    smlab = "",
    
    leftcols = c("studlab"),
    leftlabs = c("Study"),
    
    rightcols = c("effect", "ci", "w.random"),
    rightlabs = c(
      effect_label,
      "95% CI",
      "Weight"
    ),
    
    ref = 1,
    
    col.square = pal[1],
    col.square.lines = pal[1],
    
    col.random = pal[1],
    col.diamond.random = pal[2],
    col.diamond.lines.random = pal[1],
    
    col.predict = pal[4],
    
    fontsize = 14,
    fs.study = 14,
    fs.axis = 14,
    fs.heading = 14,
    fs.random = 14,
    fs.predict = 14,
    fs.hetstat = 12,
    fs.xlab = 14,
    
    ff.study = "plain",
    ff.axis = "plain",
    ff.heading = "bold",
    ff.random = "bold",
    ff.predict = "bold",
    ff.xlab = "plain",
    
    # 显示异质性统计量
    overall.hetstat = TRUE,
    print.tau2 = TRUE,
    print.I2 = TRUE,
    print.pval.Q = TRUE,
    
    # 只保留少量间距，避免异质性文字离图太远
    addrows.below.overall = 2,
    addrows.hetstat = 1,
    spacing = 1,
    
    digits = 2,
    digits.se = 2,
    digits.pval = 3
  )
  
  dev.off()
}

# ============================================================
# 11. Generate forest plots
# ============================================================

draw_forest_plot(
  meta_object = meta_LR_positive,
  output_file = file.path(
    output_dir,
    "Forest_plot_LR_positive.pdf"
  ),
  xlab = "Positive likelihood ratio",
  effect_label = "LR+"
)

draw_forest_plot(
  meta_object = meta_LR_negative,
  output_file = file.path(
    output_dir,
    "Forest_plot_LR_negative.pdf"
  ),
  xlab = "Negative likelihood ratio",
  effect_label = "LR-"
)

# ============================================================
# 12. Fagan nomogram functions
# ============================================================

probability_to_log_odds <- function(p) {
  log(p / (1 - p))
}

calculate_posttest_probability <- function(
    pretest_probability,
    likelihood_ratio
) {
  
  pretest_odds <- pretest_probability /
    (1 - pretest_probability)
  
  posttest_odds <- pretest_odds *
    likelihood_ratio
  
  posttest_odds / (1 + posttest_odds)
}

draw_probability_axis <- function(
    x_position,
    axis_title
) {
  
  probability_ticks <- c(
    0.001, 0.005, 0.01, 0.02, 0.05,
    0.10, 0.20, 0.30, 0.50, 0.70,
    0.80, 0.90, 0.95, 0.98, 0.99,
    0.995, 0.999
  )
  
  y_ticks <- probability_to_log_odds(
    probability_ticks
  )
  
  segments(
    x0 = x_position,
    y0 = min(y_ticks),
    x1 = x_position,
    y1 = max(y_ticks),
    lwd = 1.8,
    col = "#222222"
  )
  
  segments(
    x0 = x_position - 0.018,
    y0 = y_ticks,
    x1 = x_position + 0.018,
    y1 = y_ticks,
    lwd = 0.8,
    col = "#222222"
  )
  
  tick_labels <- paste0(
    formatC(
      probability_ticks * 100,
      format = "fg",
      digits = 3
    ),
    "%"
  )
  
  text(
    x = x_position - 0.030,
    y = y_ticks,
    labels = tick_labels,
    adj = 1,
    cex = 0.78,
    family = font_family
  )
  
  text(
    x = x_position,
    y = max(y_ticks) + 0.78,
    labels = axis_title,
    font = 2,
    cex = 1.05,
    family = font_family
  )
}

draw_fagan_nomogram <- function(
    pretest_probability,
    output_file
) {
  
  posttest_positive <- calculate_posttest_probability(
    pretest_probability,
    pooled_LR_positive
  )
  
  posttest_negative <- calculate_posttest_probability(
    pretest_probability,
    pooled_LR_negative
  )
  
  open_pdf(
    output_file,
    width = 14,
    height = 10
  )
  
  par(
    family = font_family,
    mar = c(5, 4, 5, 4)
  )
  
  probability_ticks <- c(
    0.001, 0.005, 0.01, 0.02, 0.05,
    0.10, 0.20, 0.30, 0.50, 0.70,
    0.80, 0.90, 0.95, 0.98, 0.99,
    0.995, 0.999
  )
  
  y_min <- min(
    probability_to_log_odds(probability_ticks)
  )
  
  y_max <- max(
    probability_to_log_odds(probability_ticks)
  )
  
  plot(
    NA,
    NA,
    type = "n",
    xlim = c(0, 1),
    ylim = c(y_min - 1.4, y_max + 1.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )
  
  x_pre <- 0.18
  x_lr <- 0.50
  x_post <- 0.82
  
  draw_probability_axis(
    x_position = x_pre,
    axis_title = "Pre-test probability"
  )
  
  draw_probability_axis(
    x_position = x_post,
    axis_title = "Post-test probability"
  )
  
  # Likelihood ratio axis
  lr_ticks <- c(
    0.01, 0.02, 0.05, 0.10, 0.20,
    0.50, 1, 2, 5, 10, 20, 50, 100
  )
  
  lr_y <- log(lr_ticks)
  
  segments(
    x0 = x_lr,
    y0 = min(lr_y),
    x1 = x_lr,
    y1 = max(lr_y),
    lwd = 1.8,
    col = "#222222"
  )
  
  segments(
    x0 = x_lr - 0.018,
    y0 = lr_y,
    x1 = x_lr + 0.018,
    y1 = lr_y,
    lwd = 0.8,
    col = "#222222"
  )
  
  text(
    x = x_lr - 0.030,
    y = lr_y,
    labels = formatC(
      lr_ticks,
      format = "fg",
      digits = 3
    ),
    adj = 1,
    cex = 0.78,
    family = font_family
  )
  
  text(
    x = x_lr,
    y = max(lr_y) + 0.78,
    labels = "Likelihood ratio",
    font = 2,
    cex = 1.05,
    family = font_family
  )
  
  title(
    main = paste0(
      "Fagan Nomogram: Pre-test Probability ",
      round(pretest_probability * 100),
      "%"
    ),
    cex.main = 1.30,
    font.main = 2,
    family = font_family
  )
  
  y_pre <- probability_to_log_odds(
    pretest_probability
  )
  
  y_lr_positive <- log(pooled_LR_positive)
  y_lr_negative <- log(pooled_LR_negative)
  
  y_post_positive <- probability_to_log_odds(
    posttest_positive
  )
  
  y_post_negative <- probability_to_log_odds(
    posttest_negative
  )
  
  # Positive test line
  lines(
    x = c(x_pre, x_lr, x_post),
    y = c(y_pre, y_lr_positive, y_post_positive),
    col = pal[1],
    lwd = 3
  )
  
  # Negative test line
  lines(
    x = c(x_pre, x_lr, x_post),
    y = c(y_pre, y_lr_negative, y_post_negative),
    col = pal[3],
    lwd = 3,
    lty = 2
  )
  
  points(
    x = c(x_pre, x_post, x_post),
    y = c(y_pre, y_post_positive, y_post_negative),
    pch = 21,
    bg = c(pal[1], pal[1], pal[3]),
    col = "white",
    cex = 1.6
  )
  
  # Move the pre-test label farther left to avoid overlap with the 50% tick
  text(
    x = x_pre - 0.105,
    y = y_pre,
    labels = paste0(
      "Pre-test: ",
      formatC(
        pretest_probability * 100,
        format = "f",
        digits = 1
      ),
      "%"
    ),
    adj = 1,
    cex = 0.88,
    family = font_family
  )
  
  text(
    x = x_post + 0.035,
    y = y_post_positive,
    labels = paste0(
      "Positive test: ",
      formatC(
        posttest_positive * 100,
        format = "f",
        digits = 1
      ),
      "%"
    ),
    adj = 0,
    cex = 0.88,
    col = pal[1],
    family = font_family
  )
  
  text(
    x = x_post + 0.035,
    y = y_post_negative,
    labels = paste0(
      "Negative test: ",
      formatC(
        posttest_negative * 100,
        format = "f",
        digits = 1
      ),
      "%"
    ),
    adj = 0,
    cex = 0.88,
    col = pal[3],
    family = font_family
  )
  
  # Do not use family in legend()
  legend(
    "bottom",
    inset = c(0, -0.16),
    legend = c(
      paste0(
        "LR+ = ",
        format_number(pooled_LR_positive),
        " (95% CI ",
        format_number(pooled_LR_positive_lower),
        "–",
        format_number(pooled_LR_positive_upper),
        ")"
      ),
      paste0(
        "LR- = ",
        format_number(pooled_LR_negative),
        " (95% CI ",
        format_number(pooled_LR_negative_lower),
        "–",
        format_number(pooled_LR_negative_upper),
        ")"
      )
    ),
    col = c(pal[1], pal[3]),
    lty = c(1, 2),
    lwd = 3,
    bty = "n",
    cex = 0.95,
    xpd = TRUE
  )
  
  dev.off()
  
  data.frame(
    Pretest_probability = pretest_probability,
    Posttest_probability_positive = posttest_positive,
    Posttest_probability_negative = posttest_negative,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 13. Generate Fagan nomograms
# ============================================================

pretest_probabilities <- c(
  0.05,
  0.20,
  0.50
)

fagan_results <- data.frame()

for (p in pretest_probabilities) {
  
  p_label <- paste0(
    round(p * 100),
    "percent"
  )
  
  one_result <- draw_fagan_nomogram(
    pretest_probability = p,
    output_file = file.path(
      output_dir,
      paste0(
        "Fagan_nomogram_pretest_probability_",
        p_label,
        ".pdf"
      )
    )
  )
  
  fagan_results <- rbind(
    fagan_results,
    one_result
  )
}

write.csv(
  fagan_results,
  file.path(
    output_dir,
    "Fagan_nomogram_posttest_probability_results.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 14. Write clinical application report
# ============================================================

report_lines <- c(
  "Clinical Application Analysis",
  "============================================================",
  "",
  paste0("Number of included studies: ", nrow(dat)),
  "",
  "Pooled positive likelihood ratio:",
  paste0(
    "LR+ = ",
    format_number(pooled_LR_positive),
    " (95% CI ",
    format_number(pooled_LR_positive_lower),
    "–",
    format_number(pooled_LR_positive_upper),
    ")"
  ),
  paste0(
    "I2 = ",
    format_number(meta_LR_positive$I2, 1),
    "%"
  ),
  paste0(
    "Heterogeneity p-value = ",
    format_p(meta_LR_positive$pval.Q)
  ),
  "",
  "Pooled negative likelihood ratio:",
  paste0(
    "LR- = ",
    format_number(pooled_LR_negative),
    " (95% CI ",
    format_number(pooled_LR_negative_lower),
    "–",
    format_number(pooled_LR_negative_upper),
    ")"
  ),
  paste0(
    "I2 = ",
    format_number(meta_LR_negative$I2, 1),
    "%"
  ),
  paste0(
    "Heterogeneity p-value = ",
    format_p(meta_LR_negative$pval.Q)
  ),
  "",
  "Fagan nomogram results:"
)

for (i in seq_len(nrow(fagan_results))) {
  
  report_lines <- c(
    report_lines,
    "",
    paste0(
      "Pre-test probability: ",
      formatC(
        fagan_results$Pretest_probability[i] * 100,
        format = "f",
        digits = 0
      ),
      "%"
    ),
    paste0(
      "Post-test probability after a positive test: ",
      formatC(
        fagan_results$Posttest_probability_positive[i] * 100,
        format = "f",
        digits = 1
      ),
      "%"
    ),
    paste0(
      "Post-test probability after a negative test: ",
      formatC(
        fagan_results$Posttest_probability_negative[i] * 100,
        format = "f",
        digits = 1
      ),
      "%"
    )
  )
}

writeLines(
  report_lines,
  file.path(
    output_dir,
    "Clinical_application_analysis_report.txt"
  ),
  useBytes = TRUE
)

# ============================================================
# 15. Console output
# ============================================================

cat("\n============================================================\n")
cat("Clinical application analysis completed.\n")
cat("============================================================\n")

cat(
  "\nPooled LR+ = ",
  format_number(pooled_LR_positive),
  " (95% CI ",
  format_number(pooled_LR_positive_lower),
  "–",
  format_number(pooled_LR_positive_upper),
  ")\n",
  sep = ""
)

cat(
  "Pooled LR- = ",
  format_number(pooled_LR_negative),
  " (95% CI ",
  format_number(pooled_LR_negative_lower),
  "–",
  format_number(pooled_LR_negative_upper),
  ")\n\n",
  sep = ""
)

print(fagan_results)

cat("\nOutput directory:\n")
cat(output_dir, "\n")



# =========================================================
# Meta-analysis of Diagnostic Likelihood Ratios
# Outcomes:
# LR+ = Sensitivity / (1 - Specificity)
# LR- = (1 - Sensitivity) / Specificity
#
# Output:
# Two separate vector PDF forest plots:
# (1) Positive likelihood ratio (LR+)
# (2) Negative likelihood ratio (LR-)
#
# Font: Arial, 14 pt
# =========================================================
# =========================================================
# Diagnostic Likelihood Ratio Meta-analysis
# Outcomes: LR+ and LR-
# Data source: TP / FP / FN / TN in main.txt
# Output: Two separate vector PDF forest plots
# Font: Arial, 14 pt
# =========================================================

# -------------------- 0. Install and load packages --------------------
pkgs <- c("meta", "dplyr")

new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(new_pkgs) > 0) {
  install.packages(
    new_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

library(meta)
library(dplyr)

# -------------------- 1. File path --------------------
work_dir <- "C:/Users/30399/Desktop/LLM"
data_file <- file.path(work_dir, "main.txt")

# If your desktop is under OneDrive, use:
# work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
# data_file <- file.path(work_dir, "main.txt")

if (!file.exists(data_file)) {
  stop("Cannot find main.txt: ", data_file)
}

# -------------------- 2. Colour palette --------------------
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

col_fixed  <- pal[2]
col_random <- pal[1]
col_study  <- pal[4]

# -------------------- 3. Read data --------------------
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
  "TP", "FP", "FN", "TN"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# -------------------- 4. Helper function --------------------
to_num <- function(x) {
  
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "-", ".", "NULL")] <- NA
  x <- gsub(",", "", x)
  
  return(suppressWarnings(as.numeric(x)))
}

# -------------------- 5. Extract 2 × 2 table --------------------
TP <- to_num(dat$TP)
FP <- to_num(dat$FP)
FN <- to_num(dat$FN)
TN <- to_num(dat$TN)

# -------------------- 6. Study labels --------------------
study_label <- paste0(
  gsub('"', "", as.character(dat$第一作者)),
  " (",
  as.character(dat$年份),
  ")"
)

study_label <- make.unique(study_label)

# -------------------- 7. Prepare LR+ / LR- data --------------------
prepare_lr_data <- function(
    TP,
    FP,
    FN,
    TN,
    Study,
    type = c("LR+", "LR-")
) {
  
  type <- match.arg(type)
  
  result <- data.frame(
    Study = Study,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    stringsAsFactors = FALSE
  )
  
  # Exclude missing or invalid 2 × 2 table data
  result <- result %>%
    filter(
      !is.na(TP),
      !is.na(FP),
      !is.na(FN),
      !is.na(TN),
      TP >= 0,
      FP >= 0,
      FN >= 0,
      TN >= 0,
      TP + FN > 0,
      FP + TN > 0
    )
  
  if (nrow(result) == 0) {
    return(result)
  }
  
  # Apply 0.5 continuity correction only if any cell is zero
  has_zero <- result$TP == 0 |
    result$FP == 0 |
    result$FN == 0 |
    result$TN == 0
  
  result$TP_cc <- result$TP
  result$FP_cc <- result$FP
  result$FN_cc <- result$FN
  result$TN_cc <- result$TN
  
  result$TP_cc[has_zero] <- result$TP_cc[has_zero] + 0.5
  result$FP_cc[has_zero] <- result$FP_cc[has_zero] + 0.5
  result$FN_cc[has_zero] <- result$FN_cc[has_zero] + 0.5
  result$TN_cc[has_zero] <- result$TN_cc[has_zero] + 0.5
  
  diseased_total <- result$TP_cc + result$FN_cc
  nondiseased_total <- result$FP_cc + result$TN_cc
  
  sensitivity <- result$TP_cc / diseased_total
  specificity <- result$TN_cc / nondiseased_total
  
  if (type == "LR+") {
    
    # LR+ = Sensitivity / (1 - Specificity)
    result$LR <- sensitivity / (1 - specificity)
    
    # SE of log(LR+)
    result$SE_logLR <- sqrt(
      1 / result$TP_cc -
        1 / diseased_total +
        1 / result$FP_cc -
        1 / nondiseased_total
    )
    
  } else {
    
    # LR- = (1 - Sensitivity) / Specificity
    result$LR <- (1 - sensitivity) / specificity
    
    # SE of log(LR-)
    result$SE_logLR <- sqrt(
      1 / result$FN_cc -
        1 / diseased_total +
        1 / result$TN_cc -
        1 / nondiseased_total
    )
  }
  
  result$logLR <- log(result$LR)
  
  # 95% CI for checking and output
  result$Lower_CI <- exp(
    result$logLR - 1.96 * result$SE_logLR
  )
  
  result$Upper_CI <- exp(
    result$logLR + 1.96 * result$SE_logLR
  )
  
  result <- result %>%
    filter(
      is.finite(LR),
      is.finite(logLR),
      is.finite(SE_logLR),
      is.finite(Lower_CI),
      is.finite(Upper_CI),
      LR > 0,
      SE_logLR > 0,
      Lower_CI > 0,
      Upper_CI > Lower_CI
    )
  
  return(result)
}

dat_lrp <- prepare_lr_data(
  TP = TP,
  FP = FP,
  FN = FN,
  TN = TN,
  Study = study_label,
  type = "LR+"
)

dat_lrn <- prepare_lr_data(
  TP = TP,
  FP = FP,
  FN = FN,
  TN = TN,
  Study = study_label,
  type = "LR-"
)

# -------------------- 8. Save included data --------------------
used_data <- bind_rows(
  dat_lrp %>% mutate(Outcome = "LR+"),
  dat_lrn %>% mutate(Outcome = "LR-")
)

write.csv(
  used_data,
  file = file.path(
    work_dir,
    "Likelihood_Ratio_Data_Used.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nIncluded studies:\n")
cat("LR+: ", nrow(dat_lrp), "\n", sep = "")
cat("LR-: ", nrow(dat_lrn), "\n", sep = "")

# -------------------- 9. Meta-analysis function --------------------
run_meta_lr <- function(lr_data, outcome_name) {
  
  if (nrow(lr_data) < 2) {
    warning(
      outcome_name,
      ": fewer than 2 eligible studies; skipped."
    )
    return(NULL)
  }
  
  meta::metagen(
    TE = lr_data$logLR,
    seTE = lr_data$SE_logLR,
    studlab = lr_data$Study,
    sm = "RR",
    common = TRUE,
    random = TRUE,
    prediction = TRUE,
    method.tau = "REML",
    title = outcome_name
  )
}

# -------------------- 10. Register Arial --------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

# -------------------- 11. Forest plot function --------------------
save_lr_forest <- function(
    meta_obj,
    lr_data,
    outcome_name,
    xlab_text,
    pdf_name
) {
  
  pdf_file <- file.path(work_dir, pdf_name)
  
  grDevices::cairo_pdf(
    filename = pdf_file,
    width = 16,
    height = max(11, 6.5 + nrow(lr_data) * 0.42),
    family = "Arial",
    onefile = TRUE
  )
  
  forest(
    meta_obj,
    
    common = TRUE,
    random = TRUE,
    prediction = TRUE,
    backtransf = TRUE,
    
    xlab = xlab_text,
    smlab = outcome_name,
    
    leftcols = c("studlab", "effect", "ci"),
    leftlabs = c("Study", "Effect", "95% CI"),
    rightcols = FALSE,
    
    # Use log scale.
    # Do NOT define xlim or at manually.
    # The plotting range will expand automatically to include all CIs.
    xlog = TRUE,
    
    digits = 2,
    
    # Font size: Arial, 14 pt
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
    ff.predict = "plain",
    
    # Colours
    col.square = col_study,
    col.square.lines = pal[1],
    
    col.common = col_fixed,
    col.random = col_random,
    col.predict = pal[3],
    
    col.diamond.common = col_fixed,
    col.diamond.random = col_random,
    col.diamond.predict = pal[3],
    
    col.diamond.lines.common = pal[1],
    col.diamond.lines.random = pal[1],
    col.diamond.lines.predict = pal[1],
    
    # Model labels
    text.common = "Fixed-effect model",
    text.random = "Random-effects model",
    text.predict = "Prediction interval",
    
    # Heterogeneity statistics
    print.tau2 = TRUE,
    print.I2 = TRUE,
    print.pval.Q = TRUE
  )
  
  dev.off()
  
  cat("Saved: ", pdf_file, "\n", sep = "")
}

# -------------------- 12. LR+ meta-analysis and forest plot --------------------
meta_lrp <- run_meta_lr(
  lr_data = dat_lrp,
  outcome_name = "Positive likelihood ratio (LR+)"
)

if (!is.null(meta_lrp)) {
  
  save_lr_forest(
    meta_obj = meta_lrp,
    lr_data = dat_lrp,
    outcome_name = "LR+",
    xlab_text = "Positive likelihood ratio (LR+)",
    pdf_name = "Meta_LR_positive_Fixed_and_Random.pdf"
  )
}

# -------------------- 13. LR- meta-analysis and forest plot --------------------
meta_lrn <- run_meta_lr(
  lr_data = dat_lrn,
  outcome_name = "Negative likelihood ratio (LR-)"
)

if (!is.null(meta_lrn)) {
  
  save_lr_forest(
    meta_obj = meta_lrn,
    lr_data = dat_lrn,
    outcome_name = "LR-",
    xlab_text = "Negative likelihood ratio (LR-)",
    pdf_name = "Meta_LR_negative_Fixed_and_Random.pdf"
  )
}

# -------------------- 14. Save results text --------------------
summary_file <- file.path(
  work_dir,
  "Likelihood_Ratio_Meta_Analysis_Results.txt"
)

sink(summary_file)

cat("Diagnostic Likelihood Ratio Meta-analysis Results\n")
cat("Fixed-effect and Random-effects Models\n")
cat("====================================================\n")

if (!is.null(meta_lrp)) {
  cat("\n\nOutcome: Positive likelihood ratio (LR+)\n")
  cat("----------------------------------------------------\n")
  print(summary(meta_lrp))
}

if (!is.null(meta_lrn)) {
  cat("\n\nOutcome: Negative likelihood ratio (LR-)\n")
  cat("----------------------------------------------------\n")
  print(summary(meta_lrn))
}

sink()

cat("\nCompleted.\n")
cat(
  "Data used: ",
  file.path(work_dir, "Likelihood_Ratio_Data_Used.csv"),
  "\n",
  sep = ""
)
cat("Results: ", summary_file, "\n", sep = "")
