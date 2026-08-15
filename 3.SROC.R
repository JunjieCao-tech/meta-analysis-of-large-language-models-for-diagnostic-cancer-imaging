# ============================================================
# SROC plot: Bivariate random-effects diagnostic meta-analysis
# Arial 14 pt; high-quality vector PDF
# ============================================================

# -------------------- 1. Install and load packages --------------------
required_packages <- c("metafor", "ggplot2", "Matrix", "scales")

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(metafor)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

# -------------------- 2. Working directory --------------------
work_dir <- "C:/Users/30399/Desktop/LLM"

if (!dir.exists(work_dir)) {
  work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
}

if (!dir.exists(work_dir)) {
  stop("Working directory was not found. Please check the folder path.")
}

input_file <- file.path(work_dir, "main.txt")

if (!file.exists(input_file)) {
  stop("main.txt was not found in the working directory.")
}

# -------------------- 3. Color palette --------------------
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

col_sroc       <- pal[1]
col_summary    <- pal[2]
col_study      <- pal[3]
col_confidence <- pal[4]
col_prediction <- pal[5]

# -------------------- 4. Arial font --------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

text_size_14 <- 14 / .pt
text_size_12 <- 12 / .pt

# -------------------- 5. Read data --------------------
dat_raw <- tryCatch(
  read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  ),
  error = function(e) {
    read.delim(
      input_file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "GB18030"
    )
  }
)

# -------------------- 6. Identify columns --------------------
clean_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", x))
}

col_names_clean <- clean_name(colnames(dat_raw))

find_column <- function(candidates, fallback_position = NULL) {
  
  candidates_clean <- clean_name(candidates)
  idx <- which(col_names_clean %in% candidates_clean)
  
  if (length(idx) > 0) {
    return(idx[1])
  }
  
  if (!is.null(fallback_position) &&
      fallback_position <= ncol(dat_raw)) {
    return(fallback_position)
  }
  
  return(NA_integer_)
}

author_col <- find_column(
  c("First author", "Author", "FirstAuthor", "Study", "StudyID"),
  fallback_position = 1
)

year_col <- find_column(
  c("Year", "Publication year", "PublicationYear"),
  fallback_position = 2
)

n_col <- find_column(
  c("N", "Total N", "TotalN", "Sample size", "SampleSize"),
  fallback_position = 11
)

tp_col <- find_column(
  c("TP", "True positive", "TruePositive"),
  fallback_position = 12
)

fp_col <- find_column(
  c("FP", "False positive", "FalsePositive"),
  fallback_position = 13
)

fn_col <- find_column(
  c("FN", "False negative", "FalseNegative"),
  fallback_position = 14
)

tn_col <- find_column(
  c("TN", "True negative", "TrueNegative"),
  fallback_position = 15
)

if (any(is.na(c(tp_col, fp_col, fn_col, tn_col)))) {
  stop("TP, FP, FN, or TN columns could not be identified.")
}

# -------------------- 7. Prepare diagnostic data --------------------
dat <- data.frame(
  Author = as.character(dat_raw[[author_col]]),
  Year   = as.character(dat_raw[[year_col]]),
  N      = suppressWarnings(as.numeric(dat_raw[[n_col]])),
  TP     = suppressWarnings(as.numeric(dat_raw[[tp_col]])),
  FP     = suppressWarnings(as.numeric(dat_raw[[fp_col]])),
  FN     = suppressWarnings(as.numeric(dat_raw[[fn_col]])),
  TN     = suppressWarnings(as.numeric(dat_raw[[tn_col]])),
  stringsAsFactors = FALSE
)

dat$Author[is.na(dat$Author) | dat$Author == ""] <- "Unknown"
dat$Year[is.na(dat$Year) | dat$Year == ""] <- ""

dat$Study <- ifelse(
  dat$Year == "",
  dat$Author,
  paste0(dat$Author, " (", dat$Year, ")")
)

dat <- dat[
  !is.na(dat$TP) &
    !is.na(dat$FP) &
    !is.na(dat$FN) &
    !is.na(dat$TN) &
    dat$TP >= 0 &
    dat$FP >= 0 &
    dat$FN >= 0 &
    dat$TN >= 0,
]

dat$Total_from_2x2 <- dat$TP + dat$FP + dat$FN + dat$TN

dat <- dat[
  dat$Total_from_2x2 > 0 &
    (dat$TP + dat$FN) > 0 &
    (dat$FP + dat$TN) > 0,
]

dat$N[is.na(dat$N) | dat$N <= 0] <- dat$Total_from_2x2[
  is.na(dat$N) | dat$N <= 0
]

if (nrow(dat) < 3) {
  stop("Fewer than three eligible studies with complete TP, FP, FN and TN data.")
}

# -------------------- 8. Continuity correction --------------------
dat$Has_zero <- dat$TP == 0 |
  dat$FP == 0 |
  dat$FN == 0 |
  dat$TN == 0

dat$TP_cc <- ifelse(dat$Has_zero, dat$TP + 0.5, dat$TP)
dat$FP_cc <- ifelse(dat$Has_zero, dat$FP + 0.5, dat$FP)
dat$FN_cc <- ifelse(dat$Has_zero, dat$FN + 0.5, dat$FN)
dat$TN_cc <- ifelse(dat$Has_zero, dat$TN + 0.5, dat$TN)

dat$Sensitivity <- dat$TP_cc / (dat$TP_cc + dat$FN_cc)
dat$Specificity <- dat$TN_cc / (dat$TN_cc + dat$FP_cc)
dat$FPR <- dat$FP_cc / (dat$FP_cc + dat$TN_cc)

dat$logit_sens <- qlogis(dat$Sensitivity)
dat$logit_fpr  <- qlogis(dat$FPR)
dat$logit_spec <- qlogis(dat$Specificity)

dat$var_sens <- 1 / dat$TP_cc + 1 / dat$FN_cc
dat$var_fpr  <- 1 / dat$FP_cc + 1 / dat$TN_cc
dat$var_spec <- 1 / dat$TN_cc + 1 / dat$FP_cc

# -------------------- 9. Fixed and random models --------------------
fit_sens_fixed <- rma.uni(
  yi = dat$logit_sens,
  vi = dat$var_sens,
  method = "FE"
)

fit_spec_fixed <- rma.uni(
  yi = dat$logit_spec,
  vi = dat$var_spec,
  method = "FE"
)

fit_sens_random <- rma.uni(
  yi = dat$logit_sens,
  vi = dat$var_sens,
  method = "REML"
)

fit_spec_random <- rma.uni(
  yi = dat$logit_spec,
  vi = dat$var_spec,
  method = "REML"
)

fixed_sens <- plogis(as.numeric(coef(fit_sens_fixed)))
fixed_spec <- plogis(as.numeric(coef(fit_spec_fixed)))

random_sens_uni <- plogis(as.numeric(coef(fit_sens_random)))
random_spec_uni <- plogis(as.numeric(coef(fit_spec_random)))

I2_sens <- fit_sens_random$I2
I2_spec <- fit_spec_random$I2

# -------------------- 10. Bivariate random-effects model --------------------
long_dat <- do.call(
  rbind,
  lapply(seq_len(nrow(dat)), function(i) {
    data.frame(
      Study = dat$Study[i],
      Outcome = factor(
        c("Sensitivity", "FPR"),
        levels = c("Sensitivity", "FPR")
      ),
      yi = c(dat$logit_sens[i], dat$logit_fpr[i]),
      stringsAsFactors = FALSE
    )
  })
)

V_list <- lapply(seq_len(nrow(dat)), function(i) {
  diag(c(dat$var_sens[i], dat$var_fpr[i]))
})

V_matrix <- as.matrix(Matrix::bdiag(V_list))

fit_bivariate <- rma.mv(
  yi = yi,
  V = V_matrix,
  mods = ~ Outcome - 1,
  random = ~ Outcome | Study,
  struct = "UN",
  method = "REML",
  data = long_dat,
  sparse = FALSE
)

# -------------------- 11. Pooled estimates --------------------
beta_hat <- as.numeric(coef(fit_bivariate))
vcov_beta <- vcov(fit_bivariate)

mu_sens <- beta_hat[1]
mu_fpr  <- beta_hat[2]

se_sens <- sqrt(vcov_beta[1, 1])
se_fpr  <- sqrt(vcov_beta[2, 2])

pooled_sens <- plogis(mu_sens)
pooled_fpr  <- plogis(mu_fpr)
pooled_spec <- 1 - pooled_fpr

pooled_sens_lcl <- plogis(mu_sens - 1.96 * se_sens)
pooled_sens_ucl <- plogis(mu_sens + 1.96 * se_sens)

pooled_fpr_lcl <- plogis(mu_fpr - 1.96 * se_fpr)
pooled_fpr_ucl <- plogis(mu_fpr + 1.96 * se_fpr)

pooled_spec_lcl <- 1 - pooled_fpr_ucl
pooled_spec_ucl <- 1 - pooled_fpr_lcl

log_DOR <- mu_sens - mu_fpr

se_log_DOR <- sqrt(
  vcov_beta[1, 1] +
    vcov_beta[2, 2] -
    2 * vcov_beta[1, 2]
)

DOR <- exp(log_DOR)
DOR_lcl <- exp(log_DOR - 1.96 * se_log_DOR)
DOR_ucl <- exp(log_DOR + 1.96 * se_log_DOR)

# -------------------- 12. SROC curve and AUC --------------------
G_matrix <- fit_bivariate$G

if (is.null(G_matrix) ||
    any(is.na(G_matrix)) ||
    nrow(G_matrix) != 2 ||
    ncol(G_matrix) != 2) {
  G_matrix <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)
}

if (abs(G_matrix[2, 2]) < 1e-10) {
  sroc_slope <- 0
} else {
  sroc_slope <- G_matrix[1, 2] / G_matrix[2, 2]
}

fpr_grid <- seq(0.001, 0.999, length.out = 1000)

sroc_dat <- data.frame(
  FPR = fpr_grid,
  Sensitivity = plogis(
    mu_sens + sroc_slope * (qlogis(fpr_grid) - mu_fpr)
  )
)

AUC <- sum(
  diff(sroc_dat$FPR) *
    (head(sroc_dat$Sensitivity, -1) +
       tail(sroc_dat$Sensitivity, -1)) / 2
)

AUC_quality <- ifelse(
  AUC >= 0.90, "Excellent diagnostic performance",
  ifelse(
    AUC >= 0.80, "Good diagnostic performance",
    ifelse(
      AUC >= 0.70, "Acceptable diagnostic performance",
      ifelse(
        AUC >= 0.60, "Limited diagnostic performance",
        "Poor diagnostic performance"
      )
    )
  )
)

# -------------------- 13. Confidence and prediction regions --------------------
ellipse_points <- function(center, covariance_matrix,
                           level = 0.95, n = 300) {
  
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(theta), sin(theta))
  
  eig <- eigen(covariance_matrix, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  
  transform_matrix <- eig$vectors %*%
    diag(sqrt(eig$values), nrow = 2) %*%
    t(eig$vectors)
  
  radius <- sqrt(qchisq(level, df = 2))
  
  points <- t(
    matrix(center, nrow = 2, ncol = n) +
      radius * transform_matrix %*% circle
  )
  
  data.frame(
    logit_sens = points[, 1],
    logit_fpr = points[, 2]
  )
}

confidence_ellipse <- ellipse_points(
  center = c(mu_sens, mu_fpr),
  covariance_matrix = vcov_beta[1:2, 1:2]
)

confidence_ellipse$Sensitivity <- plogis(confidence_ellipse$logit_sens)
confidence_ellipse$FPR <- plogis(confidence_ellipse$logit_fpr)

prediction_ellipse <- ellipse_points(
  center = c(mu_sens, mu_fpr),
  covariance_matrix = G_matrix
)

prediction_ellipse$Sensitivity <- plogis(prediction_ellipse$logit_sens)
prediction_ellipse$FPR <- plogis(prediction_ellipse$logit_fpr)

# -------------------- 14. Result labels --------------------
result_label <- paste0(
  "Bivariate random-effects model\n\n",
  "Studies: ", nrow(dat), "\n",
  "Pooled sensitivity: ",
  percent(pooled_sens, accuracy = 0.1),
  " (", percent(pooled_sens_lcl, accuracy = 0.1),
  "–", percent(pooled_sens_ucl, accuracy = 0.1), ")\n",
  "Pooled specificity: ",
  percent(pooled_spec, accuracy = 0.1),
  " (", percent(pooled_spec_lcl, accuracy = 0.1),
  "–", percent(pooled_spec_ucl, accuracy = 0.1), ")\n",
  "AUC: ", sprintf("%.3f", AUC), "\n",
  "DOR: ", sprintf("%.2f", DOR),
  " (", sprintf("%.2f", DOR_lcl),
  "–", sprintf("%.2f", DOR_ucl), ")\n",
  "I² sensitivity: ", sprintf("%.1f%%", I2_sens), "\n",
  "I² specificity: ", sprintf("%.1f%%", I2_spec), "\n\n",
  AUC_quality
)

fixed_random_label <- paste0(
  "Fixed-effect pooled estimates:\n",
  "Sensitivity = ", percent(fixed_sens, accuracy = 0.1),
  "; Specificity = ", percent(fixed_spec, accuracy = 0.1),
  "\n\n",
  "Random-effects pooled estimates:\n",
  "Sensitivity = ", percent(random_sens_uni, accuracy = 0.1),
  "; Specificity = ", percent(random_spec_uni, accuracy = 0.1)
)

# -------------------- 15. Create SROC plot --------------------
sroc_plot <- ggplot() +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey65"
  ) +
  
  geom_polygon(
    data = prediction_ellipse,
    aes(x = FPR, y = Sensitivity),
    fill = col_prediction,
    color = col_prediction,
    alpha = 0.32,
    linewidth = 0.7
  ) +
  
  geom_polygon(
    data = confidence_ellipse,
    aes(x = FPR, y = Sensitivity),
    fill = col_confidence,
    color = col_confidence,
    alpha = 0.58,
    linewidth = 0.8
  ) +
  
  geom_line(
    data = sroc_dat,
    aes(x = FPR, y = Sensitivity),
    color = col_sroc,
    linewidth = 1.5,
    lineend = "round"
  ) +
  
  geom_point(
    data = dat,
    aes(x = FPR, y = Sensitivity, size = N),
    shape = 21,
    fill = col_study,
    color = "white",
    stroke = 0.65,
    alpha = 0.92
  ) +
  
  geom_point(
    data = data.frame(FPR = pooled_fpr, Sensitivity = pooled_sens),
    aes(x = FPR, y = Sensitivity),
    shape = 23,
    size = 5.8,
    fill = col_summary,
    color = col_sroc,
    stroke = 1.1
  ) +
  
  annotate(
    "text",
    x = min(pooled_fpr + 0.04, 0.88),
    y = max(pooled_sens - 0.04, 0.08),
    label = "Pooled estimate",
    family = "Arial",
    fontface = "bold",
    size = text_size_14,
    color = col_sroc,
    hjust = 0
  ) +
  
  # Right-upper external results box
  annotate(
    "label",
    x = 1.045,
    y = 0.985,
    label = result_label,
    hjust = 0,
    vjust = 1,
    family = "Arial",
    fontface = "plain",
    size = text_size_14,
    color = "#222222",
    fill = alpha("white", 0.94),
    label.size = 0.35,
    label.r = unit(0.18, "lines"),
    label.padding = unit(0.35, "lines")
  ) +
  
  # Fixed/random results box remains in the lower-right of main plot
  annotate(
    "label",
    x = 0.975,
    y = 0.025,
    label = fixed_random_label,
    hjust = 1,
    vjust = 0,
    family = "Arial",
    fontface = "plain",
    size = text_size_12,
    color = "#222222",
    fill = alpha("white", 0.90),
    label.size = 0.35,
    label.r = unit(0.18, "lines"),
    label.padding = unit(0.32, "lines")
  ) +
  
  scale_size_continuous(
    range = c(2.8, 9),
    name = "Sample size"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1),
    expand = c(0.01, 0.01)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1),
    expand = c(0.01, 0.01)
  ) +
  
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  
  labs(
    title = "Summary Receiver Operating Characteristic (SROC) Curve",
    subtitle = "Bivariate random-effects diagnostic meta-analysis",
    x = "False-positive rate (1 − Specificity)",
    y = "Sensitivity",
    caption = paste0(
      "Solid line: SROC curve; diamond: pooled estimate; ",
      "dark shaded area: 95% confidence region; ",
      "light shaded area: 95% prediction region."
    )
  ) +
  
  theme_classic(base_family = "Arial", base_size = 14) +
  
  theme(
    text = element_text(
      family = "Arial",
      size = 14,
      color = "#222222"
    ),
    
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      family = "Arial",
      size = 14,
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      margin = margin(t = 10, r = 10, b = 10, l = 10)
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 14
    ),
    
    axis.line = element_line(
      color = "#333333",
      linewidth = 0.8
    ),
    
    axis.ticks = element_line(
      color = "#333333",
      linewidth = 0.7
    ),
    
    axis.ticks.length = unit(0.20, "cm"),
    
    # Sample-size legend moved to lower-right OUTSIDE the SROC panel
    legend.position = c(1.255, 0.205),
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 14
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 14
    ),
    
    legend.background = element_rect(
      fill = alpha("white", 0.94),
      color = "grey70",
      linewidth = 0.45
    ),
    
    legend.key = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.caption = element_text(
      family = "Arial",
      size = 12,
      hjust = 0,
      color = "grey25",
      margin = margin(t = 12)
    ),
    
    # Keep enough blank space at right for result box and sample-size legend
    plot.margin = margin(15, 390, 15, 18)
  )

# -------------------- 16. Export high-quality vector PDF --------------------
pdf_file <- file.path(
  work_dir,
  "SROC_Bivariate_Random_Effects_Meta_Analysis.pdf"
)

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 17.5,
  height = 11.5,
  family = "Arial",
  onefile = TRUE
)

print(sroc_plot)
dev.off()

# -------------------- 17. Export study data --------------------
study_data_file <- file.path(
  work_dir,
  "SROC_Analysis_Data_Used.csv"
)

write.csv(
  dat[, c(
    "Study", "N", "TP", "FP", "FN", "TN",
    "Sensitivity", "Specificity", "FPR"
  )],
  study_data_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 18. Export results --------------------
results_file <- file.path(
  work_dir,
  "SROC_Bivariate_Random_Effects_Results.txt"
)

results_text <- c(
  "SROC BIVARIATE RANDOM-EFFECTS DIAGNOSTIC META-ANALYSIS",
  "=======================================================",
  "",
  paste0("Number of included studies: ", nrow(dat)),
  "",
  "Bivariate random-effects pooled estimates:",
  paste0(
    "Sensitivity: ", sprintf("%.4f", pooled_sens),
    " (95% CI ", sprintf("%.4f", pooled_sens_lcl),
    " to ", sprintf("%.4f", pooled_sens_ucl), ")"
  ),
  paste0(
    "Specificity: ", sprintf("%.4f", pooled_spec),
    " (95% CI ", sprintf("%.4f", pooled_spec_lcl),
    " to ", sprintf("%.4f", pooled_spec_ucl), ")"
  ),
  paste0(
    "Diagnostic odds ratio: ", sprintf("%.3f", DOR),
    " (95% CI ", sprintf("%.3f", DOR_lcl),
    " to ", sprintf("%.3f", DOR_ucl), ")"
  ),
  paste0("Area under SROC curve (AUC): ", sprintf("%.4f", AUC)),
  paste0("Interpretation: ", AUC_quality),
  "",
  paste0("I-squared for sensitivity: ", sprintf("%.1f%%", I2_sens)),
  paste0("I-squared for specificity: ", sprintf("%.1f%%", I2_spec))
)

writeLines(results_text, con = results_file, useBytes = TRUE)

cat("\nSROC analysis completed.\n")
cat("PDF saved to:\n", pdf_file, "\n")