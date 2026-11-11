# Prefer the analysis-local, version-controlled R library when available.
local_analysis_lib <- file.path(Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis"), "R_packages")
if (dir.exists(local_analysis_lib)) .libPaths(c(local_analysis_lib, .libPaths()))

# ============================================================
# Sensitivity analysis for multiple observations from one patient
# Source: main.txt and Table 1.xlsx
# Output: E:/LLM/投稿版本/npj/R1/Statistical analysis/12.sensitivity analysis
# ============================================================

if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C") {
  invisible(suppressWarnings(
    Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
  ))
}

work_dir <- Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis")
Encoding(work_dir) <- "UTF-8"
main_file <- file.path(work_dir, "main.txt")
table1_file <- Sys.getenv("ANALYSIS_TABLE1_FILE", unset = file.path(dirname(work_dir), "TABLES",
                         "Table 1.xlsx"))
output_dir <- Sys.getenv(
  "ANALYSIS_OUTPUT_DIR",
  unset = file.path(work_dir, "12.sensitivity analysis")
)

if (!file.exists(main_file) || !file.exists(table1_file)) {
  stop("main.txt or the specified Table 1 workbook was not found")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("readxl", "metafor", "Matrix", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install required R package: ", pkg)
  }
}
source(file.path(work_dir, "0.Verify_Main_Data.R"), encoding = "UTF-8")
assert_main_data(main_file)

# Exact palette used in the other analysis scripts.
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")
font_family <- "Arial"

main <- read.delim(
  main_file, header = TRUE, sep = "\t", fileEncoding = "GB18030",
  check.names = FALSE, stringsAsFactors = FALSE
)
table1 <- as.data.frame(readxl::read_excel(
  table1_file, sheet = "Sheet1", skip = 1, .name_repair = "minimal"
))
required <- c("First author", "Year", "N", "Sample unit", "Unique patients")
if (!all(required %in% names(table1))) {
  stop("Table 1 is missing: ", paste(setdiff(required, names(table1)), collapse = ", "))
}
table1 <- table1[!is.na(suppressWarnings(as.integer(table1[["No"]]))), ]
if (nrow(table1) != 21L) stop("Expected 21 Table 1 rows including the red overlap row")

normalize_author <- function(x) tolower(gsub("[^a-z0-9]", "", as.character(x)))
main_key <- paste(normalize_author(main[[1]]), main[[2]], main[[11]], sep = "|")
table_key <- paste(normalize_author(table1[["First author"]]),
                   table1[["Year"]], table1[["N"]], sep = "|")
idx <- match(main_key, table_key)
if (anyNA(idx) || anyDuplicated(table_key)) {
  stop("Could not uniquely match every main.txt study to Table 1 by author, year and N")
}
t1 <- table1[idx, ]

dat <- data.frame(
  Study = paste0(as.character(main[[1]]), " (", main[[2]], ")"),
  Author = as.character(main[[1]]),
  Year = as.integer(main[[2]]),
  TP = as.numeric(main[[12]]), FP = as.numeric(main[[13]]),
  FN = as.numeric(main[[14]]), TN = as.numeric(main[[15]]),
  Sample_unit = as.character(t1[["Sample unit"]]),
  N_Table1 = as.integer(t1[["N"]]),
  Patients_reported = as.character(t1[["Unique patients"]]),
  stringsAsFactors = FALSE
)
dat$N <- rowSums(dat[, c("TP", "FP", "FN", "TN")])
if (anyNA(dat$N_Table1) || any(dat$N != dat$N_Table1)) {
  stop("A Table 1 N does not match the selected main.txt four-cell table")
}

parse_patient_count <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || grepl("^(NR|Not report)", x, ignore.case = TRUE)) {
    return(c(count = NA_real_, upper_bound = 0))
  }
  upper <- startsWith(x, "\u2264") || grepl("^(<=|at most)", x, ignore.case = TRUE)
  n <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  if (!is.finite(n) || n <= 0) stop("Unrecognized patient count: ", x)
  c(count = n, upper_bound = as.integer(upper))
}
stopifnot(
  parse_patient_count("\u22645")[["upper_bound"]] == 1,
  parse_patient_count("<=5")[["upper_bound"]] == 1,
  parse_patient_count("at most 5")[["upper_bound"]] == 1,
  parse_patient_count("5")[["upper_bound"]] == 0
)
parsed <- t(vapply(dat$Patients_reported, parse_patient_count,
                   numeric(2)))
dat$Unique_patients <- parsed[, "count"]
dat$Patient_count_is_upper_bound <- as.logical(parsed[, "upper_bound"])
if (any(!is.na(dat$Unique_patients) & dat$Unique_patients > dat$N)) {
  stop("A Table 1 unique-patient count exceeds its observation count")
}
dat$Multiplicity_lower_bound <- dat$N / dat$Unique_patients
dat$Known_repeated_observations <- !is.na(dat$Unique_patients) &
  dat$Unique_patients < dat$N
dat$Confirmed_one_per_patient <- !is.na(dat$Unique_patients) &
  dat$Unique_patients == dat$N & !dat$Patient_count_is_upper_bound
dat$Patient_count_unknown <- is.na(dat$Unique_patients)
if (sum(dat$Known_repeated_observations) < 1 ||
    sum(dat$Confirmed_one_per_patient) < 4) {
  stop("Patient-count classifications are inconsistent with Table 1")
}

# Equal-size-cluster approximation: DE = 1 + (mean observations per patient - 1) * ICC.
# For Sievert, Table 1 states at most 5 source patients; using 5 is the
# smallest defensible multiplicity and therefore a lower-bound DE.
design_effect <- function(icc) {
  de <- rep(1, nrow(dat))
  repeat_rows <- which(dat$Known_repeated_observations)
  de[repeat_rows] <- 1 +
    (dat$Multiplicity_lower_bound[repeat_rows] - 1) * icc
  de
}
for (icc in c(0.05, 0.10, 0.20, 1.00)) {
  tag <- sprintf("%03d", round(icc * 100))
  dat[[paste0("DE_ICC_", tag)]] <- design_effect(icc)
  dat[[paste0("N_eff_ICC_", tag)]] <- dat$N / design_effect(icc)
}

# The model matches script 3's logit sensitivity / logit FPR bivariate REML.
# Its zero-cell correction is applied to all four cells of an affected study.
# DE inflates BOTH within-study logit variances because class-specific patient
# multiplicities and patient-level 2x2 data were not reported.
fit_bivariate <- function(x, de, scenario) {
  if (nrow(x) < 4L) stop("Too few studies in ", scenario)
  a <- as.matrix(x[, c("TP", "FP", "FN", "TN")])
  zero <- rowSums(a == 0) > 0
  a[zero, ] <- a[zero, ] + 0.5
  sens <- log(a[, "TP"] / a[, "FN"])
  fpr <- log(a[, "FP"] / a[, "TN"])
  vs <- (1 / a[, "TP"] + 1 / a[, "FN"]) * de
  vf <- (1 / a[, "FP"] + 1 / a[, "TN"]) * de
  long <- data.frame(
    Study = factor(rep(x$Study, each = 2)),
    Outcome = factor(rep(c("Sensitivity", "FPR"), nrow(x)),
                     levels = c("Sensitivity", "FPR")),
    yi = as.vector(t(cbind(sens, fpr)))
  )
  V <- Matrix::Diagonal(x = as.vector(t(cbind(vs, vf))))
  fit <- metafor::rma.mv(
    yi = yi, V = V, mods = ~ Outcome - 1,
    random = ~ Outcome | Study, struct = "UN", method = "REML",
    data = long, sparse = FALSE
  )
  b <- as.numeric(stats::coef(fit))
  se <- sqrt(diag(stats::vcov(fit)))
  z <- stats::qnorm(0.975)
  data.frame(
    Scenario = scenario,
    Metric = c("Sensitivity", "Specificity"),
    Estimate = c(stats::plogis(b[1]), 1 - stats::plogis(b[2])),
    Lower = c(stats::plogis(b[1] - z * se[1]),
              1 - stats::plogis(b[2] + z * se[2])),
    Upper = c(stats::plogis(b[1] + z * se[1]),
              1 - stats::plogis(b[2] - z * se[2])),
    Studies = nrow(x), stringsAsFactors = FALSE
  )
}

scenario_names <- c(
  "Primary: no clustering adjustment",
  "Design effect: ICC 0.05",
  "Design effect: ICC 0.10",
  "Design effect: ICC 0.20",
  "Exclude known repeated observations",
  "Design effect: ICC 1.00 (extreme)",
  "Confirmed one observation per patient"
)
fits <- list(
  fit_bivariate(dat, rep(1, nrow(dat)), scenario_names[1]),
  fit_bivariate(dat, design_effect(0.05), scenario_names[2]),
  fit_bivariate(dat, design_effect(0.10), scenario_names[3]),
  fit_bivariate(dat, design_effect(0.20), scenario_names[4]),
  fit_bivariate(dat[!dat$Known_repeated_observations, ],
                rep(1, sum(!dat$Known_repeated_observations)), scenario_names[5]),
  fit_bivariate(dat, design_effect(1.00), scenario_names[6]),
  fit_bivariate(dat[dat$Confirmed_one_per_patient, ],
                rep(1, sum(dat$Confirmed_one_per_patient)), scenario_names[7])
)
results <- do.call(rbind, fits)
results$CI_width_pp <- 100 * (results$Upper - results$Lower)
baseline <- results[results$Scenario == scenario_names[1], ]
results$Estimate_change_pp <- 100 * (results$Estimate -
  baseline$Estimate[match(results$Metric, baseline$Metric)])
results$CI_width_change_pp <- results$CI_width_pp -
  baseline$CI_width_pp[match(results$Metric, baseline$Metric)]

# Cross-check against the updated primary SROC results (four-decimal report).
sroc_file <- file.path(work_dir, "3.SROC",
                       "SROC_Bivariate_Random_Effects_Results.txt")
if (file.exists(sroc_file)) {
  sroc_text <- readLines(sroc_file, warn = FALSE)
  for (metric in c("Sensitivity", "Specificity")) {
    line <- grep(paste0("^", metric, ":"), sroc_text, value = TRUE)
    if (length(line) != 1L) stop("Could not cross-check SROC ", metric)
    published <- as.numeric(sub(".*?: ([0-9.]+).*", "\\1", line))
    current <- baseline$Estimate[baseline$Metric == metric]
    if (abs(current - published) > 0.00005) {
      stop("Baseline clustering analysis differs from Figure 3A: ", metric)
    }
  }
}

write.csv(dat, file.path(output_dir, "Patient_Clustering_Study_Map.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(results,
          file.path(output_dir, "Patient_Clustering_Sensitivity_Results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Plot five main scenarios. The ICC=1 bound and strict one-per-patient subset
# remain in the results table; they have distinct, more restrictive meanings.
plot_names <- scenario_names[1:5]
plot_dat <- results[results$Scenario %in% plot_names, ]
plot_dat$Scenario <- factor(plot_dat$Scenario, levels = rev(plot_names))
plot_dat$Metric <- factor(plot_dat$Metric,
                          levels = c("Sensitivity", "Specificity"))
plot_dat$Label <- sprintf("%.1f%% (%.1f%%\u2013%.1f%%)",
                          100 * plot_dat$Estimate,
                          100 * plot_dat$Lower,
                          100 * plot_dat$Upper)
colours <- stats::setNames(pal, plot_names)
g <- ggplot2::ggplot(plot_dat,
                     ggplot2::aes(y = Scenario, colour = Scenario)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = Lower, xend = Upper, yend = Scenario),
    linewidth = 0.9
  ) +
  ggplot2::geom_point(ggplot2::aes(x = Estimate), size = 3.1) +
  ggplot2::geom_text(
    ggplot2::aes(x = 1.04, label = Label),
    hjust = 0, size = 14 / ggplot2::.pt, family = font_family, show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Metric, ncol = 1) +
  ggplot2::scale_colour_manual(values = colours, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
  ggplot2::labs(
    title = "Patient clustering: sensitivity of pooled diagnostic accuracy",
    subtitle = "Bivariate REML; 95% confidence intervals",
    x = "Pooled estimate", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 14, base_family = font_family) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(colour = pal[1], face = "bold", size = 14),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(colour = pal[1], face = "bold", size = 14),
    axis.text = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(colour = "#222222", size = 14),
    plot.margin = ggplot2::margin(14, 285, 14, 12, unit = "pt")
  )
pdf_file <- file.path(output_dir, "Patient_Clustering_Sensitivity.pdf")
grDevices::cairo_pdf(pdf_file, width = 13.7, height = 10.5,
                     family = font_family)
print(g)
grDevices::dev.off()

report <- c(
  "PATIENT-CLUSTERING SENSITIVITY ANALYSIS",
  paste("Primary studies:", nrow(dat)),
  paste("Known repeated-observation studies:",
        sum(dat$Known_repeated_observations)),
  paste("Studies with unknown selected-cohort patient count:",
        sum(dat$Patient_count_unknown)),
  paste("Confirmed one-observation-per-patient studies:",
        sum(dat$Confirmed_one_per_patient)),
  "Main analysis: one selected four-cell table per study, with no duplicate model outputs.",
  "Within-study observations were treated as independent in the current primary analysis.",
  "Design effect: DE = 1 + (N / unique patients - 1) * assumed ICC.",
  "For known repeated-observation studies, DE multiplies both within-study logit variances.",
  "The same assumed ICC and overall N/patient multiplicity are applied to sensitivity and FPR; class-specific cluster sizes are unavailable.",
  "For unknown patient counts, no design effect can be calculated; these studies remain unadjusted in ICC scenarios.",
  "Sievert: at most five source patients; using five yields a minimum multiplicity and a lower-bound DE.",
  "ICC values are scenarios, not estimates from the included papers.",
  "The equal-cluster-size approximation cannot recover exact patient-level covariance and can understate uncertainty when cluster sizes vary.",
  "Excluded-study scenarios also change case mix and heterogeneity; interpret them as robustness checks.",
  "Exact cluster-robust intervals require patient-linked raw classifications or class-specific cluster summaries.",
  "CI_width_change_pp compares interval widths with the unadjusted bivariate baseline.",
  "Methods: PRISMA-DTA, BMJ 2020, doi:10.1136/bmj.m2632; Genders et al., Radiology 2012, doi:10.1148/radiol.12120509; Cochrane Handbook, Chapter 23."
)
writeLines(report,
           file.path(output_dir, "Patient_Clustering_Methods_and_Limitations.txt"),
           useBytes = TRUE)
print(results[, c("Scenario", "Metric", "Studies", "Estimate", "Lower",
                  "Upper", "CI_width_change_pp")], row.names = FALSE,
      digits = 4)
cat("Outputs saved to: ", output_dir, "\n", sep = "")
