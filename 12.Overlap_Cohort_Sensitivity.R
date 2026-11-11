# Prefer the analysis-local, version-controlled R library when available.
local_analysis_lib <- file.path(Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis"), "R_packages")
if (dir.exists(local_analysis_lib)) .libPaths(c(local_analysis_lib, .libPaths()))

# ============================================================
# Sensitivity analysis for overlapping Chen cohorts in Table 2
# Primary: Chen 2024 [N=116] (32, 48, 9, 27), excluding red Chen 2025 [N=106]
# Alternative: Chen 2025 [N=106] (26, 26, 11, 43), excluding Chen 2024 [N=116]
# Both analyses contain 20 independent study rows.
# ============================================================

# Some Windows Rscript sessions start in the C locale and cannot resolve
# the Chinese directory name in the project path.
if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C") {
  invisible(suppressWarnings(
    Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936")
  ))
}

default_work_dir <- Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis")
Encoding(default_work_dir) <- "UTF-8"
work_dir <- Sys.getenv("ANALYSIS_WORK_DIR", unset = default_work_dir)
input_file <- file.path(work_dir, "table 2.txt")
main_file <- file.path(work_dir, "main.txt")
output_dir <- Sys.getenv(
  "ANALYSIS_OUTPUT_DIR",
  unset = file.path(work_dir, "12.sensitivity analysis")
)

if (!file.exists(input_file) || !file.exists(main_file)) {
  stop("Required input file is missing: table 2.txt or main.txt")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("metafor", "Matrix", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install the required R package before running: ", pkg)
  }
}

source(file.path(work_dir, "0.Verify_Main_Data.R"), encoding = "UTF-8")
assert_main_data(main_file)

# Keep the colour palette used in scripts 2-11.
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")
scenario_names <- c(
  Primary = "Primary: Chen 2024 [N=116]",
  Replacement = "Sensitivity: Chen 2025 [N=106]"
)
scenario_colours <- c(
  "Primary: Chen 2024 [N=116]" = pal[1],
  "Sensitivity: Chen 2025 [N=106]" = pal[3]
)

# The first line is a title. The second line contains the column names.
raw <- read.delim(
  input_file, skip = 1, header = TRUE, sep = "\t",
  quote = "", comment.char = "", fill = TRUE,
  fileEncoding = "GB18030", check.names = FALSE,
  stringsAsFactors = FALSE
)
required <- c("First author", "Online Year", "TP", "FP", "FN", "TN")
if (!all(required %in% names(raw))) {
  stop("table 2.txt is missing: ", paste(setdiff(required, names(raw)), collapse = ", "))
}

dat <- data.frame(
  Author = trimws(gsub('"', '', as.character(raw[["First author"]]), fixed = TRUE)),
  Year = suppressWarnings(as.integer(raw[["Online Year"]])),
  TP = suppressWarnings(as.numeric(raw[["TP"]])),
  FP = suppressWarnings(as.numeric(raw[["FP"]])),
  FN = suppressWarnings(as.numeric(raw[["FN"]])),
  TN = suppressWarnings(as.numeric(raw[["TN"]])),
  stringsAsFactors = FALSE
)
dat$Author <- gsub("†", "", dat$Author, fixed = TRUE)
# Exclude text footnotes, which have no year or four-cell table.
dat <- dat[complete.cases(dat[, c("Year", "TP", "FP", "FN", "TN")]), ]
if (nrow(dat) != 21L) stop("Expected 21 Table 2 study rows; found ", nrow(dat))
if (any(as.matrix(dat[, c("TP", "FP", "FN", "TN")]) < 0) ||
    any(as.matrix(dat[, c("TP", "FP", "FN", "TN")]) !=
        floor(as.matrix(dat[, c("TP", "FP", "FN", "TN")])))) {
  stop("A four-cell count is negative or noninteger")
}

count_key <- function(x) paste(x$TP, x$FP, x$FN, x$TN, sep = ":")
find_chen <- function(year, counts) {
  which(
    grepl("^Chen,\\s*Ziman#?$", dat$Author) & dat$Year == year &
      count_key(dat) == paste(counts, collapse = ":")
  )
}
chen_primary <- find_chen(2024L, c(32, 48, 9, 27))
chen_red <- find_chen(2025L, c(26, 26, 11, 43))
if (length(chen_primary) != 1L || length(chen_red) != 1L) {
  stop("The two overlapping Chen rows were not uniquely identified")
}

# Check that the primary selection is the same set of four-cell tables
# as main.txt. Both source files must use the same selected indicator direction.
main <- read.delim(
  main_file, header = TRUE, sep = "\t", fileEncoding = "GB18030",
  check.names = FALSE, stringsAsFactors = FALSE
)
main_counts <- data.frame(
  TP = as.numeric(main[["TP"]]), FP = as.numeric(main[["FP"]]),
  FN = as.numeric(main[["FN"]]), TN = as.numeric(main[["TN"]])
)
base_raw <- dat[-chen_red, ]
if (nrow(main_counts) != 20L ||
    !identical(sort(paste(gsub("[^a-z0-9]", "", tolower(main[[1]])), main[[2]], count_key(main_counts))),
               sort(paste(gsub("[^a-z0-9]", "", tolower(base_raw$Author)), base_raw$Year, count_key(base_raw))))) {
  stop("The primary Table 2 counts do not match main.txt; resolve inputs first")
}

# The selected four-cell table is TP=8, FP=9, FN=11, TN=20.
mehta <- which(grepl("^Mehta,", dat$Author))
if (length(mehta) != 1L || count_key(dat[mehta, ]) != "8:9:11:20") {
  stop("Mehta must be melanoma-positive: TP/FP/FN/TN = 8/9/11/20")
}

dat$Study <- paste0(gsub("#", "", dat$Author), " (", dat$Year, ")")
dat$Study[c(chen_primary, chen_red)] <- paste0(dat$Study[c(chen_primary, chen_red)], " [N=", rowSums(dat[c(chen_primary, chen_red), c("TP", "FP", "FN", "TN")]), "]")
dat$Total <- rowSums(dat[, c("TP", "FP", "FN", "TN")])
dat$Selection <- "Both scenarios"
dat$Selection[chen_primary] <- "Primary only"
dat$Selection[chen_red] <- "Replacement only (red Table 2 row)"
dat$Data_note <- "Table 2 counts unchanged"
dat$Data_note[mehta] <- "Table 2 counts synchronized"
if (any(dat$Total <= 0)) stop("A study has no observations")

primary <- dat[-chen_red, ]
replacement <- dat[-chen_primary, ]
if (nrow(primary) != 20L || nrow(replacement) != 20L ||
    anyDuplicated(primary$Study) || anyDuplicated(replacement$Study)) {
  stop("Each scenario must contain 20 uniquely labelled studies")
}

adjust_zero_cells <- function(x) {
  a <- as.matrix(x[, c("TP", "FP", "FN", "TN")])
  zero_row <- rowSums(a == 0) > 0
  a[zero_row, ] <- a[zero_row, ] + 0.5
  a
}

# Bivariate random-effects model: logit sensitivity and logit false-positive
# rate; independent conditional binomial sampling within each study;
# unstructured between-study covariance; REML estimation.
fit_bivariate <- function(x, scenario) {
  a <- adjust_zero_cells(x)
  logit_sens <- log(a[, "TP"] / a[, "FN"])
  logit_fpr <- log(a[, "FP"] / a[, "TN"])
  var_sens <- 1 / a[, "TP"] + 1 / a[, "FN"]
  var_fpr <- 1 / a[, "FP"] + 1 / a[, "TN"]
  long <- data.frame(
    Study = factor(rep(x$Study, each = 2)),
    Outcome = factor(rep(c("Sensitivity", "FPR"), nrow(x)),
                     levels = c("Sensitivity", "FPR")),
    yi = as.vector(t(cbind(logit_sens, logit_fpr)))
  )
  V <- Matrix::Diagonal(x = as.vector(t(cbind(var_sens, var_fpr))))
  fit <- metafor::rma.mv(
    yi = yi, V = V, mods = ~ Outcome - 1,
    random = ~ Outcome | Study, struct = "UN", method = "REML",
    data = long, sparse = FALSE
  )
  beta <- as.numeric(stats::coef(fit))
  se <- sqrt(diag(stats::vcov(fit)))
  z <- stats::qnorm(0.975)
  data.frame(
    Scenario = scenario, Method = "Bivariate REML",
    Metric = c("Sensitivity", "Specificity"),
    Estimate = c(stats::plogis(beta[1]), 1 - stats::plogis(beta[2])),
    Lower = c(stats::plogis(beta[1] - z * se[1]),
              1 - stats::plogis(beta[2] + z * se[2])),
    Upper = c(stats::plogis(beta[1] + z * se[1]),
              1 - stats::plogis(beta[2] - z * se[2])),
    Studies = nrow(x)
  )
}

# Supporting five-metric analysis matches the main forest-plot principle:
# logit inverse-variance random-effects meta-analysis, REML tau-squared.
fit_univariate <- function(x, scenario) {
  metric_pairs <- list(
    Sensitivity = c("TP", "FN"),
    Specificity = c("TN", "FP"),
    Accuracy = c("TP+TN", "FP+FN"),
    PPV = c("TP", "FP"),
    NPV = c("TN", "FN")
  )
  rows <- lapply(names(metric_pairs), function(metric) {
    pair <- metric_pairs[[metric]]
    # Correct a zero event/failure for this outcome only, as in script 2.
    raw <- as.matrix(x[, c("TP", "FP", "FN", "TN")])
    get_raw <- function(expr) {
      parts <- strsplit(expr, "+", fixed = TRUE)[[1]]
      rowSums(raw[, parts, drop = FALSE])
    }
    event <- get_raw(pair[1]); failure <- get_raw(pair[2])
    if (any(event + failure <= 0)) stop("Zero denominator in ", metric)
    zero <- event == 0 | failure == 0
    event[zero] <- event[zero] + 0.5
    failure[zero] <- failure[zero] + 0.5
    fit <- metafor::rma.uni(
      yi = log(event / failure),
      vi = 1 / event + 1 / failure,
      method = "REML"
    )
    beta <- as.numeric(stats::coef(fit))
    se <- sqrt(as.numeric(stats::vcov(fit)))
    z <- stats::qnorm(0.975)
    data.frame(
      Scenario = scenario, Method = "Univariate REML",
      Metric = metric,
      Estimate = stats::plogis(beta),
      Lower = stats::plogis(beta - z * se),
      Upper = stats::plogis(beta + z * se),
      Studies = nrow(x)
    )
  })
  do.call(rbind, rows)
}

results <- rbind(
  fit_bivariate(primary, scenario_names[["Primary"]]),
  fit_bivariate(replacement, scenario_names[["Replacement"]]),
  fit_univariate(primary, scenario_names[["Primary"]]),
  fit_univariate(replacement, scenario_names[["Replacement"]])
)
if (any(!is.finite(as.matrix(results[, c("Estimate", "Lower", "Upper")]))) ||
    any(results$Lower > results$Estimate) ||
    any(results$Upper < results$Estimate)) {
  stop("Invalid pooled estimate or confidence interval")
}

key <- paste(results$Method, results$Metric)
p <- results[results$Scenario == scenario_names[["Primary"]], ]
s <- results[results$Scenario == scenario_names[["Replacement"]], ]
s <- s[match(paste(p$Method, p$Metric), paste(s$Method, s$Metric)), ]
if (anyNA(s$Metric)) stop("Scenario outcomes do not align")
comparison <- data.frame(
  Method = p$Method, Metric = p$Metric, Studies_each = p$Studies,
  Primary = p$Estimate, Primary_Lower = p$Lower, Primary_Upper = p$Upper,
  Replacement = s$Estimate,
  Replacement_Lower = s$Lower, Replacement_Upper = s$Upper,
  Difference_pp = 100 * (s$Estimate - p$Estimate),
  CI_width_change_pp = 100 * ((s$Upper - s$Lower) - (p$Upper - p$Lower))
)

write.csv(dat, file.path(output_dir, "Overlap_Sensitivity_Data_Used.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(results, file.path(output_dir, "Overlap_Sensitivity_All_Estimates.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(comparison, file.path(output_dir, "Overlap_Sensitivity_Comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

make_plot <- function(method, filename, height) {
  d <- results[results$Method == method, ]
  metrics <- if (method == "Bivariate REML") {
    c("Sensitivity", "Specificity")
  } else {
    c("Sensitivity", "Specificity", "Accuracy", "PPV", "NPV")
  }
  d$y <- match(d$Metric, rev(metrics)) +
    ifelse(d$Scenario == scenario_names[["Primary"]], -0.13, 0.13)
  d$Value_label <- sprintf(
    "%.1f%% (%.1f%%\u2013%.1f%%)",
    100 * d$Estimate, 100 * d$Lower, 100 * d$Upper
  )
  g <- ggplot2::ggplot(d, ggplot2::aes(colour = Scenario)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = Lower, xend = Upper, y = y, yend = y),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Estimate, y = y), size = 3.1
    ) +
    # Numeric estimates sit in a dedicated column beyond the plotting panel.
    ggplot2::geom_text(
      ggplot2::aes(x = 1.045, y = y, label = Value_label),
      hjust = 0, size = 14 / ggplot2::.pt, family = "Arial", show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text", x = 1.045, y = length(metrics) + 0.48,
      label = "Estimate (95% CI)", hjust = 0,
      size = 14 / ggplot2::.pt, family = "Arial", fontface = "bold", colour = pal[1]
    ) +
    ggplot2::scale_colour_manual(values = scenario_colours) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(metrics), labels = rev(metrics),
      limits = c(0.5, length(metrics) + 0.65)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
    ggplot2::labs(
      title = "Overlapping Chen cohort: replacement sensitivity analysis",
      subtitle = paste(method, "| 20 studies in each scenario"),
      x = "Pooled estimate (95% CI)", y = NULL, colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(colour = pal[1], face = "bold", size = 14),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(colour = "#222222", size = 14),
      legend.text = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(12, 290, 12, 12, unit = "pt")
    )
  grDevices::cairo_pdf(file.path(output_dir, filename),
                       width = 13.5, height = height, family = "Arial")
  print(g)
  grDevices::dev.off()
}
make_plot("Bivariate REML", "Overlap_Sensitivity_Bivariate.pdf", 5.4)
make_plot("Univariate REML", "Overlap_Sensitivity_Five_Metrics.pdf", 7.2)

report <- c(
  "OVERLAPPING-COHORT SENSITIVITY ANALYSIS",
  "Data: table 2.txt (21 rows; one overlapping Chen row per scenario)",
  "Primary: Chen 2024 [N=116], TP/FP/FN/TN = 32/48/9/27; 20 studies",
  "Replacement: Chen 2025 [N=106], TP/FP/FN/TN = 26/26/11/43; 20 studies",
  "The Chen 2026 cohort remains in both scenarios.",
  "The primary Table 2 four-cell tables match main.txt exactly.",
  "Zero cells: add 0.5 to all four cells of affected studies for bivariate logit estimation.",
  "Univariate metrics: add 0.5 to a zero event/failure pair for that outcome.",
  "All reported intervals are 95% Wald intervals on the logit scale.",
  "Difference_pp = replacement minus primary in percentage points.",
  "CI_width_change_pp = replacement CI width minus primary CI width.",
  "These intervals do not account for multiple observations within a patient."
)
writeLines(report, file.path(output_dir, "Overlap_Sensitivity_Methods_and_Audit.txt"),
           useBytes = TRUE)
print(comparison, row.names = FALSE, digits = 4)
cat("Outputs saved to: ", output_dir, "\n", sep = "")
