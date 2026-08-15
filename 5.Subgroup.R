# ============================================================
# Seven subgroup meta-analysis forest plots
# Outcomes: Accuracy, NPV, PPV, Sensitivity, Specificity
# Font: Arial 14 pt
# Show BOTH fixed-effect and random-effects subgroup tests
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

# -------------------- 1. Install and load package --------------------
if (!requireNamespace("meta", quietly = TRUE)) {
  install.packages("meta", dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(meta)
})

# -------------------- 2. Paths --------------------
work_dir <- "C:/Users/30399/Desktop/LLM"

if (!dir.exists(work_dir)) {
  work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
}

if (!dir.exists(work_dir)) {
  stop("Working directory was not found.")
}

out_dir <- file.path(work_dir, "5.Subgroup")
input_file <- file.path(work_dir, "main.txt")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(input_file)) {
  stop("main.txt was not found: ", input_file)
}

# -------------------- 3. Colour palette --------------------
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

# -------------------- 4. Utility functions --------------------
safe_num <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("%", "", x)
  suppressWarnings(as.numeric(x))
}

# If the old PDF is open in Adobe/WPS, save a timestamped new PDF
get_available_pdf_name <- function(filename) {
  
  if (!file.exists(filename)) {
    return(filename)
  }
  
  deleted <- suppressWarnings(file.remove(filename))
  
  if (isTRUE(deleted) || !file.exists(filename)) {
    return(filename)
  }
  
  file_dir <- dirname(filename)
  file_base <- tools::file_path_sans_ext(basename(filename))
  file_ext <- tools::file_ext(filename)
  
  file.path(
    file_dir,
    paste0(
      file_base,
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".",
      file_ext
    )
  )
}

open_pdf <- function(filename, width, height) {
  
  filename <- get_available_pdf_name(filename)
  
  pdf_ok <- tryCatch(
    {
      grDevices::cairo_pdf(
        filename = filename,
        width = width,
        height = height,
        family = "Arial",
        pointsize = 14,
        onefile = TRUE
      )
      TRUE
    },
    error = function(e) FALSE
  )
  
  if (!pdf_ok) {
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      family = "Helvetica",
      pointsize = 14,
      onefile = TRUE
    )
  }
  
  return(filename)
}

# -------------------- 5. Read data --------------------
dat0 <- tryCatch(
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

if (is.null(dat0) || nrow(dat0) == 0) {
  dat0 <- read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    fileEncoding = "GB18030",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

if (ncol(dat0) < 15) {
  stop("main.txt has fewer than 15 columns.")
}

dat <- data.frame(
  Study       = as.character(dat0[[1]]),
  Year        = safe_num(dat0[[2]]),
  Country     = as.character(dat0[[3]]),
  Region      = as.character(dat0[[4]]),
  Organ       = as.character(dat0[[5]]),
  Modality    = as.character(dat0[[6]]),
  InputType   = as.character(dat0[[7]]),
  LLM         = as.character(dat0[[8]]),
  ModelFamily = as.character(dat0[[9]]),
  Design      = as.character(dat0[[10]]),
  N           = safe_num(dat0[[11]]),
  TP          = safe_num(dat0[[12]]),
  FP          = safe_num(dat0[[13]]),
  FN          = safe_num(dat0[[14]]),
  TN          = safe_num(dat0[[15]]),
  stringsAsFactors = FALSE
)

dat$Study <- gsub("^\"|\"$", "", trimws(dat$Study))
dat$Study[is.na(dat$Study) | dat$Study == ""] <- "Unknown"

# -------------------- 6. Seven subgroup variables --------------------

# (1) Publication year
dat$Year_group <- ifelse(
  dat$Year >= 2023 & dat$Year <= 2025,
  "2023-2025",
  ifelse(dat$Year == 2026, "2026", NA)
)

# (2) Sample size
median_n <- median(dat$N, na.rm = TRUE)

dat$Sample_size_group <- ifelse(
  dat$N <= median_n,
  paste0("≤ median (", median_n, ")"),
  paste0("> median (", median_n, ")")
)

# (3) Tumour organ
organ_text <- tolower(trimws(dat$Organ))

dat$Tumor_organ_group <- ifelse(
  grepl("thyroid|甲状腺", organ_text),
  "Thyroid",
  ifelse(
    grepl("breast|乳腺|乳房", organ_text),
    "Breast",
    "Other"
  )
)

# (4) Study design
design_text <- tolower(trimws(dat$Design))

dat$Study_design_group <- ifelse(
  grepl("prospective|前瞻", design_text),
  "Prospective",
  ifelse(
    grepl("retrospective|回顾", design_text),
    "Retrospective",
    NA
  )
)

# (5) Input type
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
      NA
    )
  )
)

# (6) Imaging modality
modality_text <- tolower(trimws(dat$Modality))

dat$Imaging_modality_group <- ifelse(
  grepl("ultrasound|超声", modality_text),
  "Ultrasound",
  "Other"
)

# (7) Country
country_text <- tolower(trimws(dat$Country))

dat$Country_group <- ifelse(
  grepl("china|中国", country_text),
  "China",
  "Other"
)

# Fixed order
dat$Year_group <- factor(
  dat$Year_group,
  levels = c("2023-2025", "2026")
)

dat$Sample_size_group <- factor(
  dat$Sample_size_group,
  levels = c(
    paste0("≤ median (", median_n, ")"),
    paste0("> median (", median_n, ")")
  )
)

dat$Tumor_organ_group <- factor(
  dat$Tumor_organ_group,
  levels = c("Thyroid", "Breast", "Other")
)

dat$Study_design_group <- factor(
  dat$Study_design_group,
  levels = c("Prospective", "Retrospective")
)

dat$Input_type_group <- factor(
  dat$Input_type_group,
  levels = c("Text", "Imaging", "Multimodal / multi-input")
)

dat$Imaging_modality_group <- factor(
  dat$Imaging_modality_group,
  levels = c("Ultrasound", "Other")
)

dat$Country_group <- factor(
  dat$Country_group,
  levels = c("China", "Other")
)

subgroups <- list(
  Year = list(variable = "Year_group", title = "Publication year"),
  Sample_size = list(variable = "Sample_size_group", title = "Sample size"),
  Tumor_organ = list(variable = "Tumor_organ_group", title = "Tumor organ"),
  Study_design = list(variable = "Study_design_group", title = "Study design"),
  Input_type = list(variable = "Input_type_group", title = "Input type"),
  Imaging_modality = list(
    variable = "Imaging_modality_group",
    title = "Imaging modality"
  ),
  Country = list(variable = "Country_group", title = "Country")
)

write.csv(
  dat,
  file.path(out_dir, "Subgroup_Meta_Analysis_Data_Used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 7. Define outcomes --------------------
make_events <- function(d, outcome) {
  
  if (outcome == "Sensitivity") {
    d$event <- d$TP
    d$total <- d$TP + d$FN
  }
  
  if (outcome == "Specificity") {
    d$event <- d$TN
    d$total <- d$TN + d$FP
  }
  
  if (outcome == "Accuracy") {
    d$event <- d$TP + d$TN
    d$total <- d$TP + d$FP + d$FN + d$TN
  }
  
  if (outcome == "PPV") {
    d$event <- d$TP
    d$total <- d$TP + d$FP
  }
  
  if (outcome == "NPV") {
    d$event <- d$TN
    d$total <- d$TN + d$FN
  }
  
  d <- d[
    is.finite(d$event) &
      is.finite(d$total) &
      d$total > 0 &
      d$event >= 0 &
      d$event <= d$total,
  ]
  
  d
}

outcomes <- c(
  "Accuracy",
  "NPV",
  "PPV",
  "Sensitivity",
  "Specificity"
)

# -------------------- 8. Forest plot function --------------------
run_forest <- function(d, outcome_name, subgroup_var, filename) {
  
  d <- make_events(d, outcome_name)
  
  d$Subgroup <- d[[subgroup_var]]
  d <- d[!is.na(d$Subgroup), ]
  
  if (nrow(d) < 2 || length(unique(d$Subgroup)) < 2) {
    return(NULL)
  }
  
  d$Subgroup <- droplevels(factor(d$Subgroup))
  
  n_studies <- nrow(d)
  n_groups <- length(levels(d$Subgroup))
  
  # Figure height includes dedicated space for the three final statistics
  plot_height <- max(
    14,
    7 + n_studies * 0.55 + n_groups * 1.35
  )
  
  meta_obj <- tryCatch(
    meta::metaprop(
      event = event,
      n = total,
      studlab = Study,
      data = d,
      sm = "PLOGIT",
      method = "GLMM",
      common = TRUE,
      random = TRUE,
      subgroup = Subgroup,
      method.random.ci = "HK",
      incr = 0.5,
      method.incr = "all"
    ),
    error = function(e) NULL
  )
  
  if (is.null(meta_obj)) {
    return(NULL)
  }
  
  actual_filename <- open_pdf(
    filename = filename,
    width = 19,
    height = plot_height
  )
  
  par(family = "Arial")
  
  meta::forest(
    meta_obj,
    
    common = TRUE,
    random = TRUE,
    overall = FALSE,
    
    overall.hetstat = TRUE,
    
    # Keep both subgroup difference tests
    test.subgroup.common = TRUE,
    test.subgroup.random = TRUE,
    
    label.common = "Common effect model",
    label.random = "Random effects model",
    
    label.test.subgroup.common =
      "Test for subgroup differences (fixed effect)",
    
    label.test.subgroup.random =
      "Test for subgroup differences (random effects)",
    
    print.subgroup.name = TRUE,
    
    # Critical settings: independent blank rows for bottom text
    addrows.hetstat = 7,
    addrows.below.overall = 4,
    
    sortvar = as.numeric(d$Subgroup),
    
    backtransf = TRUE,
    xlab = "Proportion",
    
    # Original colours and shapes unchanged
    col.square = pal[1],
    col.square.lines = pal[1],
    col.diamond = pal[2],
    col.diamond.lines = pal[1],
    col.diamond.common = pal[4],
    col.diamond.random = pal[2],
    col.subgroup = "black",
    
    # Arial, 14 pt
    fontsize = 14,
    fs.study = 14,
    fs.heading = 14,
    fs.subgroup = 14,
    fs.random = 14,
    fs.common = 14,
    fs.hetstat = 14,
    fs.test.subgroup = 14,
    
    # Keep normal study-row spacing
    spacing = 1,
    
    digits = 2
  )
  
  dev.off()
  
  list(
    meta = meta_obj,
    data = d,
    filename = actual_filename
  )
}

# -------------------- 9. Run forest plots --------------------
summary_results <- data.frame()
subgroup_p_results <- data.frame()

for (outcome_name in outcomes) {
  
  for (sg_name in names(subgroups)) {
    
    sg_var <- subgroups[[sg_name]]$variable
    sg_title <- subgroups[[sg_name]]$title
    
    output_pdf <- file.path(
      out_dir,
      paste0("Forest_", outcome_name, "_by_", sg_name, ".pdf")
    )
    
    result <- run_forest(
      d = dat,
      outcome_name = outcome_name,
      subgroup_var = sg_var,
      filename = output_pdf
    )
    
    if (is.null(result)) {
      next
    }
    
    meta_obj <- result$meta
    d_used <- result$data
    
    subgroup_p_results <- rbind(
      subgroup_p_results,
      data.frame(
        Outcome = outcome_name,
        Subgroup_type = sg_title,
        P_subgroup_common = meta_obj$pval.Q.b.common,
        P_subgroup_random = meta_obj$pval.Q.b.random,
        stringsAsFactors = FALSE
      )
    )
    
    for (lev in levels(d_used$Subgroup)) {
      
      sub_d <- d_used[d_used$Subgroup == lev, ]
      
      sub_meta <- tryCatch(
        meta::metaprop(
          event = event,
          n = total,
          studlab = Study,
          data = sub_d,
          sm = "PLOGIT",
          method = "GLMM",
          common = FALSE,
          random = TRUE,
          method.random.ci = "HK",
          incr = 0.5,
          method.incr = "all"
        ),
        error = function(e) NULL
      )
      
      if (is.null(sub_meta)) {
        next
      }
      
      summary_results <- rbind(
        summary_results,
        data.frame(
          Outcome = outcome_name,
          Subgroup_type = sg_title,
          Subgroup = as.character(lev),
          Studies = sub_meta$k,
          Pooled_random_percent = plogis(sub_meta$TE.random) * 100,
          Lower_95CI_percent = plogis(sub_meta$lower.random) * 100,
          Upper_95CI_percent = plogis(sub_meta$upper.random) * 100,
          I2_percent = sub_meta$I2,
          Tau2 = sub_meta$tau2,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# -------------------- 10. Export results --------------------
write.csv(
  summary_results,
  file.path(out_dir, "Subgroup_Meta_Analysis_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  subgroup_p_results,
  file.path(out_dir, "Subgroup_Difference_Test_Pvalues.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nSubgroup forest-plot analysis completed.\n")
cat("Output directory:\n", out_dir, "\n")













































# ============================================================
# Seven subgroup SROC plots
# Bivariate random-effects diagnostic meta-analysis
# ============================================================

# -------------------- 1. Packages --------------------
required_packages <- c(
  "metafor",
  "ggplot2",
  "Matrix",
  "scales"
)

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
  stop("Working directory was not found.")
}

input_file <- file.path(work_dir, "main.txt")

if (!file.exists(input_file)) {
  stop("main.txt was not found.")
}

output_dir <- file.path(work_dir, "5.Subgroup")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------- 3. Colour palette --------------------
pal <- c(
  "#3C3B8B",
  "#6656A6",
  "#917BBD",
  "#B19CCB",
  "#E5CBE1"
)

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

# -------------------- 6. Find columns --------------------
clean_name <- function(x) {
  tolower(
    gsub(
      "[^a-zA-Z0-9\u4e00-\u9fa5]",
      "",
      x
    )
  )
}

col_names_clean <- clean_name(colnames(dat_raw))

find_column <- function(
    candidates,
    fallback_position = NULL
) {
  
  candidates_clean <- clean_name(candidates)
  
  index <- which(
    col_names_clean %in% candidates_clean
  )
  
  if (length(index) > 0) {
    return(index[1])
  }
  
  if (
    !is.null(fallback_position) &&
    fallback_position <= ncol(dat_raw)
  ) {
    return(fallback_position)
  }
  
  return(NA_integer_)
}

author_col <- find_column(
  c(
    "第一作者",
    "First author",
    "Author",
    "FirstAuthor",
    "Study"
  ),
  fallback_position = 1
)

year_col <- find_column(
  c(
    "年份",
    "Year",
    "Publication year",
    "PublicationYear"
  ),
  fallback_position = 2
)

country_col <- find_column(
  c(
    "国家（统一）",
    "国家统一",
    "国家",
    "Country"
  ),
  fallback_position = 3
)

organ_col <- find_column(
  c(
    "肿瘤器官",
    "Tumor organ",
    "Tumour organ",
    "Organ"
  ),
  fallback_position = 4
)

imaging_col <- find_column(
  c(
    "影像方式（分组）",
    "影像方式分组",
    "影像方式",
    "Imaging modality"
  ),
  fallback_position = 5
)

input_col <- find_column(
  c(
    "输入大类",
    "Input category",
    "Input type"
  ),
  fallback_position = 6
)

design_col <- find_column(
  c(
    "研究设计（复核）",
    "研究设计复核",
    "研究设计",
    "Study design"
  ),
  fallback_position = 7
)

n_col <- find_column(
  c(
    "总样本N",
    "N",
    "Total N",
    "TotalN",
    "Sample size"
  ),
  fallback_position = 8
)

tp_col <- find_column(
  c(
    "TP",
    "True positive",
    "TruePositive"
  ),
  fallback_position = 9
)

fp_col <- find_column(
  c(
    "FP",
    "False positive",
    "FalsePositive"
  ),
  fallback_position = 10
)

fn_col <- find_column(
  c(
    "FN",
    "False negative",
    "FalseNegative"
  ),
  fallback_position = 11
)

tn_col <- find_column(
  c(
    "TN",
    "True negative",
    "TrueNegative"
  ),
  fallback_position = 12
)

required_columns <- c(
  author_col,
  year_col,
  country_col,
  organ_col,
  imaging_col,
  input_col,
  design_col,
  n_col,
  tp_col,
  fp_col,
  fn_col,
  tn_col
)

if (any(is.na(required_columns))) {
  stop(
    "One or more required columns could not be identified."
  )
}

# -------------------- 7. Prepare diagnostic data --------------------
dat <- data.frame(
  Author = as.character(dat_raw[[author_col]]),
  Year = as.character(dat_raw[[year_col]]),
  Country = as.character(dat_raw[[country_col]]),
  Organ = as.character(dat_raw[[organ_col]]),
  Imaging = as.character(dat_raw[[imaging_col]]),
  Input = as.character(dat_raw[[input_col]]),
  Design = as.character(dat_raw[[design_col]]),
  N = suppressWarnings(
    as.numeric(dat_raw[[n_col]])
  ),
  TP = suppressWarnings(
    as.numeric(dat_raw[[tp_col]])
  ),
  FP = suppressWarnings(
    as.numeric(dat_raw[[fp_col]])
  ),
  FN = suppressWarnings(
    as.numeric(dat_raw[[fn_col]])
  ),
  TN = suppressWarnings(
    as.numeric(dat_raw[[tn_col]])
  ),
  stringsAsFactors = FALSE
)

dat$Author[
  is.na(dat$Author) | dat$Author == ""
] <- "Unknown"

dat$Year[
  is.na(dat$Year) | dat$Year == ""
] <- ""

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

dat$Total_from_2x2 <- with(
  dat,
  TP + FP + FN + TN
)

dat <- dat[
  dat$Total_from_2x2 > 0 &
    (dat$TP + dat$FN) > 0 &
    (dat$FP + dat$TN) > 0,
]

dat$N[
  is.na(dat$N) | dat$N <= 0
] <- dat$Total_from_2x2[
  is.na(dat$N) | dat$N <= 0
]

if (nrow(dat) < 3) {
  stop("Fewer than three eligible studies were found.")
}

# -------------------- 8. Continuity correction --------------------
dat$Has_zero <- with(
  dat,
  TP == 0 |
    FP == 0 |
    FN == 0 |
    TN == 0
)

dat$TP_cc <- ifelse(
  dat$Has_zero,
  dat$TP + 0.5,
  dat$TP
)

dat$FP_cc <- ifelse(
  dat$Has_zero,
  dat$FP + 0.5,
  dat$FP
)

dat$FN_cc <- ifelse(
  dat$Has_zero,
  dat$FN + 0.5,
  dat$FN
)

dat$TN_cc <- ifelse(
  dat$Has_zero,
  dat$TN + 0.5,
  dat$TN
)

dat$Sensitivity <- with(
  dat,
  TP_cc / (TP_cc + FN_cc)
)

dat$Specificity <- with(
  dat,
  TN_cc / (TN_cc + FP_cc)
)

dat$FPR <- with(
  dat,
  FP_cc / (FP_cc + TN_cc)
)

dat$logit_sens <- qlogis(dat$Sensitivity)
dat$logit_fpr <- qlogis(dat$FPR)

dat$var_sens <- 1 / dat$TP_cc + 1 / dat$FN_cc
dat$var_fpr <- 1 / dat$FP_cc + 1 / dat$TN_cc

# -------------------- 9. Ellipse function --------------------
ellipse_points <- function(
    center,
    covariance_matrix,
    level = 0.95,
    n = 300
) {
  
  covariance_matrix <- as.matrix(
    covariance_matrix
  )
  
  covariance_matrix <- (
    covariance_matrix +
      t(covariance_matrix)
  ) / 2
  
  theta <- seq(
    0,
    2 * pi,
    length.out = n
  )
  
  circle <- rbind(
    cos(theta),
    sin(theta)
  )
  
  eig <- eigen(
    covariance_matrix,
    symmetric = TRUE
  )
  
  eig$values[
    eig$values < 0
  ] <- 0
  
  transform_matrix <- eig$vectors %*%
    diag(
      sqrt(eig$values),
      nrow = 2
    ) %*%
    t(eig$vectors)
  
  radius <- sqrt(
    qchisq(level, df = 2)
  )
  
  points <- t(
    matrix(
      center,
      nrow = 2,
      ncol = n
    ) +
      radius *
      transform_matrix %*%
      circle
  )
  
  data.frame(
    logit_sens = points[, 1],
    logit_fpr = points[, 2]
  )
}

# -------------------- 10. Fit one subgroup model --------------------
fit_one_subgroup <- function(
    data_sub,
    subgroup_name
) {
  
  if (nrow(data_sub) < 3) {
    return(NULL)
  }
  
  data_sub$Study_model <- paste0(
    data_sub$Study,
    "_",
    seq_len(nrow(data_sub))
  )
  
  long_dat <- do.call(
    rbind,
    lapply(
      seq_len(nrow(data_sub)),
      function(i) {
        
        data.frame(
          Study = data_sub$Study_model[i],
          Outcome = factor(
            c(
              "Sensitivity",
              "FPR"
            ),
            levels = c(
              "Sensitivity",
              "FPR"
            )
          ),
          yi = c(
            data_sub$logit_sens[i],
            data_sub$logit_fpr[i]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  V_list <- lapply(
    seq_len(nrow(data_sub)),
    function(i) {
      diag(
        c(
          data_sub$var_sens[i],
          data_sub$var_fpr[i]
        )
      )
    }
  )
  
  V_matrix <- as.matrix(
    Matrix::bdiag(V_list)
  )
  
  fit_bivariate <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ Outcome - 1,
      random = ~ Outcome | Study,
      struct = "UN",
      method = "REML",
      data = long_dat,
      sparse = FALSE
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(fit_bivariate)) {
    return(NULL)
  }
  
  beta_hat <- as.numeric(
    coef(fit_bivariate)
  )
  
  vcov_beta <- vcov(
    fit_bivariate
  )
  
  if (
    length(beta_hat) < 2 ||
    nrow(vcov_beta) < 2 ||
    any(is.na(beta_hat)) ||
    any(is.na(vcov_beta))
  ) {
    return(NULL)
  }
  
  mu_sens <- beta_hat[1]
  mu_fpr <- beta_hat[2]
  
  se_sens <- sqrt(
    max(vcov_beta[1, 1], 0)
  )
  
  se_fpr <- sqrt(
    max(vcov_beta[2, 2], 0)
  )
  
  pooled_sens <- plogis(mu_sens)
  pooled_fpr <- plogis(mu_fpr)
  pooled_spec <- 1 - pooled_fpr
  
  pooled_sens_lcl <- plogis(
    mu_sens - 1.96 * se_sens
  )
  
  pooled_sens_ucl <- plogis(
    mu_sens + 1.96 * se_sens
  )
  
  pooled_fpr_lcl <- plogis(
    mu_fpr - 1.96 * se_fpr
  )
  
  pooled_fpr_ucl <- plogis(
    mu_fpr + 1.96 * se_fpr
  )
  
  pooled_spec_lcl <- 1 - pooled_fpr_ucl
  pooled_spec_ucl <- 1 - pooled_fpr_lcl
  
  G_matrix <- fit_bivariate$G
  
  if (
    is.null(G_matrix) ||
    any(is.na(G_matrix)) ||
    nrow(G_matrix) != 2 ||
    ncol(G_matrix) != 2
  ) {
    G_matrix <- matrix(
      c(0.01, 0, 0, 0.01),
      nrow = 2
    )
  }
  
  G_matrix <- (
    G_matrix +
      t(G_matrix)
  ) / 2
  
  if (
    abs(G_matrix[2, 2]) < 1e-10
  ) {
    sroc_slope <- 0
  } else {
    sroc_slope <- G_matrix[1, 2] /
      G_matrix[2, 2]
  }
  
  fpr_grid <- seq(
    0.001,
    0.999,
    length.out = 1000
  )
  
  sroc_dat <- data.frame(
    FPR = fpr_grid,
    Sensitivity = plogis(
      mu_sens +
        sroc_slope *
        (
          qlogis(fpr_grid) -
            mu_fpr
        )
    ),
    Subgroup = subgroup_name,
    stringsAsFactors = FALSE
  )
  
  AUC <- sum(
    diff(sroc_dat$FPR) *
      (
        head(
          sroc_dat$Sensitivity,
          -1
        ) +
          tail(
            sroc_dat$Sensitivity,
            -1
          )
      ) / 2
  )
  
  confidence_ellipse <- ellipse_points(
    center = c(
      mu_sens,
      mu_fpr
    ),
    covariance_matrix = vcov_beta[1:2, 1:2]
  )
  
  confidence_ellipse$Sensitivity <- plogis(
    confidence_ellipse$logit_sens
  )
  
  confidence_ellipse$FPR <- plogis(
    confidence_ellipse$logit_fpr
  )
  
  confidence_ellipse$Subgroup <- subgroup_name
  
  prediction_ellipse <- ellipse_points(
    center = c(
      mu_sens,
      mu_fpr
    ),
    covariance_matrix = G_matrix
  )
  
  prediction_ellipse$Sensitivity <- plogis(
    prediction_ellipse$logit_sens
  )
  
  prediction_ellipse$FPR <- plogis(
    prediction_ellipse$logit_fpr
  )
  
  prediction_ellipse$Subgroup <- subgroup_name
  
  pooled_dat <- data.frame(
    FPR = pooled_fpr,
    Sensitivity = pooled_sens,
    Subgroup = subgroup_name,
    stringsAsFactors = FALSE
  )
  
  list(
    data = data_sub,
    sroc = sroc_dat,
    confidence = confidence_ellipse,
    prediction = prediction_ellipse,
    pooled = pooled_dat,
    pooled_sens = pooled_sens,
    pooled_spec = pooled_spec,
    pooled_sens_lcl = pooled_sens_lcl,
    pooled_sens_ucl = pooled_sens_ucl,
    pooled_spec_lcl = pooled_spec_lcl,
    pooled_spec_ucl = pooled_spec_ucl,
    AUC = AUC
  )
}

# -------------------- 11. Draw subgroup SROC plot --------------------
plot_subgroup_sroc <- function(
    data,
    subgroup_variable,
    subgroup_levels,
    subgroup_title,
    output_name
) {
  
  data$Subgroup <- as.character(
    data[[subgroup_variable]]
  )
  
  data <- data[
    !is.na(data$Subgroup) &
      data$Subgroup != "" &
      data$Subgroup %in% subgroup_levels,
  ]
  
  data$Subgroup <- factor(
    data$Subgroup,
    levels = subgroup_levels
  )
  
  model_list <- list()
  
  for (one_group in subgroup_levels) {
    
    one_data <- data[
      data$Subgroup == one_group,
      ,
      drop = FALSE
    ]
    
    one_model <- fit_one_subgroup(
      data_sub = one_data,
      subgroup_name = one_group
    )
    
    if (!is.null(one_model)) {
      model_list[[one_group]] <- one_model
    }
  }
  
  group_names <- subgroup_levels[
    subgroup_levels %in% names(model_list)
  ]
  
  if (length(group_names) == 0) {
    cat(
      "\nNo subgroup model could be fitted: ",
      subgroup_title,
      "\n"
    )
    return(invisible(NULL))
  }
  
  # -------------------- Colours and line types --------------------
  line_types <- c(
    "solid",
    "22",
    "11",
    "42",
    "F2"
  )
  
  group_colours <- setNames(
    pal[seq_along(group_names)],
    group_names
  )
  
  group_linetypes <- setNames(
    line_types[seq_along(group_names)],
    group_names
  )
  
  # -------------------- Data for plotting --------------------
  all_studies <- do.call(
    rbind,
    lapply(
      group_names,
      function(x) {
        z <- model_list[[x]]$data
        z$Subgroup <- factor(
          x,
          levels = group_names
        )
        z
      }
    )
  )
  
  sroc_all <- do.call(
    rbind,
    lapply(
      group_names,
      function(x) {
        z <- model_list[[x]]$sroc
        z$Subgroup <- factor(
          x,
          levels = group_names
        )
        z
      }
    )
  )
  
  confidence_all <- do.call(
    rbind,
    lapply(
      group_names,
      function(x) {
        z <- model_list[[x]]$confidence
        z$Subgroup <- factor(
          x,
          levels = group_names
        )
        z
      }
    )
  )
  
  prediction_all <- do.call(
    rbind,
    lapply(
      group_names,
      function(x) {
        z <- model_list[[x]]$prediction
        z$Subgroup <- factor(
          x,
          levels = group_names
        )
        z
      }
    )
  )
  
  pooled_all <- do.call(
    rbind,
    lapply(
      group_names,
      function(x) {
        z <- model_list[[x]]$pooled
        z$Subgroup <- factor(
          x,
          levels = group_names
        )
        z
      }
    )
  )
  
  # -------------------- Manual subgroup legend --------------------
  # 放在结果框下方，避免与结果内容重合
  n_groups <- length(group_names)
  
  legend_y <- seq(
    from = 0.36,
    to = 0.36 - 0.075 * (n_groups - 1),
    length.out = n_groups
  )
  
  subgroup_legend <- data.frame(
    Subgroup = factor(
      group_names,
      levels = group_names
    ),
    x_start = 1.045,
    x_end = 1.180,
    x_text = 1.205,
    y = legend_y,
    stringsAsFactors = FALSE
  )
  
  # -------------------- Result label --------------------
  result_label <- paste(
    vapply(
      group_names,
      function(x) {
        
        z <- model_list[[x]]
        
        paste0(
          x,
          "\n",
          "Studies: ",
          nrow(z$data),
          "\n",
          "Sensitivity: ",
          percent(
            z$pooled_sens,
            accuracy = 0.1
          ),
          " (",
          percent(
            z$pooled_sens_lcl,
            accuracy = 0.1
          ),
          "–",
          percent(
            z$pooled_sens_ucl,
            accuracy = 0.1
          ),
          ")\n",
          "Specificity: ",
          percent(
            z$pooled_spec,
            accuracy = 0.1
          ),
          " (",
          percent(
            z$pooled_spec_lcl,
            accuracy = 0.1
          ),
          "–",
          percent(
            z$pooled_spec_ucl,
            accuracy = 0.1
          ),
          ")\n",
          "AUC: ",
          sprintf(
            "%.3f",
            z$AUC
          )
        )
      },
      character(1)
    ),
    collapse = "\n\n"
  )
  
  result_label_data <- data.frame(
    x = 1.045,
    y = 0.985,
    label = result_label
  )
  
  # -------------------- Plot --------------------
  sroc_plot <- ggplot() +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.7,
      color = "grey65"
    ) +
    
    geom_polygon(
      data = prediction_all,
      aes(
        x = FPR,
        y = Sensitivity,
        group = Subgroup,
        fill = Subgroup,
        color = Subgroup
      ),
      alpha = 0.20,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    
    geom_polygon(
      data = confidence_all,
      aes(
        x = FPR,
        y = Sensitivity,
        group = Subgroup,
        fill = Subgroup,
        color = Subgroup
      ),
      alpha = 0.42,
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    
    geom_line(
      data = sroc_all,
      aes(
        x = FPR,
        y = Sensitivity,
        group = Subgroup,
        color = Subgroup,
        linetype = Subgroup
      ),
      linewidth = 1.5,
      lineend = "butt",
      show.legend = FALSE
    ) +
    
    geom_point(
      data = all_studies,
      aes(
        x = FPR,
        y = Sensitivity,
        size = N,
        fill = Subgroup
      ),
      shape = 21,
      color = "white",
      stroke = 0.65,
      alpha = 0.92
    ) +
    
    geom_point(
      data = pooled_all,
      aes(
        x = FPR,
        y = Sensitivity,
        fill = Subgroup,
        color = Subgroup
      ),
      shape = 23,
      size = 5.8,
      stroke = 1.1,
      show.legend = FALSE
    ) +
    
    # 亚组标题
    annotate(
      "text",
      x = 1.045,
      y = 0.42,
      label = subgroup_title,
      hjust = 0,
      vjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = text_size_14,
      color = "#222222"
    ) +
    
    # 亚组线型
    geom_segment(
      data = subgroup_legend,
      aes(
        x = x_start,
        xend = x_end,
        y = y,
        yend = y,
        color = Subgroup,
        linetype = Subgroup
      ),
      linewidth = 1.6,
      lineend = "butt",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    # 亚组名称
    geom_text(
      data = subgroup_legend,
      aes(
        x = x_text,
        y = y,
        label = Subgroup
      ),
      family = "Arial",
      size = text_size_14,
      hjust = 0,
      vjust = 0.5,
      color = "#222222",
      inherit.aes = FALSE
    ) +
    
    # 结果框
    geom_label(
      data = result_label_data,
      aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0,
      vjust = 1,
      family = "Arial",
      fontface = "plain",
      size = text_size_14,
      color = "#222222",
      fill = alpha("white", 0.96),
      label.size = 0.35,
      label.r = grid::unit(
        0.18,
        "lines"
      ),
      label.padding = grid::unit(
        0.35,
        "lines"
      ),
      inherit.aes = FALSE
    ) +
    
    scale_color_manual(
      values = group_colours,
      limits = group_names,
      breaks = group_names,
      guide = "none"
    ) +
    
    scale_fill_manual(
      values = group_colours,
      limits = group_names,
      breaks = group_names,
      guide = "none"
    ) +
    
    scale_linetype_manual(
      values = group_linetypes,
      limits = group_names,
      breaks = group_names,
      guide = "none"
    ) +
    
    scale_size_continuous(
      range = c(2.8, 9),
      name = "Sample size"
    ) +
    
    scale_x_continuous(
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = percent_format(
        accuracy = 1
      ),
      expand = c(
        0.01,
        0.01
      )
    ) +
    
    scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = percent_format(
        accuracy = 1
      ),
      expand = c(
        0.01,
        0.01
      )
    ) +
    
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      clip = "off"
    ) +
    
    labs(
      title = "Summary Receiver Operating Characteristic (SROC) Curve",
      subtitle = paste0(
        "Bivariate random-effects diagnostic meta-analysis: ",
        subgroup_title
      ),
      x = "False-positive rate (1 − Specificity)",
      y = "Sensitivity",
      caption = paste0(
        "Lines: subgroup-specific SROC curves; ",
        "diamonds: subgroup pooled estimates; ",
        "dark shaded areas: 95% confidence regions; ",
        "light shaded areas: 95% prediction regions."
      )
    ) +
    
    theme_classic(
      base_family = "Arial",
      base_size = 14
    ) +
    
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
        size = 14
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
      
      axis.ticks.length = grid::unit(
        0.20,
        "cm"
      ),
      
      # Sample size 图例放到左下角
      legend.position = c(
        0.05,
        0.06
      ),
      
      legend.justification = c(
        0,
        0
      ),
      
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
        fill = alpha(
          "white",
          0.94
        ),
        color = "grey70",
        linewidth = 0.45
      ),
      
      legend.key = element_rect(
        fill = "white",
        color = NA
      ),
      
      legend.key.size = grid::unit(
        0.65,
        "cm"
      ),
      
      plot.caption = element_text(
        family = "Arial",
        size = 12,
        hjust = 0,
        color = "grey25",
        margin = margin(t = 12)
      ),
      
      # 右侧保留亚组图例空间
      plot.margin = margin(
        15,
        280,
        15,
        18
      )
    ) +
    
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          fill = pal[3],
          colour = "white",
          stroke = 0.65,
          alpha = 0.92
        )
      )
    )
  
  pdf_file <- file.path(
    output_dir,
    paste0(
      "Subgroup_SROC_",
      output_name,
      ".pdf"
    )
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
  
  cat(
    "\nCompleted: ",
    subgroup_title,
    "\nPDF saved to:\n",
    pdf_file,
    "\n"
  )
  
  invisible(sroc_plot)
}

# -------------------- 12. Define subgroup variables --------------------

# (1) Publication year
year_numeric <- suppressWarnings(
  as.numeric(dat$Year)
)

dat$Year_group <- NA_character_

dat$Year_group[
  year_numeric >= 2023 &
    year_numeric <= 2025
] <- "2023–2025"

dat$Year_group[
  year_numeric == 2026
] <- "2026"

# (2) Sample size
median_N <- median(
  dat$N,
  na.rm = TRUE
)

dat$Sample_size_group <- ifelse(
  dat$N <= median_N,
  paste0(
    "≤ Median (",
    round(median_N),
    ")"
  ),
  paste0(
    "> Median (",
    round(median_N),
    ")"
  )
)

# (3) Tumour organ
dat$Organ_group <- ifelse(
  grepl(
    "thyroid|甲状腺",
    dat$Organ,
    ignore.case = TRUE
  ),
  "Thyroid",
  ifelse(
    grepl(
      "breast|乳腺",
      dat$Organ,
      ignore.case = TRUE
    ),
    "Breast",
    "Other"
  )
)

# (4) Study design
dat$Design_group <- ifelse(
  grepl(
    "prospective|前瞻",
    dat$Design,
    ignore.case = TRUE
  ),
  "Prospective",
  ifelse(
    grepl(
      "retrospective|回顾",
      dat$Design,
      ignore.case = TRUE
    ),
    "Retrospective",
    NA_character_
  )
)

# (5) Input category
dat$Input_group <- ifelse(
  grepl(
    "multi-input|multi input|多输入",
    dat$Input,
    ignore.case = TRUE
  ),
  "multi-input",
  ifelse(
    grepl(
      "multimodal|多模态",
      dat$Input,
      ignore.case = TRUE
    ),
    "Multimodal",
    ifelse(
      grepl(
        "imaging|image|影像",
        dat$Input,
        ignore.case = TRUE
      ),
      "Imaging",
      ifelse(
        grepl(
          "text|文本",
          dat$Input,
          ignore.case = TRUE
        ),
        "Text",
        NA_character_
      )
    )
  )
)

# (6) Imaging modality
dat$Imaging_group <- ifelse(
  grepl(
    "ultrasound|超声",
    dat$Imaging,
    ignore.case = TRUE
  ),
  "Ultrasound",
  "Other"
)

# (7) Country
dat$Country_group <- ifelse(
  grepl(
    "china|中国|mainland china",
    dat$Country,
    ignore.case = TRUE
  ),
  "China",
  "Other"
)

# -------------------- 13. Generate seven subgroup SROC plots --------------------

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Year_group",
  subgroup_levels = c(
    "2023–2025",
    "2026"
  ),
  subgroup_title = "Publication year",
  output_name = "Publication_year"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Sample_size_group",
  subgroup_levels = c(
    paste0(
      "≤ Median (",
      round(median_N),
      ")"
    ),
    paste0(
      "> Median (",
      round(median_N),
      ")"
    )
  ),
  subgroup_title = "Sample size",
  output_name = "Sample_size"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Organ_group",
  subgroup_levels = c(
    "Thyroid",
    "Breast",
    "Other"
  ),
  subgroup_title = "Tumour organ",
  output_name = "Tumour_organ"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Design_group",
  subgroup_levels = c(
    "Prospective",
    "Retrospective"
  ),
  subgroup_title = "Study design",
  output_name = "Study_design"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Input_group",
  subgroup_levels = c(
    "Text",
    "Imaging",
    "Multimodal",
    "multi-input"
  ),
  subgroup_title = "Input category",
  output_name = "Input_category"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Imaging_group",
  subgroup_levels = c(
    "Ultrasound",
    "Other"
  ),
  subgroup_title = "Imaging modality",
  output_name = "Imaging_modality"
)

plot_subgroup_sroc(
  data = dat,
  subgroup_variable = "Country_group",
  subgroup_levels = c(
    "China",
    "Other"
  ),
  subgroup_title = "Country",
  output_name = "Country"
)

cat(
  "\nAll subgroup SROC analyses completed.\n",
  "Output directory:\n",
  output_dir,
  "\n"
)