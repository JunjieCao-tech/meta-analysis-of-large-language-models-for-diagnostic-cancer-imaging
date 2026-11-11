# Prefer the analysis-local, version-controlled R library when available.
local_analysis_lib <- file.path(Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis"), "R_packages")
if (dir.exists(local_analysis_lib)) .libPaths(c(local_analysis_lib, .libPaths()))

# =========================================================
# QUADAS-AI publication-style plots
# Font: Arial, 14 pt
# =========================================================

# ---------- 0. Install and load packages ----------
pkgs <- c("dplyr", "tidyr", "ggplot2", "stringr")

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Required packages are missing: ", paste(missing_pkgs, collapse = ", "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ---------- 1. File path ----------
work_dir <- Sys.getenv("ANALYSIS_WORK_DIR", unset = "E:/LLM/投稿版本/npj/R1/Statistical analysis")
data_file <- file.path(work_dir, "main.txt")
output_dir <- file.path(work_dir, "1.QUADAS_AI")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 2. Read data ----------
dat <- read.delim(
  data_file,
  header = TRUE,
  sep = "\t",
  fileEncoding = "GB18030",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Select assessment columns by name. Trailing metadata are not QUADAS ratings.
rating_columns <- c("RoB_患者选择", "RoB_待评试验(LLM)", "RoB_金标准", "RoB_流程与时序",
                    "App_患者选择", "App_待评试验", "App_金标准")
required_columns <- c("第一作者", "年份", rating_columns, "Study_ID", "Analysis_role")
if (!all(required_columns %in% names(dat))) stop("Missing required assessment columns")
dat <- dat[dat$Analysis_role == "Primary", , drop = FALSE]
if (nrow(dat) != 20L || anyDuplicated(dat$Study_ID)) stop("Expected 20 distinct primary studies")
valid <- c("low", "high", "unclear")
for (field in rating_columns) {
  if (anyNA(dat[[field]]) || any(!tolower(trimws(dat[[field]])) %in% valid))
    stop("Missing or invalid QUADAS rating in column: ", field)
}
quadas <- dat[, c("第一作者", "年份", rating_columns), drop = FALSE]

names(quadas) <- c(
  "Author", "Year",
  "RoB_PatientSelection",
  "RoB_IndexTest",
  "RoB_ReferenceStandard",
  "RoB_FlowTiming",
  "App_PatientSelection",
  "App_IndexTest",
  "App_ReferenceStandard"
)

quadas <- quadas %>%
  mutate(
    Study = paste0(Author, " (", Year, ")"),
    Study = make.unique(Study)
  )

study_order <- quadas$Study
assessment_template <- data.frame(
  study = study_order,
  患者选择 = "", 指标试验_AI模型 = "", 参考标准 = "", 病例流程与时间间隔 = "",
  患者选择_适用性 = "", 指标试验_适用性 = "", 参考标准_适用性 = "",
  check.names = FALSE
)
write.csv(assessment_template, file.path(output_dir, "QUADAS_AI_assessment_template.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------- 3. Reshape and clean QUADAS-AI ratings ----------
quadas_long <- quadas %>%
  pivot_longer(
    cols = starts_with("RoB_") | starts_with("App_"),
    names_to = "DomainCode",
    values_to = "Judgement"
  ) %>%
  mutate(
    Judgement = str_to_lower(str_trim(as.character(Judgement))),
    
    Judgement = case_when(
      Judgement == "low" ~ "Low",
      Judgement == "high" ~ "High",
      Judgement == "unclear" ~ "Unclear",
      TRUE ~ NA_character_
    ),
    
    Assessment = case_when(
      startsWith(DomainCode, "RoB_") ~ "Risk of Bias",
      startsWith(DomainCode, "App_") ~ "Applicability Concerns"
    ),
    
    Domain = case_when(
      str_detect(DomainCode, "PatientSelection") ~ "Patient Selection",
      str_detect(DomainCode, "IndexTest") ~ "Index Test",
      str_detect(DomainCode, "ReferenceStandard") ~ "Reference Standard",
      str_detect(DomainCode, "FlowTiming") ~ "Flow and Timing"
    ),
    
    DomainOrder = case_when(
      str_detect(DomainCode, "PatientSelection") ~ 1,
      str_detect(DomainCode, "IndexTest") ~ 2,
      str_detect(DomainCode, "ReferenceStandard") ~ 3,
      str_detect(DomainCode, "FlowTiming") ~ 4
    ),
    
    Assessment = factor(
      Assessment,
      levels = c("Risk of Bias", "Applicability Concerns")
    ),
    
    Judgement = factor(
      Judgement,
      levels = c("High", "Unclear", "Low")
    ),
    
    Symbol = case_when(
      Judgement == "High" ~ "−",
      Judgement == "Unclear" ~ "?",
      Judgement == "Low" ~ "+"
    ),
    
    Study = factor(Study, levels = rev(study_order))
  )

# Check data quality
stopifnot(nrow(quadas_long) == 20L * 7L, !anyNA(quadas_long$Judgement))
domain_levels <- c("Patient Selection", "Index Test", "Reference Standard", "Flow and Timing")
quadas_long$Domain <- factor(quadas_long$Domain, levels = domain_levels)
write.csv(quadas, file.path(output_dir, "QUADAS_AI_study_ratings.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(table(quadas_long$Judgement, useNA = "ifany"))

# ---------- 4. Colours ----------
risk_colours <- c(
  "High" = "#3C3B8B",
  "Unclear" = "#917BBD",
  "Low" = "#E5CBE1"
)

# =========================================================
# 5. Traffic-light plot
# =========================================================
p_traffic <- ggplot(
  quadas_long,
  aes(x = Domain, y = Study)
) +
  geom_point(
    aes(fill = Judgement),
    shape = 21,
    size = 7.2,
    colour = "black",
    stroke = 0.5,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = Symbol, colour = Judgement),
    fontface = "bold",
    family = "Arial",
    size = 14 / ggplot2::.pt,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  facet_grid(
    . ~ Assessment,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = risk_colours,
    name = NULL,
    drop = FALSE
  ) +
  scale_colour_manual(values = c(High = "white", Unclear = "black", Low = "black"), guide = "none") +
  labs(
    title = "QUADAS-AI Risk of Bias and Applicability Assessment",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(
      family = "Arial", face = "bold", size = 14
    ),
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      family = "Arial",
      size = 14,
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 14,
      colour = "black"
    ),
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "Arial",
      size = 14,
      face = "bold"
    ),
    panel.spacing.x = grid::unit(1.1, "cm"),
    plot.margin = margin(15, 20, 15, 20)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 5)
    )
  )

if (interactive()) print(p_traffic)

ggsave(
  filename = file.path(
    output_dir,
    "QUADAS_AI_traffic_light_circle_Arial14.png"
  ),
  plot = p_traffic,
  width = 14,
  height = max(10, length(study_order) * 0.55),
  dpi = 600,
  bg = "white"
)

# =========================================================
# 6. Summary bar plot
# =========================================================
summary_data <- quadas_long %>%
  filter(!is.na(Judgement)) %>%
  count(Assessment, Domain, DomainOrder, Judgement) %>%
  group_by(Assessment, Domain, DomainOrder) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

stopifnot(all(summary_data %>% group_by(Assessment, Domain) %>% summarise(Total = sum(n), .groups = "drop") %>% pull(Total) == 20L))
summary_data$Domain <- factor(summary_data$Domain, levels = rev(domain_levels))
write.csv(summary_data, file.path(output_dir, "QUADAS_AI_domain_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p_summary <- ggplot(
  summary_data,
  aes(x = Percentage, y = Domain, fill = Judgement)
) +
  geom_col(
    width = 0.7,
    colour = "white",
    linewidth = 0.45
  ) +
  facet_wrap(vars(Assessment), nrow = 1, scales = "free_y") +
  scale_fill_manual(
    values = risk_colours,
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "QUADAS-AI Risk of Bias and Applicability Summary",
    subtitle = "Each domain uses the same denominator",
    x = "Percentage of Studies",
    y = NULL
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(
      family = "Arial", face = "bold", size = 14
    ),
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      family = "Arial",
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 14,
      colour = "black"
    ),
    axis.title.x = element_text(
      family = "Arial",
      size = 14
    ),
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "Arial",
      size = 14,
      face = "bold"
    ),
    panel.spacing.x = grid::unit(1.1, "cm"),
    plot.margin = margin(15, 20, 15, 20)
  )

if (interactive()) print(p_summary)

ggsave(
  filename = file.path(
    output_dir,
    "QUADAS_AI_summary_bar_Arial14.png"
  ),
  plot = p_summary,
  width = 14,
  height = 7,
  dpi = 600,
  bg = "white"
)

cat("\nDone. Figures saved to:\n")
cat(file.path(output_dir, "QUADAS_AI_traffic_light_circle_Arial14.png"), "\n")
cat(file.path(output_dir, "QUADAS_AI_summary_bar_Arial14.png"), "\n")
# Vector outputs support assembly and verification of Supplementary Figure 1.
ggsave(file.path(output_dir, "QUADAS_AI_summary_bar_Arial14.pdf"), p_summary,
       width = 14, height = 7, device = cairo_pdf, bg = "white")
ggsave(file.path(output_dir, "QUADAS_AI_traffic_light_circle_Arial14.pdf"), p_traffic,
       width = 14, height = max(10, length(study_order) * 0.55), device = cairo_pdf, bg = "white")
