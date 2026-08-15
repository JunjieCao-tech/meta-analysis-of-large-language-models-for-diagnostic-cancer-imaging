# =========================================================
# QUADAS-AI publication-style plots
# Font: Arial, 14 pt
# =========================================================

# ---------- 0. Install and load packages ----------
pkgs <- c("dplyr", "tidyr", "ggplot2", "stringr")

new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ---------- 1. File path ----------
work_dir <- "C:/Users/30399/Desktop/LLM"
data_file <- file.path(work_dir, "main.txt")

# If your Desktop is under OneDrive, change the path to:
# work_dir <- "C:/Users/30399/OneDrive/Desktop/LLM"
# data_file <- file.path(work_dir, "main.txt")

# ---------- 2. Read data ----------
dat <- read.delim(
  data_file,
  header = TRUE,
  sep = "\t",
  fileEncoding = "GB18030",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Extract: Author + Year + last 7 QUADAS-AI columns
quadas <- dat[, c(1, 2, (ncol(dat) - 6):ncol(dat))]

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
print(table(quadas_long$Judgement, useNA = "ifany"))

# ---------- 4. Colours ----------
risk_colours <- c(
  "High" = "#D40000",
  "Unclear" = "#E5E500",
  "Low" = "#00B500"
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
    aes(label = Symbol),
    fontface = "bold",
    size = 4.5,
    colour = "black",
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
  labs(
    title = "QUADAS-AI Risk of Bias Assessment",
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

print(p_traffic)

ggsave(
  filename = file.path(
    work_dir,
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

p_summary <- ggplot(
  summary_data,
  aes(x = Percentage, y = Domain, fill = Judgement)
) +
  geom_col(
    width = 0.7,
    colour = "white",
    linewidth = 0.45
  ) +
  facet_grid(
    . ~ Assessment,
    scales = "free_y",
    space = "free_y"
  ) +
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
    title = "QUADAS-AI Risk of Bias Summary",
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

print(p_summary)

ggsave(
  filename = file.path(
    work_dir,
    "QUADAS_AI_summary_bar_Arial14.png"
  ),
  plot = p_summary,
  width = 14,
  height = 7,
  dpi = 600,
  bg = "white"
)

cat("\nDone. Figures saved to:\n")
cat(file.path(work_dir, "QUADAS_AI_traffic_light_circle_Arial14.png"), "\n")
cat(file.path(work_dir, "QUADAS_AI_summary_bar_Arial14.png"), "\n")