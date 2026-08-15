# ============================================================
# 阈值效应检验：Spearman 相关系数
# logit(Sensitivity) 与 logit(1 - Specificity) 的相关性
# ============================================================

# -------------------- 1. 安装并加载 R 包 --------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
}

library(ggplot2)

# -------------------- 2. 设置输入与输出路径 --------------------
input_file <- "C:/Users/30399/Desktop/LLM/main.txt"
output_dir <- "C:/Users/30399/Desktop/LLM/6.Threshold"

if (!file.exists(input_file)) {
  stop("未找到原始数据文件：", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 固定紫色色板
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

# -------------------- 3. 安全读取 main.txt --------------------
# 自动尝试常用编码；main.txt 为制表符分隔文件
read_txt_safely <- function(file_path) {
  
  encodings_to_try <- c("GB18030", "GBK", "UTF-8", "latin1")
  
  for (enc in encodings_to_try) {
    
    dat_try <- tryCatch({
      
      con <- file(file_path, open = "r", encoding = enc)
      raw_lines <- readLines(con, warn = FALSE)
      close(con)
      
      read.table(
        text = raw_lines,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        fill = TRUE,
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = c("", "NA", "N/A")
      )
      
    }, error = function(e) NULL)
    
    if (!is.null(dat_try) && nrow(dat_try) >= 3 && ncol(dat_try) >= 15) {
      message("成功读取 main.txt；使用编码：", enc)
      return(dat_try)
    }
  }
  
  stop("无法读取 main.txt。请确认文件未被 Excel 占用，且为制表符分隔文本。")
}

dat <- read_txt_safely(input_file)

# -------------------- 4. 提取四格表数据 --------------------
# 第1列：第一作者
# 第2列：年份
# 第11列：总样本N
# 第12列：TP
# 第13列：FP
# 第14列：FN
# 第15列：TN

threshold_dat <- data.frame(
  Author = as.character(dat[[1]]),
  Year   = as.character(dat[[2]]),
  TotalN = suppressWarnings(as.numeric(dat[[11]])),
  TP     = suppressWarnings(as.numeric(dat[[12]])),
  FP     = suppressWarnings(as.numeric(dat[[13]])),
  FN     = suppressWarnings(as.numeric(dat[[14]])),
  TN     = suppressWarnings(as.numeric(dat[[15]])),
  stringsAsFactors = FALSE
)

# 删除四格表缺失或不合理研究
threshold_dat <- threshold_dat[
  complete.cases(threshold_dat[, c("TP", "FP", "FN", "TN")]) &
    threshold_dat$TP >= 0 &
    threshold_dat$FP >= 0 &
    threshold_dat$FN >= 0 &
    threshold_dat$TN >= 0 &
    (threshold_dat$TP + threshold_dat$FN) > 0 &
    (threshold_dat$FP + threshold_dat$TN) > 0,
  ,
  drop = FALSE
]

if (nrow(threshold_dat) < 3) {
  stop("有效研究数量少于 3 项，不能进行阈值效应检验。")
}

threshold_dat$Study <- paste0(
  trimws(threshold_dat$Author),
  " (",
  trimws(threshold_dat$Year),
  ")"
)

# -------------------- 5. 连续性校正 --------------------
# 若任何四格表单元格为 0，全部单元格加 0.5
# 这样可以避免 Sensitivity 或 Specificity 等于 0 或 1，
# 从而避免 logit 值为无穷大。

threshold_dat$Zero_cell <- with(
  threshold_dat,
  TP == 0 | FP == 0 | FN == 0 | TN == 0
)

threshold_dat$TP_adj <- ifelse(
  threshold_dat$Zero_cell,
  threshold_dat$TP + 0.5,
  threshold_dat$TP
)

threshold_dat$FP_adj <- ifelse(
  threshold_dat$Zero_cell,
  threshold_dat$FP + 0.5,
  threshold_dat$FP
)

threshold_dat$FN_adj <- ifelse(
  threshold_dat$Zero_cell,
  threshold_dat$FN + 0.5,
  threshold_dat$FN
)

threshold_dat$TN_adj <- ifelse(
  threshold_dat$Zero_cell,
  threshold_dat$TN + 0.5,
  threshold_dat$TN
)

# -------------------- 6. 计算 Sensitivity、Specificity 及 logit 值 --------------------
threshold_dat$Sensitivity <- with(
  threshold_dat,
  TP_adj / (TP_adj + FN_adj)
)

threshold_dat$Specificity <- with(
  threshold_dat,
  TN_adj / (TN_adj + FP_adj)
)

threshold_dat$False_Positive_Rate <- 1 - threshold_dat$Specificity

# logit(p) = log[p / (1 - p)]
threshold_dat$logit_Sensitivity <- with(
  threshold_dat,
  log(Sensitivity / (1 - Sensitivity))
)

threshold_dat$logit_FPR <- with(
  threshold_dat,
  log(False_Positive_Rate / (1 - False_Positive_Rate))
)

# 删除不可计算值
threshold_dat <- threshold_dat[
  is.finite(threshold_dat$Sensitivity) &
    is.finite(threshold_dat$Specificity) &
    is.finite(threshold_dat$False_Positive_Rate) &
    is.finite(threshold_dat$logit_Sensitivity) &
    is.finite(threshold_dat$logit_FPR),
  ,
  drop = FALSE
]

k <- nrow(threshold_dat)

if (k < 3) {
  stop("可用于阈值效应检验的研究数量少于 3 项。")
}

# -------------------- 7. Spearman 阈值效应检验 --------------------
# H0：logit(Sensitivity) 与 logit(1 - Specificity) 无相关性
# 若 Spearman rho 为正且 P < 0.05，通常提示存在阈值效应。

threshold_test <- cor.test(
  x = threshold_dat$logit_Sensitivity,
  y = threshold_dat$logit_FPR,
  method = "spearman",
  exact = FALSE,
  alternative = "two.sided"
)

rho <- unname(threshold_test$estimate)
p_value <- threshold_test$p.value

p_text <- ifelse(
  p_value < 0.0001,
  "P < 0.0001",
  paste0("P = ", sprintf("%.4f", p_value))
)

interpretation <- if (rho > 0 && p_value < 0.05) {
  "Spearman 相关系数为正且 P < 0.05，提示存在统计学显著的阈值效应。"
} else {
  "Spearman 相关性未达到统计学显著性，未发现明确的阈值效应证据。"
}

# -------------------- 8. 保存阈值效应分析数据 --------------------
output_data <- threshold_dat[, c(
  "Study", "Author", "Year", "TotalN",
  "TP", "FP", "FN", "TN",
  "Zero_cell",
  "TP_adj", "FP_adj", "FN_adj", "TN_adj",
  "Sensitivity", "Specificity", "False_Positive_Rate",
  "logit_Sensitivity", "logit_FPR"
)]

write.csv(
  output_data,
  file = file.path(output_dir, "Threshold_Effect_Analysis_Data_Used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 9. 保存阈值效应检验结果 --------------------
result_file <- file.path(
  output_dir,
  "Threshold_Effect_Spearman_Test_Results.txt"
)

sink(result_file)

cat("============================================================\n")
cat("Threshold Effect Test: Spearman Correlation Analysis\n")
cat("============================================================\n\n")

cat("Input file: ", input_file, "\n", sep = "")
cat("Number of included studies (k): ", k, "\n\n", sep = "")

cat("Method:\n")
cat("Spearman correlation between logit(Sensitivity) and ")
cat("logit(1 - Specificity).\n")
cat("If any cell in the 2 x 2 table equals zero, 0.5 is added ")
cat("to all four cells.\n\n")

cat("Results:\n")
cat("Spearman correlation coefficient (rho) = ", sprintf("%.4f", rho), "\n", sep = "")
cat("P value = ", format.pval(p_value, digits = 4, eps = 0.0001), "\n\n", sep = "")

cat("Interpretation:\n")
cat(interpretation, "\n\n")

cat("Interpretation rule:\n")
cat("A positive Spearman correlation with P < 0.05 is generally ")
cat("considered evidence of a threshold effect.\n")

sink()

# -------------------- 10. 绘制阈值效应散点图 --------------------
# 仅用于趋势展示；统计推断以 Spearman 检验结果为准。

linear_fit <- lm(
  logit_Sensitivity ~ logit_FPR,
  data = threshold_dat
)

annotation_text <- paste0(
  "Studies = ", k,
  "\nSpearman rho = ", sprintf("%.3f", rho),
  "\n", p_text
)

threshold_plot <- ggplot(
  threshold_dat,
  aes(x = logit_FPR, y = logit_Sensitivity)
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = pal[1],
    fill = pal[4],
    linewidth = 1.05,
    alpha = 0.40
  ) +
  geom_point(
    aes(size = TotalN),
    shape = 21,
    fill = pal[3],
    colour = "white",
    stroke = 0.65,
    alpha = 0.92
  ) +
  annotate(
    "label",
    x = max(threshold_dat$logit_FPR),
    y = max(threshold_dat$logit_Sensitivity),
    label = annotation_text,
    hjust = 1.05,
    vjust = 1.15,
    family = "Arial",
    size = 4.2,
    colour = pal[1],
    fill = "white",
    label.size = 0.35,
    label.r = grid::unit(0.15, "lines")
  ) +
  scale_size_continuous(
    name = "Sample size",
    range = c(3.5, 9)
  ) +
  labs(
    title = "Threshold Effect Assessment",
    subtitle = "Spearman correlation between logit(Sensitivity) and logit(1 − Specificity)",
    x = expression("logit(1 " - " Specificity)"),
    y = expression("logit(Sensitivity)")
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1],
      size = 16,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      colour = pal[2],
      size = 12,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title = element_text(
      face = "bold",
      colour = pal[1]
    ),
    axis.text = element_text(colour = "black"),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      colour = pal[1],
      size = 13
    ),
    legend.text = element_text(
      colour = "black",
      size = 12
    ),
    legend.background = element_rect(
      fill = "white",
      colour = pal[4],
      linewidth = 0.5
    ),
    legend.key = element_rect(
      fill = "white",
      colour = "white"
    ),
    legend.spacing.y = grid::unit(0.25, "cm"),
    legend.margin = margin(8, 10, 8, 10),
    plot.margin = margin(12, 22, 12, 12)
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        shape = 21,
        fill = pal[3],
        colour = "white",
        stroke = 0.65,
        alpha = 0.92
      )
    )
  )

# -------------------- 11. 保存矢量 PDF 图 --------------------
ggsave(
  filename = file.path(output_dir, "Threshold_Effect_Scatter_Plot.pdf"),
  plot = threshold_plot,
  width = 11,
  height = 7,
  units = "in",
  device = grDevices::cairo_pdf,
  family = "Arial"
)

# -------------------- 12. 输出完成提示 --------------------
cat("\n============================================================\n")
cat("阈值效应检验完成。\n")
cat("输出路径：", output_dir, "\n")
cat("Spearman rho：", sprintf("%.4f", rho), "\n")
cat("P value：", format.pval(p_value, digits = 4, eps = 0.0001), "\n")
cat("============================================================\n")