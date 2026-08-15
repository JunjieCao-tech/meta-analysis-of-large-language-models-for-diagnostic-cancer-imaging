# ============================================================
# Deeks Funnel Plot + Deeks Asymmetry Test
# 适用于：GBK / GB18030 编码、制表符分隔的 main.txt
# ============================================================

# -------------------- 1. 安装并加载包 --------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
}
library(ggplot2)

# -------------------- 2. 设置路径 --------------------
input_file <- "C:/Users/30399/Desktop/LLM/main.txt"
output_dir <- "C:/Users/30399/Desktop/LLM/7.Deeks"

if (!file.exists(input_file)) {
  stop("未找到输入文件：", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 固定紫色色板
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

# -------------------- 3. 安全读取 main.txt --------------------
# 原始文件是：制表符分隔（tab-separated），且可能为 GBK / GB18030 编码
read_txt_safely <- function(file_path) {
  
  encodings_to_try <- c("GB18030", "GBK", "UTF-8", "latin1")
  
  for (enc in encodings_to_try) {
    
    dat_try <- tryCatch({
      
      con <- file(file_path, open = "r", encoding = enc)
      on.exit(close(con), add = TRUE)
      
      raw_lines <- readLines(con, warn = FALSE)
      
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
    
    # 数据至少应包含 TP、FP、FN、TN 所在的前 15 列
    if (!is.null(dat_try) && nrow(dat_try) >= 3 && ncol(dat_try) >= 15) {
      message("成功读取文件；使用编码：", enc)
      return(dat_try)
    }
  }
  
  stop(
    "无法读取 main.txt。请确认文件没有被 Excel 占用，",
    "并确认其为制表符分隔的文本文件。"
  )
}

dat <- read_txt_safely(input_file)

# -------------------- 4. 按列位置提取 Deeks 所需变量 --------------------
# main.txt 中的固定列位置：
# 第1列：第一作者
# 第2列：年份
# 第11列：总样本N
# 第12列：TP
# 第13列：FP
# 第14列：FN
# 第15列：TN

deeks_dat <- data.frame(
  Author = as.character(dat[[1]]),
  Year   = as.character(dat[[2]]),
  TotalN = suppressWarnings(as.numeric(dat[[11]])),
  TP     = suppressWarnings(as.numeric(dat[[12]])),
  FP     = suppressWarnings(as.numeric(dat[[13]])),
  FN     = suppressWarnings(as.numeric(dat[[14]])),
  TN     = suppressWarnings(as.numeric(dat[[15]])),
  stringsAsFactors = FALSE
)

# 删除四格表数据不完整或不合理的研究
deeks_dat <- deeks_dat[
  complete.cases(deeks_dat[, c("TP", "FP", "FN", "TN")]) &
    deeks_dat$TP >= 0 &
    deeks_dat$FP >= 0 &
    deeks_dat$FN >= 0 &
    deeks_dat$TN >= 0,
  ,
  drop = FALSE
]

if (nrow(deeks_dat) < 3) {
  stop("有效研究数不足，无法进行 Deeks 不对称检验。")
}

# 研究标签
deeks_dat$Study <- paste0(trimws(deeks_dat$Author), " (", deeks_dat$Year, ")")

# -------------------- 5. 零单元格连续性校正 --------------------
# 只要四格表中任一格为 0，四格均加 0.5
deeks_dat$Zero_cell <- with(
  deeks_dat,
  TP == 0 | FP == 0 | FN == 0 | TN == 0
)

deeks_dat$TP_adj <- ifelse(deeks_dat$Zero_cell, deeks_dat$TP + 0.5, deeks_dat$TP)
deeks_dat$FP_adj <- ifelse(deeks_dat$Zero_cell, deeks_dat$FP + 0.5, deeks_dat$FP)
deeks_dat$FN_adj <- ifelse(deeks_dat$Zero_cell, deeks_dat$FN + 0.5, deeks_dat$FN)
deeks_dat$TN_adj <- ifelse(deeks_dat$Zero_cell, deeks_dat$TN + 0.5, deeks_dat$TN)

# -------------------- 6. 计算 Deeks 分析数据 --------------------
# DOR = (TP × TN) / (FP × FN)
deeks_dat$DOR <- with(
  deeks_dat,
  (TP_adj * TN_adj) / (FP_adj * FN_adj)
)

deeks_dat$lnDOR <- log(deeks_dat$DOR)

# ESS = 4 / [1/(TP + FN) + 1/(FP + TN)]
deeks_dat$ESS <- with(
  deeks_dat,
  4 / (1 / (TP_adj + FN_adj) + 1 / (FP_adj + TN_adj))
)

# Deeks 漏斗图横坐标
deeks_dat$Inv_sqrt_ESS <- 1 / sqrt(deeks_dat$ESS)

# 删除异常值
deeks_dat <- deeks_dat[
  is.finite(deeks_dat$lnDOR) &
    is.finite(deeks_dat$ESS) &
    deeks_dat$ESS > 0 &
    is.finite(deeks_dat$Inv_sqrt_ESS),
  ,
  drop = FALSE
]

k <- nrow(deeks_dat)

if (k < 3) {
  stop("可用于 Deeks 检验的研究不足 3 项。")
}

# -------------------- 7. 标准 Deeks 加权回归 --------------------
# ln(DOR) ~ 1/sqrt(ESS)
# 权重为 ESS
# Deeks 不对称性检验：检验斜率是否为 0

deeks_fit <- lm(
  lnDOR ~ Inv_sqrt_ESS,
  data = deeks_dat,
  weights = ESS
)

fit_sum <- summary(deeks_fit)
coef_tab <- fit_sum$coefficients

intercept_est <- coef_tab["(Intercept)", "Estimate"]
intercept_se  <- coef_tab["(Intercept)", "Std. Error"]
intercept_t   <- coef_tab["(Intercept)", "t value"]
intercept_p   <- coef_tab["(Intercept)", "Pr(>|t|)"]

slope_est <- coef_tab["Inv_sqrt_ESS", "Estimate"]
slope_se  <- coef_tab["Inv_sqrt_ESS", "Std. Error"]
slope_t   <- coef_tab["Inv_sqrt_ESS", "t value"]
slope_p   <- coef_tab["Inv_sqrt_ESS", "Pr(>|t|)"]

df_residual <- deeks_fit$df.residual
t_critical <- qt(0.975, df = df_residual)

slope_ci_lower <- slope_est - t_critical * slope_se
slope_ci_upper <- slope_est + t_critical * slope_se

# -------------------- 8. 输出分析数据 --------------------
analysis_data <- deeks_dat[, c(
  "Study", "Author", "Year", "TotalN",
  "TP", "FP", "FN", "TN",
  "Zero_cell",
  "TP_adj", "FP_adj", "FN_adj", "TN_adj",
  "DOR", "lnDOR", "ESS", "Inv_sqrt_ESS"
)]

write.csv(
  analysis_data,
  file = file.path(output_dir, "Deeks_Analysis_Data_Used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 9. 输出检验结果 --------------------
if (slope_p < 0.05) {
  interpretation <- paste0(
    "Deeks 不对称检验提示存在统计学显著的漏斗图不对称性，",
    "可能存在发表偏倚。"
  )
} else {
  interpretation <- paste0(
    "Deeks 不对称检验未显示统计学显著的漏斗图不对称性；",
    "但不能完全排除发表偏倚。"
  )
}

result_file <- file.path(output_dir, "Deeks_Asymmetry_Test_Results.txt")

sink(result_file)

cat("============================================================\n")
cat("Deeks Funnel Plot and Deeks Asymmetry Test\n")
cat("============================================================\n\n")

cat("Input file:", input_file, "\n")
cat("Number of included studies (k):", k, "\n\n")

cat("Method:\n")
cat("Weighted linear regression of ln(DOR) against 1/sqrt(ESS).\n")
cat("Weight = effective sample size (ESS).\n")
cat("If any cell in the 2 x 2 table equals zero, 0.5 is added to all cells.\n\n")

cat("Deeks asymmetry test (test of slope = 0):\n")
cat("Slope =", sprintf("%.4f", slope_est), "\n")
cat("SE =", sprintf("%.4f", slope_se), "\n")
cat(
  "95% CI = [",
  sprintf("%.4f", slope_ci_lower),
  ", ",
  sprintf("%.4f", slope_ci_upper),
  "]\n",
  sep = ""
)
cat("t =", sprintf("%.4f", slope_t), "\n")
cat("df =", df_residual, "\n")
cat(
  "P value for asymmetry =",
  format.pval(slope_p, digits = 4, eps = 0.0001),
  "\n\n"
)

cat("Intercept (for regression description only):\n")
cat("Intercept =", sprintf("%.4f", intercept_est), "\n")
cat("SE =", sprintf("%.4f", intercept_se), "\n")
cat("t =", sprintf("%.4f", intercept_t), "\n")
cat("P =", format.pval(intercept_p, digits = 4, eps = 0.0001), "\n\n")

cat("Conclusion:\n")
cat(interpretation, "\n")

sink()

# -------------------- 10. 回归线与 95% CI 数据 --------------------
x_seq <- seq(
  min(deeks_dat$Inv_sqrt_ESS),
  max(deeks_dat$Inv_sqrt_ESS),
  length.out = 300
)

pred_dat <- data.frame(Inv_sqrt_ESS = x_seq)

pred <- predict(
  deeks_fit,
  newdata = pred_dat,
  interval = "confidence",
  level = 0.95
)

pred_dat$fit <- pred[, "fit"]
pred_dat$lwr <- pred[, "lwr"]
pred_dat$upr <- pred[, "upr"]

# -------------------- 11. 绘制 Deeks 漏斗图 --------------------
p_text <- ifelse(
  slope_p < 0.0001,
  "P < 0.0001",
  paste0("P = ", sprintf("%.4f", slope_p))
)

annotation_text <- paste0(
  "Studies = ", k,
  "\nSlope = ", sprintf("%.3f", slope_est),
  "\n95% CI: ", sprintf("%.3f", slope_ci_lower),
  " to ", sprintf("%.3f", slope_ci_upper),
  "\nt = ", sprintf("%.3f", slope_t),
  ", df = ", df_residual,
  "\n", p_text
)

deeks_plot <- ggplot(
  deeks_dat,
  aes(x = Inv_sqrt_ESS, y = lnDOR)
) +
  geom_ribbon(
    data = pred_dat,
    aes(x = Inv_sqrt_ESS, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = pal[4],
    alpha = 0.38
  ) +
  geom_line(
    data = pred_dat,
    aes(x = Inv_sqrt_ESS, y = fit),
    inherit.aes = FALSE,
    colour = pal[1],
    linewidth = 1.1
  ) +
  geom_point(
    aes(size = ESS),
    shape = 21,
    fill = pal[3],
    colour = "white",
    stroke = 0.65,
    alpha = 0.92
  ) +
  annotate(
    "label",
    x = max(deeks_dat$Inv_sqrt_ESS) * 0.98,
    y = max(deeks_dat$lnDOR),
    label = annotation_text,
    hjust = 1,
    vjust = 1,
    family = "Arial",
    size = 4.2,
    colour = pal[1],
    fill = "white",
    label.size = 0.35,
    label.r = grid::unit(0.15, "lines")
  ) +
  scale_size_continuous(
    name = "Effective sample size",
    range = c(3.5, 10)
  ) +
  labs(
    title = "Deeks Funnel Plot",
    subtitle = "Weighted regression of ln(DOR) against 1/sqrt(ESS)",
    x = expression("1/" * sqrt(ESS)),
    y = expression(ln(DOR))
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
    legend.title = element_text(
      face = "bold",
      colour = pal[1]
    ),
    legend.position = c(0.10, 0.16),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      colour = pal[4],
      linewidth = 0.4
    ),
    plot.margin = margin(12, 18, 12, 12)
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

# -------------------- 12. 保存矢量 PDF --------------------
ggsave(
  filename = file.path(output_dir, "Deeks_Funnel_Plot.pdf"),
  plot = deeks_plot,
  width = 9,
  height = 7,
  units = "in",
  device = grDevices::cairo_pdf,
  family = "Arial"
)

# -------------------- 13. 提示完成 --------------------
cat("\n============================================================\n")
cat("Deeks 分析完成。\n")
cat("输出路径：", output_dir, "\n")
cat(
  "P value for asymmetry：",
  format.pval(slope_p, digits = 4, eps = 0.0001),
  "\n"
)
cat("============================================================\n")