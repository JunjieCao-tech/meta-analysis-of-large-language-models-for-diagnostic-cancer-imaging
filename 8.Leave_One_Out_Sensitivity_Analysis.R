# ============================================================
# 逐一剔除研究敏感性分析（Leave-one-out）
# 指标：Sensitivity、Specificity、Accuracy、NPV、PPV
# 不计算 AUC
# ============================================================

# -------------------- 1. 安装并加载 R 包 --------------------
if (!requireNamespace("metafor", quietly = TRUE)) {
  install.packages("metafor", dependencies = TRUE)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
}

library(metafor)
library(ggplot2)

# -------------------- 2. 设置路径 --------------------
input_file <- "C:/Users/30399/Desktop/LLM/main.txt"
output_dir <- "C:/Users/30399/Desktop/LLM/8.Leave_One_Out_Sensitivity_Analysis"

if (!file.exists(input_file)) {
  stop("未找到原始文件：", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 通用字体，避免 Arial 字体报错
font_family <- "sans"

# 紫色色板
pal <- c("#3C3B8B", "#6656A6", "#917BBD", "#B19CCB", "#E5CBE1")

# 预设“影响性研究”判定阈值
effect_change_threshold <- 0.05
i2_drop_threshold <- 25

# -------------------- 3. 安全读取 main.txt --------------------
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
    
    if (!is.null(dat_try) && nrow(dat_try) >= 4 && ncol(dat_try) >= 15) {
      message("成功读取 main.txt；使用编码：", enc)
      return(dat_try)
    }
  }
  
  stop("无法读取 main.txt。请确认文件为制表符分隔文本，且未被 Excel 占用。")
}

dat <- read_txt_safely(input_file)

# -------------------- 4. 提取四格表数据 --------------------
meta_dat <- data.frame(
  Author = as.character(dat[[1]]),
  Year = as.character(dat[[2]]),
  TotalN = suppressWarnings(as.numeric(dat[[11]])),
  TP = suppressWarnings(as.numeric(dat[[12]])),
  FP = suppressWarnings(as.numeric(dat[[13]])),
  FN = suppressWarnings(as.numeric(dat[[14]])),
  TN = suppressWarnings(as.numeric(dat[[15]])),
  stringsAsFactors = FALSE
)

meta_dat <- meta_dat[
  complete.cases(meta_dat[, c("TP", "FP", "FN", "TN")]) &
    meta_dat$TP >= 0 &
    meta_dat$FP >= 0 &
    meta_dat$FN >= 0 &
    meta_dat$TN >= 0 &
    (meta_dat$TP + meta_dat$FN) > 0 &
    (meta_dat$FP + meta_dat$TN) > 0 &
    (meta_dat$TP + meta_dat$FP) > 0 &
    (meta_dat$TN + meta_dat$FN) > 0,
  ,
  drop = FALSE
]

if (nrow(meta_dat) < 4) {
  stop("有效研究数少于 4 项，无法进行逐一剔除敏感性分析。")
}

meta_dat$Study <- paste0(
  trimws(meta_dat$Author),
  " (",
  trimws(meta_dat$Year),
  ")"
)

meta_dat$Study <- make.unique(meta_dat$Study, sep = "_")

# -------------------- 5. 连续性校正 --------------------
meta_dat$Zero_cell <- with(
  meta_dat,
  TP == 0 | FP == 0 | FN == 0 | TN == 0
)

meta_dat$TP_adj <- ifelse(meta_dat$Zero_cell, meta_dat$TP + 0.5, meta_dat$TP)
meta_dat$FP_adj <- ifelse(meta_dat$Zero_cell, meta_dat$FP + 0.5, meta_dat$FP)
meta_dat$FN_adj <- ifelse(meta_dat$Zero_cell, meta_dat$FN + 0.5, meta_dat$FN)
meta_dat$TN_adj <- ifelse(meta_dat$Zero_cell, meta_dat$TN + 0.5, meta_dat$TN)

# -------------------- 6. 计算各研究诊断指标 --------------------
meta_dat$Sensitivity <- with(meta_dat, TP_adj / (TP_adj + FN_adj))
meta_dat$Specificity <- with(meta_dat, TN_adj / (TN_adj + FP_adj))
meta_dat$Accuracy <- with(
  meta_dat,
  (TP_adj + TN_adj) / (TP_adj + FP_adj + FN_adj + TN_adj)
)
meta_dat$PPV <- with(meta_dat, TP_adj / (TP_adj + FP_adj))
meta_dat$NPV <- with(meta_dat, TN_adj / (TN_adj + FN_adj))

# -------------------- 7. 随机效应比例 Meta 分析函数 --------------------
run_proportion_meta <- function(events, totals) {
  
  proportion <- events / totals
  
  yi <- log(proportion / (1 - proportion))
  vi <- 1 / events + 1 / (totals - events)
  
  fit <- rma(
    yi = yi,
    vi = vi,
    method = "REML"
  )
  
  pooled_logit <- as.numeric(fit$b[1, 1])
  pooled_se <- as.numeric(fit$se[1])
  
  ci_lower_logit <- pooled_logit - 1.96 * pooled_se
  ci_upper_logit <- pooled_logit + 1.96 * pooled_se
  
  list(
    pooled = plogis(pooled_logit),
    lower = plogis(ci_lower_logit),
    upper = plogis(ci_upper_logit),
    I2 = as.numeric(fit$I2),
    tau2 = as.numeric(fit$tau2),
    Q = as.numeric(fit$QE),
    Q_p = as.numeric(fit$QEp)
  )
}

# -------------------- 8. 完整分析函数 --------------------
run_full_meta_analysis <- function(dat_sub) {
  
  sen_fit <- run_proportion_meta(
    events = dat_sub$TP_adj,
    totals = dat_sub$TP_adj + dat_sub$FN_adj
  )
  
  spe_fit <- run_proportion_meta(
    events = dat_sub$TN_adj,
    totals = dat_sub$TN_adj + dat_sub$FP_adj
  )
  
  acc_fit <- run_proportion_meta(
    events = dat_sub$TP_adj + dat_sub$TN_adj,
    totals = dat_sub$TP_adj + dat_sub$FP_adj +
      dat_sub$FN_adj + dat_sub$TN_adj
  )
  
  ppv_fit <- run_proportion_meta(
    events = dat_sub$TP_adj,
    totals = dat_sub$TP_adj + dat_sub$FP_adj
  )
  
  npv_fit <- run_proportion_meta(
    events = dat_sub$TN_adj,
    totals = dat_sub$TN_adj + dat_sub$FN_adj
  )
  
  data.frame(
    k = nrow(dat_sub),
    
    Sen = sen_fit$pooled,
    Sen_Lower = sen_fit$lower,
    Sen_Upper = sen_fit$upper,
    Sen_I2 = sen_fit$I2,
    Sen_tau2 = sen_fit$tau2,
    Sen_Q_p = sen_fit$Q_p,
    
    Spe = spe_fit$pooled,
    Spe_Lower = spe_fit$lower,
    Spe_Upper = spe_fit$upper,
    Spe_I2 = spe_fit$I2,
    Spe_tau2 = spe_fit$tau2,
    Spe_Q_p = spe_fit$Q_p,
    
    Accuracy = acc_fit$pooled,
    Accuracy_Lower = acc_fit$lower,
    Accuracy_Upper = acc_fit$upper,
    Accuracy_I2 = acc_fit$I2,
    Accuracy_tau2 = acc_fit$tau2,
    Accuracy_Q_p = acc_fit$Q_p,
    
    PPV = ppv_fit$pooled,
    PPV_Lower = ppv_fit$lower,
    PPV_Upper = ppv_fit$upper,
    PPV_I2 = ppv_fit$I2,
    PPV_tau2 = ppv_fit$tau2,
    PPV_Q_p = ppv_fit$Q_p,
    
    NPV = npv_fit$pooled,
    NPV_Lower = npv_fit$lower,
    NPV_Upper = npv_fit$upper,
    NPV_I2 = npv_fit$I2,
    NPV_tau2 = npv_fit$tau2,
    NPV_Q_p = npv_fit$Q_p,
    
    stringsAsFactors = FALSE
  )
}

# -------------------- 9. 主分析 --------------------
cat("正在进行主分析...\n")

overall_result <- run_full_meta_analysis(meta_dat)
overall_result$Analysis <- "Overall"
overall_result$Excluded_Study <- "None"

overall_result <- overall_result[, c(
  "Analysis", "Excluded_Study", "k",
  
  "Sen", "Sen_Lower", "Sen_Upper", "Sen_I2", "Sen_tau2", "Sen_Q_p",
  "Spe", "Spe_Lower", "Spe_Upper", "Spe_I2", "Spe_tau2", "Spe_Q_p",
  "Accuracy", "Accuracy_Lower", "Accuracy_Upper",
  "Accuracy_I2", "Accuracy_tau2", "Accuracy_Q_p",
  "PPV", "PPV_Lower", "PPV_Upper", "PPV_I2", "PPV_tau2", "PPV_Q_p",
  "NPV", "NPV_Lower", "NPV_Upper", "NPV_I2", "NPV_tau2", "NPV_Q_p"
)]

# -------------------- 10. 逐一剔除单篇研究 --------------------
cat("正在逐一剔除研究...\n")

loo_list <- vector("list", nrow(meta_dat))

for (i in seq_len(nrow(meta_dat))) {
  
  excluded_study <- meta_dat$Study[i]
  dat_loo <- meta_dat[-i, , drop = FALSE]
  
  cat("剔除：", excluded_study, "\n")
  
  loo_result <- tryCatch(
    run_full_meta_analysis(dat_loo),
    error = function(e) {
      data.frame(
        k = nrow(dat_loo),
        
        Sen = NA_real_, Sen_Lower = NA_real_, Sen_Upper = NA_real_,
        Sen_I2 = NA_real_, Sen_tau2 = NA_real_, Sen_Q_p = NA_real_,
        
        Spe = NA_real_, Spe_Lower = NA_real_, Spe_Upper = NA_real_,
        Spe_I2 = NA_real_, Spe_tau2 = NA_real_, Spe_Q_p = NA_real_,
        
        Accuracy = NA_real_, Accuracy_Lower = NA_real_,
        Accuracy_Upper = NA_real_, Accuracy_I2 = NA_real_,
        Accuracy_tau2 = NA_real_, Accuracy_Q_p = NA_real_,
        
        PPV = NA_real_, PPV_Lower = NA_real_, PPV_Upper = NA_real_,
        PPV_I2 = NA_real_, PPV_tau2 = NA_real_, PPV_Q_p = NA_real_,
        
        NPV = NA_real_, NPV_Lower = NA_real_, NPV_Upper = NA_real_,
        NPV_I2 = NA_real_, NPV_tau2 = NA_real_, NPV_Q_p = NA_real_
      )
    }
  )
  
  loo_result$Analysis <- "Leave-one-out"
  loo_result$Excluded_Study <- excluded_study
  loo_result <- loo_result[, names(overall_result)]
  
  loo_list[[i]] <- loo_result
}

loo_results <- do.call(rbind, loo_list)

# -------------------- 11. 计算五项指标变化 --------------------
indicators <- c("Sen", "Spe", "Accuracy", "PPV", "NPV")

for (indicator in indicators) {
  
  loo_results[[paste0(indicator, "_Change")]] <-
    loo_results[[indicator]] - overall_result[[indicator]]
  
  loo_results[[paste0(indicator, "_I2_Drop")]] <-
    overall_result[[paste0(indicator, "_I2")]] -
    loo_results[[paste0(indicator, "_I2")]]
  
  loo_results[[paste0(indicator, "_Material_Change")]] <-
    abs(loo_results[[paste0(indicator, "_Change")]]) >=
    effect_change_threshold
  
  loo_results[[paste0(indicator, "_I2_Material_Drop")]] <-
    loo_results[[paste0(indicator, "_I2_Drop")]] >=
    i2_drop_threshold
}

# -------------------- 12. 标记潜在影响性研究 --------------------
loo_results$Potential_Influential_Study <- FALSE
loo_results$Flag_Reason <- ""

for (i in seq_len(nrow(loo_results))) {
  
  reasons <- c()
  
  for (indicator in indicators) {
    
    effect_flag <- loo_results[[paste0(indicator, "_Material_Change")]][i]
    i2_flag <- loo_results[[paste0(indicator, "_I2_Material_Drop")]][i]
    
    if (isTRUE(effect_flag)) {
      reasons <- c(reasons, paste0(indicator, " change >= 0.05"))
    }
    
    if (isTRUE(i2_flag)) {
      reasons <- c(reasons, paste0(indicator, " I2 drop >= 25%"))
    }
  }
  
  if (length(reasons) > 0) {
    loo_results$Potential_Influential_Study[i] <- TRUE
    loo_results$Flag_Reason[i] <- paste(reasons, collapse = "; ")
  } else {
    loo_results$Flag_Reason[i] <- "No material influence detected"
  }
}

# -------------------- 13. 保存结果 --------------------
write.csv(
  meta_dat,
  file = file.path(output_dir, "Sensitivity_Analysis_Data_Used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  overall_result,
  file = file.path(output_dir, "Overall_Meta_Analysis_Results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  loo_results,
  file = file.path(output_dir, "Leave_One_Out_Sensitivity_Analysis_Results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

influential_studies <- loo_results[
  loo_results$Potential_Influential_Study,
  ,
  drop = FALSE
]

write.csv(
  influential_studies,
  file = file.path(output_dir, "Potentially_Influential_Studies.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------- 14. 输出文字报告 --------------------
report_file <- file.path(
  output_dir,
  "Leave_One_Out_Sensitivity_Analysis_Report.txt"
)

sink(report_file)

cat("============================================================\n")
cat("Leave-one-out Sensitivity Analysis\n")
cat("Indicators: Sensitivity, Specificity, Accuracy, PPV, NPV\n")
cat("============================================================\n\n")

cat("Number of studies:", nrow(meta_dat), "\n\n")

cat("Overall pooled results:\n")
cat("Sensitivity =", sprintf("%.4f", overall_result$Sen), "\n")
cat("Specificity =", sprintf("%.4f", overall_result$Spe), "\n")
cat("Accuracy =", sprintf("%.4f", overall_result$Accuracy), "\n")
cat("PPV =", sprintf("%.4f", overall_result$PPV), "\n")
cat("NPV =", sprintf("%.4f", overall_result$NPV), "\n\n")

cat("Overall heterogeneity (I2):\n")
cat("Sensitivity I2 =", sprintf("%.2f", overall_result$Sen_I2), "%\n")
cat("Specificity I2 =", sprintf("%.2f", overall_result$Spe_I2), "%\n")
cat("Accuracy I2 =", sprintf("%.2f", overall_result$Accuracy_I2), "%\n")
cat("PPV I2 =", sprintf("%.2f", overall_result$PPV_I2), "%\n")
cat("NPV I2 =", sprintf("%.2f", overall_result$NPV_I2), "%\n\n")

cat("Flag criteria:\n")
cat("- Absolute pooled estimate change >= 0.05\n")
cat("- I2 reduction >= 25%\n\n")

if (nrow(influential_studies) == 0) {
  cat("Conclusion:\n")
  cat("No study met the predefined influence criteria.\n")
  cat("The pooled diagnostic performance estimates remained stable after sequential exclusion.\n")
} else {
  cat("Potentially influential studies:\n\n")
  
  for (i in seq_len(nrow(influential_studies))) {
    cat("Excluded study:", influential_studies$Excluded_Study[i], "\n")
    cat("Reason:", influential_studies$Flag_Reason[i], "\n\n")
  }
  
  cat("Conclusion:\n")
  cat("The listed studies may contribute substantially to heterogeneity or pooled-effect variation.\n")
  cat("They should be discussed according to population, threshold, reference standard,\n")
  cat("study design, and risk-of-bias characteristics; they should not be excluded solely\n")
  cat("because they influence the pooled estimate.\n")
}

sink()

# -------------------- 15. 绘图数据：五项合并效应值 --------------------
plot_dat <- loo_results

plot_dat$Excluded_Study <- factor(
  plot_dat$Excluded_Study,
  levels = rev(plot_dat$Excluded_Study)
)

effect_plot_dat <- rbind(
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Sensitivity",
    Estimate = plot_dat$Sen,
    Overall = overall_result$Sen
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Specificity",
    Estimate = plot_dat$Spe,
    Overall = overall_result$Spe
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Accuracy",
    Estimate = plot_dat$Accuracy,
    Overall = overall_result$Accuracy
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "PPV",
    Estimate = plot_dat$PPV,
    Overall = overall_result$PPV
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "NPV",
    Estimate = plot_dat$NPV,
    Overall = overall_result$NPV
  )
)

effect_plot <- ggplot(
  effect_plot_dat,
  aes(x = Estimate, y = Excluded_Study)
) +
  geom_vline(
    aes(xintercept = Overall),
    colour = pal[2],
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    fill = pal[3],
    colour = "white",
    stroke = 0.6
  ) +
  facet_wrap(~ Indicator, scales = "free_x", nrow = 1) +
  
  # 仅修改横坐标：每个分面约显示 3 个刻度，保留 3 位小数
  scale_x_continuous(
    n.breaks = 3,
    labels = scales::label_number(accuracy = 0.001)
  ) +
  
  labs(
    title = "Leave-one-out Sensitivity Analysis",
    subtitle = "Dashed line indicates the overall pooled estimate",
    x = "Pooled estimate after excluding one study",
    y = "Excluded study"
  ) +
  theme_classic(base_size = 13, base_family = font_family) +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1],
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      colour = pal[2],
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      colour = pal[1]
    ),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 10),
    strip.background = element_rect(
      fill = pal[5],
      colour = pal[4]
    ),
    strip.text = element_text(
      face = "bold",
      colour = pal[1]
    )
  )

# -------------------- 16. 绘图数据：五项 I² --------------------
i2_plot_dat <- rbind(
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Sensitivity I2",
    I2 = plot_dat$Sen_I2,
    Overall_I2 = overall_result$Sen_I2
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Specificity I2",
    I2 = plot_dat$Spe_I2,
    Overall_I2 = overall_result$Spe_I2
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "Accuracy I2",
    I2 = plot_dat$Accuracy_I2,
    Overall_I2 = overall_result$Accuracy_I2
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "PPV I2",
    I2 = plot_dat$PPV_I2,
    Overall_I2 = overall_result$PPV_I2
  ),
  data.frame(
    Excluded_Study = plot_dat$Excluded_Study,
    Indicator = "NPV I2",
    I2 = plot_dat$NPV_I2,
    Overall_I2 = overall_result$NPV_I2
  )
)

i2_plot <- ggplot(
  i2_plot_dat,
  aes(x = I2, y = Excluded_Study)
) +
  geom_vline(
    aes(xintercept = Overall_I2),
    colour = pal[2],
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    fill = pal[3],
    colour = "white",
    stroke = 0.6
  ) +
  facet_wrap(~ Indicator, scales = "free_x", nrow = 1) +
  
  # 同样避免 I² 横坐标刻度文字重叠
  scale_x_continuous(
    n.breaks = 3,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  
  labs(
    title = "Heterogeneity After Excluding Each Study",
    subtitle = "Dashed line indicates the overall I2",
    x = "I2 (%) after excluding one study",
    y = "Excluded study"
  ) +
  theme_classic(base_size = 13, base_family = font_family) +
  theme(
    plot.title = element_text(
      face = "bold",
      colour = pal[1],
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      colour = pal[2],
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      colour = pal[1]
    ),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 10),
    strip.background = element_rect(
      fill = pal[5],
      colour = pal[4]
    ),
    strip.text = element_text(
      face = "bold",
      colour = pal[1]
    )
  )

# -------------------- 17. 保存两页矢量 PDF --------------------
pdf_file <- file.path(
  output_dir,
  "Leave_One_Out_Sensitivity_Analysis_Plots.pdf"
)

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 16,
  height = 9
)

print(effect_plot)
print(i2_plot)

dev.off()

# -------------------- 18. 完成提示 --------------------
cat("\n============================================================\n")
cat("逐一剔除研究敏感性分析完成。\n")
cat("输出目录：", output_dir, "\n")
cat("纳入研究数：", nrow(meta_dat), "\n")
cat("潜在影响性研究数：", nrow(influential_studies), "\n")
cat("============================================================\n")