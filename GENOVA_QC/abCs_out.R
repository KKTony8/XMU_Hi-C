# 10:19 2025/7/2 tested by wky
# 整个workflow的思路就是利用了A/B区室的划分就是利用主成分分析PCA来定义，用PC1的分数来划分区室
# 流程： 对每个样本的染色体按分辨率区间bin来划分，计算得到相关的矩阵；再通过距离标准化得到O/E矩阵；对矩阵进行PCA分析提取PC1分数作为区室得分
# 例如：  
    #对于样本1，相关矩阵的行（或列）作为变量，进行PCA。假设得到的PC1为：
        # bin1: -0.5, bin2: -0.4, bin3: 0.0, bin4: 0.4, bin5: 0.5
    # 对于样本2，假设PC1为：
        # bin1: -0.6, bin2: -0.5, bin3: 0.1, bin4: 0.4, bin5: 0.6
    # 得到的结果Cs(compartment score)_out是这样格式：
    #| chrom | start   | end     | bin | SRR5050 | XM12  |
    #|-------|---------|---------|-----|---------|-------|
    #| chr1  | 0       | 500000  | 1   | 0.5     | 0.3   |
    #| chr1  | 500000  | 1000000 | 2   | -0.2    | -0.1  |
    #| chr1  | 1000000 | 1500000 | 3   | -0.7    | -0.6  |

# !/usr/bin/env Rscript

# plot_compartment_chr1.R
# 对 SRR5005050 和 XM-QGY-MME-13_L7_1 两个样本在 chr1 上做 A/B compartment 分数绘图
library(GENOVA)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 参数设置
cool_dir    <- "/data5/Wuky/2XMU_T7-HiC/MME12scc/cool_500kb"
resolution  <- 500000
chr_to_plot <- "chr1"

# 更新文件路径和样本名称
file1 <- file.path(cool_dir, "SRR5005050.allValidPairs_500000.cool")
file2 <- file.path(cool_dir, "XM-QGY-MME-13_L7_1.allValidPairs_500000.cool")
name1 <- "SRR50"  # 修改样本名称为SRR50
name2 <- "XM13"   # 保持XM13不变

# 2. 加载 .cool 数据
hic1 <- load_contacts(file1, resolution=resolution, sample_name=name1,
                      balancing=TRUE, legacy=FALSE)
hic2 <- load_contacts(file2, resolution=resolution, sample_name=name2,
                      balancing=TRUE, legacy=FALSE)

# 3. 计算 compartment score（PC1）
CS_out <- compartment_score(list(hic1, hic2), ev=1, ref=NULL)

# 4. 绘图：在 chr1 上可视化
#visualise(
  #CS_out,
  #chr          = chr_to_plot,
  #track_type   = "compartment",
  #track_labels = c(name1, name2),
  #ylim_comp    = c(-2, 2),
  #title        = paste0("A/B Compartment (PC1) on ", chr_to_plot)
#)

# 得到的cs_out是这样的格式
    #| chrom | start   | end     | bin | SRR5050 | XM12  |
    #|-------|---------|---------|-----|---------|-------|
    #| chr1  | 0       | 500000  | 1   | 0.5     | 0.3   |
    #| chr1  | 500000  | 1000000 | 2   | -0.2    | -0.1  |
    #| chr1  | 1000000 | 1500000 | 3   | -0.7    | -0.6  | 
# 而我想得到的plot_df是这样的格式，于是进行5的处理
    # | bin | Mb  | Sample   | PC1  |
    # |-----|-----|----------|------|
    # | 1   | 0   | SRR5050  | 0.5  |
    # | 1   | 0   | XM12     | 0.3  |
    # | 2   | 0.5 | SRR5050  | -0.2 |
    # | 3   | 1.0 | SRR5050  | -0.7 |
# 1. 提取 chr1 上的 PC1 并构建长格式数据框
df <- CS_out$compart_scores %>%
  filter(chrom == chr_to_plot) %>%
  transmute(
    bin,
    Mb = start / 1e6,
    SRR50 = .data[[name1]],  # 使用变量名引用
    XM13 = .data[[name2]]    # 使用变量名引用
  )

plot_df <- df %>%
  pivot_longer(
    cols      = c("SRR50", "XM13"),
    names_to  = "Sample",
    values_to = "PC1"
  )

# 2. 用 ggplot2 绘图
p <- ggplot(plot_df, aes(x = Mb, y = PC1, color = Sample)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(
    values = c(SRR50 = "darkorange", XM13 = "steelblue"),  # 更新颜色
    labels = c(SRR50 = "SRR5005050", XM13 = "XM-QGY-MME-13_L7_1")  # 更新标签
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = seq(0, max(plot_df$Mb, na.rm=TRUE), by = 50)) +
  labs(
    title = "A/B Compartment (PC1) on chr1",
    x     = "Location chr1 (Mb)",
    y     = "Compartment score (PC1)",
    color = "Sample"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# 3. 保存图像
output_file <- paste0("compartment_", chr_to_plot, "_SRR50_vs_XM13.png")
ggsave(output_file,
       plot   = p,
       width  = 10,
       height = 4,
       dpi    = 300)

cat(paste0("Plot saved as: ", output_file, "\n"))

### 计算ab区室，即PC1分数相关性系数的散点图，这里先对chr1计算

# 用 Pearson 可以告诉你“两条曲线的振幅变化是否严格同步”
# 用 Spearman 可以告诉你“总体上大区室分（PC1 大）/小区室分（PC1 小）的位置趋势是否一致”
# 1. 提取 chr1 compartment 得分
df_chr1 <- CS_out$compart_scores %>%
  filter(chrom == "chr1") %>%
  transmute(
    SRR50 = .data[["SRR50"]],
    XM13 = .data[["XM13"]]
  ) %>%
  filter(!is.na(SRR50) & !is.na(XM13))

# 2. 计算相关性
pearson_corr  <- cor(df_chr1$SRR50, df_chr1$XM13, method = "pearson")
spearman_corr <- cor(df_chr1$SRR50, df_chr1$XM13, method = "spearman")

# 3. 创建更专业的散点图（不使用ggpubr）
library(ggplot2)

# 计算回归线数据
regression_line <- lm(XM13 ~ SRR50, data = df_chr1)
slope <- coef(regression_line)[2]
intercept <- coef(regression_line)[1]

# 创建更专业的散点图
p <- ggplot(df_chr1, aes(x = SRR50, y = XM13)) +
  # 使用六边形分箱图展示密度
  geom_hex(bins = 80, alpha = 0.9) +
  
  # 添加回归线
  geom_abline(slope = slope, intercept = intercept, 
              color = "red", linewidth = 1.2) +
  
  # 添加对角线参考线
  geom_abline(slope = 1, intercept = 0, 
              color = "gray40", linetype = "dashed", linewidth = 0.7) +
  
  # 添加相关系数标注
  annotate("text", x = min(df_chr1$SRR50, na.rm = TRUE), 
           y = max(df_chr1$XM13, na.rm = TRUE) * 0.95,
           label = paste0("Pearson r = ", round(pearson_corr, 3)),
           hjust = 0, vjust = 1, size = 5, color = "black") +
  
  annotate("text", x = min(df_chr1$SRR50, na.rm = TRUE), 
           y = max(df_chr1$XM13, na.rm = TRUE) * 0.88,
           label = paste0("Spearman ρ = ", round(spearman_corr, 3)),
           hjust = 0, vjust = 1, size = 5, color = "black") +
  
  # 添加数据点数量标注
  annotate("text", x = min(df_chr1$SRR50, na.rm = TRUE), 
           y = max(df_chr1$XM13, na.rm = TRUE) * 0.81,
           label = paste0("n = ", nrow(df_chr1), " bins"),
           hjust = 0, vjust = 1, size = 4.5, color = "gray30") +
  
  # 颜色方案
  scale_fill_viridis_c(option = "magma", trans = "log10", 
                       name = "Log10(count)") +
  
  # 主题设置
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray92"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = c(0.95, 0.05),  # 右下角
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "gray50"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    aspect.ratio = 1  # 保持正方形比例
  ) +
  
  # 标签和标题
  labs(
    title = "Compartment Score Correlation on chr1",
    subtitle = "SRR5005050 vs XM-QGY-MME-13_L7_1",
    x = "SRR5005050 PC1",
    y = "XM-QGY-MME-13_L7_1 PC1"
  )

# 4. 保存高质量图像
output_file <- "compart_corr_chr1_SRR50_vs_XM13.png"
ggsave(output_file, plot = p, width = 8, height = 8, dpi = 300, bg = "white")

# 5. 显示图表和统计信息
print(p)
cat(paste0(
  "chr1 Compartment Correlation:\n",
  "  Pearson r = ", round(pearson_corr, 4), "\n",
  "  Spearman ρ = ", round(spearman_corr, 4), "\n",
  "  Data points: ", nrow(df_chr1), "\n",
  "Plot saved as: ", output_file, "\n"
))

### 计算全部染色体chr的PC1相关性的散点图

# 1. 提取所有染色体的 compartment 得分
all_chr_data <- CS_out$compart_scores %>%
  filter(!is.na(.data[["SRR50"]]) & !is.na(.data[["XM13"]])) %>%
  select(SRR50 = .data[["SRR50"]], XM13 = .data[["XM13"]])

# 2. 计算全局相关性
pearson_corr  <- cor(all_chr_data$SRR50, all_chr_data$XM13, method = "pearson")
spearman_corr <- cor(all_chr_data$SRR50, all_chr_data$XM13, method = "spearman")
n_points <- nrow(all_chr_data)

# 3. 计算回归线数据
regression_line <- lm(XM13 ~ SRR50, data = all_chr_data)
slope <- coef(regression_line)[2]
intercept <- coef(regression_line)[1]

# 4. 创建专业级全局相关性图表
p <- ggplot(all_chr_data, aes(x = SRR50, y = XM13)) +
  # 使用六边形分箱图展示密度
  geom_hex(bins = 100, alpha = 0.85) +  # 适当减少分箱数量提高性能
  
  # 添加回归线
  geom_abline(slope = slope, intercept = intercept, 
              color = "red", linewidth = 1.3) +
  
  # 添加对角线参考线
  geom_abline(slope = 1, intercept = 0, 
              color = "gray40", linetype = "dashed", linewidth = 0.8) +
  
  # 添加相关系数标注
  annotate("text", x = min(all_chr_data$SRR50, na.rm = TRUE), 
           y = max(all_chr_data$XM13, na.rm = TRUE) * 0.95,
           label = paste0("Pearson r = ", round(pearson_corr, 3)),
           hjust = 0, vjust = 1, size = 5.5, color = "black", fontface = "bold") +
  
  annotate("text", x = min(all_chr_data$SRR50, na.rm = TRUE), 
           y = max(all_chr_data$XM13, na.rm = TRUE) * 0.88,
           label = paste0("Spearman ρ = ", round(spearman_corr, 3)),
           hjust = 0, vjust = 1, size = 5.5, color = "black", fontface = "bold") +
  
  # 添加数据点数量标注
  annotate("text", x = min(all_chr_data$SRR50, na.rm = TRUE), 
           y = max(all_chr_data$XM13, na.rm = TRUE) * 0.81,
           label = paste0("n = ", format(n_points, big.mark = ","), " bins"),
           hjust = 0, vjust = 1, size = 5, color = "gray30") +
  
  # 颜色方案
  scale_fill_viridis_c(
    option = "inferno", 
    trans = "log10", 
    name = "Log10(count)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  # 主题设置
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = c(0.95, 0.05),  # 右下角
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "gray50"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15)),
    plot.subtitle = element_text(hjust = 0.5, size = 16, color = "gray40", margin = margin(b = 20)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "gray30"),
    aspect.ratio = 1  # 保持正方形比例
  ) +
  
  # 标签和标题
  labs(
    title = "Genome-wide Compartment Score Correlation",
    subtitle = "SRR5005050 vs XM-QGY-MME-13_L7_1",
    x = "SRR5005050 PC1",
    y = "XM-QGY-MME-13_L7_1 PC1",
    caption = "Data includes all autosomal chromosomes"
  ) +
  
  # 坐标轴美化
  coord_cartesian(clip = "off")  # 确保标注不被裁剪

# 5. 保存高质量图像
output_file <- "compart_corr_all_chromosomes_SRR50_vs_XM13.png"
ggsave(output_file, plot = p, width = 10, height = 9, dpi = 350, bg = "white")

# 6. 显示统计信息
cat(paste0(
  "Global Compartment Correlation (All Chromosomes):\n",
  "  Pearson r = ", round(pearson_corr, 4), "\n",
  "  Spearman ρ = ", round(spearman_corr, 4), "\n",
  "  Data points: ", format(n_points, big.mark = ","), "\n",
  "Plot saved as: ", output_file, "\n"
))

### 计算每个染色体的Pearson和Spearman相关系数，然后用条形图展示。
    
    # 以下是步骤：
    # a. 从CS_out中提取所有染色体的compartment scores。
    # b. 按染色体分组，计算每个染色体的Pearson相关系数（以及可选的Spearman）。
    # c. 绘制条形图，x轴为染色体，y轴为相关系数。

# 1. 提取所有染色体的 compartment 得分
all_chr_data <- CS_out$compart_scores %>%
  filter(!is.na(.data[["SRR50"]]) & !is.na(.data[["XM13"]])) %>%
  select(chrom, bin, SRR50 = .data[["SRR50"]], XM13 = .data[["XM13"]])

# 2. 按染色体分组计算相关性
corr_by_chr <- all_chr_data %>%
  group_by(chrom) %>%
  filter(n() >= 10) %>%  # 只保留有足够数据点的染色体
  summarise(
    Pearson = cor(SRR50, XM13, method = "pearson"),
    Spearman = cor(SRR50, XM13, method = "spearman"),
    .groups = "drop"
  ) %>%
  mutate(chrom = factor(chrom, levels = unique(chrom)))  # 保持染色体原始顺序

# 3. 转换为长格式用于绘图
corr_long <- corr_by_chr %>%
  pivot_longer(
    cols = c(Pearson, Spearman),
    names_to = "Correlation_Type",
    values_to = "Correlation"
  )

# 4. 创建汇总图表
library(ggplot2)
p_summary <- ggplot(corr_long, aes(x = chrom, y = Correlation, fill = Correlation_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Correlation, 3)), 
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(Pearson = "#4E79A7", Spearman = "#F28E2B")) +
  labs(
    title = "Compartment Score Correlation by Chromosome",
    subtitle = "SRR5005050 vs XM-QGY-MME-13_L7_1",
    x = "Chromosome",
    y = "Correlation Coefficient",
    fill = "Correlation Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
  ) +
  ylim(min(0, min(corr_long$Correlation)) * 1.1, 1.1)  # 设置合适的Y轴范围

# 5. 保存汇总图表
ggsave("compart_corr_summary_by_chromosome.png", p_summary, width = 12, height = 7, dpi = 300)
cat("Summary plot of correlations by chromosome saved as: compart_corr_summary_by_chromosome.png\n")

