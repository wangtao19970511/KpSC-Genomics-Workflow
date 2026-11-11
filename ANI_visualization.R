# ===================== 1. 加载必要的R包 =====================
# 如果未安装，请先运行：install.packages(c("tidyverse", "scales", "ggrepel"))
library(tidyverse)      # 包含 ggplot2、dplyr、readr 等
library(scales)         # 坐标轴格式化
library(ggrepel)        # 可选，用于标注特殊样本

# ===================== 2. 读取fastANI结果 =====================
# 请将 "data/fastani_output.txt" 替换为您本地实际的fastANI结果文件路径
ani_raw <- read_tsv("data/fastani_output.txt",
                    col_names = c("Query", "Reference", "ANI", "Matches", "Total"),
                    col_types = "ccddd")

# ===================== 3. 数据处理 =====================
ani_data <- ani_raw %>%
  mutate(Sample = basename(Query)) %>%          # 提取文件名（不含路径）
  arrange(ANI) %>%                               # 按ANI升序排列
  mutate(
    Rank = row_number(),                         # 生成序号（1~N）
    Category = case_when(
      ANI >= 95 ~ "≥95% (同种)",
      ANI >= 90 ~ "90-95% (同属)",
      TRUE      ~ "<90% (不同属)"
    ),
    Category = factor(Category, levels = c("≥95% (同种)", "90-95% (同属)", "<90% (不同属)"))
  )

# 计算统计信息用于副标题
total_n <- nrow(ani_data)
n_ge95 <- sum(ani_data$ANI >= 95)
n_90_95 <- sum(ani_data$ANI >= 90 & ani_data$ANI < 95)
n_lt90 <- sum(ani_data$ANI < 90)
median_ani <- median(ani_data$ANI)

# ===================== 4. 设置坐标轴范围 =====================
# 根据数据自动调整，也可手动指定
x_min <- floor(min(ani_data$ANI) * 10) / 10   # 向下取整到0.1
x_max <- ceiling(max(ani_data$ANI) * 10) / 10 # 向上取整到0.1
# 如果数据集中在较窄区间，可适当扩展边界
if (x_max - x_min < 2) {
  x_min <- x_min - 0.2
  x_max <- x_max + 0.2
}
# 设置刻度间隔（根据范围智能选择）
if (x_max - x_min <= 3) {
  break_step <- 0.5
} else if (x_max - x_min <= 5) {
  break_step <- 1
} else {
  break_step <- 2
}
x_breaks <- seq(floor(x_min), ceiling(x_max), by = break_step)

# ===================== 5. 绘制优化后的抖动点图 =====================
p <- ggplot(ani_data, aes(x = ANI, y = Rank, color = Category)) +
  # 抖动点：垂直抖动0.8，透明度0.6，固定随机种子保证可重复
  geom_point(size = 1.8, alpha = 0.6,
             position = position_jitter(width = 0, height = 0.8, seed = 123)) +
  # 95% 物种界定阈值线（红色虚线）
  geom_vline(xintercept = 95, linetype = "dashed", color = "#d73027", size = 0.8) +
  # 中位线（灰色实线，可选）
  geom_vline(xintercept = median_ani, linetype = "solid",
             color = "gray40", size = 0.5, alpha = 0.8) +
  # 自定义颜色（色盲友好型配色）
  scale_color_manual(
    values = c("≥95% (同种)" = "#d73027",
               "90-95% (同属)" = "#4575b4",
               "<90% (不同属)" = "#313695"),
    name = NULL  # 隐藏图例标题（若需显示可改为 "ANI 分组"）
  ) +
  # 横坐标设置
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = x_breaks,
    expand = c(0.02, 0),               # 两端添加2%边距
    labels = number_format(accuracy = 0.1)  # 保留一位小数
  ) +
  # 纵坐标：样本序号
  scale_y_continuous(
    limits = c(0, total_n + 5),         # 留出顶端空间
    breaks = seq(0, total_n, by = 50),
    expand = c(0.01, 0)
  ) +
  # 坐标轴标签与标题
  labs(
    x = "Average Nucleotide Identity (%)",
    y = "Isolates (ordered by increasing ANI)",
    title = "Average Nucleotide Identity (ANI) to Reference Strain",
    subtitle = paste0("Total isolates: ", total_n,
                      "  |  ≥95%: ", n_ge95,
                      "  |  90-95%: ", n_90_95,
                      "  |  <90%: ", n_lt90,
                      "  |  Median ANI: ", round(median_ani, 2), "%"),
    caption = "Red dashed line: 95% species cutoff"
  ) +
  # 经典主题，无网格线，加边框
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    # 图例放在图内右下角
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.key = element_blank(),
    legend.text = element_text(size = 9),
    # 标题
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0, margin = margin(b = 10)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  # 图例中点放大，更易看清
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# 显示图形
print(p)

# ===================== 6. 保存图片 =====================
# 保存为PDF（矢量图，适合投稿）
ggsave("ANI_plot_final.pdf", plot = p, width = 7, height = 9, units = "in")

# 保存为PNG（位图，300 DPI，适合预览）
ggsave("ANI_plot_final.png", plot = p, width = 7, height = 9, units = "in", dpi = 300)

# 保存为SVG（可缩放矢量图，适合网页或再编辑）
ggsave("ANI_plot_final.svg", plot = p, width = 7, height = 9, units = "in")