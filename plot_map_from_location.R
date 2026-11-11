# ==================== 基于经纬度/地名的全球分布统计与绘图（修正版）====================
# 输入：地区总计.txt（包含 GCA 编号和经纬度/地点描述）
# 输出：国家统计表、全球分布地图（PNG、PDF、SVG）

# -------------------- 1. 加载必要的包 --------------------
required_packages <- c("maps", "ggplot2", "dplyr", "viridis", "svglite", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. 读取数据 --------------------
# 假设文件位于当前工作目录
data <- read.delim("地区总计.txt", header = TRUE, stringsAsFactors = FALSE)
cat("读取到", nrow(data), "条记录。\n")

# -------------------- 3. 定义解析函数 --------------------
# 将坐标字符串转换为十进制经纬度（返回 list(lat, lon)）
parse_coords <- function(coord_str) {
  # 如果是 NA 或空字符串，直接返回 NULL
  if (is.na(coord_str) || trimws(coord_str) == "") {
    return(NULL)
  }
  
  coord_str <- trimws(coord_str)
  
  # 情况1：带方向字母的格式，如 "27.939 N 120.695 E"
  # 使用 tryCatch 防止 str_detect 在 NA 上出错（但我们已经检查过）
  if (str_detect(coord_str, "[NSEW]")) {
    parts <- unlist(str_split(coord_str, "\\s+"))
    if (length(parts) == 4) {
      lat_val <- as.numeric(parts[1])
      lat_dir <- parts[2]
      lon_val <- as.numeric(parts[3])
      lon_dir <- parts[4]
      # 检查数值是否有效
      if (!is.na(lat_val) && !is.na(lon_val)) {
        lat <- ifelse(lat_dir == "S", -lat_val, lat_val)
        lon <- ifelse(lon_dir == "W", -lon_val, lon_val)
        return(list(lat = lat, lon = lon))
      }
    }
  }
  
  # 情况2：斜杠分隔的可能带负号格式，如 "-37.8462/144.9947"
  if (str_detect(coord_str, "[/]")) {
    parts <- unlist(str_split(coord_str, "/"))
    if (length(parts) == 2) {
      lat <- as.numeric(parts[1])
      lon <- as.numeric(parts[2])
      if (!is.na(lat) && !is.na(lon)) {
        return(list(lat = lat, lon = lon))
      }
    }
  }
  
  # 其他无法解析的返回 NULL
  return(NULL)
}

# 将地点字符串标准化为国家名称
get_country_from_location <- function(loc) {
  # 处理无效输入
  if (is.na(loc) || trimws(loc) %in% c("NA", "Not found", "", "0")) {
    return(NA_character_)
  }
  
  loc <- trimws(loc)
  
  # 尝试解析坐标
  coords <- parse_coords(loc)
  if (!is.null(coords)) {
    # 用 maps 包查询国家（注意：经度 x，纬度 y）
    country <- maps::map.where(x = coords$lon, y = coords$lat)
    if (length(country) > 0 && !is.na(country)) {
      # 多个结果时取第一个（例如边界可能返回两个）
      return(country[1])
    } else {
      return(NA_character_)
    }
  }
  
  # 处理冒号格式：取冒号前的部分
  if (str_detect(loc, ":")) {
    country <- trimws(str_split(loc, ":")[[1]][1])
    return(country)
  }
  
  # 直接作为国家名返回
  return(loc)
}

# -------------------- 4. 解析所有地点并获取国家 --------------------
raw_countries <- sapply(data$Latitude.and.longitude, get_country_from_location)
names(raw_countries) <- NULL

# 查看未能解析的条目（调试用）
unresolved <- data[is.na(raw_countries), ]
if (nrow(unresolved) > 0) {
  cat("以下地点未能解析为国家，将被忽略：\n")
  print(unresolved)
}

# -------------------- 5. 国家名称标准化（适配 maps 包） --------------------
# 定义映射表
country_map <- c(
  "United States"     = "USA",
  "USA"               = "USA",
  "U.S.A."            = "USA",
  "United Kingdom"    = "UK",
  "UK"                = "UK",
  "Republic of Korea" = "South Korea",
  "South Korea"       = "South Korea",
  "Hong Kong"         = "China",   # 香港归入中国
  "Taiwan"            = "China"     # 台湾归入中国
)

standardize_country <- function(cnt) {
  if (is.na(cnt)) return(NA_character_)
  # 去除多余空格
  cnt <- trimws(cnt)
  # 应用映射
  if (cnt %in% names(country_map)) {
    return(country_map[cnt])
  }
  # 否则原样返回（需确保 maps 中存在）
  return(cnt)
}

countries_std <- sapply(raw_countries, standardize_country)

# -------------------- 6. 统计每个国家的数量 --------------------
country_counts <- table(countries_std, useNA = "ifany")
cat("\n国家统计（包括 NA）：\n")
print(country_counts)

# 去除 NA 的统计（NA 代表未能解析的条目）
valid_counts <- country_counts[!is.na(names(country_counts))]
cat("\n有效国家数量：", sum(valid_counts), "株\n")

# 转换为数据框，便于后续绘图
country_df <- data.frame(
  Country_Std = names(valid_counts),
  genomes = as.numeric(valid_counts),
  stringsAsFactors = FALSE
)

# 检查总数量是否与原始记录数一致（忽略 NA 后）
cat("原始总记录数：", nrow(data), "\n")
cat("有效解析数 ：", sum(valid_counts), "\n")

# -------------------- 7. 准备地图数据 --------------------
world_map <- ggplot2::map_data("world")

# 合并数据
plot_data <- world_map %>%
  left_join(country_df, by = c("region" = "Country_Std")) %>%
  mutate(
    genomes = ifelse(is.na(genomes), 0, genomes),
    log_genomes = ifelse(genomes > 0, log10(genomes), NA)
  )

# -------------------- 8. 绘制地图（色盲友好配色）--------------------
p <- ggplot() +
  geom_polygon(
    data = plot_data,
    aes(x = long, y = lat, group = group, fill = log_genomes),
    color = "white", size = 0.2
  ) +
  scale_fill_viridis_c(
    name = "Genome Count\n(log10 scale)",
    option = "viridis",
    na.value = "grey90",
    breaks = log10(c(1, 3, 10, 30, 100)),
    labels = c("1", "3", "10", "30", "100")
  ) +
  coord_fixed(ratio = 1.3) +
  labs(
    title = "Global Distribution of Sequenced KpSC Genomes",
    subtitle = paste0(
      "Total genomes: ", sum(valid_counts),
      " | Countries/Regions with data: ", nrow(country_df)
    ),
    caption = "Data source: NCBI GenBank | Visualization: R (maps + ggplot2)"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    plot.caption = element_text(color = "grey60", hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

print(p)

# -------------------- 9. 保存地图 --------------------
output_dir <- "./figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# PNG
png_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.png")
ggsave(png_file, plot = p, width = 14, height = 7, dpi = 600)
cat("\n✓ PNG 地图已保存至：", png_file, "\n")

# PDF
pdf_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.pdf")
ggsave(pdf_file, plot = p, width = 14, height = 7, device = cairo_pdf)
cat("✓ PDF 地图已保存至：", pdf_file, "\n")

# SVG
svg_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.svg")
ggsave(svg_file, plot = p, width = 14, height = 7, device = "svg")
cat("✓ SVG 地图已保存至：", svg_file, "\n")

# -------------------- 10. 输出国家统计表（可选）--------------------
write.table(country_df, file = "country_counts.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("国家统计表已保存至 country_counts.txt\n")