# ==================== Global Distribution Analysis and Mapping (Revised) ====================
# Input: location_total.txt (contains GCA IDs and latitude/longitude or location descriptions)
# Output: Country counts table, global distribution map (PNG, PDF, SVG)

# -------------------- 1. Load required packages --------------------
required_packages <- c("maps", "ggplot2", "dplyr", "viridis", "svglite", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. Read data --------------------
# Assuming the file is in the current working directory
# Ensure your input file is named "location_total.txt"
data <- read.delim("location_total.txt", header = TRUE, stringsAsFactors = FALSE)
cat("Read", nrow(data), "records.\n")

# -------------------- 3. Define parsing functions --------------------
# Convert coordinate string to decimal latitude/longitude (returns list(lat, lon))
parse_coords <- function(coord_str) {
  # If NA or empty string, return NULL
  if (is.na(coord_str) || trimws(coord_str) == "") {
    return(NULL)
  }
  
  coord_str <- trimws(coord_str)
  
  # Case 1: Format with direction letters, e.g., "27.939 N 120.695 E"
  if (str_detect(coord_str, "[NSEW]")) {
    parts <- unlist(str_split(coord_str, "\\s+"))
    if (length(parts) == 4) {
      lat_val <- as.numeric(parts[1])
      lat_dir <- parts[2]
      lon_val <- as.numeric(parts[3])
      lon_dir <- parts[4]
      if (!is.na(lat_val) && !is.na(lon_val)) {
        lat <- ifelse(lat_dir == "S", -lat_val, lat_val)
        lon <- ifelse(lon_dir == "W", -lon_val, lon_val)
        return(list(lat = lat, lon = lon))
      }
    }
  }
  
  # Case 2: Slash-separated negative format, e.g., "-37.8462/144.9947"
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
  
  return(NULL)
}

# Standardize location string to country name
get_country_from_location <- function(loc) {
  if (is.na(loc) || trimws(loc) %in% c("NA", "Not found", "", "0")) {
    return(NA_character_)
  }
  
  loc <- trimws(loc)
  
  coords <- parse_coords(loc)
  if (!is.null(coords)) {
    country <- maps::map.where(x = coords$lon, y = coords$lat)
    if (length(country) > 0 && !is.na(country)) {
      return(country[1])
    } else {
      return(NA_character_)
    }
  }
  
  # Handle colon format: take the part before the colon
  if (str_detect(loc, ":")) {
    country <- trimws(str_split(loc, ":")[[1]][1])
    return(country)
  }
  
  return(loc)
}

# -------------------- 4. Parse all locations and get countries --------------------
raw_countries <- sapply(data$Latitude.and.longitude, get_country_from_location)
names(raw_countries) <- NULL

# Check unresolved entries (for debugging)
unresolved <- data[is.na(raw_countries), ]
if (nrow(unresolved) > 0) {
  cat("The following locations could not be resolved to a country and will be ignored:\n")
  print(unresolved)
}

# -------------------- 5. Standardize country names (Adapt for maps package) --------------------
country_map <- c(
  "United States"     = "USA",
  "USA"               = "USA",
  "U.S.A."            = "USA",
  "United Kingdom"    = "UK",
  "UK"                = "UK",
  "Republic of Korea" = "South Korea",
  "South Korea"       = "South Korea",
  "Hong Kong"         = "China",   
  "Taiwan"            = "China"     
)

standardize_country <- function(cnt) {
  if (is.na(cnt)) return(NA_character_)
  cnt <- trimws(cnt)
  if (cnt %in% names(country_map)) {
    return(country_map[cnt])
  }
  return(cnt)
}

countries_std <- sapply(raw_countries, standardize_country)

# -------------------- 6. Count genomes per country --------------------
country_counts <- table(countries_std, useNA = "ifany")
cat("\nCountry counts (including NA):\n")
print(country_counts)

valid_counts <- country_counts[!is.na(names(country_counts))]
cat("\nTotal valid countries/regions:", sum(valid_counts), "genomes\n")

country_df <- data.frame(
  Country_Std = names(valid_counts),
  genomes = as.numeric(valid_counts),
  stringsAsFactors = FALSE
)

cat("Total raw records:", nrow(data), "\n")
cat("Valid records parsed:", sum(valid_counts), "\n")

# -------------------- 7. Prepare map data --------------------
world_map <- ggplot2::map_data("world")

plot_data <- world_map %>%
  left_join(country_df, by = c("region" = "Country_Std")) %>%
  mutate(
    genomes = ifelse(is.na(genomes), 0, genomes),
    log_genomes = ifelse(genomes > 0, log10(genomes), NA)
  )

# -------------------- 8. Plot map (colorblind-friendly palette) --------------------
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

# -------------------- 9. Save map --------------------
output_dir <- "./figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

png_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.png")
ggsave(png_file, plot = p, width = 14, height = 7, dpi = 600)
cat("\n[OK] PNG map saved to:", png_file, "\n")

pdf_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.pdf")
ggsave(pdf_file, plot = p, width = 14, height = 7, device = cairo_pdf)
cat("[OK] PDF map saved to:", pdf_file, "\n")

svg_file <- file.path(output_dir, "KpSC_Global_Distribution_Map.svg")
ggsave(svg_file, plot = p, width = 14, height = 7, device = "svg")
cat("[OK] SVG map saved to:", svg_file, "\n")

# -------------------- 10. Output country counts table --------------------
write.table(country_df, file = "country_counts.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Country count table saved to country_counts.txt\n")