# ====================================================================
# Klebsiella Core-genome ML Tree - Rigorous Version (Advanced Gem Palette)
# Filename: plot_circular_tree.R
# ====================================================================

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(tidyr)
library(ape)
library(ggnewscale)

# Ensure your working directory is set correctly
# setwd("your/working/directory") 

# Read input files
tree <- read.tree("Core_ML_tree.treefile")
metadata <- read.table("species_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annotation <- read.table("annotation_cleaned.txt", header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, quote = "")

# Add reference strain info manually if missing
reference_info <- data.frame(
  Assembly.Accession = "reference",
  Organism.Name = "Klebsiella pneumoniae subsp. pneumoniae",
  stringsAsFactors = FALSE
)
metadata <- rbind(metadata, reference_info)

reference_anno <- data.frame(
  Assembly.Accession = "reference",
  Organism.Name = "Klebsiella pneumoniae subsp. pneumoniae",
  Collection.date = NA, Host = NA, Host_category = NA,
  stringsAsFactors = FALSE
)
annotation <- rbind(annotation, reference_anno)

# =========================================================================
# 🌟 Core modification: Advanced Gem Color Palette
# =========================================================================

# 1. Sector background colors: Low saturation soft tones to highlight tree branches
species_colors <- c(
  "Klebsiella africana" = "#F5B0CB",                                   # Soft Pink
  "Klebsiella quasivariicola" = "#EDBB99",                             # Light Terracotta
  "Klebsiella quasipneumoniae subsp. similipneumoniae" = "#F5CBA7",    # Light Peach
  "Klebsiella quasipneumoniae subsp. quasipneumoniae" = "#D2B4DE",     # Soft Lilac
  "Klebsiella variicola subsp. variicola" = "#AED6F1",                 # Glacier Blue
  "Klebsiella variicola subsp. tropica" = "#A2D9CE",                   # Soft Mint
  "Klebsiella pneumoniae subsp. pneumoniae" = "#F5B7B1"                # Rose Pink
)

# 2. Tip species point colors: High saturation dark colors for texture
tip_colors <- c(
  "Klebsiella africana" = "#D11586",                                   # Deep Magenta
  "Klebsiella quasivariicola" = "#D35400",                             # Rust Orange
  "Klebsiella quasipneumoniae subsp. similipneumoniae" = "#B9770E",    # Bronze Gold
  "Klebsiella quasipneumoniae subsp. quasipneumoniae" = "#6C3483",     # Amethyst
  "Klebsiella variicola subsp. variicola" = "#1F618D",                 # Sapphire
  "Klebsiella variicola subsp. tropica" = "#117A65",                   # Emerald
  "Klebsiella pneumoniae subsp. pneumoniae" = "#C0392B"                # Ruby Red
)
# =========================================================================

# Note: Outer ring year and host colors remain unchanged for visual consistency
date_colors <- c(
  "1994-2013" = "#4A90E2",
  "2014-2016" = "#50C878",
  "2017-2019" = "#F5A623",
  "2020-2022" = "#E94B3C"
)

host_colors <- c(
  "Homo sapiens" = "#E74C3C",
  "Animal" = "#F39C12",
  "Plant" = "#27AE60",
  "Environment" = "#3498DB"
)

annotation$Year <- suppressWarnings(as.numeric(annotation$Collection.date))
annotation$Year_group <- cut(annotation$Year,
                             breaks = c(-Inf, 2013, 2016, 2019, 2022, Inf),
                             labels = c("1994-2013", "2014-2016", "2017-2019", "2020-2022", "2023+"),
                             right = TRUE)

tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
tip_data <- left_join(tip_data, metadata %>% select(Assembly.Accession, Organism.Name), 
                     by = c("label" = "Assembly.Accession"))
colnames(tip_data)[2] <- "species"
tip_data <- left_join(tip_data, 
                     annotation %>% select(Assembly.Accession, Collection.date, 
                                          Host_category, Year, Year_group),
                     by = c("label" = "Assembly.Accession"))

tip_data$date_display <- ifelse(!is.na(tip_data$Year_group) & 
                                tip_data$Year_group != "2023+",
                                as.character(tip_data$Year_group), NA)
tip_data$host_display <- ifelse(
  !is.na(tip_data$Host_category) &
  tip_data$Host_category %in% names(host_colors),
  tip_data$Host_category, NA
)

# Calculate true tree depth
max_dist <- max(node.depth.edgelength(tree))

layer1_pos <- max_dist + max_dist * 0.08
layer2_pos <- max_dist + max_dist * 0.15
xlim_max <- max_dist + max_dist * 0.25

# Plot tree with actual branch lengths
p <- ggtree(tree, layout = "circular", color = "grey30", size = 0.3) %<+% tip_data +
  xlim(0, xlim_max)

# Sector background
species_groups <- tip_data %>% filter(!is.na(species)) %>% group_by(species) %>% summarise(tips = list(label), .groups = 'drop')
group_list <- setNames(species_groups$tips, species_groups$species)

for(sp in names(group_list)) {
  tips <- group_list[[sp]]
  if(length(tips) > 1) { node_id <- getMRCA(tree, tips) } else { node_id <- which(tree$tip.label == tips[1]) }
  if(!is.null(node_id) && !is.na(node_id)) {
    p <- p + geom_hilight(node = node_id, fill = species_colors[sp],
                         alpha = 0.3, to.bottom = TRUE, extendto = max_dist)
  }
}

# Draw connecting lines: link tips to a unified outer ring
p <- p + geom_tiplab(aes(label = ""), align = TRUE, linetype = "dotted", size = 0.2, color = "grey60", linesize = 0.3)

# Actual tip species points: allow overlap, add transparency
p <- p + 
  geom_tippoint(aes(fill = species), size = 2.5, shape = 21, 
                color = "black", stroke = 0.2, alpha = 0.8, show.legend = TRUE) +
  scale_fill_manual(
    values = tip_colors, name = "Species", na.value = "grey50",
    guide = guide_legend(override.aes = list(size = 4, shape = 21, color = "black", stroke = 0.3, alpha = 1), order = 1)
  )

# First annotation layer: Collection Date
p <- p + 
  new_scale_fill() +
  geom_tippoint(
    aes(fill = date_display, subset = !is.na(date_display), x = layer1_pos),
    size = 2.5, shape = 24, color = "black", stroke = 0.2, show.legend = TRUE
  ) +
  scale_fill_manual(
    values = date_colors, name = "Collection Year", na.translate = FALSE, breaks = names(date_colors),
    guide = guide_legend(override.aes = list(size = 4, shape = 24, color = "black", stroke = 0.3), order = 2)
  )

# Second annotation layer: Host Category
p <- p + 
  new_scale_fill() +
  geom_tippoint(
    aes(fill = host_display, subset = !is.na(host_display), x = layer2_pos),
    size = 2.5, shape = 22, color = "black", stroke = 0.2, show.legend = TRUE
  ) +
  scale_fill_manual(
    values = host_colors, name = "Host Category", na.translate = FALSE, breaks = names(host_colors),
    guide = guide_legend(override.aes = list(size = 4, shape = 22, color = "black", stroke = 0.3), order = 3)
  )

p <- p + 
  theme(
    legend.position = "right", legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9), legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.4, "cm"), legend.box = "vertical",
    legend.box.spacing = unit(0.5, "cm"), plot.margin = margin(10, 10, 10, 10)
  )

# Save plots
ggsave("klebsiella_ML_tree_rigorous.pdf", plot = p, width = 24, height = 24, units = "in", dpi = 300)
ggsave("klebsiella_ML_tree_rigorous.png", plot = p, width = 24, height = 24, units = "in", dpi = 300)
ggsave("klebsiella_ML_tree_rigorous.svg", plot = p, width = 24, height = 24, units = "in")

print(p)