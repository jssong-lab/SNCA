
############################################################################
## Author  : Porter B. Howland                                            ##
## Email   : pbh2@illinois.edu                                            ##
## Date    : September 7, 2025                                            ##
## Updated :                                                              ##
##                                                                        ##
## Description :                                                          ##
##                                                                        ##
## Creates a heatmap using genes in selected GO categories from the       ##
## analysis of data in the IGH experiment.                                ##
##                                                                        ##
############################################################################


library(scales)
library(pheatmap)
library(gridExtra)
library(grid)
library(gtable)


## read in all genes for used GO categories
df <- read.csv("path_to_gene_data.csv")
rownames(df) <- df$gene
df <- df[, !(names(df) %in% c("gene"))]
colnames(df) <- paste("patient", 1:4)


## colormap adjustments ##

# True limits
MIN <- min(df)
MAX <- max(df)

ABS_MAX <- abs(MIN)  # Symmetric upper limit for palette, not for data

# Total number of colors (more = smoother)
n <- 100

# Equal-length negative and positive palettes based on ABS_MAX
n_each <- n / 2
neg_colors <- colorRampPalette(c("red", "orange", "white"))(n_each)
pos_full_colors <- colorRampPalette(c("white", "blue"))(n_each)

# For the positive portion, we truncate the gradient to just the MAX range
# Compute how much of the positive gradient we want
truncate_ratio <- MAX / ABS_MAX
n_pos <- max(1, round(n_each * truncate_ratio))  # at least 1 color

# Use only the first n_pos colors (light blues) from the full blue gradient
pos_colors <- pos_full_colors[1:n_pos]

# Combine palettes
my_palette <- c(neg_colors, pos_colors)

# Define breaks: from MIN to MAX, length = colors + 1
my_breaks <- seq(MIN, MAX, length.out = length(my_palette) + 1)

## end colormap adjustment ##

p1 <- pheatmap(
  df[c(1:29),],
  color = my_palette,
  breaks = my_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 24,
  cellheight = 12,
  legend = FALSE,
  silent = TRUE
)

p2 <- pheatmap(
  df[c(30:length(df)),],
  color = my_palette,
  breaks = my_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 24,
  cellheight = 12,
  legend = FALSE,
  silent = TRUE
)

## extract the legend
plegend <- pheatmap(
  matrix(seq(MIN, MAX, length.out = n), ncol = 1),
  colors = my_palette,
  breaks = my_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  silent = TRUE
)

# Extract grobs
g1 <- p1$gtable
g2 <- p2$gtable
g_legend <- plegend$gtable$grobs[[which(sapply(plegend$gtable$grobs, function(x) x$name) == "legend")]]

  
plegend$gtable

# Arrange plots: two heatmaps and a single legend
grid.newpage()
grid.arrange(g1, g2, g_legend, ncol = 3, widths = c(4, 4, 1))
