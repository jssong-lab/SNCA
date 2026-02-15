
############################################################################
## Author  : Porter B. Howland                                            ##
## Email   : pbh2@illinois.edu                                            ##
## Date    : March 6, 2025                                                ##
## Updated : June 4, 2025                                                 ##
##                                                                        ##
## Description :                                                          ##
##                                                                        ##
##   Creates heatmaps of highest magnitude average log2 fold change genes ##
##   within a gene ontology category, sorted by control pseudobulk        ##
##   expression.  Unbiased clustering is performed on the full category   ##
##   and used for the displayed dendrogram. Updated column annotation to  ##
##   match colors in Seurat.                                              ##
############################################################################


library(Seurat)
library(pheatmap)
library(ggplot2)
library(grid)
library(scales)


###################################
## Global Parameters (constants) ##
###################################

## 'groups' correspond to the Seurat objects 'orig.ident' idents.

top.n <- 15
nblocks <- 10
groups <- c("k14-control", "a53t-hyper", "a53t-hypo")
labels <- c("control", "hyper", "hypo")
pseudobulk.threshold <- 1e-6


##########################
## Function Definitions ##
##########################

## --- coarse.grain
##
## Takes a single list of length N, and divides it uniformly into
## nblocks, and computes the average for each block.  If nblocks does
## not divide N, the final block can be smaller than N / nblocks.
## Returns a list length nblocks containing the averaged values.
##
coarse.grain <- function(xs) {
  
  block.stride <- ceiling(length(xs) / nblocks)
  
  block.avgs <- list()
  for (i in seq(1, length(xs), by=block.stride)) {
    block.avgs <- c(block.avgs, mean(xs[c(i:(i+block.stride-1))], na.rm = TRUE))
  }
  
  return(block.avgs)
}


## --- make.block
##
## Given a Seurat object and a list of genes, returns a matrix of
## those genes coarse-grained into nblocks for each experimental
## group.  Genes are listed in rows and groups in chunks of nblocks in
## columns.
##
make.block <- function(sobj, genes) {
  
  xs <- list()
  for (group in groups) {
    x <- FetchData(sobj, vars = genes, layer = paste("data", group, sep='.'))
    xs <- rbind(xs, sapply(x, coarse.grain))
  }
  
  # re-shape into a matrix
  mat <- t(matrix(unlist(xs),
                  nrow = length(groups) * nblocks,
                  ncol = length(genes)))
  
  # assign row and column names
  rownames(mat) <- genes
  colnames(mat) <- paste(sapply(labels, \(x) rep(x, nblocks)), c(1:nblocks), sep = '-')
  
  return(mat)
}



## --- rank.order
##
## Given a dataframe df of numerical values, it returns a copy of df
## with the value replaced by the scaled rank (i.e., the location in
## the list if all values were sorted, divided by the total number of
## values).  The original dataframe df is unaltered.
##
rank.order <- function(df) {
  
  ranks <- data.frame(x=c(), y=c(), val=c())
  for (x in rownames(df)) {
    for (y in colnames(df)) {
      ranks <- rbind(ranks, data.frame(x=x, y=y, val=df[x, y]))
    }
  }

  ranks <- ranks[order(ranks$val),]
  ranks$val <- 1:dim(ranks)[[1]] / dim(ranks)[[1]]
  
  ranked.df <- data.table::copy(df)
  for (i in 1:nrow(ranks)) {
    ranked.df[ranks[i, 1], ranks[i, 2]] <- ranks[i, 3]
  }
  
  return(ranked.df)
}


## Anonymous functions - assign column annotations and re-order
## dataframe
mrepN <- \(xs) unlist(lapply(xs, \(x) rep(x, nblocks)))
order_by_control <- \(df) rev(order(apply(df[,paste0("control-", 1:10)], 1, mean)))


##########
## Main ##
##########


## Open the Seurat object and select only the melanocyte subcluster
## (here k = 4).

## Load the Seurat object (sobj) from the R Data Structure file
# sobj <- readRDS("../preprocessed_counts_2000_vf.rds")
sobj <- readRDS("your/path/to/saved.rds")
Idents(sobj) <- "seurat_clusters"
sobj <- subset(sobj, ident = 4)
Idents(sobj) <- "orig.ident"


## Collect file names where gene ontology (GO) data is stored.  These
## are currated manually.  Output file names are derived from the
## input category names.
go.dir <- file.path(".", "merged-categories")
filenames <- list.files(go.dir, pattern="*.txt", full.names = TRUE)
outfiles <- paste(
  sapply(
    strsplit(filenames, '/'),
    \(x) strsplit(tail(x, 1), '.', fixed=TRUE)[[1]])[1,],
  "pdf", sep='.')


## Read in adjusted p-values as a tab-separated file and re-order the
## columns to display in the desired order in the heatmap.  This data
## currated (via Python scripts) from the output files of Seurat's
## FindMarkers() function for DEGs within the melanocyte cluster
## between experimental groups.
##
## Note: this file is just a reformatting of data (up to binarization)
## contained in Seurat's differential expression output.
##
df.pvals <- read.table(
  file = "merged_adj_pvals.tsv",
  sep = '\t',
  header = TRUE,
  stringsAsFactors = TRUE,
  row.names = 1)
df.pvals <- df.pvals[,c("hypo.vs.hyper", "hypo.vs.control", "hyper.vs.control")]


## Filter gene list by pseudobulk expression.  The threshold is set
## low in order to include lowly expressed genes but exclude those
## with zero base expression in the given cluster.
pb <- as.data.frame(PseudobulkExpression(sobj,
                                         features = Features(sobj),
                                         layer = "data")$RNA)
pb <- pb[apply(pb, 1, max) > pseudobulk.threshold,]
all.genes <- rownames(pb)

column_colors <- hue_pal()(3)

## Annotation colors for statistical significance which is passed to
## the pheatmap function.
ann_colors <- list(
  groups = c("control" = column_colors[1], "hyper" = column_colors[2], "hypo" = column_colors[3]),
  hyper.vs.control = c("white", "black"),
  hypo.vs.control = c("white", "black"),
  hypo.vs.hyper = c("white", "black")
)

## Sorted avg_log2FC lists for up and down regulated genes.  These
## lists are currated from the DEGs files.
gs.up.sorted <- readLines("genes_hyper_or_hypo_up_sorted.txt")
gs.down.sorted <- readLines("genes_hyper_or_hypo_down_sorted.txt")
gs.all.sig <- readLines("genes_all_sig_03.txt")

## GO category names (for titles).
categories <- c("Pigmentation", "Neurodegeneration", "Lysosome",
                "Autophagy", "Senescence", "Unfolded Protein Response", "Mitochondrion")

## make pheatmaps for each category
for (i in c(1:length(filenames))) {
  
  ## Excluding genes with zero expression (will cause an error), make
  ## coarse-grained block for genes in the GO category.
  fs <- intersect(all.genes, readLines(filenames[[i]]))
  df <- make.block(sobj, fs)

  ## Build dataframe for column annotations (colors experimental
  ## groups at top of heatmap below dendrogram).
  df.annotation.cols <- as.data.frame(mrepN(labels),
                                      row.names = colnames(df))
  colnames(df.annotation.cols) <- c("groups")
  
  ## Hierarchical cluster columns using *all* genes in the GO category
  ## (unbiased).  The heatmap is not shown but saved so that the
  ## dendrogram can be display on the filtered figure.
  hm <- pheatmap(df, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "average", silent = TRUE)
  
  ## subset the dataframe rows by taking only those statistically
  ## significant genes in the category.
  df <- df[rownames(df) %in% gs.all.sig,]

  ## Split the dataframe into down and up regulated genes, and sort by
  ## magnitude of avg_log2FC.  gs.{down,up}.sorted lists read as
  ## input, and constructed by processing the DEGs files.
  df.down <- df[gs.down.sorted[gs.down.sorted %in% rownames(df)],]
  df.up <- df[gs.up.sorted[gs.up.sorted %in% rownames(df)],]
  
  ## Take top n (if available) avg log2 FC genes (down / up).
  num.down <- min(top.n, nrow(df.down))
  num.up <- min(top.n, nrow(df.up))
  
  df.down <- head(df.down, num.down)
  df.up <- head(df.up, num.up)
  
  ## Sort down by descending control pseudobulk among those selected
  ## genes; up by descending as well.
  df.down <- df.down[order_by_control(df.down),]
  df.up <- df.up[order_by_control(df.up),]

  ## Bind both dataframes together (down on top, up on bottom) and
  ## rank order.
  df <- rbind(df.down, df.up)
  df <- rank.order(df)

  max.val <- max(df)
  
  ## show subset of rows on new pheatmap
  pheatmap(df,
           # main=categories[[i]],
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           annotation_legend_row = FALSE,
           annotation_legend = FALSE,
           cluster_rows = FALSE,
           gaps_row = c(num.down),  # rep to adjust amount of spacing
           cluster_cols = hm$tree_col,
           annotation_col = df.annotation.cols,
           annotation_colors = ann_colors,
           cutree_cols = 3,
           annotation_row = df.pvals[rownames(df),],
           show_colnames = FALSE,
           color = colorRampPalette(c(rgb(23/255, 127/255, 230/255), "cyan", "green", "yellow", "red"))(200),
           breaks = seq(0, max.val, max.val / (200 - 1)),
           filename=paste0("hm_down_up_desc_FC_sorted_ranked_unannotated_", outfiles[[i]]),
           cellheight = 12,
           cellwidth = 12,
           width = 8,
           height = 7
           )
}


