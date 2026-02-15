
############################################################################
## Author  : Porter B. Howland                                            ##
## Email   : pbh2@illinois.edu                                            ##
## Date    : August 31, 2024                                              ##
## Updated : September 5, 2024                                            ##
##                                                                        ##
## Description :                                                          ##
##                                                                        ##
##     Using the counts matrices ouputted from the 10x genomics pipeline, ##
##     read in multiple datasets and merge them.  Standard preprocessing  ##
##     is applied (i.e., discarding cells with high mitochondria counts,  ##
##     small number of features, etc.).  Various postprocessing is then   ##
##     done including PCA, clustering, and differential expression.  The  ##
##     bulk of the computation serialized to an *.rds file so that it can ##
##     simply be loaded again with for quicker downstream analysis.       ##
############################################################################


## load packages required for plotting (patchwork, ggplot2) and
## analysis (dplyr, Seurat)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)


## File for serialized results of umap/clustering and output directory
## name (to keep workspace tidy).  The RDS filename should probably match
## the one in the else block of the !file.exists(rds_filename) statement.
## This is purely to save repeating expensive computation.
rds_filename <- "name_for_serialized_structure.rds"
outdir <- "directory_for_all_figures_and_output_files"


## Check to see if the directory already exists and warn user...
## ... we don't want to clobber pre-existing results.
if (dir.exists(outdir)) {
    print(paste("The output directory '", outdir,
                "' already exists!  Exiting...", sep = ''))
    stop(1)
} else {
    dir.create(outdir)
}



## Check if a prepared RDS file exists already and load it to save
## expensive computation.
if (!file.exists(rds_filename)) {

    ## These blocks are specific to the file structure on our local
    ## cluster and will need to be adapted to take as intput the output
    ## files from the 10x genomics pipeline.
    
    ## the order of projects, (directory) identifiers, and prefixes should
    ## all coincide
    projects <- list("k14-control", "a53t-hyper", "a53t-hypo")
    identifiers <- list("H7TNMDRX5", "H7KWJDRX5", "H7KTGDRX5")
    prefixes <- list("control", "hyper", "hypo")
    n <- length(projects)

    ## thresholds for count selection - see documentation for
    ## 'CreateSeuratObject'
    min_cells <- 3
    min_features <- 200

    ## construct filepaths (these will be different in each experiment)
    filepaths = list()
    for (i in 1:n) {
        filepaths[[i]] <- file.path(
            "/", "home", "groups", "song", "songlab2", "shared", "NicholasTheodosakisData",
            "scRNA-seq4", "Analysis", "cell-ranger-outputs",
            paste("count", identifiers[[i]], "no_XY_with_SNCA", sep = '_'),
            "outs", "filtered_feature_bc_matrix"
        )
    }    
    
    ## read in (sparse) counts matrices and create seurat objects
    ##sobs = rep(list(NULL), n)
    sobs = list()    
    for (i in 1:n) {
        sobs[[i]] <- CreateSeuratObject(
            counts = Read10X(data.dir = filepaths[[i]]),
            project = projects[[i]],
            min.cells = min_cells,
            min.features = min_features
        )
    }

    ## merge datasets together
    df <- merge(
        x = sobs[[1]],
        y = sobs[-1],
        add.cell.ids = prefixes,
        project = "merged_samples"
    )

    ## perform normalization and data scaling
    df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
    df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
    df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
    df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
    df <- ScaleData(df, features = rownames(df))

    ## PCA and clustering
    df <- RunPCA(df, features = VariableFeatures(object = df))
    df <- RunUMAP(df, dims = 1:30)
    df <- FindNeighbors(df, dims = 1:30)
    df <- FindClusters(df, resolution = 0.1)

    ## save dataframe (R Data Structure) to file
    saveRDS(df, file = rds_filename)
    
} else {

    ## read the RDS file if it is already available
    df <- readRDS(rds_filename)

}


################################################################
## additional postprocessing (comment / uncomment as desired) ##
################################################################


## ** NOTE **
##
## The Idents() method sets the identity class for the Seurat object
## (dataframe, df).  This is used by other Seurat functions for
## plotting and analysis, usually through the variables ident.1 and
## ident.2, but sometimes only implicitely.  At risk of being verbose,
## each independent task below has the correct identity class
## explicitely set at the onset.
##
## The available identity classes correspond to the columns of the
## Seurat object, which can be found by, e.g., 'head(df)'.  The values
## taken by any particular identity class can be found by calling
## 'levels(df)' with Idents() set to the desired class.  The number of
## cells/reads in each level of a class can be found using
## table(df$orig.ident), or replacing 'orig.ident' with the desired
## identity class label.
##
## You can also use the 'group.by' argument in many functions rather
## than explicitely changing the identity class.
##


## differential expression of cluster k to others
Idents(df) <- "seurat_clusters"
joined_df <- JoinLayers(df)
for (k in as.integer(levels(df))) {
    sink(file = file.path(outdir, paste("de_cluster_", k, "_vs_all_else.txt", sep = '')))
    print(FindMarkers(joined_df, ident.1 = k))
    sink()
}
remove(joined_df)

## plot clusters
Idents(df) <- "seurat_clusters"
DimPlot(df, reduction = "umap", alpha = 0.2, label = TRUE)
ggsave(file.path(outdir, "cluster.png"))
ggsave(file.path(outdir, "cluster.pdf"))


##
## zoom in on melanocyte cluster (specific to SNCA mouse experiment)
##
Idents(df) <- "seurat_clusters"
df_mel <- subset(df, idents = 3)

Idents(df_mel) <- "orig.ident"
DimPlot(df_mel, reduction = "umap", alpha = 0.2, label = FALSE) + xlim(c(3, 10)) + ylim(c(-15, -7.5)) + NoLegend()
ggsave("melanocyte_cluster.png")

FeaturePlot(df_mel, features = "SNCA", alpha = 0.2) + xlim(c(3, 10)) + ylim(c(-15, -7.5)) + NoLegend()
ggsave("melanocyte_SNCA.png")



## dotplot with specific features
Idents(df) <- "seurat_clusters"
fs = c("Cdh1", "Krt1", "Krt5", "Krt10", "Krt14", "Trp63", "Dsc1",
       "Ivl", "Krt2", "Krt17", "Sostdc1", "Fst", "Cdk1", "Ccnb1",
       "Ccnb2", "Ccna2", "Bub1", "Mitf", "Tyr", "Dct", "Mlana",
       "Tyrp1", "Pax3", "Kit", "Sox10", "Snca", "SNCA", "Ptprc",
       "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1", "Cd7", "Cd28",
       "Cd44", "Cd48", "Cd69", "Cd82", "Foxp3", "Il2ra", "Ctla2a",
       "Ctla2b", "Ctla4", "Trdc", "Trdv4", "Trgv2", "Lag3", "Cd14",
       "Cd83", "Cd80", "Cd86", "Cd207", "Ccr2", "Il16", "Il18",
       "Cd274", "Cd200", "Cd74", "Prex2", "Plcb4", "Mef2c", "Cadm1",
       "H2-K1", "H2-Q7", "H2-D1", "H2-Aa", "H2-Ab1", "H2-Eb1",
       "H2-M2")
DotPlot(df, features = fs) + RotatedAxis()
ggsave(file.path(outdir, "dotplot.png"), units = "in", width = 18, height = 7)
ggsave(file.path(outdir, "dotplot.pdf"), units = "in", width = 18, height = 7)
remove(fs)


## all samples plot
Idents(df) <- "orig.ident"
DimPlot(df, reduction = "umap", alpha = 0.2)
ggsave(file.path(outdir, "umap_samples.png"))
ggsave(file.path(outdir, "umap_samples.pdf"))



## make feature plots across all samples
Idents(df) <- "orig.ident"
fs <- c("Mitf", "Tyr", "Dct", "Mlana", "Tyrp1", "Pax3", "Kit", "Sox10", "Snca",
        "Pmel", "S100b", "S100a1")
for (f in fs) {
    FeaturePlot(df, features = f, alpha = 0.2)
    ggsave(file.path(outdir, paste("feature_", f, ".png" , sep = '')))
    ggsave(file.path(outdir, paste("feature_", f, ".pdf" , sep = '')))
}

FeaturePlot(df, features = "SNCA", alpha = 0.2)
ggsave(file.path(outdir, paste("feature_SNCA_artificial.png" , sep = '')))
ggsave(file.path(outdir, paste("feature_SNCA_artificial.pdf" , sep = '')))

sink(file = file.path(outdir, "average_expressions.txt"))
print(AverageExpression(df, features = c("Snca", "SNCA")))
sink()


## find differentially expressed genes between samples
Idents(df) <- "orig.ident"
joined_df <- JoinLayers(df)

sink(file = file.path(outdir, "de_hyper_vs_control.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, "de_hypo_vs_control.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, "de_hypo_vs_hyper.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "a53t-hyper"))
sink()

remove(joined_df)


## find differentially expressed genes between samples limited to a
## single cluster

for (k in as.integer(levels(df$seurat_clusters))) {
   
    Idents(df) <- "seurat_clusters"
    cluster_k <- subset(df, idents = k)
    Idents(cluster_k) <- "orig.ident"
    joined_df <- JoinLayers(cluster_k)

    sink(file = file.path(outdir, paste("de_cluster", k, "hyper_vs_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hypo_vs_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hypo_vs_hyper.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "a53t-hyper"))
    sink()
}


## find differentially expressed genes between samples limited to a
## single cluster (melanocyte)  -- all filters removed

for (k in as.integer(levels(df$seurat_clusters))) {
    
    Idents(df) <- "seurat_clusters"
    cluster_k <- subset(df, idents = k)
    Idents(cluster_k) <- "orig.ident"
    joined_df <- JoinLayers(cluster_k)

    sink(file = file.path(outdir, paste("de_cluster", k, "hyper_vs_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control",
                      logfc.threshold = 0.0, min.pct = 0.0))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hypo_vs_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control",
                      logfc.threshold = 0.0, min.pct = 0.0))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hypo_vs_hyper.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "a53t-hyper",
                      logfc.threshold = 0.0, min.pct = 0.0))
    sink()
}


## cluster k complement   
Idents(df) <- "seurat_clusters"
k <- 3
cluster_k <- subset(df, idents = (0:9)[-4]) # everything except cluster k = 3
Idents(cluster_k) <- "orig.ident"
joined_df <- JoinLayers(cluster_k)

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hyper_to_control.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hypo_to_control.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hyper_to_hypo.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "a53t-hypo"))
sink()


## make violin plots between samples limited to melanocyte cluster (cluster k=3)
genes = c("Mlana", "Sox10", "S100b", "Pmel", "Tyr", "Mitf", "Bace1",
          "Bace2", "Actb", "Gapdh", "Pgk1")


for (gene in genes) {
    VlnPlot(joined_df, features = gene) + NoLegend()
    ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".png", sep = '')))
    ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".pdf", sep = '')))
}

remove(k, cluster_k, joined_df)



########  A FEW QUICK PLOTS


## make expression level plots for MITF and some target genes
Idents(df) <- "seurat_clusters"
k <- 3
cluster_k <- subset(df, idents = k) # cluster 3 *only*
Idents(cluster_k) <- "orig.ident"
joined_df <- JoinLayers(cluster_k)

genes <- c("Trpm1", "Mitf", "Tyr", "Pmel", "Oca2", "Ptgds", "Prnp",
           "Dbi", "Rps24", "Apoc1", "Actb")

for (gene in genes) {
    VlnPlot(joined_df, features = gene) + NoLegend()
    ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".png", sep = '')))
    RidgePlot(object = joined_df, features = gene) + NoLegend()
    ggsave(file.path(outdir, paste("rdgplot_cluster_", k, "_", gene, ".png", sep = '')))
}


VlnPlot(joined_df, features = gene) + NoLegend()
ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".png", sep = '')))

RidgePlot(object = joined_df, features = gene) + NoLegend()
ggsave(file.path(outdir, paste("rdgplot_cluster_", k, "_", gene, ".png", sep = '')))

remove(k, cluster_k, joined_df)


## make heatmap of single cell expression levels in melanocyte cluster
Idents(df) <- "seurat_clusters"
k <- 3
cluster_k <- subset(df, idents = k) # cluster 3 *only*
Idents(cluster_k) <- "orig.ident"
joined_df <- JoinLayers(cluster_k)


################################################################
## additional postprocessing (comment / uncomment as desired) ##
################################################################


## ** NOTE **
##
## The Idents() method sets the identity class for the Seurat object
## (dataframe, df).  This is used by other Seurat functions for
## plotting and analysis, usually through the variables ident.1 and
## ident.2, but sometimes only implicitely.  At risk of being verbose,
## each independent task below has the correct identity class
## explicitely set at the onset.
##
## The available identity classes correspond to the columns of the
## Seurat object, which can be found by, e.g., 'head(df)'.  The values
## taken by any particular identity class can be found by calling
## 'levels(df)' with Idents() set to the desired class.  The number of
## cells/reads in each level of a class can be found using
## table(df$orig.ident), or replacing 'orig.ident' with the desired
## identity class label.
##
## You can also use the 'group.by' argument in many functions rather
## than explicitely changing the identity class.
##


## plot clusters
Idents(df) <- "seurat_clusters"
DimPlot(df, reduction = "umap", alpha = 0.2, label = TRUE)
ggsave(file.path(outdir, "cluster.png"))
ggsave(file.path(outdir, "cluster.pdf"))


## dotplot with specific features
Idents(df) <- "seurat_clusters"
fs = c("Cdh1", "Krt1", "Krt5", "Krt10", "Krt14", "Trp63", "Dsc1", "Ivl",
       "Krt2", "Krt17", "Sostdc1", "Fst", "Cdk1", "Ccnb1", "Ccnb2", "Ccna2",
       "Bub1", "Mitf", "Tyr", "Dct", "Mlana", "Tyrp1", "Pax3", "Kit", "Sox10",
       "Snca", "Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1", "Cd7",
       "Cd28", "Cd44", "Cd48", "Cd69", "Cd82", "Foxp3", "Il2ra", "Ctla2a",
       "Ctla2b", "Ctla4", "Trdc", "Trdv4", "Trgv2", "Lag3", "Cd14", "Cd83",
       "Cd80", "Cd86", "Cd207", "Ccr2", "Il16", "Il18", "Cd274", "Cd200",
       "Cd74", "Prex2", "Plcb4", "Mef2c", "Cadm1", "H2-K1", "H2-Q7", "H2-D1",
       "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-M2")
DotPlot(df, features = fs) + RotatedAxis()
ggsave(file.path(outdir, "dotplot.png"), units = "in", width = 18, height = 7)
ggsave(file.path(outdir, "dotplot.pdf"), units = "in", width = 18, height = 7)
remove(fs)


## all samples plot
Idents(df) <- "orig.ident"
DimPlot(df, reduction = "umap", alpha = 0.2)
ggsave(file.path(outdir, "umap_samples.png"))
ggsave(file.path(outdir, "umap_samples.pdf"))


## differential expression of cluster k to others
Idents(df) <- "seurat_clusters"
joined_df <- JoinLayers(df)
for (k in as.integer(levels(df))) {
    sink(file = file.path(outdir, paste("de_cluster_", k, ".txt", sep = '')))
    print(FindMarkers(joined_df, ident.1 = k))
    sink()
}
remove(joined_df)


## make feature plots across all samples
Idents(df) <- "orig.ident"
fs <- c("Mitf", "Tyr", "Dct", "Mlana", "Tyrp1", "Pax3", "Kit", "Sox10", "Snca",
        "Pmel", "S100b", "S100a1")
for (f in fs) {
    FeaturePlot(df, features = f, alpha = 0.2)
    ggsave(file.path(outdir, paste("feature_", f, ".png" , sep = '')))
    ggsave(file.path(outdir, paste("feature_", f, ".pdf" , sep = '')))
}



## find differentially expressed genes between samples
Idents(df) <- "orig.ident"
joined_df <- JoinLayers(df)

sink(file = file.path(outdir, "de_hyper_to_control.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, "de_hypo_to_control.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, "de_hyper_to_hypo.txt"))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "a53t-hypo"))
sink()

remove(joined_df)


## find differentially expressed genes between samples limited to a
## single cluster

for (k in as.integer(levels(df$seurat_clusters))) {
   
    Idents(df) <- "seurat_clusters"
    cluster_k <- subset(df, idents = k)
    Idents(cluster_k) <- "orig.ident"
    joined_df <- JoinLayers(cluster_k)

    sink(file = file.path(outdir, paste("de_cluster", k, "hyper_to_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hypo_to_control.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
    sink()

    sink(file = file.path(outdir, paste("de_cluster", k, "hyper_to_hypo.txt", sep = '_')))
    print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "a53t-hypo"))
    sink()
}


## cluster k complement   
Idents(df) <- "seurat_clusters"
k <- 3
cluster_k <- subset(df, idents = (0:9)[-4]) # everything except cluster k = 3
Idents(cluster_k) <- "orig.ident"
joined_df <- JoinLayers(cluster_k)

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hyper_to_control.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hypo_to_control.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hypo", ident.2 = "k14-control"))
sink()

sink(file = file.path(outdir, paste("de_comp_cluster", k, "hyper_to_hypo.txt", sep = '_')))
print(FindMarkers(joined_df, ident.1 = "a53t-hyper", ident.2 = "a53t-hypo"))
sink()


## make violin plots between samples limited to melanocyte cluster (cluster k=3)
genes = c("Mlana", "Sox10", "S100b", "Pmel", "Tyr", "Mitf", "Bace1",
          "Bace2", "Actb", "Gapdh", "Pgk1")


for (gene in genes) {
    VlnPlot(joined_df, features = gene) + NoLegend()
    ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".png", sep = '')))
    ggsave(file.path(outdir, paste("vlnplot_cluster_", k, "_", gene, ".pdf", sep = '')))
}

remove(k, cluster_k, joined_df)





## average expression of biomarkers scatter plots

Idents(df) <- "seurat_clusters"
sdf <- JoinLayers(subset(df, idents = 3))  # subset melanocyte cluster
Idents(sdf) <- "orig.ident"

melanocyte.markers <- c("Mitf", "Mlana", "Tyr", "Dct", "Tyrp1", "Pax3", "Kit", "Sox10")
keratinocyte.markers <- c("Cdh1", "Krt1", "Krt5", "Krt14", "Trp63")

mean.markers <- function(df, markers) {
    acc <- FetchData(df, vars = markers[[1]])
    for (marker in tail(markers, -1)) {
        acc <- acc + FetchData(df, vars = marker)
    }
    acc <- acc / length(markers)
}

sdf <- AddMetaData(sdf,
                   metadata = mean.markers(sdf, melanocyte.markers),
                   col.name = "melanocyte.markers")

sdf <- AddMetaData(sdf,
                   metadata = mean.markers(sdf, keratinocyte.markers),
                   col.name = "keratinocyte.markers")

FeatureScatter(sdf, feature1 = "melanocyte.markers", feature2 = "keratinocyte.markers")
ggsave(file.path(outdir, "biomarker-scatter-plot.png")

