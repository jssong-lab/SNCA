
############################################################################
## Author  : Porter B. Howland                                            ##
## Email   : pbh2@illinois.edu                                            ##
## Date    : August 31, 2024                                              ##
## Updated : March 31, 2025
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

##
## useful links:
##     https://satijalab.org/seurat/articles/integration_introduction.html
##

library(dplyr)
library(Seurat)
library(scales)
library(patchwork)
library(stringdist)
library(metap)


## file for serialized results of umap/clustering and output directory
## name (to keep workspace tidy)
rds_filename <- "integrated_sample_0_2.rds"
outdir <- "merged_keratinocyte_analysis"

## parameters for downstream analysis:
num_pca_dims <- 20
clustering_resolution <- 0.2
subclustering_resolution <- 0.1


## check to see if the directory already exists and warn user...
if (dir.exists(outdir)) {
    print(paste("The output directory '", outdir,
                "' already exists!  Exiting...", sep = ''))
    stop(1)
} else {
    dir.create(outdir)
}


## check if a prepared RDS file exists already
if (!file.exists(rds_filename)) {

    ## Again, the directory structure and file locations are specific to
    ## the cluster we are using; this section should be adjusted to read
    ## in the output of 10x pipeline.
    
    ## the order of projects, (directory) identifiers, and prefixes should
    ## all coincide -- maybe change the name of this vector?
    ids <- list("IGH01-lesional", "IGH01-perilesional",
                "IGH02-lesional", "IGH02-perilesional",
                "IGH03-lesional", "IGH03-perilesional",
                "IGH04-lesional", "IGH04-perilesional")
    n <- length(ids)

    ## thresholds for count selection - see documentation for
    ## 'CreateSeuratObject'
    min_cells <- 3
    min_features <- 200

    ## construct filepaths (these will be different in each experiment)
    ## -- location at:
    ## -- /home/groups/song/songlab2/shared/NicholasTheodosakisData/scRNA-seq6/Analysis
    filepaths = list()
    for (i in 1:n) {
        filepaths[[i]] <- file.path(
            "/", "home", "groups", "song", "songlab2", "shared", "NicholasTheodosakisData",
            "scRNA-seq6", "Analysis", paste("count", ids[[i]], sep = '_'),
            "outs", "filtered_feature_bc_matrix"
        )
    }    
    
    ## read in (sparse) counts matrices and create seurat objects
    sobs = list()    
    for (i in 1:n) {
        sobs[[i]] <- CreateSeuratObject(
            counts = Read10X(data.dir = filepaths[[i]]),
            project = ids[[i]],
            min.cells = min_cells,
            min.features = min_features
        )
    }

    ## add patient and sample (type) labels
    label <- 1
    for (i in seq(1, 8, by=2)) {
        sobs[[i]]$patient <- as.character(label)
        sobs[[i+1]]$patient <- as.character(label)
        sobs[[i]]$sample <- "lesional"
        sobs[[i+1]]$sample <- "perilesional"
        label <- label + 1
    }
    remove(label)

    ## merge datasets together
    df <- merge(
        x = sobs[[1]],
        y = sobs[-1],
        add.cell.ids = ids,
        project = "all_samples_merged"
    )

    ## run some basic quality control to check data
    for (i in 1:n) {
        ## compute mitochondria content for each sample
        sobs[[i]][["percent.mt"]] <- PercentageFeatureSet(sobs[[i]], pattern = "^MT-")
        
        ## violin plots
        VlnPlot(sobs[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        ggsave(file.path(outdir, paste("qc_vln_plots_", ids[[i]], ".png", sep='')))
        ggsave(file.path(outdir, paste("qc_vln_plots_", ids[[i]], ".pdf", sep='')))

        ## scatter plots for mitochondria counts
        FeatureScatter(sobs[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
        ggsave(file.path(outdir, paste("qc_percent_MT_plots_", ids[[i]], ".png", sep='')))
        ggsave(file.path(outdir, paste("qc_percent_MT_plots_", ids[[i]], ".pdf", sep='')))

        ## scatter plots for features
        FeatureScatter(sobs[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        ggsave(file.path(outdir, paste("qc_nFeature_plots_", ids[[i]], ".png", sep='')))
        ggsave(file.path(outdir, paste("qc_nFeature_plots_", ids[[i]], ".pdf", sep='')))

    }

    ## perform normalization and data scaling
    df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
    df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
    df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
    df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
    df <- ScaleData(df, features = rownames(df))

    ## split based of off patient
    df <- JoinLayers(df)
    df[["RNA"]] <- split(df[["RNA"]], f = df$patient)
    df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)  # Why do I need this and the line above?  Check this again later...

    ## PCA and clustering, unintegrated
    df <- RunPCA(df, reduction.name = "patient_unintegrated_pca")
    df <- FindNeighbors(df, reduction = "patient_unintegrated_pca", dims = 1:num_pca_dims)
    df <- FindClusters(df,
                       resolution = clustering_resolution,
                       cluster.name = "patient_unintegrated_clusters")
    df <- RunUMAP(df, dims = 1:num_pca_dims,
                  reduction = "patient_unintegrated_pca",
                  reduction.name = "patient_unintegrated_umap")
    
    ## integrate layers, cluster, and UMAP
    df <- IntegrateLayers(object = df, method = CCAIntegration,
                          orig.reduction = "patient_unintegrated_pca",
                          new.reduction = "patient_integrated.cca",
                          verbose = FALSE)
    df <- JoinLayers(df)
    df <- FindNeighbors(df,
                        reduction = "patient_integrated.cca",
                        dims = 1:num_pca_dims,
                        graph.name = "patient_integrated_snn")
    df <- FindClusters(df,
                       resolution = clustering_resolution,
                       graph.name = "patient_integrated_snn",
                       cluster.name = "patient_integrated_clusters")
    df <- RunUMAP(df, dims = 1:num_pca_dims,
                  reduction = "patient_integrated.cca",
                  reduction.name = "patient_integrated_umap")
    
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


##
## find number of cells in each idents
##
sink(file = file.path(outdir, "idents_tables.txt"))
Idents(df) <- "patient_integrated_clusters"
print(table(df$patient_integrated_clusters))
print(table(df$sample))
print(table(df$patient))
for (k in levels(df$patient_integrated_clusters)) {
    print(paste("## cluster", k))
    print(table(subset(df, idents = k)$patient))
    print(table(subset(df, idents = k)$sample))
}
sink()


## re-order the levels for the clustering data
df$patient_integrated_clusters <- factor(df$patient_integrated_clusters,
                                         levels(df$patient_integrated_clusters)[order(as.numeric(levels(df$patient_integrated_clusters)))])


## remove clusters with <= 5 cells as these are likely outliers
Idents(df) <- "patient_integrated_clusters"
cells_to_remove <- c()
cluster_table <- table(df$patient_integrated_clusters)
for (k in names(cluster_table)) {
    if (cluster_table[k] <= 5) {
        cells_to_remove <- append(cells_to_remove,
                                  Cells(df)[which(Idents(df) == k)])
    }
}

df <- subset(df, cells = cells_to_remove, invert = TRUE)
df$patient_integrated_clusters <- droplevels(df$patient_integrated_clusters)
df$patient_unintegrated_clusters <- droplevels(df$patient_unintegrated_clusters)


##
## now just check DE genes between clusters, with lesional and perilesional combined,
## all patients combined (for cell-type identification)
##

Idents(df) <- "patient_integrated_clusters"
for (cluster in levels(df)) {
    sink(file = file.path(outdir, paste0("de_cluster_", cluster, "_to_all_else.txt")))
    print(FindMarkers(df, ident.1 = cluster))
    sink()
}

## check spinous keratinocytes and other mystery cell types...
sink(file = file.path(outdir, paste0("de_cluster_12_to_1.txt")))
print(FindMarkers(df, ident.1 = 12, ident.2 = 1))
sink()


## elbow plots for PCA dims determination
ElbowPlot(df, ndims = 30, reduction = "patient_unintegrated_pca")
ggsave(file.path(outdir, "elbow_plot_unintegrated.png"))
ggsave(file.path(outdir, "elbow_plot_unintegrated.pdf"))

## plot clusters by patient (1,2,3,4) and sample (lesional/perilesional)
DimPlot(df, reduction = "patient_integrated_umap", alpha = 0.2, group.by = "patient", label = TRUE)
ggsave(file.path(outdir, "dimplot_patient.png"))
ggsave(file.path(outdir, "dimplot_patient.pdf"))

DimPlot(df, reduction = "patient_integrated_umap", alpha = 0.2, group.by = "sample", label = TRUE)
ggsave(file.path(outdir, "dimplot_sample.png"))
ggsave(file.path(outdir, "dimplot_sample.pdf"))

Idents(df) <- "seurat_clusters"
DimPlot(df, reduction = "patient_integrated_umap", alpha = 0.2, label = TRUE)
ggsave(file.path(outdir, "dimplot_seurat_clusters.png"))
ggsave(file.path(outdir, "dimplot_seurat_clusters.pdf"))

Idents(df) <- "orig.ident"
DimPlot(df, reduction = "patient_integrated_umap", alpha = 0.2)
ggsave(file.path(outdir, "dimplot_orig_ident.png"))
ggsave(file.path(outdir, "dimplot_orig_ident.pdf"))

features to look for and plot
fs <- c("CDH1", "KRT1", "KRT5", "KRT10", "KRT14", "TP63", "DSC1", "IVL",
       "KRT2", "KRT17", "SOSTDC1", "FST", "CDK1", "CCNB1", "CCNB2", "CCNA2",
       "BUB1", "MITF", "TYR", "DCT", "MLANA", "TYRP1", "PAX3", "KIT", "SOX10",
       "SNCA", "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B1", "CD7",
       "CD28", "CD44", "CD48", "CD69", "CD82", "FOXP3", "IL2RA", "CTLA2A",
       "CTLA2B", "CTLA4", "TRDC", "TRDV4", "TRGV2", "LAG3", "CD14", "CD83",
       "CD80", "CD86", "CD207", "CCR2", "IL16", "IL18", "CD274", "CD200",
       "CD74", "PREX2", "PLCB4", "MEF2C", "CADM1", "H2-K1", "H2-Q7", "H2-D1",
       "H2-AA", "H2-AB1", "H2-EB1", "H2-M2")

all_features <- Features(df)
mask <- fs %in% all_features

## make a dotplot of all found features
Idents(df) <- "patient_integrated_clusters"
DotPlot(df, features = fs[mask], split.by = "patient", cols = hue_pal()(4)) +
    RotatedAxis() + xlab("Gene") + ylab("Cluster_Patient") +
    ggtitle("Patient integrated clusters: all patients")
ggsave(file.path(outdir, "dotplot_all_patients.png"), units = "in", width = 18, height = 12)
ggsave(file.path(outdir, "dotplot_all_patients.pdf"), units = "in", width = 18, height = 12)

## dotplot separately for each patient
Idents(df) <- "patient"
for (k in levels(Idents(df))) {
    DotPlot(subset(df, subset = patient == k),
            features = fs[mask],
            group.by = "patient_integrated_clusters") + RotatedAxis() +
        xlab("Gene") + ylab("Cluster") +
        ggtitle(paste("Patient integrated clusters: patient", k))
    ggsave(file.path(outdir, paste0("dotplot_patient_", k, ".png")), units = "in", width = 18, height = 7)
    ggsave(file.path(outdir, paste0("dotplot_patient_", k, ".pdf")), units = "in", width = 18, height = 7)
}

## make feature plots across all samples with only found features
Idents(df) <- "orig.ident"
for (f in fs[mask]) {
    FeaturePlot(df, features = f, alpha = 0.2, reduction = "patient_integrated_umap")
    ggsave(file.path(outdir, paste("feature_", f, ".png" , sep = '')))
    ggsave(file.path(outdir, paste("feature_", f, ".pdf" , sep = '')))
}

## look for similar features (typos?), note in log file, and make
## feature plots using the closest matches
## NOTE: this requires the stringdist package
closest_matches <- all_features[amatch(fs[!mask], all_features, maxDist = Inf)]
print(data.frame(fs[!mask], closest_matches))

for (f in closest_matches) {
    FeaturePlot(df, features = f, alpha = 0.2, reduction = "patient_integrated_umap")
    ggsave(file.path(outdir, paste("feature_", f, ".png" , sep = '')))
    ggsave(file.path(outdir, paste("feature_", f, ".pdf" , sep = '')))
}

remove(fs, all_features, mask, closest_matches)


## ## run differential expression analysis integrated by patients for
## ## each cluster using default meta-analysis methods for combination
## ## of p-values; p-values will be combined manually later so this
## ## portion of output is essentially ignored

## Idents(df) <- "patient_integrated_clusters"
## for (k in as.integer(levels(df$patient_integrated_clusters))) {    

##     cluster_k <- JoinLayers(subset(df, idents = k))
##     Idents(cluster_k) <- "sample"

##     sink(file = file.path(outdir, paste("de_cluster", k, "fisher_integrated_lesional_vs_perilesional.txt", sep = '_')))
##     print(FindConservedMarkers(cluster_k,
##                                ident.1 = "lesional",
##                                ident.2 = "perilesional",
##                                grouping.var = "patient"))
##     sink()

## }


##
## subcluster within the melanocyte cluster
##

k <- 2  # melanocyte cluster
df <- FindSubCluster(object = df,
                     cluster = k,
                     graph.name = "patient_integrated_snn",
                     subcluster.name = "melanocyte_subclustering",
                     resolution = subclustering_resolution)

## manually re-order the levels - this is only to make the legend of
## figure below look nicer
df@meta.data$melanocyte_subclustering <- factor(df@meta.data$melanocyte_subclustering,
                                                levels = c("0", "1", "2_0", "2_1", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

## adjust the color map so 2_0 and 2_1 contrasting colors
colors <- hue_pal()(16)
colors[4] <- colors[9]
colors <- colors[-9]

## plot melanocyte subclustering
DimPlot(df,
        reduction = "patient_integrated_umap",
        alpha = 0.3,
        group.by = "melanocyte_subclustering",
        cols = colors,
        label = TRUE)
ggsave(file.path(outdir, "dimplot_melanocyte_subcluster.png"))
ggsave(file.path(outdir, "dimplot_melanocyte_subcluster.pdf"))

## updated dotplot features
fs <- c("CDH1", "KRT1", "KRT5", "KRT10", "KRT14", "TP63", "DSC1",
        "IVL", "KRT2", "MYO5B", "KRT17", "SOSTDC1", "FST", "CDK1",
        "CCNB1", "CCNB2", "CCNA2", "BUB1", "MITF", "TYR", "DCT",
        "MLANA", "TYRP1", "PAX3", "KIT", "SOX10", "SNCA", "PTPRC",
        "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B1", "CD7", "CD28",
        "CD44", "CD48", "CD69", "CD82", "FOXP3", "IL2RA", "CTLA2A",
        "CTLA2B", "CTLA4", "TRDC", "TRDV4", "TRGV2", "CD14", "CD83",
        "CD86", "CD207", "SELE", "VWF", "ACE", "CD34", "COL1A1",
        "COL1A2", "FBLN1", "IL16", "IL18", "CD200", "CD74", "PREX2",
        "PLCB4", "MEF2C", "CADM1", "H2-K1", "H2-Q7", "H2-D1", "H2-AA",
        "H2-AB1", "H2-EB1", "H2-M2", "LYVE1", "PROX1", "PDPN")

all_features <- Features(df)
mask <- fs %in% all_features

## dotplot with melanocyte subclusters
Idents(df) <- "melanocyte_subclustering"
DotPlot(df, features = fs[mask]) +  RotatedAxis() + xlab("Gene") + ylab("Cluster")
ggsave(file.path(outdir, "dotplot_melanocyte_subcluster_markers.png"), units = "in", width = 18, height = 7)
ggsave(file.path(outdir, "dotplot_melanocyte_subcluster_markers.pdf"), units = "in", width = 18, height = 7)


##
## now find the proportions of cells in each cluster for each patient and sample
##

## cast sample and patient metadata columns to factors
df@meta.data$patient <- as.factor(df@meta.data$patient)
df@meta.data$sample <- as.factor(df@meta.data$sample)

df.counts <- data.frame(cluster = c(), patient = c(), sample = c(), count = c())
for (cluster in levels(df@meta.data$melanocyte_subclustering)) {
    for (p in levels(df@meta.data$patient)) {
        for (s in levels(df@meta.data$sample)) {
            count <- df@meta.data %>%
                filter(melanocyte_subclustering == cluster) %>%
                filter(patient == p) %>%
                filter(sample == s) %>%
                nrow
            df.counts <- dplyr::bind_rows(
                                    df.counts,
                                    data.frame(cluster = c(cluster),
                                               patient = c(p),
                                               sample = c(s),
                                               count = c(count))
                                )
        }
    }
}

## Save an RDS file with only the melanocyte subcluster;
## since there is a huge amount of data for this experiment,
## R-studio has trouble reading in a full processed RDS file
## for all clusters.  Saving only the subcluster of interest
## (melanocyte) mitigates this issue.
##
Idents(df) <- "melanocyte_subclustering"
saveRDS(subset(df, idents = "2_0"), file = "strict_melanocyte_cluster_only_IGH.rds")



## print number of cells in subclustering
sink(file = file.path(outdir, "idents_subclustering.txt"))
print(table(df$melanocyte_subclustering))
sink()


## Further analysis for subcluster.  The scatter plot of mean
## melanocyte and keratinocyte markers is generated here.

## markers to use for feature scatter plots and GMM subclustering
melanocyte.markers <- c("MITF", "TYR", "DCT", "MLANA", "TYRP1", "PAX3", "KIT", "SOX10", "SNCA")
keratinocyte.markers <- c("CDH1", "KRT1", "KRT2", "KRT5", "KRT10", "KRT14", "KRT17", "TP63")
all.markers <- append(melanocyte.markers, keratinocyte.markers)


## DotPlot(df, features = all.markers, idents = sub_levels(df$melanocyte_subclustering, k)) +
DotPlot(df, features = fs[mask]) +
    RotatedAxis() + xlab("Gene") + ylab("Cluster") + ggtitle("Subclustering of melanocytes")
ggsave(file.path(outdir, "dotplot_melanocyte_subcluster.png"), units = "in", width = 18, height = 7)
ggsave(file.path(outdir, "dotplot_melanocyte_subcluster.pdf"), units = "in", width = 18, height = 7)


## plot melanocyte subclustering
DimPlot(df, reduction = "patient_integrated_umap", alpha = 0.3, label = TRUE)
ggsave(file.path(outdir, "dimplot_melanocyte_subcluster.png"))
ggsave(file.path(outdir, "dimplot_melanocyte_subcluster.pdf"))


Idents(df) <- "patient_integrated_clusters"
sdf <- JoinLayers(subset(df, idents = k))  # subset melanocyte cluster
Idents(sdf) <- "patient"


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

df@meta.data[rownames(sdf@meta.data), "melanocyte.markers"] <- sdf@meta.data[, "melanocyte.markers"]
df@meta.data[rownames(sdf@meta.data), "keratinocyte.markers"] <- sdf@meta.data[, "keratinocyte.markers"]


## feature scatters of melanocytes split by patient (without groups)
FeatureScatter(sdf, feature1 = "melanocyte.markers", feature2 = "keratinocyte.markers", split.by = "patient")
ggsave(file.path(outdir, "biomarkers_patient_split.png"), units = "in", width = 17, height = 4)
ggsave(file.path(outdir, "biomarkers_patient_split.pdf"), units = "in", width = 17, height = 4)


## biomarker scatter plots using Seurat's FindSubCluster()
FeatureScatter(sdf, feature1 = "melanocyte.markers", feature2 = "keratinocyte.markers", split.by = "patient", group.by = "melanocyte_subclustering")
ggsave(file.path(outdir, "biomarkers_patient_split_subcluster_group.png"), units = "in", width = 17, height = 4)
ggsave(file.path(outdir, "biomarkers_patient_split_subcluster_group.pdf"), units = "in", width = 17, height = 4)

FeatureScatter(sdf, feature1 = "melanocyte.markers", feature2 = "keratinocyte.markers", group.by = "melanocyte_subclustering")
ggsave(file.path(outdir, "biomarkers_all_subcluster_group.png"))
ggsave(file.path(outdir, "biomarkers_all_subcluster_group.pdf"))


##
## now find conserved markers using the melanocyte marker dominant portion of melanocyte cells
##

melanocyte.strict <- JoinLayers(subset(df, subset = melanocyte_subclustering == "2_0" ))
Idents(melanocyte.strict) <- "sample"

sink(file = file.path(outdir, "melanocyte_subclustering_integrated_lesional_vs_perilesional.txt"))
print(FindConservedMarkers(melanocyte.strict,
                           ident.1 = "lesional",
                           ident.2 = "perilesional",
                           grouping.var = "patient",
                           meta.method = metap::sumlog))
sink()

