library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SoupX)
library(decontX)
library(DropletUtils)


# RunSoupX function
RunSoupX <- function(filter_matrix, raw_matrix, samplename){
    # SoupX
    # @filter_matrix: filtered_feature_bc_matrix.h5 of cellragner output
    # @raw_matrix: raw_feature_bc_matrix.h5 of cellranger output
    # @samplename: sample names
    filter.matrix <- Read10X_h5(filter_matrix)
    raw.matrix <- Read10X_h5(raw_matrix)

    srat <- CreateSeuratObject(counts = filter.matrix, project = samplename)

    srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(srat)
    srat <- ScaleData(srat, features = all.genes)

    srat <- RunPCA(srat, features = VariableFeatures(object = srat))
    srat    <- FindNeighbors(srat, dims = pcs.num, verbose = F)
    srat    <- FindClusters(srat, resolution = res, verbose = T)
    srat <- RunUMAP(srat, dims = pcs.num)
    
    meta <- srat@meta.data
    umap <- srat@reductions$umap@cell.embeddings

    soup.channel <- SoupChannel(raw.matrix, filter.matrix)
    soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
    soup.channel  <- setDR(soup.channel, umap)
    soup.channel  <- autoEstCont(soup.channel)
    saveRDS(soup.channel, paste0(samplename,"_SoupX_obj.rds"))
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
    DropletUtils:::write10xCounts("./SoupXCounts", adj.matrix)

    obj <- CreateSeuratObject(counts = adj.matrix, project = samplename, min.cells = 3, min.features = 200)
    return(obj)
    #sc = load10X(cellranger_outs_path)
    #sc = autoEstCont(sc)
    #out = adjustCounts(sc,roundToInt=TRUE)
    #DropletUtils:::write10xCounts("./SoupXCounts", out)
}

# RunDecontX function
RunDecontX <- function(filter_path, raw_path, samplename, group){
    # @filter_path : the "filtered_feature_bc_matrix" directory of cellranger outs.
    # @raw_path: the "raw_feature_bc_matrix" directory of cellranger outs.
    # @samplename: the input sample names. aways is file name.
    # @group: the group information of samples.
    filter_counts <- Read10X(filter_path)
    sce.filter <- SingleCellExperiment(list(counts = filter_counts))

    raw_counts <- Read10X(raw_path)
    sce.raw <- SingleCellExperiment(list(counts = filter_counts))

    decontX_results <- decontX(sce.filter, background = sce.raw)

    saveRDS(decontX_results, paste0(samplename,"_decontX_results.rds"))

    prop <- plotDecontXContamination(decontX_results)
    ggsave(paste0(samplename, "_decontX_contamination.pdf"), plot = prop,  width = 8, height=6)

    obj <- CreateSeuratObject(counts = round(decontXcounts(decontX_results)), project = samplename, min.cells = 0, min.features = 0)
    obj$group <- group
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

    obj$contamination <- sce@metadata$decontX$contamination
    obj$decontX_label  <- ifelse(obj$contamination < 0.2, "Normal", "Contaminated")

    return(obj)
}

# RunDoubletFinder function
RunDoubletFinder <- function(obj) {
    # @obj : the input Seurat object.
    # DoubletFinder
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res <- paramSweep(obj, PCs = pcs.num, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    optimalpk <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]))

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- obj$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    obj <- doubletFinder(obj, PCs = pcs.num, pN = 0.25, pK = optimalpk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(obj@meta.data)[grep("DF.classifications", colnames(obj@meta.data))] <- "DF.classifications"
    #pANN_value <- paste0("pANN_","0.25_",optimalpk,"_",nExp_poi)
    #obj <- doubletFinder(obj, PCs = pcs.num, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = FALSE)

    return(obj)
}


sample_info <- read.table("../sampleinfo.txt", sep = '\t', header = T)
# samplename matrixpath group

samplename <- basename(getwd())
group <- sample_info[sample_info$samplename == samplename,]$group

#samplename = sample_info$samplename
#matrixpath = sample_info$matrixpath
samplename <- basename(getwd())
group <- sample_info[sample_info$samplename == samplename,]$group


filter_matrix = "./outs/filtered_feature_bc_matrix.h5"
raw_matrix = "./outs/raw_feature_bc_matrix.h5"

#Run SoupX
obj <- RunSoupX(filter_matrix, raw_matrix, samplename)
obj$group <- group


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# Normalization and scale the raw data
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
# Run pca
obj <- RunPCA(obj,features = VariableFeatures(object = obj))
ElbowPlot(obj, ndims = 50)

# Finder cluter and run umap
obj <- FindNeighbors(obj, dims=pcs.num)
obj <- FindClusters(obj, resolution=res)
obj <- RunUMAP(obj, dims=pcs.num)

obj <- RunDoubletFinder(obj)
saveRDS(obj, paste0(samplename, "soupx_doubletfiner.out.rds"))


