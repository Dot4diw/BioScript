library(Seurat)
library(DoubletFinder)
library(hdf5r)
library(dplyr)
library(magrittr)
library(patchwork)
library(SoupX)
library(DropletUtils)


# sample info content
samplename	sample_type	group	sampling_time	raw_matrixpath	filtered_matrixpath	velocyto_loompath
S1	CTRL	CTRL	0m	/data/S1/cellranger/outs/raw_feature_bc_matrix.h5	/data/S1/cellranger/outs/filtered_feature_bc_matrix.h5	/data/S1/cellranger/velocyto/cellranger_count_tdTomato.loom
S2	CTRL	CTRL	2m	/data/S1/cellranger/outs/raw_feature_bc_matrix.h5	/data/S1/cellranger/outs/filtered_feature_bc_matrix.h5	/data/S1/cellranger/velocyto/cellranger_count_tdTomato.loom

sampleinfo <- read.table("./sample_info.txt", header=T, sep = "\t")

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

    srat <- ScaleData(srat, features = VariableFeatures(object = srat))
    srat <- RunPCA(srat, features = VariableFeatures(object = srat))

    # Find significant PCs
    stdev <- srat[["pca"]]@stdev
    percent_stdev <- (stdev/sum(stdev)) * 100
    cumulative <- cumsum(percent_stdev)
    co1 <- which(cumulative > 90 & percent_stdev < 5)[1] 
    co2 <- sort(which((percent_stdev[1:length(percent_stdev) - 1] - 
                         percent_stdev[2:length(percent_stdev)]) > 0.1), 
                decreasing = T)[1] + 1
    min_pc <- min(co1, co2)
    
    pcs.num <- 1:min_pc
    print(paste0("The significant PCs used in the next analysis : ", min_pc))
    srat    <- FindNeighbors(srat, dims = pcs.num, verbose = FALSE)
    srat    <- FindClusters(srat, resolution = res, verbose = FALSE)
    srat <- RunUMAP(srat, dims = pcs.num,verbose = FALSE)
    
    meta <- srat@meta.data
    umap <- srat@reductions$umap@cell.embeddings

    soup.channel <- SoupChannel(raw.matrix, filter.matrix)
    soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
    soup.channel  <- setDR(soup.channel, umap)
    soup.channel  <- autoEstCont(soup.channel)
    saveRDS(soup.channel, paste0(samplename,"_SoupX_obj.rds"))
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

    DropletUtils:::write10xCounts(paste0("./SoupXCounts/",samplename), adj.matrix)

    obj <- CreateSeuratObject(counts = adj.matrix, project = samplename, min.cells = 3, min.features = 500)
    return(obj)
}


# RunDoubletFinder function
RunDoubletFinder <- function(obj, preset.doublet.rate=NULL) {
    # @obj : the input Seurat object is filtered(Normalizae & Scale & RunPCA & FindNeighbor & RunUMAP & FindClusters .
    # DoubletFinder
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    recover.cell_num <- dim(obj)[2]
    
    # Get the multiplet rate if not provided
    if(is.null(preset.doublet.rate)){
        print('preset.doublet.rate not provided.......')
        print('estimating multiplet rate from cells in dataset......')
    
        # 10X multiplet rates table
        #https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
        multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.024, 0.032, 0.040, 0.048, 0.056, 0.064, 0.072, 0.080),
                                      'Loaded_cells' = c(825, 1650, 3300, 4950, 6600, 8250, 9900, 11550, 13200, 14850, 16500),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
        print("======================Multiplet Rates of 10x Single Cell======================")
        print(multiplet_rates_10x)
        print("==============================================================================")
        
        preset.doublet.rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < recover.cell_num) %>% 
            dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
            dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
        
        print( paste0("cell number: ", recover.cell_num, " ==> preset doublet rate is : ", preset.doublet.rate) )
    }
    
    sweep.res <- paramSweep(obj, PCs = pcs.num, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    optimalpk <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]))

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- obj$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(preset.doublet.rate*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    obj <- doubletFinder(obj, PCs = pcs.num, pN = 0.25, pK = optimalpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    pANN_value <- paste0("pANN_","0.25_", optimalpk,"_", nExp_poi)
    obj <- doubletFinder(obj, PCs = pcs.num, pN = 0.25, pK = optimalpk, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = FALSE)
    
    doubletfinder_col <- paste0("DF.classifications_0.25_", optimalpk, "_", nExp_poi.adj)
    colnames(obj@meta.data)[grep(doubletfinder_col, colnames(obj@meta.data))] <- "doublet_finder"
  
    return(obj)
}


res = 1

obj_list <- list()

for ( samplename in sampleinfo$samplename)
{
    # Read data and create seurat object
    filter_matrix <- sampleinfo[sampleinfo$samplename == samplename,]$filtered_matrixpath
    raw_matrix <- sampleinfo[sampleinfo$samplename == samplename,]$raw_matrixpath
    
    sampletype <- sampleinfo[sampleinfo$samplename == samplename,]$sample_type
    group <- sampleinfo[sampleinfo$samplename == samplename,]$group
    sampling_time <- sampleinfo[sampleinfo$samplename == samplename,]$sampling_time

    obj_list[[samplename]] <- RunSoupX(filter_matrix, raw_matrix, samplename)

    # Add sampl info
    obj_list[[samplename]][["percent.mt"]] <- PercentageFeatureSet(obj_list[[samplename]], pattern = "^mt-")
    obj_list[[samplename]][["percent.ribo"]] <- PercentageFeatureSet(obj_list[[samplename]], pattern = "^Rp[sl]")
    obj_list[[samplename]][["percent.hb"]] <- PercentageFeatureSet(obj_list[[samplename]], pattern = "^Hb[^(p)]")


    obj_list[[samplename]][["sample_name"]] <- samplename
    obj_list[[samplename]][["sample_type"]] <- sampletype
    obj_list[[samplename]][["group"]] <- group
    obj_list[[samplename]][["sampling_time"]] <- sampling_time

    # Filter cells
    obj_list[[samplename]] <- subset(obj_list[[samplename]], nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)

    # Normalize and Scale
    obj_list[[samplename]] <- NormalizeData(obj_list[[samplename]], normalization.method = "LogNormalize", scale.factor = 10000)
    obj_list[[samplename]] <- FindVariableFeatures(obj_list[[samplename]], selection.method = "vst", nfeatures = 2000)
    
    #all.genes <- rownames(obj_list[[samplename]])
    obj_list[[samplename]] <- ScaleData(obj_list[[samplename]], features = VariableFeatures(object = obj_list[[samplename]]))
    
    # Run PCA 
    obj_list[[samplename]] <- RunPCA(obj_list[[samplename]], features = VariableFeatures(object = obj_list[[samplename]]), verbose = FALSE)

    # Find significant PCs
    stdev <- obj_list[[samplename]][["pca"]]@stdev
    percent_stdev <- (stdev/sum(stdev)) * 100
    cumulative <- cumsum(percent_stdev)
    co1 <- which(cumulative > 90 & percent_stdev < 5)[1] 
    co2 <- sort(which((percent_stdev[1:length(percent_stdev) - 1] - 
                         percent_stdev[2:length(percent_stdev)]) > 0.1), 
                decreasing = T)[1] + 1
    min_pc <- min(co1, co2)
    
    pcs.num <- 1:min_pc
    print(paste0("The significant PCs used in the next analysis : ", min_pc))

    # Findcluster and run umap
    obj_list[[samplename]] <- FindNeighbors(obj_list[[samplename]], dims=pcs.num, verbose = FALSE)
    obj_list[[samplename]] <- FindClusters(obj_list[[samplename]], resolution=res, verbose = FALSE)
    obj_list[[samplename]] <- RunUMAP(obj_list[[samplename]], dims=pcs.num, verbose = FALSE)

    # Run DoubletFinder
    obj_list[[samplename]] <- RunDoubletFinder(obj_list[[samplename]], preset.doublet.rate=NULL)
    saveRDS(obj_list[[samplename]], paste0("./soupx_doubletfinder_result/",samplename, "_after_soupx_and_doubletfinder.rds"))
}

