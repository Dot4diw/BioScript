library(Seurat)
library(monocle)
library(ggplot2)
library(convert2anndata)
library(anndata)



obj <- readRDS("./seurat.rds")
cds <- readRDS("./seurat.rds_monocle2.rds")

pseudo_coord <- as.data.frame(t(cds@reducedDimS))
colnames(pseudo_coord) <- c("Component_1","Component_2","Component_3")
head(pseudo_coord)
#obj@reductions$umap <- CreateDimReducObject(embedding = as.matrix(pseudo_coord), key = "Component_")
obj@reductions$monocle2 <- CreateDimReducObject(embedding = as.matrix(pseudo_coord), key = "Component_")

DimPlot(object = obj, reduction = "monocle2")

obj$Size_Factor <- pData(cds)$Size_Factor
obj$num_genes_expressed <- pData(cds)$num_genes_expressed
obj$Pseudotime <- pData(cds)$Pseudotime
obj$State <- pData(cds)$State

# Convert to SingleCellExperiment if necessary
sce <- convert_seurat_to_sce(obj)
# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)
# Save the AnnData object
write_h5ad(ad, "scanpy_merger_monocle2.h5ad")
