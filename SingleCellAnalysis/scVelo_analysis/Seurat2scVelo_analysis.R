library(Seurat)
library(qs)

obj <- qread("./seurat_analysis_obj.qs")
#obj <- JoinLayers(obj)

# Joinlayers
obj_join = obj_td

# save metadata table:
obj_join$barcodes <- colnames(obj_join)
obj_join$UMAP_1 <- obj_join@reductions$umap@cell.embeddings[,1]
obj_join$UMAP_2 <- obj_join@reductions$umap@cell.embeddings[,2]

write.csv(obj_join@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(obj_join, assay='RNA', layer = 'counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(obj_join@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
write.csv(obj_join@reductions$harmony@cell.embeddings, file='harmony.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
