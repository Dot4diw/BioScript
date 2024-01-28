#!/usr/bin/bash
set -x
Rscript prepatring_cellphonedb_inputfile.r
# The Content of `prepatring_cellphonedb_inputfile.r`
###  #!/usr/bin/R
###  
###  library(Seurat)
###  
###  scRNA.obj <- readRDS("../../celltype_annotation_by_manual3.rds")
###  
###  count_raw <- scRNA.obj@assays$RNA@data
###  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
###  
###  #write.table(as.matrix(scRNA.obj@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
###  write.table(count_norm, 'cellphonedb_count.txt', sep='\t', quote=F)
###  
###  
###  meta_data <- cbind(rownames(scRNA.obj@meta.data), scRNA.obj@meta.data[,'celltype_manual', drop=F])
###  meta_data <- as.matrix(meta_data)
###  meta_data[is.na(meta_data)] = "Unkown" #  The celltype cant't NA
###  
###  write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

# cellphonedb_meta.txt and cellphonedb_count.txt are output files of prepatring_cellphonedb_inputfile.r
# Install the cellphonedb by pip or conda.
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt \
    --threads=8 \
    --counts-data=gene_name \
    --output-format csv

cellphonedb plot dot_plot \
    --means-path ./out/means.csv \
    --pvalues-path ./out/pvalues.csv

cellphonedb plot heatmap_plot cellphonedb_meta.txt \
    --means-path ./out/means.csv \
    --pvalues-path ./out/pvalues.csv