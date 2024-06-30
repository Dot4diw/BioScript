# Install the edgeR package.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("edgeR", quietly = TRUE))
    BiocManager::install("edgeR")

if (!requireNamespace('EnhancedVolcano', quietly = TRUE))
    BiocManager::install('EnhancedVolcano')



# load library 
library("edgeR")
library("EnhancedVolcano")

# read the matrix data 
count_matrix <- as.matrix(read.csv("df_sc.csv", row.names = "gene"))

n_count <- 3 # normal sample count
t_count <- 3 # tumor sample count
# Create DGEList data class
sample_info <- colnames(count_matrix)
group_info <- c(rep("ctl", n_count), rep("trt", t_count))  # 3 ctl samples, 3 trt samples
group <- factor(group_info)

dge <- DGEList(counts = count_matrix, group = group)

# Filter out the genes with low counts
keep <- filterByExpr(y = dge)
# Recalculating the library sizes (keep.lib.sizes=FALSE) 
# for each sample is recommended following the filtering step. 
# After filtering, library sizes will be slightly changed for each sample.
dge <- dge[keep, , keep.lib.sizes=FALSE]


# Normalization and effective library sizes
# The normalization is performed using the TMM (Trimmed Mean of M-values) between-sample normalization method. 
dge <- calcNormFactors(object = dge, method = 'TMM')

# Model fitting and estimating dispersions
design <- model.matrix(~group)
dge <- estimateDisp(y = dge, design = design)

# Testing for differential gene expression
et <- exactTest(object = dge)

# topTags() function is useful to extract the table with adjusted p values (FDR). 
# The output table is ordered by p values.
# As the comparison of groups is trt-ctr, 
# the positive log fold change represents the gene is 
# more highly expressed in the trt condition as compared 
# to the ctr condition and vice versa.
top_degs <-  topTags(object = et, n = "Inf")
diff_res <- as.data.frame(top_degs)

# Export differential gene expression analysis table to CSV file.
write.csv(diff_res, file="condition_treatment_vs_control_dges.csv")

# Volcano plot
volcano_plot <- EnhancedVolcano(diff_res,
                                lab = rownames(diff_res),
                                title = "edgeR results",
                                subtitle = "Differential expression",
                                x = 'logFC',
                                y = 'PValue',
                                legendPosition = 'right',
                                legendLabSize = 12,
                                legendIconSize = 5.0)
pdf("edgeR_differential_expression_volcano_plot.pdf", width = 8, height = 6)
volcano_plot
dev.off()