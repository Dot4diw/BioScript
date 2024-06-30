# Install the packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!requireNamespace('EnhancedVolcano', quietly = TRUE))
    BiocManager::install('EnhancedVolcano')


# load library 
library("DESeq2")
library("EnhancedVolcano")

# Read the matrix data.
count_matrix <- as.matrix(read.csv("df_sc.csv", row.names = "gene"))

n_count <- 3 # normal sample count
t_count <- 3 # tumor sample count

# creat the sample information (metadata).
coldata <- data.frame(
   sample = colnames(count_matrix),
   condition = c(rep("control", n_count), rep("infected", t_count)), 
   row.names = "sample" )

coldata$condition <- as.factor(coldata$condition)

# Do Differential Gene Expressoin Analysis.
if ( all(rownames(coldata) == colnames(count_matrix)) ){
    # construct DESeqDataSet for DGE analysis
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, 
                              design = ~ condition)

    # remove the genes which have < 10 reads 
    # (this can vary based on research goal) in total across all the samples.
    dds <- dds[rowSums(counts(dds)) >= 10,]

    # set control condition as reference
    dds$condition <- relevel(dds$condition, ref = "control")

    dds <- DESeq(dds)

    dge_result <- results(dds)
    dge_result_order <- dge_result[order(dge_result$padj),]
    diff_res <- as.data.frame(dge_result_order)

    write.csv(diff_res, file="condition_treatment_vs_control_dges.csv")

    # Get the normalized counts
    normalized_counts <- counts(dds, normalized=TRUE)
    write.csv(as.data.frame(normalized_counts), file="counts_matrix_normalized.csv")

    # Volcano plot
    volcano_plot <- EnhancedVolcano(diff_res,
                                    lab = rownames(diff_res),
                                    title = "edgeR results",
                                    subtitle = "Differential expression",
                                    x = 'log2FoldChange',
                                    y = 'pvalue',
                                    legendPosition = 'right',
                                    legendLabSize = 12,
                                    legendIconSize = 5.0)
    pdf("DESeq2_differential_expression_volcano_plot.pdf", width = 8, height = 6)
    volcano_plot
    dev.off()
} else {
    print("The names of the columns in the counting matrix are not in the same order as the sample names.")
}
