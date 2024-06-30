library(Seurat)
library(ggplot2)
scRNA <- read("scRNA.rds")
cell_percent <- as.data.frame(prop.table(table(scRNA$annotation)))

colnames(cell_percent) <- c("Celltype","Frequence")
cell_percent <- cell_percent[sort(cell_percent$Frequence,decreasing=F,index.return=TRUE)$ix,]
cell_percent$Celltype <- factor(cell_percent$Celltype, levels = cell_percent$Celltype)
cell_percent

options(repr.plot.width=8, repr.plot.height=6)
ggplot(cell_percent,aes(x = Frequence, y = Celltype,fill=Frequence)) + geom_bar(stat="identity", width=0.8) + theme_bw() + 
    scale_fill_gradient(low = "#47df89", high = "#c74a14") + 
    theme(axis.text.y = element_text(size = 12, face = "bold",color="black"),
          axis.text.x = element_text(size = 12, face = "bold",color="black"),
          axis.title = element_text(size=14,face="bold"),
          legend.title=element_text(size=12,face="bold",color="black"),
          legend.text=element_text(size=12,face="bold",color="black"),
          plot.title = element_text(hjust = 0.5,face="bold",size=14),
          panel.border = element_rect(fill=NA,color="black", size=1.8)
         )

# multiple sample
cell_percent <- as.data.frame(prop.table(table(scRNA$project, scRNA$annotation)))
colnames(cell_percent) <- c("platform","Celltype","Frequence")
cell_percent

options(repr.plot.width=8, repr.plot.height=10)

ggplot(cell_percent, aes(fill = Celltype, y = Frequence, x = platform)) + geom_bar(position = "fill", stat = "identity") +
    ggtitle("The proportion of all cells in each plateform.") +
    theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
    theme(axis.text.y = element_text(size = 12,color="black"),
          axis.text.x = element_text(size = 12,color="black"),
          axis.title = element_text(size=14),
          legend.title=element_text(size=12,color="black"),
          legend.text=element_text(size=12,color="black"),
          plot.title = element_text(hjust = 0.5,size=14),
          panel.border = element_rect(fill=NA,color="black", size=1.8)) +
    scale_fill_manual(values = c("#D5C7E3",
                                 "#81958F",
                                 "#7B96DF",
                                 "#7DE4DB",
                                 "#C1E297",
                                 "#CC37E7",
                                 "#D6E64E",
                                 "#8066CC",
                                 "#D45083",
                                 "#DDB64C",
                                 "#D3E2D3",
                                 "#D95ED3",
                                 "#DAC190",
                                 "#6DE5A1",
                                 "#DB97DC",
                                 "#DC7659",
                                 "#663DDA",
                                 "#6FBEDE",
                                 "#7CE758"))
