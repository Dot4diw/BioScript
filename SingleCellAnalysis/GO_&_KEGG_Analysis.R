rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(dplyr)
library(GseaVis)

# Go Plot
goplot <- function(go_obj, topn=20,reg="UP") 
{
  ego <- as.data.frame(go_obj)
  top_ego <- ego %>% group_by(ONTOLOGY) %>% arrange(ONTOLOGY, p.adjust) %>% slice(1:topn)
  # RichRation = Fold Enrichment
  top_ego$RichRatio <- sapply(strsplit(top_ego$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])) / sapply(strsplit(top_ego$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  top_ego <- top_ego[order(top_ego$RichRatio),]
  top_ego$Description<-factor(top_ego$Description,levels = top_ego$Description)
  
  plot <- ggplot(top_ego,aes(x=Description,y=RichRatio)) + 
    geom_point(aes(size=Count,color= -log10(p.adjust))) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    coord_flip() + facet_grid(ONTOLOGY~., scale="free",labeller = labeller(ONTOLOGY = c("BP" = "Biological Process", 
                                                                                        "CC" = "Cellular Component", 
                                                                                        "MF" = "Molecular Function"))
    ) + ggtitle(paste0(i,": ",reg, "\nTOP GO TERMS")) + 
    scale_colour_gradient(low="yellow",high="red",space = "Lab") + theme_bw() +
    theme(strip.text.y = element_text(size = 10, colour = "#000000")) + 
    theme(strip.background.y = element_rect(fill = "#FFFFFF", colour = "#000000"), plot.title = element_text(size = 8, hjust = 0))
  return(plot)
}

# KEGG plot
keggplot <- function(kegg_obj, topn=30, reg="UP")
{
  ekegg <- as.data.frame(kegg_obj)
  ekegg <- ekegg[order(ekegg$p.adjust), ]
  top_kegg <- head(ekegg,topn)
  top_kegg$RichRatio <- sapply(strsplit(top_kegg$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])) / sapply(strsplit(top_kegg$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  top_kegg <- top_kegg[order(top_kegg$RichRatio),]
  top_kegg$Description<-factor(top_kegg$Description,levels = top_kegg$Description)
  
  plot <- ggplot(top_kegg,aes(x=Description,y=RichRatio)) + 
    geom_point(aes(size=Count,color= -log10(p.adjust))) + ggtitle(paste0(i, ":",reg,"\nTOP KEGG PATHWAYS")) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    coord_flip() + scale_colour_gradient(low="yellow",high="red",space = "Lab") + theme_bw() + theme(plot.title = element_text(size = 8, hjust = 0))
  return(plot)
}

# KEGG plot
reactomplot <- function(kegg_obj, topn=30, reg="UP")
{
  ekegg <- as.data.frame(kegg_obj)
  ekegg <- ekegg[order(ekegg$p.adjust), ]
  top_kegg <- head(ekegg,topn)
  top_kegg$RichRatio <- sapply(strsplit(top_kegg$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])) / sapply(strsplit(top_kegg$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  top_kegg <- top_kegg[order(top_kegg$RichRatio),]
  top_kegg$Description<-factor(top_kegg$Description,levels = top_kegg$Description)
  
  plot <- ggplot(top_kegg,aes(x=Description,y=RichRatio)) + 
    geom_point(aes(size=Count,color= -log10(p.adjust))) + ggtitle(paste0(i, ":",reg,"\nTOP ReactomPA")) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    coord_flip() + scale_colour_gradient(low="yellow",high="red",space = "Lab") + theme_bw() + theme(plot.title = element_text(size = 8, hjust = 0))
  return(plot)
}


deg <- read.csv("deg.csv", row.names = 1)
head(deg)

# DO GO&KEGG
# # Convert the gene Symbol to Entrez gene ID
orgdb = "org.Mm.eg.db"
organ = "mmu"

i = "Samplename"
reg = "UP"

geneid <- bitr(deg$gene,
               fromType = 'SYMBOL',
               toType   = c('ENTREZID'),
               OrgDb    = orgdb)

ego <- enrichGO(gene          = geneid$ENTREZID,
                OrgDb         = orgdb,
                keyType       = 'ENTREZID',
                ont           = "All",
                pAdjustMethod = "BH",
                readable      = TRUE,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)

ekegg <- enrichKEGG(gene          = geneid$ENTREZID,
                    keyType       = 'kegg',
                    organism      = organ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)


ekegg <- setReadable(ekegg, OrgDb = orgdb, keyType="ENTREZID")
ekegg@result$Description <- str_remove(ekegg@result$Description, " - Mus musculus \\(house mouse\\)")

write.table(ego, file = paste0(i,"_", reg,  "_GO.txt"), sep="\t", row.names =F, quote = F) 
write.table(ekegg, file =  paste0(i,"_",reg,"_KEGG.txt"), sep="\t", row.names =F, quote = F)

saveRDS(ego, file=paste0(i,"_",reg, "_GO.rds"))
saveRDS(ekegg, file=paste0(i,"_",reg, "_KEGG.rds"))


plot1 <- goplot(ego,topn = 10, reg=reg)
#plot1 <- dotplot(ego, split="ONTOLOGY",showCategory = 20,title = 'GO Ontology Dotplot',label_format=50) + 
#  facet_grid(ONTOLOGY~., scale="free") + 
#  theme(strip.text.y = element_text(size = 15, colour = "#000000")) + 
#  theme(strip.background.y = element_rect(fill = "#FFFFFF", colour = "#000000"))

ggsave(file = paste0(i, "_",reg, "_GO.pdf"), plot = plot1, width = 7, height = 6 )
plot1
rm(plot1)

plot2 <- keggplot(ekegg, topn = 30, reg=reg)
plot2
#plot2 <- dotplot(ekegg, showCategory = 30,title = 'KEGG Dotplot',label_format=50) + scale_color_continuous(low="#6D2884", high="#D099F2")
ggsave(file = paste0(i,"_",reg, "_KEGG.pdf"), plot = plot2, width = 7, height = 6)

rm(plot2)
