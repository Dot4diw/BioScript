library(ggplot2)
library(tidyverse)
library(ggrepel)
library(dplyr)

df <- read.table("diff_combined_4classSample.txt",header = T)

df$label <- ifelse(df$p_val_adj<0.001,"adjust P-val<0.001","adjust P-val>=0.001")

group_by(df,vs_info) %>% top_n(n = 50, wt = abs(avg_log2FC)) -> top10


vs_levels <- unique(df$vs_info)

df$vs_info <- factor(df$vs_info,levels=vs_levels)

bar_value_p <- c()
bar_value_n <- c()
for (i in 1:length(unique(df$vs_info)))
{
    bar_value_p <- append(bar_value_p, max(df[(df$vs_info == vs_levels[[i]]),2]))
    bar_value_n <- append(bar_value_n, min(df[(df$vs_info == vs_levels[[i]]),2]))
    
}
dfbar_p <-data.frame(x=vs_levels, y=bar_value_p)
dfbar_n <-data.frame(x=vs_levels, y=bar_value_n)


dfcol<-data.frame(x=vs_levels, y=0, label=vs_levels)
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F")


p <- ggplot() + geom_jitter(data = df,
                            aes(x = vs_info, y = avg_log2FC, color = label),
                            size = 2,
                            width =0.4) +
    geom_col(data = dfbar_p,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar_n,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6) +
    geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F) +
    geom_text_repel(data=top10,
                    aes(x=vs_info,y=avg_log2FC,label=gene),
                    force = 10,
                    min.segment.length = Inf,
                    #vjust = -0.5,hjust = 0,
                    point.size = NA)

p + theme_bw() + theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12, color="black", vjust = 1, hjust=1, angle = 45), 
    axis.text.y = element_text(size = 12, ,color="black"),
    axis.title = element_text(size=12,),
    legend.title = element_text(size=12, color="black"),
    legend.text = element_text(size=12, color="black"),
    plot.title = element_text(hjust = 0.5,size=12)) 
