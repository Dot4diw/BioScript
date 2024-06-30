library(tidyr)
library(reshape2)
library(ggplot2)
library(ggsignif)

# Input data format
###########################################################################################
#TCPA    TCPB    TCPC    TCPD    TCPE    TCPF
#113      89      49      49      101.5   76.5
#82       86      51      51      98.5    73.5
#125      102     52      49      114.5   89.5
###########################################################################################
#Store the data in the above format in a txt file

# set work directory
rm(list = ls())
setwd(".\\")
data <- read.table(file = 'test.txt', sep = "\t", header = T, row.names = NULL)

df <- melt(data, variable.name = "Group", value.name ="Values") 

#compared_list = list(c("Before", "After")) # Define two groups that require hypothesis testing
# Setting the fonts
windowsFonts(A=windowsFont("Times New Roman"),
    B=windowsFont("Arial"))

plot <-  ggplot(df, aes(Group, Values))+geom_boxplot(aes(fill=Group),size=1, color="black") + 
               stat_boxplot(geom = "errorbar", width=0.2, color = "black") + # errorbar 
               #scale_y_continuous(breaks=seq(0, range(df$Values)[2] + 0.2*range(df$Values)[2]) ,range(df$Values)[2]/10)
               #geom_signif(comparisons = compared_list, test = t.test, map_signif_level=T, size = 1, textsize = 5) + # add signif of target
               geom_point(aes(fill=Group), shape=21, size=2, alpha = 1/5, colour = "black", fill="black") +
               #geom_jitter(aes(fill=group), width=0.5, shape=21, size=4) + # 
               theme_bw() + 
               #scale_fill_manual(values=c("#8064a2","#f8a968")) + # Custom the box fill colors
               theme(#panel.grid.major = element_blank(), 
                     #panel.grid.minor = element_blank(),
                     panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     #legend.title=element_text(size=14,face="bold",color="black"),
                     text = element_text(family = "B", face = "bold"),
                     axis.text.x = element_text(size = 14, face = "bold",color="black",vjust = 1, hjust=1, angle = 45), 
                     axis.text.y = element_text(size = 14, face = "bold",color="black"),
                     axis.title = element_text(size=14,face="bold"),
                     legend.title = element_text(size=14,face="bold",color="black"),
                     legend.text = element_text(size=14,face="bold",color="black"),
                     plot.title = element_text(hjust = 0.5,face="bold",size=14))

png("boxplot.png",units="in",width=6, height=6,res=600)
plot
dev.off()
