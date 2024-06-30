library(ggplot2)
library(Rmisc)
library("RColorBrewer")
library(ggpubr)

data <- ToothGrowth
data$dose <- as.factor(data$dose)
data_summary <- summarySE(data, measurevar = "len", groupvars = c("supp","dose"))

## head(data)
##    len supp dose
## 1  4.2   VC  0.5
## 2 11.5   VC  0.5
## 3  7.3   VC  0.5
## 4  5.8   VC  0.5
## 5  6.4   VC  0.5
## 6 10.0   VC  0.5

p <- ggplot(data, aes(x = dose, y = len, fill = supp)) +
     geom_dotplot(binaxis = 'y', stackdir = 'center', 
               position = position_dodge(0.8)) +
     geom_errorbar(data = data_summary, aes(ymin = len - sd, ymax = len + sd), 
                width = 0.6, size = 0.5, position = position_dodge(0.8)) + 
     stat_summary(fun.y=mean, geom="crossbar", size=0.25, width = 0.7, 
                position = position_dodge(0.8), show.legend = F) +
     stat_compare_means(aes(group=supp), label = "p.signif", show.legend=FALSE) +
     scale_fill_manual(values = c("#d39200","#08519c") ) + 
     labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length") +
     theme_bw() + 
     theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(size=10),
          legend.key.size=unit(1.0,'cm'))

pdf("dotplot.pdf",width=6,height=5)
p
dev.off()
