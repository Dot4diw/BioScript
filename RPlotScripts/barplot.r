library(ggplot2)

df <- data.frame(category=c('A', 'B', 'C', 'D', 'E'),
                 value=c(12, 17, 30, 22, 19),
                 sd=c(4, 5, 7, 4, 2))

fill_color=c("#357ebd","#8c57a2","#cbcacb","#165f83","#f39b7f")

p <- ggplot(data=df, 
            aes(x=category,y=value)) + theme_bw()  +
            geom_bar(stat="identity",
            color="black", 
            fill=fill_color,
            size=0.5) + geom_errorbar(aes(x=category, ymin=value-sd, ymax=value+sd), width=0.2,size=0.5) + 
            theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
p

p2<-ggbarplot(ToothGrowth, x="dose", y="len", add = "mean_se", color = "black", fill="supp",
    palette = "jco", position = position_dodge(0.8))+
    stat_compare_means(aes(group=supp), label = "p.signif", label.y = 29)
p2
