# install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library(ggVennDiagram)
library(ggvenn)

genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350),
  )
# The test file content as shown below
#  AA       BB     CC       DD
# Gene1     Gene11  Gene4   Gene10
# Gene2     Gene12  Gene5   Gene11
# Gene3     Gene13  Gene6   Gene12
# Gene4     Gene14  Gene7   Gene13
# Gene5     Gene15  Gene8   Gene14
# Gene6     Gene16  Gene20  Gene15
# Gene7     Gene17  Gene21  Gene16
# Gene8     Gene3   Gene22  Gene31
# Gene9     Gene4   Gene23  Gene32
# Gene10    Gene5   Gene24  Gene33
# Gene11    Gene9   Gene25  Gene34
# Gene12    Gene1           Gene35
# Gene13    Gene2           Gene36
# Gene14                    Gene37
#                           Gene38
#                           Gene39
#                           Gene40
#                           Gene41
#                           Gene42
#                           Gene43
# 

aa <- read.table("test.txt", sep = '\t', header = TRUE)
xx <- as.list(aa)

#如何去除R的list中纯""的元素

for (i in names(xx))
{
		if ("" %in% xx[[i]]) 
        {
            xx[[i]] = xx[[i]][-which(xx[[i]]=="")]
        }
}

overlap_list <- process_region_data(Venn(xx))

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

