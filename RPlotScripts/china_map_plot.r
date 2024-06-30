# Install Packages
install.packages("mapchina")
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("xmc811/mapchina", ref = "dev")

# Import Packages
library(mapchina)
library(sf)
library(tidyverse)
library(ggplot2)

sf_use_s2(FALSE)
# Read my data
mydata <- read.csv("mymapdata.csv")
#  Name_Province        Count
#  新疆维吾尔自治区     29
#  西藏自治区   0
#  青海省       4
#  甘肃省       85
#  宁夏回族自治区       35
#  陕西省       60
#  四川省       635
#  重庆市       91
#  湖北省       179
#  河南省       38
#  山西省       59
#  河北省       45
#  山东省       188
#  北京市       11
#  天津市       1
#  内蒙古自治区 37
#  辽宁省       20
#  吉林省       112
#  黑龙江省     21
#  江苏省       30
#  上海市       1
#  安徽省       159
#  浙江省       179
#  福建省       227
#  江西省       348
#  湖南省       1683
#  贵州省       1326
#  广西壮族自治区       213
#  广东省       153
#  海南省       57
#  云南省       4010
#  香港特别行政区       0
#  台湾省       0
#  澳门特别行政区       0

# Extract the Province data.
df <- china %>%
        filter(Code_Province %in% as.character(unique(china$Code_Province)))
        
df <- df %>%
        group_by(Name_Province) %>%
        summarise(geometry = st_union(geometry))


df$Count <- rep(0, nrow(df))

# Add our data into the Province Map data.
for (i in 1:length(df$Name_Province))
{
  print(subset(mydata, mydata$Name_Province == df[i,]$Name_Province)$Count)
  df[i,]$Count <- subset(mydata, mydata$Name_Province == df[i,]$Name_Province)$Count
}

ggplot(data = df) +
    geom_sf(aes(fill = Count)) + theme_bw() + scale_fill_gradient(
        low="#fff5f0",high = "red") +  geom_sf_label(aes(label = Count))
ggplot(data = df) +
    geom_sf(aes(fill = Count)) + theme_void() + scale_fill_gradient(
        low="#fff5f0",high = "red") +  geom_sf_text(aes(label = Count))
