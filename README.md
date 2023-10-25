# qqpcr
library(devtools)  

devtools::install_github('liuyuchenlab/qqpcr')  

library(qqpcr)  

PCR <-read.csv('test.csv')#只能用这个命名和文件格式  

result <- qqpcr(PCR,'Gapdh')#Gapdh可以换成其他基因，这一步可以得到relative和三个图片，图片在plot里没有保存，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记。  

write.csv(result,'Gapdh_relative.csv')#获得relative expression level。  

