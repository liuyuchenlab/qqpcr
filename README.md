# qqpcr
library(devtools)
devtools::install_github('liuyuchenlab/qqpcr')
library(qqpcr)
PCR <-read.csv('test.csv')#只能用这个命名和文件格式
result <- qqpcr(PCR,'Gapdh')#Gapdh可以换成其他基因
write.csv(result,'Gapdh_relative.csv')#获得relative expression level。
