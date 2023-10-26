# qqpcr！ 
###### 在学习qPCR实验时写的，可以一键获得相对表达量，可将结果拿到其他软件绘图
###### 适用于两个组每组两个样本，对照组在上


library(devtools)  

devtools::install_github('liuyuchenlab/qqpcr')  

library(qqpcr)  

###### 只能用PCR命名和示例文件的格式  

PCR <-read.csv('test.csv')  

###### Gapdh可以换成其他基因，这一步可以得到相对表达量和三个图片：  

###### 图片在plot里没有保存，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记。  

result <- qqpcr(PCR,'Gapdh')  

###### 保存相对表达量  

write.csv(result,'Gapdh_relative.csv') 

