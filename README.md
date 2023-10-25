# qqpcr！ 
###### 在学习qPCR实验时写的，可以便捷获得relative expression方法，可将结果拿到其他软件绘图
###### 这个代码应该只能用于两个组每组两个样本，对照组放上面。


library(devtools)  

devtools::install_github('liuyuchenlab/qqpcr')  

library(qqpcr)  

PCR <-read.csv('test.csv')#只能用这个命名和文件格式  

###### Gapdh可以换成其他基因，这一步可以得到relative和三个图片：  

###### 图片在plot里没有保存，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记。  

result <- qqpcr(PCR,'Gapdh')  

###### 获得relative expression level。 

write.csv(result,'Gapdh_relative.csv') 

