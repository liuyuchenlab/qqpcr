# qqpcr！ 
###### 在学习qPCR实验时写的，可以一键获得相对表达量，可将结果拿到其他软件绘图
###### 适用于两个组每组两个样本，对照组在上


library(devtools)  

devtools::install_github('liuyuchenlab/qqpcr')  

###### 如果连接失败，尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

library(qqpcr)  

###### 只能用PCR命名or not  

PCR <-read.csv('test.csv')  

###### Gapdh可以换成其他基因，这一步可以得到相对表达量和三个图片：  

###### 图片在plot里没有保存，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记。  

result <- qqpcr(PCR,'Gapdh')  

##### 柱形图  
![image](https://github.com/liuyuchenlab/qqpcr/assets/131852185/0d421abe-6cc9-47fc-be97-1622fe80c0c9)  

##### p值柱形图  
![image](https://github.com/liuyuchenlab/qqpcr/assets/131852185/df40ac97-74ea-4642-87c4-be2d18c1851e)  

##### 显著性标记柱形图  

![image](https://github.com/liuyuchenlab/qqpcr/assets/131852185/303d40f6-083b-4fe0-993e-3954ed8609cb)  



###### 保存相对表达量  

write.csv(result,'Gapdh_relative.csv') 

