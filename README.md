# qqpcr！ 
###### 在学习qPCR实验时写的，可以一键获得相对表达量
###### 最适用于两个组每组两个样本，因为我学的时候就这样，多了估计也行

install.packages("devtools")

library(devtools)  

###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

devtools::install_github('liuyuchenlab/qqpcr')  


library(qqpcr)  

###### 只能用PCR命名or not  
###### 内置数据
![image](https://github.com/user-attachments/assets/0cf6077a-b55e-4271-8b1f-2c25b360efe4)


PCR <- read.csv('test.csv')  

###### Gapdh可以换成其他基因，这一步可以得到相对表达量和三个图片：  

###### 图片在plot里没有保存，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记 
###### p值和显著性标记的位置在分组超过3个时会比较乱
###### p值计算使用近似T检验，但我感觉算得不太对，可能是因为示例数据样本量每组只有两个，这边还是建议只拿relative数据

result <- qqpcr(PCR, reference_gene = "Gapdh", control_group = "control")  

##### 柱形图  
![image](https://github.com/user-attachments/assets/e6645e0a-d706-4de7-af60-32a649079b49)
  

##### p值柱形图  
![image](https://github.com/user-attachments/assets/df12fb80-c903-488c-8c86-b77b85bd4205)


##### 显著性标记柱形图  

![image](https://github.com/user-attachments/assets/cbad41dc-86df-4f37-bfaf-5dff88143e38)

###### 函数运行时会自动保存行列转换的CT值、相对表达值、统计数据和p值
![image](https://github.com/user-attachments/assets/949b6eb3-e931-4b80-86fd-720bb36115e5)


##### 也可以添加自定义颜色和截断

custom_colors <- c("#FF0000", "#00FF00", "#0000FF")  # 用户自定义的颜色

result <- qqpcr(PCR, 'Gapdh', 'control', breaks = list(c(4, 5), c(14, 15), c(20, 21)), custom_colors = custom_colors)

###### 也可以自定义保存

![image](https://github.com/user-attachments/assets/44836f5f-9995-40b2-b1ba-1ef4b036154c)





