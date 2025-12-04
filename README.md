# qqpcr！ 
###### 在学习qPCR实验时写的，可以一键获得相对表达量
##### 安装

```
install.packages("devtools")
library(devtools)
```

###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

```
devtools::install_github('liuyuchenlab/qqpcr')  
library(qqpcr)  
```

###### 内置数据格式
```
data(PCR)
```
<img width="693" height="549" alt="image" src="https://github.com/user-attachments/assets/b36a27be-41ac-4fff-85ac-6dad5dfab318" />


###### 读取数据（只能用PCR命名or not）  
```
#PCR <- read.csv('test.csv')
PCR <- read.xlsx("test.xlsx") #推荐此方式
```
###### Gapdh可以换成其他基因，这一步可以得到相对表达量和三个图片：  

###### 图片在plot里，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记 
###### p值计算使用T检验得到的padj，可以使用relative数据去其他知名软件计算并绘图
```
results <- qqpcr(PCR, reference_gene = "Gapdh", control_group = "control")  
```
##### 柱形图  
<img width="936" height="512" alt="image" src="https://github.com/user-attachments/assets/cbf53aaf-25f2-4a1e-a4b5-cc01ae56aa52" />

##### p值柱形图  
<img width="936" height="512" alt="image" src="https://github.com/user-attachments/assets/25d4ecaf-0f4e-4eb1-8277-5042444b635f" />

##### 显著性标记柱形图  
<img width="936" height="512" alt="image" src="https://github.com/user-attachments/assets/21e99ab7-44c6-45c1-b05d-aef4be241ecd" />

###### 函数运行时会自动保存行列转换的CT值、相对表达值、统计数据和p值还有一张basic图片到自动创建的文件夹中，文件夹名为当天的日期。
![image](https://github.com/user-attachments/assets/5f8052c7-84cc-4bd9-92d3-ca7012261c1b)

##### 默认只有16个颜色，也可以添加自定义颜色和截断
```
custom_colors <- c("#FF0000", "#00FF00", "#0000FF")  # 用户自定义的颜色

result <- qqpcr(PCR, 'Gapdh', 'control', 
                breaks = list(c(2,4), c(14, 20)), #设置截断位置
                custom_colors = custom_colors  #设置自定义颜色
                )
```
##### 示例展示
<img width="936" height="512" alt="image" src="https://github.com/user-attachments/assets/5e1882e0-fc71-4747-a5b7-1543d83e23d5" />

##### 也可以自定义保存
![image](https://github.com/user-attachments/assets/b661ec38-b609-431c-b799-658d3f0577bc)

### 快去试试吧！






