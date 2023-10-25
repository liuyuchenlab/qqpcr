#####qqpcr!
# Hello, world!
qqpcr <- function(PCR,reference_gene){
  suppressMessages(library(ggpubr))
  suppressMessages(library(psych))
  suppressMessages(library(tidyverse))
  suppressMessages(library(dplyr))
  all_col <- colnames(PCR)
  elements_to_remove <- c(reference_gene, "group")
  target_genes <- all_col[!(all_col%in% elements_to_remove)]#目的基因
  # 创建空向量来存储结果
  mrna_results <- vector("list", length(target_genes))
  # 循环处理每个目标基因的数据
  for (i in seq_along(target_genes)) {
    # 计算∆Ct（目标基因Ct - 内参基因Ct）
    dct <- PCR[[target_genes[i]]] - PCR[[reference_gene]]
    # 计算∆∆Ct（∆Ct - 对照组Ct均值）
    ddct <- dct - mean(dct[1:2])
    # 计算mRNA相对表达量
    mrna <- 2^-ddct
    # 将结果存储在列表中
    mrna_results[[i]] <- mrna
  }
  # 将结果添加回原始数据框中
  for (i in seq_along(target_genes)) {
    col_name <- paste(target_genes[i], "mrna", sep = "_")
    PCR[[col_name]] <- mrna_results[[i]]
  }
  # 输出包含多基因qPCR数据分析结果的数据框
  print(PCR)
  # 计算∆Ct、∆∆Ct和mRNA相对表达量，并创建新的PCR_long数据框
  PCR_long <- data.frame()
  for (i in seq_along(target_genes)) {
    dct <- PCR[[target_genes[i]]] - PCR[[reference_gene]]
    ddct <- dct - mean(dct[1:2])
    mrna <- 2^-ddct
    col_name <- paste(target_genes[i], "mrna", sep = "_")
    PCR_long <- rbind(PCR_long, data.frame(Gene = target_genes[i], Relative_mRNA = mrna, group = PCR$group))
  }
  print(PCR_long)
  # 绘制多基因的柱状图，按分组着色
  p <- ggbarplot(PCR_long,
                 'Gene',  # 将 'Gene' 替换为实际的基因名列名
                 'Relative_mRNA',  # 将 'Relative_mRNA' 替换为实际的相对mRNA表达量列名
                 fill = 'group',
                 color = 'group',
                 palette = "jco",
                 add = "mean_sd",
                 position = position_dodge(width = 0.8), # 不堆积，按分组分开
                 xlab = F, ylab = 'Relative mRNA expression', legend = 'right',
                 ggtheme = theme_bw())

  print(p)
  #截断
  #suppressMessages(library(ggbreak))
  #p+scale_y_break(c(1.5, 15))
  ###统计检验
  # 创建空数据框来存储结果
  result_df <- data.frame(Gene = character(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

  # 循环处理每个目标基因的数据
  for (i in seq_along(target_genes)) {
    # 使用dplyr进行分组，并进行t检验，但先检查每个组是否有足够的数据点
    t_test_result <- PCR %>%
      filter(!is.na(.[[target_genes[i]]])) %>%  # 过滤掉缺失值
      group_by(group) %>%
      summarize(p_value = t.test(.[[target_genes[i]]])$p.value)

    # 将结果添加到结果数据框中
    result_df <- rbind(result_df, data.frame(Gene = target_genes[i], p_value = t_test_result$p_value))
  }

  # 将p值添加到pcr_long数据框中
  PCR_long_with_p <- PCR_long %>%
    left_join(result_df, by = "Gene",relationship = "many-to-many") %>%
    mutate(p_value = ifelse(!is.na(p_value), p_value, ""))
  # 绘制多基因的柱状图，并标记p值
  p <- ggbarplot(PCR_long_with_p,
                 'Gene',
                 'Relative_mRNA',
                 fill = 'group',
                 color = 'group',
                 palette = "jco",
                 add = "mean_sd",
                 position = position_dodge(width = 0.8),
                 xlab = F, ylab = 'Relative mRNA expression', legend = 'right',
                 ggtheme = theme_bw())
  #  geom_jitter(color='black',size=2)
  p.format <- p + stat_compare_means(aes(group = group), label = "p.format", label.y = max(PCR_long$Relative_mRNA)*1.1,method = "t.test")
  p.signif <- p + stat_compare_means(aes(group = group), label = "p.signif", label.y = max(PCR_long$Relative_mRNA)*1.1,method = 't.test')
  print(p.format)
  print(p.signif)
  return(PCR)
}
#ok
