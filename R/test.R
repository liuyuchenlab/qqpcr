# 初始化环境
renv::init()
# 如果包为使用时必须的，则需要设置 type = "Imports"
usethis::use_package(package = "ggpubr", type = "Imports")
usethis::use_package(package = "psych", type = "Imports")
usethis::use_package(package = "tidyverse", type = "depends")
usethis::use_package(package = "ggbreak", type = "Imports")
usethis::use_package(package = "dplyr", type = "depends")
usethis::use_package(package = "rstatix", type = "depends")
# 保存当前所用的包环境
renv::snapshot()
#创建 README
usethis::use_readme_rmd()
# 通过 RStudio 的 File > New File > R Script 也一样
file.create("R/qqpcr.R")
#生成内置数据
PCR <-read.csv("test.csv")
usethis::use_data(PCR)
#自定义函数#
qqpcr <- function(PCR, reference_gene, control_group, breaks = NULL, custom_colors = NULL) {
  suppressMessages(library(ggpubr))
  suppressMessages(library(psych))
  suppressMessages(library(tidyverse))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggbreak))
  suppressMessages(library(rstatix))
  # 检查 control_group 是否在 PCR$group 中
  if (!(control_group %in% unique(PCR$group))) {
    stop(paste("Error: The control group", control_group, "does not exist in the PCR data. Please check the group names."))
  }

  PCR_Ct <- PCR
  all_col <- colnames(PCR)
  elements_to_remove <- c(reference_gene, "group")
  target_genes <- all_col[!(all_col %in% elements_to_remove)]

  # 创建空向量来存储结果
  results <- vector("list", length(target_genes))

  # 循环处理每个目标基因的数据
  for (i in seq_along(target_genes)) {
    dct <- PCR[[target_genes[i]]] - PCR[[reference_gene]]
    control_mean <- mean(dct[PCR$group == control_group], na.rm = TRUE)
    ddct <- dct - control_mean
    relative <- 2^-ddct
    col_name <- target_genes[i]
    PCR[[col_name]] <- relative
    results[[i]] <- data.frame(Gene = target_genes[i], Relative = relative, group = PCR$group)
  }

  # 创建两个数据框：一个包含原始Ct值，另一个包含计算的相对表达量
  PCR_relative <- PCR[, c("group", target_genes)]

  # 将所有基因的数据合并到一个长格式数据框中
  PCR_long <- bind_rows(results)

  # 确保分组和基因的顺序
  PCR_long <- PCR_long %>%
    mutate(
      Gene = factor(Gene, levels = target_genes),  # 确保基因按照输入顺序显示
      group = factor(group, levels = unique(PCR$group))  # 确保分组按照输入顺序显示
    )

  # 计算每个基因每个组的均值、标准误差等统计信息
  summary_stats <- PCR_long %>%
    group_by(Gene, group) %>%
    summarize(
      mean_relative = mean(Relative, na.rm = TRUE),
      sd_relative = sd(Relative, na.rm = TRUE),
      var_relative = var(Relative, na.rm = TRUE),  # 方差
      n = n(),
      se_relative = sd_relative / sqrt(n),
      upper_ci = mean_relative + se_relative,
      lower_ci = mean_relative - se_relative,
      .groups = 'drop'
    )
  # 统计检验：仅保留非list列，避免write.csv报错
  stat.test <- PCR_long %>%
    group_by(Gene) %>%
    t_test(Relative~group, ref.group = control_group) %>%
    adjust_pvalue() %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Gene") %>%
    # 筛选出需要的列（移除list类型列）
    select(Gene, group1, group2, p, p.adj, p.adj.signif, y.position, x, xmin, xmax)
  # 找到最高的 upper_ci
  max_upper_ci <- max(summary_stats$upper_ci, na.rm = TRUE)



  # 默认颜色设置
  default_colors <- c("#B22222","#4A70A9","#FBCE6A","#77AE43", "#EDB021","#acc", "#c34a56", "#c7b363", "#2f4858", "#8c564b", "#ff7f0e", "#845ec2", "#d65db1", "#2ca02c", "#d62728", "#00c9a7")

  # 使用用户自定义颜色（如果提供了）或默认颜色
  group_colors <- if (is.null(custom_colors)) {
    setNames(default_colors, unique(PCR$group))
  } else {
    setNames(custom_colors, unique(PCR$group))
  }

  # 创建基础柱状图
  p_basic <- ggbarplot(PCR_long,
                       x = 'Gene',
                       y = 'Relative',
                       fill = 'group',
                       color = 'group',
                       palette =  group_colors,
                       add = "mean_sd",
                       add.params = list(width = 0.4),
                       position = position_dodge(width = 0.8),
                       xlab = "", ylab = 'Relative expression', legend = 'right',
                       ggtheme = theme_classic()) +
    scale_y_continuous(limits = c(0, max_upper_ci * 1.6), expand = c(0, 0)) + # 设置 y 轴范围，并增加空白以显示标记
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key.height = unit(1, "lines"),
      legend.key.spacing.y = unit(0.5, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.ticks = element_line(size = 0.8),
      axis.ticks.length = unit(2, "mm"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "mm") ,
      axis.title.y = element_text(size = 15, face = "bold", vjust = 4),
      axis.line.x = element_line(linetype = 1, color = "black", size = 0.8),
      axis.line.y = element_line(linetype = 1, color = "black", size = 0.8),
      plot.title = element_text(size = 15, color = "black",  hjust = 0.5,vjust = 2)) +
    labs(title = paste("Relative to", control_group, "(", reference_gene, ")")) +
    guides(fill = guide_legend(byrow = TRUE))

  # 如果提供了 breaks 参数
  if (!is.null(breaks) && length(breaks) > 0) {
    # 确保 breaks 是一个合适的列表
    if (!is.list(breaks) || !all(sapply(breaks, function(x) length(x) == 2 && is.numeric(x)))) {
      stop("The 'breaks' parameter should be a list of numeric vectors with exactly two elements each.")
    }
    for (break_set in breaks) {
      p_basic <- p_basic + scale_y_break(breaks = break_set, scales = "free", space = .5) +
        theme(
          axis.line.y.right = element_blank(),  # 去掉右边的y轴线
          axis.ticks.y.right = element_blank(),  # 去掉右边的y轴刻度线
          axis.text.y.right = element_blank()    # 去掉右边的y轴标签
        )
    }
  }

  # 对照组不显示的偏移调整
  dodge <- position_dodge(width = 0.8)

  #添加统计信息
  p_signif <- p_basic + stat_pvalue_manual(stat.test,label="p.adj.signif",hide.ns=F,#label="p"orlabel="p.adj"
                                           tip.length=0.01,label.size=4)
  p_padj <- p_basic + stat_pvalue_manual(stat.test,label="p.adj",hide.ns=F,#label="p"orlabel="p.adj"
                                         tip.length=0.01,label.size=4)
  # 打印图形
  print(p_basic)
  print(p_padj)
  print(p_signif)
  # 数据格式转换
  PCR_Ct <- t(PCR_Ct)
  PCR_relative <- t(PCR_relative)

  # 移除列名
  colnames(PCR_Ct) <- NULL
  colnames(PCR_relative) <- NULL

  # 获取当前日期
  current_date <- format(Sys.Date(), "%Y-%m-%d")

  # 创建日期文件夹
  dir.create(current_date, showWarnings = FALSE)

  # 创建带日期和内参基因名的文件名
  ct_filename <- file.path(current_date, paste0(control_group, "_", reference_gene, "_PCR_Ct_", current_date, ".csv"))
  relative_filename <- file.path(current_date, paste0(control_group, "_", reference_gene, "_PCR_relative_", current_date, ".csv"))

  # 创建带日期的p值文件名
  p_value_filename <- file.path(current_date, paste0(control_group, "_", reference_gene, "_p_values_", current_date, ".csv"))


  # 保存 p 值数据到 CSV 文件
  stat.test <- as.data.frame(stat.test)
  write.csv(stat.test, file = p_value_filename, row.names = FALSE)
  # 确保列名不被写入文件中
  write.csv(PCR_Ct, file = ct_filename, row.names = TRUE)
  write.csv(PCR_relative, file = relative_filename, row.names = TRUE)

  # 创建带日期的统计信息文件名
  summary_stats_filename <- file.path(current_date, paste0(control_group, "_", reference_gene, "_summary_stats_", current_date, ".csv"))
  # 保存 summary_stats 数据到 CSV 文件
  write.csv(summary_stats, file = summary_stats_filename, row.names = FALSE)

  # 保存图形
  ggsave(file.path(current_date, paste0(control_group, "_", reference_gene, "_Relative_expression_", current_date, ".pdf")), plot = p_basic, width = 10, height = 6)

  # 返回包含Ct值、相对表达量、统计信息和图形的列表
  return(list(PCR_Ct = PCR_Ct,
              PCR_relative = PCR_relative,
              p_value = stat.test,
              summary_stats = summary_stats,  # 返回统计信息
              p_basic = p_basic,
              p_with_p = p_padj,
              p_with_signif = p_signif))
}


# 示例数据读取和函数调用
# PCR <- read.csv("/mnt/data/your_file_name.csv")
# result <- qqpcr(PCR, reference_gene = "Gapdh", control_group = "control")
# print(result$p)
# print(result$p_value_plot)
# print(result$p_signif_plot)

# 示例数据读取和函数调用
# PCR <- read.csv("/mnt/data/your_file_name.csv")
# result <- qqpcr(PCR, reference_gene = "Gapdh", control_group = "control")
data(PCR)
#PCR <-read.csv("1.csv")
results  <- qqpcr(PCR,'Rps18',"treat")
custom_colors <- c("#FF0000", "#00FF00", "#0000FF")  # 用户自定义的颜色
# 示例 breaks 参数定义了多个截断
result <- qqpcr(PCR, 'Gapdh', 'control', breaks = list(c(5,6), c(14, 15), c(20, 21)), custom_colors = custom_colors)
