#####qqpcr!
# Hello, world!
#' Title
#'
#' @param PCR PCR_name
#' @param reference_gene reference_gene
#' @param control_group control_group
#' @param breaks breaks
#' @param custom_colors custom_colors
#'
#' @return relative expression and plots
#' @export
#'
#' @examples
#' data(PCR)
#' results  <- qqpcr(PCR,'Gapdh',"control")
qqpcr <- function(PCR, reference_gene, control_group, breaks = NULL, custom_colors = NULL) {
  suppressMessages(library(ggpubr))
  suppressMessages(library(psych))
  suppressMessages(library(tidyverse))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggbreak))
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

  # 计算每个基因每个组的均值、标准误差
  summary_stats <- PCR_long %>%
    group_by(Gene, group) %>%
    summarize(
      mean_relative = mean(Relative, na.rm = TRUE),
      sd_relative = sd(Relative, na.rm = TRUE),
      n = n(),
      se_relative = sd_relative / sqrt(n),
      upper_ci = mean_relative + se_relative,
      lower_ci = mean_relative - se_relative,
      .groups = 'drop'
    )

  # 找到最高的 upper_ci
  max_upper_ci <- max(summary_stats$upper_ci, na.rm = TRUE)

  # 计算每个目标基因的 p 值
  result_df <- data.frame(Gene = character(),
                          Comparison = character(),
                          group = character(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

  for (i in seq_along(target_genes)) {
    for (group in unique(PCR$group)) {
      if (group != control_group) {
        t_test_result <- t.test(PCR[[target_genes[i]]][PCR$group == control_group], PCR[[target_genes[i]]][PCR$group == group])
        result_df <- rbind(result_df, data.frame(Gene = target_genes[i], Comparison = paste(control_group, "vs", group), group = group, p_value = t_test_result$p.value))
      }
    }
  }

  # 合并 p 值和显著性标记数据
  p_value_data <- summary_stats %>%
    left_join(result_df, by = c("Gene", "group")) %>%
    mutate(
      label_y_p = upper_ci + 0.1,  # p 值的位置
      label_y_signif = upper_ci + 0.2,  # 显著性标记的位置
      label_p = ifelse(!is.na(p_value), paste0("p = ", format(p_value, digits = 2)), ""),
      label_signif = case_when(
        !is.na(p_value) & p_value < 0.001 ~ "***",
        !is.na(p_value) & p_value < 0.01 ~ "**",
        !is.na(p_value) & p_value < 0.05 ~ "*",
        !is.na(p_value) ~ "ns",
        TRUE ~ ""
      ),
      label_combined = ifelse(!is.na(p_value), paste(label_p, label_signif, sep = "\n"), "")
    )

  # 默认颜色设置
  default_colors <- c("#008ccc", "#c7b363","#c34a36", "#8c564b",  "#845ec2", "#d65db1", "#ff7f0e", "#2ca02c", "#d62728","#2f4858", "#00c9a7")

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
                       palette = group_colors,
                       add = "mean_sd",
                       position = position_dodge(width = 0.8),
                       xlab = "", ylab = 'Relative expression', legend = 'right',
                       ggtheme = theme_bw()) +
    scale_y_continuous(limits = c(0, max_upper_ci * 1.2), expand = c(0, 0))  # 设置 y 轴范围

  # 如果提供了 breaks 参数
  if (!is.null(breaks) && length(breaks) > 0) {
    # 确保 breaks 是一个合适的列表
    if (!is.list(breaks) || !all(sapply(breaks, function(x) length(x) == 2 && is.numeric(x)))) {
      stop("The 'breaks' parameter should be a list of numeric vectors with exactly two elements each.")
    }
    for (break_set in breaks) {
      p_basic <- p_basic + scale_y_break(breaks = break_set)
    }
  }

  # 对照组不显示的偏移调整
  dodge <- position_dodge(width = 0.8)

  # 添加 p 值到柱子上方，使用分组颜色
  p_with_p <- p_basic +
    geom_text(data = p_value_data %>% filter(group != control_group),
              aes(x = Gene, y = label_y_p, label = label_p, color = group),
              position = dodge, vjust = -0.5, size = 3, show.legend = FALSE)

  # 添加显著性标记到柱子上方，使用分组颜色
  p_with_signif <- p_basic +
    geom_text(data = p_value_data %>% filter(group != control_group),
              aes(x = Gene, y = label_y_signif, label = label_signif, color = group),
              position = dodge, vjust = -0.5, size = 3, show.legend = FALSE)

  # 打印图形
  print(p_basic)
  print(p_with_p)
  print(p_with_signif)
  #数据格式转换
  PCR_Ct <- t(PCR_Ct)
  PCR_relative <- t(PCR_relative)
  # 移除列名
  colnames(PCR_Ct) <- NULL
  colnames(PCR_relative) <- NULL
  # 获取当前日期
  current_date <- format(Sys.Date(), "%Y-%m-%d")

  # 创建带日期和内参基因名的文件名
  ct_filename <- paste0("PCR_Ct_", current_date, ".csv")
  relative_filename <- paste0("PCR_relative_", reference_gene,"_" , current_date ,".csv")
  # 创建带日期的p值文件名
  p_value_filename <- paste0("p_values_", "relative_to_",reference_gene,"_" ,current_date, ".csv")
  result_df <- result_df[,-2]
  # 保存 p 值数据到 CSV 文件
  write.csv(result_df, file = p_value_filename, row.names = F)
  # 确保列名不被写入文件中
  write.csv(PCR_Ct, file = ct_filename, row.names = T)
  write.csv(PCR_relative, file = relative_filename, row.names = T)
  return(list(PCR_Ct = PCR_Ct, PCR_relative = PCR_relative, p_basic = p_basic, p_with_p = p_with_p, p_with_signif = p_with_signif))
}
#ok
