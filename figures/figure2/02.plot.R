# number barplot
{
  dtLeadSNP = fread("resultCombined_leadSNP_cis.tsv", sep='\t', quote=FALSE)
  leadPairs = dtLeadSNP[, .(cellType, pairs)] %>%
    unique() %>%
    .[, .N, by=cellType] %>%
    setnames('N', "leadPairsNum")
  leadSNP = dtLeadSNP[, .(cellType, snps)] %>%
    unique() %>%
    .[, .N, by=cellType] %>%
    setnames('N', "leadSNPNum")
  leadTag = dtLeadSNP[, .(cellType, APATag)] %>%
    unique() %>%
    .[, .N, by=cellType] %>%
    setnames('N', "leadTagNum")
  dataPlot = sigPairs %>%
    merge.data.table(leadPairs, by="cellType") %>%
    merge.data.table(leadSNP, by="cellType") %>%
    merge.data.table(leadTag, by="cellType") %>%
    setorder(leadSNPNum) %>%
    .[, cellType:=factor(cellType, levels=cellType)]

  celltype_colors <- as.data.frame(celltype_colors)
  celltype_colors$cellType <- rownames(celltype_colors)
  dataPlot <- merge(dataPlot, celltype_colors, by = "cellType", all.x = TRUE, sort = FALSE)
  dataPlot$Color <- dataPlot$celltype_colors
  dataPlot <- dataPlot %>%
    arrange(leadTagNum) %>%
    mutate(cellType = factor(cellType, levels = cellType))
  
  p1 <- ggplot(dataPlot, aes(x = cellType, y = leadTagNum)) +
    geom_bar(stat = "identity", fill = dataPlot$Color) +
    labs(title = "Number of APAs (3'UTR)", x = "", y = "") +
    scale_y_reverse(expand = c(0, 0), limits = c(4200, 0), breaks = seq(0, 4000, by = 1000)) +
    scale_x_discrete(position = "top", expand = expansion(add = c(1, 1))) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),  # 标题水平居中
      axis.text = element_text(size = 14, color = 'black'),
      axis.title = element_text(size = 14, color = 'black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank()
    )
  print(p1)
  
  p2 <- ggplot(dataPlot, aes(x = cellType, y = leadSNPNum)) +
    geom_bar(stat = "identity", fill = dataPlot$Color) +
    labs(title = "Number of 3'aSNP", x = "", y = "") +
    scale_x_discrete(expand = expansion(add = c(1, 1))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 6400), breaks = seq(0, 6000, by = 2000)) +
    theme_minimal() +
    coord_flip() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 14, color = 'black'),
      axis.title = element_text(size = 14, color = 'black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank()
    )
  print(p2)
  
  pdf("PAS.num.pdf", width = 4, height = 8)
  print(p1)
  dev.off()
  
  pdf("SNP.num.pdf", width = 4, height = 8)
  print(p2)
  dev.off()
}
