# slingshot
{
  seuratB <- readRDS("seuratB.rds")
  curve_1 = slingshot::slingCurves(seuratB@misc$slingshot)[[1]]
  lineData = curve_1$s[curve_1$ord, 1:2]
  colnames(lineData) = c("PHATE_1", "PHATE_2")
  p <- dittoSeq::dittoDimPlot(seuratB, var="pseudoTime", reduction.use="phate",
                              size=0.2,
                              do.raster=TRUE, raster.dpi=450) +
    geom_path(aes(PHATE_1, PHATE_2), data=as.data.frame(lineData), linewidth=0.75) +
    labs(title="", x="PHATE 1", y="PHATE 2") +
    scale_color_viridis_c(name="Pseudotime") +
    theme_classic() +
    theme(
      plot.title   = element_text(color = "black"),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      axis.text    = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 14, color = "black"),
      legend.text  = element_text(size = 12, color = "black")
    ) +
    coord_fixed(ratio = 1)
  
  pdf("slingshot.pdf", width = 5, height = 4)
  print(p)
  dev.off()
}

# dynamic PAS heatmap
{
  dynamic_dir <- "Dynamic"
  
  pas_mean <- fread("mean.data.txt")
  res_trend <- fread("dynamic_PAS_lm.res")
  res_trend <- subset(res_trend, dynamic_type != "Static")
  table(res_trend$dynamic_type)
  sigModel <- setNames(res_trend[["dynamic_type"]], res_trend[["APATag"]])
  
  dataPlot = copy(pas_mean) %>%
    column_to_rownames("APATag") %>%
    as.matrix() %>%
    .[names(sigModel), ] %>%
    t() %>%
    scale() %>%
    t()
  
  method_vec <- sigModel[rownames(dataPlot)]
  method_vec <- factor(
    method_vec,
    levels = c("Linear", "Quadratic", "Linear & Quadratic")
  )
  method_col <- c(
    "Linear" = "#1f78b4",
    "Quadratic" = "#e31a1c",
    "Linear & Quadratic" = "#ff7f00"
  )
  
  row_ha <- rowAnnotation(
    Method = method_vec,
    col = list(Method = method_col),
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )
  
  colorHeatmap <- list(low = "#5E3C99", mid = "white", high = "#BA2E1F")
  p <- Heatmap(
    dataPlot,
    col=circlize::colorRamp2(c(min(dataPlot, na.rm=T), 0, max(dataPlot, na.rm=T)),
                             c(colorHeatmap[["low"]], colorHeatmap[["mid"]], colorHeatmap[["high"]])),
    cluster_rows = TRUE,
    cluster_row_slices = FALSE, 
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_names = TRUE,
    row_names_rot  = 90,
    column_title = "Dynamic PAS usage (3'UTR)",
    row_title = NULL,
    heatmap_legend_param = list(title = "PAS usage"),
    left_annotation = row_ha,
    row_split = method_vec,
    gap = unit(2, "mm") 
  )
  
  pdf("dynamic.PAS.heatmap.pdf", width = 4, height = 6)
  print(p)
  dev.off()
}
