celltypeV <- c("CD4_TCM", "CD4_Naive", "CD8_TEM", "NK", "CD8_Naive", "CD4_TEM", "CD8_TCM", "T_Reg", "gdT",
               "MAIT", "CD4_CTL", "NK_R", "B_IN", "B_Mem", "DC", "Mono_C", "Mono_NC")
celltype_colors <- c(CD4_Naive = '#5E2750', CD4_TCM = '#7A1F5C', CD4_TEM = '#BD5D7C', CD4_CTL = '#FF9A9B',
                     CD8_Naive = '#1A4174', CD8_TCM = '#456FA5', CD8_TEM = '#6F9CD6', T_Reg = '#2F5256',
                     gdT = '#3B8E8B', MAIT = '#7DCCA6',  NK = '#1C7C1C', NK_R = '#33B733',B_Mem = '#A77533',
                     B_IN = '#D3B04F', Mono_C = '#D84F19', Mono_NC = '#F88A31', DC = '#C1272D')

diseaseV <- c("MS.GCST003566", "PSO.GCST90243956", "SLE.GCST003156", "T1DM.GCST90014023", 
              "CD.ibdgc", "IBD.ibdgc", "RA.jenger", "Celiac.phecode557")
disease_colors <- c(CD = "#ED6A5A", Celiac = "#F5A623", IBD = "#A56C85", MS = "#DA8A67", 
                    PSO = "#8EB78F", RA = "#4FA092", SLE = "#6DAEDB", T1DM = "#4E6C91")

# Figure 1B
{
  hg38 <- readDNAStringSet("Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  names(hg38) <- sub(" .*", "", names(hg38)) 
  
  all_pas <- fread("metadata_intronic_and_tandem_revised.txt")
  all_pas <- all_pas[, c("peak", "strand", "seqnames","start", "end", "Intron.exon_location", "hg38_Position", "peak_short_name")]
  all_pas[, position := as.integer(tstrsplit(peak_short_name, "_")[[2]])]
  
  pos <- all_pas$position
  chr <- all_pas$seqnames
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = pmax(pos-50, 1), end = pos+50),
                strand = all_pas$strand)
  seqlevels(gr) <- sub("^chr", "", seqlevels(gr))
  seqnames(gr) <- sub("^chr", "", seqnames(gr))
  
  seqs <- getSeq(hg38, gr)
  all_pas$read <- as.character(seqs)
  
  get_base_freq <- function(df) {
    seq_matrix <- do.call(rbind, strsplit(df$read, split = ""))
    seq_df <- as.data.frame(seq_matrix, stringsAsFactors = FALSE)
    colnames(seq_df) <- -50:50
    base_freq <- apply(seq_df, 2, function(x) {
      table(factor(x, levels = c("A","C","G","T"))) / length(x)
    })
    base_freq_df <- as.data.frame(t(base_freq))
    base_freq_df$pos <- as.numeric(rownames(base_freq_df))
    base_freq_long <- base_freq_df %>%
      pivot_longer(cols = c("A","C","G","T"), names_to = "base", values_to = "freq")
    base_freq_long$PAS_group <- unique(df$Intron.exon_location)
    return(base_freq_long)
  }
  
  utr_freq <- get_base_freq(subset(all_pas, Intron.exon_location == "3'_most_exon"))
  intron_freq <- get_base_freq(subset(all_pas, Intron.exon_location == "Intron"))
  utr_freq$PAS_group <- "3'UTR-APA"
  intron_freq$PAS_group <- "IPA"
  freq_all <- rbind(utr_freq, intron_freq)

  p <- ggplot(freq_all, aes(x = pos, y = freq, color = base)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.2, size = 0.7) +
    scale_color_manual(values = c("A"="#00BD2F", "C"="#3600E7", "G"="#F9A932", "T"="#E01412")) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_wrap(~ PAS_group, ncol = 1, scales = "free_y", strip.position = "top") +
    labs(x = "Distance to polyA sites (nt)", 
         y = "Nucleotide proportion", 
         color = "Base") + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.ticks = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_blank(),
      strip.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "grey90", color = "black")
    )
  
  pdf("Fig1B.pdf", width = 5, height = 6)
  print(p)
  dev.off()
  
}

# Figure 1F
{
  pas_info <- fread("metadata_intronic_and_tandem_revised.txt")
  pas_info <- pas_info[, c("seqnames", "hg38_Position", "peak_short_name", "Ensemble_ID")]
  pas_info <- pas_info[, APATag := paste(Ensemble_ID, seqnames, hg38_Position, sep = ".")]
  head(pas_info)
  
  apa_mat <- readRDS("PAScount.rds")
  meta <- fread("metaData.list")
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$cell
  
  common_cells <- intersect(colnames(apa_mat), rownames(meta))
  apa_mat <- apa_mat[, common_cells]
  meta <- meta[common_cells, ]
  
  seurat_obj <- CreateSeuratObject(counts = apa_mat, meta.data = meta)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj10 <- RunUMAP(seurat_obj, dims = 1:10)
  
  p <- DimPlot(seurat_obj10, group.by = "cellTypePaper", label = TRUE, cols = celltype_colors, repel = TRUE, 
               label.size = 6) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      axis.ticks = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.position = "none"
    )
  pdf("Fig1F.pdf", width = 5, height = 6)
  print(p)
  dev.off()
}

# Figure 2A
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
  
  pdf("Fig2A.PAS.pdf", width = 4, height = 8)
  print(p1)
  dev.off()
  
  pdf("Fig2A.SNP.pdf", width = 4, height = 8)
  print(p2)
  dev.off()
}

# Figure 3A
{
  out <- readRDS("fastqtl_to_mash_output/sample.mash.rds")
  maxb <- out$strong.b
  maxz <- out$strong.z
  out <- readRDS("mashr_flashr_workflow_output/sample.mash.EZ.posterior.rds")
  pm.mash <- out$PosteriorMean
  lfsr.all <- out$lfsr
  standard.error <- maxb/maxz
  pm.mash.beta <- pm.mash*standard.error
  
  thresh=0.05
  pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
  lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]
  shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
  colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(maxz)
  for(i in 1:ncol(lfsr.mash)){
    for(j in 1:ncol(lfsr.mash)){
      sig.row=which(lfsr.mash[,i]<thresh)
      sig.col=which(lfsr.mash[,j]<thresh)
      a=(union(sig.row,sig.col))
      quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
      shared.fold.size[i,j] = mean(quotient > 0, na.rm = TRUE)
    }
  }
  
  rownames(shared.fold.size) <- sub("\\.tsv$", "", rownames(shared.fold.size))
  colnames(shared.fold.size) <- sub("\\.tsv$", "", colnames(shared.fold.size))
  lat <- shared.fold.size
  
  celltypes_order <- c("Mono_NC", "Mono_C", "DC", "B_Mem", "B_IN", "NK_R", "CD4_CTL", "MAIT",
                       "gdT", "T_Reg", "CD8_TCM", "CD4_TEM", "CD8_Naive", "NK", "CD8_TEM",
                       "CD4_Naive", "CD4_TCM")
  celltypes_order <- rev(celltypes_order)
  lat <- lat[celltypes_order, celltypes_order]
  lat[lower.tri(lat)] <- NA
  
  clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                 "#E0F3F8","#91BFDB","#4575B4")))(100)
  
  celltype_colors <- celltype_colors[celltypes_order]
  
  annotation_col <- data.frame(cellType = factor(colnames(lat), levels = celltypes_order))
  rownames(annotation_col) <- colnames(lat)
  
  annotation_colors_col <- list()
  annotation_colors_col[[1]] <- celltype_colors[celltypes_order]
  names(annotation_colors_col) <- NULL
  
  p <- pheatmap(lat,
                color = clrs,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                display_numbers = FALSE,
                fontsize_row = 12,
                fontsize_col = 10,
                angle_col = 45,
                na_col = "white",
                border_color = "white",
                cellwidth = 16,
                cellheight = 16, 
                labels_col = rep("", ncol(lat)),  
                annotation_row = annotation_col,
                annotation_col = annotation_col,
                annotation_colors = list(cellType = celltype_colors), 
                annotation_names_row = FALSE,
                annotation_names_col = FALSE,
                annotation_legend = FALSE
  )
  
  pdf("Fig.3A.pdf", width = 8, height = 8)
  print(p)
  dev.off()
}

# Figure 3B
{
  out <- readRDS("fastqtl_to_mash_output/sample.mash.rds")
  maxb <- out$strong.b
  maxz <- out$strong.z
  out <- readRDS("mashr_flashr_workflow_output/sample.mash.EZ.posterior.rds")
  pm.mash <- out$PosteriorMean
  lfsr.all <- out$lfsr
  standard.error <- maxb/maxz
  pm.mash.beta <- pm.mash*standard.error
  
  thresh <- 0.05
  pm.mash.beta <- pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
  lfsr.mash <- lfsr.all[rowSums(lfsr.all<0.05)>0,]
  shared.fold.size <- matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
  colnames(shared.fold.size) <- rownames(shared.fold.size) <- colnames(maxz)
  for (i in 1:ncol(lfsr.mash))
    for (j in 1:ncol(lfsr.mash)) {
      sig.row=which(lfsr.mash[,i]<thresh)
      sig.col=which(lfsr.mash[,j]<thresh)
      a=(union(sig.row,sig.col))
      quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
      shared.fold.size[i,j] = mean(quotient > 0.5 & quotient < 2, na.rm = TRUE)
    }
  
  rownames(shared.fold.size) <- sub("\\.tsv$", "", rownames(shared.fold.size))
  colnames(shared.fold.size) <- sub("\\.tsv$", "", colnames(shared.fold.size))
  lat <- shared.fold.size
  
  celltypes_order <- c("Mono_NC", "Mono_C", "DC", "B_Mem", "B_IN", "NK_R", "CD4_CTL", "MAIT",
                       "gdT", "T_Reg", "CD8_TCM", "CD4_TEM", "CD8_Naive", "NK", "CD8_TEM",
                       "CD4_Naive", "CD4_TCM")
  celltypes_order <- rev(celltypes_order)
  lat <- lat[celltypes_order, celltypes_order]
  lat[lower.tri(lat)] <- NA
  
  clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                 "#E0F3F8","#91BFDB","#4575B4")))(100)
  
  celltype_colors <- celltype_colors[celltypes_order]
  
  annotation_col <- data.frame(cellType = factor(colnames(lat), levels = celltypes_order))
  rownames(annotation_col) <- colnames(lat)
  
  annotation_colors_col <- list()
  annotation_colors_col[[1]] <- celltype_colors[celltypes_order]
  names(annotation_colors_col) <- NULL
  
  p <- pheatmap(lat,
                color = clrs,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                display_numbers = FALSE,
                fontsize_row = 12,
                fontsize_col = 10,
                angle_col = 45,
                na_col = "white",
                border_color = "white",
                cellwidth = 16,
                cellheight = 16, 
                labels_col = rep("", ncol(lat)),  
                annotation_row = annotation_col,
                annotation_col = annotation_col,
                annotation_colors = list(cellType = celltype_colors), 
                annotation_names_row = FALSE,
                annotation_names_col = FALSE,
                annotation_legend = FALSE
  )

  pdf("Fig3B.pdf", width = 8, height = 8)
  print(p)
  dev.off()
}

# Figure 4A
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
  
  pdf("Fig4A.pdf", width = 5, height = 4)
  print(p)
  dev.off()
}

# Figure 4H
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
  
  pdf("/data1/wangwenhui/WYH/Figure.V2/Fig.5/04.dyn.PAS.pdf", width = 4, height = 6)
  print(p)
  dev.off()
}

# Figure 5E
{
  sig_coloc_pt <- "coloc.res.all.txt"
  coloc <- fread(sig_coloc_pt)
  coloc$pairs <- paste0(coloc$lead_snp, "@", coloc$gene)
  coloc_sig <- coloc[PP.H4 > 0.75]
  coloc_sig <- coloc_sig[disease %in% diseaseV & celltype %in% cellTypeV]
  coloc_sig$disease <- sub("\\..*$", "", coloc_sig$disease)
  setnames(coloc_sig, old = "gene", new = "APATag")
  
  diseaseV <- c("CD", "Celiac", "IBD", "MS", "PSO", "RA", "SLE", "T1DM")
  
  apa_info_pt <- "tagRef.txt"
  apa_info <- fread(apa_info_pt)
  apa_info <- apa_info[, c("idPAS", "Ensemble_ID", "Gene_Symbol")]
  setnames(apa_info, old = "idPAS", new = "APATag")
  apa_info[, idPAS := paste0(Gene_Symbol, sub("^[^.]+", "", APATag))]
  
  coloc_sig <- merge(coloc_sig, apa_info, by = "APATag", all.x = TRUE)
  head(coloc_sig)
  head(coloc_sig)
  length(unique(coloc_sig$APATag))
  length(unique(coloc_sig$Gene_Symbol))
 
  coloc_counts <- coloc_sig %>%
    group_by(Gene_Symbol, celltype, disease) %>%
    summarise(n_coloc = n(), .groups = "drop")
  
  coloc_wide <- coloc_counts %>%
    tidyr::pivot_wider(
      names_from = disease,
      values_from = n_coloc,
      values_fill = 0
    )
  coloc_wide <- coloc_wide %>%
    rowwise() %>%
    mutate(total_coloc = sum(c_across(where(is.numeric)))) %>%
    ungroup()
  coloc_wide <- coloc_wide %>%
    mutate(y = as.numeric(factor(celltype, levels = rev(cellTypeV))),
           x = as.numeric(factor(Gene_Symbol, levels = rev(unique(coloc_sig$Gene_Symbol)))))
  
  p <- ggplot() +
    geom_scatterpie(
      aes(x = x, y = y, r = 0.5),
      data = coloc_wide,
      cols = diseaseV,
      color = NA
    ) +
    scale_fill_manual(name = "Disease", values = disease_colors) +
    scale_y_continuous(breaks = seq_along(cellTypeV),
                       labels = rev(cellTypeV)) +
    scale_x_continuous(breaks = 1:length(unique(coloc_sig$Gene_Symbol)),
                       labels = rev(unique(coloc_sig$Gene_Symbol)),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = NULL, y = NULL, title = "3'aQTL") +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", size = 0.5),
      axis.text.x = element_text(size = 14, color = "black", angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black", vjust = 0.5, hjust = 1),
      axis.title = element_text(size = 16),
      axis.ticks = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      plot.title = element_text(size = 16, color = "black", hjust = 0.5),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(color = "black", size = 12),
      legend.text  = element_text(color = "black", size = 12)
    ) +
    coord_fixed() +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )
  
  pdf("Fig5E.pdf", width = 15, height = 8)
  print(p)
  dev.off()
}

# Figure 6A
{
  apa_info_pt <- "tagRef.tsv"
  apa_info <- fread(apa_info_pt)
  apa_info <- apa_info[, c("idPAS", "Ensemble_ID", "Gene_Symbol")]
  apa_info[, pos := {
    parts <- tstrsplit(idPAS, "\\.")
    paste(parts[[2]], parts[[3]], sep = ".")
  }]
  apa_info[, idNew := paste0(Gene_Symbol, ".", pos)]
  head(apa_info)
  
  sig_pt <- "COLOC.SMR.txt"
  sig <- fread(sig_pt)
  sig <- sig[disease %in% diseaseV & cellType %in% cellTypeV]
  sig[, disease := sub("\\..*$", "", disease)]
  min(sig$p_HEIDI.smr, na.rm = TRUE)
  sig[which.min(p_HEIDI.smr)]
  sig[, SMR := fifelse(!is.na(pAdjust.smr) & !is.na(p_HEIDI.smr) & pAdjust.smr < 0.05 & p_HEIDI.smr > 0.01, "sig.smr", "")]
  table(sig$SMR)
  sig <- subset(sig, COLOC == "sig.coloc" & SMR == "sig.smr")
  
  apa_count <- sig[
    , .N, by = .(Gene_Symbol, APATag)
  ]
  top_apa <- apa_count[
    , .SD[which.max(N)], by = Gene_Symbol
  ]
  table(top_apa$Gene_Symbol)
  sig_top <- sig[
    top_apa, on = .(Gene_Symbol, APATag)
  ]
  sig_apas <- unique(sig$APATag)
  unique(sig$Gene_Symbol)
  
  length(unique(sig$APATag))
  length(unique(sig$Gene_Symbol))
  
  diseaseV <- c("CD", "Celiac", "IBD", "MS", "PSO", "RA", "SLE", "T1DM")
  smr_all_pt <- "smr.sig.dat"
  smr_all <- fread(smr_all_pt)
  smr_all <- smr_all[, c("disease", "celltype", "probeID", "b_SMR")]
  colnames(smr_all) <- c("disease", "cellType", "APATag", "b_SMR")
  smr_sub <- smr_all[APATag %in% sig_apas]
  smr_sub[, disease := sub("\\..*$", "", disease)]
  
  map_dt <- apa_info[, .(idPAS, idNew)]
  setkey(map_dt, idPAS)
  setkey(smr_sub, APATag)
  
  smr_sub <- map_dt[smr_sub]
  smr_sub[, APATag := idNew]
  smr_sub[, idNew := NULL]
  head(smr_sub)

  sig_apas_new <- c("ZFP36L1.chr14.68788632", "TUFM.chr16.28842411", "SF3A3.chr1.37957513", 
                    "ORMDL3.chr17.39922275", "ERAP2.chr5.96918770", "ARPC2.chr2.218254348")
  smr_sub <- smr_sub[smr_sub$APATag %in% sig_apas_new, ]
  full_grid <- CJ(
    disease = diseaseV,
    cellType = cellTypeV,
    APATag = sig_apas_new,
    unique = TRUE
  )
  
  setkey(full_grid, disease, cellType, APATag)
  setkey(smr_sub,  disease, cellType, APATag)
  
  smr_sub_full <- smr_sub[full_grid]
  smr_sub_full[, disease := factor(disease, levels = diseaseV)]
  smr_sub_full[, cellType := factor(cellType, levels = cellTypeV)]
  smr_sub_full[, APATag := factor(APATag, levels = rev(sig_apas_new))]
  
  setorder(smr_sub_full, disease, cellType, APATag)
  smr_sub_full[, geneName := tstrsplit(APATag, "\\.", fixed = FALSE)[[1]]]
  smr_sub_full[
    , b_SMR.new := fifelse(
      b_SMR > 1,  1,
      fifelse(b_SMR < -1, -1, b_SMR)
    )
  ]
  
  message("Number of NA b_SMR: ", sum(is.na(smr_sub_full$b_SMR)))
  summary(smr_sub_full$b_SMR)
  
  cols <- brewer.pal(11, "RdBu")
  colors <- c(cols[2], cols[6], cols[10])
  values <- scales::rescale(c(-1, 0, 1))
  
  p <- ggplot(smr_sub_full, aes(x = disease, y = APATag, fill = b_SMR.new)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = colors,
      values = values,
      limits = c(-1, 1),
      breaks = seq(-1, 1, 1),
      oob = scales::squish,
      na.value = "grey90", 
      name = "b_SMR",
      guide = guide_colorbar(
        ticks = FALSE
      )
    ) +
    facet_grid(. ~ cellType) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_blank(),
      axis.ticks = element_line(size = 0.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      strip.text = element_text(size = 12),
      strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),
      panel.grid = element_blank(),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    coord_fixed(ratio = 1)
  p
  
  pdf("Fig6A.pdf", width = 17, height = 5)
  print(p)
  dev.off()
}


