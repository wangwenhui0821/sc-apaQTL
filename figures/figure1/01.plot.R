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

# Nucleotide proportion
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
  
  pdf("Nucleotide.proportion.pdf", width = 5, height = 6)
  print(p)
  dev.off()
}

# example pas
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
  pdf("example.pas.pdf", width = 5, height = 6)
  print(p)
  dev.off()
}
