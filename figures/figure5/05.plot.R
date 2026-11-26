# bubble plot
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
  
  pdf("bubble.plot.pdf", width = 15, height = 8)
  print(p)
  dev.off()
}
