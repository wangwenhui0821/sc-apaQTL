# coloc & smr
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
  
  pdf("coloc.smr.pdf", width = 17, height = 5)
  print(p)
  dev.off()
}
