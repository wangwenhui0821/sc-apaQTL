# mash sign
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
  
  pdf("mash.sign.pdf", width = 8, height = 8)
  print(p)
  dev.off()
}

# mash magnitude
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

  pdf("mash.magnitude.pdf", width = 8, height = 8)
  print(p)
  dev.off()
}
