make_pedersen_heatmap <- function(
  filename,
  mod_mod_cor_l,
  pheno_rowmod_cor_l,
  colsidecolors=NULL,
  dev_height=10,
  dev_width=10, 
  dendro_show="both",
  symm = FALSE,
  keyvalname,
  revC,
  pheno_p_cutoff = 0.05,
  labRow=NULL,
  labCol=NULL,
  dev="pdf",
  margins = c(23, 33),
  row_modality_color = NULL,
  col_label_color = NULL,
  modality_bar=NULL,
  deepSplit=4,
  minClusterSize=3
){
  
  if (any(rownames(mod_mod_cor_l$est_m) != rownames(pheno_rowmod_cor_l$est_m))){
    stop("rows of mod_mod_cor_l do not match rows of pheno_rowmod_cor_l")
  }
  

  
  ### cluster rows of matrix
  hdist <- dist(mod_mod_cor_l$est_m)
  hc <- hclust(hdist, method = "average")
  dend = as.dendrogram(hc)
  
  ### cluster columns of matrix
  hdist2 <- dist(t(mod_mod_cor_l$est_m))
  hc2 <- hclust(hdist2, method = "average")
  dend2 =as.dendrogram(hc2)
  phenobar <- make_phenobar(pheno_p= pheno_rowmod_cor_l$pval_m, 
                            pheno_est=pheno_rowmod_cor_l$est_m,
                            p_cutoff = pheno_p_cutoff)
  plotmat <- mod_mod_cor_l$est_m
  stars <- mod_mod_cor_l$p_star_m
  
  # define superclusters
  
  if (col_label_color == "supercluster"){
    dynamic_mods <- dynamicTreeCut::cutreeHybrid(dendro =hc2, 
                                                 distM = as.matrix(hdist2),
                                                 deepSplit = deepSplit, pamRespectsDendro = TRUE,
                                                 minClusterSize = minClusterSize)
    col_modality_color <- labels2colors(dynamic_mods$labels)
  }
  
  ### make heatmap with all the selected KEGG modules
  switch(dev,
         "pdf" = pdf (filename, height = dev_height, width = dev_width),
         "png" =  png (filename, height = dev_height, width = dev_width))
  # pdf (filename, height = dev_height, width = dev_width)
  hm_args = list(plotmat,
                 Rowv = dend, 
                 Colv = dend2, 
                 dendrogram = dendro_show,
                 col = bluered (100),
                 symbreaks = T,
                 revC = FALSE,
                 key = T,
                 symkey = T,
                 symm = symm,
                 keysize = 1,
                 KeyValueName = keyvalname,
                 lhei = c(0.6, 4),
                 cellnote = stars, 
                 notecol = "black",
                 notecex = 1.1,
                 trace = "none",
                 labRow = labRow,
                 labCol = labCol,
                 RowSideColors = t(phenobar[rownames(plotmat),]), 
                 side.height.fraction = 0.35,    
                 NumColSideColors = dim(phenobar)[2], 
                 margins = margins,
                 cexRow = 1,
                 cexCol = 1,
                 modality_bar=modality_bar
  )
  
  if (!is.null(colsidecolors)){
    hm_args <- c(hm_args, list(ColSideColors=colsidecolors))
  }
  if (!is.null(row_modality_color)){
    hm_args <- c(hm_args, list(row_modality_color =  row_modality_color[hc$order]))
  }
  if (!is.null(col_modality_color)){
    hm_args <- c(hm_args, list(col_modality_color =  col_modality_color[hc2$order]))
  }
  do.call(heatmap.3, hm_args)
  
  dev.off ()
  # return(list(col_hc = hc2, col_dist = hdist2))
}
