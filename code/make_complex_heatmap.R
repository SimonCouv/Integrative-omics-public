make_complex_heatmap <- function(
  mat, # pw correlations, list of matrices n_clust x n_clust
  pheno_l,  # cluster-pheno associations, list of matrices n_clust x n_phenotypes
  annotation_df = NULL,
  annotationdf_colors_list = NULL, # named (=df names) list of named (=single annotation levels/extremes) vectors
  supercluster = "calc",
  deepSplit = 4,
  minClusterSize = 3,
  col_labels=colnames(mat),
  show_annotation_legends = NULL,
  row_labels=rownames(mat),
  core_col=NULL,
  core_heatmap_legend_title = "heatmap",
  pheno_p_cutoff = 0.05,
  supercluster_pheno_f,
  km=1,
  clustersizes = NULL,
  hc_method = "average",
  y,  # outcomes, n_subject x n_phenotypes
  me   # cluster PC1s, n_subject x n_clust
){
  
  library(WGCNA)

  
  # protect against drop to vector when subsetting and add index
  if (!is.null(annotation_df)){
    annotation_df <- as_tibble(annotation_df)
    annotation_df$feature <- row_labels
  }
  
  
  # filter out grey module
  # browser()
  grey_index <- grep(row_labels, pattern = "grey")
  if (length(grey_index)>0){
    pheno_l <- pheno_l %>% map(~.x[-grey_index,])
    mat <- mat[-grey_index, -grey_index]
    col_labels <- col_labels[-grey_index]
    row_labels <- row_labels[-grey_index]
    me <- me[,-grey_index]
    if (!is.null(annotation_df))
      annotation_df <- annotation_df[-grey_index,]
    if (!is.null(clustersizes))
      clustersizes <- clustersizes[-grey_index]
  }
  
  
  ### cluster rows of matrix
  hdist <- dist(mat)
  hc <- hclust(hdist, method = hc_method)
  dend <- as.dendrogram(hc)
  
  ### cluster columns of matrix
  hdist2 <- dist(t(mat))
  hc2 <- hclust(hdist2, method = hc_method)
  dend2 <- as.dendrogram(hc2)
  # phenobar <- make_phenobar(pheno_p= pheno_rowmod_cor_l$pval_m, 
  #                           pheno_est=pheno_rowmod_cor_l$est_m,
  #                           p_cutoff = pheno_p_cutoff)
  

  
  
  # default color definition for annotation_df
  if (is.null(annotationdf_colors_list)){
    annotationdf_colors_list <- map(annotation_df, ~define_colors(.x, breaks = c(0,1)))
  }

    
  
  # define supercluster
  if (supercluster=="calc"){
    dynamic_mods <- dynamicTreeCut::cutreeHybrid(dendro =hc, 
                                                 distM = as.matrix(hdist),
                                                 deepSplit = deepSplit, pamRespectsDendro = TRUE,
                                                 minClusterSize = minClusterSize)
    # col_modality_color <- labels2colors(dynamic_mods$labels)
    # dend_colored <- dendextend::color_branches(dend, clusters = dynamic_mods$labels[hc2$order])
    # annotation_df <- annotation_df[hc$order,]
    
    # rename superclusters in order of their appearance in the heatmap
    old_names_heatmap_order <- unique(dynamic_mods$labels[hc2$order])
    # remap <- data.frame(old_name_by_order = heatmap_supercluster_order,
    #                     order = 1:n_distinct(dynamic_mods$labels))
    supercluster_labels <- c()
    for (i in seq_along(dynamic_mods$labels)) supercluster_labels[i] <- which(old_names_heatmap_order == dynamic_mods$labels[i])
    #test: success
    # supercluster_labels[hc$order]
    # browser()
    # supercluster_labels_long <- paste0("cluster_", supercluster_labels)
    
    # supercluster PC1 representativeness (var explained)
    supercluster_me <- moduleEigengenes(
      me, 
      colors =  supercluster_labels, 
      nPC=10, impute = FALSE,
      excludeGrey = FALSE, subHubs = FALSE, 
      scale = FALSE)
    
    supercluster_rep <- supercluster_me[['varExplained']] %>% 
      .[1,] %>% 
      t(.) %>% 
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(supercluster = str_replace(rowname, "X", "")) %>%
      dplyr::select(`supercluster PC1 var expl.` = `1`,
                    supercluster)
    
    # add to annotation_df
    annotation_df$supercluster = paste0(supercluster_labels)
    annotation_df <- annotation_df %>%
      left_join(supercluster_rep)
    # annotationdf_colors_list$supercluster <- define_colors(annotation_df$supercluster, breaks = c(0,1)) %>% .[order(names(.))]
    # annotationdf_colors_list[['supercluster']] <- define_colors(annotation_df$supercluster[order(hc$order)], breaks = c(0,1))
    # annotationdf_colors_list[['supercluster']] <- structure(ggplot_qual_colors(n_distinct(dynamic_mods$labels)), names=unique(dynamic_mods$labels))
    annotationdf_colors_list <- map(annotation_df, ~define_colors(.x, breaks = c(0,1)))
    annotationdf_colors_list$supercluster <- define_colors(annotation_df$supercluster, breaks = c(0,1)) %>% .[order(names(.))]
    annotation_df <- annotation_df %>% dplyr::select(supercluster, `supercluster PC1 var expl.`, everything())
    if (is.null(show_annotation_legends))
      show_annotation_legends <- c(TRUE, FALSE, rep(TRUE, ncol(annotation_df)-2))
    else if (show_annotation_legends == "all")
      show_annotation_legends <- rep(TRUE, ncol(annotation_df))
    
    # make top_annotation
    top_annotation <- HeatmapAnnotation(
      supercluster=annotation_df$supercluster,
      col = list(supercluster=structure(ggplot_qual_colors(n_distinct(supercluster_labels)), names=unique(paste0(supercluster_labels)))),
      show_legend = FALSE
      # show_annotation_name = c(supercluster=FALSE)
      )
    
    # browser()
    # make supercluster pheno bar
    supercluster_me_pheno_l <- do.call(supercluster_pheno_f, list(supercluster_me$eigengenes, y)) %>%
      .[1:3] %>%
      map(~.x %>% 
            as_tibble() %>% 
            dplyr::select(CVD =  matches("lrm_cvd_[^(sa)]"),
                          `CVD, age+sex adj` = matches("lrm_cvd_sa_"),
                          `CVD, age+sex+statin adj` = matches("lrm_cvd_sas_")) %>%
            as.matrix(.))
    # duplicate according to cluster-supercluster membership
    supercluster_me_pheno_l <- supercluster_me_pheno_l %>% map(~.x[supercluster_labels,])
    
    
    super_pheno_df <- make_phenobar(pheno_p = supercluster_me_pheno_l$pval_m, 
                              pheno_est = supercluster_me_pheno_l$est_m,
                              p_cutoff=pheno_p_cutoff) %>% 
      as_tibble() %>%
      mutate_all(~recode(.,"grey"="non-signif.", "darkred" = "positive", "darkblue"="negative"))
    
    super_pheno_colors <- c("non-signif."="grey",  "positive"="darkred", "negative" = "darkgreen")
    
    super_pheno_mat <- as.matrix(super_pheno_df)
    rownames(super_pheno_mat) <- row_labels
    super_pheno_hm <- Heatmap(matrix=super_pheno_mat, 
                        col =super_pheno_colors, 
                        row_labels = row_labels, 
                        column_labels = colnames(super_pheno_mat),
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        name = "CVD association_superclust",
                        show_heatmap_legend = FALSE)
    
  }
  
  # define phenotype annotations
  pheno_df <- make_phenobar(pheno_p = pheno_l$pval_m, pheno_est = pheno_l$est_m, p_cutoff=pheno_p_cutoff) %>% as_tibble() %>%
    mutate_all(~recode(.,"grey"="non-signif.", "darkred" = "positive", "darkblue"="negative"))
  pheno_colors <- c("non-signif."="grey",  "positive"="darkred", "negative" = "darkgreen")
  pheno_mat <- as.matrix(pheno_df)
  rownames(pheno_mat) <- row_labels
  pheno_hm <- Heatmap(matrix=pheno_mat, 
                      col =pheno_colors, 
                      row_labels = row_labels, 
                      column_labels = colnames(pheno_mat),
                      cluster_rows = FALSE, 
                      cluster_columns = FALSE,
                      name = "CVD association")

  # browser()
  if (is.null(core_col)){
    core_col <- define_colors(as.vector(mat), breaks = c(-1,0,1), quant_colors = c("blue", "white", "red"))
  }
  
  hm_args <- list(
    matrix = mat,
    cluster_columns = rev(dend2),
    cluster_rows = dend,
    col=core_col,
    row_labels = row_labels,
    column_labels = col_labels,
    name = core_heatmap_legend_title,
    column_names_rot = 55,
    column_names_side = "top",
    column_dend_side = "bottom",
    show_row_dend=FALSE
    # km=km
  )

  # browser()
  if (! is.null(annotation_df)){
    supercluster_id_hm <- HeatmapAnnotation(n = anno_text(annotation_df$supercluster), which = "row", show_annotation_name = TRUE)
    cluster_names_hm <- HeatmapAnnotation(n = anno_text(row_labels), which = "row", show_annotation_name = TRUE)
    anno_hm <- HeatmapAnnotation(df = as.data.frame(dplyr::select(annotation_df, -feature), stringsAsFactors=FALSE),
                                 which = "row", col=annotationdf_colors_list, show_legend = show_annotation_legends)
    # anno_hm_ext <-  HeatmapAnnotation(df = as.data.frame(annotation_df, stringsAsFactors=FALSE),
    #                                   size = anno_text(clustersizes),
    #                                   which = "row", col=annotationdf_colors_list, 
    #                                   show_legend = show_annotation_legends)
    if (!is.null(clustersizes)){
      size_hm <- HeatmapAnnotation(size =  anno_barplot(clustersizes), which = "row", show_annotation_name = c(size=TRUE))
      hm_list <- super_pheno_hm + pheno_hm + do.call(Heatmap, c(hm_args, list(top_annotation=top_annotation))) + anno_hm + size_hm + cluster_names_hm
      full_hm <- draw(hm_list, ht_gap = unit(c(2,1,4,1,1), "mm"), main_heatmap = 3)
    } else {
      hm_list <- super_pheno_hm + pheno_hm + do.call(Heatmap, c(hm_args, list(top_annotation=top_annotation))) + anno_hm + cluster_names_hm
      full_hm <- draw(hm_list, ht_gap = unit(c(2,1,4,1), "mm"), main_heatmap = 3)
    }
    
    # visualise cluster IDs, to map returned supercluster_me eigengene order
    # hm_list <- super_pheno_hm + pheno_hm + do.call(Heatmap, c(hm_args, list(top_annotation=top_annotation))) + supercluster_id_hm + anno_hm + size_hm + cluster_names_hm
    # full_hm <- draw(hm_list, ht_gap = unit(c(2,1,4,1,1,1), "mm"), main_heatmap = 3)

    
    
    # core_hm <- do.call(Heatmap, c(hm_args, list(right_annotation=anno_hm_ext, top_annotation=top_annotation)))  # preserves feature names (draw(hm+anno_hm) does not)
  } else {
    full_hm <- do.call(Heatmap, hm_args)
  }
  
  return(list(hm=full_hm, row_order = hc$order, supercluster_me = supercluster_me, supercluster_colorcode = annotationdf_colors_list$supercluster))
  # return(list(hm=full_hm, row_order = hc$order, supercluster_me = ret_supercluster_me, supercluster_colorcode = annotationdf_colors_list$supercluster))
  
  # DO NOT USE, labeling incorrect
  # ret_supercluster_me <- structure(supercluster_me$eigengenes, names = names(annotationdf_colors_list$supercluster))
  # ret_supercluster_me <- ret_supercluster_me[,rev(unique(as.character(dynamic_mods$labels)[hc$order]))]   # ! reordering is based on (character) column names, not indices!, uncomment visualise cluster IDs to comfirm when in doubt 
  # #order supercluster_me columns to match order (left-right) in which they occur in the heatmmap
  # return(list(hm=full_hm, row_order = hc$order, supercluster_me = ret_supercluster_me, supercluster_colorcode = annotationdf_colors_list$supercluster))
}
