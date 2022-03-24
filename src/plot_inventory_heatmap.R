get_sample_inventory <- function(inventory, unique = TRUE) {
  if (unique == TRUE) {
    sample_inventory <- inventory %>%
      dplyr::select(
        patient_id,
        sample_type,
        tumor_type,
        tumor_supersite,
        tumor_site,
        tumor_subsite,
        therapy,
        specimen_lab_medicine_drawn_date_timepoint,
        technique,
        is_site_matched
      ) %>%
      dplyr::mutate(technique_status = 1) %>%
      distinct(
        patient_id,
        sample_type,
        tumor_supersite,
        tumor_site,
        tumor_subsite,
        therapy,
        specimen_lab_medicine_drawn_date_timepoint,
        technique,
        .keep_all = TRUE
      ) %>%
      spread(technique, technique_status, fill = 0)
  } else {
    sample_inventory <- inventory %>%
      dplyr::select(
        patient_id,
        sample_type,
        tumor_type,
        tumor_supersite,
        tumor_site,
        tumor_subsite,
        therapy,
        specimen_lab_medicine_drawn_date_timepoint,
        technique,
        is_site_matched
      ) %>%
      dplyr::mutate(technique_status = 1) %>%
      group_by(
        patient_id,
        sample_type,
        tumor_type,
        tumor_supersite,
        tumor_site,
        tumor_subsite,
        therapy,
        specimen_lab_medicine_drawn_date_timepoint,
        technique,
        is_site_matched
      ) %>%
      summarise(n_samples = n()) %>%
      ungroup() %>%
      spread(technique, n_samples, fill = 0)
  }
  
}

# plot_inventory_heatmap <-
#   function(inventory_mat,
#            bottom_annotation,
#            left_annotation,
#            column_order,
#            row_split = "technique",
#            column_split = "consensus_signature",
#            column_split_shared = NULL) {
#     patient_id_col = clrs$patient_id
#     tumor_supersite_col = list("Site" = clrs$tumor_supersite[!names(clrs$tumor_supersite)%in%c("UQ","Unknown")])
#     signature_col = list("Mutational\nsignature" = clrs$consensus_signature[names(clrs$consensus_signature)!="NA"])
#     chemo_intent_col = list("Setting" = clrs$chemo_intent)
#     pathology_stage_col = list("Stage" = clrs$pathology_stage)
#     age_col_fun = circlize::colorRamp2(seq(30, 90, 20), brewer.pal(4, "YlOrRd"))
#     
#     row_split = left_annotation[[row_split]]
#     column_split = bottom_annotation[[column_split]]
#     
#     column_gap = unit(4, "mm")
#     
#     signature_colors = plyr::mapvalues(
#       bottom_annotation$consensus_signature,
#       from = names(signature_col$`Mutational\nsignature`),
#       to = signature_col$`Mutational\nsignature`
#     ) %>% as.character
#     site_colors = plyr::mapvalues(
#       left_annotation$tumor_supersite,
#       from = names(tumor_supersite_col$`Site`),
#       to = tumor_supersite_col$`Site`
#     ) %>% as.character
#     
#     column_title = "MSK SPECTRUM"
#     
#     top_df = as_tibble(colSums(inventory_mat), rownames = "patient_id") %>%
#       column_to_rownames(var = "patient_id") %>%
#       as.matrix
#     
#     top_annotation <- ComplexHeatmap::columnAnnotation(
#       "# samples\nper patient" = anno_barplot(
#         top_df,
#         gp = gpar(col = signature_colors, fill = signature_colors),
#         # axis_param = list(side = "right"),
#         border = FALSE
#       ),
#       annotation_name_gp = list(fontsize = 10),
#       annotation_width = unit(c(1, 4), "cm")
#     )
#     
#     bottom_df = bottom_annotation %>%
#       dplyr::select(
#         patient_age,
#         gyn_diagnosis_chemo_intent_description,
#         gyn_diagnosis_figo_stage,
#         consensus_signature
#       ) %>%
#       as.data.frame
#     colnames(bottom_df) <- c("Age", "Setting", "Stage", "Signature")
#     
#     bottom_annotation <- ComplexHeatmap::columnAnnotation(
#       df = bottom_df,
#       col = list(
#         "Age" = age_col_fun,
#         "Setting" = clrs$chemo_intent[names(clrs$chemo_intent) != "Salvage"],
#         "Stage" = clrs$pathology_stage,
#         "Signature" = clrs$consensus_signature[names(clrs$consensus_signature) != "NA"]
#       ),
#       show_annotation_name = TRUE,
#       annotation_name_side = "left",
#       annotation_name_gp = list(fontsize = 10),
#       # annotation_legend_param = list(
#       #     signature = list(title = "Signature")),
#       simple_anno_size = unit(0.35, "cm")
#       # na_col = "white", gp = gpar(col = "white")
#     )
#     
#     left_df = left_annotation %>%
#       dplyr::select(tumor_supersite) %>%
#       as.data.frame
#     colnames(left_df) <- c("Site")
#     
#     left_annotation <- ComplexHeatmap::rowAnnotation(
#       df = left_df,
#       col = list("Site" = clrs$tumor_supersite[!names(clrs$tumor_supersite) %in% c("UQ","Unknown")]),
#       show_annotation_name = TRUE,
#       annotation_name_side = "bottom",
#       annotation_name_gp = list(fontsize = 10),
#       simple_anno_size = unit(0.35, "cm")
#     )
#     
#     site_count = rowSums(inventory_mat)
#     
#     right_annotation <- ComplexHeatmap::rowAnnotation(
#       "test" = anno_text(
#         site_count, 
#         gp = gpar(fontsize = 8)),
#       "# samples\nper site" = anno_barplot(
#         site_count,
#         gp = gpar(col = site_colors, fill = site_colors),
#         border = FALSE),
#       annotation_name_gp = list(fontsize = 10)
#     )
#     
#     at = seq(0, 2, by = 1)
#     col_fun = circlize::colorRamp2(c(0, 1, 2),
#                                    c("grey90", "grey70", "grey50"))
#     
#     heatmap_legend_param = list(
#       at = at,
#       color_bar = "discrete",
#       legend_gp = gpar(fill = col_fun)
#     )
#     
#     ht <- ComplexHeatmap::Heatmap(
#       inventory_mat,
#       name = "# samples",
#       heatmap_legend_param = heatmap_legend_param,
#       col = col_fun,
#       show_row_dend = FALSE,
#       show_column_dend = FALSE,
#       show_row_names = FALSE,
#       show_column_names = TRUE,
#       border = FALSE,
#       # top_annotation = top_annotation,
#       # bottom_annotation = bottom_annotation,
#       left_annotation = left_annotation,
#       right_annotation = right_annotation,
#       row_order = rownames(inventory_mat),
#       column_title = "Patient",
#       column_title_side = "bottom",
#       row_names_side = "left",
#       column_names_side = "bottom",
#       row_split = row_split,
#       column_order = column_order,
#       cluster_row_slices = FALSE,
#       cluster_column_slices = FALSE,
#       cluster_rows = FALSE,
#       cluster_columns = FALSE,
#       row_title_gp = gpar(fontsize = 10),
#       row_title_rot = 0,
#       column_title_gp = gpar(fontsize = 12),
#       column_names_gp = gpar(fontsize = 10),
#       rect_gp = gpar(col = "white", lwd = 0.5)
#     )
#     
#     return(ht)
#   }


plot_blood_inventory_heatmap <-
  function(inventory_mat,
           bottom_annotation,
           left_annotation,
           column_order,
           row_split = "technique",
           column_split = "consensus_signature",
           column_split_shared = NULL) {
    patient_id_col = clrs$patient_id
    timepoint_col = list("Timepoint" = clrs$timepoint)
    signature_col = list("Mutational\nsignature" = clrs$consensus_signature[names(clrs$consensus_signature)!="NA"])
    chemo_intent_col = list("Setting" = clrs$chemo_intent)
    pathology_stage_col = list("Stage" = clrs$pathology_stage)
    age_col_fun = circlize::colorRamp2(seq(30, 90, 20), brewer.pal(4, "YlOrRd"))
    
    row_split = left_annotation[[row_split]]
    column_split = bottom_annotation[[column_split]]
    
    column_gap = unit(4, "mm")
    
    signature_colors = plyr::mapvalues(
      bottom_annotation$consensus_signature,
      from = names(signature_col$`Mutational\nsignature`),
      to = signature_col$`Mutational\nsignature`
    ) %>% as.character
    timepoint_colors = plyr::mapvalues(
      left_annotation$specimen_lab_medicine_drawn_date_timepoint,
      from = names(timepoint_col$`Timepoint`),
      to = timepoint_col$`Timepoint`
    ) %>% as.character
    
    column_title = "MSK SPECTRUM"
    
    top_df = as_tibble(colSums(inventory_mat), rownames = "patient_id") %>%
      column_to_rownames(var = "patient_id") %>%
      as.matrix
    
    top_annotation <- ComplexHeatmap::columnAnnotation(
      "# samples\nper patient" = anno_barplot(
        top_df,
        gp = gpar(col = signature_colors, fill = signature_colors),
        # axis_param = list(side = "right"),
        border = FALSE
      ),
      annotation_name_gp = list(fontsize = 10),
      annotation_width = unit(c(1, 4), "cm")
    )
    
    bottom_df = bottom_annotation %>%
      dplyr::select(
        patient_age,
        gyn_diagnosis_chemo_intent_description,
        gyn_diagnosis_figo_stage,
        consensus_signature
      ) %>%
      as.data.frame
    colnames(bottom_df) <- c("Age", "Setting", "Stage", "Signature")
    
    bottom_annotation <- ComplexHeatmap::columnAnnotation(
      df = bottom_df,
      col = list(
        "Age" = age_col_fun,
        "Setting" = clrs$chemo_intent[names(clrs$chemo_intent) != "Salvage"],
        "Stage" = clrs$pathology_stage,
        "Signature" = clrs$consensus_signature[names(clrs$consensus_signature) != "NA"]
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "left",
      annotation_name_gp = list(fontsize = 10),
      # annotation_legend_param = list(
      #     signature = list(title = "Signature")),
      simple_anno_size = unit(0.35, "cm")
      # na_col = "white", gp = gpar(col = "white")
    )
    
    left_df = left_annotation %>%
      dplyr::select(specimen_lab_medicine_drawn_date_timepoint) %>%
      as.data.frame
    colnames(left_df) <- c("Timepoint")
    
    left_annotation <- ComplexHeatmap::rowAnnotation(
      df = left_df,
      col = timepoint_col,
      # col = list("Site" = clrs$tumor_supersite[!names(clrs$tumor_supersite) %in% c("UQ")]),
      show_annotation_name = FALSE,
      annotation_name_side = "bottom",
      annotation_name_gp = list(fontsize = 10),
      simple_anno_size = unit(0.35, "cm")
    )
    
    site_count = rowSums(inventory_mat)
    
    right_annotation <- ComplexHeatmap::rowAnnotation(
      "site" = anno_text(
        site_count, 
        gp = gpar(fontsize = 8)),
      "# samples\nper site" = anno_barplot(
        site_count,
        gp = gpar(col = timepoint_colors, fill = timepoint_colors),
        border = FALSE),
      annotation_name_gp = list(fontsize = 10)
    )
    
    at = seq(0, 2, by = 1)
    col_fun = circlize::colorRamp2(c(0, 1, 2),
                                   c("grey90", "grey70", "grey50"))
    
    heatmap_legend_param = list(
      at = at,
      color_bar = "discrete",
      legend_gp = gpar(fill = col_fun)
    )
    
    ht <- ComplexHeatmap::Heatmap(
      inventory_mat,
      name = "# samples",
      heatmap_legend_param = heatmap_legend_param,
      col = col_fun,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      show_row_names = FALSE,
      show_column_names = TRUE,
      border = FALSE,
      # top_annotation = top_annotation,
      # bottom_annotation = bottom_annotation,
      left_annotation = left_annotation,
      right_annotation = right_annotation,
      row_order = rownames(inventory_mat),
      column_title = "Patient",
      column_title_side = "bottom",
      row_names_side = "left",
      column_names_side = "bottom",
      row_split = row_split,
      column_order = column_order,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 10),
      row_title_rot = 0,
      column_title_gp = gpar(fontsize = 12),
      column_names_gp = gpar(fontsize = 10),
      rect_gp = gpar(col = "white", lwd = 0.5)
      # layer_fun = function(j, i, x, y, w, h, fill) {
      #   # restore_matrix() is explained after this chunk of code
      #   ind_mat = restore_matrix(j, i, x, y)
      #   for(ir in seq_len(nrow(ind_mat))) {
      #     # start from the second column
      #     for(ic in seq_len(ncol(ind_mat))[-1]) {
      #       ind1 = ind_mat[ir, ic-1] # previous column
      #       ind2 = ind_mat[ir, ic]   # current column
      #       v1 = inventory_mat[i[ind1], j[ind1]]
      #       v2 = inventory_mat[i[ind2], j[ind2]]
      #       if(v1 * v2 > 0) { # if they have the same sign
      #         grid.segments(x[ind1], y[ind1], x[ind2], y[ind2],
      #                       gp = gpar(lwd = 2))
      #         grid.points(x[c(ind1, ind2)], y[c(ind1, ind2)], 
      #                     pch = 16, size = unit(4, "mm"))
      #       }
      #     }
      #   }
      # }
      # layer_fun = function(j, i, x, y, w, h, fill) {
      #   # restore_matrix() is explained after this chunk of code
      #   ind_mat = restore_matrix(j, i, x, y)
      #   for(ir in seq_len(ncol(ind_mat))) {
      #     # start from the second row
      #     for(ic in seq_len(nrow(ind_mat))[-1]) {
      #       ind1 = ind_mat[ir-1, ic] # previous row
      #       ind2 = ind_mat[ir, ic]   # current row
      #       v1 = inventory_mat[i[ind1], j[ind1]]
      #       v2 = inventory_mat[i[ind2], j[ind2]]
      #       # print(v1)
      #       # print(v2)
      #       if(v1 == v2) { # if they have the same sign
      #         grid.segments(x[ind1], y[ind1], x[ind2], y[ind2], gp = gpar(lwd = 2))
      #         grid.points(x[c(ind1, ind2)], y[c(ind1, ind2)], pch = 16, size = unit(4, "mm"))
      #       }
      #     }
      #   }
      # }
    )
    
    return(ht)
  }

plot_tissue_inventory_heatmap <-
  function(inventory_mat,
           bottom_annotation,
           left_annotation,
           column_order,
           row_split = "technique",
           column_split = "consensus_signature",
           column_split_shared = NULL) {
    patient_id_col = clrs$patient_id
    tumor_supersite_col = list("Site" = clrs$tumor_supersite[!names(clrs$tumor_supersite)%in%c("UQ","Unknown")])
    signature_col = list("Mutational\nsignature" = clrs$consensus_signature[names(clrs$consensus_signature)!="NA"])
    chemo_intent_col = list("Setting" = clrs$chemo_intent)
    pathology_stage_col = list("Stage" = clrs$pathology_stage)
    age_col_fun = circlize::colorRamp2(seq(30, 90, 20), brewer.pal(4, "YlOrRd"))
    
    row_split = left_annotation[[row_split]]
    column_split = bottom_annotation[[column_split]]
    
    column_gap = unit(4, "mm")
    
    signature_colors = plyr::mapvalues(
      bottom_annotation$consensus_signature,
      from = names(signature_col$`Mutational\nsignature`),
      to = signature_col$`Mutational\nsignature`
    ) %>% as.character
    site_colors = plyr::mapvalues(
      left_annotation$tumor_supersite,
      from = names(tumor_supersite_col$`Site`),
      to = tumor_supersite_col$`Site`
    ) %>% as.character
    
    column_title = "MSK SPECTRUM"
    
    top_df <- as_tibble(colSums(inventory_mat), rownames = "patient_id") %>%
      column_to_rownames(var = "patient_id") %>%
      as.matrix
    
    top_annotation <- ComplexHeatmap::columnAnnotation(
      "# samples\nper patient" = anno_barplot(
        top_df,
        gp = gpar(col = signature_colors, fill = signature_colors),
        # axis_param = list(side = "right"),
        border = FALSE
      ),
      annotation_name_gp = list(fontsize = 10),
      annotation_width = unit(c(1, 4), "cm")
    )
    
    bottom_df <- bottom_annotation %>%
      dplyr::select(
        patient_age,
        gyn_diagnosis_chemo_intent_description,
        gyn_diagnosis_figo_stage,
        consensus_signature
      ) %>%
      as.data.frame
    colnames(bottom_df) <- c("Age", "Setting", "Stage", "Signature")
    
    bottom_annotation <- ComplexHeatmap::columnAnnotation(
      df = bottom_df,
      col = list(
        "Age" = age_col_fun,
        "Setting" = clrs$chemo_intent[names(clrs$chemo_intent) != "Salvage"],
        "Stage" = clrs$pathology_stage,
        "Signature" = clrs$consensus_signature[names(clrs$consensus_signature) != "NA"]
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "left",
      annotation_name_gp = list(fontsize = 10),
      # annotation_legend_param = list(
      #     signature = list(title = "Signature")),
      simple_anno_size = unit(0.35, "cm")
      # na_col = "white", gp = gpar(col = "white")
    )
    
    left_df <- left_annotation %>%
      dplyr::select(tumor_supersite, sample_type) %>%
      as.data.frame
    colnames(left_df) <- c("Site", "Sample type")
    
    left_annotation <- ComplexHeatmap::rowAnnotation(
      # df = left_df,
      # col = list("Site" = clrs$tumor_supersite[!names(clrs$tumor_supersite) %in% c("UQ","Unknown")]),
      site_sample_type = anno_simple(
        left_df$`Site`,
        col = clrs$tumor_supersite[!names(clrs$tumor_supersite) %in% c("UQ","Unknown")],
        # gp = gpar(col = site_colors, fill = site_colors),
        pch = ifelse(left_df$`Sample type` == "Tumor", 1, 2),
        # width = unit(0.35, "cm"),
        pt_size = unit(0.25, "cm")
      ),
      # 
      show_annotation_name = FALSE,
      annotation_name_side = "bottom",
      annotation_name_gp = list(fontsize = 10),
      simple_anno_size = unit(0.35, "cm")
    )
    
    site_count = rowSums(inventory_mat)
    
    right_annotation <- ComplexHeatmap::rowAnnotation(
      "site_count" = anno_text(
        site_count, 
        gp = gpar(fontsize = 8)),
      "# samples\nper site" = anno_barplot(
        site_count,
        gp = gpar(col = site_colors, fill = site_colors),
        border = FALSE),
      show_annotation_name = FALSE,
      annotation_name_side = "bottom",
      annotation_name_gp = list(fontsize = 10)
    )
    
    at = seq(0, 2, by = 1)
    col_fun = circlize::colorRamp2(c(0, 1, 2),
                                   c("grey90", "grey70", "grey50"))
    
    heatmap_legend_param = list(
      at = at,
      color_bar = "discrete",
      legend_gp = gpar(fill = col_fun)
    )
    
    ht <- ComplexHeatmap::Heatmap(
      inventory_mat,
      name = "# samples",
      heatmap_legend_param = heatmap_legend_param,
      col = col_fun,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      show_row_names = FALSE,
      show_column_names = TRUE,
      border = FALSE,
      # top_annotation = top_annotation,
      # bottom_annotation = bottom_annotation,
      left_annotation = left_annotation,
      right_annotation = right_annotation,
      row_order = rownames(inventory_mat),
      column_title = "Patient",
      column_title_side = "bottom",
      row_names_side = "left",
      column_names_side = "bottom",
      row_split = row_split,
      column_order = column_order,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 10),
      row_title_rot = 0,
      column_title_gp = gpar(fontsize = 12),
      column_names_gp = gpar(fontsize = 10),
      rect_gp = gpar(col = "white", lwd = 0.5)
    )
    
    return(ht)
  }