Cluster_Highlight_Image_DimPlot <- function(
    seurat_object,
    fov,
    split_by,
    highlight_col = "grey",
    background_col = "navy",
    cluster_level,
    sample = "all", 
    cr = T,
    size = 0.5,
    isPdf = F,
    w = 8,
    h = 11,
    col = 4
){
  
  # Initialize an empty list to store the generated color vectors
  color_vectors <- list()
  
  clusters <- unlist(unique(seurat_object[[cluster_level]]))
  
  # Initialize an empty list to store the generated color vectors
  color_vectors <- list()
  
  # Use a for loop to generate the desired number of color vectors
  for (i in 1:length(clusters)) {
    # Create a vector filled with blue
    current_vector <- rep(background_col, length = length(clusters))
    
    # Replace one blue with grey based on the loop iteration
    current_vector[i] <- highlight_col
    
    # Add the current vector to the list of color vectors
    color_vectors[[i]] <- current_vector
  }  
  
  
  images <- list()
  
  if(isPdf){
    fn <- paste0(getwd(), "/05_plots/",Sys.Date(), "_", sample, "_", cluster_level, "_clusterImageHighlight.pdf")
    # grDevices::cairo_pdf(file=fn, width = w, height = h, onefile=T)
    CairoPDF(file=fn, width = w, height = h, onefile=T)
  }
  # img <- list()
  for(i in 1:length(clusters)){
    print(ImageDimPlot(seurat_object,
                       fov=fov,
                       group.by = cluster_level,
                       cols = color_vectors[[i]],
                       size = size,
                       crop=cr,
                       combine=T) & theme(legend.position = 'none') & ggtitle(paste0("Highlighted Cluster: ",levels(clusters)[i])))
    # wrap_plots(img, ncol = col) + plot_layout(guides = "collect")
  }
  
  if(isPdf){
    dev.off()
  }
}


#' @importFrom magrittr %>% %<>%
NULL

# TODO:
# - edit args

#' Intermediate solution to \code{subset()}:
#' subset FOVs/centroids if selected cells are NOT found in each FOV
#' NOTE: some code parts and args are taken from SeuratObject

#' Function params/args:
#' @param object An S4 object or A \code{FOV} object
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param Update.slots If to update slots of an object
#' @param Update.object If to update final object, default to TRUE.
#' @param ... Arguments passed to \code{subset()} and other methods


subset_opt <- function(
    object = NULL, 
    subset = NULL, 
    cells = NULL, 
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    ...)
{
  
  if (Update.slots) { 
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }
  
  message("Cloing object..")
  obj_subset <- object
  
  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) { 
    cells <- Cells(obj_subset)[cells]
  }
  
  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }
  
  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset, 
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }
  
  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else { 
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq, 
             function(i) { 
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
             }) %>% unlist
  }
  
  if (all(cells_check)) { 
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells, 
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                     cells = cells, 
                     idents = idents, 
                     features = features, 
                     ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }
    
  } else { 
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <- 
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells, 
                       idents = idents, 
                       features = features, 
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n", 
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] } 
    
    # subset final object
    message("..subset final object")
    obj_subset %<>% 
      base::subset(cells = cells,
                   idents = idents,
                   features = features, 
                   ...)
  }
  
  if (Update.object && !class(obj_subset) == "FOV") { 
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }
  
  message("Object is ready!")
  return(obj_subset)
  
}

plot_stat_meta <- function(dataset, 
                      plot_type, 
                      clusterLabel = "seurat_clusters",
                      group_by = "sample",
                      pal_setup = 'Set2',
                      plot_ratio = 1,
                      text_size = 10,
                      tilt_text = FALSE) {
  
  if (is.data.frame(pal_setup)) {
    pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
  } else {
    pal <- pal_setup
  }
  
  stat <- tibble::tibble(group = dataset[[group_by]][[1]], 
                         cluster = dataset[[clusterLabel]][[1]])
  stat %<>%
    group_by(.data$group, 
             .data$cluster) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  ncolors <- if (plot_type == 'prop_fill') {
    length(unique(dataset[[clusterLabel]][[1]]))
  } else {
    length(unique(dataset[[group_by]][[1]]))
  }
  colors <- set_colors(pal, ncolors)
  
  thm <- theme(aspect.ratio = plot_ratio,
               legend.title = element_text(size = text_size),
               legend.text = element_text(size = text_size),
               axis.title = element_text(size = text_size),
               axis.text = element_text(size = text_size),
               axis.title.x = element_blank()
  ) + theme_bw()
  
  thm2 <- theme(legend.position = "none")
  thm3 <- theme(axis.text.x = element_text(angle = 45, 
                                           vjust = 0.5))
  
  switch(plot_type,
         group_count = stat %>%
           group_by(.data$group) %>%
           summarise(sum(n)) %>%
           ggplot(aes(x = .data$group, 
                      y = .data$`sum(n)`)) +
           geom_col(aes(fill = .data$group)) +
           geom_text(aes(label = .data$`sum(n)`), 
                     vjust = -0.5, 
                     size = text_size * 0.35) +
           scale_fill_manual(values = colors) + 
           labs(x = group_by, y = "Number of Cells") + 
           thm + thm2 + if (tilt_text) {thm3},
         
         prop_fill = ggplot(stat) + 
           geom_bar(aes(x = .data$group, 
                        y = .data$freq, 
                        fill = .data$cluster), 
                    position = "fill", 
                    stat = "identity") +
           scale_y_continuous(labels = scales::percent) +
           scale_fill_manual(values = colors, 
                             name = "Cluster") +
           labs(x = group_by, y = "Proportion") + 
           thm + if (tilt_text) {thm3},
         
         prop_multi = stat %>%
           mutate(freq = round(.data$freq, 3)) %>%
           ggplot() + 
           geom_bar(aes(x = .data$group,
                        y = .data$freq, 
                        fill = .data$group), 
                    stat = "identity") +
           geom_text(aes(x = .data$group, 
                         y = .data$freq, 
                         label = scales::percent(.data$freq)), 
                     vjust = -0.5, 
                     size = text_size * 0.35) +
           scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                              labels = scales::percent_format()) +
           facet_wrap(~ cluster, 
                      ncol = 4, 
                      scales = "free") +
           scale_fill_manual(values = colors, 
                             name = "Group") +
           labs(x = NULL, 
                y = "Proportion") + 
           theme(strip.text.x = element_text(size = text_size)) + 
           thm + thm2 + if (tilt_text) {thm3},
         
         stop("Unknown plot type")
  )
}

set_colors <- function(pal, n) {
  
  if (all(pal %in% rownames(brewer.pal.info))) {
    num <- c()
    for (i in seq(length(pal))) {
      num[i] <- brewer.pal.info[pal[i],][[1]]
    }
    full_pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
  } else if (all(are_colors(pal))) {
    full_pal <- pal
  } else {
    stop('Incorrect palette setup. Please input valid RColorBrewer palette names or color names.')
  }
  
  if (n <= length(full_pal)) {
    return(full_pal[1:n])
  } else {
    warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.", 
            immediate. = TRUE)
    return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
  }
}

are_colors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

getBulkExpression <- function(seu,
                              cellType,
                              donor,
                              assays=NULL,
                              disease,
                              diseaseLevels){
  bulk <- AggregateExpression(seu, 
                              return.seurat = T, 
                              slot = "counts", 
                              assays = assays, 
                              group.by = c(cellType,
                                           donor,
                                           disease))
  
  bulk$celltype <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
  bulk$donor <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
  bulk$disease <- sapply(strsplit(Cells(bulk), split = "_"), "[", 3)
  bulk$disease <- factor(x = bulk$disease, levels = diseaseLevels)
  return(bulk)
}


plot_stat <- function(dataset, 
                      plot_type, 
                      group_by = "sample",
                      pal_setup = 'Set2',
                      clusterLabel,
                      plot_ratio = 1,
                      text_size = 10,
                      tilt_text = FALSE) {
  
  # if (is.data.frame(pal_setup)) {
  #   pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
  # } else {
  #   pal <- pal_setup
  # }
  
  stat <- tibble::tibble(group = dataset[[group_by]][[1]], 
                         cluster = dataset[[clusterLabel]][[1]])
  stat %<>%
    group_by(.data$group, 
             .data$cluster) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  ncolors <- if (plot_type == 'prop_fill') {
    length(unique(dataset[[clusterLabel]][[1]]))
  } else {
    length(unique(dataset[[group_by]][[1]]))
  }
  
  colors <- set_colors(pal, ncolors)
  
  thm <- theme(aspect.ratio = plot_ratio,
               legend.title = element_text(size = text_size),
               legend.text = element_text(size = text_size),
               axis.title = element_text(size = text_size),
               axis.text = element_text(size = text_size),
               axis.title.x = element_blank()
  ) + theme_bw()
  
  thm2 <- theme(legend.position = "none")
  thm3 <- theme(axis.text.x = element_text(angle = 45, 
                                           vjust = 0.5))
  
  switch(plot_type,
         group_count = stat %>%
           group_by(.data$group) %>%
           summarise(sum(n)) %>%
           ggplot(aes(x = .data$group, 
                      y = .data$`sum(n)`)) +
           geom_col(aes(fill = .data$group)) +
           geom_text(aes(label = .data$`sum(n)`), 
                     vjust = -0.5, 
                     size = text_size * 0.35) +
           scale_fill_manual(values = colors) + 
           labs(x = group_by, y = "Number of Cells") + 
           thm + thm2 + if (tilt_text) {thm3},
         
         prop_fill = ggplot(stat) + 
           geom_bar(aes(x = .data$group, 
                        y = .data$freq, 
                        fill = .data$cluster), 
                    position = "fill", 
                    stat = "identity") +
           scale_y_continuous(labels = scales::percent) +
           scale_fill_manual(values = colors, 
                             name = "Cluster") +
           labs(x = group_by, y = "Proportion") + 
           thm + if (tilt_text) {thm3},
         
         prop_multi = stat %>%
           mutate(freq = round(.data$freq, 3)) %>%
           ggplot() + 
           geom_bar(aes(x = .data$group,
                        y = .data$freq, 
                        fill = .data$group), 
                    stat = "identity") +
           geom_text(aes(x = .data$group, 
                         y = .data$freq, 
                         label = scales::percent(.data$freq)), 
                     vjust = -0.5, 
                     size = text_size * 0.35) +
           scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                              labels = scales::percent_format()) +
           facet_wrap(~ cluster, 
                      ncol = 4, 
                      scales = "free") +
           scale_fill_manual(values = colors, 
                             name = "Group") +
           labs(x = NULL, 
                y = "Proportion") + 
           theme(strip.text.x = element_text(size = text_size)) + 
           thm + thm2 + if (tilt_text) {thm3},
         
         stop("Unknown plot type")
  )
}

getProportions <- function(seu_obj,
                           sample_id,
                           annotation,
                           loc,
                           status){
  
  # get a list of all cells and their corresponding identifiers (usually a patient ID)
  sample <- seu_obj[[sample_id]]
  # generate a factor of sample, then add barcodes
  s <- sample[,]
  names(s) <- rownames(sample)
  # get a list of all cells and their corresponding annotation, convert to factor and add barcodes
  annot <- seu_obj[[annotation]]
  a <- annot[,]
  names(a) <- rownames(annot)
  # get a list of all cells and their corresponding location, convert to factor and add barcodes
  loca <- seu_obj[[loc]]
  l <- loca[,]
  names(l) <- rownames(loca)
  
  

  cell_counts <- table(s, a)
  cell_counts.loc <- table(s, l)
  
  # Calculate total cells per original ID
  total_cells_per_orig_id <- rowSums(cell_counts)
  
  
  # Get unique combinations of sample_id and annotation
  unique_combinations <- unique(cbind(seu_obj[[sample_id]], seu_obj[[annotation]], seu_obj[[status]], seu_obj[[loc]]))
  unique_combinations[] <- lapply(unique_combinations, as.character)
  
  # Initialize an empty data frame
  result_table <- data.frame(
    sample = character(),
    status = character(),
    loc = character(),
    annotation = character(),
    cells = numeric(),
    total.cells = numeric(),
    proportion = numeric(),
    total.loc = numeric(),
    proportion.loc = numeric(),
    # area.loc = numeric(),
    # area.prop = numeric(),
    stringsAsFactors = TRUE
  )
  
  # Iterate through each unique combination and populate the data frame
  for (i in 1:nrow(unique_combinations)) {
    current_orig_id <- unique_combinations[i, 1]
    current_annotation <- unique_combinations[i, 2]
    current_status <- unique_combinations[i, 3]
    current_loc <- unique_combinations[i, 4]
    
    # Get cell count and calculate proportion
    current_count <- cell_counts[current_orig_id, current_annotation]
    current_proportion <- current_count / total_cells_per_orig_id[current_orig_id]
    current_loc_total <- cell_counts.loc[current_orig_id,current_loc]
    
    current_proportion.loc <- current_count / current_loc_total
    
    # Add new row to the data frame with total cells
    result_table <- rbind.data.frame(result_table, c(current_orig_id,
                                                     current_annotation,
                                                     current_status,
                                                     current_loc,
                                                     current_count,
                                                     total_cells_per_orig_id[current_orig_id],
                                                     current_proportion,
                                                     current_loc_total,
                                                     current_proportion.loc))
    colnames <- c(sample_id,
                  annotation,
                  status,
                  loc,
                  "count",
                  "total.all",
                  "prop.all",
                  "total.loc",
                  "prop.loc")
    setnames(result_table, old = names(result_table), new = colnames)
  }
  
  result_table <- result_table[order(result_table[[annotation]], result_table[[status]]), ]
  rownames(result_table) <- NULL
  
  # Select the last 5 columns and convert to numeric
  result_table[, -c(1:4)] <- lapply(result_table[, -c(1:4)], as.numeric)
  
  return(result_table)
  
}

getCellsPerArea <- function(seu_obj,
                            sample_id = "pt.id",
                            section_id = "orig.ident",
                            annotation,
                            loc,
                            status,
                            area){
  
  
  # Get unique combinations of sample_id and annotation
  unique_combinations <- unique(cbind(seu_obj[[sample_id]], seu_obj[[annotation]], seu_obj[[status]], seu_obj[[loc]], seu_obj[[section_id]]))
  unique_combinations[] <- lapply(unique_combinations, as.character)
  
  # Initialize an empty data frame
  result_table <- data.frame(
    sample = character(),
    status = character(),
    loc = character(),
    annotation = character(),
    cells = numeric(),
    total.cells = numeric(),
    area.loc = numeric(),
    area.prop = numeric(),
    stringsAsFactors = TRUE
  )
  
  meta <- seu_obj@meta.data
  
  # Iterate through each unique combination and populate the data frame
  for (i in 1:nrow(unique_combinations)) {
    current_sample_id <- unique_combinations[i, 1]
    current_section_id <- unique_combinations[i, 5]
    current_annotation <- unique_combinations[i, 2]
    current_status <- unique_combinations[i, 3]
    current_loc <- unique_combinations[i, 4]
    current_meta <- meta[meta[[annotation]] == current_annotation & meta[[section_id]] == current_section_id, ]
    loc_area <- unique(current_meta$total_area_mm2)/1000 # divide by 1k because I calculated total area wrong - fix once hpc is back up
    
    total_nuc_area <- sum(current_meta$nucleus_area)/1000000 # convert from um2 to mm2
    total_cells <- length(current_meta$orig.ident)
    
    
    
    # Add new row to the data frame with total cells
    result_table <- rbind.data.frame(result_table, c(current_section_id,
                                                     current_sample_id,
                                                     current_annotation,
                                                     current_status,
                                                     current_loc,
                                                     loc_area,
                                                     total_nuc_area,
                                                     total_cells))
    colnames <- c(section_id,
                  sample_id,
                  annotation,
                  status,
                  loc,
                  "region_area",
                  "total_nuc_area",
                  "total_cells")
    setnames(result_table, old = names(result_table), new = colnames)
  }
  
  result_table <- result_table[order(result_table[[annotation]], result_table[[status]]), ]
  rownames(result_table) <- NULL
  
  # Select the last 3 columns
  last_cols <- names(result_table)[(ncol(result_table) - 2):ncol(result_table)]  
  
  # Convert the last 3 columns to numeric (coerce errors to NA)
  result_table[, last_cols] <- lapply(result_table[, last_cols], function(x) as.numeric(x, na.rm = TRUE))
  
  result_table <- result_table %>%
    group_by_at(c(sample_id,annotation)) %>%
    mutate(total_nuc_area_pt = sum(total_nuc_area),
           region_area_pt = sum(region_area),
           total_cells_pt = sum(total_cells))
  
  result_table <- result_table[!duplicated(result_table[c(2,3)]),]
  
  # Get cell count and calculate proportion
  result_table$cells_per_region_mm2 <- result_table$total_cells_pt/result_table$region_area_pt
  result_table$nuc_area_per_region_mm2 <- result_table$total_nuc_area_pt/result_table$region_area_pt
  
  
  return(result_table)
  
}

propFacetPlot <- function(prop_data,
                          x_var,
                          y_var,
                          fill_var,
                          add_var,
                          facet_var,
                          stat_test){
  p_load(ggpubr, rstatix)
  stat.test <- result_table %>%
    group_by(names(.[2])) %>%
    t_test(names(.[9]) ~ names(.[3])) %>%
    adjust_pvalue(method = stat_test) %>%
    add_significance()
  stat.test 
  
  bxp <- ggboxplot(
    result_table, x = x_var, y = y_var, fill = fill_var, 
    facet.by = facet_var, scales = "free", add = add_var
  )
  # Make facet and add p-values
  stat.test <- stat.test %>% add_xy_position(x = x_var)
  bxp +  
    stat_pvalue_manual(
      stat.test, hide.ns = TRUE,
      label = "{p.adj}{p.adj.signif}"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    theme_ggprism_mod()
}

cropRegion <- function(
    orig.dir,
    new.dir,
    cropFile
){
  R.utils::copyDirectory(orig.dir, new.dir)
  # R.utils::copyDirectory("/data/williamsdrw/20240405__165613__Human_gingiva2_VT/HV137A/", "/data/williamsdrw/20240405__165613__Human_gingiva2_VT/HV137A_2/")
  # 
  # cropFile <- fread("/data/williamsdrw/20240405__165613__Human_gingiva2_VT/Selection_1_coordinates.csv")
  
  crop <- fread(cropFile)  
  xmax <- max(crop$X)
  xmin <- min(crop$X)
  ymax <- max(crop$Y)
  ymin <- min(crop$Y)
  cb <- fread(paste0(new.dir, "cell_boundaries.csv.gz"))
  
  cb_new <-  cb[ cb$vertex_x > xmin &
                   cb$vertex_x < xmax &
                   cb$vertex_y > ymin &
                   cb$vertex_y < ymax, ]
  
  unique_cells <- length(unique(cb_new$cell_id))
  
  c <- fread(paste0(new.dir, "cells.csv.gz"))
  
  c_new <- c[c$x_centroid > xmin &
               c$x_centroid < xmax &
               c$y_centroid > ymin &
               c$y_centroid < ymax, ]
  
  t <- fread(paste0(new.dir, "transcripts.csv.gz"))
  
  t_new <- t[t$x_location > xmin &
               t$x_location < xmax &
               t$y_location > ymin &
               t$y_location < ymax, ]
  
  barcode <- fread(paste0(new.dir, "cell_feature_matrix/barcodes.tsv.gz"), header = F)
  rem <- barcode$V1 %in% c_new$cell_id
  bar <- barcode[barcode$V1 %in% c_new$cell_id,]
  
  m <- Matrix::readMM(paste0(new.dir, "cell_feature_matrix/matrix.mtx.gz"))
  m <- m[,rem]
  
  fwrite(c_new, file = paste0(new.dir, "cells.csv.gz"))
  fwrite(t_new, file = paste0(new.dir, "transcripts.csv.gz"))
  fwrite(cb_new, file = paste0(new.dir, "cell_boundaries.csv.gz"))
  fwrite(bar, file = paste0(new.dir, "cell_feature_matrix/barcodes.tsv.gz"), col.names = F)
  Matrix::writeMM(m, file = paste0(new.dir, "cell_feature_matrix/matrix.mtx.gz"))
  
}


cropRegion_complex <- function(
    orig.dir,
    new.dir,
    cropFile
){
  R.utils::copyDirectory(orig.dir, new.dir)
  
  files <- list.files(new.dir, full.names = TRUE)
  
  # Filter files not ending with .gz using regular expressions
  files_to_delete <- files[!grepl(".gz$", files)]
  
  # Check if there are any files to delete
  if (length(files_to_delete) > 0) {
    # Shell command to delete (use with caution!)
    shell_cmd <- paste0("rm -rf ", files_to_delete)
    shell_cmd <- shell_cmd[c(1:5,7:18)]
    for(i in seq_along(shell_cmd)){
      system(shell_cmd[i], intern = TRUE)
    }
    message(paste0("Deleted ", length(shell_cmd), " files- and directories."))
  } else {
    message("No files found for deletion (all files end with .gz).")
  }
  
  
  crop <- fread(cropFile)  
  cb <- fread(paste0(new.dir, "cell_boundaries.csv.gz"))
  sf_df <- st_as_sf(cb, coords = c("vertex_x", "vertex_y"))
  
  c_hull <- 
    st_as_sf(crop, coords = c("X", "Y")) %>% 
    st_combine() %>% st_convex_hull()
  
  crop_df <-
    st_intersection(sf_df, c_hull)
  
  cb <- fread(paste0(new.dir, "cell_boundaries.csv.gz"))
  cb_new <- cb[cb$cell_id %in% unique(crop_df$cell_id),]
  c <- fread(paste0(new.dir, "cells.csv.gz"))
  c_new <- c[c$cell_id %in% unique(crop_df$cell_id),]
  t <- fread(paste0(new.dir, "transcripts.csv.gz"))
  t_new <- t[t$cell_id %in% unique(crop_df$cell_id),]
  
  
  barcode <- fread(paste0(new.dir, "cell_feature_matrix/barcodes.tsv.gz"), header = F)
  rem <- barcode$V1 %in% c_new$cell_id
  bar <- barcode[barcode$V1 %in% c_new$cell_id,]
  m <- Matrix::readMM(paste0(new.dir, "cell_feature_matrix/matrix.mtx.gz"))
  m <- m[,rem]
  
  fwrite(c_new, file = paste0(new.dir, "cells.csv.gz"))
  fwrite(t_new, file = paste0(new.dir, "transcripts.csv.gz"))
  fwrite(cb_new, file = paste0(new.dir, "cell_boundaries.csv.gz"))
  fwrite(bar, file = paste0(new.dir, "cell_feature_matrix/barcodes.tsv.gz"), col.names = F)
  Matrix::writeMM(m, file = paste0(new.dir, "cell_feature_matrix/matrix.mtx.gz"))
}

# Function to read and combine dataframes
obtain_10x_areas <- function(files) {
  # Initialize an empty dataframe
  combined_df <- data.frame()
  
  # Loop through each file
  for (file in files) {
    # Read the data using read.csv with gz argument
    df <- read.csv(gzFile(file), stringsAsFactors = FALSE)
    # Append the data to the combined dataframe
    combined_df <- rbind(combined_df, df)
  }
  
  # Return the combined dataframe
  return(combined_df)
}