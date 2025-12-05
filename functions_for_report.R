# functions_for_report.R

# Setup some directories in this scope so that things are saved properly

# .report_cfg <- new.env(parent = emptyenv())
# 
# # called once from the Rmd to set paths
# set_report_paths <- function(base.path = NULL,
#                              ggplot.dir   = NULL,
#                              dataset.path  = NULL,
#                              out.dir   = NULL) {
#   if (!is.null(base.path)) .report_cfg$base.path <- base.path
#   if (!is.null(ggplot.dir))   .report_cfg$ggplot.dir   <- ggplot.dir
#   if (!is.null(dataset.path))  .report_cfg$dataset.path  <- dataset.path
#   if (!is.null(out.dir))   .report_cfg$out.dir   <- out.dir
# }
# 
# # helper if you want it
# get_report_cfg <- function() .report_cfg


# Helper function that formats filenames. Allows for adding
# suffix to config for alternate versions of analysis (e.g. 
# different clustering resolutions, filters, etc)
suffix_path <- function(base_dir, filename, suffix = "") { 
	# cfg <- get_report_cfg() 
	# base_dir <- cfg$base.dir

	root <- tools::file_path_sans_ext(filename) 
	ext  <- tools::file_ext(filename) 
	fn   <- if (nzchar(suffix)) paste0(root, suffix, if (nzchar(ext)) paste0(".", ext) else "") else filename 
	file.path(base_dir, fn)
}

# Helper function to make the sample names pretty and workable
sanitize_sample_names <- function(formatted.samples){

  ids <- formatted.samples
  
  # Use sub instead of gsub in the first call. gsub replaces all, while sub 
  # only replaces the first instance. The first dash is just there as formatting
  # in the config file. This allows for names with dashes in them.
  ids <- sub("^-", "", ids)
  ids <- gsub("\ -", "\ ", ids)
  ids <- gsub("\"", "", ids)
  ids <- strsplit(ids, "[[:space:]]")
  ids <- unlist(ids)
  
  
  ids

}


# This got complicated as I was learning about how ggsave saves images (i.e. 
# be default it saves huge images) and how the pdf rendering presents those
# images (i.e. it makes them huge if they're save huge)
save_png_plot <- function(p, filename, ggplot.dir,
                          width = 7, height = 5, dpi = 150) {

	# cfg <- get_report_cfg()
	# dir <- cfg$ggplot.dir 
	
	if (!dir.exists(ggplot.dir)) dir.create(ggplot.dir, recursive = TRUE, showWarnings = FALSE) 
	out <- file.path(ggplot.dir,  filename) 
	print(out)
	
	# works for ggplot or any grid grob via cowplot::ggdraw 
	if (inherits(p, "ggplot")) { 
		ggplot2::ggsave(out, plot = p, width = width, height = height, 
				units = "in", dpi = dpi, limitsize = TRUE) 
	} else if (!is.null(p$gtable)) { # pheatmap object 
		png(out, width = width, height = height, units = "in", res = dpi) 
		grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off() 
	} else if (inherits(p, "grob") || inherits(p, "gTree")) { 
		png(out, width = width, height = height, units = "in", res = dpi) 
		grid::grid.newpage(); grid::grid.draw(p); dev.off() 
	} else { 
		# last resort: try plotting it 
		png(out, width = width, height = height, units = "in", res = dpi) 
		print(p); dev.off() 
	} 
	invisible(out) 
}


# Filter out the metrics we don't want, and format everything nicely
make_qc_table <- function(ids, dataset.path){
  
	d10x.metrics <- lapply(ids, function(sample.name){ 
				       metrics.path.cleaned <- file.path(dataset.path, sample.name, "metrics_summary.csv") 
				       metrics.path.raw     <- file.path(dataset.path, sample.name, "outs/metrics_summary.csv") 
				       tryCatch({ 
					       if (file.exists(metrics.path.cleaned)) { 
						       read.csv(metrics.path.cleaned, colClasses = "character") 
					       } else if (file.exists(metrics.path.raw)) { 
						       read.csv(metrics.path.raw, colClasses = "character") 
					       } else { 
						       data.frame() 
					       } 
				       }, error = function(cond){ 
					       message(paste0("Error loading the sample ", sample.name, ": ", conditionMessage(cond))) 
					       data.frame() 
				       }) 
				})

	# --- Make them a data.frame --- 
	experiment.metrics <- do.call("rbind", d10x.metrics) 
	rownames(experiment.metrics) <- ids 
	df <- if (inherits(experiment.metrics, "data.frame")) { 
		experiment.metrics 
	} else { 
		as.data.frame(do.call(rbind, experiment.metrics), stringsAsFactors = FALSE) 
	} 

	# --- Keep only desired metrics --- 
	keep <- c("Estimated.Number.of.Cells", 
		  "Mean.Reads.per.Cell", 
		  "Median.Genes.per.Cell", 
		  "Valid.Barcodes") 
	df_sub <- df[, keep, drop = FALSE] 
	
	# --- Prettify metric names --- 
	pretty_names <- function(x) { 
		s <- gsub("\\.", " ", x) 
		s <- tolower(s) 
		substr(s, 1, 1) <- toupper(substr(s, 1, 1)) 
		s 
	} 
	colnames(df_sub) <- pretty_names(colnames(df_sub))
  
 
	# --- Transpose; put 'metric' first --- 
	qc.table <- t(df_sub) 
	qc.table <- as.data.frame(qc.table, stringsAsFactors = FALSE) 
	qc.table$metric <- rownames(qc.table) 
	row.names(qc.table) <- NULL 
	qc.table <- qc.table[, c(ncol(qc.table), 1:(ncol(qc.table)-1))]
  
  	# --- Order metrics (optional) --- 
	metric_order <- c("Estimated number of cells", 
			  "Mean reads per cell", 
			  "Median genes per cell", 
			  "Valid barcodes") 
	qc.table <- qc.table %>% 
		mutate(metric = factor(metric, levels = metric_order)) %>% 
		arrange(metric) %>% 
		mutate(metric = as.character(metric))
	# return 
	qc.table
}


# This is mostly used for putting the sample names in the right order so things are consistant in the plots
make_metric_table_from_list <- function(seurat.list, 
					metric, 
					sample_levels = NULL, 
					sort_alpha = TRUE) { 
	df <- do.call(rbind, lapply(names(seurat.list), function(s) { 
					    so <- seurat.list[[s]] 
					    stopifnot(metric %in% colnames(so[[]])) 
					    data.frame(Sample = s, value = so[[metric, drop = TRUE]], 
						       row.names = colnames(so), check.names = FALSE) 
				       })) 
	names(df)[2] <- metric 

	if (!is.null(sample_levels)) { 
		df$Sample <- factor(df$Sample, levels = sample_levels) 
	} else if (sort_alpha) { 
		df$Sample <- factor(df$Sample, levels = sort(unique(df$Sample))) 
	} else { 
		df$Sample <- factor(df$Sample, levels = unique(df$Sample))
  }
  
  df
}


make_metric_table <- function(obj, sample_col, metric, sample_levels = NULL, sort_alpha = TRUE) {
  stopifnot(metric %in% colnames(obj[[]]))
  out <- data.frame(
    Sample = obj[[sample_col, drop = TRUE]],
    value  = obj[[metric, drop = TRUE]],
    row.names = colnames(obj),
    check.names = FALSE
  )
  names(out)[2] <- metric
  
  if (!is.null(sample_levels)) {
    out$Sample <- factor(as.character(out$Sample), levels = sample_levels)
  } else if (sort_alpha) {
    out$Sample <- factor(as.character(out$Sample), levels = sort(unique(as.character(out$Sample))))
  } else {
    out$Sample <- factor(as.character(out$Sample), levels = unique(as.character(out$Sample)))
  }
  out
}


################################################################################
# Takes a DF and prints out a nice formatted table to the rendered PDF
################################################################################

render_formatted_table <- function(input.df) {
  # ---- format numbers (no manual LaTeX) ----
  num_only <- function(x) as.numeric(gsub("[^0-9.]", "", x))
  
  fmt_one_cell <- function(metric, val_chr){
    if (is.na(val_chr) || val_chr == "") return(NA_character_)
    if (identical(metric, "Valid barcodes")){
      v <- num_only(val_chr) / 100
      if (is.na(v)) return(NA_character_)
      percent(v, accuracy = 0.1)                # e.g. "93.8%"
    } else {
      v <- num_only(val_chr)
      if (is.na(v)) return(NA_character_)
      comma(v)                                   # e.g. "24,047"
    }
  }
  
  for (j in 2:ncol(input.df)) {
    input.df[[j]] <- mapply(fmt_one_cell, input.df$metric, input.df[[j]])
  }
  input.df[is.na(input.df)] <- ""
  input.df[] <- lapply(input.df, as.character)
  
  # ---- render (escape TRUE) ----
  input.df %>%
    rename(Metric = metric) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          align    = c("l", rep("r", ncol(input.df)-1)),
          caption  = "Sequencing / QC metrics by sample",
          # escape = TRUE by default; keep it that way!
          longtable = FALSE) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      font_size = 12                                # body font
    ) %>%
    column_spec(1, bold = TRUE) %>%                 # bold Metric column
    add_header_above(c(" " = 1, "Samples" = ncol(input.df)-1)) %>%
    row_spec(0, angle = 45, font_size = 8)          # rotate/shrink header row (sample names)
}


# helper to locate 10x paths
read_10x_any <- function(sample.name, dataset.path) {
	# cfg <- get_report_cfg()
	# dataset.path <- cfg$dataset.path

	if (dir.exists(file.path(dataset.path, sample.name, "outs"))){ 
		# tenx.h5.path     <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix.h5") 
		tenx.matrix.path <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix") 
	} else { 
		# tenx.h5.path     <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix.h5") 
		tenx.matrix.path <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix") 
	} 

	# if (isTRUE(import.h5) && file.exists(tenx.h5.path)) { 
	# 	Read10X_h5(tenx.h5.path) 
	# } else { 
		Read10X(tenx.matrix.path)
#   }
}

# Remove Y chromosome genes first
read_10x_no_Y_chr <- function(sample.name,
                              y_genes_file = "y_genes_symbols.txt",
                              verbose = TRUE) {
  # 1) Resolve paths exactly as in read_10x_any
  if (dir.exists(file.path(dataset.path, sample.name, "outs"))) {
    tenx.h5.path     <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path, sample.name, "outs/raw_matrix")
  } else {
    tenx.h5.path     <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path, sample.name, "raw_matrix")
  }
  
  # 2) Load Y-chromosome gene names (symbols)
  if (!file.exists(y_genes_file)) {
    stop("Y gene list file not found: ", y_genes_file)
  }
  y_genes <- readLines(y_genes_file)
  y_genes <- unique(y_genes[y_genes != ""])
  
  if (verbose) {
    message("Loaded ", length(y_genes), " Y-chromosome gene symbols from: ", y_genes_file)
  }
  
  # 3) Read 10x data (h5 if requested and present)
  if (isTRUE(import.h5) && file.exists(tenx.h5.path)) {
    mat <- Read10X_h5(tenx.h5.path)
  } else {
    mat <- Read10X(tenx.matrix.path)
  }
  
  # 4) Remove Y genes from the matrix (or list of matrices)
  drop_y <- function(m) {
    if (is.null(rownames(m))) {
      warning("Matrix has no rownames; cannot drop Y genes reliably.")
      return(m)
    }
    is_y <- rownames(m) %in% y_genes
    n_y  <- sum(is_y)
    if (verbose) {
      message("  Found ", n_y, " Y-chromosome genes in this feature set; removing them.")
    }
    m[!is_y, , drop = FALSE]
  }
  
  if (is.list(mat)) {
    # multi-modal object: apply to each
    mat <- lapply(mat, drop_y)
  } else {
    mat <- drop_y(mat)
  }
  
  return(mat)
}
read_pipseeker_any <- function(sample.name) {
  if (dir.exists(file.path(dataset.path, sample.name, "outs"))){
    tenx.h5.path     <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path, sample.name, "outs/raw_matrix")
  } else {
    tenx.h5.path     <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path,sample.name, "raw_matrix") # for pipseeker
    
    # tenx.matrix.path <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix")
  }
  if (isTRUE(import.h5) && file.exists(tenx.h5.path)) {
    Read10X_h5(tenx.h5.path)
  } else {
    Read10X(tenx.matrix.path)
  }
}


vln_boxplot <- function(data, x, y, title) {
  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)

  p
}

shaded_vln_boxplot <- function(data, 
                               x, 
                               y, 
                               title, 
                               show.outliers = TRUE,
                               this_floor = -Inf, 
                               this_ceiling = Inf){
  ggplot(data, aes(x = {{x}}, y = {{y}})) +
    # shaded bands (put first so they sit behind the geoms)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = this_ceiling, ymax = Inf,
             fill = "red", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = this_floor,
             fill = "red", alpha = 0.08) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_hline(yintercept = c(this_floor, this_ceiling),
               color = "red", linetype = "dashed", linewidth = 0.6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
}


feat_plots_top_genes <- function(cluster.num,
                                 seurat.obj,
                                 marker.genes.df,
                                 filename.suffix,
                                 ggplot.dir,
                                 assay.used = "SCT",
                                 genes_per_page = 6,
                                 fig_width = 9,     # inches
                                 fig_height = 10,    # inches
                                 dpi = 150,
                                 umap_pt_size = 0.02,
                                 vln_pt_size = 0.05,
                                 base_text = 10) {
  
	# cfg <- get_report_cfg()
	# ggplot.dir <- cfg$ggplot.dir
  if (!dir.exists(ggplot.dir)) dir.create(ggplot.dir, recursive = TRUE, showWarnings = FALSE)
  
  # paths
  feat.plot.path <- file.path(ggplot.dir,
                              paste0("clustering_", 
                                     filename.suffix, 
                                     "_marker_gene_featPlot_cl_", 
                                     cluster.num, 
                                     "_", assay.used, ".png"))
  vln.plot.path  <- file.path(ggplot.dir,
                              paste0("clustering_", 
                                     filename.suffix, 
                                     "_marker_gene_vlnPlot_cl_", 
                                     cluster.num, 
                                     "_", assay.used, ".png"))
  
  # gene list
  curr.markers <- marker.genes.df[marker.genes.df$cluster == cluster.num, ]
  genes <- na.omit(curr.markers$gene)
  genes <- genes[seq_len(min(length(genes), genes_per_page))]
  
  # ---------- FEATURE PLOTS (rasterized) ----------
  fp_list <- Seurat::FeaturePlot(
    seurat.obj,
    features   = genes,
    reduction  = "umap.postint",
    cols       = c("grey85", "navy"),
    ncol       = 3,
    pt.size    = umap_pt_size,
    label.size = base_text,
    combine    = FALSE
    # raster     = TRUE,        # <- key for file size & speed
    # raster.dpi = c(dpi, dpi)
  )
  
  fp_list <- lapply(fp_list, function(p) {
    p +
      theme(
        axis.text.x = element_text(size = base_text, angle = 45, hjust = 1),
        axis.text.y = element_text(size = base_text),
        plot.title  = element_text(size = base_text + 1)
      )
  })
  
  fp_grid <- cowplot::plot_grid(plotlist = fp_list, ncol = 3)
  title_g <- cowplot::ggdraw() + cowplot::draw_label(
    paste0("Cluster ", cluster.num, " — top marker features"),
    fontface = "bold", size = base_text + 2
  )
  fp_final <- cowplot::plot_grid(title_g, fp_grid, ncol = 1, rel_heights = c(0.12, 1))
  
  ggplot2::ggsave(
    filename = feat.plot.path, plot = fp_final,
    width = fig_width, height = fig_height, units = "in",
    dpi = dpi, limitsize = TRUE
  )
  
  # ---------- VIOLIN PLOTS ----------
  vp_list <- Seurat::VlnPlot(
    object   = seurat.obj,
    features = genes,
    pt.size  = vln_pt_size,
    combine  = FALSE
  )
  
  vp_list <- lapply(vp_list, function(p) {
    p +
      theme(
        axis.text.x = element_text(size = base_text, angle = 45, hjust = 1),
        axis.text.y = element_text(size = base_text),
        plot.title  = element_text(size = base_text + 1)
      )
  })
  
  vp_grid  <- cowplot::plot_grid(plotlist = vp_list, ncol = 3)
  title_v  <- cowplot::ggdraw() + cowplot::draw_label(
    paste0("Cluster ", cluster.num, " — top marker violins"),
    fontface = "bold", size = base_text + 2
  )
  vp_final <- cowplot::plot_grid(title_v, vp_grid, ncol = 1, rel_heights = c(0.12, 1))
  
  ggplot2::ggsave(
    filename = vln.plot.path, plot = vp_final,
    width = fig_width, height = fig_height, units = "in",
    dpi = dpi, limitsize = TRUE
  )
  
  invisible(list(feature_plot = feat.plot.path, violin_plot = vln.plot.path))
}
