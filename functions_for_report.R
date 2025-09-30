suffix_path <- function(base_dir, filename, suffix = "") {
  root <- tools::file_path_sans_ext(filename)
  ext  <- tools::file_ext(filename)
  fn   <- if (nzchar(suffix)) paste0(root, suffix, if (nzchar(ext)) paste0(".", ext) else "") else filename
  file.path(base_dir, fn)
}

save_png_plot <- function(plot, filename, this.width=10, this.height=6){
  save.path <- file.path(ggplot.directory, paste0("part1", part1.suffix, "_", filename))
  ggsave(save.path, plot, width=this.width, height=this.height)
}


# helper to locate 10x paths
read_10x_any <- function(sample.name) {
  if (dir.exists(file.path(dataset.path, sample.name, "outs"))){
    tenx.h5.path     <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path, sample.name, "outs/raw_feature_bc_matrix")
  } else {
    tenx.h5.path     <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix.h5")
    tenx.matrix.path <- file.path(dataset.path, sample.name, "raw_feature_bc_matrix")
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

shaded_vln_boxplot <- function(data, x, y, title, floor = -Inf, ceiling = Inf){
  ggplot(data, aes(x = {{x}}, y = {{y}})) +
    # shaded bands (put first so they sit behind the geoms)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ceiling, ymax = Inf,
             fill = "red", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = floor,
             fill = "red", alpha = 0.08) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_hline(yintercept = c(floor, ceiling),
               color = "red", linetype = "dashed", linewidth = 0.6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
}


feat_plots_top_genes <- function(cluster.num, 
                                 seurat.obj, 
                                 marker.genes.df, 
                                 ggplot.dir, 
                                 filename.suffix){  
  
  saved.fig.width = 25
  saved.fig.height = 20
  
  # print(paste0("Cluster number ", cluster.num, ", width: ", saved.fig.width, 
  #              ", height: ", saved.fig.height))
  curr.markers <- marker.genes.df[marker.genes.df$cluster == cluster.num,] 
  num.markers.to.plot <- 6
  # print(head(curr.markers))

  feat.plot.path <- file.path(ggplot.dir, 
                              paste0("clustering_", filename.suffix, 
                                     "_marker_gene_featPlot_cl_", cluster.num, ".png"))
  
  vln.plot.path <- file.path(ggplot.dir, 
                             paste0("clustering_", filename.suffix, 
                                    "_marker_gene_vlnPlot_cl_", cluster.num, ".png"))

  feat.plot.list <- FeaturePlot(
    seurat.obj,
    curr.markers$gene[1:num.markers.to.plot],
    reduction = "umap.postint",
    cols = c("lightgrey", "blue"),
    ncol = 3,
    pt.size = 2,
    label.size = 30)

  # Fix the axis labels and main title sizes for the feature plots
  for (i in 1:num.markers.to.plot){
    feat.plot.list[[i]] <- feat.plot.list[[i]] + theme(axis.text.x = element_text(angle=45, hjust=1, size=30),
                                                       axis.text.y = element_text(size = 30),
                                                       plot.title = element_text(size=30))
  }

  plot.title <- ggdraw() + draw_label(paste0("Cluster ", cluster.num), fontface = 'bold', size = 30)
  print(cowplot::plot_grid(plot.title, feat.plot.list, ncol = 1, rel_heights = c(0.1, 1)))


  ggsave(feat.plot.path, width=saved.fig.width, height=saved.fig.height)


  vln.plot.list <- VlnPlot(object = seurat.obj,
                      features = curr.markers$gene[1:num.markers.to.plot],
                      pt.size = 0.05)

  # Fix the axis labels and main title sizes for the violin plots
  for (i in 1:num.markers.to.plot){
    vln.plot.list[[i]] <- vln.plot.list[[i]] + theme(axis.text.x = element_text(angle=45, hjust=1, size=30),
                                                       axis.text.y = element_text(size = 30),
                                                       plot.title = element_text(size=30))
  }
  print(cowplot::plot_grid(plot.title, vln.plot.list, ncol = 1, rel_heights = c(0.1, 1)))


  ggsave(vln.plot.path, width=saved.fig.width, height=saved.fig.height)


}
