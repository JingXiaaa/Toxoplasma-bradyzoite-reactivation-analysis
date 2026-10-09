library(Seurat)
library(ggplot2)

highlight_genes <- c(
  "TGME49-200385",
  "TGME49-311100",
  "TGME49-238440"
)

genes_to_plot <- c(
  "TGME49-259020",
  "TGME49-291040",
  "TGME49-268860",
  "TGME49-233460",
  "TGME49-233480",
  "TGME49-268850",
  highlight_genes
)

samples_to_use <- c("EW_0h", "EW_4h", "EW_8h", "EW_24h")

gene_colors <- c(
  "TGME49-233480" = "#E7533C",
  "TGME49-268850" = "#F29B7F",
  "TGME49-233460" = "#FFB41D",
  "TGME49-259020" = "#4DBBD5",
  "TGME49-291040" = "#3C5488",
  "TGME49-268860" = "#8391B4",
  "TGME49-200385" = "#4C752A",
  "TGME49-311100" = "#A3B18A",
  "TGME49-238440" = "black"
)

plot_gene_splines <- function(clusters_to_use) {
  expr <- GetAssayData(
    merge_all,
    assay = "RNA",
    layer = "data"
  )

  spline_data <- lapply(genes_to_plot, function(gene) {
    df <- data.frame(
      pseudotime = merge_all$pseudotime_scaled,
      expression = expr[gene, ],
      sample = merge_all$orig.sample,
      cluster = merge_all$subcluster
    )

    df <- df[
      complete.cases(df) &
        df$cluster %in% clusters_to_use &
        df$sample %in% samples_to_use,
    ]

    if (nrow(df) < 10) return(NULL)

    fit <- smooth.spline(
      x = df$pseudotime,
      y = df$expression,
      spar = 1
    )

    data.frame(
      pseudotime = fit$x,
      expression = fit$y,
      gene = gene,
      highlighted = gene %in% highlight_genes
    )
  })

  spline_data <- do.call(rbind, spline_data)

  ggplot(
    spline_data,
    aes(
      x = pseudotime,
      y = expression,
      color = gene,
      group = gene
    )
  ) +
    geom_line(
      data = subset(spline_data, !highlighted),
      linetype = "dashed",
      alpha = 0.6,
      linewidth = 2
    ) +
    geom_line(
      data = subset(spline_data, highlighted),
      linewidth = 2
    ) +
    scale_color_manual(values = gene_colors) +
    theme_classic() +
    labs(
      x = "Pseudotime",
      y = "Normalized expression",
      color = "Gene ID"
    )
}

# Figure S2C
plot_gene_splines(c("6", "13", "7", "5", "9"))

# Figure S2D
plot_gene_splines(c("6", "13", "7", "19", "18", "21b"))
