library(Seurat)
library(readxl)
library(ggplot2)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)

# Figure 4A
ref <- UpdateSeuratObject(readRDS("./S.O.tg_RH_boothroyd.rds"))
orthologs <- read_excel("./TGGT1_ME49_orthologs.xlsx")
gene_map <- setNames(orthologs$TGME49, orthologs$TGGT1)

counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
norm <- GetAssayData(ref, assay = "RNA", layer = "data")

rename_genes <- function(x, map) {
  new_names <- map[rownames(x)]
  new_names[is.na(new_names)] <- rownames(x)[is.na(new_names)]
  rownames(x) <- new_names
  x
}

counts <- rename_genes(counts, gene_map)
norm <- rename_genes(norm, gene_map)

ref[["RNA_renamed"]] <- CreateAssay5Object(counts = counts, data = norm)
DefaultAssay(ref) <- "RNA_renamed"

common_genes <- intersect(rownames(ref), rownames(merge_all))

anchors <- FindTransferAnchors(
  reference = ref,
  query = merge_all,
  features = common_genes,
  dims = 1:2
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = ref$phase,
  dims = 1:2
)

merge_all <- AddMetaData(merge_all, predictions)

DimPlot(
  merge_all,
  reduction = "UMAP_Loupe",
  group.by = "predicted.id",
  pt.size = 0.01
) +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("G1.a", "G1.b", "S", "M", "C")
  ) +
  theme_classic()


# Figure 4G
sce <- as.SingleCellExperiment(merge_all)

sce <- slingshot(
  sce,
  clusterLabels = "RNA_snn_res.1_updated",
  reducedDim = "UMAP",
  start.clus = "6"
)

pt <- slingPseudotime(sce)[, 1]
umap <- reducedDim(sce, "UMAP")

plot_df <- data.frame(
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  pseudotime = pt
)

curves <- slingCurves(sce)
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  x <- as.data.frame(curves[[i]]$s)
  colnames(x) <- c("UMAP_1", "UMAP_2")
  x$lineage <- i
  x
}))

pt_colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = pseudotime), size = 0.8) +
  scale_color_gradientn(colors = pt_colors, na.value = "grey80") +
  geom_path(
    data = curve_df,
    aes(UMAP_1, UMAP_2, group = lineage),
    inherit.aes = FALSE,
    linewidth = 0.75
  ) +
  coord_fixed() +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Pseudotime")


# Figures 4H and 4I
genes_to_plot <- c(
  "TGME49-259020", "TGME49-291040", "TGME49-268860",
  "TGME49-233460", "TGME49-233480", "TGME49-268850"
)

samples_to_plot <- c("EW_0h", "EW_4h", "EW_8h", "EW_24h")

gene_colors <- c(
  "TGME49-233480" = "#E7533C",
  "TGME49-268850" = "#F29B7F",
  "TGME49-233460" = "#FFB41D",
  "TGME49-259020" = "#4DBBD5",
  "TGME49-291040" = "#3C5488",
  "TGME49-268860" = "#8391B4"
)

plot_pseudotime_genes <- function(clusters) {
  expr <- GetAssayData(merge_all, assay = "RNA", layer = "data")

  spline_data <- lapply(genes_to_plot, function(gene) {
    df <- data.frame(
      pseudotime = merge_all$pseudotime_scaled,
      expression = expr[gene, ],
      sample = merge_all$orig.sample,
      cluster = merge_all$subcluster
    )

    df <- df[
      complete.cases(df) &
        df$cluster %in% clusters &
        df$sample %in% samples_to_plot,
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
      gene = gene
    )
  })

  spline_data <- do.call(rbind, spline_data)

  ggplot(
    spline_data,
    aes(pseudotime, expression, color = gene, group = gene)
  ) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = gene_colors) +
    theme_classic() +
    labs(
      x = "Pseudotime",
      y = "Normalized expression",
      color = "Gene ID"
    )
}

# Figure 4H
plot_pseudotime_genes(c("6", "13", "7", "5", "9"))

# Figure 4I
plot_pseudotime_genes(c("6", "13", "7", "19", "18", "21b"))
