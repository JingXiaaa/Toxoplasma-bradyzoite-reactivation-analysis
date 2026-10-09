library(Seurat)
library(scCustomize)
library(readxl)
library(pheatmap)

sample_dirs <- c(
  EW_0h = "./Cellranger/count/EW_0h/outs/filtered_feature_bc_matrix/",
  EW_4h = "./Cellranger/count/EW_4h/outs/filtered_feature_bc_matrix/",
  EW_non_invasive = "./Cellranger/count/EW_non-invasive/outs/filtered_feature_bc_matrix/",
  EW_8h = "./Cellranger/count/EW_8h/outs/filtered_feature_bc_matrix/",
  EW_24h = "./Cellranger/count/EW_24h/outs/filtered_feature_bc_matrix/"
)

seurat_list <- lapply(names(sample_dirs), function(sample_name) {
  counts <- Read10X(data.dir = sample_dirs[[sample_name]])
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj <- subset(
    obj,
    subset = nFeature_RNA > 150 & nCount_RNA > 250 & nCount_RNA < 5000
  )
  obj$orig.sample <- sample_name
  obj
})

names(seurat_list) <- names(sample_dirs)

merge_all <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(sample_dirs),
  project = "merged"
)

merge_all <- NormalizeData(merge_all)
merge_all <- FindVariableFeatures(merge_all, selection.method = "vst", nfeatures = 2000)
merge_all <- ScaleData(merge_all, features = rownames(merge_all))
merge_all <- RunPCA(merge_all)
merge_all <- FindNeighbors(merge_all, reduction = "pca", dims = 1:20)
merge_all <- FindClusters(merge_all, resolution = 1.0)
merge_all <- RunUMAP(
  merge_all,
  reduction = "pca",
  dims = 1:20,
  min.dist = 0.3,
  n.neighbors = 30
)

# Figure 3B
DimPlot(merge_all, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Figure 3C
DimPlot(merge_all, reduction = "umap", group.by = "orig.sample")

# Figure 3D
merge_all.markers <- FindAllMarkers(
  merge_all,
  group.by = "seurat_clusters",
  logfc.threshold = 1,
  only.pos = FALSE
)

brady_genes_to_plot <- unique(read_excel("./markers_clusters.xlsx", sheet = "B clusters")$gene)
tachy_genes_to_plot <- unique(read_excel("./markers_clusters.xlsx", sheet = "T clusters")$gene)
all_genes_to_plot <- c(brady_genes_to_plot, tachy_genes_to_plot)

cluster_means <- AverageExpression(
  merge_all,
  features = all_genes_to_plot,
  group.by = "seurat_clusters",
  layer = "scale.data"
)$RNA

final_cluster_order <- c(
  "g16", "g12", "g8", "g14", "g6", "g10", "g11", "g24",
  "g23", "g1", "g13", "g17", "g3", "g2", "g4", "g7",
  "g5", "g9", "g20", "g15", "g19", "g18", "g21", "g22"
)

mat_subset <- cluster_means[, final_cluster_order]

brady_peak <- apply(mat_subset[brady_genes_to_plot, , drop = FALSE], 1, which.max)
tachy_peak <- apply(mat_subset[tachy_genes_to_plot, , drop = FALSE], 1, which.max)

ordered_brady_genes <- brady_genes_to_plot[order(brady_peak)]
ordered_tachy_genes <- tachy_genes_to_plot[order(tachy_peak)]
final_ordered_genes <- c(ordered_brady_genes, ordered_tachy_genes)
mat_final <- mat_subset[final_ordered_genes, , drop = FALSE]

row_annot <- data.frame(
  Type = c(
    rep("Bradyzoite", length(brady_genes_to_plot)),
    rep("Tachyzoite", length(tachy_genes_to_plot))
  )
)
rownames(row_annot) <- c(brady_genes_to_plot, tachy_genes_to_plot)
row_annot <- row_annot[final_ordered_genes, , drop = FALSE]

pheatmap(
  mat_final,
  color = colorRampPalette(c("#003566", "white", "#b71c2c"))(99),
  breaks = seq(-3, 3, length.out = 100),
  annotation_row = row_annot,
  annotation_colors = list(Type = c(Bradyzoite = "#1f77b4", Tachyzoite = "#ff7f0e")),
  gaps_row = length(ordered_brady_genes),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  border_color = NA
)

# Figure 3E
FeaturePlot_scCustom(
  seurat_object = merge_all,
  features = "TGME49-259020",
  colors_use = viridis_inferno_dark_high,
  pt.size = 0.1,
  label = TRUE,
  alpha_exp = 1
) + Blank_Theme()