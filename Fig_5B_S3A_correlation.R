library(Seurat)
library(Hmisc)

genes_A <- read.csv("./TFs.csv", stringsAsFactors = FALSE)[[1]]
genes_B <- read.csv("./tachy-markers.csv", stringsAsFactors = FALSE)[[1]]

clusters_to_use <- c("6", "13", "7", "19", "18", "21b")
samples_to_use <- c("EW_0h", "EW_4h", "EW_8h", "EW_24h")

obj <- subset(
  merge_all,
  subset = subcluster %in% clusters_to_use &
           orig.sample %in% samples_to_use
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

genes_A <- intersect(genes_A, rownames(expr))
genes_B <- intersect(genes_B, rownames(expr))

combined_expr <- rbind(
  expr[genes_A, , drop = FALSE],
  expr[genes_B, , drop = FALSE]
)

cor_result <- rcorr(
  t(as.matrix(combined_expr)),
  type = "spearman"
)

# Correlation coefficients and P values between the two gene sets
final_r_values <- cor_result$r[genes_A, genes_B, drop = FALSE]
final_p_values <- cor_result$P[genes_A, genes_B, drop = FALSE]
