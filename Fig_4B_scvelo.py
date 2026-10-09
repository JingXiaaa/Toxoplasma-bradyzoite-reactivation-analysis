from pathlib import Path
import pandas as pd
import scanpy as sc
import scvelo as scv

samples = ["EW_0h", "EW_4h", "EW_non_invasive", "EW_8h", "EW_24h"]

cellranger_dir = Path("./Cellranger")
input_dir = Path("./velocity_input")
output_dir = Path("./velocity_output")
output_dir.mkdir(exist_ok=True)

cluster_files = {
    "EW_0h": "clusters_0h.csv",
    "EW_4h": "clusters_4h.csv",
    "EW_non_invasive": "clusters_non_invasive.csv",
    "EW_8h": "clusters_8h.csv",
    "EW_24h": "clusters_24h.csv",
}

umap_files = {
    "EW_0h": "umap_0h.csv",
    "EW_4h": "umap_4h.csv",
    "EW_non_invasive": "umap_non_invasive.csv",
    "EW_8h": "umap_8h.csv",
    "EW_24h": "umap_24h.csv",
}


def load_sample(sample):
    clusters = pd.read_csv(
        input_dir / cluster_files[sample],
        index_col=0
    )

    umap = pd.read_csv(
        input_dir / umap_files[sample],
        index_col=0
    )

    count_dir = cellranger_dir / f"count_{sample}" / "outs" / "filtered_feature_bc_matrix"
    loom_file = cellranger_dir / f"count_{sample}" / "velocyto" / f"count_{sample}.loom"

    adata = sc.read_10x_mtx(
        count_dir,
        var_names="gene_symbols"
    )

    # Retain cells used in the final analysis.
    keep = adata.obs_names.intersection(clusters.index)
    adata = adata[keep].copy()

    cluster_values = clusters.reindex(adata.obs_names)
    adata.obs["Clusters"] = cluster_values.iloc[:, 0].values
    adata.obsm["X_umap"] = umap.reindex(adata.obs_names).to_numpy()

    velocity_data = scv.read(loom_file)
    adata = scv.utils.merge(adata, velocity_data)

    return adata


for sample in samples:
    adata = load_sample(sample)

    scv.pp.normalize_per_cell(adata, enforce=True)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=20, n_neighbors=30)

    scv.tl.recover_dynamics(adata, n_jobs=4)
    scv.tl.velocity(adata, mode="dynamical", n_jobs=4)
    scv.tl.velocity_graph(adata, n_jobs=4)

    adata.write(
        output_dir / f"{sample}_velocity_graph.h5ad",
        compression="gzip"
    )
