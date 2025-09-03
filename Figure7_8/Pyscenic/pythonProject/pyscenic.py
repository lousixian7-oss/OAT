import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# ----------------------------
# Scanpy global settings
# ----------------------------
sc.settings.verbosity = 3  # verbosity levels: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()  # print versions of key dependencies
sc.set_figure_params(dpi=300, figsize=(5,5), dpi_save=600)  # high-res figures

# ----------------------------
# Function: Convert Seurat-exported files into an AnnData object
# ----------------------------
def seurat_to_adata(counts,      # path to counts.mtx (exported from Seurat)
                    meta,        # path to metadata.csv (cell metadata from Seurat)
                    gene_name,   # path to gene_names.csv
                    pca,         # path to pca.csv
                    reduction1,  # column name for UMAP/TSNE dimension 1 (e.g., "UMAP_1")
                    reduction2): # column name for UMAP/TSNE dimension 2 (e.g., "UMAP_2")

    # Load expression matrix (genes x cells)
    X = io.mmread(counts)

    # Create AnnData object (cells x genes, so transpose needed)
    adata = anndata.AnnData(X=X.transpose().tocsr())

    # Load cell metadata
    cell_meta = pd.read_csv(meta)

    # Load gene names
    with open(gene_name, 'r') as f:
        gene_names = f.read().splitlines()

    # Assign cell metadata (obs) and gene metadata (var)
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']  # set cell barcodes as index
    adata.var.index = gene_names            # set gene names as index

    # Load PCA embeddings
    pca = pd.read_csv(pca)
    pca.index = adata.obs.index
    adata.obsm['X_pca'] = pca.to_numpy()

    # Load UMAP/TSNE embeddings
    adata.obsm['X_umap'] = np.vstack(
        (adata.obs[reduction1].to_numpy(),
         adata.obs[reduction2].to_numpy())
    ).T

    return adata

# ----------------------------
# Apply the function to load Seurat-exported files
# ----------------------------
sce_test = seurat_to_adata(
    counts='D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.2.SCENIC_all/seurat_counts.mtx',
    meta='D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.2.SCENIC_all/seurat_metadata.csv',
    gene_name='D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.2.SCENIC_all/seurat_gene_names.csv',
    pca='D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.2.SCENIC_all/seurat_pca.csv',
    reduction1='UMAP_1',
    reduction2='UMAP_2'
)

# Print summary of the AnnData object to confirm successful conversion
print(sce_test)

# ----------------------------
# Visualization: UMAP plots
# ----------------------------
# Example usage:
# sc.pl.umap(sce_test, color=['annotation_myo'])
# sc.pl.umap(sce_test, color='VEGFA')

sc.pl.umap(
    sce_test, color='anntation',
    add_outline=True,
    legend_loc='on data',
    legend_fontsize=4,
    legend_fontoutline=1,
    frameon=False,
    title='clustering of cells',
    palette='Set1',
    size=8
)

# ----------------------------
# Save AnnData object for downstream Scanpy analysis
# ----------------------------
sce_test.write('sce_test.h5ad')
