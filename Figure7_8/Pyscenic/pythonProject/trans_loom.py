import os
import loompy as lp
import numpy as np
import scanpy as sc

# ----------------------------
# Step 1. Load expression matrix
# ----------------------------
# Read expression matrix exported from R/Seurat (genes x cells, CSV format)
# The matrix should have genes as rows and cells as columns
x = sc.read_csv("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.2.SCENIC_all/counts.csv")

# ----------------------------
# Step 2. Define row (gene) and column (cell) attributes
# ----------------------------
# Row attributes: assign gene names
row_attrs = {"Gene": np.array(x.var_names)}

# Column attributes: assign cell IDs
col_attrs = {"CellID": np.array(x.obs_names)}

# ----------------------------
# Step 3. Create a loom file
# ----------------------------
# Transpose the expression matrix (cells x genes → required by loom)
# Save as "sce.loom" for downstream SCENIC analysis
lp.create("sce.loom", x.X.transpose(), row_attrs, col_attrs)
