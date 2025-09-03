README:

This repository contains the complete workflow for Mendelian Randomization (MR), Colocalization, and single-cell RNA-seq integration in druggable genes in periodontitis. It covers data sources (eQTLGen, FINNGEN, DGIdb, Finan et al.), MR and coloc analyses (Fig2A–B), PheWAS validation (Fig2D), functional enrichment (Fig2H), and downstream single-cell analyses (clustering, GSEA/UCell, pseudotime, metabolism, CellChat, SCENIC). Scripts, input data, and results are organized in a modular structure for reproducibility and direct mapping to manuscript figures.
————————————————————————
### Mendelian Randomization and Colocalization
#### H1.0 Data sources

**Exposure**

Druggable genes: DGIdb + Finan et al. → 5,883 genes

Cis-eQTLs: eQTLGen (31,684 samples, 16,987 genes)

After intersection: 3,275 exposure genes

**Outcome**

FINNGEN R11 (2024-06-24, gingivitis/periodontitis)

N=406,876 (118,404 cases; 288,472 controls)

#### H2.0 Software and dependencies

- R ≥ 4.3.0

- Key packages: TwoSampleMR, MRPRESSO, coloc, clusterProfiler, org.Hs.eg.db, locuscomparer

- Optional: MetaForest (forest plot), plink (LD clumping)

#### H3.0 Directory layout
- data/ (eQTLGen, FINNGEN, OAT.txt, DPP8.txt)
- scripts/ (Fig2A.MR_function.R, Fig2B.coloc.R, plot_PheWAS.R, enrichment_GO_KEGG.R)
- results/ (mendelian test/, coloc/, phewas/, enrichment/)

#### H4.0 Analysis workflow

**Step 1. Mendelian Randomization (Fig2A)**

Input: Fig2A.eqtlgen-lite.txt

LD clumping → Fig2A.MR.clump

Run Fig2A.MR_function.R (get_F, MR-PRESSO, Steiger, power)

Run Fig2A.MR.R → results in mendelian test/

Visualization: Fig2A.MetaForest.R (forest plot)

**Step 2. Colocalization (Fig2B)**

Input: cis-eQTLs (e.g., OAT.txt, DPP8.txt) + FINNGEN GWAS

Run Fig2B.coloc.R with default priors (p1=1e−4, p2=1e−4, p12=1e−5)

Filter PPH4 > 0.9

Visualization: Fig2B.coloc_print.R / locuscomparer → OAT.pdf, DPP8.pdf

**Step 3. PheWAS (Fig2D)**

Data: AstraZeneca PheWAS portal (UK Biobank, 450k samples, >17k traits)

Threshold: 2×10⁻⁹

Run plot_PheWAS.R → plot/OAT_binary_*_1.png, _2.png

**Step 4. Enrichment (Fig2H)**

GO: enrichGO(ont=BP/CC/MF, p<0.05, q<0.05)

KEGG: enrichKEGG("hsa")

Output: GO.txt, KEGG.txt, and dot/bar plots

Highlight key pathways (e.g., Arginine and proline metabolism)

#### H5.0 Expected Outputs 

**Fig2A.MR**

In: Fig2A.eqtlgen-lite.txt, Fig2A.MR.periodontal.out

Out (per exposure):

`dat.xlsx`, `res.xlsx`, `MR-PRESSO.xlsx`

`<exp>_scatter.pdf`, `_forest.pdf`, `_funnel.pdf`, `_leave_one_out.pdf`

Summary: mendelian testres.csv (includes FDR on IVW p-values)

**Fig2B.coloc

In: cis-eQTL region files (OAT.txt, DPP8.txt), FINNGEN GWAS

Out: coloc_data_result1.csv (PPH0–PPH4), coloc_data_result2.csv (per-SNP), OAT.pdf, DPP8.pdf

**Fig2D.plot_PheWAS

In: AZ portal CSV (e.g., OAT_binary_*.csv), traits.xlsx

Out: plot/<title>_<type>_1.png, _2.png, and <title_name>_<type>.xlsx

**Fig2H.GO_KEGG

In: core_gene_set.csv (SYMBOLs)

Out: GO.txt, KEGG.txt, and the bar/bubble PDFs above

### Single-cell RNA-seq Preprocessing and Annotation Workflow
#### H1.0 Input data（GSE171213）
Path: `D:\scRNAseq\06.Proj_processing\Sen_yazhouyan\H1. harmony\GSE171213\`
- **HC samples (Healthy controls)**
    - `HC1`: GSM5220920_HC1.cellname.list.txt.gz, GSM5220920_HC1.counts.tsv.gz
    - `HC2`: GSM5220921_HC2.cellname.list.txt.gz, GSM5220921_HC2.counts.tsv.gz
    - `HC3`: GSM5220922_HC3.cellname.list.txt.gz, GSM5220922_HC3.counts.tsv.gz
    - `HC4`: GSM5220923_HC4.cellname.list.txt.gz, GSM5220923_HC4.counts.tsv.gz
- **PD samples (Patients)**
    - `PD1`: GSM5220924_PD1.cellname.list.txt.gz, GSM5220924_PD1.counts.tsv.gz
    - `PD2`: GSM5220925_PD2.cellname.list.txt.gz, GSM5220925_PD2.counts.tsv.gz
    - `PD3`: GSM5220926_PD3.cellname.list.txt.gz, GSM5220926_PD3.counts.tsv.gz
    - `PD4`: GSM5220927_PD4.cellname.list.txt.gz, GSM5220927_PD4.counts.tsv.gz
    - `PD5`: GSM5220928_PD5.cellname.list.txt.gz, GSM5220928_PD5.counts.tsv.gz
#### H2.0 Code execution order
- **Note:** Specific parameter settings (e.g., filtering thresholds, number of PCs, clustering resolutions) should be referred to in the **Methods section of the manuscript** and the **accompanying code**.

**Step 1. Data Preprocess**
- Load raw count matrices (`counts.tsv.gz`) and cell barcode mapping (`cellname.list.txt.gz`).
- Replace cell indices (C1, C2, …) with actual barcodes.
- Construct Seurat objects (`CreateSeuratObject`) with QC filters:
    - min.cells = 3
    - min.features = 200
- Save Seurat object list → `scRNAlist.rds`.

**Step 2. QC and Merge**
- Run `scRNAdataPreProcessing()` for each object:
    - Filter cells by nFeature_RNA, nCount_RNA, %mito.
    - Add cell cycle scoring.
    - Run DoubletFinder (`DF.classify`) and DecontX for ambient RNA removal.
- Merge all objects (`merge`).
- Remove doublets and low-RNA cells.
- Save merged object → `preprocessed.rds`.

**Step 3. Harmony Integration and Clustering**
- Normalize (LogNormalize), find variable features, scale data.
- Run PCA → ElbowPlot → Harmony batch correction.
- Run UMAP, neighbors, clusters (`FindClusters`, res = 0.1–1.5).
- Assign disease status (`HC` vs `PD`).
- Visualize UMAP plots (cluster-level and sample-level).
- Save integrated object → `preprocessed.rds`.
    

**Step 4. Cell Type Annotation**
- Identify cluster markers using **COSG** (`cosg`).
- Validate expression of marker genes (`FeaturePlot`, `DotPlot`).
- Manually map clusters to cell types (T cells, B cells, Fibroblasts, etc.).
- Subset to major immune/structural populations.
- Save annotated object → `annotation_finished.rds`.

**Step 5. Visualization (SCP & scplotter)**
- **UMAP** plots with `CellDimPlot`.
- **Cell proportion plots** with `CellStatPlot` (stack, ring, pie).
- **Feature expression plots** (`FeatureDimPlot`, `FeatureStatPlot`) for OAT/DPP8.
- **DotPlot** of canonical marker genes.

#### H3.0 Expected output
- **UMAP** 
- **Cell proportion plots** 
- **Feature expression plots** 
- **DotPlot**

### GSEA + UCell
- To identify enriched pathways at the **cluster level** using GSEA on average expression profiles, and validate key pathways at the **single-cell level** using UCell scoring.  
- This approach compensates for the limitations of GSEA in scRNA-seq by combining **bulk-like analysis (GSEA)** and **single-cell robustness (UCell)**.
#### H1.0 Input data
- Preprocessed and annotated Seurat object

#### H2.0 Code execution order
1. **Average expression per cluster**
    - `AverageExpression()` → compute cluster-level expression matrix.
    - Select top 1000 variable genes.
    - Visualize correlation heatmap between clusters.
2. **Prepare gene sets**
    - Use `msigdbr` to download MSigDB **Reactome (C2:CP:REACTOME)** pathways.
3. **Run GSEA for each cluster**
    - Sort average expression values (descending).
    - Filter genes with expression > 0.1.
    - Run `GSEA()` with Reactome TERM2GENE.
    - Save enrichment results per cluster.
4. **Combine enrichment scores**
    - Collect `enrichmentScore` across clusters into a matrix (`es.max`).
    - Plot heatmap of enrichment scores across pathways.
5. **Identify top pathways per cluster**
    - Calculate fold-change = (cluster-specific ES − median of others).
    - Select **top 5 enriched pathways per cluster**.
    - Plot clustered heatmap of selected pathways.
6. **Validate with UCell**
    - Extract key Reactome pathways (e.g. ECM organization, ECM degradation, MET signaling).
    - Run `AddModuleScore_UCell()` at single-cell resolution.
    - Visualize pathway activity on UMAP using `FeaturePlot()`.
#### H3.0 Expected output
- Cluster-level results (GSEA)
- Pathway-level validation (UCell)


### Monocle2 + CytoTRACE Pseudotime Analysis
- o reconstruct fibroblast differentiation trajectories using **Monocle2** and refine the root state assignment with **CytoTRACE**.  Monocle2 builds a pseudotime trajectory based on ordering genes, while CytoTRACE provides an independent stemness score to help identify the most primitive (stem-like) cell states.
#### H1.0 Input data
- Seurat fibroblast subset object
- Conversion to Monocle2 CDS for trajectory inference
- Gene expression matrix (`RNA@counts`) used for CytoTRACE scoring
#### H2.0 Code execution order
1. **Convert Seurat → Monocle2**
    - `seurat_to_monocle()` with RNA counts.
    - Normalize with `estimateSizeFactors()` and `estimateDispersions()`.
    - Gene filtering: keep genes expressed in ≥10 cells.
    - Add UMI counts and remove extreme cells.
2. **Identify ordering genes**
    - Differential expression test (`differentialGeneTest`) across cell groups (`annotation3`).
    - Select significant genes (q-value < 0.01) as ordering genes.
    - Apply `setOrderingFilter()` and visualize with `plot_ordering_genes()`.
3. **Trajectory construction (Monocle2)**
    - Dimension reduction with `reduceDimension(..., DDRTree)`.
    - Order cells along pseudotime with `orderCells()`.
    - Save trajectory object (`cds_DGT`).
    - Visualization:
        - `plot_cell_trajectory()` colored by annotation, state, pseudotime, disease status.
        - Gene-specific pseudotime plots (e.g., OAT).
4. **CytoTRACE scoring**
    - Input = raw counts matrix from `obj.fib`.
    - Run `CytoTRACE()` to compute stemness scores.
    - Add scores to Seurat metadata and Monocle2 CDS.
    - Identify root state = state with **highest average CytoTRACE score**.
    - Re-order cells in trajectory with `orderCells(root_state = ...)`.
5. **Visualization (integrated)**
    - Pseudotime trajectory colored by CytoTRACE score (Spectral palette).
    - UMAP visualization of CytoTRACE scores in Seurat (`FeatureDimPlot`).
    - Export plots (PDF): pseudotime, OAT dynamics, CytoTRACE gradients.
#### H3.0 Expected output
- Monocle2 trajectory plots
- Gene expression dynamics
- CytoTRACE results


### IrGSEA + ScMetabolism
- To investigate fibroblast metabolic pathway activity by integrating **irGSEA** scoring (AUCell, UCell, singscore, ssGSEA) with curated **KEGG** and **Reactome metabolic gene sets** (from `scMetabolism`).  
- This allows identification of metabolism-related pathway differences across fibroblast subgroups (e.g., OAT⁺ vs OAT⁻ Fib), and correlation of pathway activity with specific genes such as OAT. 
#### H1.0 Input data
- Seurat fibroblast subset object
- **Metabolic gene sets**:
	- KEGG: `KEGG_metabolism_nc.gmt`
	- Reactome: `REACTOME_metabolism.gmt`  
	    (provided in `scMetabolism` package data folder)
#### H2.0 Code execution order
- **Load metabolic gene sets**
    - Import KEGG and Reactome metabolic pathways via `read.gmt()`.
    - Convert to proper list format for `irGSEA`.
- **Run irGSEA scoring**
    - Apply `irGSEA.score()` with methods: `AUCell`, `UCell`, `singscore`, `ssgsea`.
    - Run separately for KEGG and Reactome gene sets.
- **Integrate and visualize results**
    - `irGSEA.integrate()` → combine results at cluster level (`annotation2`).
    - Visualization:
        - **Heatmap** of top 50 pathways.
        - **Bubble plot** of pathway enrichment.
        - **Density scatterplot** (UMAP-based pathway activity, e.g. glycolysis).
- **Correlation analysis (OAT⁺ fibroblasts)**
    - Subset `OAT+ Fib`.
    - Extract expression of genes of interest (`OAT`, `DPP8`).
    - Extract pathway activity scores (e.g., `N-Glycan biosynthesis`, `Glutamate/glutamine metabolism`).
    - Merge gene expression with pathway scores → correlation analysis.
    - Visualization:
        - Scatterplot with regression line.
        - Pearson correlation coefficient + p-value.
        - Marginal density plots along X/Y axis.
- **Save results**
    - Save integrated irGSEA results (`.rds`).
    - Export heatmap, bubble plot, and correlation plots (PNG/PDF).
#### H3.0 Expected output
- Pathway-level enrichment
- Correlation analysis

### pseudo-bulk analysis
- This section describes the pseudo-bulk strategy to compare OAT+ vs OAT− fibroblasts by aggregating single-cell counts to the sample level, performing differential expression (DESeq2), and running pathway enrichment (fgsea / ssGSEA).

#### H1.0 Input data
- `obj.fib`: Seurat object containing fibroblasts, annotated with OAT status (`OAT+ Fib` / `OAT- Fib`).
- Metadata includes:
    - `orig.ident` → sample ID (e.g., HC1–HC4, PD1–PD5).
    - `annotation2` → OAT subgroup (OAT+ / OAT−).
- Gene set files for enrichment analysis:
    - `KEGG_metabolism_nc.gmt` (KEGG metabolic pathways).
    - `REACTOME_metabolism.gmt` (Reactome metabolic pathways).

#### H2.0 Code execution order
- **Prepare pseudo-bulk counts**
    - Concatenate `sample ID + OAT status` into `sample_group`.
    - Use `AggregateExpression()` to sum raw counts at the sample-group level.
- **Build metadata**
    - Create metadata (`meta`) with columns: `sample`, `OAT_status`.
    - Assign `OAT_pos` / `OAT_neg` labels for downstream DE.
- **Differential expression with DESeq2**
    - Construct DESeq2 object (`dds`).
    - Filter lowly expressed genes (`rowSums >= 10`).
    - Run DE analysis (`DESeq()`).
    - Extract significant DEGs and classify into **up / down / stable**.
- **Volcano plot of DEGs**
    - Plot `log2FoldChange` vs `-log10(pvalue)`.
    - Highlight up- and downregulated genes.
- **Pathway enrichment with fgsea**
    - Rank genes by log2FC.
    - Run fgsea on KEGG/Reactome pathways.
    - Extract top enriched pathways.
- **Barplot of enriched pathways**
    - Select top 15 significant pathways.
    - Plot NES values, grouped by OAT+ or OAT− enrichment.
- **ssGSEA validation**
    - Run ssGSEA (`gsva`) across pseudo-bulk samples.
    - Compare pathway activity (e.g., "Arginine and proline metabolism") between OAT+ and OAT− groups with boxplots.
- **Statistical testing & heatmap**
    - Perform Mann-Whitney tests for each pathway.
    - Select top 10 significant pathways.
    - Plot heatmap of ssGSEA scores, annotated by group (OAT+ / OAT−).
- **Styled fgsea scatterplot**
    - Plot NES for all pathways.
    - Scale point size by set size, color by `-log10(pvalue)`.
    - Highlight key metabolic pathways (e.g., arginine, glycosaminoglycan, linoleic acid metabolism).

#### H3.0 Expected output
- **Pseudo-bulk counts matrix** (genes × sample_groups)
- **DESeq2 results table** of DEGs between OAT+ and OAT− fibroblasts.
- **ssGSEA validation results**: pathway activity per sample with boxplots.
### Cell-Cell Communication Analysis (CellChat)
- This section describes the workflow to compare cell-cell communication between **healthy controls (HC)** and **patients (PD)** using the **CellChat** R package.
#### H1.0 Input data
- **`obj`**: annotated Seurat object containing all cells with `diseaseStatus` (HC/PD) and `annotation` metadata.
- **`obj.fib`**: fibroblast-only Seurat object with refined clustering (e.g., OAT+ vs OAT− fibroblasts).

#### H2.0 Code execution order
- Step 1. Split Seurat object into HC and PD groups
- Step 2. Create CellChat objects for HC and PD
- Step 3. Run CellChat preprocessing (HC and PD separately)
- Step 4. Save and reload objects
- Step 5. Merge CellChat objects for comparison
- Step 6. Visualization of global interaction networks
- Step 7. Pathway-level comparison
- Step 8. Bubble plots and chord diagrams
- Step 9. Gene expression validation
- Step 10. Correlation analysis 

#### H3.0 Expected output
- **RDS files**
    - `results/cellchat.HC.rds`
    - `results/cellchat.PD.rds`
- **Figures**
    - `cellchat.p1`: Barplots of interaction counts and weights (HC vs PD)
    - `cellchat.p2`: Heatmaps of communication networks
    - `cellchat.p3.1-3.2`: Differential interaction networks
    - `cellchat.p4`: Compare signaling pathway activities between HC and PD
    - `cellchat.p5.1-5.2`: Scatter plots showing signaling role analysis
    - `cellchat.p6.1-6.4`: Chord diagrams for specific signaling pathways in PD
    - `cellchat.p7.1-7.4`: Violin plots showing gene expression of pathway genes
    - `cellchat.p8.1-8.4`: Network role analysis for selected pathways
    - `cellchat.p12-14`: Correlation scatterplot with marginal histograms


### SCENIC (R & Python Analysis)
#### H1.0 Input data
- **From Seurat (R export)**
    - `seurat_counts.mtx` → Raw count matrix (genes × cells, sparse format)
    - `seurat_gene_names.csv` → Gene name list
    - `seurat_metadata.csv` → Cell metadata (e.g., cluster, sample ID, disease status)
    - `seurat_pca.csv` → PCA embeddings
    - `counts.csv` → Alternative full expression matrix (dense format, genes × cells)
        
- **Reference databases (for PySCENIC)**
    - `hs_hgnc_tfs.txt` → List of human transcription factors
    - `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` → Motif annotation
#### H2.0 Code execution order
- **R preprocessing (SCENIC.R)**
    - Export Seurat object into SCENIC-compatible format (`counts.mtx`, metadata, gene names, PCA).
    - Example:
        `seurat_to_adata(obj.fib, Dimension="UMAP", path=".../H9. SCENIC/")`
- **Convert to Loom format (Python, `trans_loom.py`)**
    `x = sc.read_csv("counts.csv") row_attrs = {"Gene": np.array(x.var_names)} col_attrs = {"CellID": np.array(x.obs_names)} lp.create("sce.loom", x.X.transpose(), row_attrs, col_attrs)`
- **Run PySCENIC pipeline (pyscenic.py)**
    - Step 1: Build co-expression modules (GRN inference).
    - Step 2: Motif enrichment pruning (using TF + motif databases).
    - Step 3: Regulon scoring → output `sce_SCENIC.loom`.
- **Post-analysis (R + Python visualization)**
    - In R (`SCENIC` package): calculate **RSS scores**, visualize regulon activity, extract TF–target networks.
    - In Python (`scanpy`):
        `sce = sc.read("sce_test.h5ad") sc.pl.umap(sce, color="annotation", legend_loc="on data", size=8)`
#### H3.0 Expected output
1. **Loom / H5AD objects**
    - `sce.loom` → raw loom file (counts + metadata).
    - `sce_SCENIC.loom` → processed loom file with regulon AUC activity
    - `sce_test.h5ad` → AnnData object for downstream Scanpy visualization.
2. **Regulon results**
    - `sce.adj.csv` → adjacency matrix (gene–gene co-expression edges).
    - `sce.regulons.csv` / `regulons.csv` → TF–target mapping and regulon scores.