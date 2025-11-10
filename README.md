# Single-Cell RNA-seq Analysis Pipeline
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  
A comprehensive collection of R and Python scripts for single-cell RNA sequencing (scRNA-seq) data analysis, from raw data processing to advanced downstream analyses.

## Overview

This repository contains a complete workflow for analyzing single-cell RNA-seq data, including quality control, normalization, clustering, cell type annotation, trajectory analysis, cell-cell communication, and various enrichment analyses.

## Dependencies

### R Packages
- **Core**: `Seurat`, `data.table`, `stringr`, `tibble`
- **Visualization**: `ggplot2`, `patchwork`, `pheatmap`, `ComplexHeatmap`, `ggheatmap`
- **Analysis Tools**: `DoubletFinder`, `harmony`, `SingleR`, `monocle`, `CellChat`, `SCENIC`, `copykat`
- **Enrichment**: `clusterProfiler`, `GSVA`, `GSEABase`, `AUCell`
- **Utilities**: `clustree`, `Rogue`, `dplyr`, `tidyverse`

### Python Packages
- `pyscenic==0.12.1`
- `loompy`
- `scanpy==1.10.0`
- `numpy==1.23.5`
- `pandas==1.5.3`

### Command-line Tools
- `CellRanger`
- `prefetch` (SRA Toolkit)
- `parallel-fastq-dump`

## Pipeline Workflow

### 1. Data Acquisition and Preprocessing (`1.sh`)
**Purpose**: Download SRA data and perform alignment with CellRanger

**Key steps**:
- Download raw sequencing data from SRA using `prefetch`
- Convert SRA files to FASTQ format using `parallel-fastq-dump`
- Merge multiple sequencing runs for the same sample
- Align reads to reference genome using CellRanger
- Generate gene expression matrices

**Usage**: Modify SRR accession numbers and file paths, then run the shell script.

---

### 2. Data Import (`3.R`)
**Purpose**: Read single-cell data from various file formats

**Supports**:
- 10X Genomics format (matrix.mtx, genes.tsv, barcodes.tsv)
- H5 files
- R data files (RDS/RData)
- Plain text files (TXT/CSV)

**Key functions**:
- `Read10X()` - Read 10X format data
- `Read10X_h5()` - Read H5 files
- `CreateSeuratObject()` - Create Seurat objects
- `merge()` - Combine multiple samples
- `JoinLayers()` - Join data layers

---

### 3. Quality Control (`4.R`)
**Purpose**: Filter low-quality cells and visualize QC metrics

**QC metrics**:
- Number of detected genes per cell (`nFeature_RNA`)
- Total UMI counts per cell (`nCount_RNA`)
- Mitochondrial gene percentage (`percent.mt`)
- Hemoglobin gene percentage (`percent.HB`)

**Output**: Violin plots before/after QC, filtered Seurat object

---

### 4. Doublet Detection (`5.R`)
**Purpose**: Identify and remove doublets using DoubletFinder

**Key steps**:
- Normalization, feature selection, and scaling
- PCA and UMAP dimensionality reduction
- Optimal pK parameter selection via parameter sweep
- Homotypic doublet proportion estimation
- Doublet classification

**Output**: Doublet scores and classifications added to metadata

---

### 5. Cell Cycle Scoring (`6.R`)
**Purpose**: Score cells for cell cycle phase

**Features**:
- Assign S phase and G2M phase scores
- Classify cells into G1, S, or G2M phases
- Visualize cell cycle distribution

**Output**: Cell cycle phase annotations in metadata

---

### 6. Normalization and Scaling (`7.R`)
**Purpose**: Normalize, identify variable features, and scale data

**Methods**:
- **LogNormalize**: Standard log-normalization
- **FindVariableFeatures**: Identify highly variable genes
- **ScaleData**: Z-score normalization with optional regression
- **SCTransform**: Advanced normalization method

**Key parameters**:
- Regress out cell cycle effects
- Scale factor: 10,000
- Variable features: 2,000-3,000 genes

---

### 7. Dimensionality Reduction and Batch Correction (`8.R`)
**Purpose**: PCA, batch effect correction, and clustering

**Key steps**:
- Principal Component Analysis (PCA)
- Harmony batch effect correction
- Optimal resolution selection using clustree
- UMAP and t-SNE visualization
- Neighborhood graph construction and clustering

**Output**: Dimensionally reduced data with cluster assignments

---

### 8. Automated Cell Type Annotation (`9.R`)
**Purpose**: Automated cell type annotation using SingleR

**Features**:
- Reference-based cell type prediction
- Uses Human Primary Cell Atlas or custom references
- Calculates Spearman correlation between query and reference
- Provides confidence scores for predictions

**Output**: Cell type predictions added to metadata

---

### 9. Manual Cell Type Annotation (`10.R`)
**Purpose**: Manual annotation using marker genes

**Key steps**:
- Find marker genes for each cluster (`FindAllMarkers`)
- Visualize known cell type markers (FeaturePlot, DotPlot, VlnPlot)
- Iterative manual annotation based on marker expression
- Generate heatmaps of top marker genes
- Cell proportion analysis across samples

**Markers used**:
- Immune: PTPRC
- Epithelial: EPCAM, KRT8, KRT18, KRT19
- Stromal: PECAM1, VWF, ACTA2
- T cells: CD3D, CD3E, CD8A, CD4
- Myeloid: CD14, CD68, CD163

---

### 10. Trajectory Analysis (`11.R`)
**Purpose**: Pseudo-temporal ordering with Monocle 2

**Key steps**:
- Create CellDataSet object
- Select ordering genes (high variance)
- Dimensionality reduction with DDRTree
- Order cells along trajectory
- Identify branch-dependent genes (BEAM analysis)
- Find pseudotime-dependent genes

**Output**: 
- Trajectory plots colored by state, pseudotime, cell type
- Pseudotime-dependent gene lists
- Branch-specific gene expression

---

### 11. Cell-Cell Communication (`12.R`)
**Purpose**: Infer cell-cell interactions using CellChat

**Key steps**:
- Load ligand-receptor database (CellChatDB)
- Identify overexpressed ligands and receptors
- Compute communication probability
- Infer signaling pathways
- Network visualization (circle plots, chord diagrams, heatmaps)

**Output**:
- Cell-cell communication networks
- Signaling pathway activity
- Ligand-receptor pair interactions
- Sender/receiver/mediator/influencer roles

---

### 12. Transcription Factor Activity Analysis (`13_1.R`, `create_loom.py`, `run_pySCENIC.sh`)
**Purpose**: Identify active transcription factor regulons using SCENIC

**Workflow**:
1. **R**: Export expression matrix from Seurat
2. **Python** (`create_loom.py`): Convert to loom format
3. **Shell** (`run_pySCENIC.sh`): Run pySCENIC pipeline
   - GRN inference (GRNBoost2)
   - Regulon prediction (cisTarget)
   - AUCell scoring
4. **R**: Import results, calculate RSS, differential analysis

**Output**: 
- Cell type-specific transcription factors
- Regulon activity scores
- Binary regulon activity

---

### 13. Data Integration (`14.R`)
**Purpose**: Compare multiple batch correction methods

**Methods tested**:
- CCA (Canonical Correlation Analysis)
- RPCA (Reciprocal PCA)
- JointPCA
- Harmony
- FastMNN (optional)

**Features**: Works with both RNA and SCT assays

---

### 14. Tumor Cell Identification (`15.R`)
**Purpose**: Distinguish malignant from normal cells using CopyKAT

**Key steps**:
- Extract tumor and reference cells
- Infer copy number variations
- Classify cells as diploid, aneuploid, or not defined
- Visualize results on UMAP

**Output**: Malignancy predictions per cell

---

### 15. Gene Set Enrichment Analysis - ssGSEA (`16.R`)
**Purpose**: Score cells for custom gene set activity

**Method**: Single-sample GSEA (ssGSEA)

**Key steps**:
- Load custom gene sets from Excel file
- Run ssGSEA on expression matrix
- Normalize scores
- Add to Seurat metadata
- Compare between groups (wilcox test)

**Output**: Heatmaps, UMAP plots, boxplots of pathway scores

---

### 16. Gene Set Enrichment Analysis - AUCell (`17.R`)
**Purpose**: Alternative gene set scoring using AUCell

**Method**: Area Under the Curve (AUCell)

**Advantages**: 
- Better for sparse data
- Faster computation
- Threshold-based binary classification

**Output**: Similar to ssGSEA (scores, plots, comparisons)

---

### 17. GO and KEGG Enrichment (`18.R`)
**Purpose**: Functional enrichment analysis of marker genes

**Features**:
- Gene Ontology (GO) enrichment (BP, CC, MF)
- KEGG pathway enrichment
- Symbol to Entrez ID conversion
- Multiple visualization types (bar plots, bubble plots)

**Analysis**: Performed separately for each cell type cluster

---

### 18. Gene Set Enrichment Analysis - GSEA (`19.R`)
**Purpose**: Pre-ranked GSEA analysis

**Features**:
- Ranks genes by log fold change
- Tests enrichment of pre-defined gene sets
- Identifies upregulated and downregulated pathways
- Generates enrichment plots

**Gene sets**: KEGG pathways (customizable)

---

## Directory Structure

```
.
├── 1.sh                    # Data download and CellRanger
├── 3.R                     # Data import
├── 4.R                     # Quality control
├── 5.R                     # Doublet detection
├── 6.R                     # Cell cycle scoring
├── 7.R                     # Normalization
├── 8.R                     # Batch correction and clustering
├── 9.R                     # SingleR annotation
├── 10.R                    # Manual annotation
├── 11.R                    # Trajectory analysis
├── 12.R                    # Cell-cell communication
├── 13_1.R                  # SCENIC (R portion)
├── 14.R                    # Data integration
├── 15.R                    # CopyKAT
├── 16.R                    # ssGSEA
├── 17.R                    # AUCell
├── 18.R                    # GO/KEGG enrichment
├── 19.R                    # GSEA
├── create_loom.py          # SCENIC loom creation
└── run_pySCENIC.sh         # SCENIC pipeline
```

## Usage Notes

1. **Sequential execution**: Scripts are numbered to indicate the typical analysis order
2. **File paths**: Update `setwd()` commands and file paths to match your system
3. **Species**: Most scripts are configured for human data (use `org.Hs.eg.db`, `hg38`, etc.). For mouse data, change to `org.Mm.eg.db`, `mm10`, etc.
4. **Computational resources**: Adjust `cores`, `workers`, and `threads` parameters based on available resources
5. **Intermediate files**: Each script saves `.rds` files that can be loaded by subsequent scripts

## Typical Workflow

A standard analysis typically follows this order:

```
1.sh → 3.R → 4.R → 5.R → 6.R → 7.R → 8.R → 9.R/10.R → downstream analyses (11.R-19.R)
```

Downstream analyses (scripts 11-19) can be run independently based on research questions.



## Citation

If you use these scripts, please cite the relevant tools:
- **Seurat**: Hao et al., Cell 2021
- **DoubletFinder**: McGinnis et al., Cell Systems 2019
- **Harmony**: Korsunsky et al., Nature Methods 2019
- **SingleR**: Aran et al., Nature Immunology 2019
- **Monocle**: Trapnell et al., Nature Biotechnology 2014
- **CellChat**: Jin et al., Nature Communications 2021
- **SCENIC**: Aibar et al., Nature Methods 2017
- **CopyKAT**: Gao et al., Nature Biotechnology 2021



## Contact

For questions or issues, please open an issue on GitHub.
