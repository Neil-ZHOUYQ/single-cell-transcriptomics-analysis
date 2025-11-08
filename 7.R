#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#Directory
setwd("E:/learning_materials/sc/7Normalize/7Normalize")

#Read
pbmc = readRDS("E:/learning_materials/sc/6CellCycle/6CellCycle/3pbmc_CellCycle.rds")

# Assumes that each cell initially contains the same number of RNA molecules.
# For each cells, the counts of each gene: (gene count UMI/total UMIs) * scale_factor, then ln(1+x)
# Results are stored in pbmc@assays$RNA@layers$data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)



# Since single-cell matrices are sparse, many gene expression values are almost 0
# Helps to highlight biological signals in single-cell datasets
# Find the top 2000 highly variable genes
# vst: variance-stabilizing transformation
# results are saved in pbmc@assays$RNA@meta.data

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# ScaleData() 1. regress out 2. scaling
# The effect of cell cycle on gene expression is confounding variable but not technical noise. 
# model of cell cycle genes expression: lm(Gene_A ~ S.Score + G2M.Score)
# Residual = observed - predicted, residual is true gene expression signal without pollution
# Different genes' expression values varies in different scale (0-8, 0-2), but they share same importance in regulation 
# For each gene, use (X-mean)/sd(X) to scale and normalize the expression matrix

# Shift the expression of each gene so that the mean expression across cells is 0
# Scale the expression of each gene so that the variance across cells is 1
# (X - mean)/sd(X)
# Give equal weight in downstream analyses.
# Remove the effect of cell cycle here
# Results are stored in pbmc@assays$RNA@layers$scale.data
pbmc <- ScaleData(pbmc,vars.to.regress = c("S.Score", "G2M.Score"))
#pbmc <- ScaleData(pbmc,vars.to.regress = c("S.Score", "G2M.Score"),features = rownames(pbmc)) #scale all the genes not only high variance 2000 

# SCT, equivalent to replacing the above three functions: NormalizeData, FindVariable, ScaleData
# Find 3000 highly variable genes

# In conventional analysis, using a small number of PCs can focus on key biological differences
# while not introducing more technical differences, which is a conservative approach.
# It loses some biological difference information but is relatively safe in conventional methods.
# But SCT's normalization and standardization are done well,
# inputting more PCs can extract more biological differences, while ensuring no technical errors are introduced.
# SCT believes: the newly added 1000 genes contain subtle biological differences not detected before.
# Moreover, even if all genes are used for downstream analysis, the results are similar to SCT's results

# Results are stored in pbmc@assays$SCT
pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score"))

#Default matrix
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "RNA"
DefaultAssay(pbmc)

#Subsequent analysis will continue using SCT as an example
DefaultAssay(pbmc) = "SCT"
DefaultAssay(pbmc)

#Save
saveRDS(pbmc,"4pbmc_Normalize.rds")