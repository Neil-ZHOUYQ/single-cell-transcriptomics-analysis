#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#Directory
setwd("E/learning_materials/sc/6CellCycle/6CellCycle")

#Read
pbmc = readRDS("E:/learning_materials/sc/5double/5double/2pbmc_double.rds")

#Cell cycle scoring
pbmc <- NormalizeData(pbmc)

#Get G2M phase related genes
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))

#Get S phase related genes
s_genes <- cc.genes$s.genes  
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))

#Cell cycle phase scoring
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)

colnames(pbmc@meta.data)
table(pbmc$Phase)

#Plot
DimPlot(pbmc,group.by = "Phase",reduction = "tsne")

#Save
saveRDS(pbmc,"3pbmc_CellCycle.rds")
