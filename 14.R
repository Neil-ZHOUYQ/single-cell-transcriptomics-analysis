#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(monocle)
library(tidyverse)
library(CellChat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(foreach)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggheatmap)
library(reshape2)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(SeuratWrappers)
library(SeuratData)
options(future.globals.maxSize = 1e9)

#Directory
setwd("F:/4scRNA/1code/14integrated")

#Read
pbmc = readRDS("F:/4scRNA/1code/7Normalize/4pbmc_Normalize.rds")

#Split
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)

# #Join
# pbmc = JoinLayers(pbmc)

#Normalize
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc,vars.to.regress = c("S.Score", "G2M.Score"))
# pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score"))

#PCA
pbmc <- RunPCA(pbmc, verbose = F)

#Set PCs
pcs = 1:40

#CCA
pbmc <- IntegrateLayers(
  object = pbmc, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE)

#RPCA
pbmc <- IntegrateLayers(
  object = pbmc, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE)

#JointPCA
pbmc <- IntegrateLayers(
  object = pbmc, method = JointPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.JointPCA",
  verbose = FALSE)

#Harmony
pbmc <- IntegrateLayers(
  object = pbmc, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE)

# #FastMNN
# Sys.setenv(LANGUAGE = "en")
# pbmc <- IntegrateLayers(
#   object = pbmc, method = FastMNNIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.mnn",
#   verbose = FALSE)


# #SCT
# #CCA
# pbmc <- IntegrateLayers(
#   object = pbmc,method = CCAIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE)
# #
# #RPCA
# pbmc <- IntegrateLayers(
#   object = pbmc, method = RPCAIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.rpca",
#   verbose = FALSE)
# #
# #Harmony
# pbmc <- IntegrateLayers(
#   object = pbmc, method = HarmonyIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.harmony",
#   verbose = FALSE)
# #
# #FastMNN
# # Sys.setenv(LANGUAGE = "en")
# # pbmc <- IntegrateLayers(
# #   object = pbmc, method = FastMNNIntegration,
# #   normalization.method = "SCT",new.reduction = "integrated.mnn",
# #   verbose = FALSE)

#Dimensionality reduction clustering
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = pcs)
pbmc <- FindClusters(pbmc, resolution = 1, cluster.name = "cca_clusters")
pbmc <- RunUMAP(pbmc, reduction = "integrated.cca", dims = pcs, reduction.name = "umap.cca")
pbmc <- RunTSNE(pbmc, reduction = "integrated.cca", dims = pcs, reduction.name = "tsne.cca")

#Plot
pdf(file="1cca.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "umap.cca", label = T)
dev.off()

pdf(file="2cca.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap.cca",label = F,group.by = "orig.ident")
dev.off()

pdf(file="3cca.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap.cca",label = F,group.by = "Is_Double")
dev.off()

pdf(file="4cca.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "tsne.cca", label = T)
dev.off()

pdf(file="5cca.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne.cca",label = F,group.by = "orig.ident")
dev.off()

pdf(file="6cca.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne.cca",label = F,group.by = "Is_Double")
dev.off()

DefaultAssay(pbmc)
#Join Layers
pbmc = JoinLayers(pbmc)

#Save
saveRDS(pbmc,"14pbmc_UMPA.TSNE.rds")
