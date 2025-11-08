#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(Rogue)
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
library(copykat)

#Directory
setwd("E:/learning_materials/sc/15copykat/15copykat")

#Read
pbmc <- readRDS("E:/learning_materials/sc/9SingleR/9SingleR/6pbmc_SingleR.rds")

#Run on a single sample
pbmc = subset(pbmc, orig.ident %in% "2")
table(pbmc@meta.data$SingleR)

#Extract tumor-derived cells and some normal reference cells
pbmc1 <- subset(pbmc, SingleR %in% 'Epithelial_cells')
pbmc2 <- subset(pbmc, SingleR %in% c("T_cells","CMP","Macrophage"))                      # customize reference normal cells group

#Merge
pbmc <- merge(pbmc1, pbmc2)    # merge metadata, but still have pbmc@assays$RNA@layers$counts.pbmc1, pbmc@assays$RNA@layers$counts.pbmc2

#Join
DefaultAssay(pbmc) = "RNA"
pbmc = JoinLayers(pbmc)

#Extract expression matrix for analysis
counts <- as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "counts"))

#Set normal reference cells
ref <- colnames(pbmc2)

#Run
res <- copykat(rawmat=counts,ngene.chr=5,norm.cell.names=ref,sam.name="all",n.cores=10)
#No reference cells
#res <- copykat(rawmat=counts,ngene.chr=5,sam.name="all",n.cores=10)                    # automatic find normal cells  

#Save results
saveRDS(res,file="copykat.res.rds")

#Import copykat prediction results into seurat object
pbmc <- readRDS("F:/4scRNA/1code/9SingleR/6pbmc_SingleR.rds")

#Read in
malignant <- read.delim("all_copykat_prediction.txt")

#Convert to data frame
malignant <- data.frame(copykat.pred = malignant$copykat.pred, row.names = malignant$cell.names)

#Add to meta.data
pbmc <- AddMetaData(pbmc, metadata = malignant)
table(pbmc$copykat.pred)

#Modify
pbmc$copykat <- recode(pbmc@meta.data$copykat.pred,
                       "not.defined" = "diploid")
table(pbmc$copykat)

#Plot
p1 <- DimPlot(pbmc, group.by = "copykat", cols = c("red", "blue", "gray50"))
ggsave("2copykat_res.png", p1, width = 14, height = 7)

#Save
saveRDS(pbmc, file = "15pbmc_copykat.rds")
