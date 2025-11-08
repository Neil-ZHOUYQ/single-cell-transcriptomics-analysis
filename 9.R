install.packages("BiocManager")
BiocManager::install("SingleR")

n#Load
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
#SingleR version 2.4.1

#Directory
setwd("E:/learning_materials/sc/9SingleR/9SingleR")

#Read
pbmc = readRDS("E:/learning_materials/sc/8harmony/8harmony/5pbmc_UMPA.TSNE.rds")

#Load reference database
load("ref_Human_all.RData")

#Get expression matrix
#testdata = GetAssayData(object = pbmc@assays$RNA, layer = "data")
testdata = GetAssayData(object = pbmc@assays$SCT, layer = "data")

#Get clusters
clusters <- pbmc@meta.data$seurat_clusters

#Main and fine cell populations of the reference database
table(ref_Human_all@colData@listData[["label.main"]])
table(ref_Human_all@colData@listData[["label.fine"]])

#Run singleR
# only use marker genes
# calculate the Spearman's Rank Correlation of each cluster and cell type
cellpred <- SingleR(test = testdata, ref = ref_Human_all, clusters = clusters, assay.type.test = "logcounts", 
          labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")

#Get the cell type for each cluster
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#Add to the seurat.metadata object
pbmc@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
 pbmc@meta.data[which(pbmc$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}

#Plot
p = plotScoreHeatmap(cellpred)
ggsave("1.pdf", p, width = 12, height = 5)

p1 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "umap")
p <- p1 | p2
ggsave("2.pdf", p, width = 12, height = 5)

#Save
saveRDS(pbmc,"6pbmc_SingleR.rds")