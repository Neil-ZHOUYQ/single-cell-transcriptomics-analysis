install.packages("Rogue")
install.packages("clustree")
install.packages("harmony")

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

# PCA linear dimensionality-reduction 
# standards of selecting PCs
# PCA accumulated contribution > 90%
# PC's standard deviation contribution < 5%
# continuous PCs have difference less than 0.1%






#Directory
setwd("E:/learning_materials/sc/8harmony/8harmony")

#Read
pbmc = readRDS("E:/learning_materials/sc/7Normalize/7Normalize/4pbmc_Normalize.rds")

#Set active matrix and grouping information
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "SCT"
DefaultAssay(pbmc)

table(Idents(pbmc))
Idents(pbmc) = "orig.ident"
table(Idents(pbmc))

#Top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10) # SCTransform found variable features in previous step and saved in pbmc@assays$SCT@var.feature

#Plot
pdf(file="1.pdf",width=7,height=6)
VariableFeaturePlot(object = pbmc)
dev.off()

pdf(file="2.pdf",width=7,height=6)
LabelPoints(plot = VariableFeaturePlot(object = pbmc), points = top10, repel = TRUE)  # repel = True, making sure gene names don't overlap
dev.off()

#PCA
pbmc <- RunPCA(pbmc, verbose = F)   # verbose = F , don't print on console, pbmc@reductions$pca

#Principal component analysis plot
pdf(file="3.pdf",width=7,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#Plot genes associated with each PCA component
pdf(file="4.pdf",width=10,height=9)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)      # visualize the weights of 20 genes in pc from 1 to 4
dev.off()

#Principal component analysis heatmap
pdf(file="5.pdf",width=10,height=9)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2) # 250 top-scored and bottom-scored cells; 15 top-weighed and bottom-weighted genes 
dev.off()

#Select appropriate PCs
#Cumulative contribution of principal components > 90%, choose the elbow point
pdf(file="6.pdf",width=7,height=6)
ElbowPlot(pbmc, ndims = 50)                 # The decrease of standard variation represents the variance PCs explain
dev.off()

#Determine the percentage of each PC  
pct <- pbmc [["pca"]]@stdev / sum( pbmc [["pca"]]@stdev) * 100  # The contribution of standard variation each PC contributes
pct

#Calculate the cumulative percentage of each PC                      # Cumulative sum
cumu <- cumsum(pct)
cumu

#Set PCs
pcs = 1:40



# harmony Batch Effect Correction
# each iteration, pull the centroid of sample A and sample B to global centroid of that cell type
pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20, dims.use = pcs)

table(pbmc@meta.data$orig.ident)

#Select appropriate resolution
#Run through resolutions from 0.1-2           
seq = seq(0.1,2,by=0.1)
pbmc <- FindNeighbors(pbmc, dims = pcs) 
for (res in seq){
 pbmc = FindClusters(pbmc, resolution = res)
}

#Plot
p1 = clustree(pbmc,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("SCT_sun_res.png", p, width = 30, height = 14)       # select the resolution which gives clear enough clustering

#Dimensionality reduction clustering
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 1)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
#Dimensionality reduction clustering
#pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = pcs) %>% FindClusters(resolution = 1)   # ?

#pbmc <- RunUMAP(pbmc, reduction = "pca", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "pca")

colnames(pbmc@meta.data)
#Plot
pdf(file="7.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "umap", label = T)      #resolution = 1, 23 clusters
dev.off() 

pdf(file="8.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

pdf(file="9.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

pdf(file="10.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "tsne", label = T)
dev.off()

pdf(file="11.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne",label = F,group.by = "orig.ident")
dev.off()

pdf(file="12.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()

#Save
saveRDS(pbmc,"5pbmc_UMPA.TSNE.rds")