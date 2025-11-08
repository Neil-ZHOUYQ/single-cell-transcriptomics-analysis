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

#Directory
setwd("E:/learning_materials/sc/10annotation/10annotation")

#Read
pbmc = readRDS("E:/learning_materials/sc/9SingleR/9SingleR/6pbmc_SingleR.rds")

DefaultAssay(pbmc)
#DefaultAssay(pbmc) = "RNA"
table(Idents(pbmc))
#Idents(pbmc) = "nCount_RNA"

#Find markers
pbmc.markers1 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1)
write.csv(pbmc.markers1,file="markers.1.SCT.csv")

DefaultAssay(pbmc) = "RNA"
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1)
write.csv(pbmc.markers2,file="markers.2.RNA.csv")
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "RNA"

#Create marker set
markers <- c("PTPRC", #immune
      "EPCAM", #epithelial
      "MME","PECAM1") #stromal

#Plot
p <- FeaturePlot(pbmc, features = markers, ncol = 3,reduction = "umap")
ggsave("1_1.pdf", p, width = 15, height = 10)

p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("2_1.pdf", p, width = 14, height = 6)

p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("3_1.pdf", p, width = 14, height = 6)

pdf(file="4_1.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "umap", label = T)
dev.off()

#First time
#immune epithelial stromal
# manually annotation
pbmc$celltype.1 <- recode(pbmc@meta.data$seurat_clusters,
              "0" = "epithelial",
              "1" = "immune",
              "2" = "epithelial",
              "3" = "3",
              "4" = "epithelial",
              "5" = "epithelial",
              "6" = "stromal",
              "7" = "epithelial",
              "8" = "epithelial",
              "9" = "epithelial",
              "10" = "immune",
              "11" = "epithelial",
              "12" = "epithelial",
              "13" = "immune",
              "14" = "immune",
              "15" = "epithelial",
              "16" = "stromal",
              "17" = "stromal",
              "18" = "stromal",
              "19" = "19",
V              "20" = "20",
              "21" = "epithelial")

#Check the subtypes of each cell population
table(pbmc@meta.data$celltype.1)

#Change color palette
Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
      '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
      '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
      '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
      '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
      '#968175')

#Plot
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.1",cols = Biocols)
ggsave("5_1.pdf", p, width = 7, height = 6)







#Second time
#Marker gene annotation
# get the marker genes from publications
markers <- c("ACTA2", #Fibroblast
      "PECAM1","VWF","ENG", #Endothelial
      "CMA1","MS4A2","TPSAB1","TPSB2", #Mast
      "AR","KRT19","KRT18","KRT8", #Luminal
      "TP63","KRT14","KRT5", #Basal/intermediate
      "LYZ","FCGR3A","CSF1R","CD68","CD163","CD14","UCHL1","HAVCR2", #Monolytic
      "PDCD1","CTLA4","CD8A","SELL","PTPRC","CD4","BTLA","IL2RA","IL7R","CCR7","CD28","CD27","SLAMF1","DPP4","CD7","CD2","CD3G","CD3E","CD3D") #T

#Plot
# p <- FeaturePlot(pbmc, features = markers, ncol = 3)
# ggsave("6.pdf", p, width = 15, height = 32)
p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("7_1.pdf", p, width = 14, height = 6)
p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("8_1.pdf", p, width = 14, height = 6)

#Fibroblast Endothelia Mast Luminal Basal/intermediate Monolytic T
pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
              "0" = "Luminal",
              "1" = "Mast",
              "2" = "Luminal",
              "3" = "Fibroblast",
              "4" = "Luminal",
              "5" = "5",
              "6" = "Endothelial",
              "7" = "Luminal",
              "8" = "8",
              "9" = "9",
              "10" = "Monolytic",
              "11" = "11",
              "12" = "12",
              "13" = "T",
              "14" = "Mast",
              "15" = "15",
              "16" = "Endothelial",
              "17" = "Fibroblast_17",
              "18" = "Endothelial",
              "19" = "Basal/interrmediate",
              "20" = "20",
              "21" = "Luminal")

#Check the subtypes of each cell population
table(pbmc@meta.data$celltype.main)
#Plot
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
ggsave("9_1.pdf", p, width = 7, height = 6)

p1 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.1")
ggsave("10_1.pdf", p1|p2, width = 15, height = 6)







#Third time
#Fibroblast Endothelia Mast Luminal Basal/intermediate Monolytic T
pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
              "0" = "Luminal",
              "1" = "Mast",
              "2" = "Luminal",
         _      "3" = "Fibroblast",
              "4" = "Luminal",
              "5" = "Luminal",
              "6" = "Endothelial",
              "7" = "Luminal",
              "8" = "Luminal",
              "9" = "Luminal",
              "10" = "Monolytic",
              "11" = "Luminal",
              "12" = "Luminal",
              "13" = "T",
              "14" = "Mast",
V              "15" = "Luminal",
          _     "16" = "Endothelial",
              "17" = "Fibroblast_17",
              "18" = "Endothelial",
       _        "19" = "Basal/interrmediate",
              "20" = "20",
              "21" = "Luminal")

#Check the subtypes of each cell population
table(pbmc@meta.data$celltype.main)

#Plot
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
ggsave("11_1.pdf", p, width = 7, height = 6)

p1 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.1")
ggsave("12_1.pdf", p1|p2, width = 15, height = 6)

p1 = DimPlot(pbmc, reduction = "tsne", label = T,group.by = "celltype.main")
p2 = DimPlot(pbmc, reduction = "tsne", label = T,group.by = "celltype.1")
ggsave("13_1.pdf", p1|p2, width = 15, height = 6)






#Fourth time
#Marker gene annotation
markers <- c("ACTA2", #Fibroblast
      "PECAM1","VWF","ENG", #Endothelial
      "MS4A2","TPSAB1","TPSB2", #Mast
      "AR","KRT19","KRT18","KRT8", #Luminal
      "TP63","KRT14","KRT5", #Basal/intermediate
      "LYZ","FCGR3A","CSF1R","CD68","CD163","CD14", #Monolytic
      "DPP4","CD7","CD2","CD3G","CD3E","CD3D") #T

#Plot
p <- DotPlot(pbmc, features = markers,group.by = "celltype.main") + RotatedAxis()
ggsave("14_1.pdf", p, width = 14, height = 6)

#Check the subtypes of each cell population
table(pbmc@meta.data$celltype.main)









#Grouping
table(Idents(pbmc))
Idents(pbmc) = "celltype.main"
table(Idents(pbmc))

#Matrix
DefaultAssay(pbmc)

#Find markers for each cell type
pbmc.markers3 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.3)
write.csv(pbmc.markers3,file="markers.celltype.RNA.csv")

#Extract top 10 markers
top10 <- pbmc.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Scale data for these genes
markers = as.data.frame(top10[,"gene"])
pbmc <- ScaleData(pbmc, features = as.character(unique(markers$gene)))

#Plot
p = DoHeatmap(pbmc,
     features = as.character(unique(markers$gene)),
     group.by = "celltype.main")
ggsave("15_1.pdf", p, width = 10, height = 9)

#Sample cells from each cell subtype (number of cells = smallest cell subtype)
allCells = names(Idents(pbmc))
allType = levels(Idents(pbmc))
choose_Cells = unlist(lapply(allType, function(x){
 cgCells = allCells[Idents(pbmc)== x ]
 cg=sample(cgCells,min(table(pbmc@meta.data$celltype.main)))
 cg
}))

#Extract
cg_sce = pbmc[, allCells %in% choose_Cells]
table(Idents(cg_sce))

#Plot
p = DoHeatmap(cg_sce,
       features = as.character(unique(markers$gene)),
       group.by = "celltype.main")
ggsave("16_1.pdf", p, width = 10, height = 9)




#Stacked bar chart
cell.prop<-as.data.frame(prop.table(table(pbmc@meta.data$celltype.main, pbmc@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

p = ggplot(cell.prop,aes(group,proportion,fill=cluster))+
 geom_bar(stat="identity",position="fill")+
 ggtitle("")+
 theme_bw()+
 theme(axis.ticks.length=unit(0.5,'cm'))+
 guides(fill=guide_legend(title=NULL))
ggsave("17_1.pdf", p, width = 10, height = 9)

#Save
saveRDS(pbmc,"7pbmc_celltype.rds")



# Clustering annotation results are not good, go back and re-adjust resolution
# Extract cell types for subdivision subtype definition
pbmc.T = pbmc[, Idents(pbmc) %in% c("T")]
# To further annotate the T cells 
# remove cell cycle effect, remove the effect only in T cells
# NormalizeData()
# FindVariableFeatures() to select high variable genes
# ScaleData according to T cells only
# run PCA to bulid new PCs
# RunHarmony()  remove batch effect only in T cells
# RunUMAP(), FindNeighbors(), FindClusters()