
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)


#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#Directory
setwd("E:/learning_materials/sc/5double/5double")

#Read
pbmc = readRDS("E:/learning_materials/sc/4QC/4QC/1pbmc_qc.rds")

#%>% (forward-pipe operator) is the most commonly used operator,
#it passes the data or expression on the left side to the function call or expression on the right side for execution,
#allowing continuous operations like a chain.



#Data normalization, finding variable features, data scaling
pbmc = pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

#PCA, UMAP, finding neighbors, clustering

pbmc = pbmc %>% 
  RunPCA() %>%    
  RunUMAP(dims = 1:30) %>%  
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.1)

#pbmc@reductions$pca
#pbmc@reductions$umap
#pbmc@reductions$tsne
#pbmc@graphs$RNA_snn
#pbmc@meta.data$seurat_clusters



DimPlot(pbmc, reduction = "umap", label = TRUE)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident")
# Plot the expression of the "CD3D" gene
# The darker the color (e.g., yellow/red), the higher the CD3D expression in that cell
FeaturePlot(pbmc, reduction = "umap", features = "CD3D")





# randomly pick 2 cells to sythesize "artificial doublets"
# pANN(Proportion of artificial nearest neighbors):Given a neighborhood size of a cell, calculate the proportion of artificial cells in it
# pk value: neighborhood size
# pk value is perfect when real singlets have low pANN and artificial doublets have high pANN. At this time, a Bimodal distribution of pANN value is shown,
# BCmetric: meature the Bimodal Distribution. Higher the BCmetric, better the distribution and so is the pk value.


#First, obtain the optimal pK value
#pK represents the neighborhood size
sweep.res.list <- paramSweep(pbmc, PCs = 1:30, sct = FALSE)   # try different pk value, calculate the pANN of each pk
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)     # calculate BCmetric(bimodality coefficient)
bcmvn <- find.pK(sweep.stats)                                 # bcmvn is a form of BCmetric with paired pk value
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()




#Then estimate the proportion of homotypic doublets in the doublet population
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)   # according to cell types in metadata, the proportion of homologous doublets in all doublets are estimated
print(homotypic.prop)


#The doublet proportion is around 7%
nExp_poi <- round(0.07*nrow(pbmc@meta.data))   # 7% is from 10X statistics, the prediction of all doublets proportion    
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # predicted heterologous cells number


# The number of simulated artificial doublets. Different values have little effect on the identification results, default is 0.25
# give total number of heterologous cells and neighbourhood size, label the cells that are considered as heterologous
pbmc <- doubletFinder(pbmc, PCs = 1:30, 
                      pN = 0.25,    # the proportion of artificial doublets in all cells
                      pK = pk_best,  # best performance pk value
                      nExp = nExp_poi.adj, #  predicted heterologous cells number
sct = FALSE)


#Change column names to "Double_score" and "Is_Double"
colnames(pbmc@meta.data)

colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))-1] <- "Double_score"
colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))] <- "Is_Double"

#View DoubletFinder analysis results
head(pbmc@meta.data[, c("Double_score", "Is_Double")])

#Plot tsne plot of DoubletFinder classification
DimPlot(pbmc, reduction = "tsne", group.by = "Is_Double")

#Plot violin plot of doublet classification
VlnPlot(pbmc, group.by = "Is_Double", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0, ncol = 2)

#Filter non-singlet data
#pbmc <- subset(pbmc, Is_Double == "Singlet")

#Save
saveRDS(pbmc,"2pbmc_double.rds")