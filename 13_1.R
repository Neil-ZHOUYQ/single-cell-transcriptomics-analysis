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

#Directory
setwd("E:/4scRNA/13pySCENIC")

#Read in
pbmc <- readRDS("pbmc_Fib.rds")

#Get expression matrix
matrix1 = as.matrix(pbmc@assays$RNA@counts)
dim(matrix1)
#Randomly select 500 cells to speed up runtime
# matrix1 <- matrix1[,sample(colnames(matrix1), 500)]
# dim(matrix1)

#Transpose, export
write.csv(t(matrix1),file = "pbmc_Fib_exp.csv")



# run run_pySCENIC.sh



#The following content runs in R
setwd("E:/4scRNA/13pySCENIC")

#Read in
sce_SCENIC <- open_loom("sce_SCENIC.loom")

#regulons (each regulon is 1 TF and its target genes, like a gene set)
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")

#regulons list
regulons <- regulonsToGeneLists(regulons_incidMat)

#AUCell matrix
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')

#Read in scRNA
sco <- readRDS("pbmc_Fib.rds")
colnames(sco@meta.data)

#RSS analysis, calculate regulon specificity score, find cell type-specific transcription factors
rss <- calcRSS(regulonAUC, sco$fib.celltype.main)
rss <- na.omit(rss)

#Scatter plot
rssPlot <- plotRSS(rss,zThreshold = 1,cluster_columns = FALSE,
                   order_rows = TRUE,
                   thr=0.01,
                   varName = "cellType",
                   col.low = '#330066',
                   col.mid = '#66CC66',
                   col.high = '#FFCC33')
pdf(file="1rssPlot.pdf", width=8, height=20)
print(rssPlot)
dev.off()

#Heatmap
rss_data <- rssPlot$plot$data
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
col_ann <- data.frame(group= colnames(rss_data))#Column annotation
rownames(col_ann) <- colnames(rss_data)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC","#330066","#FFCC33")
names(groupcol) <- colnames(rss_data)
col <- list(group=groupcol)
text_columns <- sample(colnames(rss_data),0)#Do not show column names
p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
pdf(file="2heatmap.pdf", width=8, height=20)
print(p)
dev.off()

#Plot
topN=5
plot <- list()
for(i in colnames(rss)){
  df.i <- data.frame(SpecificityScore=rss[,i], labels=rownames(rss))
  df.i <- arrange(df.i, desc(SpecificityScore)) %>% as.data.frame()
  df.i$Regulons <- 1:nrow(df.i)
  df.i$color <- ifelse(df.i$Regulons <= topN, "red", "gray")
  df.i$labels[df.i$Regulons > topN] <- NA
  p <- ggplot(df.i, aes(Regulons, SpecificityScore)) +
    geom_point(size=3, color = df.i$color) + 
    geom_text_repel(aes(label=labels), size=4) + 
    scale_x_continuous(limits = c(0, 200)) +
    ggtitle(i) + xlab("Regulons Rank") + ylab("Specificity Score") + 
    theme_bw(base_size = 12) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
  plot[[i]] <- p
}
p <- wrap_plots(plot, ncol = 3)
ggsave("3RSS.pdf", p, width = 9, height = 8, limitsize = F)

#Calculate activity threshold
regulonAUC = assay(regulonAUC)
bin.T <- AUCell_exploreThresholds(regulonAUC,
                                  smallestPopPercent=0.25,
                                  assignCells=TRUE, 
                                  plotHist=FALSE,
                                  verbose=FALSE)

#Binary transformation
regulonBin <- lapply(rownames(regulonAUC), function(reg){
  as.numeric(colnames(regulonAUC) %in% bin.T[[reg]][["assignment"]])
})
regulonBin <- do.call("rbind", regulonBin)
dimnames(regulonBin) <- list(rownames(regulonAUC), colnames(regulonAUC))

#Differential analysis, Seurat
sco[["Regulon"]] <- CreateAssayObject(counts = regulonAUC)
sco[["binRegulon"]] <- CreateAssayObject(counts = regulonBin)
DefaultAssay(sco) <- "Regulon"
sco <- ScaleData(sco, features = rownames(sco))
Idents(sco) <- "fib.celltype.main"
deg <- FindAllMarkers(sco, only.pos = T, logfc.threshold = 0)
top <- group_by(deg, cluster) %>% top_n(10, avg_log2FC) %>% pull(gene) %>% unique()
p <- DoHeatmap(sco, features = top, label = F)
ggsave("4DEG.pdf", p, width = 12, height = 6.5, limitsize = F)

#Regulon overview
rn  = "OSR1(+)"
tf = "OSR1"
gp = "fib.celltype.main"
p1 <- DimPlot(sco, reduction = "umap", group.by = "fib.celltype.main", label = T) + NoLegend()
DefaultAssay(sco) <- "binRegulon"
p2 <- FeaturePlot(sco, reduction = "umap", features = rn) + ggtitle(paste0(rn,"_binRAS"))
DefaultAssay(sco) <- "Regulon"
p3 <- FeaturePlot(sco, reduction = "umap", features = rn) + ggtitle(paste0(rn,"_RAS"))
DefaultAssay(sco) <- "RNA"
p4 <- FeaturePlot(sco, reduction = "umap", features = tf) + ggtitle(paste0(tf,"_Expression"))
df1 <- data.frame(cluster=sco@meta.data[,gp,drop=T], auc=sco@assays$Regulon@counts[rn,], row.names = colnames(sco))
p5 <- ggboxplot(df1, x='cluster', y='auc', fill='cluster', bxp.errorbar = T, outlier.shape = NA) +
  ggtitle(paste0(rn,"_RAS")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5)) + RotatedAxis() 
df2 <- data.frame(cluster=sco@meta.data[,gp,drop=T], expr=sco@assays$RNA@data[tf,], row.names = colnames(sco))
p6 <- ggboxplot(df2, x='cluster', y='expr', fill='cluster', bxp.errorbar =T, outlier.shape = NA) +
  ggtitle(paste0(tf,"_Expression")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5)) + RotatedAxis() 
p <- (p1|p3|p4)/(p2|p5|p6)
ggsave(paste0(tf,"_overview.pdf"), p, width = 16, height = 9)