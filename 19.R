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
library(copykat)
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringi)  
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')

#Directory
setwd("F:/4scRNA/1code/19GSEA")

#Read
pbmc = readRDS("F:/4scRNA/1code/10annotation/7pbmc_celltype.rds")

#Matrix
DefaultAssay(pbmc)

#Find markers for each cell type
markers3 <- FindAllMarkers(pbmc)

#/ is the directory operator
markers3$cluster = gsub("Basal/intermediate","Basal.intermediate",markers3$cluster)

#Save
write.csv(markers3,file="markers.celltype.RNA.all.csv")

#Read gene set file
gmt=read.gmt("c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
gmt[,1] = gsub("KEGG_", "", gmt[,1])

#According to cluster, sequentially extract genes from each cluster for GSEA
for (cluster2 in unique(markers3$cluster)) {
  
  #Extract markers for each cell subtype
  rt = subset(markers3,cluster == cluster2)
  
  #Sort by logFC
  rt = rt[order(rt[,"avg_log2FC"],decreasing=T),]
  
  #Get the logFC for each gene
  logFC = as.vector(rt[,"avg_log2FC"])
  names(logFC) = as.vector(rt[,"gene"])
  logFC[1:10]
  
  #Perform GSEA enrichment analysis on the ranked genes
  kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
  kkTab=as.data.frame(kk)
  kkTab=kkTab[kkTab$p.adjust<0.05,]
  write.table(kkTab,file=paste0(cluster2,".GSEA.result.KEGG.txt"),sep="\t",quote=F,row.names = F)
  
  #Plots for enrichment in the experimental group
  #Number of pathways to display
  termNum=5     
  kkUp=kkTab[kkTab$NES>0,]
  if (!is.na(kkUp[1,1])) {
    if(nrow(kkUp)>=termNum){showTerm=row.names(kkUp)[1:termNum]
    }else{
      showTerm=row.names(kkUp)[1:nrow(kkUp)]}
    gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in ",cluster2," up"), pvalue_table = T)
    pdf(file=paste0(cluster2,".up.GSEA.KEGG.pdf"), width=14, height=11)
    print(gseaplot)
    dev.off()
  }
  
  #Plots for enrichment in the normal group
  #Number of pathways to display
  termNum=5      
  kkDown=kkTab[kkTab$NES<0,]
  if (!is.na(kkDown[1,1])) {
    if(nrow(kkDown)>=termNum){showTerm=row.names(kkDown)[1:termNum]
    }else{
      showTerm=row.names(kkDown)[1:nrow(kkDown)]}
    gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in ",cluster2," down"))
    pdf(file=paste0(cluster2,".down.GSEA.KEGG.pdf"), width=14, height=11)
    print(gseaplot)
    dev.off()
  }
  
}