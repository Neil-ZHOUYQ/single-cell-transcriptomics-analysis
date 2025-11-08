# install.packages("devtools")
devtools::install_github("aertslab/SCENIC")
devtools::install_github("rcastelo/GSVA")
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install RcisTarget (which you are missing) and AUCell (the third step of SCENIC)
BiocManager::install(c("org.Hs.eg.db"))
BiocManager::install("enrichplot", force = TRUE)
install.packages('GOplot')

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
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringi)  
library(GOplot)
library(tidyverse)

R.utils::setOption("clusterProfiler.download.method",'auto') #automatically choose the right method to download when clusterProfiler is running

#Directory
setwd("E:/learning_materials/sc/18GOKEGG/18GOKEGG")

#Read in
markers <- read_csv("markers.celltype.RNA.csv")

#/ is the directory operator
markers$cluster = gsub("Basal/intermediate","Basal.intermediate",markers$cluster)

#Filter
#markers <- subset(markers, p_val < 0.01 & avg_log2FC>1)

#According to cluster, sequentially extract genes from each cluster for GO, KEGG
for (cluster2 in unique(markers$cluster)) {
  
  #Get the gene names for each cluster
  input_gene = as.vector(subset(markers,cluster == cluster2)[,1])[[1]]
  
  #Convert gene symbols to gene ids
  #org.Hs.eg is for human species
  #https://www.jianshu.com/p/84e70566a6c6
  entrezIDs=mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA) #transform input_genes to entrezID through dictionary org.HS.egSYMOL2EG
  entrezIDs=as.character(entrezIDs)
  #Remove genes with NA gene id
  gene=entrezIDs[entrezIDs!="NA"]
  #Remove multiple IDs
  gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  
  #Filter conditions
  pvalueFilter=0.05
  qvalueFilter=1       
  
  if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}
  
  #GO
  kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
  GO=as.data.frame(kk)
  GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
  
  #Save
  write.table(GO, file=paste0(cluster2,".GO.txt"), sep="\t", quote=F, row.names = F)

  #Number to display
  showNum=10
  
  #Bar plot
  pdf(file=paste0(cluster2,".GObarplot.pdf"), width=10, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bar)
  dev.off()
  
  #Bubble plot
  pdf(file=paste0(cluster2,".GObubble.pdf"), width=10, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bub)
  dev.off()
  
  #KEGG
  kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
  KEGG=as.data.frame(kk)
  KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(input_gene[match(strsplit(x,"/")[[1]],as.character(entrezIDs))],collapse="/")))
  KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
  
  #Save
  write.table(KEGG, file=paste0(cluster2,".KEGG.txt"), sep="\t", quote=F, row.names = F)
  
  #Define the number of pathways to display
  showNum=20
  
  #Bar plot
  pdf(file=paste0(cluster2,".KEGGbarplot.pdf"), width=9, height=7)
  bub = barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
  print(bub)
  dev.off()
  
  #Bubble plot
  pdf(file=paste0(cluster2,".KEGGbubble.pdf"), width = 9, height = 7)
  bub = dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
  print(bub)
  dev.off()
}
