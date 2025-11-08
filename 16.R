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

#Directory
setwd("F:/4scRNA/1code/16ssGSEA")

#Read in
pbmc <- readRDS("F:/4scRNA/1code/9SingleR/6pbmc_SingleR.rds")
DefaultAssay(pbmc) <- "RNA"
#Normalize
pbmc <- NormalizeData(pbmc)
dim(pbmc)

#Get expression matrix
expr = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "data"))
#Get average expression
#expr2 = as.matrix(AverageExpression(pbmc, assays = "RNA", layer = "data", group.by = "Is_Double")[[1]])
dim(expr)

#Filter
expr <- expr[rowSums(expr)>0,]
dim(expr)

#Read in gene sets, put the gene sets you want to score into genesets.xlsx
#The second column in xlsx needs to be a description of the gene set, if not, set to na
geneSets2 = read_excel("genesets.xlsx",col_names = F)

#Export as gmt file
write.table(geneSets2,file = "my_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

#Read gmt file
geneSets=getGmt("my_genesets.gmt", geneIdType=SymbolIdentifier())

#Analyze method=c("gsva", "ssgsea", "zscore", "plage")
gsvaResult=gsva(expr, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

#If the following error occurs, it is due to the GSVA package being updated to the new version 1.52
#Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').
#Use the following command instead of the one above
#It is no longer possible to download and install the old version from either github or bioconductor, to use the old version, you can only copy it from an R environment that has the old version of the GSVA package

#Analyze
#gsvapar <- gsvaParam(expr, geneSets,kcdf='Gaussian',absRanking=TRUE)
#gsvaResult <- gsva(gsvapar)

#Normalize the score
normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)

#Export
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="ssgseaOut.txt", sep="\t", quote=F, col.names=F)

#Add to meta.data
pbmc <- AddMetaData(pbmc, metadata = t(gsvaResult))

#Plot
p <- FeaturePlot(pbmc, features = rownames(gsvaResult),reduction = "umap", ncol = 3)
ggsave("1umap.pdf", p, width = 15, height = 5)
p2 <- DimPlot(pbmc, group.by = "Is_Double", label = T,reduction = "umap")
ggsave("2double.pdf", p2, width = 6, height = 5)

#Get group information
group <- data.frame(group1 = pbmc$Is_Double, row.names = colnames(pbmc))

#Sort by group
gsvaResult1 = gsvaResult[,group == "Singlet"]
gsvaResult2 = gsvaResult[,group == "Doublet"]
gsvaResult = cbind(gsvaResult1,gsvaResult2)
gsvaResult2 = gsvaResult

#Get the number of samples in each group
conNum=length(group[group=="Singlet"])
treatNum=length(group[group=="Doublet"])
Type=c(rep(1,conNum), rep(2,treatNum))
Type2 <- data.frame(group1 = c(rep("Singlet",conNum),rep("Doublet",treatNum)), 
                    row.names = colnames(gsvaResult2))

#Merge
gsvaResult = cbind(t(gsvaResult),Type)

#Heatmap
pdf(file="3heatmap.pdf", width=8, height=3)
pheatmap(gsvaResult2, 
         annotation=Type2, 
         color = colorRampPalette(c(rep("DodgerBlue1",5), "white", rep("Firebrick2",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#Extract significantly different pathways, if no need to filter significant pathways, change 0.05 below to 1
sigGene=c()
for(i in colnames(gsvaResult)[1:(ncol(gsvaResult)-1)]){
  test=wilcox.test(gsvaResult[,i] ~ gsvaResult[,"Type"])
  pvalue=test$p.value
  if(pvalue<1){ #can be modified
    sigGene=c(sigGene, i)
  }
}

#Extract the above significant pathways, remove non-significant ones
hmExp=gsvaResult[,sigGene]
Type=c(rep("Singlet",conNum),rep("Doublet",treatNum))
names(Type)=rownames(gsvaResult)
Type=as.data.frame(Type)
hmExp = t(hmExp)

#ggboxplot plot
#Convert data to ggplot2 input file
hmExp2 = t(hmExp)
hmExp2 = cbind(hmExp2,Type)
rt=melt(hmExp2,id.vars=c("Type"))
colnames(rt)=c("Type","Genesets","Expression")

#Set comparison groups
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=c("Singlet","Doublet"))

#Draw boxplot
boxplot=ggboxplot(rt, x="Genesets", y="Expression", fill="Type",
                  xlab="",
                  ylab="Score",#Modify Y-axis title name
                  legend.title="Type",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+#Modify color
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
      panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
       legend.title = element_text(face ="bold.italic",size = 13)
  )

#Output image
pdf(file="4ggboxplot.pdf", width=8, height=6)
print(boxplot)
dev.off()