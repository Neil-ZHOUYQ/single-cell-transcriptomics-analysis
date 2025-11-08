#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

#Directory
setwd("E:/learning_materials/sc/4QC/4QC")

#Must use fread, other functions will take a long time
pbmc <-fread("GSE141445/data.raw.matrix.txt", sep="\t")
pbmc[1:5,1:5]

#Directly convert to row names
pbmc <- 	column_to_rownames(pbmc,"V1")
pbmc[1:5,1:5]

#Create seurat object
pbmc <- CreateSeuratObject(pbmc, min.features = 200, min.cells = 3)

#Join Layers
pbmc = JoinLayers(pbmc)

#View
table(pbmc@meta.data$orig.ident)

table(str_split(colnames(pbmc),'-',simplify = T)[,2])
#Add or modify meta.data information
pbmc <- AddMetaData(object = pbmc, 
                    metadata = str_split(colnames(pbmc),'-',simplify = T)[,2],   
                    col.name = "orig.ident") 

table(pbmc@meta.data$orig.ident)

#%in% checks if elements of the first vector are in the second vector, returns a boolean value.
table(pbmc@meta.data$orig.ident %in% c("1","2","3"))

dim(pbmc)

#Extract
pbmc = subset(pbmc,orig.ident %in% c("1","2","3"))

dim(pbmc)

#View
table(pbmc@meta.data$orig.ident)

#Quality control
#Mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Red blood cell percentage
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 

#View correlation
FeatureScatter(pbmc, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

#View QC metrics
#Set plotting elements
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
#Violin plot before QC
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin
#Save
ggsave("1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8) 
dim(pbmc)



#Set QC metrics
quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
#plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
quantile(pbmc$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.9, 1, 0.01))
#plots[[2]] + geom_hline(yintercept = 22000)
quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
#plots[[3]] + geom_hline(yintercept = 20)
quantile(pbmc$percent.HB, seq(0.9, 1, 0.01))
#plots[[4]] + geom_hline(yintercept = 1)

#Set QC standards
#Genes
minGene=300
maxGene=10000
#counts
minUMI=600
#Mitochondria
pctMT=10
#Blood cells
pctHB=1

#Data QC and draw violin plot
pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
violin

#Save
ggsave("2vlnplot_after_qc.pdf", plot = violin, width = 14, height = 8) 
dim(pbmc)

#Save
saveRDS(pbmc,"1pbmc_qc.rds")

#Read
pbmc = readRDS("1pbmc_qc.rds")