if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
           'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
           'SummarizedExperiment', 'batchelor', 'HDF5Array',
           'terra', 'ggrastr'))
BiocManager::install("monocle")
BiocManager::install("HSMMSingleCell")


#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(Rogue)
library(clustree)  #?
library(harmony)
library(SingleR)
library(dplyr)
library(monocle)
library(tidyverse)



#Directory
setwd("E:/learning_materials/sc/11monocle/11monocle")

#Read in
sco <- readRDS("sco.brca.cd4Tcell.rds")

#Plot
DimPlot(sco, group.by = "celltype_rename", label = T) + DimPlot(sco, group.by = "seurat_clusters", label = T)

#Rename
sco$celltype <- sco$celltype_rename
colnames(sco@meta.data)

#Get expression matrix
data <- GetAssayData(sco, assay = "RNA", slot = "counts")
#Get cell annotation information, phenotype information
pd <- new('AnnotatedDataFrame', data = sco@meta.data[,c(1,5,15,17)])  # monocle requires metadata to be in AnnotatedDataFrame
#Get gene names
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Create a new cell dataset object
dim(sco)



# create monocle object: cellDataSet
mycds <- newCellDataSet(as.matrix(data),
            #If data is too large, convert to sparse matrix
            #as(as.matrix(data),"sparseMatrix"),
            phenoData = pd,
            featureData = fd,
            lowerDetectionLimit = 0.5,
            expressionFamily = negbinomial.size() )# data is original UMI counts and required negative binomial modal to handle




#Data preprocessing, estimate size factors and dispersion, similar to normalization, standardization
#Takes a long time to run, can modify the number of cores
mycds <- estimateSizeFactors(mycds)      # step1    
mycds <- estimateDispersions(mycds, cores=8)  # step2

### Select ordering genes, around 2000
# select high variant genes
disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 1 * dispersion_fit) %>% pull(gene_id) %>% as.character()

#Mark these genes
mycds <- setOrderingFilter(mycds, order.genes) #only using order.genes in the following steps

#Set ordering genes
plot_ordering_genes(mycds)
p <- plot_ordering_genes(mycds)
ggsave("0OrderGenes_1.pdf", p, width = 8, height = 6)



### Dimensionality reduction ordering
#residualModelFormulaStr 
#Subtract the influence of "uninteresting" sources of variation to reduce their impact on clustering.
mycds <- reduceDimension(mycds, max_components = 2, reduction_method = 'DDRTree', 
            residualModelFormulaStr = "~orig.ident")
# residualModelFormulaStr = "~orig.ident", remove batch effect. build a model: gene expression ~ orig.ident
# algorithm DDRTree, dimensionality reduction. build the tree which is trajectory



#Takes a long time to run
mycds <- orderCells(mycds)
#If an error occurs
#Error :
# !nei() was deprecated in igraph 2.1.0 and is now defunct.
# please use.nei() instead.
#Run rlang::last trace() to see where the error occurred.

#Need to downgrade igraph to version 1.5.1
#First uninstall this package 
# remove.packages("igraph")
#Then install the specified version 
# packageurl = 'https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.5.1.tar.gz'
# install.packages(packageurl, repos = NULL, type = 'source')
#If it still fails to install, watch 2.R package installation video, install locally.

#Error
#Error in if (class(projection) != "matrix") projection <- as.matrix(projection) :
#Run
#trace('project2MST', edit = T, where = asNamespace("monocle"))
#Find
#if (elass(projection) != "matrix")
#Delete and save

# Can manually set the starting point, root
# mycds <- orderCells(mycds,root_state = 5) #must run after: mycds <- orderCells(mycds)!!!

# Visualize results
# naive (Tn), central memory (Tcm), effector memory (Tem)
# regulatory (tregs)，T helper (Th)




# State
p1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("1Trajectory_State_1.pdf", plot = p1, width = 10, height = 6.5)
# Pseudotime
p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("2Trajectory_Pseudotime_1.pdf", plot = p2, width = 10, height = 6.5)
# Celltype
p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("3Trajectory_Celltype2_1.pdf", plot = p3, width = 10, height = 6.5)
# orig.ident
p4 <- plot_cell_trajectory(mycds, color_by = "orig.ident")
ggsave("4Trajectory_Sample_1.pdf", plot = p4, width = 10, height = 6.5)
# Dendrogram
p5 <- plot_complex_cell_trajectory(mycds, x = 1, y = 2,
                 color_by = "celltype")
ggsave("5Trajectory_dendrogram_1.pdf", plot = p5, width = 10, height = 6.5)
# Cell density plot
p6 <- ggplot(pData(mycds),aes(Pseudotime,colour = celltype,fill = celltype)) +
 geom_density(bw = 0.5, size = 1, alpha = 0.5)+theme_classic()
ggsave("6Trajectory_Density_1.pdf", plot = p6, width = 10, height = 6.5)
# Expression change of specified genes
genes = c(order.genes)[1:4]
p1 = plot_genes_in_pseudotime(mycds[genes],color_by="State")
p2 = plot_genes_in_pseudotime(mycds[genes],color_by="celltype")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by="Pseudotime")
ggsave("7Trajectory_Pseudotime_1.pdf", plot = p1|p2|p3, width = 10, height = 6.5)

p1 = plot_genes_jitter(mycds[genes],grouping="State",color_by="State")
p2 = plot_genes_violin(mycds[genes],grouping="State",color_by="State")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by="State")
ggsave("8Trajectory_jitter_1.pdf", plot = p1|p2|p3, width = 10, height = 6.5)

pData(mycds)$SAT1 = log2(exprs(mycds)["SAT1",] +1)
p1 = plot_cell_trajectory(mycds,color_by = "SAT1") + 
 scale_color_continuous(type = "viridis")

pData(mycds)$KLRB1 = log2(exprs(mycds)["KLRB1",] +1)
p2 = plot_cell_trajectory(mycds,color_by = "KLRB1") + 
 scale_color_continuous(type = "viridis")

pData(mycds)$HSPA8 = log2(exprs(mycds)["HSPA8",] +1)
p3 = plot_cell_trajectory(mycds,color_by = "HSPA8") + 
 scale_color_continuous(type = "viridis")

pData(mycds)$G3BP2 = log2(exprs(mycds)["G3BP2",] +1)
p4 = plot_cell_trajectory(mycds,color_by = "G3BP2") + 
 scale_color_continuous(type = "viridis")

ggsave("9Trajectory_Expression_1.pdf",plot = p1|p2|p3|p4,width=14,height=6.5)

#Save to seurat object
pdata <- Biobase::pData(mycds)
sco <- AddMetaData(sco, metadata = pdata[,c("Pseudotime","State")])
saveRDS(sco, file = "10sco.pseudotime.rds")





#Find pseudotime differential genes, using monocle method
# does the gene expression change as the change of psedotime changes?
Time_diff <- differentialGeneTest(mycds, cores = 10,
                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff, "11Time_diff_all.csv", row.names = F)
#Plot, sort by qval, select top 100
Time_genes <- Time_diff[order(Time_diff$qval), "gene_short_name"][1:100]
#num_clusters cluster by row, how many clusters
p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3, 
              show_rownames=T, return_heatmap=T)
ggsave("12Time_heatmap_1.pdf", p, width = 5, height = 10)
#Save
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "13Time_diff_sig_1.csv", row.names = F)

# Single-cell trajectory "branch" analysis
# Find genes associated with branch points
# BEAM analysis, used to find genes regulated in a branch-dependent manner.
beam_res <- BEAM(mycds, branch_point = 1, cores = 10,   # Focus on branch point one
        progenitor_method = "duplicate")
write.csv(beam_res, "14BEAM_all_1.csv", row.names = F)
#Top 100
BEAM_genes <- beam_res[order(beam_res$qval), "gene_short_name"][1:100]
p <- plot_genes_branched_heatmap(
 mycds[BEAM_genes,], branch_point = 1, num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave("15BEAM_heatmap_1.pdf", p$ph_res, width = 6.5, height = 10)
#Save
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- beam_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig, "16BEAM_sig_1.csv", row.names = F)

