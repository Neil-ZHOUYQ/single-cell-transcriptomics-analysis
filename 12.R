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
#CellChat version 1.6.1

#Directory
setwd("E:/learning_materials/sc/12CellChat/12CellChat")

#Read in
pbmc <- readRDS("sco.brca.cd4Tcell.rds")

#Create CellChat object
table(pbmc@meta.data$celltype_rename)
table(pbmc@meta.data$orig.ident)

#Note, the first input is the expression matrix, you can try several formats
#pbmc@assays$RNA@data
#pbmc@assays$SCT@data
#GetAssayData(object = pbmc@assays$RNA, layer = "data")
#GetAssayData(object = pbmc@assays$SCT, layer = "data")

cellchat <- createCellChat(pbmc@assays$RNA@data,meta=pbmc@meta.data, # build CellChat object, analyzing with group information
             group.by="celltype_rename")

#Set ligand-receptor interaction database, use CellChatDB.mouse for mouse
CellChatDB <- CellChatDB.human   # human ligand-receptor database

#Autocrine/Paracrine signaling interactions
#Extracellular matrix (ECM) receptor interactions
#Cell-cell contact interactions
showDatabaseCategory( CellChatDB)

#Add reference database
#Can choose a subset of the database for analysis
cellchat@DB
#CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB





#Extract cell communication signaling genes
cellchat <- subsetData(cellchat)

#Multithreading
future::plan("multisession", workers = 12)

#Identify over-expressed ligand, receptor genes
cellchat <- identifyOverExpressedGenes(cellchat)

#Identify over-expressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)

#Data correction (optional)
#Smooth gene expression values based on those defined in a highly credible experimentally validated PPI network
cellchat <- projectData(cellchat, PPI.human)

#Calculate communication probability and infer cellchat network, if using non-corrected data, raw.use = T
cellchat <- computeCommunProb(cellchat, raw.use = F)  # CORE COMPUTING STEP 1: compute the communication probability of each ligand-receptor pairs
cellchat <- filterCommunication(cellchat, min.cells = 10) # filter the little expressed communication pairs

#Infer cell-cell communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)     # CORE COMPUTING STEP 2: enrich the ligand-receptor pairs to their signaling pathways 

#Extract cell interaction results
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "1Gene.csv", row.names = F) # gene level

#Signaling pathway level
df.netP <- subsetCommunication(cellchat, slot.name = "netP") # pathway level 
write.csv(df.netP, "2Pathway.csv", row.names = F)

#Calculate integrated cell communication network
cellchat <- aggregateNet(cellchat)

#Save results
saveRDS(cellchat, file = "3cellchat.rds")

cellchat = readRDS("3cellchat.rds")
#Plot
#Total
groupSize <- as.numeric(table(cellchat@idents))
pdf("4NetVisual_overview_all.pdf", width = 8, height = 6)
par(xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
        label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
        label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#Separate
pdf("5NetVisual_overview_split.pdf", width = 6, height = 5)
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
 mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
 mat2[i, ] <- mat[i, ]
 par(xpd = TRUE)
 netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
         edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#By signaling pathway, take the first 5 as an example
mypathways <- cellchat@netP$pathways
mypathways <- mypathways[1:5]
mypathways

#Signaling pathway
pdf("6NetVisual_pathways_circle.pdf", width = 6, height = 5)
for(pathways.show in mypathways){
  par(xpd = TRUE)
 netVisual_aggregate(cellchat,signaling = pathways.show,layout = "circle")
}
dev.off()

#Chord diagram
pdf("7NetVisual_pathways_chord.pdf", width = 10, height = 8)
for(pathways.show in mypathways){
 netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
}
dev.off()

#Heatmap
pdf("8NetVisual_pathways_heatmap.pdf", width = 10, height = 7)
for(pathways.show in mypathways){
  par(xpd = TRUE)
 p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
 plot(p)
}
dev.off()

#Ligand-receptor within signaling pathway
dir.create("9Pathways")
for(pathways.show in mypathways){
 pdf(paste0("9Pathways/", pathways.show, ".pdf"), width = 8, height = 6.5)
 # Display ligand-receptor contribution
 netAnalysis_contribution(cellchat, signaling = pathways.show)
 pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
 for(LR.show in pairLR){
  # Network plot shows cell-cell ligand-receptor interaction
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  }
 for(LR.show in pairLR){
  # Chord diagram shows cell-cell ligand-receptor interaction
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
 }
 dev.off()
}

levels(cellchat@idents)

#Show all cell-cell ligand-receptor interactions
p <- netVisual_bubble(cellchat, sources.use = 1:length(levels(cellchat@idents)), 
           targets.use = 1:length(levels(cellchat@idents)), remove.isolate = FALSE)
ggsave("10CCI_all.pdf", p, width = 8, height = 20, limitsize = F)

#Specify cell-cell ligand-receptor interactions
p <- netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4:5, remove.isolate = FALSE)
ggsave("11CCI_subcell.pdf", p, width = 8, height = 15, limitsize = F)

#Specify cell and signaling pathway ligand-receptor interactions
p <- netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4:5, signaling = c("CCL","TNF"), remove.isolate = FALSE)
ggsave("12CCI_subcell_subpathway.pdf", p, width = 8, height = 6, limitsize = F)

#Specify cell and signaling pathway ligand-receptor interactions
pairLR.use <- c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1B")
pairLR.use <- data.frame(interaction_name = pairLR.use)
p <- netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4:5, pairLR.use = pairLR.use, remove.isolate = FALSE)
ggsave("13CCI_subcell_subLR.pdf", p, width = 8, height = 6, limitsize = F)

#Ligand, receptor gene expression
p <- plotGeneExpression(cellchat, signaling = "CCL")
ggsave("14GeneExpression_violin_sig.pdf", p, width = 10, height = 9, limitsize = F)
p <- plotGeneExpression(cellchat, signaling = "CCL", enriched.only = FALSE)
ggsave("15GeneExpression_violin_all.pdf", p, width = 10, height = 9, limitsize = F)

#Cell communication network system analysis
#Calculate multiple network central measures for each cell group
#Identify dominant senders, receivers, mediators, and influencers in cell-cell communication networks at any time.
cellchat2 <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("16SignalingRole.pdf", width = 6, height = 4.5)
for(pathways.show in mypathways){
 netAnalysis_signalingRole_network(cellchat2,signaling=pathways.show,
                  width=8,height=2.5,font.size=10)
}
dev.off()