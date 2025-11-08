#Load
library(Seurat)
library(data.table)
library(stringr)
library(tibble)


#Working directory
setwd("E:/learning_materials/sc/3Read/3Read")


#####1.matrix.mtx, genes.tsv and barcodes.tsv####

#Directory
setwd("E:/learning_materials/sc/3Read/3Read/1GSE234527")

# Get a list of all sample files in the data folder
samples <- list.files("seurat/")
samples

# Create an empty list
seurat_list <- list()

#Read data and create Seurat object
#Delete cells with less than 200 genes expressed, and genes expressed in less than 3 cells
for (sample in samples) {
  #File path
  data.path <- paste0("seurat/", sample)
  
  #Read 10x data
  seurat_data <- Read10X(data.dir = data.path)

  #Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  #Add to list
  seurat_list <- append(seurat_list, seurat_obj)
}

#Merge
seurat_combined <- merge(x = seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

old_seurat_combined <- seurat_combined

meta <- seurat_combined@meta.data

batch1 <- seurat_combined@assays$RNA@layers$counts.GSM7470392

#Join Layers
pbmc = JoinLayers(seurat_combined)

combined <- pbmc@assays$RNA@layers$counts
020308combined1 <- pbmc[["RNA"]][["counts"]]
a = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "counts")[1:20,1:20])

#####2.H5 format#####

#Directory
setwd("F:/4scRNA/1code/3/2GSE199866")

#Get file names
fs=list.files(pattern = '.h5')
fs

#Read in sequentially and create Seurat objects
sceList = lapply(fs, function(x){
  a=Read10X_h5(x)
  p=str_split(x,'_',simplify = T)[,1]
  sce <- CreateSeuratObject(a,project = p ,min.features = 200,min.cells = 3)
})

#Get sample GSM number
folders = substr(fs,1,10)
folders

#Use the merge function to combine
sce.big <- merge(sceList[[1]], 
                y = sceList[-1], 
                add.cell.ids = folders)

#Join Layers
pbmc = JoinLayers(sce.big)

#####3.R data files (RDS/RDATA files)####

#Directory
setwd("F:/4scRNA/1code/3/3")

#Read RDATA file
load(file ="1pbmc.RData")

#Read RDS file
pbmc2 = readRDS("1pbmc.rds")


#####4.TXT or CSV####

#Directory
setwd("E:/learning_materials/sc_bili/3Read/3Read/4/GSE153935")

#Read in, use fread, because the matrix is large and the format is messy
#Must use fread, otherwise other functions will take a long time
pbmc <-fread("GSE153935_TLDS_AllCells.txt", sep="\t")
pbmc[1:5,1:5]

#Directly convert to row names
pbmc <-  column_to_rownames(pbmc,"V1")
pbmc[1:5,1:5]

#Create seurat object
pbmc <- CreateSeuratObject(pbmc, min.features = 300, min.cells = 3)

#Join Layers
pbmc = JoinLayers(pbmc)





#Directory
setwd("E:/learning_materials/sc_bili/3Read/3Read/4/GSE165722")

#Get sample locations
fs2=list.files(pattern = '.txt')
fs2
fs3=list.files(pattern = '.tsv')
fs3

#Get sample names
folders = substr(fs2,1,10) #substr(string, start, stop)
folders



#Read in
sceList = list()
for (i in 1:length(fs2)) {
  #Read in
  abc123 = as.data.frame(fread(fs2[i]))
  abc456 = as.data.frame(fread(fs3[i]))
  #Convert the gene column to row names, gene names
  abc456 <-  column_to_rownames(abc456,"gene")
  #Add matrix column names, cell names
  colnames(abc456) = abc123[,1]
  #Create seurat object
  sce <- CreateSeuratObject(abc456,project = folders[i],min.features = 300, min.cells = 3)
  #Put into the list sequentially
  sceList[i] = sce
}


#Merge
sce.big <- merge(sceList[[1]], 
                 y = sceList[-1], 
                 add.cell.ids = folders)

#Join Layers
pbmc = JoinLayers(sce.big)



#Add cell annotation information

#Get cell annotation information
abc789 = pbmc@meta.data 

#Export
write.table(data.frame(ID=rownames(abc789),abc789),file="meta.txt", sep="\t", quote=F, row.names = F,col.names = T)

#Read in
meta = fread("meta.xlsx")

#Convert the first column to row names
meta <-  column_to_rownames(meta,"ID")

#Add meta.data information, the order of cell names must be consistent
pbmc <- AddMetaData(object = pbmc, 
                    metadata = meta,   
                    col.name = c("group1","group2")) 


#When inconsistent, sort according to the names in the seurat object
meta = meta[colnames(pbmc),]
