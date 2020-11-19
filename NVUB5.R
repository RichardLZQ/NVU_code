# Reanalysis the EC batch5 data for Glp1r project
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(devtools)
library(ggpmisc)
library(ggpubr)
source_gist("524eade46135f6348140") # Load the add formular to plot function
#______________________________________________________________________________________________________________________________#
# 1. Data import####

old5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/old5/")
young5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/young5/")
oldEx5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/oldEx5/")

old5 <- CreateSeuratObject(counts  = old5.data, project = "old5",min.cells = 5)
young5 <- CreateSeuratObject(counts  = young5.data, project = "young5",min.cells = 5)
oldEx5 <- CreateSeuratObject(counts = oldEx5.data, project = 'oldEx5', min.cells = 5)

old5 <- RenameCells(old5, add.cell.id = "B5O") 
young5 <- RenameCells(young5, add.cell.id = "B5Y")
oldEx5 <- RenameCells(oldEx5, add.cell.id = "B5E") 


ovy5 <- merge(x = old5,y = list(oldEx5,young5) ,project = "glp1r")

ovy5 <- PercentageFeatureSet(ovy5, pattern = "^mt-", col.name = "percent.mt")
ovy5 <- subset(x = ovy5, subset = quantile(ovy5$nCount_RNA,0.95,na.rm = T)  > nCount_RNA & nCount_RNA > quantile(ovy5$nCount_RNA,0.05,na.rm = T) & percent.mt < 20 )

ovy5 <- SCTransform(ovy5,verbose = T)

ovy5 <- RunPCA(ovy5)
ovy5 <- RunUMAP(ovy5,dims = 1:30,umap.method = "umap-learn")
ovy5 <- FindNeighbors(ovy5,dims = 1:30)
ovy5 <- FindClusters(ovy5,resolution = 1.5)
ovy5 <- FindClusters(ovy5,resolution = 1.4)
ovy5 <- FindClusters(ovy5,resolution = 1.2)
#______________________________________________________________________________________________________________________________#
# 3. Celltype identify ####
Idents(ovy5) <- "SCT_snn_res.1.4"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("21")))<-"PC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("15","27")))<-"SMC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("1","2","4","14","23")))<-"EC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("6","11","7","8","12","19","29")))<-"MG"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("3","0","5","38")))<-"AC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("9","31")))<-"OPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("40")))<-"OLG"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("24")))<-"imNeur"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("46","48")))<-"mNeur"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("33")))<-"NRP"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("36")))<-"EPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("16","20","30")))<-"CPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("32")))<-"Hb_EC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("22")))<-"MAC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("52")))<-"TNC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("45")))<-"MNC"
ovy5$Celltype <- Idents(ovy5)
ovy5 <- subset(ovy5,idents = c("PC","SMC","EC","MG","AC","OPC","NRP","OLG","EPC","TNC","imNeur","mNeur","CPC","Hb_EC","MAC","MNC"))

ovy5.less <- subset(ovy5, cells= sample(colnames(ovy5),6000)) #Downsample to 6000 cells
#______________________________________________________________________________________________________________________________#
# 4. DEG calculation ####
Idents(ovy5) <- "Celltype"
PC.data <- subset(ovy5,idents = c("PC"))
SMC.data <- subset(ovy5,idents = c("SMC"))
EC.data <- subset(ovy5,idents = c("EC"))
MG.data <- subset(ovy5,idents = c("MG"))
AC.data <- subset(ovy5,idents = c("AC"))
OPC.data <- subset(ovy5,idents = c("OPC"))
OLG.data <- subset(ovy5,idents = c("OLG"))
imNeur.data <- subset(ovy5,idents = c("imNeur"))
mNeur.data <- subset(ovy5,idents = c("mNeur"))
CPC.data <- subset(ovy5,idents = c("CPC"))
Hb.data <- subset(ovy5,idents = c("Hb_EC"))
MAC.data <- subset(ovy5,idents = c("MAC"))
MNC.data <- subset(ovy5,idents = c("MNC"))

celltype.list <- list(PC.data,SMC.data,EC.data,MG.data,AC.data,OPC.data,OLG.data,
                      imNeur.data,mNeur.data,CPC.data,Hb.data,MAC.data,MNC.data)

ovy.deg <- list() # old vs young
for (i in 1:length(celltype.list)) {
  Idents(celltype.list[[i]]) <- "orig.ident"
  ovy.deg[[i]] <- FindMarkers(celltype.list[[i]],ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(ovy.deg) <- c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","CPC","Hb_EC","MAC","MNC")


evo.deg <- list() # old with exenatide vs old
for (i in 1:length(celltype.list)) {
  Idents(celltype.list[[i]]) <- "orig.ident"
  evo.deg[[i]] <- FindMarkers(celltype.list[[i]],ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)
}

names(evo.deg) <- c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","CPC","Hb_EC","MAC","MNC")


# Export the deg list

deg.expt <- function(obj,path){
  for (i in 1:length(obj)) {
    write_csv(as_tibble(obj[[i]],rownames="Gene"), paste0(path,names(obj)[i],"_deg_all.csv"))
  }
  
}


deg.path <- "~/bioinfo/Richard/NVUB5/Text/DEG_list/OVY/"
deg.expt(ovy.deg,deg.path)

deg.path <- "~/bioinfo/Richard/NVUB5/Text/DEG_list/EVO/"
deg.expt(evo.deg,deg.path)