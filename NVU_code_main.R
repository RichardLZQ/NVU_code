# Single cell data analysis pipeline
library(Seurat)
library(ggsci)
library(gprofiler2)
library(DOSE)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)


#______________________________________________________________________________________________________________________________#
# 1. Data import ####

old1.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_1/old1/outs/filtered_feature_bc_matrix/")
young1.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_1/young1/outs/filtered_feature_bc_matrix/")
oldex1.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_1/oldEx1/outs/filtered_feature_bc_matrix/")

old2.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_2/old2/outs/filtered_feature_bc_matrix/")
young2.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_2/young2/outs/filtered_feature_bc_matrix/")
oldex2.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_2/oldEx2/outs/filtered_feature_bc_matrix/")

old4.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_4/old4/outs/filtered_feature_bc_matrix/")
young4.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_4/young4/outs/filtered_feature_bc_matrix/")
oldex4.data <- Read10X(data.dir ="/media/khlab/Treasure/Seq/Ex4/Batch_4/oldEx4/outs/filtered_feature_bc_matrix/")



old1 <- CreateSeuratObject(counts  = old1.data, project = "old1",min.cells = 3)
old2 <- CreateSeuratObject(counts  = old2.data, project = "old2",min.cells = 3)
old4 <- CreateSeuratObject(counts  = old4.data, project = "old4",min.cells = 3)

young1 <- CreateSeuratObject(counts  = young1.data, project = "young1",min.cells = 3)
young2 <- CreateSeuratObject(counts  = young2.data, project = "young2",min.cells = 3)
young4 <- CreateSeuratObject(counts  = young4.data, project = "young4",min.cells = 3)

oldex1 <- CreateSeuratObject(counts  = oldex1.data, project = "oldex1",min.cells = 3)
oldex2 <- CreateSeuratObject(counts  = oldex2.data, project = "oldex2",min.cells = 3)
oldex4 <- CreateSeuratObject(counts  = oldex4.data, project = "oldex4",min.cells = 3)

old1 <- RenameCells(old1, add.cell.id = "B1") # Adding the batch label to cell names incase duplicated barcode
old2 <- RenameCells(old2, add.cell.id = "B2")
old4 <- RenameCells(old4, add.cell.id = "B4")

young1 <- RenameCells(young1, add.cell.id = "B1")
young2 <- RenameCells(young2, add.cell.id = "B2")
young4 <- RenameCells(young4, add.cell.id = "B4")

oldex1 <- RenameCells(oldex1, add.cell.id = "B1")
oldex2 <- RenameCells(oldex2, add.cell.id = "B2")
oldex4 <- RenameCells(oldex4, add.cell.id = "B4")

#______________________________________________________________________________________________________________________________#
## 2. Prepare use CCA based method to remove batch effect ####

Batch1 <- merge(x = old1, y=list(young1,oldex1))
Batch2 <- merge(x =old2, y=list(young2,oldex2))
Batch4 <- merge(x=old4, y=list(young4,oldex4))

ovy.list <- list(batch1= Batch1, batch2=Batch2, batch4=Batch4)

for (i in 1:length(x = ovy.list)) {
  ovy.list[[i]] <- NormalizeData(object = ovy.list[[i]], verbose = FALSE)
  ovy.list[[i]] <- FindVariableFeatures(object = ovy.list[[i]], 
                                        selection.method = "vst", nfeatures = 4000, verbose = FALSE)
}

ovy.anchors <- FindIntegrationAnchors(object.list = ovy.list, dims = 1:30,anchor.features = 4000)
ovy.integrated <- IntegrateData(anchorset = ovy.anchors, dims = 1:30)

ovy.integrated <- ScaleData(object = ovy.integrated,  vars.to.regress = c("nCount_RNA")) # Adding UMI and percentage of mito genes as regression variables

ovy.integrated <- RunPCA(object = ovy.integrated, npcs = 60, verbose = FALSE)
ovy.integrated <- RunTSNE(object = ovy.integrated,dims = 1:50)


#ovy.integrated <- JackStraw(object = ovy.integrated, num.replicate = 100,dims = 50)
#ovy.integrated <- ScoreJackStraw(object = ovy.integrated, dims = 1:50)

#JackStrawPlot(object = ovy.integrated, dims = 1:50) # Check the results

ovy.integrated <- FindNeighbors(object = ovy.integrated, dims = 1:50)

ovy.integrated <- FindClusters(object = ovy.integrated, resolution = 1.6)

#______________________________________________________________________________________________________________________________#
## 3. Adding sample and batch info to data ####

Idents(ovy.integrated) <- "orig.ident"

Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("old1","young1","oldex1")))<-"1" # Adding batch info to object
Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("old2","young2","oldex2")))<-"2"
Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("old4","young4","oldex4")))<-"4"
ovy.integrated$Batch <- Idents(ovy.integrated)

Idents(ovy.integrated) <- "orig.ident"

Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("old1","old2","old4")))<-"old" # Adding group info to object
Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("young1","young2","young4")))<-"young"
Idents(ovy.integrated, cells = WhichCells(object = ovy.integrated, idents = c("oldex1","oldex2","oldex4")))<-"oldex"
ovy.integrated$Age <- Idents(ovy.integrated)

#______________________________________________________________________________________________________________________________#
## 4.QC ####

DefaultAssay(object = ovy.integrated) <- "RNA"

mito.features <- grep(pattern = "^mt-", x = rownames(x = ovy.integrated), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ovy.integrated, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ovy.integrated, slot = 'counts'))

ribo.feature <- grep(pattern = c("^Rpl|^Rps"), x = rownames(x = ovy.integrated), value = TRUE)
percent.ribo <- Matrix::colSums(x = GetAssayData(object = ovy.integrated, slot = 'counts')[ribo.feature, ]) / Matrix::colSums(x = GetAssayData(object = ovy.integrated, slot = 'counts'))

ovy.integrated[['percent.mito']] <- percent.mito # Adding to meta data
ovy.integrated[['percent.ribo']] <- percent.ribo

VlnPlot(object = ovy.integrated, features = c("nCount_RNA","nFeature_RNA"),pt.size = 0,ncol = 1)
#FeatureScatter(object = ovy.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito") # Generate some QC plots
#FeatureScatter(object = ovy.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#FeatureScatter(object = ovy.integrated, feature1 = "percent.mito", feature2 = "percent.ribo")

DefaultAssay(object = ovy.integrated) <- "integrated"
ovy.clean <- subset(x = ovy.integrated, subset =quantile(ovy.integrated$nFeature_RNA,0.95,na.rm = T) > nFeature_RNA & nFeature_RNA > quantile(ovy.integrated$nFeature_RNA,0.05,na.rm = T) & quantile(ovy.integrated$nCount_RNA,0.95,na.rm = T)  > nCount_RNA & nCount_RNA > quantile(ovy.integrated$nCount_RNA,0.05,na.rm = T) & percent.mito < 0.20 )

#______________________________________________________________________________________________________________________________#
## 5.Cell type ####

Idents(ovy.clean) <- "integrated_snn_res.1.6"
nFeature <- ovy.clean$nFeature_RNA %>% tibble::enframe(name = "BC",value = "nFeature")
clust <- ovy.clean$integrated_snn_res.1.6 %>% tibble::enframe(name = "BC",value = "cluster")
nFeature <- inner_join(nFeature, clust, by="BC")
low_quality <- nFeature %>% group_by(cluster) %>% summarise(nF_clust=mean(nFeature)) %>% filter(nF_clust<500) #20,21,36 was removed here
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("19")))<-"PC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("12","35")))<-"SMC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("0","2","8","22")))<-"EC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("1","3","9","28","29")))<-"MG"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("4","5","6","25")))<-"AC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("31","11")))<-"OPC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("32")))<-"NRP"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("39","52")))<-"OLG"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("24")))<-"imNeur"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("34")))<-"mNeur"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("23")))<-"EPC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("10","14","18")))<-"CPC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("30","43")))<-"Hb_EC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("17")))<-"MAC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("38")))<-"TNC"
Idents(ovy.clean, cells = WhichCells(object = ovy.clean, idents = c("33")))<-"MNC"
ovy.clean$Celltype <- Idents(ovy.clean)
ovy.clean <- subset(ovy.clean, idents =c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","NRP","EPC","CPC","Hb_EC","MAC","TNC","MNC"))
DimPlot(ovy.clean, group.by = "Celltype",label = T) #Check cell types


#______________________________________________________________________________________________________________________________#
## 6.Subtype ####

DefaultAssay(object = ovy.clean) <- "RNA"
Idents(ovy.clean) <- "Celltype"

EC <- subset(ovy.clean,idents = "EC")
SMC <- subset(ovy.clean, idents = "SMC")

EC <- as.SingleCellExperiment(EC)
SMC <- as.SingleCellExperiment(SMC)
saveRDS(EC,"./inter_data/EC.data")
saveRDS(SMC,"./inter_data/SMC.data")

EC.type <- readRDS("./cellassign/EC_results.rds")
SMC.type <- readRDS("./cellassign/SMC_results.rds")
EC.tmp<- as.character(EC.type)
names(EC.tmp)<- names(EC.type)

SMC.tmp <- as.character(SMC.type)
names(SMC.tmp) <- names(SMC.type)
ovy.clean$Subtype <- c(EC.tmp, SMC.tmp) # Send the results back to Seurat object

Idents(ovy.clean) <- "Celltype"

Celltype <- ovy.clean$Celltype
#______________________________________________________________________________________________________________________________#
## 7.Redo the normalization and scale ####
Batch1 <- merge(x = old1, y=list(young1,oldex1))
Batch2 <- merge(x =old2, y=list(young2,oldex2))
Batch4 <- merge(x=old4, y=list(young4,oldex4))

Batch1 <- subset(Batch1, cells =  intersect(colnames(Batch1), substr(names(Celltype),1,19))) # Only keep cells can be identified
Batch2 <- subset(Batch2, cells =  intersect(colnames(Batch2), substr(names(Celltype),1,19)))
Batch4 <- subset(Batch4, cells =  intersect(colnames(Batch4), substr(names(Celltype),1,19)))

ovy.list <- list(batch2=Batch2, batch4=Batch4)

ovy.final <- merge(x=Batch1, y=ovy.list)

ovy.final <- NormalizeData(ovy.final, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

ovy.final <- FindVariableFeatures(object = ovy.final, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

DefaultAssay(object = ovy.final) <- "RNA"

mito.features <- grep(pattern = "^mt-", x = rownames(x = ovy.final), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ovy.final, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ovy.final, slot = 'counts'))

ribo.feature <- grep(pattern = c("^Rpl|^Rps"), x = rownames(x = ovy.final), value = TRUE)
percent.ribo <- Matrix::colSums(x = GetAssayData(object = ovy.final, slot = 'counts')[ribo.feature, ]) / Matrix::colSums(x = GetAssayData(object = ovy.final, slot = 'counts'))

ovy.final[['percent.mito']] <- percent.mito # Adding to meta data
ovy.final[['percent.ribo']] <- percent.ribo

ovy.final <- ScaleData(object = ovy.final,  vars.to.regress = c("nCount_RNA", "percent.mito")) # Adding UMI and percentage of mito genes as

##7.1 Remarker the cells 

Idents(ovy.final) <- "orig.ident"

Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("old1","young1","oldex1")))<-"1" # Adding batch info to object
Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("old2","young2","oldex2")))<-"2"
Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("old4","young4","oldex4")))<-"4"
ovy.final$Batch <- Idents(ovy.final)

Idents(ovy.final) <- "orig.ident"

Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("old1","old2","old4")))<-"old" # Adding group info to object
Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("young1","young2","young4")))<-"young"
Idents(ovy.final, cells = WhichCells(object = ovy.final, idents = c("oldex1","oldex2","oldex4")))<-"oldex"
ovy.final$Age <- Idents(ovy.final)

ovy.final$Celltype <- ovy.clean$Celltype
ovy.final$Celltype<- factor(ovy.final$Celltype, levels = c("MNC","TNC","MAC","Hb_EC","CPC","EPC","NRP","mNeur","imNeur","OLG", "AC","MG", "EC","SMC","PC","OPC"))

ovy.final$Subtype <- ovy.clean$Subtype
ovy.final$Subtype<- factor(ovy.final$Subtype, levels = c( "A1","A2","aaSMC","aSMC","AV","Cap", "V", "VCap", "vSMC","EC_unclassified","SMC_unclassified"))

#7.2 Redo tSNE and PCA

ovy.final <-  RunPCA(object = ovy.final, npcs = 50, verbose = T)
ovy.final <- RunTSNE(object = ovy.final,dims = 1:30)
ovy.final <- FindNeighbors(object = ovy.final, dims = 1:30)
# ovy.final <- FindClusters(object = ovy.final, resolution = 0.6)
# ovy.final <- FindClusters(object = ovy.final, resolution = 0.4)
ovy.final <- FindClusters(object = ovy.final, resolution = 1.6)

#______________________________________________________________________________________________________________________________#
# 8.Find degs ####

### 8.1 EC and subtypes

Idents(ovy.final) <- "Subtype"

A1.data <- subset(ovy.final, idents = "A1") # Extract data into individual obj ovy.finalts
A2.data <- subset(ovy.final, idents = "A2")
V.data <- subset(ovy.final, idents = "V")
AV.data <- subset(ovy.final, idents = "AV")
Cap.data <- subset(ovy.final, idents = "Cap")
VCap.data <- subset(ovy.final, idents = "VCap")
Idents(ovy.final) <- "Celltype"
EC.data <- subset(ovy.final,idents = "EC")

Idents(A1.data) <- "Age" # Set the idents back to age
Idents(A2.data) <- "Age"
Idents(V.data) <- "Age"
Idents(AV.data) <- "Age"
Idents(Cap.data) <- "Age"
Idents(VCap.data) <- "Age"
Idents(EC.data) <- "Age"

EC.ove.diff <- FindMarkers(object = EC.data, ident.1 = "old", ident.2 = "oldex",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
A1.ove.diff <- FindMarkers(object = A1.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
A2.ove.diff <- FindMarkers(object = A2.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
V.ove.diff <- FindMarkers(object = V.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
AV.ove.diff <- FindMarkers(object = AV.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
Cap.ove.diff <- FindMarkers(object = Cap.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
VCap.ove.diff <- FindMarkers(object = VCap.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")
EC.ove.dif <- list(EC=EC.ove.diff, A1=A1.ove.diff,A2=A2.ove.diff,Cap=Cap.ove.diff,VCap=VCap.ove.diff,V=V.ove.diff, AV=AV.ove.diff)

EC.ovy.diff <- FindMarkers(object = EC.data, ident.1 = "old", ident.2 = "young",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
A1.ovy.diff <- FindMarkers(object = A1.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
A2.ovy.diff <- FindMarkers(object = A2.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
V.ovy.diff <- FindMarkers(object = V.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
AV.ovy.diff <- FindMarkers(object = AV.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
Cap.ovy.diff <- FindMarkers(object = Cap.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
VCap.ovy.diff <- FindMarkers(object = VCap.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
EC.ovy.dif <- list(EC=EC.ovy.diff, A1=A1.ovy.diff,A2=A2.ovy.diff,Cap=Cap.ovy.diff, VCap=VCap.ovy.diff,V=V.ovy.diff, AV=AV.ovy.diff)


### 8.2 SMCs
Idents(ovy.final) <- "Celltype"
SMC.data <-subset(ovy.final,idents = "SMC")
Idents(ovy.final) <- "Subtype"
aaSMC.data <- subset(ovy.final,idents = "aaSMC")
aSMC.data <- subset(ovy.final,idents = "aSMC")
vSMC.data <- subset(ovy.final,idents = "vSMC")

Idents(SMC.data) <-"Age"
Idents(aaSMC.data) <-"Age"
Idents(aSMC.data) <-"Age"
Idents(vSMC.data) <-"Age"

SMC.ovy.diff <- FindMarkers(object = SMC.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
aSMC.ovy.diff <- FindMarkers(object = aSMC.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
aaSMC.ovy.diff <- FindMarkers(object = aaSMC.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")
vSMC.ovy.diff <- FindMarkers(object = vSMC.data, ident.1 = "old", ident.2 = "young", logfc.threshold =0, assay = "RNA",test.use = "MAST")

SMC.ovy.dif <- list(SMC=SMC.ovy.diff, aSMC= aSMC.ovy.diff, aaSMC= aaSMC.ovy.diff, vSMC= vSMC.ovy.diff)
SMC.ove.diff <- FindMarkers(object = SMC.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold =0, assay = "RNA",test.use = "MAST")


### 8.3 PCs

Idents(ovy.final) <- "Celltype"
PC.data <- subset(ovy.final, idents = "PC")
Idents(PC.data) <- "Age"
PC.ovy.diff <- FindMarkers(object = PC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

PC.ove.diff <- FindMarkers(object = PC.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

### 8.4 ACs

Idents(ovy.final) <- "Celltype"
AC.data <- subset(ovy.final, idents = "AC")
Idents(AC.data) <- "Age"
AC.ovy.diff <- FindMarkers(object = AC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

AC.ove.diff <- FindMarkers(object = AC.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

### 8.5 MGs

Idents(ovy.final) <- "Celltype"
MG.data <- subset(ovy.final, idents = "MG")
Idents(MG.data) <- "Age"
MG.ovy.diff <- FindMarkers(object = MG.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
MG.ove.diff <- FindMarkers(object = MG.data, ident.1 = "old", ident.2 = "oldex", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

### 8.6 All other celltypes
Idents(ovy.final) <- "Celltype"
MNC.data <- subset(ovy.final, idents = "MNC")
TNC.data <- subset(ovy.final, idents = "TNC")
MAC.data <- subset(ovy.final, idents = "MAC")
Hb_EC.data <- subset(ovy.final, idents = "Hb_EC")
CPC.data <- subset(ovy.final, idents = "CPC")
EPC.data <- subset(ovy.final, idents = "EPC")
NRP.data <- subset(ovy.final, idents = "NRP")
mNeur.data <- subset(ovy.final, idents = "mNeur")
imNeur.data <- subset(ovy.final, idents = "imNeur")

Idents(MNC.data) <- "Age"
Idents(TNC.data) <- "Age"
Idents(MAC.data ) <- "Age"
Idents(Hb_EC.data) <- "Age"
Idents(CPC.data) <- "Age"
Idents(EPC.data) <- "Age"
Idents(NRP.data) <- "Age"
Idents(mNeur.data) <- "Age"
Idents(imNeur.data) <- "Age"

TNC.ovy.diff <- FindMarkers(object = TNC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
MNC.ovy.diff <- FindMarkers(object = MNC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
MAC.ovy.diff <- FindMarkers(object = MAC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
Hb_EC.ovy.diff <- FindMarkers(object = Hb_EC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
CPC.ovy.diff <- FindMarkers(object = CPC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
EPC.ovy.diff <- FindMarkers(object = EPC.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
NRP.ovy.diff <- FindMarkers(object = NRP.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
mNeur.ovy.diff <- FindMarkers(object = mNeur.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")
imNeur.ovy.diff <- FindMarkers(object = imNeur.data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, assay = "RNA",test.use = "MAST")

othertype.deg.list <- list(SMC=SMC.ovy.diff,PC=PC.ovy.diff,AC=AC.ovy.diff,MG=MG.ovy.diff,
                           MNC=MNC.ovy.diff,TNC=TNC.ovy.diff,MAC=MAC.ovy.diff,Hb_EC=Hb_EC.ovy.diff,
                           CPC=CPC.ovy.diff,EPC=EPC.ovy.diff,NRP=NRP.ovy.diff,mNeur=mNeur.ovy.diff,
                           imNeur=imNeur.ovy.diff)
#
othertype.deg.ovy <- list(SMC=SMC.ovy.diff,PC=PC.ovy.diff,AC=AC.ovy.diff,MG=MG.ovy.diff)
othertype.deg.ove <- list(SMC=SMC.ove.diff,PC=PC.ove.diff,AC=AC.ove.diff,MG=MG.ove.diff)
#Save the degs to folders

degsv <- function(deglist,fc){
  for(i in 1:length(deglist)){
    tmp <- deglist[[i]]
    tmp.up <- tmp[(tmp$p_val_adj<0.05)&((tmp$avg_logFC)>fc),]
    tmp.down <- tmp[(tmp$p_val_adj<0.05)&((tmp$avg_logFC)<(-fc)),]
    write.csv(tmp.up,file = paste0("./DEGs/",fc,"/",names(deglist)[i],"_",fc,"_up.csv"))
    write.csv(tmp.down,file = paste0("./DEGs/",fc,"/",names(deglist)[i],"_",fc,"_down.csv"))
  }
}

degsv_all <- function(deglist,fc){
  for(i in 1:length(deglist)){
    tmp <- deglist[[i]]
    tmp.up <- tmp[(tmp$p_val_adj<0.05)&(abs(tmp$avg_logFC)>fc),]
    write.csv(tmp.up,file = paste0("./DEGs/SMCPC/",names(deglist)[i],".csv"))
  }
}

#Compare the DEGs between batch_5 and batch_1-4
library(devtools)
source_gist("524eade46135f6348140")
degbat <- function(x,y){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(x)) {
    df1 <- x[[i]]
    df2 <- y[[i]]
    its.genes <- intersect(rownames(df1),rownames(df2))
    sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
    sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
    sig=sig1+sig2
    df <- tibble(batch_4=df1[its.genes,"avg_logFC"],batch_5=df2[its.genes,"avg_logFC"], 
                 sig=sig ,row.names = its.genes)
    df <- df %>% filter(sig>1) # Filter the unsignificant genes
    new.levels <- c("batch_4_sig","batch_5_sig","both_sig")
    df$sig <- factor(new.levels[df$sig],levels = c("batch_4_sig","batch_5_sig","both_sig"))
    #write.csv(df,file=paste0("./bt5/",names(x)[i],".csv"))
    df <- df %>%  filter(sig=="both_sig")
    write.csv(df,file=paste0("./bt5/",names(x)[i],".csv"))
    formula <- y ~ x
    print(ggplot(df, aes(x=batch_4,y=batch_5))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
            stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T)+
            geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
            ggtitle(paste0(names(x)[i],"_Subtype"))+theme_bw()+
            stat_fit_glance(method = "lm", 
                            method.args = list(formula = formula),
                            label.x = "right",
                            label.y = "bottom",
                            aes(label = paste("italic(P)*\"-value = \"*", 
                                              signif(..p.value.., digits = 4), sep = "")),
                            parse = TRUE)+
            geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 
pdf("Batch_coor_sig.pdf")
degbat(EC.ovy.bt4,EC.ovy.bt5)
dev.off()

