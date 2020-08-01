# Scriptes to generate the NVU paper figures
#Mypalette <- c("#C0344D","#66898D","#2D7A8E","#D0CC83","#D36F4C","#649FB1")
Mypalette <-  c("#E64B35FF","#F39B7FFF","#00A087FF","#227c91","#3C5488FF","#8491B4FF")
tSNE.palette <- c("#F8766D","#E68613","#CD9600","#ABA300","#7CAE00","#0CB702","#00BE67","#00C19A","#00BFC4","#00B8E7","#00A9FF","#8494FF","#C77CFF","#ED68ED","#FF61CC","#FF68A1")

library(Seurat)
library(tidyverse)
library(ggsci)
library(ggpmisc)
library(ggrepel)
library(venn)
library(reshape2)
library(readxl)


library(devtools)
source_gist("524eade46135f6348140") # Load the add formular to plot function

# Figure.1 ####

#F1.a Cell group tSNE####
pdf("./NVU_figure/F1/Figure_1_Age_tSNE.pdf")
Idents(ovy.clean) <- "Age"
DimPlot(ovy.clean,reduction = 'tsne',group.by = 'Age',cells = WhichCells(object = ovy.clean,idents = c('old','young')), label = F,cols = c("#cc6670","#3993ac"))
dev.off()

#F1.b Cell group tSNE####
pdf("./NVU_figure/F1/Figure_1_Main_cell_types.pdf")
Idents(ovy.clean) <- "Age"
ovy.oy <- subset(ovy.clean,cells = WhichCells(object = ovy.clean,idents = c('old','young')))
ovy.oy$Celltype <- factor(ovy.oy$Celltype,levels = c("Hb_EC","TNC","MAC","EPC","PC","EC","MNC","SMC","imNeur","CPC","AC","OLG","MG","mNeur","OPC","NRP"))
Idents(ovy.oy) <- "Celltype"
DimPlot(ovy.oy,label = T,group.by = "Celltype",coord.fixed = T,label.size = 5)
dev.off()

#F1.c EC subtypes tSNE####
Idents(ovy.oy) <- "Celltype"
ovy.oyec <- subset(ovy.oy,idents="EC")
ovy.oyless <- subset(ovy.oy, cells= sample(colnames(ovy.oyec),6000)) #Downsample to 6000 cells
ovy.oyless$Subtype <- factor(ovy.oyless$Subtype, levels = c("A1","A2","Cap","VCap","V","AV","Others"))
pdf("./NVU_figure/F1/Figure_1_Subtypes_tSNE.pdf",width = 11,height = 8)
for (i in 1:6) {
  Idents(ovy.oyless) <- "Subtype"
  tmp <- ovy.oyless$Subtype
  sub.name<- levels(ovy.oyless)[i]
  tmp[!grepl(paste0("^",levels(ovy.oyless)[i],"$"),ovy.oyless$Subtype)] <- "Others"
  ovy.oyless$singlesubtype <- tmp
  Idents(ovy.oyless) <- "singlesubtype"
  print(DimPlot(ovy.oyless,pt.size = 0.8,order=c(sub.name,"Others"),cols = c("#e6e6e6",Mypalette[i]),na.value = "#e6e6e6"))
}
dev.off()

#F1.d EC subtypes markers tSNE####
sub.marker <- c("Bmx", "Vegfc", "Mfsd2a", "Slc38a5", "Nr2f2", "Vwf")
pdf("./NVU_figure/F1/Figure_1_Subtype_markers.pdf",width = 11,height = 8)
for (i in 1:length(sub.marker)){
  print(FeaturePlot(ovy.oyless,features = sub.marker[i] ,cols =c("#e6e6e6","#0d0d5e"),pt.size = 0.8,order=T,coord.fixed = F)+expand_limits(colour = seq(0, 5, by = 1))
  )
}
dev.off()

#F1.e Subtypes cell proportion####
Idents(ovy.clean) <- "Age"
ovy.old <- subset(ovy.clean,cells = WhichCells(object = ovy.clean,idents = c('old')))
ovy.young <- subset(ovy.clean,cells = WhichCells(object = ovy.clean,idents = c('young')))
subtype.old <- as.data.frame(table(ovy.old$Subtype)) %>% as_tibble() %>% filter(Var1 %in% c("A1","A2","Cap","VCap","V","AV"))
subtype.old <- subtype.old[c(1,2,4,6,5,3),] %>% mutate(prop= round(Freq/sum(.$Freq)*100,3)) %>% mutate(lab.ypos = round(cumsum(prop) - 0.5*prop,2))
subtype.young <- as.data.frame(table(ovy.young$Subtype)) %>% as_tibble() %>% filter(Var1 %in% c("A1","A2","Cap","VCap","V","AV"))
subtype.young <- subtype.young[c(1,2,4,6,5,3),] %>% mutate(prop= round(Freq/sum(.$Freq)*100,3)) %>% mutate(lab.ypos = round(cumsum(prop) - 0.5*prop,2))
subtype.old$Var1 <- c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC")
subtype.young$Var1 <- c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC")
subtype.old$Var1 <- factor(subtype.old$Var1, levels = c("CapEC","vCapEC","vEC","avEC","aEC1","aEC2"))
subtype.young$Var1 <- factor(subtype.young$Var1, levels = c("CapEC","vCapEC","vEC","avEC","aEC1","aEC2"))
pdf("./NVU_figure/F1/Figure_1_Subtype_proportion.pdf")
ggplot(subtype.young, aes(x = 2, y = prop, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0,direction = -1)+
  #geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  scale_fill_manual(values = c("#00A087FF","#227c91","#3C5488FF","#8491B4FF","#E64B35FF","#F39B7FFF"))+
  theme_void()+
  ggtitle("Young")+
  xlim(0.5, 2.5)

ggplot(subtype.old, aes(x = 2, y = prop, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0,direction = -1)+
  #geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  scale_fill_manual(values = c("#00A087FF","#227c91","#3C5488FF","#8491B4FF","#E64B35FF","#F39B7FFF"))+
  theme_void()+
  ggtitle("Old")+
  xlim(0.5, 2.5)
dev.off()

###F1.f Subtype DEG venn diagram####
pdf("./NVU_figure/F1/Figure_1_DEG_Venn_up_down.pdf")
venn.list <- lapply(EC.ovy.dif[2:6], function(x){rownames(x)[(x$p_val_adj<0.05 & abs(x$avg_logFC)>0.1)]})
venn.list <- venn.list[c(3,4,5,1,2)]
venn(venn.list,zcolor = Mypalette[c(3,4,5,1,2)] , cexil = 1, cexsn = 0.8,borders = T,opacity = 0.4)

venn.list <- lapply(EC.ovy.dif[2:6], function(x){rownames(x)[(x$p_val_adj<0.05 & (x$avg_logFC)<(-0.1))]})
venn.list <- venn.list[c(3,4,5,1,2)]
venn(venn.list,zcolor = Mypalette[c(3,4,5,1,2)], cexil = 1, cexsn = 0.8,borders = T,opacity = 0.4)

dev.off()

###F1.g Subtype DEG quantification####
deglist <- EC.ovy.dif

deglist.up <- lapply(deglist[2:7], function(x){x <- x[(x$p_val_adj <0.05&x$avg_logFC>0.1),]}) 
deglist.down <- lapply(deglist[2:7], function(x){x <- x[(x$p_val_adj <0.05&x$avg_logFC<(-0.1)),]}) 

deglist.up <- sapply(deglist.up, nrow)
deglist.down <- sapply(deglist.down, nrow)
deglist.down <- deglist.down*-1

deglist.up <- as_tibble(deglist.up,rownames="Subtype")
deglist.down <- as_tibble(deglist.down, rownames="Subtype")
deglist.num <- inner_join(deglist.up, deglist.down, by="Subtype")

new.levels <- c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC")
old.levels <- c("A1","A2","Cap","VCap","V","AV")


deglist.num$Subtype <- factor(new.levels[1:6],levels = rev(new.levels))

pdf("./NVU_figure/F1/Figure_1_DEG_barplot.pdf",paper = "a4r")
ggplot(data = deglist.num)+geom_bar(aes(x=Subtype,y=value.x),position = "stack", stat="identity",fill="#E64B35",width =  0.75)+
  geom_bar(aes(x=Subtype,y=value.y),position = "stack", stat="identity",fill="steelblue",width = 0.75)+
  scale_y_continuous(breaks=seq(-100,250,50), limits=c(-100,250)) +coord_flip()
dev.off()

###F1.h Shared DEG heatmap ####
pdf("./NVU_figure/F1/Figure_1_DEG_heatmap_shared.pdf")
deglist <- EC.ovy.dif[2:7]
# UP
deg.sig <- lapply(deglist, function(x){x <- x[x$p_val_adj <0.05,]})
updeg.sig <- lapply(deg.sig, function(x){x <- x[x$avg_logFC >FCth,]})
upgene.names <- c() 
for (i in 1:length(deg.sig)) {
  upgene.names <- c(upgene.names, rownames(updeg.sig[[i]]))
} # Extract unique gene names
upgene.names <- upgene.names[!duplicated(upgene.names)] # Remove duplicated gene names
updeg.inx <- matrix(0,nrow = length(upgene.names), ncol = length(deg.sig),dimnames = (list(upgene.names, names(deg.sig)))) # Initialize the index matrix

for (i in 1:length(updeg.sig)){
  updeg.inx[rownames(updeg.sig[[i]]),i] <- 1
} # Mark the gene as significant in index matrix
updeg.ind <- updeg.inx
indsum <- apply(updeg.ind,1,sum)
updeg.ind <- cbind(updeg.ind, sum=indsum)
updeg.ind <- updeg.ind[order(updeg.ind[,ncol(updeg.ind)], decreasing = T),]
# 3. Do heatmap on filtered table and fill with logFC.
updeg.ind <- updeg.ind[updeg.ind[,ncol(updeg.ind)]>threshold,]
updeg.FC <- updeg.ind[,-ncol(updeg.ind)] # Initialize another matrix to store logFC value
for (i in 1:ncol(updeg.FC)){
  updeg.FC[updeg.ind[,i]==1,i] <- updeg.sig[[i]][rownames(updeg.ind)[updeg.ind[,i]==1],"avg_logFC"]
  updeg.FC[!updeg.ind[,i]==1,i] <- NA
}
updeg.FC <- updeg.FC[nrow(updeg.FC):1,]
updeg.FC <- reshape2::melt(updeg.FC)
ggplot(updeg.FC, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(low = "white",high = "#E64B35",na.value = "#d8cfcd",labels=seq(0,1.5, by = 0.5),breaks=seq(0,1.5, by = 0.5),limits=c(0,1.5))+
  scale_x_discrete(position = "top")
#Down
downdeg.sig <- lapply(deg.sig, function(x){x <- x[x$avg_logFC < -FCth,]})
downgene.names <- c() 
for (i in 1:length(deg.sig)) {
  downgene.names <- c(downgene.names, rownames(downdeg.sig[[i]]))
} # Extract unique gene names
downgene.names <- downgene.names[!duplicated(downgene.names)] # Remove duplicated gene names
downdeg.inx <- matrix(0,nrow = length(downgene.names), ncol = length(deg.sig),dimnames = (list(downgene.names, names(deg.sig)))) # Initialize the index matrix

for (i in 1:length(downdeg.sig)){
  downdeg.inx[rownames(downdeg.sig[[i]]),i] <- 1
} # Mark the gene as significant in index matrix
# 2. find the intersective genes among all gene set by row sum. Set a threshold: the row sum larger than 3.
downdeg.ind <- downdeg.inx
indsum <- apply(downdeg.ind,1,sum)
downdeg.ind <- cbind(downdeg.ind, sum=indsum)
downdeg.ind <- downdeg.ind[order(downdeg.ind[,ncol(downdeg.ind)], decreasing = T),]
# 3. Do heatmap on filtered table and fill with logFC.
downdeg.ind <- downdeg.ind[downdeg.ind[,ncol(downdeg.ind)]>threshold,]
downdeg.FC <- downdeg.ind[,-ncol(downdeg.ind)] # Initialize another matrix to store logFC value
for (i in 1:ncol(downdeg.FC)){
  downdeg.FC[downdeg.ind[,i]==1,i] <- downdeg.sig[[i]][rownames(downdeg.ind)[downdeg.ind[,i]==1],"avg_logFC"]
  downdeg.FC[!downdeg.ind[,i]==1,i] <- NA
}

downdeg.FC <- reshape2::melt(downdeg.FC)
ggplot(downdeg.FC, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(high  = "white",low = "steelblue",na.value = "#d8cfcd",labels=seq(-1,0, by = 0.5),breaks=seq(-1,0, by = 0.5),limits=c(-1,0)) +
  scale_x_discrete(position = "bottom")

###F1.i Unique DEG heatmap ####
deglist <- EC.ovy.dif[2:7]
# UP markers
deg.sig <- lapply(deglist, function(x){x <- x[x$p_val_adj <0.05,]})
updeg.sig <- lapply(deg.sig, function(x){x <- x[x$avg_logFC >FCth,]})
upgene.names <- c() 
for (i in 1:length(deg.sig)) {
  upgene.names <- c(upgene.names, rownames(updeg.sig[[i]]))
} # Extract unique gene names
upgene.names <- upgene.names[!duplicated(upgene.names)] # Remove duplicated gene names
updeg.inx <- matrix(0,nrow = length(upgene.names), ncol = length(deg.sig),dimnames = (list(upgene.names, names(deg.sig)))) # Initialize the index matrix

for (i in 1:length(updeg.sig)){
  updeg.inx[rownames(updeg.sig[[i]]),i] <- 1
} # Mark the gene as significant in index matrix
updeg.ind <- updeg.inx
indsum <- apply(updeg.ind,1,sum)
updeg.ind <- cbind(updeg.ind, sum=indsum)
updeg.ind <- updeg.ind[order(updeg.ind[,ncol(updeg.ind)], decreasing = T),]
updeg.ind <- updeg.ind[updeg.ind[,ncol(updeg.ind)]==1,]
updeg.sum <- colSums(updeg.ind[,-ncol(updeg.ind)])
# 2. Order this list by logFC value and use top3 most representative genes to do heatmap. 
upmarker.FC <- updeg.ind[,-ncol(updeg.ind)] # Copy a new matrix
for (i in 1:ncol(upmarker.FC)){
  upmarker.FC[updeg.ind[,i]==1,i] <- updeg.sig[[i]][rownames(updeg.ind)[updeg.ind[,i]==1],"avg_logFC"]
  upmarker.FC[!updeg.ind[,i]==1,i] <- NA
} # Get FC value from input data list
upmarker.FC <- as.data.frame(upmarker.FC)
upmarker.FCs <- data.frame()
expr <- paste0("upmarker.FC[with(upmarker.FC,order(",paste(colnames(upmarker.FC),collapse = ','),",decreasing = T)),]") # Reorder the data frame column by column
upmarker.FC <- eval(parse(text = expr))
for (i in 1:6){
  if (sum(updeg.ind[,i]==1)>topn){
    upmarker.FCs <- rbind(upmarker.FCs, upmarker.FC[updeg.ind[,i]==1,][1:topn,])}
  else{upmarker.FCs <- rbind(upmarker.FCs, upmarker.FC[updeg.ind[,i]==1,])}
}
# for (i in 1:ncol(upmarker.FCs)){
#   intergenes <- intersect(rownames(upmarker.FCs)[is.na(upmarker.FCs[,i])], rownames(deg.sig[[i]]))
#   if (length(intergenes)>0){
#     upmarker.FCs[intergenes,i] <-deg.sig[[i]][intergenes,"avg_logFC"]
#   }
# }
upmarker.FCs <- upmarker.FCs[nrow(upmarker.FCs):1,]
upmarker.FCs <- melt(as.matrix(upmarker.FCs[,1:5]))
ggplot(upmarker.FCs, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(low = "white",high = "#E64B35",na.value = "#d8cfcd",labels=seq(0,0.6, by = 0.2),breaks=seq(0,0.6, by = 0.2),limits=c(0,0.6))+
  scale_x_discrete(position = "top") 
# Down markers
downdeg.sig <- lapply(deg.sig, function(x){x <- x[x$avg_logFC < -FCth,]})
downgene.names <- c() 
for (i in 1:length(deg.sig)) {
  downgene.names <- c(downgene.names, rownames(downdeg.sig[[i]]))
} # Extract unique gene names
downgene.names <- downgene.names[!duplicated(downgene.names)] # Remove duplicated gene names
downdeg.inx <- matrix(0,nrow = length(downgene.names), ncol = length(deg.sig),dimnames = (list(downgene.names, names(deg.sig)))) # Initialize the index matrix

for (i in 1:length(downdeg.sig)){
  downdeg.inx[rownames(downdeg.sig[[i]]),i] <- 1
} # Mark the gene as significant in index matrix
# 2. find the intersective genes among all gene set by row sum. Set a threshold: the row sum larger than 3.
downdeg.ind <- downdeg.inx 
indsum <- apply(downdeg.ind,1,sum)
downdeg.ind <- cbind(downdeg.ind, sum=indsum)
downdeg.ind <- downdeg.ind[order(downdeg.ind[,ncol(downdeg.ind)], decreasing = T),]
downdeg.ind <- downdeg.ind[downdeg.ind[,ncol(downdeg.ind)]==1,]
downdeg.sum <- colSums(downdeg.ind[,-ncol(downdeg.ind)])
downmarker.FC <- downdeg.ind[,-ncol(downdeg.ind)] # Copy a new matrix
for (i in 1:ncol(downmarker.FC)){
  downmarker.FC[downdeg.ind[,i]==1,i] <- downdeg.sig[[i]][rownames(downdeg.ind)[downdeg.ind[,i]==1],"avg_logFC"]
  downmarker.FC[!downdeg.ind[,i]==1,i] <- NA
}
downmarker.FC <- as.data.frame(downmarker.FC)
downmarker.FCs <- data.frame()
expr <- paste0("downmarker.FC[with(downmarker.FC,order(",paste(colnames(downmarker.FC),collapse = ','),")),]") # Reorder the data frame column by column
downmarker.FC <- eval(parse(text = expr))
for (i in 1:ncol(downmarker.FC)){
  if (sum(downdeg.ind[,i]==1)>topn){
    downmarker.FCs <- rbind(downmarker.FCs, downmarker.FC[downdeg.ind[,i]==1,][1:topn,])}
  else{downmarker.FCs <- rbind(downmarker.FCs, downmarker.FC[downdeg.ind[,i]==1,])}
}
# for (i in 1:ncol(downmarker.FCs)){
#   intergenes <- intersect(rownames(downmarker.FCs)[is.na(downmarker.FCs[,i])], rownames(deg.sig[[i]]))
#   if (length(intergenes)>0){
#     downmarker.FCs[intergenes,i] <-deg.sig[[i]][intergenes,"avg_logFC"]
#   }
# }
downmarker.FCs <- melt(as.matrix(downmarker.FCs[,1:ncol(downmarker.FC)-1]))
ggplot(downmarker.FCs, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(high  = "white",low = "steelblue",na.value = "#d8cfcd",labels=seq(-0.6,0, by = 0.2),breaks=seq(-0.6,0, by = 0.2),limits=c(-0.6,0))#+scale_x_discrete(position = "top")
dev.off()
#_________________________________####
# Figure.2 ####

###F2.a Pathway enrichment results####
allmat.score <- read_csv("./inter_data/cleaned_score.csv") #Load the manually curated list

allmat.score$`SuperPath Name` <- factor(allmat.score$`SuperPath Name`,levels = rev(allmat.score$`SuperPath Name`))
allmat.score <- allmat.score[nrow(allmat.score):1,]
finalmat.up.score <- allmat.score[,c(2:7)]
colnames(finalmat.up.score)[2:6] <- c("aEC1","aEC2","CapEC","vCapEC","vEC")
finalmat.down.score <- allmat.score[,c(2,8:12)]
colnames(finalmat.down.score)[2:6] <- c("aEC1","aEC2","CapEC","vCapEC","vEC")
finalmat.up.score <- gather(finalmat.up.score, Celltype, expr,-"SuperPath Name")
finalmat.down.score <- gather(finalmat.down.score, Celltype, expr,-"SuperPath Name")

allmat.overlap <- read_csv("./inter_data/cleaned_overlap_new.csv") #Load the manually corrected list

allmat.overlap$`SuperPath Name` <- factor(allmat.overlap$`SuperPath Name`,levels = rev(allmat.overlap$`SuperPath Name`))
allmat.overlap <- allmat.overlap[nrow(allmat.overlap):1,]
finalmat.up.overlap <- allmat.overlap[,c(2:7)]
colnames(finalmat.up.overlap)[2:6] <- c("aEC1","aEC2","CapEC","vCapEC","vEC")
finalmat.down.overlap <- allmat.overlap[,c(2,8:12)]
colnames(finalmat.down.overlap)[2:6] <- c("aEC1","aEC2","CapEC","vCapEC","vEC")
finalmat.up.overlap <- gather(finalmat.up.overlap, Celltype, overlap,-"SuperPath Name")
finalmat.down.overlap <- gather(finalmat.down.overlap, Celltype, overlap,-"SuperPath Name")

finalmat.up <- inner_join(finalmat.up.score,finalmat.up.overlap,by=c("SuperPath Name","Celltype"))
finalmat.down <- inner_join(finalmat.down.score,finalmat.down.overlap,by=c("SuperPath Name","Celltype"))

pdf("./NVU_figure/F2/Figure_2_DEG_heatmap_shared.pdf",paper = "a4r")
ggplot(finalmat.up, aes(x=Celltype, y=finalmat.up$`SuperPath Name`,color=expr,size=overlap)) + 
  geom_point()+scale_x_discrete(position = "top")+scale_size_continuous(range = c(1,15),limits = c(1,20))+
  scale_color_gradient(high  = "black", low = "#E64B35",na.value = "#d8cfcd")+
  expand_limits(colour = seq(0, 30, by = 5))
ggplot(finalmat.down, aes(x=Celltype, y=finalmat.up$`SuperPath Name`,color=expr,size=overlap)) + 
  geom_point()+scale_size_continuous(range = c(1,15),limits = c(1,20))+
  scale_color_gradient(high  = "black", low = "steelblue",na.value = "#d8cfcd")+
  expand_limits(colour = seq(0, 30, by = 5))
dev.off()

###F2.b The background genes in pathway enrichment results####

#1. Read the discussed genes 
overlap.gene <- read_csv("./inter_data/Tobediscussed.csv",col_names = T)
overlap.gene <- overlap.gene %>% distinct()
#2. Load the pval and lnFC and reshape them to fit ggplot2
ovy.overlap.pval <- lapply(seq(2,7), function(i){
  x <- EC.ovy.dif[[i]]
  name.tmp <- names(EC.ovy.dif)[i]
  x <- as_tibble(x,rownames="Gene") %>% filter(Gene %in% overlap.gene$Gene) %>% select(Gene,!! name.tmp:=p_val_adj) 
})

ovy.overlap.pval <- purrr::reduce(ovy.overlap.pval,full_join,by="Gene")  %>% gather(key = "Subtype",value = "p_val_adj",-"Gene")

ovy.overlap.FC <- lapply(seq(2,7), function(i){
  x <- EC.ovy.dif[[i]]
  name.tmp <- names(EC.ovy.dif)[i]
  x <- as_tibble(x,rownames="Gene") %>% filter(Gene %in% overlap.gene$Gene) %>% select(Gene,!!name.tmp:=avg_logFC) 
})

ovy.overlap.FC <- purrr::reduce(ovy.overlap.FC,full_join,by="Gene")  %>% gather(key = "Subtype",value = "avg_logFC",-"Gene")

ovy.overlap.mat <- inner_join(ovy.overlap.FC, ovy.overlap.pval, by=c("Gene","Subtype"))

ovy.overlap.mat$Gene <- factor(ovy.overlap.mat$Gene, levels = rev(overlap.gene$Gene))
ovy.overlap.mat$Subtype <- as.factor(ovy.overlap.mat$Subtype)
new.levels <- c("aEC1","aEC2","avEC","CapEC","vEC","vCapEC")
names(new.levels) <- levels(ovy.overlap.mat$Subtype)
ovy.overlap.mat$Subtype <- factor(new.levels[ovy.overlap.mat$Subtype],levels = c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC"))
ovy.overlap.mat$p_val_adj <- (-1)*log2(ovy.overlap.mat$p_val_adj)
#3. Plot data
# ggplot(ovy.overlap.mat, aes(x=Subtype, y=Gene,color=avg_logFC,size=p_val_adj)) + 
#   geom_point()+scale_x_discrete(position = "top")+scale_size_continuous(range = c(1,12))+
#   scale_color_gradient2(high  = "#E64B35", low = "steelblue",mid = "white",na.value = "#d8cfcd")
#expand_limits(colour = seq(0, 30, by = 5))
# Plot seperately
ovy.overlap.mat.pos <- ovy.overlap.mat %>% filter(avg_logFC>0&Subtype %in% c("aEC1","aEC2","CapEC","vCapEC","vEC") & p_val_adj>4.322) %>% 
  mutate(logFC=ifelse(avg_logFC>0.5,0.5,avg_logFC), `-logP_val`=ifelse(p_val_adj>40,40,p_val_adj)) 
ovy.overlap.mat.neg <- ovy.overlap.mat %>% filter(avg_logFC<0&Subtype %in% c("aEC1","aEC2","CapEC","vCapEC","vEC") & p_val_adj>4.322) %>% 
  mutate(logFC=ifelse(avg_logFC<(-0.5),-0.5,avg_logFC), `-logP_val`=ifelse(p_val_adj>40,40,p_val_adj))

pdf("Figure_2_Pathway_genes_dotplot.pdf")
  ggplot(ovy.overlap.mat.pos, aes(x=Subtype, y=Gene,color=logFC,size=`-logP_val`)) + 
    geom_point()+scale_x_discrete(position = "top")+scale_size_continuous(range = c(1,12))+
    scale_color_gradient(high  = "#E64B35", low = "white",na.value = "#d8cfcd",limits=c(0,0.5))+theme_classic()+expand_limits(size=c(0,40))
  ggplot(ovy.overlap.mat.neg, aes(x=Subtype, y=Gene,color=logFC,size=`-logP_val`)) + 
    geom_point()+scale_x_discrete(position = "top")+scale_size_continuous(range = c(1,12))+
    scale_color_gradient(high  = "white", low = "#2b506e",na.value = "#d8cfcd",limits=c(-0.5,0))+theme_classic()+expand_limits(size=c(0,40))
dev.off()

###F2.c to e Not from sequencing data####

#_________________________________####
# Figure.3 ####


###F3.a GWAS Disease over-representation test dotplot####
dict2 <- readRDS("./dicts/dict2.rds")
deg.sig <- lapply(EC.ovy.dif[2:7], function(x){x <- x[(x$p_val_adj <0.05&abs(x$avg_logFC)>FCth),]})
dislist <- dict2[overlap] # Convert back to mmu PS: The overlap data is calculated by enrich_brain.R script (another individual analysis).
brain.mat.pval <- matrix(data = 0, nrow = 5, ncol = length(dislist), dimnames = list(names(deg.sig)[1:5],dislist))
for (i in 1:5){
  brain.mat.pval[i,] <- deg.sig[[i]][dislist,"p_val_adj"]
}
brain.mat.pval <- melt(brain.mat.pval)
brain.mat.pval <- brain.mat.pval %>% mutate(pval=(-log2(value))) %>% select(Var1,Var2,pval)
brain.mat.FC <- matrix(data = 0, nrow = 5, ncol = length(dislist), dimnames = list(names(deg.sig)[1:5],dislist))
for (i in 1:5){
  brain.mat.FC[i,] <- deg.sig[[i]][dislist,"avg_logFC"]
}
brain.mat.FC <- melt(brain.mat.FC)
colnames(brain.mat.FC)[3] <- "LnFC"
brain.mat <- inner_join(brain.mat.FC, brain.mat.pval, by=c("Var1","Var2"))
brain.mat$Var1 <- factor(brain.mat$Var1, levels = rev(c("A1","A2","Cap","VCap","V")))

pdf("./NVU_figure/F3/Figure_3_Brain_Disease_heatmap.pdf",width = 11.69,height = 5)
ggplot(brain.mat, aes(x=Var2, y=Var1,color=LnFC,size=pval)) +
  geom_point()+scale_size_continuous(range = c(1,10))+
  scale_color_gradient2(mid = "white", low = "steelblue",high="#E64B35",na.value = "#d8cfcd",breaks=seq(-0.4,0.6,0.2),
                        limits=c(-0.4,0.6),labels=seq(-0.4, 0.6, by = 0.2))+
  theme_classic()+expand_limits(size=c(0,50))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

###F3.b Human AD brain DEG compared to mouse EC DEG####

#The allen brain data process and comparsion with mouse data plot was finished by Allen_data.R

###F3.c Human AD brain DEG boxplot####

#The boxplot was finished by Allen_boxplot.R including the genes we discussed in paper and other concordant genes.

#_________________________________####
# Figure.4 ####

###F4.d The reversalability of exenatide on different gene sets ####
degdot <- function(x,y){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(x)) {
    df1 <- x[[i]]
    df2 <- y[[i]]
    its.genes <- intersect(rownames(df1),rownames(df2))
    sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
    sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
    sig=sig1+sig2
    df <- tibble(evo=-df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
                 sig=sig ,row.names = its.genes)
    df <- df %>% filter(sig>0) # Filter the unsignificant genes
    new.levels <- c("evo_sig","ovy_sig","both_sig")
    df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
    formula <- y ~ x
    print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
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
            expand_limits(x=c(-1, 1.5), y=c(-0.6,0.6))+
            geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 

pdf("./NVU_figure/F4/Figure_4_Ex_Subtypes.pdf")
degdot(EC.ove.dif,EC.ovy.dif)
dev.off()

degdot_enrich_Cap <- function(x,y,target){ # Compare the EVO and OVY effect based on a target gene list
  target <- target[!duplicated(target)]
  df1 <- x$Cap
  df2 <- y$Cap
  df1 <- df1[target,]
  df2 <- df2[target,]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=-df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>1) # Filter the unsignificant genes
  new.levels <- c("evo_sig","ovy_sig","both_sig")
  df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
  formula <- y ~ x
  #write.csv(df,file=paste0(names(x)[i],".csv"))
  print(ggplot(df, aes(x=ovy,y=evo,label=row.names))+ geom_point(size=3,alpha = 0.7,colour="#325d81")+
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,ypos = 0.4)+geom_text_repel()+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0("CapEC","_Subtype"))+theme_bw()+
          stat_fit_glance(method = "lm", 
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*", 
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          expand_limits(x=c(-0.5, 0.5), y=c(-0.3,0.6))+ 
          scale_y_continuous(breaks=c(-0.3,0,0.3,0.6), limits=c(-0.3,0.6))+
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0))
} 

pdf("./NVU_figure/F4/Figure_4_Ex_Pathway_BBB.pdf")
degdot_enrich_Cap(x = EC.ove.dif,y=EC.ovy.dif,target = BBB.genes)
dev.off()

pdf("./NVU_figure/F4/Figure_4_Ex_Pathway_IMM.pdf")
degdot_enrich_Cap(x = EC.ove.dif,y=EC.ovy.dif,target = IMM.genes)
dev.off()

pdf("./NVU_figure/F4/Figure_4_Ex_Pathway_ENG.pdf")
degdot_enrich_Cap(x = EC.ove.dif,y=EC.ovy.dif,target = ENG.genes)
dev.off()

###F4.e1 The reversalability of exenatide on AD relevent gene sets ####

pdf("./NVU_figure/F4/Figure_4_Ex_AD.pdf")
degdot_enrich_Cap(x = EC.ove.dif,y=EC.ovy.dif,target = AD.genes)
dev.off()

###F4.e2 The reversalability of exenatide on AD and GTEX concordent genes ####

degdot_enrich_same_dir <- function(x,y,target){ # Compare the EVO and OVY effect based on a target gene list
  df1 <- x$EC
  df2 <- y$EC
  df1 <- df1[target$Gene,]
  df2 <- df2[target$Gene,]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=-df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  new.levels <- c("evo_sig","ovy_sig","both_sig")
  df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
  df <- inner_join(df,target,by=c("row.names"="Gene"))
  formula <- y ~ x
  df[33,"row.names"] <- "C1orf54"
  print(ggplot(df, aes(x=ovy,y=evo,label=row.names))+ geom_point(size=3,alpha = 0.8,aes(colour=Source))+
          stat_smooth_func(geom="text",method="lm",hjust=0,parse=T,ypos = 0.4)+geom_text_repel()+theme_bw()+
          expand_limits(x=c(-0.6, 0.6), y=c(-0.3,0.6))+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+ggtitle("Same_direction")+
          stat_fit_glance(method = "lm", 
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*", 
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  
} 
#PS.Criteria 1. hsa_sig=T, 2. mouse <- 0.05+0.1 3.Same direction 4.LogFC>0.7 in EC cells 5.Mmu use CapEC as FC value

#GTEX data
same.di <- bind_rows((final.data.gt %>% filter(sig_hsa==TRUE) %>% filter(mmuFC>0 & hsaFC >0)),(final.data.gt %>% filter(sig_hsa==TRUE) %>% filter(mmuFC<0 & hsaFC <0)))
same.di <- same.di %>% filter(sig_EC==T) %>% .[["Gene"]]
gpro.conv <- read_csv("./dicts/gpro_gtex.csv")
same.di <- c("CAVIN2","ADGRF5", same.di) # Update some symbol names to lastet version
same.gtex <- gpro.conv %>% filter(ortholog_name %in% same.di) %>% select(input) %>% .[['input']]

#Allen AD data
same.di <- bind_rows((final.data.al %>% filter(sig_hsa==TRUE) %>% filter(mmuFC>0 & hsaFC >0)),(final.data.al %>% filter(sig_hsa==TRUE) %>% filter(mmuFC<0 & hsaFC <0)))
same.di <- same.di %>% filter(sig_mmu==T) %>% .[["Gene"]]
same.allen <- gpro.conv %>% filter(ortholog_name %in% same.di) %>% select(input) %>% .[['input']]
same.all <- tibble(Gene=c(same.gtex, same.allen)) %>% mutate(Source=ifelse(Gene %in% same.gtex, "GTEx","Allen"))

pdf("./NVU_figure/F4/Figure_4_Ex_Same_dir_All.pdf")
degdot_enrich_same_dir(x = EC.ove.dif,y=EC.ovy.dif,target = same.all)
dev.off()

#_________________________________####
# Figure.S1 ####

###S1.a Violin plot for major cell types####

Idents(ovy.oy) <- "Celltype"
pdf("./NVU_figure/S1/Figure_1S_Cluster_markers.pdf")
VlnPlot(object = ovy.oy, features = c("Kcnj8"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Acta2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Cldn5"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Ctss"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Ntsr2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Syt1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Cdk1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Pdgfra"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Cldn11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Sox11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Ccdc153"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Ttr"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Alas2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Pf4"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Rax"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy.oy, features = c("Plac8"),pt.size = 0,ncol = 1,y.max = 8)
dev.off()

###S1.b Heatmap for EC markers expression####

EC.data <- subset(ovy.final,idents = "EC")
EC.file <- "./cellassign/EC.csv" # EC markers
EC.marker <- list() # Prepare markers
EC.tmp <-read.csv(file = EC.file,header = T)
for (i in (1:6)){
  EC.marker[i] <- strsplit(as.character(EC.tmp$gene.symbols[i]),"\\s+")
}

names(EC.marker) <- c("A1","A2","V","AV","Cap","VCap")

EC.marker <- EC.marker[c(1,2,5,6,3,4)]
remgenes <- function (x) {
  if (length(x) > 50) {
    short.genes <- x[1:50]
    
  }
  else {
    short.genes <- x
  }
  return(short.genes)
} #Remove some genes from list to decrease computating work
EC.marker <- lapply(EC.marker, remgenes)

EC.marker <- lapply(EC.marker, function(x) x[x %in% rownames(EC.data)]) # Remove the non-exist gene names

subtypeexpr <-function (x, smarker) {
  mean(unlist(FetchData(EC.data, vars = x,cells = WhichCells(object = EC.data, idents = smarker),slot = "data")))
}
Idents(EC.data) <- "Subtype"
EC.sig.genes <-c()
for (i in 1:6) { # Reorder the gene list by expression level
  #EC.marker <- EC.marker[intersect(rownames(EC.data), rownames(EC.marker)),]
  marker.list <- EC.marker[[i]]
  ECmarker.exp <- lapply((marker.list),subtypeexpr,smarker = names(EC.marker)[i]) 
  
  marker.scored <- data.frame(unlist(ECmarker.exp), row.names = (marker.list))
  
  marker.sorted <- marker.scored[order(marker.scored$unlist.ECmarker.exp.,decreasing = T),,drop=FALSE]
  
  EC.sig.genes <-c(EC.sig.genes,rownames(marker.sorted)[1:17])
}
EC.data <- ScaleData(EC.data,features = EC.sig.genes)
Idents(EC.data,cells=WhichCells(EC.data,idents = "AV")) <- "AV"
Idents(EC.data,cells=WhichCells(EC.data,idents = "V")) <- "V"
Idents(EC.data,cells=WhichCells(EC.data,idents = "VCap")) <- "VCap"
Idents(EC.data,cells=WhichCells(EC.data,idents = "Cap")) <- "Cap"
Idents(EC.data,cells=WhichCells(EC.data,idents = "A2")) <- "A2"
Idents(EC.data,cells=WhichCells(EC.data,idents = "A1")) <- "A1"
EC.data$ECSubtype <- Idents(EC.data)
EC.data <- SubsetData(EC.data,ident.remove = "EC_unstable")
heat.obj<- DoHeatmap(EC.data, features = EC.sig.genes, group.by = "ECSubtype",label = F,draw.lines = F,slot = "scale.data")


###S1.c Subtype marker gene expression level between old and young####

expr.mat <- FetchData(ovy.clean,vars = c("Mfsd2a","Tfrc","Bmx","Vegfc","Nr2f2","Slc38a5","Vwf","Vcam1"))
expr.mat <- as_tibble(expr.mat,rownames = "cell.bc")
Subtype.marker <- as_tibble(as.data.frame(ovy.clean$Subtype),rownames = "cell.bc")
Age.marker <- as_tibble(as.data.frame(ovy.clean$Age),rownames = "cell.bc")
expr.mat <- inner_join(Subtype.marker,expr.mat,by="cell.bc")
expr.mat <- inner_join(Age.marker, expr.mat, by="cell.bc")

expr.mat <- expr.mat %>% filter(`ovy.clean$Age` %in% c("old","young") & `ovy.clean$Subtype` %in% c("A1","A2","Cap","VCap","V","AV"))
expr.mat <- expr.mat %>% group_by(`ovy.clean$Age`,`ovy.clean$Subtype`) %>% summarise_at(.funs = mean,.vars = vars(colnames(expr.mat)[4:11]))
expr.mat <- gather(expr.mat, marker, expr,-c(`ovy.clean$Age`,`ovy.clean$Subtype`))
expr.mat <- mutate(expr.mat, area=case_when(marker=="Bmx" | marker=="Vegfc" ~ "Artery",
                                            marker=="Mfsd2a" | marker=="Tfrc" ~ "Capillary",
                                            marker=="Nr2f2"| marker=="Slc38a5" ~ "Vein",
                                            marker=="Vwf" | marker=="Vcam1" ~ "Arterial_vein"))
expr.mat$`ovy.clean$Subtype` <- factor(expr.mat$`ovy.clean$Subtype`, levels = c("A1","A2","Cap","VCap","V","AV"))
expr.mat$area <- factor(expr.mat$area, levels = c("Artery","Capillary","Vein","Arterial_vein"))
#expr.mat <- expr.mat %>% group_by(marker, `ovy.clean$Subtype`) %>% summarise(expr =mean(expr))
pdf("./NVU_figure/S1/Figure_S1_Markers_level.pdf")
ggplot(data = expr.mat,aes(y = expr, x = `ovy.clean$Subtype`, group=marker)) + 
  geom_line(aes(colour = marker),stat="identity",size=1)+ 
  geom_point(aes(colour = marker),size=1.5) + scale_color_npg()+
  theme_bw()+facet_grid(area ~ `ovy.clean$Age`)
dev.off()
#_________________________________####
# Figure.S2 ####

###S2.a No. of DEGs vesus No. of cells####
DEGs.no <- lapply(EC.ovy.dif, function(x) dim(x[x$p_val_adj < 0.05 & abs(x$avg_logFC)>0.1,])[1]) # Extract sig DEGs number first
DEGs.no <- unlist(DEGs.no)
DEGs.no <- enframe(DEGs.no,name = 'Celltype',value = 'DEGno')

DEGs.cell.young <- subtype.young[,1:2] %>% bind_rows(tibble(Var1='EC',Freq=sum(subtype.young$Freq)))
DEGs.cell.old <- subtype.old[,1:2] %>% bind_rows(tibble(Var1='EC',Freq=sum(subtype.old$Freq)))

DEGs.cell.young <- left_join(DEGs.cell.young,DEGs.no,by=c('Var1'='Celltype'))
DEGs.cell.old <- left_join(DEGs.cell.old,DEGs.no,by=c('Var1'='Celltype'))

DEGs.cell.all <- tibble(Var1=DEGs.cell.old$Var1,Freq=DEGs.cell.old$Freq+DEGs.cell.young$Freq,DEGno=DEGs.cell.old$DEGno) #Construct a total number

pdf(file = './NVU_figure/S2/DEGs_vs_Cell_number_All.pdf')
ggplot(DEGs.cell.all[1:6,],aes(x=Freq,y=DEGno,label=Var1))+
  geom_point()+stat_smooth(method = 'lm',se = FALSE)+stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,ypos = 0.4)+
  geom_text_repel()+theme_classic()+
  stat_fit_glance(method = "lm",method.args = list(formula = formula),label.x = "right",
                  label.y = "bottom",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE)
dev.off()

###S2.b Venn diagram between pooled EC DEG and merged subtype DEGs####
tt.EC.DEG <- as_tibble(EC.ovy.dif[[1]],rownames = "Gene")
tt.EC.DEG <- tt.EC.DEG %>% filter((p_val_adj<0.05)&abs(avg_logFC)>0.1) %>% select(Gene)

mg.EC.DEG <- lapply(EC.ovy.dif[2:7], as_tibble,rownames="Gene")
mg.EC.DEG <- lapply(mg.EC.DEG, function(x)filter(x,(p_val_adj<0.05)&abs(avg_logFC)>0.1))
mg.EC.DEG <- lapply(mg.EC.DEG, function(x)x$Gene)
mg.EC.DEG <- unlist(mg.EC.DEG)
mg.EC.DEG <- mg.EC.DEG[!duplicated(mg.EC.DEG)]

DEG.list <- list(`Pooled EC`=tt.EC.DEG$Gene,Merged=mg.EC.DEG)

pdf("./NVU_figure/S2/Venn.pdf")
venn(DEG.list, zcolor = c("#F9FA9B","#FF7777") ,cexil = 1, cexsn = 0.8,borders = T,opacity = 0.4)
dev.off()

###S2.c Correlation between pooled EC and EC subtypes####
inter.DEG <- intersect(tt.EC.DEG$Gene,mg.EC.DEG)

tt.EC.list <- as_tibble(EC.ovy.dif[[1]],rownames = "Gene") %>% filter(Gene %in% inter.DEG) %>% select(Gene,avg_logFC)

mg.EC.list <- lapply(EC.ovy.dif[2:7], as_tibble,rownames="Gene")
mg.EC.list <- lapply(1:6,function(x){
  mg.EC.list[[x]] %>% mutate(Subtype=names(mg.EC.list)[x])
})

mg.EC.list <- do.call(bind_rows, mg.EC.list)
mg.EC.list <- mg.EC.list %>% filter(Gene %in% inter.DEG) %>%filter((p_val_adj<0.05)&abs(avg_logFC)>0.1) %>% 
  select(Gene,avg_logFC,Subtype)
mg.EC.list <- mg.EC.list %>% rename(mg_lnFC=avg_logFC)
tt.EC.list <- tt.EC.list %>% rename(tt_lnFC =avg_logFC)

EC.list <- merge(tt.EC.list,mg.EC.list,by="Gene")
EC.list$Subtype <- factor(EC.list$Subtype,levels = c("A1","A2","Cap","VCap","V","AV"))

pdf("./NVU_figure/S2/Cor_TtEC&mergeEC.pdf")
ggplot(EC.list, aes(x=mg_lnFC,y=tt_lnFC))+ geom_point(size=2,alpha = 0.8,aes(colour=Subtype))+
  geom_segment(aes(x = -1, y = -1, xend = 1.5, yend = 1.5))+
  scale_color_manual(values = c("#E64B35FF","#F39B7FFF","#00A087FF","#227c91","#3C5488FF","#8491B4FF"))+theme_bw()+
  geom_hline(yintercept = 0) +geom_vline(xintercept = 0)+coord_fixed()
dev.off()
#_________________________________####
# Figure.S3 ####

###S3. The Vcam1 gene in avEC####

pdf("./NVU_figure/S3/Vcam1.pdf")
VCAM1.AV.dif <- FindMarkers(AV.data,ident.1 = "old",ident.2 = "young",features = "Vcam1",logfc.threshold = 0,test.use = 't')
VlnPlot(AV.data,features = "Vcam1",idents = c("young","old"))
VCAM1.pos.av <- subset(AV.data,subset = Vcam1 > 0)
VCAM1.AV.pos.dif <- FindMarkers(VCAM1.pos.av,ident.1 = "old",ident.2 = "young",features = "Vcam1",logfc.threshold = 0,test.use = 't')
VlnPlot(VCAM1.pos.av,features = "Vcam1",idents = c("young","old"))
dev.off()

#_________________________________####
# Figure.S5 ####

###S5.a The comparision analsyis combined with Mayo dataset####
#5. Mayo RNAseq dataset
#5.1 Mayo DEG with Allen data
Mayo.AD.DEG <- read_tsv('../NVU_revision/Out_source/Mayo_TCX/MayoRNAseq_RNAseq_TCX_ADvsCON_Comprehensive.txt') #Load the Mayo DEG result directly
Mayo.AD.DEG.sig <- Mayo.AD.DEG %>% filter(Dx.qValue<0.05)
Mayo.AD.DEG.sig.enrich.ovysig <- Mayo.AD.DEG.sig %>% filter(GeneName %in% gpro.conv$ortholog_name)

# Check the concordance
final.data.al <- final.data.al %>% select(Gene, mmuFC, hsaFC, Enrich_logFC, sig_mmu) %>% mutate(Source="Allen")

Mayo.final <- left_join(Mayo.AD.DEG.sig.enrich.ovysig,mmu.data,by=c('GeneName'='Symbol')) %>% 
  filter(Ind==TRUE) %>% select(GeneName,avg_logFC,Dx.Beta,Enrich_logFC,Ind) %>% mutate(Source="Mayo")

colnames(Mayo.final) <- colnames(final.data.al)
Mayo.final <- bind_rows(Mayo.final, final.data.al)

pdf("./NVU_figure/S5/Figure_S5_Allen_Mayo.pdf")
ggplot(Mayo.final, aes(y=hsaFC,x=mmuFC,label=Gene))+ 
  geom_point(size=4,alpha = 0.8,aes(colour=Enrich_logFC))+ggrepel::geom_text_repel(size=4)+
  scale_color_gradient(low="#71a0c8",high = "black")+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  expand_limits(x=c(-0.5, 0.5), y=c(-0.3,0.6))+scale_y_continuous(breaks = c(-0.3,-0.15,0,0.15,0.3,0.45,0.6))+
  theme_classic()+coord_fixed(ratio = 1.2)
dev.off()

###S5.b The concordant gene in AD human brain ####

Mayo.patient <- read_csv('../NVU_revision/Out_source/Mayo_TCX/MayoRNAseq_RNAseq_TCX_covariates.csv')
Mayo.mat <- read_tsv('../NVU_revision/Out_source/Mayo_TCX/MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv',)
Mayo.4box <- Mayo.mat %>% filter(ensembl_id %in% c('ENSG00000164741','ENSG00000137962','ENSG00000140443','ENSG00000140575'))

Mayo.4box <- tibble(Patient= colnames(Mayo.4box), DLC1=t(Mayo.4box[1,]),ARHGAP29=t(Mayo.4box[2,]),IGF1R=t(Mayo.4box[3,]),IQGAP1=t(Mayo.4box[4,]))
Mayo.4box <- Mayo.4box[-1,]

Mayo.ctrl <- Mayo.patient %>% filter(Diagnosis=="Control") %>% select(SampleID)
Mayo.AD <- Mayo.patient %>% filter(Diagnosis=="AD") %>% select(SampleID)

Mayo.4box <- Mayo.4box %>% mutate(AD=ifelse(Patient %in% Mayo.ctrl$SampleID, 'Ctrl','AD'))

Mayo.4box$AD <- factor(Mayo.4box$AD,levels = c('Ctrl','AD'))
Mayo.4box$DLC1 <- as.numeric(Mayo.4box$DLC1)
Mayo.4box$ARHGAP29 <- as.numeric(Mayo.4box$ARHGAP29)
Mayo.4box$IGF1R <- as.numeric(Mayo.4box$IGF1R)
Mayo.4box$IQGAP1 <- as.numeric(Mayo.4box$IQGAP1)


pdf(file = './NVU_figure/S5/Mayo_box.pdf')

g <- ggboxplot(Mayo.4box, x = "AD", y = "DLC1",
               palette = "npg",add = "jitter",color="black",title = 'DLC1',
               legend = "right",add.params = list(color = "grey"),
               width = 0.4,ylim=c(0,120))+scale_y_continuous(breaks = c(0,30,60,90,120))
print(g+ stat_compare_means(method = "t.test"))

g <- ggboxplot(Mayo.4box, x = "AD", y = "ARHGAP29",
               palette = "npg",add = "jitter",color="black",title = 'ARHGAP29',
               legend = "right",add.params = list(color = "grey"),
               width = 0.4,ylim=c(0,150))+scale_y_continuous(breaks = c(0,30,60,90,120,150))
print(g+ stat_compare_means(method = "t.test"))

g <- ggboxplot(Mayo.4box, x = "AD", y = "IGF1R",
               palette = "npg",add = "jitter",color="black",title = 'IGF1R',
               legend = "right",add.params = list(color = "grey"),
               width = 0.4,ylim=c(0,280))+scale_y_continuous(breaks = c(0,70,140,210,280))
print(g+ stat_compare_means(method = "t.test"))

g <- ggboxplot(Mayo.4box, x = "AD", y = "IQGAP1",
               palette = "npg",add = "jitter",color="black",title = 'IQGAP1',
               legend = "right",add.params = list(color = "grey"),
               width = 0.4,ylim=c(0,120))+scale_y_continuous(breaks = c(0,30,60,90,120))
print(g+ stat_compare_means(method = "t.test"))
dev.off()

#_________________________________####
# Figure.S7 ####

###S7.a The correlation between GTEX DEG and mouse DEG in EC enriched genes ####

#Please check GTEx_data.R to finish analysis

###S7.b The slope of regression line between GTEX DEG and mouse DEG in EC enriched genes ####

#Please check Update_Gtex.R to finish analysis

###S7.c The boxplot of concordant and discordant genes in GTEX dataset ####

#Please check Gtex_boxplot.R to finish analysis

###S7.d The regression result of selected genes ####

#Please check Update_Gtex.R to finish analysis

#_________________________________####
# Figure.S8 ####

###S8. The exenatide reversal effect on other EC subtypes####

pdf("./NVU_figure/S8/Figure_8_Ex_Subtypes.pdf")
degdot(EC.ove.dif,EC.ovy.dif)
dev.off()

#_________________________________####
# Figure.S9 ####

###S9.c The tSNE plot of selected MG markers

Idents(ovy.clean) <- "Celltype"
MG.plt.dat <- subset(ovy.clean,idents="MG")
Idents(MG.plt.dat) <- "seurat_clusters"
MG.plt.dat <- subset(MG.plt.dat,idents=c("1","3","9","29"))
Idents(MG.plt.dat) <- "Age"

old.bc <- MG.plt.dat$Age %>% enframe(name = 'cellbc',value = 'age') %>% filter(age=="old") 
old.bc <- old.bc$cellbc
young.bc <- MG.plt.dat$Age %>% enframe(name = 'cellbc',value = 'age') %>% filter(age=="young") 
young.bc <- sample(young.bc$cellbc,2839)
oldex.bc <- MG.plt.dat$Age %>% enframe(name = 'cellbc',value = 'age') %>% filter(age=="oldex") 
oldex.bc <- sample(oldex.bc$cellbc,2839)

MG.old <- subset(MG.plt.dat,cells=old.bc)
MG.young <- subset(MG.plt.dat,cells=young.bc)
MG.oldex <- subset(MG.plt.dat,cells=oldex.bc)

mg.marker <- c("Apoe","Ccl6","Cd9","Timp2")
pdf("./NVU_figure/S9/microglia_tSNE.pdf",width = 5,height = 15)
for (i in 1:length(mg.marker)) {
  old.plt <- FeaturePlot(MG.old,features = mg.marker[i],min.cutoff = 0, max.cutoff = 4)+
    scale_y_continuous(limits = c(-10,25),breaks=seq(-10,25,5))+scale_x_continuous(limits = c(-27,3),breaks=seq(-25,0,5))
  
  young.plt <- FeaturePlot(MG.young,features = mg.marker[i],min.cutoff = 0, max.cutoff = 4)+
    scale_y_continuous(limits = c(-10,25),breaks=seq(-10,25,5))+scale_x_continuous(limits = c(-27,3),breaks=seq(-25,0,5))
  
  oldex.plt <- FeaturePlot(MG.oldex,features = mg.marker[i],min.cutoff = 0, max.cutoff = 4)+
    scale_y_continuous(limits = c(-10,25),breaks=seq(-10,25,5))+scale_x_continuous(limits = c(-27,3),breaks=seq(-25,0,5))
  col.plt <- cowplot::plot_grid(young.plt,old.plt,oldex.plt,ncol=1)
  print(col.plt)
}
dev.off()

###S9.d The expression level of selected genes in agiang and exenatide treated brains

mg.marker <- read_csv('~/inter_data/partial_markers_final.csv')
mg.marker <- mg.marker[!duplicated(mg.marker$Gene),]



mgmak.ovy <- as_tibble(MG.ovy.diff,rownames = "Gene")
mgmak.ove <- as_tibble(MG.ove.diff,rownames = "Gene")
mgmak.ove$avg_logFC <- mgmak.ove$avg_logFC*-1

mgmak.ovy <- inner_join(mgmak.ovy,mg.marker,by="Gene")
mgmak.ove <- inner_join(mgmak.ove,mg.marker,by="Gene")

mgmak.ovy <- mgmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
mgmak.ove <- mgmak.ove %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

mgmak.dat <- bind_rows(mgmak.ovy,mgmak.ove)
mgmak.dat <- mgmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-40),40,(-1*log2pval)))
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC>0.4,0.4,mgmak.dat$avg_logFC)
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC<(-0.2),(-0.2),mgmak.dat$avg_logFC)

mgmak.dat$Type <- factor(mgmak.dat$Type,levels = c("OVY","EVO"))
mgmak.dat$Gene <- factor(mgmak.dat$Gene,levels = c("Ptprc","Cd68","Aif1","Timp2","Lgals3bp","Cd52","Cd9","Ccl6","Ccl4","Ccl3","Ccl2","Trem2","Apoe"))
mgmak.dat$Gene <- fct_recode(mgmak.dat$Gene,Iba1="Aif1",Cd45="Ptprc")

pdf("./NVU_figure/S9/microglia_markers.pdf")
ggplot(mgmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,40))+ 
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.2,0.41))
dev.off()