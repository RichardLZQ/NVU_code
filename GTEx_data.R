library(dplyr)
library(stringr)
library(ggpubr)
library(readr)

#1. Load the data downloaded from https://www.gtexportal.org/home/datasets
Ph <- read.table("~/bioinfo/Gtex_test/Data/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",header = T,sep = "\t")  #Load phenotype and annotation info
Attr <- read.table("~/bioinfo/Gtex_test/Data/GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

options(stringsAsFactors = F)
GTEx=read.table('~/bioinfo/Gtex_test/Data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz',header = T,sep = '\t',skip = 2)

#2. Extract the cortex data 

ID.list <- Attr[grepl("Brain",Attr$SMTS),]

ID.list <- ID.list[!(grepl("Cerebellum",ID.list$SMTSD) | grepl("cervical",ID.list$SMTSD)),1]
ID.list <- as.character(lapply(ID.list, function(x) gsub('[-]','.',x)))

Cortex.data <-cbind(GTEx[,1:2],GTEx[,colnames(GTEx) %in% ID.list]) #Keep cortex data first

#3. Prepare mmu data 
# EC.enrich <- FindMarkers(object = ovy.final, ident.1 = "EC",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
# saveRDS(EC.enrich,file = "./inter_data/EC_enriched_genes.rds")
EC.enrich <- readRDS("./inter_data/EC_enriched_genes.rds") #Use this to accelerate the analysis
EC.enrich.sig <- EC.enrich[(EC.enrich$avg_logFC)>0.7,] # Define the really enriched genes
EC.enrich.sig <- as_tibble(EC.enrich.sig,rownames = "Gene")
mer.list <- vector()
for (i in 2: length(EC.ovy.dif)){ #Keep the sig DEGs
  tmp <- EC.ovy.dif[[i]]
  tmp <- tmp[(tmp$p_val_adj <0.05 & abs(tmp$avg_logFC)>0.1),]
  mer.list <- c(mer.list, rownames(tmp))
} # Merged list
mer.list <- mer.list[!duplicated(mer.list)]
mer.list <- mer.list[mer.list %in% rownames(EC.ovy.dif[[1]])] 
mer.tibble <- tibble(Gene=mer.list) 
ec.enrich.ind <- tibble(Ind=(mer.list %in% EC.enrich.sig$Gene), Symbol = mer.list)
ec.enrich.FC <- left_join(mer.tibble,EC.enrich.sig,by="Gene") %>% rename(Symbol=Gene,Enrich_logFC=avg_logFC) %>% select(Symbol,Enrich_logFC)
mmu.data <- EC.ovy.dif[[1]][mer.list,]

#4. Load the convert dic

gpro.conv <- read_csv("./dicts/gpro_gtex.csv")
gpro.conv <- gpro.conv %>% filter(ortholog_name!="N/A" & input %in% mer.list) #
gpro.conv <- gpro.conv[!duplicated(gpro.conv$input),]

mmu.data <- mmu.data[gpro.conv$input,] #Clean wrong genes
mmu.data <- as_tibble(mmu.data,rownames = NA) %>% tibble::rownames_to_column(var = "Symbol") %>% mutate(ENSG=gpro.conv$ortholog_ensg)

Cortex.data$Name <- str_sub(Cortex.data$Name,1,15)
mmu.data$ENSG[!mmu.data$ENSG %in% Cortex.data$Name] # Confirm all converted genes are in Gtex database

mmu.data <- mmu.data %>% filter(ENSG!="ENSG00000273749"& ENSG!="ENSG00000283586")
Cortex.data <- Cortex.data %>% filter(Name %in% mmu.data$ENSG) %>% tibble::column_to_rownames(var = "Name") 
Cortex.data <- Cortex.data[mmu.data$ENSG,]

# 5. Add the patient info 

ID.list2 <- colnames(Cortex.data)[2:length(colnames(Cortex.data))]

Age.list <- Ph
Age.list$SUBJID <- as.character(lapply(Age.list$SUBJID, function(x) gsub('[-]','.',x)))
Age.ind <- unlist(lapply(Age.list$SUBJID, function(x) sum(grepl(x, ID.list2))>0))
Age.list <- Age.list[Age.ind,]
rownames(Age.list) <- Age.list$SUBJID

Age <- as.character(lapply(ID.list2,function(x) Age.list$AGE[str_detect(x,rownames(Age.list))]))
Age <- c("NA",Age)
Age <- as.factor(Age)
new.levels <- c("Young","Young","Old","Old","NA")
Age <- factor(new.levels[Age],levels = c("Young","Old","NA"))
names(Age) <- colnames(Cortex.data)

Cortex.t <- t(Cortex.data)
Cortex.t <- data.frame(Age, Cortex.t)

# 6. Calculate the change p-value (qvalue) + boxplot
p.table <- tibble()
for (i in 1:nrow(Cortex.data)){
  data <- (data.frame(Cortex.t[,c("Age",rownames(Cortex.data)[i])]))
  gene.name <- data[1,2]
  data <- data[2:nrow(data),]
  colnames(data) <- c('V1','V2')
  data$V2 <- unlist(lapply(data$V2,as.numeric ))
  old.mean <- data %>% filter(V1=="Old") %>% summarise(mean=mean(V2))
  young.mean <- data %>% filter(V1=="Young") %>% summarise(mean=mean(V2))
  tmp<- compare_means(V2~V1, data = data, method = "t.test",paired = F)
  tmp <- mutate(tmp, logFC = log(as.numeric(old.mean/young.mean)),Gene = gene.name)
  p.table <- bind_rows(p.table, tmp)
  if (i%%1000==0){print(i)}
}

p.table <- mutate(p.table, p.adj = p.adjust(p.table$p, method = "fdr", n = length(p.table$p)))

p.table.sig <- p.table %>% filter(p.adj<0.05)
mmu.data <- inner_join(x = mmu.data, y = ec.enrich.FC, by= "Symbol")
mmu.data <- inner_join(x = mmu.data, y = ec.enrich.ind, by= "Symbol")

final.data.gt <- tibble(Gene=p.table$Gene,mmuFC=mmu.data$avg_logFC, hsaFC=p.table$logFC, sig_hsa= ifelse(p.table$Gene %in% p.table.sig$Gene,TRUE,FALSE), Enrich_logFC= mmu.data$Enrich_logFC,sig_EC= mmu.data$Ind)

final.data.gt <- final.data.gt %>% filter((Enrich_logFC)>0.7 & sig_hsa==TRUE) #Final.data object is calculated by Gtex.final script

pdf("../NVU_figure/F4/Figure_4_Gtex_coordinate_0909.pdf")
ggplot(final.data.gt, aes(y=hsaFC,x=mmuFC,label=Gene))+ 
  geom_point(size=4,alpha = 0.8,aes(colour=Enrich_logFC))+ggrepel::geom_text_repel(size=4)+
  scale_color_gradient(low="#71a0c8",high = "black",)+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  expand_limits(x=c(-0.6, 0.6), y=c(-0.4,1.2))+scale_y_continuous(breaks = c(-0.4, 0,0.4,0.8,1.2))+
  theme_classic()+coord_fixed(ratio = 0.75)
dev.off()