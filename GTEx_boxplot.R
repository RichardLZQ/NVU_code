#Gtex boxplot script####

library(dplyr)
library(stringr)
library(ggpubr)

#1. Load the data
Ph <- read.table("~/bioinfo/Gtex_test/Data/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",header = T,sep = "\t")  #Load phenotype and annotation info
Attr <- read.table("~/bioinfo/Gtex_test/Data/GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

options(stringsAsFactors = F)
GTEx=read.table('~/bioinfo/Gtex_test/Data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz',header = T,sep = '\t',skip = 2)


overlap.aging <- read_xlsx("./inter_data/Boxplot_genes/Boxpolot_genes.xlsx")  %>% filter(Source=="GTEx") %>% select(mmu_genes)
#2. Extract the cortex data 

ID.list <- Attr[grepl("Brain",Attr$SMTS),]

ID.list <- ID.list[!(grepl("Cerebellum",ID.list$SMTSD) | grepl("cervical",ID.list$SMTSD)),1]
#ID.list <- ID.list[(grepl("Cortex",ID.list$SMTSD)),1]
ID.list <- as.character(lapply(ID.list, function(x) gsub('[-]','.',x)))

Cortex.data <- cbind(GTEx[,1:2],GTEx[,colnames(GTEx) %in% ID.list]) #Keep cortex data first

#3. Load the convert dic

gpro.conv <- read_csv("~/bioinfo/Gtex_test/gpro_gtex.csv")

gpro.conv <- gpro.conv %>% filter(input %in% overlap.aging$mmu_genes) #
gpro.conv <- gpro.conv[!duplicated(gpro.conv$input),]
Cortex.data <- Cortex.data %>% filter(Description %in% gpro.conv$ortholog_name)

#4. Add the patient info 

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
colnames(Cortex.t)[2:ncol(Cortex.t)] <- Cortex.t[2,2:ncol(Cortex.t)]

#5. Plot boxplot
setwd("../NVU_figure/F4/Gtex/")
for (i in 1:length(Cortex.data$Description)){
  data <- (data.frame(Cortex.t[,c("Age",Cortex.data$Description[i])]))
  gene.name <- data[2,2]
  data <- data[3:nrow(data),]
  colnames(data) <- c('V1','V2')
  #data$V1 <- unlist(lapply(data$V1,as.numeric ))
  data$V2 <- unlist(lapply(data$V2,as.numeric ))
  uplimit <- ceiling(max(data$V2)/10)*10
  pdf(paste0(gene.name,".pdf"))
  g <- ggboxplot(data, x = "V1", y = "V2",
                 palette = "npg",add = "jitter",color="black",title = gene.name,
                 legend = "right",add.params = list(color = "grey"),
                 width = 0.4,ylim=c(0,uplimit))
  print(g+ stat_compare_means(method = "t.test"))
  dev.off()
}
setwd("../../../NVU_code/")
