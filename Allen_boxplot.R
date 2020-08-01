# Script to generate boxplot on allen brain AD dataset


#1. Load the data
allen.data<- read_csv("~/bioinfo/Allen_data/raw/fpkm_table_normalized.csv")
sample.info <- read_csv("~/bioinfo/Allen_data/raw/columns-samples.csv")
gene.annotation <- read_csv("~/bioinfo/Allen_data/raw/rows-genes.csv")
patient.info <- read_csv("~/bioinfo/Allen_data/raw/DonorInformation.csv")

# 2. Prepare the annotation relationship
sample.region <- sample.info #%>% filter(structure_acronym=="PCx" )#| structure_acronym=="TCx"

ctrl <- patient.info %>% filter(ever_tbi_w_loc=="N" & dsm_iv_clinical_diagnosis == "No Dementia") %>% filter(donor_id %in% sample.region$donor_id)

AD <- patient.info %>% filter(ever_tbi_w_loc=="N" & dsm_iv_clinical_diagnosis == "Alzheimer's Disease Type") %>% filter(donor_id %in% sample.region$donor_id)

ctrl.sd <- sample.region %>% filter(donor_id %in% ctrl$donor_id)

AD.sd <- sample.region %>% filter(donor_id %in% AD$donor_id)

allen.ad <- allen.data[,(colnames(allen.data)%in%ctrl.sd$rnaseq_profile_id | colnames(allen.data)%in% AD.sd$rnaseq_profile_id)]  

allen.ad <- bind_cols(Gene=gene.annotation$gene_symbol, as_tibble(allen.ad))
# 3. Prepare the target gene list
mer.list <- vector()
for (i in 2: length(EC.ovy.dif)){
  tmp <- EC.ovy.dif[[i]]
  tmp <- tmp[(tmp$p_val_adj <0.05& abs(tmp$avg_logFC)>0.1),]
  mer.list <- c(mer.list, rownames(tmp))
}
mer.list <- mer.list[!duplicated(mer.list)]
mer.list <- mer.list[mer.list %in% rownames(EC.ovy.dif[[1]])] 
mer.tibble <- tibble(Gene=mer.list) 
ec.enrich.ind <- tibble(Ind=(mer.list %in% EC.enrich.sig$Gene), Input = mer.list)
ec.enrich.FC <- left_join(mer.tibble,EC.enrich.sig,by="Gene") %>% rename(Input=Gene,Enrich_logFC=avg_logFC) %>% select(Input,Enrich_logFC)
mmu.data <- EC.ovy.dif[[1]][mer.list,]


gpro.conv <- read_csv("~/bioinfo/Allen_data/gpro.csv")

gpro.conv <- gpro.conv %>% filter(ortholog_name!="N/A" & input %in% mer.list) #

gpro.conv <- gpro.conv[!duplicated(gpro.conv$input),]


mmu.data <- mmu.data[gpro.conv$input,]
mmu.data <- as_tibble(mmu.data,rownames = NA) %>% tibble::rownames_to_column(var = "Input") %>%  mutate(Symbol =gpro.conv$ortholog_name)
rownames(mmu.data) <- gpro.conv$ortholog_name

# 4. Transform allen data and combined with old data
allen.ad <- allen.ad %>% tibble::column_to_rownames(var = "Gene")
mmu.data <- mmu.data[intersect(rownames(mmu.data),rownames(allen.ad)),]
allen.ad <- allen.ad[mmu.data$Symbol,]

AD.status<- sapply(colnames(allen.ad), function(x) ifelse(x %in% ctrl.sd$rnaseq_profile_id, "non-AD","AD"))

allen.ad <- rbind(AD.status, allen.ad)
rownames(allen.ad)[1] <- "AD_stat" # Ad AD status into allen dataset

allen.t <- t(allen.ad)

# Load target genes
library(readxl)
overlap.aging <- c("SLC2A1","MFSD2A","IFITM3","ATOX1","NOSTRIN")

for (i in 1:length(overlap.aging)) {
  data <- data.frame(allen.t[,c("AD_stat", overlap.aging[i])])
  colnames(data) <- c('V1','V2')
  data$V1 <- unlist(lapply(data$V1,as.factor ))
  data$V2 <- unlist(lapply(data$V2,as.numeric ))
  gene.name <- overlap.aging[i]
  uplimit <- ceiling(max(data$V2)/10)*10
  pdf(paste0(gene.name,".pdf"))
  g <- ggboxplot(data, x = "V1", y = "V2",
                 palette = "npg",add = "jitter",color="black",title = gene.name,
                 legend = "right",add.params = list(color = "grey"),
                 width = 0.4,ylim=c(0,uplimit))
  print(g+ stat_compare_means(method = "t.test"))
  dev.off()
}

