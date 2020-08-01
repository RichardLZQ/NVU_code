#1. Load the data downloaed from https://aging.brain-map.org/api/v2/well_known_file_download/502999992
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
#EC.enrich <- FindMarkers(object = ovy.final, ident.1 = "EC",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
EC.enrich <- readRDS("~/bioinfo/NVU_code/inter_data/EC_enriched_genes.rds")
EC.enrich.sig <- EC.enrich[(EC.enrich$avg_logFC)>0.7,] # Define the really enriched genes
EC.enrich.sig <- as_tibble(EC.enrich.sig,rownames = "Gene")
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


gpro.conv <- read_csv("~/bioinfo/NVU_code/dicts/gpro_allen.csv")

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

p.table <- tibble()
for (i in 1:nrow(mmu.data)) {
  data <- data.frame(allen.t[,c("AD_stat", mmu.data$Symbol[i])])
  colnames(data) <- c('V1','V2')
  data$V1 <- unlist(lapply(data$V1,as.factor ))
  data$V2 <- unlist(lapply(data$V2,as.numeric ))
  logFC <- data %>% group_by(V1) %>% summarise(mean=mean(V2)) %>% arrange(V1)
  tmp<- compare_means(V2~V1, data = data, method = "t.test",paired = F)
  tmp <- mutate(tmp, logFC = as.numeric(log(logFC[2,2]/logFC[1,2])),Gene = mmu.data$Symbol[i])
  p.table <- bind_rows(p.table, tmp)
}
p.table <- mutate(p.table, p.adj =p.adjust(p.table$p, method = "fdr", n = length(p.table$p)))

mmu.data <- mmu.data[mmu.data$Symbol %in% p.table$Gene,]
mmu.data <- inner_join(x = mmu.data, y = ec.enrich.FC, by= "Input")


mmu.data <- inner_join(x = mmu.data, y = ec.enrich.ind, by= "Input")

p.table.sig <- p.table %>% filter(p.adj<0.05)

final.data <- tibble(Gene=p.table$Gene,mmuFC=mmu.data$avg_logFC, hsaFC=p.table$logFC, sig_hsa = ifelse(p.table$Gene %in% p.table.sig$Gene,TRUE,FALSE), Enrich_logFC = mmu.data$Enrich_logFC,sig_mmu = mmu.data$Ind)

final.data.al <- final.data %>% filter((Enrich_logFC)>0.7 & sig_hsa==TRUE) #Final.data object is calculated by Gtex.final script
#final.data.p <- mutate(final.data.p, glabel =  ifelse(sig_hsa==TRUE,Gene,""))
pdf("../NVU_figure/F4/Figure_4_Allen_coordinate_0903.pdf")
ggplot(final.data.al, aes(y=hsaFC,x=mmuFC,label=Gene))+ 
  geom_point(size=4,alpha = 0.8,aes(colour=Enrich_logFC))+ggrepel::geom_text_repel(size=4)+
  scale_color_gradient(low="#71a0c8",high = "black")+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  expand_limits(x=c(-0.5, 0.5), y=c(-0.3,0.3))+scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))+
  theme_classic()+coord_fixed(ratio = 1.2)
#xlim(-0.5, 0.5)+ylim(-0.4,1.2)
dev.off()