# The code generating NVUB5 project figures

# a. The cell type umap
# a1. UMAP of selected celltypes
Idents(ovy5) <- "Celltype"
ovy5.less <- subset(ovy5, idents = c("AC","OPC","MG","MAC","OLG","PC","SMC","EC")) 
ovy5.less <- subset(ovy5, cells= sample(colnames(ovy5.less),6000)) #Downsample to 6000 cells
ovy5.less$Celltype <- droplevels(ovy5.less$Celltype)
ovy5.less$Celltype <- factor(ovy5.less$Celltype, levels = c("EC","OLG","OPC","AC","MG","MAC","SMC","PC"))
DimPlot(ovy5.less,group.by = "Celltype",label=T)

# a2. UMAP seperated by sample ids
Idents(ovy5.less) <- "orig.ident"
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "old5"),pt.size = 0.2,sizes.highlight = 0.3)
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "oldEx5"),pt.size = 0.2,sizes.highlight = 0.3)
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "young5"),pt.size = 0.2,sizes.highlight = 0.3)

# a3. The cell number
Orig.num <- enframe(ovy5$orig.ident,"cellbc","orig")
celltype.num <- enframe(ovy5$Celltype,"cellbc","celltype")

# b. The reversal effect coordinate plot

degdot.final <- function(x,y,ct.list){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(ct.list)) {
    ct <- ct.list[i]
    df1 <- x[[ct]]
    df2 <- y[[ct]]
    its.genes <- intersect(rownames(df1),rownames(df2))
    sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
    sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
    sig=sig1+sig2
    df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
                 sig=sig ,row.names = its.genes)
    df <- df %>% filter(sig>0) # Filter the unsignificant genes
    new.levels <- c("evo_sig","ovy_sig","both_sig")
    df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
    #df <- df %>% filter(sig!="evo_sig")
    formula <- y ~ x
    print(nrow(df))
    return(nrow(df))
    # print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
    #         stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
    #         geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
    #         ggtitle(paste0(ct.list[i],"_Celltype"))+theme_bw()+
    #         stat_fit_glance(method = "lm",
    #                         method.args = list(formula = formula),
    #                         label.x = "right",
    #                         label.y = "bottom",
    #                         aes(label = paste("italic(P)*\"-value = \"*",
    #                                           signif(..p.value.., digits = 4), sep = "")),
    #                         parse = TRUE)+
    #         scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
    #         geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 

ct.list <- c("PC","SMC","EC","MG","AC","OPC","OLG","MAC")
pdf("~/bioinfo/Richard/NVUB5/Figures/Figure/b. Cooridnate/allcell4patent.pdf")
degdot.final(evo.deg,ovy.deg,ct.list = ct.list)
dev.off()


# c. The slope and reversal propotion line plot
slope.list <- c()
cor.dataframe <- data.frame()
for (i in 1:length(ct.list)) { # Calculate the linear regression and retrive slope 
  ct <- ct.list[i]
  df1 <- x[[ct]]
  df2 <- y[[ct]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  new.levels <- c("evo_sig","ovy_sig","both_sig")
  df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
  #df <- df %>% filter(sig!="evo_sig")
  #df <- df %>% filter(ovy>-1 &ovy<1) %>% filter(evo> -0.75 &evo<0.75) #To make sure manually calculated number is same with automatical results
  linearMod <- lm(evo ~ ovy, data=df)
  cor.dataframe <- rbind(cor.dataframe,t(as.data.frame(linearMod$coefficients)))
  slope.list <- c(slope.list,linearMod$coefficients[2])
}

rownames(cor.dataframe) <-  c("PC","SMC","EC","MG","AC","OPC","OLG","MAC")
names(slope.list) <-  c("PC","SMC","EC","MG","AC","OPC","OLG","MAC")
slope.list <- enframe(slope.list, name = "Celltype","Slope")
slope.list <- slope.list[c(5,6,4,8,7,2,1,3),]
line.data <- left_join(slope.list,rev.porp,by="Celltype")
line.data$Celltype <- factor(line.data$Celltype,levels = line.data$Celltype)
#line.data$Slope <- line.data$Slope*-1
#line.data$Slope[9] <- 0 # Force negative to 0
#line.data$bk_remove[9] <- 0
line.data$group[4] <- "Glia"
line.data <- line.data[,-5]

line.data2 <- gather(line.data, "Type","Value",-c("Celltype","group")) 
#line.data2[9:,4] <- line.data2[9:24,4]*0.5 # Reshapre the data to fit "overlap" plot
line.data2$group2 <- c(pull(line.data2,group)[1:8],paste0(pull(line.data2,group)[1:8],"2"))
#line.data2$Value[1:8] <- ((line.data2$Value[1:8]-0.1)/0.9)*0.3+0.7
line.data2$Value[1:8] <- line.data2$Value[1:8]/2
line.data2$Value[9:16] <-  line.data2$Value[9:16]-0.5

ggplot(data=line.data2, aes(x=Celltype, y=Value, group=group2))+
  geom_line(aes(colour=group,linetype=Type))+
  geom_point(aes(colour=group,shape=Type))+scale_y_continuous(
    name = "Pct of DEGs reversed", limits = c(0,0.5),
    sec.axis = sec_axis( trans=~./2, name="Negative slope",breaks = seq(-1,-0.2,0.2))# Add a second axis and specify its features
  ) 

# d. The bubble plot of selected genes

# d1. Microglia
mg.marker <- read_csv("./Outsource/Selected_gene/MG_final.csv")
marker.type <- mg.marker$Type
mg.marker <- mg.marker[!duplicated(mg.marker$Gene),]

mgmak.ovy <- as_tibble(ovy.deg$MG,rownames = "Gene")
mgmak.evo <- as_tibble(evo.deg$MG,rownames = "Gene")

mgmak.ovy <- inner_join(mgmak.ovy,mg.marker,by = "Gene")
mgmak.evo <- inner_join(mgmak.evo,mg.marker,by = "Gene")

mgmak.ovy <- mgmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
mgmak.evo <- mgmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

mgmak.dat <- bind_rows(mgmak.ovy,mgmak.evo)
mgmak.dat$p_val_adj <- ifelse(mgmak.dat$p_val_adj==0,1e-320,mgmak.dat$p_val_adj)
mgmak.dat <- mgmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC>0.3,0.3,mgmak.dat$avg_logFC)
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC<(-0.3),(-0.3),mgmak.dat$avg_logFC)

mgmak.dat$Type <- factor(mgmak.dat$Type,levels = c("OVY","EVO"))
mgmak.dat$Gene <- factor(mgmak.dat$Gene,levels = rev(mg.marker$Gene))

ggplot(mgmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+ 
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))


#d2. Astrocyte

ac.marker <- read_csv("./Outsource/Selected_gene/AC_final.csv")
ac.marker <- ac.marker[!duplicated(ac.marker$Gene),]

acmak.ovy <- as_tibble(ovy.deg$AC,rownames = "Gene")
acmak.evo <- as_tibble(evo.deg$AC,rownames = "Gene")

acmak.ovy <- inner_join(acmak.ovy,ac.marker,by = "Gene")
acmak.evo <- inner_join(acmak.evo,ac.marker,by = "Gene")

acmak.ovy <- acmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
acmak.evo <- acmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

acmak.dat <- bind_rows(acmak.ovy,acmak.evo)
acmak.dat$p_val_adj <- ifelse(acmak.dat$p_val_adj==0,1e-320,acmak.dat$p_val_adj)
acmak.dat <- acmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
acmak.dat$avg_logFC <- ifelse(acmak.dat$avg_logFC<(-0.3),(-0.3),acmak.dat$avg_logFC)
acmak.dat$avg_logFC <- ifelse(acmak.dat$avg_logFC>(0.3),(0.3),acmak.dat$avg_logFC)

acmak.dat$Type <- factor(acmak.dat$Type,levels = c("OVY","EVO"))
acmak.dat$Gene <- factor(acmak.dat$Gene, levels = rev(ac.marker$Gene))

ggplot(acmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))


#d3. SMC

smc.marker <- read_csv("./Outsource/Selected_gene/SMC_final.csv")
smc.marker <- smc.marker[!duplicated(smc.marker$Gene),]

smcmak.ovy <- as_tibble(ovy.deg$SMC,rownames = "Gene")
smcmak.evo <- as_tibble(evo.deg$SMC,rownames = "Gene")

smcmak.ovy <- inner_join(smcmak.ovy,smc.marker,by = "Gene")
smcmak.evo <- inner_join(smcmak.evo,smc.marker,by = "Gene")

smcmak.ovy <- smcmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
smcmak.evo <- smcmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

smcmak.dat <- bind_rows(smcmak.ovy,smcmak.evo)
smcmak.dat <- smcmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
smcmak.dat$avg_logFC <- ifelse(smcmak.dat$avg_logFC>0.3,0.3,smcmak.dat$avg_logFC)
smcmak.dat$avg_logFC <- ifelse(smcmak.dat$avg_logFC<(-0.3),(-0.3),smcmak.dat$avg_logFC)

smcmak.dat$Type <- factor(smcmak.dat$Type,levels = c("OVY","EVO"))
smcmak.dat$Gene <- factor(smcmak.dat$Gene, levels = rev(smc.marker$Gene))

ggplot(smcmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))


# Seperate by type

mg.marker <- read_csv("./Outsource/Selected_gene/MG_marker.csv")
marker.type <- mg.marker$Type[!duplicated(mg.marker$Type)]
mg.marker <- mg.marker %>% filter(Type==marker.type[3])
mg.marker <- mg.marker[!duplicated(mg.marker$Gene),]

mgmak.ovy <- as_tibble(ovy.deg$MG,rownames = "Gene")
mgmak.evo <- as_tibble(evo.deg$MG,rownames = "Gene")

mgmak.ovy <- inner_join(mgmak.ovy,mg.marker,by = "Gene")
mgmak.evo <- inner_join(mgmak.evo,mg.marker,by = "Gene")

mgmak.ovy <- mgmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
mgmak.evo <- mgmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

mgmak.dat <- bind_rows(mgmak.ovy,mgmak.evo)
mgmak.dat$p_val_adj <- ifelse(mgmak.dat$p_val_adj==0,1e-320,mgmak.dat$p_val_adj)
mgmak.dat <- mgmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-40),(-1*log2pval),(-1*log2pval)))
#mgmak.dat <- mgmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-40),40,(-1*log2pval)))
# mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC>0.4,0.4,mgmak.dat$avg_logFC)
# mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC<(-0.2),(-0.2),mgmak.dat$avg_logFC)

mgmak.dat$Type <- factor(mgmak.dat$Type,levels = c("OVY","EVO"))
mgmak.dat$Gene <- factor(mgmak.dat$Gene,levels = rev(mg.marker$Gene))

ggplot(mgmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+ 
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab")



# 3 in 1 plots
ac.marker <- read_csv("./Outsource/Selected_gene/AC_final.csv")
mg.marker <- read_csv("./Outsource/Selected_gene/MG_final.csv")
smc.marker <- read_csv("./Outsource/Selected_gene/SMC_final.csv")

thr.marker <- bind_rows(ac.marker[,c(1,2)],mg.marker, smc.marker)

ac.marker <- ac.marker[!duplicated(ac.marker$Gene),]
acmak.ovy <- as_tibble(ovy.deg$AC,rownames = "Gene")
acmak.evo <- as_tibble(evo.deg$AC,rownames = "Gene")
acmak.ovy <- inner_join(acmak.ovy,ac.marker,by = "Gene")
acmak.evo <- inner_join(acmak.evo,ac.marker,by = "Gene")
acmak.ovy <- acmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="AC_OVY") %>% filter(p_val_adj!="NA")
acmak.evo <- acmak.evo %>% select(c(1,3,6)) %>% mutate(Type="AC_EVO") %>% filter(p_val_adj!="NA")
acmak.dat <- bind_rows(acmak.ovy,acmak.evo)

mg.marker <- mg.marker[!duplicated(mg.marker$Gene),]
mgmak.ovy <- as_tibble(ovy.deg$MG,rownames = "Gene")
mgmak.evo <- as_tibble(evo.deg$MG,rownames = "Gene")
mgmak.ovy <- inner_join(mgmak.ovy,mg.marker,by = "Gene")
mgmak.evo <- inner_join(mgmak.evo,mg.marker,by = "Gene")
mgmak.ovy <- mgmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="MG_OVY") %>% filter(p_val_adj!="NA")
mgmak.evo <- mgmak.evo %>% select(c(1,3,6)) %>% mutate(Type="MG_EVO") %>% filter(p_val_adj!="NA")
mgmak.dat <- bind_rows(mgmak.ovy,mgmak.evo)

smcmak.ovy <- as_tibble(ovy.deg$SMC,rownames = "Gene")
smcmak.evo <- as_tibble(evo.deg$SMC,rownames = "Gene")
smcmak.ovy <- inner_join(smcmak.ovy,smc.marker,by = "Gene")
smcmak.evo <- inner_join(smcmak.evo,smc.marker,by = "Gene")
smcmak.ovy <- smcmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="SMC_OVY") %>% filter(p_val_adj!="NA")
smcmak.evo <- smcmak.evo %>% select(c(1,3,6)) %>% mutate(Type="SMC_EVO") %>% filter(p_val_adj!="NA")
smcmak.dat <- bind_rows(smcmak.ovy,smcmak.evo)

thr.dat <- bind_rows(acmak.dat,mgmak.dat,smcmak.dat)
thr.dat <- thr.dat %>% mutate(`log2pval` = log2(p_val_adj))
thr.dat <- thr.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-300),300,(-1*log2pval)))
thr.dat$avg_logFC <- ifelse(thr.dat$avg_logFC>0.4,0.4,thr.dat$avg_logFC)
thr.dat$avg_logFC <- ifelse(thr.dat$avg_logFC<  -0.4,-0.4,thr.dat$avg_logFC)

thr.dat$Type <- factor(thr.dat$Type,levels = c("AC_OVY","AC_EVO","MG_OVY","MG_EVO","SMC_OVY","SMC_EVO"))
thr.dat$Gene <- factor(thr.dat$Gene,levels = rev(thr.marker$Gene[!duplicated(thr.marker$Gene)]))

ggplot(thr.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,5),limits = c(0,300))+
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.41,0.41))
