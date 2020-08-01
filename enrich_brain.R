# enrichBrain is an over reprentative analysis function finished by KHlab to fullfill the requirement that doing analysis focusing on brain diseases
  library(purrr)
  library(xlsx)
  library(stringr)
  library(purrr)
  
  # 1. Impoort DEGlist and brain diseases list
  FCthreshold=0.1
  deglist=EC.ovy.dif
  braindislist="./dicts/GWASSummary_Richard_Final(less).xlsx"
  deglist$Merge <- do.call(rbind, deglist[2:7]) # Load the 
  
  brain <- loadWorkbook(braindislist)
  brain <- getSheets(brain)
  disnames <- names(brain) # Load diseases names
  brain <- list()
  for (i in 2:length(disnames)) {
    tmp <- read.xlsx(braindislist, i,header = F ,stringsAsFactors=F)
    tmp <- unlist(tmp)# Add second columns to list
    tmp <- as.vector(na.omit(tmp))
    brain[[i-1]] <- tmp
  }
  names(brain) <- disnames[2:length(disnames)] 
  
  # 2. Doing hyper geometric analysis
  
  results <- list()
  
  for(i in 1:length(deglist)) {
    if (i ==length(deglist) ) { MergeDEGs <- deglist[[length(deglist)]]
    MergeDEGs <- rownames(MergeDEGs[(MergeDEGs$p_val_adj < 0.05 & abs(MergeDEGs$avg_logFC)>FCthreshold),])
    MergeDEGs <- strsplit(MergeDEGs,split = "\\.")
    DEGs <- rapply(MergeDEGs, function(x) return(x[2]))
    }
    else {DEGs <- rownames(deglist[[i]][(deglist[[i]]$p_val_adj < 0.05 & abs(deglist[[i]]$avg_logFC)>FCthreshold),])} # Extract gene names with signiciant
    DEGs.hsa<- mmudic(DEGs)
    subres <- data.frame()
    for (j in 1:length(brain)) {
      
      overlap <- intersect(DEGs.hsa, brain[[j]])
      q = length(overlap)
      m = length(brain[[j]])
      n = 19746 - m #The total detected gene number
      k = length(DEGs.hsa)
      pval <- phyper(q = q-1, m = m, n = n,k = k, lower.tail = FALSE)
      
      tmpdf <- DataFrame(p_value = pval, overlap_gene =paste(overlap,collapse = "/"),row.names = names(brain)[j])
      
      subres <- rbind(subres, tmpdf)
      
      
    }
    subres$adj_p <- p.adjust(subres$p_value, method = "fdr")
    results[[i]] <- subres
    
  
  names(results) <- names(deglist)
  # 3. Output the reuslts list 
  overlap <- lapply(results[2:5], function(x){x@listData[["overlap_gene"]]})
  
  overlap <- as.character(unlist(overlap))
  
  overlap <- overlap[nchar(overlap)>0]
  
  overlap <- paste(overlap, collapse = '/')
  
  overlap <- unlist(strsplit(overlap,split = '/'))
  
  overlap <- overlap[!duplicated(overlap)]
  
  dis.ind<- sapply(overlap, function(y)detect_index(brain, function(x) sum(str_detect(x,paste0("^",y,"$")))>=1))
  overlap.ind <- tibble(Gene=overlap, Disease=names(brain)[dis.ind])
  overlap.ind$Disease <- factor(overlap.ind$Disease,levels = c("SVD_WMH","Stroke_SVD","PD","AD","ALS","FTD","MSA"))
  overlap.ind <- arrange(overlap.ind,Disease)
  overlap <- overlap.ind$Gene
  
 
  return(overlap)
  
}
dict <- readRDS("./dicts/dict.rds")
mmudic <- function(DEGs){
  res <- dict[DEGs]
  ind <- sapply(dict[DEGs], function(x)length(x)>0)
  res <- res[ind]
  res <- as.character(unlist(res))
  return(res)
}
