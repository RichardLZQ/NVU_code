# Cellassign method with Tensorflow platform
library(SingleCellExperiment)
library(cellassign)
library(tensorflow)
library(scran)


#0. Load the data ####
#EC <- readRDS("./inter_data/EC.data") # Load the EC and SMC data
EC <- computeSumFactors(EC) # Calc the size factor
EC.s <- sizeFactors(EC)

#1. Loading the cell markers. ####

EC.file <- "./cellassign/EC.csv"

EC.marker <- list() # Prepare markers
EC.tmp <-read.csv(file = EC.file,header = T)
for (i in (1:6)){
  EC.marker[i] <- strsplit(as.character(EC.tmp$gene.symbols[i]),"\\s+")
}

names(EC.marker) <- c("A1","A2","V","AV","Cap","VCap")

# 2. Prepare the markers ####

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

EC.marker <- lapply(EC.marker, function(x) x[x %in% rownames(EC)]) # Remove the non-exist gene names

EC.marker <- marker_list_to_mat(EC.marker,include_other = F) # Transform marker to required matrix

# 3. Run the model
set.seed(2019)                # set R rng seed
reticulate::py_set_seed(2019) # set python rng seed

EC.fit <- cellassign(exprs_obj = t(as.matrix(SummarizedExperiment::assay(EC[rownames(EC.marker)], "counts"))), 
                     marker_gene_info = EC.marker, 
                     s = EC.s, 
                     learning_rate = 1e-2, 
                     shrinkage = TRUE,
                     verbose = T,
                     num_runs = 3)

Sys.time()

# For SMC data
#SMC <- readRDS("./inter_data/SMC.data")
SMC <- computeSumFactors(SMC)
SMC.s <- sizeFactors(SMC)
SMC.marker <- readRDS("./cellassign/smc_marker.rds")

SMC.marker <- lapply(SMC.marker, function(x) x[x %in% rownames(SMC)])

SMC.marker <- marker_list_to_mat(SMC.marker,include_other = F) # Transform marker to required


SMC <- as.matrix(SummarizedExperiment::assay(SMC[rownames(SMC.marker)], "counts"))

SMC.s <- SMC.s[!(colSums(SMC)==0)]
SMC <- SMC[,!(colSums(SMC)==0)]

set.seed(2019)                # set R rng seed
reticulate::py_set_seed(2019) # set python rng seed
SMC.fit <- cellassign(exprs_obj = t(SMC),
                      marker_gene_info = SMC.marker,
                      s = SMC.s,
                      learning_rate = 1e-2,
                      shrinkage = TRUE,
                      verbose = FALSE,
                      num_runs=3
)
Sys.time()

# 4. Summarize the results and return to main code
EC.type <- EC.fit$cell_type
names(EC.type) <- colnames(EC)

SMC.type <- SMC.fit$cell_type
names(SMC.type) <- colnames(SMC)

saveRDS(EC.type,"./assign/EC_results.rds")
saveRDS(SMC.type,"./assign/SMC_results.rds")
