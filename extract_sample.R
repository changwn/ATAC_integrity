
setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac')
load('C:/Users/wnchang/Documents/F/PhD_Research/Functions8.1.RData')  #edit TCGA colnames function

sample_info <- read.csv('DataS1_si.csv', header = T, skip=37)

ATAC_sample_id <- as.character(unique(sample_info[, 'submitter_id']) )
ATAC_sample_id <- gsub('-', '_', ATAC_sample_id)
#ATAC_sample_id <- as.character(unique(sample_info[, 'Tissue_Barcode']) )
ATAC_cancer <- as.character(unique(sample_info[,'cohort']))

sample_info[, 'submitter_id'] <- gsub('-','_',as.character(sample_info[, 'submitter_id']))  

ATAC_TCGA_common <- list()
cancer_str <- c()
pointer <- 1
#tmp <- c()
for(i in 1:length(ATAC_cancer))
{
  if(ATAC_cancer[i] %in% c('CHOL','MESO', 'THCA', 'UCEC') ) next
  print(ATAC_cancer[i])
  data_chr <- paste("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/","TCGA-",ATAC_cancer[i],"_FPKM_T.RData",sep='' )
  load(data_chr)
  
  TCGA_exist <- edit_colnames3(colnames(data_t))
  overlap_n <- length(intersect(ATAC_sample_id, TCGA_exist))
  print(overlap_n)
  
  if(length(intersect(ATAC_sample_id, TCGA_exist)) > 0)
  {
    
    #tmp <- data.frame(cancer=rep(ATAC_cancer[i], overlap_n), sample=intersect(ATAC_sample_id, TCGA_exist))
    #ATAC_TCGA_common <- rbind(ATAC_TCGA_common, tmp)
    common_sp <- intersect(ATAC_sample_id, TCGA_exist)
    ATAC_id <- list()
    for(j in 1:length(common_sp)){
      ATAC_id[[j]] <- as.character(sample_info[sample_info[,'submitter_id']==common_sp[j], 1])
    }
    names(ATAC_id) <- common_sp
    cancer_str <- c(cancer_str, ATAC_cancer[i])
    ATAC_TCGA_common[[pointer]] <- ATAC_id
    pointer <- pointer+1
  }
  
  
}
names(ATAC_TCGA_common) <- cancer_str



save(list=c('ATAC_sample_id','ATAC_cancer','ATAC_TCGA_common'), file='ATAC_extract_02.RData')





