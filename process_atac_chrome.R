
setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac')
load('./Homo_sapiens.GRCh38_91_info_CHR.RDAta')

table(Homo_sapiens.GRCh38_91_info[,5])
Homo_prot_code <- Homo_sapiens.GRCh38_91_info[which(Homo_sapiens.GRCh38_91_info[, 5] == 'protein_coding'), ]

Homo_prot_code[, 9] <- as.character(Homo_prot_code[, 9])

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts')
file_norm <- list.files()

atac_norm_tmp <- read.table(file_norm[1], header=T)


for(ii in 1:nrow(atac_norm_tmp))
{

	chr_tmp <- as.character(atac_norm_tmp[ii, 1])
	start_tmp <- atac_norm_tmp[ii, 2]
	end_tmp <- atac_norm_tmp[ii, 3]

	found_count <- judge_in_range(chr_tmp, start_tmp, end_tmp, Homo_prot_code)


}


#-------------function------------------
judge_in_range <- function(chr_tmp, start_tmp, end_tmp, Homo_prot_code)
{
	chr_tmp <- unlist(strsplit(chr_tmp, 'r', fixed=T))[[2]]
	Homo_prot_code_tmp <- Homo_prot_code[which(Homo_prot_code[, 9] == chr_tmp), ]
	found_flag <- 0
	for(i in 1:nrow(Homo_prot_code_tmp))
	{
		if(Homo_prot_code_tmp[i, 6]<=start_tmp && end_tmp<=Homo_prot_code_tmp[i, 7])
		{
			found_flag <- found_flag + 1
		}

	}

	return(found_flag)
}


# setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_PeakCalls')
# file_peak <- list.files()
# atac_peak_tmp <- read.csv(file_peak[1], header=T, sep='\t')














