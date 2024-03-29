
setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac')
load('./Homo_sapiens.GRCh38_91_info_CHR.RDAta')

table(Homo_sapiens.GRCh38_91_info[,5])
Homo_prot_code <- Homo_sapiens.GRCh38_91_info[which(Homo_sapiens.GRCh38_91_info[, 5] == 'protein_coding'), ]

Homo_prot_code[, 9] <- as.character(Homo_prot_code[, 9])

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts')
file_norm <- list.files()

extract_atac <- list()
cancer_str <- c()
extract_2more <- list()
for(i in 1:length(file_norm))
{
	setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts')
	atac_norm_tmp <- read.table(file_norm[i], header=T)
	cancer_str <- c(cancer_str, unlist(strsplit(file_norm[i], '_', fixed=T))[[1]])
	print(cancer_str[i])
	

	found1_count <- 0
	found2_count <- 0
	storage_info <- data.frame()
	more_than2_storage <- c()
	k <- 1
	#for(ii in 1:1000)
	for(ii in 1:nrow(atac_norm_tmp))
	{
		if(ii %% 3000 == 0) print(ii)
		#print(ii)
		chr_tmp <- as.character(atac_norm_tmp[ii, 1])
		start_tmp <- atac_norm_tmp[ii, 2]
		end_tmp <- atac_norm_tmp[ii, 3]

		found_stat <- judge_in_range(chr_tmp, start_tmp, end_tmp, Homo_prot_code)
		found_count <- found_stat[[1]]
		if(found_count > 0) 
		{
			#print(paste(atac_norm_tmp[ii, 4], '_found the_', found_count, sep=''))
			if(found_count == 1){
				found1_count <- found1_count + 1
				storage_info[k, 1] <- ii
				storage_info[k, 2] <- atac_norm_tmp[ii, 4]
				storage_info[k, 3] <- found_stat[[2]]
				storage_info[k, 4] <- found_stat[[3]]
				storage_info[k, 5] <- 1
			}
			if(found_count == 2){
				found2_count <- found2_count + 1
				#print(paste('ii : ', ii, sep=''))
				#print(paste('k = ', k, sep=''))
				storage_info[k, 1] <- ii
				storage_info[k, 2] <- atac_norm_tmp[ii, 4]
				storage_info[k, 3] <- combine_with_semicolon(found_stat[[2]])
				storage_info[k, 4] <- combine_with_semicolon(found_stat[[3]])
				storage_info[k, 5] <- 2
			}
			if(found_count > 2){
				# #print('There is the situation that found more than 2 !!!')
				# more_than2_tmp <- matrix(c(ii, found_count), 1,2)
				# more_than2_storage <- rbind(more_than2_storage, more_than2_tmp)
			}

			k <- k + 1
		}

	}

	atac_norm_extract <- atac_norm_tmp[storage_info[,1], ]
	atac_norm_extract <- cbind(storage_info[,4], atac_norm_extract)
	atac_norm_extract <- cbind(storage_info[,5], atac_norm_extract)
	colnames(atac_norm_extract)[1] <- 'correspond_gene_number'
	colnames(atac_norm_extract)[2] <- 'gene_name'

	extract_atac[[i]] <- atac_norm_extract
	extract_2more[[i]] <- more_than2_storage

}
names(extract_atac) <- cancer_str
names(extract_2more) <- cancer_str

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac')
save(extract_atac, file='extract_atac_log2norm.RData')
save(extract_2more, file='extract_2more_info.RData')

#-------------------------------------

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts')
file_norm <- list.files()

extract_atac_pro <- list()
cancer_str <- c()
#extract_2more <- list()
for(i in 1:length(file_norm))
{
	setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts')
	atac_norm_tmp <- read.table(file_norm[i], header=T)
	cancer_str <- c(cancer_str, unlist(strsplit(file_norm[i], '_', fixed=T))[[1]])
	print(cancer_str[i])
	

	found1_count <- 0
	#found2_count <- 0
	storage_info <- data.frame()
	#more_than2_storage <- c()
	k <- 1
	#for(ii in 1:1000)
	for(ii in 1:nrow(atac_norm_tmp))
	{
		if(ii %% 3000 == 0) print(ii)
		#print(ii)
		chr_tmp <- as.character(atac_norm_tmp[ii, 1])
		start_tmp <- atac_norm_tmp[ii, 2]
		end_tmp <- atac_norm_tmp[ii, 3]

		found_stat <- judge_in_promotor(chr_tmp, start_tmp, end_tmp, Homo_prot_code)
		found_count <- found_stat[[1]]
		if(found_count > 0) 
		{
			#print(paste(atac_norm_tmp[ii, 4], '_found the_', found_count, sep=''))
			if(found_count == 1){
				found1_count <- found1_count + 1
				storage_info[k, 1] <- ii
				storage_info[k, 2] <- atac_norm_tmp[ii, 4]
				storage_info[k, 3] <- found_stat[[2]]
				storage_info[k, 4] <- found_stat[[3]]
				storage_info[k, 5] <- 1
			}else{
				storage_info[k, 1] <- ii
				storage_info[k, 2] <- atac_norm_tmp[ii, 4]
				storage_info[k, 3] <- combine_with_semicolon_N( found_stat[[2]] )
				storage_info[k, 4] <- combine_with_semicolon_N( found_stat[[3]])
				storage_info[k, 5] <- found_count
			}

			# if(found_count == 2){
			# 	found2_count <- found2_count + 1
			# 	#print(paste('ii : ', ii, sep=''))
			# 	#print(paste('k = ', k, sep=''))
			# 	storage_info[k, 1] <- ii
			# 	storage_info[k, 2] <- atac_norm_tmp[ii, 4]
			# 	storage_info[k, 3] <- combine_with_semicolon(found_stat[[2]])
			# 	storage_info[k, 4] <- combine_with_semicolon(found_stat[[3]])
			# 	storage_info[k, 5] <- 2
			# }
			# if(found_count > 2){
			# 	# #print('There is the situation that found more than 2 !!!')
			# 	# more_than2_tmp <- matrix(c(ii, found_count), 1,2)
			# 	# more_than2_storage <- rbind(more_than2_storage, more_than2_tmp)
			# }

			k <- k + 1
		}

	}

	atac_norm_extract <- atac_norm_tmp[storage_info[,1], ]
	atac_norm_extract <- cbind(storage_info[,4], atac_norm_extract)
	atac_norm_extract <- cbind(storage_info[,5], atac_norm_extract)
	colnames(atac_norm_extract)[1] <- 'correspond_gene_number'
	colnames(atac_norm_extract)[2] <- 'gene_name'

	extract_atac_pro[[i]] <- atac_norm_extract

}
names(extract_atac_pro) <- cancer_str

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac')
save(extract_atac_pro, file='extract_atac_log2norm_promoter.RData')



for(kk in 1:length(extract_atac_pro))
{
	print(dim(extract_atac_pro[[kk]]))
}

#----------------



#-------------function------------------
judge_in_range <- function(chr_tmp, start_tmp, end_tmp, Homo_prot_code)
{
	chr_tmp <- unlist(strsplit(chr_tmp, 'r', fixed=T))[[2]]
	Homo_prot_code_tmp <- Homo_prot_code[which(Homo_prot_code[, 9] == chr_tmp), ]
	found_flag <- 0
	found_loca <- c()
	found_gene <- c()
	for(i in 1:nrow(Homo_prot_code_tmp))
	{
		if(Homo_prot_code_tmp[i, 6]<=start_tmp && end_tmp<=Homo_prot_code_tmp[i, 7])
		{
			found_flag <- found_flag + 1
			found_loca <- c(found_loca, i)
			found_gene <- c(found_gene, Homo_prot_code_tmp[i, 3])
		}

	}

	return(list(found_flag, found_loca, found_gene))
}

judge_in_promotor <- function(chr_tmp, start_tmp, end_tmp, Homo_prot_code)
{
  chr_tmp <- unlist(strsplit(chr_tmp, 'r', fixed=T))[[2]]
  Homo_prot_code_tmp <- Homo_prot_code[which(Homo_prot_code[, 9] == chr_tmp), ]
  
  found_flag <- 0
  found_loca <- c()
  found_gene <- c()
  found_direction <- c()
  for(i in 1:nrow(Homo_prot_code_tmp))
  {
    if(Homo_prot_code_tmp[i, 8] == '+')
    {
      if( (Homo_prot_code_tmp[i, 6]-800) <= start_tmp && end_tmp <= Homo_prot_code_tmp[i, 6])
      {
        found_flag <- found_flag + 1
        found_loca <- c(found_loca, i)
        found_gene <- c(found_gene, Homo_prot_code_tmp[i, 3])
        found_direction <- c(found_direction, '+')
      }
    }else{
      if(Homo_prot_code_tmp[i, 7] <= start_tmp && end_tmp <= (Homo_prot_code_tmp[i, 7]+800) )
      #if(Homo_prot_code_tmp[i, 7] = end_tmp && end_tmp <= (Homo_prot_code_tmp[i, 7]+600) )
      #if(Homo_prot_code_tmp[i, 7] = start_tmp && start_tmp <= (Homo_prot_code_tmp[i, 7]+600) )
      {
        found_flag <- found_flag + 1
        found_loca <- c(found_loca, i)
        found_gene <- c(found_gene, Homo_prot_code_tmp[i, 3])
        found_direction <- c(found_direction, '-')
      }
    }
    
  }
  
  return(list(found_flag, found_loca, found_gene, found_direction))
}

#combine the 2 element into 1 element
combine_with_semicolon <- function(my.obj)
{
	combine_result <- paste0(my.obj[1], sep=';', my.obj[2])

	return(combine_result)

}

#combine the N element into 1 element
combine_with_semicolon_N <- function(my.obj)
{
	len1 <- length(my.obj)
	combine_result <- my.obj[1]
	for(i in 2:len1)
	{
		combine_result <- paste0(combine_result, sep=';', my.obj[i])

	}

	return(combine_result)

}


# setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_08_08_atac/TCGA-ATAC_Cancer_Type-specific_PeakCalls')
# file_peak <- list.files()
# atac_peak_tmp <- read.csv(file_peak[1], header=T, sep='\t')














