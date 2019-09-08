
#----------load R
# qsub -I -q interactive -l nodes=1:ppn=6,vmem=128gb,walltime=03:00:00 -m ae 
# module load r/3.6.0


setwd('/N/dc2/projects/zhangclab/wennan/TCGA_Methydata')
file_here <- list.files(pattern='450')

i = 1

load(file_here[1])
cancer_str <- unlist(strsplit(file_here[i], split='.',fixed=T))[[1]]


setwd("/N/dc2/projects/zhangclab/wennan/Methy_ictd")
#load("/N/dc2/projects/zhangclab/wennan/TCGA_Methydata/BRCA.methy450.RData")

#extract gene name
rownames(datat) <- gsub('^.*\\|', '', rownames(datat))

###############################
#loading ICTD
############################
mainDir <- '/N/dc2/projects/zhangclab/wennan/ICTD_source_file'
setwd(mainDir)
#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header.r")
source("./Loading_header_new_IM.r")
#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header_new_Brain.r")
source("./All_new_functions.r")
source("./Step_3_linking_graph_function.r")
source("./R2_all_function.R")
source("./temp_func.r")
source("./temp_func3.r")
source("./temp_func3_new.r")
source("./01-23NMF/NMF_functions_new.r")
source("./Step2plus_Celltype_marker_Hierarchical_CTES.r")
source("./Step2plus_Celltype_marker_Hierarchical_CTES2.r")

library(NMF)
set.seed(123456)
source("./02-01-NMF_function/nmf.library.R")
source("./02-01-NMF_function/ini.R")
source("./02-01-NMF_function/Step_3_linking_graph_function.r")
source("./02-01-NMF_function/R4RR_selection_functions.r")
source("./02-01-NMF_function/NMF_functions_new.r")
source("./02-01-NMF_function/NMF_method1_test_version.r")
source("./02-01-NMF_function/NMF_method1_test_version2.r")
source("./02-01-NMF_function/NMF_method1_test_version3.r")

source("./debug_20190217.R")  #debug "compute_CompRowspace_NN_uni" function



#-------------run ictd
source("/N/dc2/projects/zhangclab/wennan/ATAC_ictd/ICTD_ATAC_seq_functions.r")


	data_c3<-datat[which(rownames(datat)%in%rownames(immune_cell_uni_table0_GS)==1),]
	data_c3 <- na.omit(data_c3)

	data_t <- data_c3
    ccc<-apply(data_t,1,var)
	cutc<-quantile(ccc,0.3)
	tg_ids<-which(ccc>cutc)
	data_t0<-data_t[tg_ids,]
	tg_dgenes<-names(sort(table(rownames(data_t0)),decreasing=T)[1:200])
	tg_d_ids<-c()
	for(i in 1:length(tg_dgenes))
	{
	tg_d_ids<-c(tg_d_ids,which(rownames(data_t0)==tg_dgenes[i]))
	}
	data_t0<-data_t0[-tg_d_ids,]
	#use data_t0​


	data01<-normalize_data2(data_t0)
    data_CORS_cancer<-data01
    data_ccc<-data_CORS_cancer

    nn<-rownames(data_ccc)
    nn<-paste(nn,1:length(nn),sep="_")
    rownames(data_ccc)<-nn

    tg_R<-paste(cancer_str,"_Methy_R1.RData",sep="")
    list_c1<-MRHCA_IM_compute_MR2(data_CORS_cancer=data_ccc,IM_id_list)
    MR_IM_result_new_c<-MRHCA_IM_compute_full_pub_new(data_CORS_cancer=data_ccc,list_c=list_c1,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
    list_new_c2<-Process_MR_IM_result_new2(MR_IM_result_new_c,tg_key_c=tg_key0,cell_type_enrich_cut=0.4,cor_cut0=0.75,num_cut=6,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
    R1_filter_step1_results_new<-R1_list_filtering_step1_new2(list_new_c2,data_CORS_cancer=data_ccc,max_cut=10,cutn0=7,cut10=0.75,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)#cut10=0.8 for RNA-seq, cut10=0.75 for microarray, cut10=0.85 for scRNA-seq simulated data
    tg_R1_lists<-R1_filter_step1_results_new[[4]]

	print(length(tg_R1_lists))

    setwd('/N/dc2/projects/zhangclab/wennan/Methy_ictd/R1_result')
    getwd()
    #save tmp
    #save(data_t, file="BRCA_Methy_tumor.RData")
    save(data_CORS_cancer,list_new_c2,R1_filter_step1_results_new,tg_R1_lists,file=tg_R)




