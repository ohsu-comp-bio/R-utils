###############################################################################
# This script contains a variety of functions to facilitate microarray, proteomics, ChIP-seq, and RNA-seq analyses
# 
# Authors: Mark Fisher, Jessica Minnier, Theresa Lusardi, Julja Burchard
# April-August, 2016
###############################################################################

require(ggplot2)
require(RColorBrewer)
require(NMF)
require(qvalue)

############################
##Utilities for gs.wrapper##
############################

Remove_gene_lists_of_length_one_from_setlists = function(setlist){
  #If one of the signature sets from among the list that you're passing to gs.wrapper has length <= 1, gs.wrapper will error out. This function filters those sets out.
  #setlist is the list of gene signature lists that you pass to gw.wrapper
  final_list = list()
  for (i in 1:length(setlist)){
    if(length(setlist[[i]])>1){
      final_list[[length(final_list)+1]] = setlist[[i]]
      names(final_list)[length(final_list)] = names(setlist)[i]
    }
  }
  return(final_list)
}

add_rank_ratio_column = function(pathway_annotation_dt){
  tmp_dt = pathway_annotation_dt[order(Condx, Source, E.Value)]
  final_dt = NULL
  for (i in unique(tmp_dt[["Condx"]])){
    for(j in unique(tmp_dt[["Source"]])){
      sub_dt = tmp_dt[tmp_dt[["Source"]]==j & tmp_dt[["Condx"]]==i,]
      if (nrow(sub_dt)>0){
        rank_order = sort(tmp_dt[tmp_dt[["Source"]]==j & tmp_dt[["Condx"]]==i,][["E.Value"]])
        ranks = c(1:length(rank_order))
        contributing_dt = as.data.table(cbind(sub_dt,ranks))
        final_dt = rbind(final_dt, contributing_dt)
      }
    }
  }
  colnames(final_dt) = c(colnames(pathway_annotation_dt), "Rank")
  return (final_dt)
}

add_enrichment_ratio_column = function(pathway_annotation_dt){
  pathway_annotation_dt[,enrichment_ratio := (Overlap/Sig.Size)/(Set.Size/Uni.Size)]
  return(pathway_annotation_dt)
}

###############################################
##Directly relevant to all abundance projects##
###############################################

get_mat_of_greater_than_ave_sd_probesets = function(E_lr_av_all_data_frame){
  #takes a data frame of abundance to av all ratios with proper rownames and returns a subset of rows where the standard dev. in ratio is greater than global standard dev.
  sd_by_row = sapply(1:nrow(E_lr_av_all_data_frame),function(x){sd(E_lr_av_all_data_frame[x,],na.rm=T)})
  sd_by_row = as.matrix(sd_by_row)
  rownames(sd_by_row)=rownames(E_lr_av_all_data_frame)
  global_sd = sd(as.matrix(E_lr_av_all_data_frame))
  sd_greater_log_mat = sd_by_row>global_sd
  E_lr_av_changing_genes = E_lr_av_all_data_frame[sd_greater_log_mat[,1],]
  return(E_lr_av_changing_genes)
}


make_heatmap_versatile = function(matrix_with_samples_for_sig_IDs, sample_string_vec, gsub_remove_str='', annotation_mat=NULL, annotate_with_gene_names = F, ID_colname=NULL, Symbol_colname=NULL, save_image_file=F, string_to_lead_file_name_with = "", c.lim=quantile(matrix_with_samples_for_sig_IDs[,sample_string_vec], probs=.95,na.rm=T)){
  #Makes a heat map using the NMF aheatmap function, using only the columns of a matrix in sample_string_vec
  #matrix_with_samples_for_sig_IDs is a matrix containing all rows with with you want to make a heatmap. It can contain some columns you don't want, because you specify which columns you want in sample_string_vec. Assumes that matrix_with_samples_for_sig_IDs already has rownames of ID (e.g., transcript cluster ID).
  #sample_string_vec is a vector of column names that you want to use to contstruct the heatmap. The most permissive thing you can do is sample_string_vec=colnames(matrix_with_samples_for_sig_IDs).
  #gsub_remove_str: let's say that your column names are very long, and you don't want to label the columns of your heat map with them. You would provide a regular expression pattern in gsub_remove_str, and it will REMOVE things that match that pattern when naming the columns in the heat map.
  #annotation_mat is a matrix containing at the very least a column of IDs (can be longer than rownames(matrix_with_samples_for_sig_IDs)), and a column of something else (e.g., gene symbols) that match those IDs. To be used in case, say, you want to annotate the rows of the heat map with gene symbols instead of transcript cluster IDs.
  #annotate_with_gene_names a boolean deciding whether you want to annotate the heat map using annotation_mat
  #ID_colname is the column name in annotation_mat that contains the IDs of the same type as in rownames(matrix_with_samples_for_sig_IDs), e.g. transcript cluster ID
  #Symbol_colname is the column name in annotation_mat that contains the values (e.g, gene symbols) that the user actually wants to annotate the heat map with
  #save_image_file is a boolean specifying whether 600 dpi png files should be saved or not
  #string_to_lead_file_name_with is a string the represents the first part of the file name to be saved. Can also be used to specify the destination path of the save file if not current working directory.
  #c.lim are the limits in which most of the color variation of the heat map will occur. The min. and max values will also be included in the heat map color range, but everything between them and the c.lim range will be the darkest two colors.
  #Ann_col
  subset_of_mat_of_interest = matrix_with_samples_for_sig_IDs[,sample_string_vec]
  class(subset_of_mat_of_interest) = "numeric"
  #c.lim = quantile(subset_of_mat_of_interest, probs=.95,na.rm=T)
  if(class(subset_of_mat_of_interest) != "matrix"){
    subset_of_mat_of_interest = t(as.matrix(subset_of_mat_of_interest))
  }
  if(ncol(subset_of_mat_of_interest)>2 & nrow(subset_of_mat_of_interest)>1){
    if(annotate_with_gene_names ==F){
      tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_remove_str, '',sample_string_vec),labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
      print(tmp)
      indices_in_heat_map_order = rev(tmp$rowInd)
      ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
      if(save_image_file==T){
        png(filename=paste0(string_to_lead_file_name_with,"_heatmap.png"),width=5,height=5.4,units="in",res=600)
        tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_remove_str, '',sample_string_vec),labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
        dev.off()
      }
      return(list(ordered_probeset_ids))
    } else{
      #then annotate_with_gene_names =T
      ##Find the gene names that correpond to rownames of subset_of_mat_of_interest
      gn_names = vec_of_gene_symbols_given_vec_of_IDs(vec_of_IDs=rownames(subset_of_mat_of_interest), master_mat=annotation_mat, gene_symbol_colname_in_master=Symbol_colname, ID_colname_in_master=ID_colname)
      tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_remove_str, '',sample_string_vec),labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = gn_names)
      print(tmp)
      indices_in_heat_map_order = rev(tmp$rowInd)
      ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
      ordered_gene_ids = vec_of_gene_symbols_given_vec_of_IDs(vec_of_IDs=ordered_probeset_ids, master_mat=annotation_mat, gene_symbol_colname_in_master=Symbol_colname, ID_colname_in_master=ID_colname)
      if(save_image_file==T){
        png(filename=paste0(string_to_lead_file_name_with,"_heatmap.png"),width=5,height=5.4,units="in",res=600)
        tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_remove_str, '',sample_string_vec),labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = gn_names)
        dev.off()
      }
      return(list(ordered_probeset_ids, ordered_gene_ids))
    }
  } else{
    return("")
  }
}

gene_pcaplot <- function(exprdat,sampleid,groupdat=NULL,colorfactor=NULL,shapefactor=NULL,
                         plot_sampleids=TRUE, pcnum=1:2, plottitle = "PCA Plot") {
  #borrowed from Jessica Minnier, 2016
  #adapted from DESeq2:::plotPCA.DESeqTransform
  pca <- prcomp(t(exprdat))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if(is.null(groupdat)) groupdat = data.frame("group"=rep(1,ncol(exprdat)))
  intgroup = colnames(groupdat)
  allgroup <- if (length(intgroup) > 1) {
    factor(apply(groupdat, 1, paste, collapse = ":"))
  }else{allgroup <- intgroup}
  d <- data.frame(PC1 = pca$x[, pcnum[1]], PC2 = pca$x[, pcnum[2]], uniquegroup = allgroup, 
                  groupdat, name = sampleid)
  percentVar <- round(100 * percentVar)
  if(is.null(colorfactor)) {d$color=as.factor(1)}else{
    colnames(d)[colnames(d)==colorfactor] <- "color"}
  if(is.null(shapefactor)) {d$shape=as.factor(1)}else{
    colnames(d)[colnames(d)==shapefactor] <- "shape"
  }
  if(identical(shapefactor,colorfactor)) {d$shape = d$color}
  p <- ggplot(d, aes(PC1, PC2, color=color, shape=shape, size=3)) 
  if(plot_sampleids) {
    p <- p + geom_text(aes(label=name,size=10))
  }else{
    p <- p + geom_point()
  }
  
  if(!is.null(colorfactor)) {
    p <- p + guides(color=guide_legend(title=colorfactor))
  }else {
    p <- p + guides(color = "none")
  }
  if(!is.null(shapefactor)) {
    p <- p + guides(shape=guide_legend(title=shapefactor))
  }else{
    p <- p + guides(shape = "none")
  }
  p <- p + guides(size= "none") + theme_bw() + 
    xlab(paste0("PC",pcnum[1],": ",percentVar[pcnum[1]],"% variance")) +
    ylab(paste0("PC",pcnum[2],": ",percentVar[pcnum[2]],"% variance")) + ggtitle(plottitle)
  
  return(p)
}

#####################################################################
##Utility functions likely also useful for all project applications##
#####################################################################

rename_nonunique_cols = function(data_mat_dt_df){
  #Renames second, third, etc. occurences of a name in colnames(data_mat_dt_df) by concatenating an "_#n" at the end, where #n is the nth occurence of that name detected.
  #data_mat_dt_df is a matrix, data frame, or data table with colnames()
  unique_colnames = unique(colnames(data_mat_dt_df))
  for (cn in unique_colnames){
    cols_with_same_name = which(colnames(data_mat_dt_df)==cn)
    if(length(cols_with_same_name)>1){
      for(i in 2:length(cols_with_same_name)){
        colnames(data_mat_dt_df)[cols_with_same_name[i]] = paste(colnames(data_mat_dt_df)[cols_with_same_name[i]], i, sep="_")
      }
    }
  }
  return(data_mat_dt_df)
}


Generate_mean_by_treatment_matrix = function(Normalized_abundance_list, vec_of_treatments){
  #Function to average abundance values across treatment and return a list of matrices with same rownames and only one col. per unique treatment type.
  #Normalized_abundance_list is a list of abundance matrices, data frames, or data tables (with names() of the list populated). Each matrix/data frame/data table has rownames() and colnames()
  #vec_of_treatments is a vector with a treatment assigned to each sample/column. E.g., if .CEL files were called 001.CEL, 002.CEL, 003.CEL, and 004.CEL, and 001 and 002 are treatment x and 003 and 004 are treatment y, vector would be c("x", "x", "y", "y").
  #Can deal with treatments that only have one column of data/one sample
  #Returns a list of matrices in the same order as the original list. Each colname with be a treatment name from unique(vec_of_treatments) concatenated to, "_Mn"
  final_mat_list =list()
  final_mat = NULL
  for (i in 1:length(Normalized_abundance_list)){
    for (a in 1:length(unique(vec_of_treatments))){
      cols_of_interest = which(vec_of_treatments==unique(vec_of_treatments)[a])
      if (length(cols_of_interest)==1){
        final_mat = cbind(final_mat, as.data.frame(Normalized_abundance_list[[i]])[,cols_of_interest])
      } else{
        final_mat = cbind(final_mat, rowMeans(as.data.frame(Normalized_abundance_list[[i]])[,cols_of_interest]))
      }
      colnames(final_mat)[a] =paste0(unique(vec_of_treatments)[a],"_Mn")
    }
    rownames(final_mat) = rownames(Normalized_abundance_list[[i]])
    final_mat_list[[i]] = final_mat
    names(final_mat_list)[i] = names(Normalized_abundance_list)[i]
    rm(final_mat)
    final_mat=NULL
  }
  return(final_mat_list)
}


Rm_quantile_of_very_different_intensities = function(lin_scale_abundance_mat, quant_val){
  #Removes the quant_val...th quantile of the most differing abundance metric (e.g., probe intensity) across an expression matrix; returns the remainder of the matrix
  #lin_scale_abundance_mat is a matrix of linear scale abundance data
  #quant_val is the cutoff you'd like to use (so, only retain rows where the diff is less than the quant_val the quantile of the data)
  #returns a linear-scale matrix subset of lin_scale_abundance_mat
  log_v_mat = log2(lin_scale_abundance_mat)
  z_tmp = sapply(1:nrow(log_v_mat), FUN=function(x){diff(range(log_v_mat[x,], na.rm=T))})
  log_vec_for_include = z_tmp<quantile(z_tmp,probs=quant_val) #make argument
  final_mat = lin_scale_abundance_mat[log_vec_for_include,]
  rownames(final_mat) = rownames(lin_scale_abundance_mat)[log_vec_for_include]
  return(lin_scale_abundance_mat[log_vec_for_include,])
}


match_up_row_order = function(dt_1, obj_2, is_second_arg_a_list=F){
  #this function re-orders the rownames of obj_2 (see below) to match those of dt_1
  #dt_1 is a data table, data.frame, or matrix (the latter two not yet tested) by which the other data table(s) are to be standardized/rendered in the same order. Must have rownames values.
  #obj_2 can either be a data.table/frame/matrix (the latter two not yet tested) or a list of these
  #is_second_arg_a_list specifices whether obj_2 is a list or not
  if(is_second_arg_a_list==T){
    for(ls_len in 1:length(obj_2)){
      if(nrow(dt_1)!=nrow(obj_2[[ls_len]])){
        print (paste0("Item #",ls_len," does not have the same row number as the first argument"))
        return(NULL)
      }
      order_wrt_dt_1 = match(rownames(dt_1),rownames(obj_2[[ls_len]]))
      obj_2[[ls_len]] = obj_2[[ls_len]][order_wrt_dt_1,]
      rownames(obj_2[[ls_len]]) = rownames(dt_1)
    }
    return(obj_2)
  } else{
    if(nrow(dt_1)!=nrow(obj_2)){
      print (paste0("Item #",ls_len," does not have the same row number as the first argument"))
      return(NULL)
    }
    order_wrt_dt_1 = match(rownames(dt_1),rownames(obj_2))
    obj_2 = obj_2[order_wrt_dt_1,]
    rownames(obj_2) = rownames(dt_1)
    return(obj_2)
  }
}



vec_of_gene_symbols_given_vec_of_IDs = function(vec_of_IDs, master_mat, gene_symbol_colname_in_master, ID_colname_in_master){
  #returns a vector of strings that map to the input vector of strings in a given matrix
  #vec_of_IDs is a vector of strings
  #master_mat is the matrix containing at least two columns, one of which is a superset of vec_of_IDs and the other of which is a superset of the matched strings
  #gene_symbol_colname_in_master is the column name in master_mat corresponding to the target strings that you want as output
  #ID_colname_in_master is the column name in master_mat corresponding to the input vector of strings
  master_mat_match = match(vec_of_IDs,master_mat[,ID_colname_in_master])
  idx1=which(!is.na(master_mat_match))
  master_mat_match=master_mat_match[idx1]
  return(master_mat[master_mat_match, gene_symbol_colname_in_master])
}

########################
##Regression functions##
########################

get_pairwise_log_ratio = function(cntrl, other_treatment_str, E_treatment_wise_intensities){
  cntrl_col = E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""), colnames(E_treatment_wise_intensities))]
  other_treatment_col = E_treatment_wise_intensities[, grep(paste("^",other_treatment_str,"$", sep=""), colnames(E_treatment_wise_intensities))]
  vec_of_lfc = other_treatment_col - cntrl_col
  return(vec_of_lfc)
}

generate_lr_p_beta_q_lfc_aveIntensities = function(cntrl, other_treatment_str, E_treatment_wise_intensities, E, tx, lm.obj, dpath, filename_str){
  all_stats_list = get_lm_obj_stats (lm.obj, cntrl, other_treatment_str, dpath, filename_str)
  ##begin output matrix with pairwise log ratios/log fold change
  output_matrix = get_pairwise_log_ratio(cntrl, other_treatment_str, E_treatment_wise_intensities) 
  output_matrix = as.matrix(output_matrix)
  ##add betas
  output_matrix = cbind(output_matrix, all_stats_list[[1]][,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(all_stats_list[[1]]))])
  ##add p-vals and q-vals
  output_matrix= cbind(output_matrix, all_stats_list[[2]][,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(all_stats_list[[2]]))], all_stats_list[[3]]$qvalues)
  ##add treatmentwise average intensities
  output_matrix= cbind(output_matrix, E_treatment_wise_intensities[,grep(paste("^",other_treatment_str,"$",sep=""),colnames(E_treatment_wise_intensities))], E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""),colnames(E_treatment_wise_intensities))])
  ##add column names
  colnames(output_matrix)=c(paste(other_treatment_str,"lr",cntrl,sep="_"), paste("beta", other_treatment_str,"vs",cntrl, sep="_"), paste("p_val", other_treatment_str,"vs",cntrl, sep="_"), paste("q_val", other_treatment_str,"vs",cntrl, sep="_"), paste("avg_intensity",other_treatment_str,sep="_"), paste("avg_intensity", cntrl, sep="_"))
  rownames(output_matrix) = rownames(E)
  return(output_matrix)
}

generate_lr_p_q_aveIntensities_wilcox = function(cntrl, other_treatment_str, E_treatment_wise_intensities, E, tx, all_stats_list){
  ##begin output matrix with pairwise log ratios/log fold change
  output_matrix = get_pairwise_log_ratio(cntrl, other_treatment_str, E_treatment_wise_intensities) 
  output_matrix = as.matrix(output_matrix)
  
  ##add p-vals and q-vals
  output_matrix= cbind(output_matrix, all_stats_list[[1]], all_stats_list[[2]]$qvalues)
  
  ##add treatmentwise average intensities
  output_matrix= cbind(output_matrix, E_treatment_wise_intensities[,grep(paste("^",other_treatment_str,"$",sep=""),colnames(E_treatment_wise_intensities))], E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""),colnames(E_treatment_wise_intensities))])
  
  ##add column names
  colnames(output_matrix)=c(paste(other_treatment_str,"lr",cntrl,sep="_"), paste("p_val", other_treatment_str,"vs",cntrl, sep="_"), paste("q_val", other_treatment_str,"vs",cntrl, sep="_"), paste("avg_intensity",other_treatment_str,sep="_"), paste("avg_intensity", cntrl, sep="_"))
  rownames(output_matrix) = rownames(E)
  
  return(output_matrix)
}

lm_factor_any_number_of_covariates_dummy_contrasts = function(vec_of_covariate_strings, vec_of_control_strings, phenotype_dt, exprs_matrix){
  #Return an object of type lm() using dummy contrasts
  #Could be made less hacky using factor(tx$Treatment_text, levels=c("Onex","Twox")) instead of adding AA in front of control
  #A generalization of subset_and_lm_by_pair_with_second_treatment_covariate and subset_and_lm_by_pair_with_MDS_covariates with any number of covariates
  #vec_of_covariate_strings is a vector of strings corresponding to colnames(phenotype_dt) of interest to include as covariates
  #vec_of_control_strings is a vector of strings, each element corresponding to the control treatment in the column referred to in the corresponding element in vec_of_covariate_strings
  #phenotype_dt is a data table with the each row corresponding to phenotype of of a sample. Name of phenotype is the colname() of that column. Example below
  # Inv_sample_name Treatment_text Treatment_text_2 Treatment_text_3 Treatment_text_4 Treatment_text_5
  #1                9           Twox              Low          Jan2716            X0077               OS
  #2               10           Onex             High          Oct1215            X1308               OD
  #3               11           Twox              Low           Feb116            X0108               OD
  #4               12           Twox             High          Oct1215            X1308               OS
  #exprs_matrix is any matrix of abundance data with rownames() and colnames()
  
  covariate_list =list()
  for (i in 1:length(vec_of_covariate_strings)){
    col_of_interest = which(colnames(phenotype_dt)==vec_of_covariate_strings[i])
    #assign(paste0("TX_",i), phenotype_dt$vec_of_covariate_strings[i])
    covariate_list[[i]] = phenotype_dt[[col_of_interest]]
    covariate_list[[i]] [covariate_list[[i]]==vec_of_control_strings[i]] = paste0("AA", vec_of_control_strings[i])
    covariate_list[[i]] = as.factor(unlist(covariate_list[[i]], use.names = F))
    names(covariate_list)[i] = paste0("TX_",i)
  }
  master_str = "lm(t(exprs_matrix)~"
  for (i in 1:length(covariate_list)){
    new_str = paste0("unlist(covariate_list[[",i,"]])")
    master_str = paste0(master_str, new_str)
    if (i<length(covariate_list)){ #add a plus sign to everything but the last one
      master_str = paste0(master_str, "+")
    } 
  }
  master_str = paste0(master_str,")")
  return(eval(parse(text=master_str)))
}

get_lm_obj_stats_any_num_covariates_dummy = function (lm.obj, vec_of_other_treatment_strings, vec_of_control_strings, dpath, filename_str){
  #WARNING: can't remember whether this still needs troubleshooting
  #Returns a list with element 1 being a matrix of regression coefficients (betas)
  #Element 2 is a matrix of p-values
  #Element 3 is a list produced by the qvalue() function
  #lm.obj is an object created by an lm() call to a matrix. The colnames() of that matrix were temporarily assigned strings where treatment corresponded to the particular sample (e.g., 1, 2, 3, 4 were reassigned High Low High Low) before being passed to one of the lm
  #vec_of_other_treatment_strings is a vector of characters ...
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  q_val_list = list()
  for(other_treatment_string_cnt in 1: length(vec_of_other_treatment_strings)){
    p_vals_we_care_about = as.vector(p_vals[,grep(paste("^","unlist",vec_of_other_treatment_strings[other_treatment_string_cnt],"$",sep=""),colnames(p_vals))])
    #Report when there was an NA p-value, which has potential to gum up downstream works
    #Need to standardize what to do here
    if(sum(is.na(p_vals_we_care_about))>0){
      rows_with_NAs = which(is.na(p_vals_we_care_about))
      print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
      #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
      #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
      #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
    }
    plot_p_val_hist(p_vals_we_care_about, cntrl_str, other_treatment_str, dpath, filename_str)
    q_val_list[[other_treatment_string_cnt]] = qvalue(p_vals_we_care_about)
    names(q_val_list)[other_treatment_string_cnt] = paste0("qvalue_of_", vec_of_other_treatment_strings[other_treatment_string_cnt])
    #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
    print(plot(q_val_list[[other_treatment_string_cnt]], rng=c(0,0.8)))
    #dev.off()
  }
  return(list(betas, p_vals, q_val_list))
}

get_wilcox_based_qval_vec = function (wilcox_pval_vec, cntrl_str, other_treatment_str, dpath){
  q_vals = qvalue(wilcox_pval_vec);
  #png(filename=paste(dpath,"wilcox_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144);
  print(plot(q_vals));
  #dev.off();
  return(q_vals);
}

plot_p_val_hist = function (pval_vec, cntrl_str, other_treatment_str, dpath, filename_str){
  #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_pDist.png", sep=""),width=5,height=5.4,units="in",res=144)
  print(hist(pval_vec, breaks = round(length(pval_vec)/150), main = paste(filename_str, other_treatment_str, "vs.", cntrl_str, sep=" "), xlab="P-value"))
  #dev.off()
}

wilcox_loop = function (cntrl_str, other_treatment_str, tx, exprs_matrix){
  wilcox_obj_ls = NULL
  for (x in 1:nrow(exprs_matrix)){
    w_obj = wilcox.test(exprs_matrix[x,grep(cntrl_str,colnames(exprs_matrix))],exprs_matrix[x,grep(other_treatment_str,colnames(exprs_matrix))])$p.value
    wilcox_obj_ls = c(wilcox_obj_ls, w_obj)
    names(wilcox_obj_ls)[x] = paste(other_treatment_str, "vs", cntrl_str, rownames(exprs_matrix)[x], sep="_")
  }
  return (wilcox_obj_ls)
}


#######################################################################################
##Directly relevant to all abundance projects, but more useful with some minor tweaks##
#######################################################################################

MA_plot_two_samples_same_treatment = function(Abundance_mat, vec_of_treatment_categories, yrng=NULL , xrngc=range(Abundance_mat, na.rm=T), gsub_str='A_11H1.*') {
  #Function that selects two samples at random from the same treatment in a matrix and generates MA plots with them
  #Likely a good opportunity here to wrap Theresa's MA plot function within this
  #Abundance_mat is a matrix with abundance values. Has rownames() and colnames().
  #vec_of_treatment_categories is a vector of strings of treatment categories. Same length as length(colnames(Abundance_mat)) and corresponds to colnames(Abundance_mat). vec_of_treatment_categories[1] describes the phenotype of colnames(Abundance_mat)[1], etc.
  #yrng can be a vector of two numbers specifying a range for ylim for the MA plot. Otherwise, it will be the 99.9% quantile of the y-values in the log ratios of the two samples of interest.
  #xrngc is a numerical vector of length two specifying xlim for the MA plot.
  #gsub_str: let's say that your column names are very long, and you don't want to label the MA plot with them. You would provide a regular expression pattern in gsub_str, and it will REMOVE things that match that pattern when labeling axes in the MA plot.
  uniq_tmnt = unique(vec_of_treatment_categories)
  for(uniq_cntr in 1:length(uniq_tmnt)){
    cols_of_interest = which(vec_of_treatment_categories %in% uniq_tmnt[uniq_cntr])
    if(length(cols_of_interest)>1){
      rand_two_to_get = sample(cols_of_interest,2)
      xe = rowMeans(Abundance_mat[,rand_two_to_get])
      ye = Abundance_mat[,rand_two_to_get[1]] - Abundance_mat[,rand_two_to_get[2]]
      if(is.null(yrng)){
        yrng = c(-quantile(ye, probs=.999,na.rm=T),quantile(ye, probs=.999,na.rm=T))
      }
      plot(x=xe, y=ye, xlab=paste0('Mn. lg2 ', gsub(gsub_str,'',colnames(Abundance_mat)[rand_two_to_get[1]]), " and ", gsub(gsub_str,'',colnames(Abundance_mat)[rand_two_to_get[2]])), ylab=paste0('Lg ratio ', gsub(gsub_str,'',colnames(Abundance_mat)[rand_two_to_get[1]]), "/", gsub(gsub_str,'',colnames(Abundance_mat)[rand_two_to_get[2]])), cex.lab=1.5, pch='.', col='blue', ylim=yrng, main=uniq_tmnt[uniq_cntr])
      abline(h=0,col=grey(.5),lty="22")
    } else{
      paste0("Not enough replicates for ", uniq_tmnt[uniq_cntr],", silly goose!")
    }
  }
}

make_MA_plot = function(Expression_list, vec_of_treatments, cntrl_str){
  ##Function to plot succession of MA plots with one treatment category vs. another as the "M" axis, with ncol = number of normalization algorithms and nrow=non-control treatments in array set.
  #Inputs are the list of expression matrices (each element represents a particular normalization algorithm), a vector with a treatment assigned to each array. E.g., if .CEL files were called 001.CEL, 002.CEL, 003.CEL, and 004.CEL, and 001 and 002 are treatment x and 003 and 004 are treatment y, vector would be c("x", "x", "y", "y"), and the string of the control treatment's name (e.g., "Mock" or "Control", or perhaps "X" in our example).
  #This could probably be split into two functions-one that generates log ratios to control for all non-control treatments, and one that plots MA plots.
  Treatment_wise_mean_expression_list = Generate_mean_by_treatment_matrix(Expression_list, vec_of_treatments)
  log_ratio_list = list()
  ylims_list = list()
  for(i in 1:length(Treatment_wise_mean_expression_list)){
    cntrl_col = which(colnames(Treatment_wise_mean_expression_list[[i]])==paste0(cntrl_str,"_Mn"))
    log_ratio_mat = NULL
    treatmnt_nm_vec = NULL
    for (j in 1:ncol(Treatment_wise_mean_expression_list[[i]])){
      if (j != cntrl_col){
        log_ratio_mat = cbind(log_ratio_mat, Treatment_wise_mean_expression_list[[i]][,j]-Treatment_wise_mean_expression_list[[i]][,cntrl_col])
        colnames(log_ratio_mat)[ncol(log_ratio_mat)] = paste0(colnames(Treatment_wise_mean_expression_list[[i]])[j], "_vs_", colnames(Treatment_wise_mean_expression_list[[i]])[cntrl_col])
        treatmnt_nm_vec = c(treatmnt_nm_vec,colnames(Treatment_wise_mean_expression_list[[i]])[j])
      }
    }
    rownames(log_ratio_mat) = rownames(Expression_list[[i]])
    log_ratio_list[[i]] = log_ratio_mat
    names(log_ratio_list)[i] = names(Expression_list)[i]
    ylims_list[[i]] = quantile(abs(log_ratio_mat),probs=0.9999, na.rm=T)
    names(ylims_list)[i] = names(Expression_list)[i]
    rm(log_ratio_mat)
  }
  #rs = ncol(Treatment_wise_mean_expression_list[[i]])-1 #number of non-control treatments
  #cs = length(Treatment_wise_mean_expression_list) # number of normalization types
  #par(mfcol = c(rs,cs), cex = 1.2, mar=c(4,4,0.5,0.5), oma=c(1,1,1,1), mgp=c(2,.6,0))
  for(e in 1:length(Treatment_wise_mean_expression_list)){
    for(i in 1:length(treatmnt_nm_vec)){
      xe = rowMeans(cbind(Treatment_wise_mean_expression_list[[e]][,paste0(cntrl_str,"_Mn")], Treatment_wise_mean_expression_list[[e]][,treatmnt_nm_vec[i]]))
      ye = log_ratio_list[[e]][,i]
      plot(x=xe, y=ye, xlab=paste0('Mn lg2 int, cntrl + ', treatmnt_nm_vec[i], ": ", names(Treatment_wise_mean_expression_list)[e]), ylab=colnames(log_ratio_list[[e]])[i], cex.lab=1.5, pch='.', col='blue', ylim=c(-ylims_list[[e]], ylims_list[[e]]))
      abline(h=0,col=grey(.5),lty="22")
    }
  }
  #dev.off()
}

#################################
##Microarray specific functions##
#################################

Noramlization_wrap = function(wd_path_str=params$proj_path, array_type_str=params$array_type){
  setwd(wd_path_str)
  eCELs = list.celfiles('.',full.names=T)
  switch(array_type_str,
         "PrimeView_Human" = {
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Loess##
           normalized_loess_Data_obj_precursor=normalize.loess(exprs(Data_obj))
           colnames(normalized_loess_Data_obj_precursor) = colnames(Expression_rma_normalized_Data_obj)
           loess_Data_obj = Data_obj
           loess_Data_obj = `exprs<-`(loess_Data_obj, normalized_loess_Data_obj_precursor)
           sampleNames(phenoData(loess_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(loess_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           normalized_loess_Data_obj = computeExprSet(loess_Data_obj, pmcorrect.method="pmonly", summary.method="medianpolish")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj - E_av_all_loess
           
           ##Combine##
           Expression_lr_to_av_all_list = list(E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("RMA_quantile", "LOESS")
           Expression_all_list = list(Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("RMA_quantile", "LOESS")
           
           pckg_name="primeView_human"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj, Expression_non_normalized_affy_obj))
         },
         "HG-U133_Plus_2" ={
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##MAS5##
           mas5_normalized_Data_obj=mas5(Data_obj, normalize=T,sc=500, analysis="absolute")
           Expression_mas5_normalized_Data_obj = log2(exprs(mas5_normalized_Data_obj)) 
           E_av_all_mas5 = rowMeans(Expression_mas5_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_mas5 = Expression_mas5_normalized_Data_obj - E_av_all_mas5
           qc_Data_obj = qc(Data_obj)
           percent_p = percent.present(qc_Data_obj)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Loess##
           normalized_loess_Data_obj=expresso(Data_obj,bgcorrect.method="none", normalize.method="loess", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj- E_av_all_loess
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Combine##
           Expression_lr_to_av_all_list = list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           Expression_all_list = list(Expression_mas5_normalized_Data_obj, Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           
           library(hgu133plus2.db)
           pckg_name="hgu133plus2"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj,Expression_non_normalized_affy_obj))
         },
         "Mouse430_2" ={
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##MAS5##
           mas5_normalized_Data_obj=mas5(Data_obj, normalize=T,sc=500, analysis="absolute")
           Expression_mas5_normalized_Data_obj = log2(exprs(mas5_normalized_Data_obj))
           E_av_all_mas5 = rowMeans(Expression_mas5_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_mas5 = Expression_mas5_normalized_Data_obj - E_av_all_mas5
           qc_Data_obj = qc(Data_obj)
           percent_p = percent.present(qc_Data_obj)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Loess##
           normalized_loess_Data_obj=expresso(Data_obj,bgcorrect.method="none", normalize.method="loess", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj- E_av_all_loess
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Combine##
           Expression_all_list = list(Expression_mas5_normalized_Data_obj, Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           Expression_lr_to_av_all_list = list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("MAS 5", "RMA quantile", "LOESS")
           
           library(mouse4302.db)
           pckg_name="mouse4302"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj,Expression_non_normalized_affy_obj))
         },
         "HTA-2_0" ={ 
           Data_obj = read.celfiles(eCELs)
           
           ##MAS5##
           #No mismatch probes for ST arrays, so MAS5 not a thing
           
           ##Unnormalized##
           #summarize with the multimappers
           non_norm_Data_obj_precursor = exprs(Data_obj)
           non_norm_Data_obj_with_multi = Summarize_by_some_custom_ID(non_norm_Data_obj_precursor, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #Remove wildly variable probes and summarize
           non_norm_Data_obj_precursor_no_offenders = Rm_quantile_of_very_different_intensities(non_norm_Data_obj_precursor,0.99999) #linear scale
           probe_level_abundance_list = list(non_norm_Data_obj_precursor_no_offenders)
           names(probe_level_abundance_list)[length(probe_level_abundance_list)] = "Non_norm_no_offenders"
           
           #summarize that
           non_norm_Data_obj_with_multi_no_offenders = Summarize_by_some_custom_ID(non_norm_Data_obj_precursor_no_offenders, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #summarize non_norm_Data_obj_precursor_no_offenders (no offenders!) without the multimappers
           idx_in_f2_with_no_multi = match(non_multiple_mapper_probes,file2_probeset$probe_id)
           idx1=which(!is.na(idx_in_f2_with_no_multi)) #None are na...
           idx_in_f2_with_no_multi=idx_in_f2_with_no_multi[idx1]
           gene_probeset_name_no_multi = file2_probeset$gene_probeset_name[idx_in_f2_with_no_multi]
           non_normalized_Data_obj_no_multi_no_offenders =  Summarize_by_some_custom_ID(non_norm_Data_obj_precursor_no_offenders, non_multiple_mapper_probes, gene_probeset_name_no_multi)
           
           #patch those missing from the latter with the former using rbind; patched will remove  wild probes for now
           missing_multi_only_non_norm = non_norm_Data_obj_with_multi_no_offenders[!(non_norm_Data_obj_with_multi_no_offenders$probeset_id %chin% non_normalized_Data_obj_no_multi_no_offenders$probeset_id),]
           rownames(missing_multi_only_non_norm) = missing_multi_only_non_norm$probeset_id
           Expression_minimally_multi_mapping_no_norm_Data_obj = rbind(non_normalized_Data_obj_no_multi_no_offenders, missing_multi_only_non_norm)
           
           #Give the Expression_minimally_multi_mapping_no_norm_Data_obj rows probeset_id names
           rnm = Expression_minimally_multi_mapping_no_norm_Data_obj[["probeset_id"]]
           Expression_minimally_multi_mapping_no_norm_Data_obj = Expression_minimally_multi_mapping_no_norm_Data_obj[,-ncol(Expression_minimally_multi_mapping_no_norm_Data_obj), with=F]
           rownames(Expression_minimally_multi_mapping_no_norm_Data_obj) = rnm
           Expression_non_normalized_affy_obj = Expression_minimally_multi_mapping_no_norm_Data_obj
           
           #Give the non_norm_Data_obj_with_multi_no_offenders rows probeset_id names and remove the column that currently contains that info.
           rnm2 = non_norm_Data_obj_with_multi_no_offenders[["probeset_id"]]
           non_norm_Data_obj_with_multi_no_offenders = non_norm_Data_obj_with_multi_no_offenders[,-ncol(non_norm_Data_obj_with_multi_no_offenders), with=F]
           rownames(non_norm_Data_obj_with_multi_no_offenders) = rnm2
           
           E_av_all_no_norm = rowMeans(Expression_minimally_multi_mapping_no_norm_Data_obj, na.rm=TRUE)
           E_lr_av_all_no_norm = Expression_minimally_multi_mapping_no_norm_Data_obj- E_av_all_no_norm
           
           ##RMA##
           normalized_rma_Data_obj_precursor=preprocessCore:::normalize.quantiles(exprs(Data_obj))
           rownames(normalized_rma_Data_obj_precursor) = rownames(exprs(Data_obj))
           colnames(normalized_rma_Data_obj_precursor) = colnames(exprs(Data_obj))
           normalized_rma_Data_obj = Summarize_by_some_custom_ID(normalized_rma_Data_obj_precursor, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #Remove wildly variable probes
           normalized_rma_Data_obj_precursor_no_offenders = Rm_quantile_of_very_different_intensities(normalized_rma_Data_obj_precursor,0.99999)
           probe_level_abundance_list[[length(probe_level_abundance_list)+1]] = normalized_rma_Data_obj_precursor_no_offenders
           names(probe_level_abundance_list)[length(probe_level_abundance_list)] = "RMA_probes_no_offenders"
           
           #summarize that
           normalized_rma_Data_obj_no_offenders = Summarize_by_some_custom_ID(normalized_matrix_with_rownames=normalized_rma_Data_obj_precursor_no_offenders, probe_ID_vec=file2_probeset$probe_id, custom_ID_vec=file2_probeset$gene_probeset_name)
           
           #summarize normalized_rma_Data_obj_no_offenders (no offenders!) without the multimappers
           idx_in_f2_with_no_multi = match(non_multiple_mapper_probes, file2_probeset$probe_id)
           idx1=which(!is.na(idx_in_f2_with_no_multi))
           idx_in_f2_with_no_multi=idx_in_f2_with_no_multi[idx1]
           gene_probeset_name_no_multi = file2_probeset$gene_probeset_name[idx_in_f2_with_no_multi]
           normalized_rma_Data_obj_no_offenders_no_multi =  Summarize_by_some_custom_ID(normalized_rma_Data_obj_precursor_no_offenders, non_multiple_mapper_probes, gene_probeset_name_no_multi)
           
           #patch those missing from the latter with the former using rbind; removed wildly variable probes for now
           missing_multi_only_rma_norm = normalized_rma_Data_obj_no_offenders[!(normalized_rma_Data_obj_no_offenders$probeset_id %chin% normalized_rma_Data_obj_no_offenders_no_multi$probeset_id),]
           rownames(missing_multi_only_rma_norm) = missing_multi_only_rma_norm$probeset_id
           Expression_minimally_multi_mapping_rma_Data_obj = rbind(normalized_rma_Data_obj_no_offenders_no_multi, missing_multi_only_rma_norm)
           
           #Give the Expression_minimally_multi_mapping_rma_Data_obj rows probeset_id names  and remove the column that currently contains that info.
           rnm = Expression_minimally_multi_mapping_rma_Data_obj[["probeset_id"]]
           Expression_minimally_multi_mapping_rma_Data_obj = Expression_minimally_multi_mapping_rma_Data_obj[,-ncol(Expression_minimally_multi_mapping_rma_Data_obj), with=F]
           rownames(Expression_minimally_multi_mapping_rma_Data_obj) = rnm
           
           #Give the normalized_rma_Data_obj_no_offenders rows probeset_id names and remove the column that currently contains that info.
           rnm2 = normalized_rma_Data_obj_no_offenders[["probeset_id"]]
           normalized_rma_Data_obj_no_offenders = normalized_rma_Data_obj_no_offenders[,-ncol(normalized_rma_Data_obj_no_offenders), with=F]
           rownames(normalized_rma_Data_obj_no_offenders) = rnm2
           
           E_av_all_rma = rowMeans(Expression_minimally_multi_mapping_rma_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = as.data.frame(Expression_minimally_multi_mapping_rma_Data_obj)- E_av_all_no_norm
           rownames(E_lr_av_all_rma) = rownames(Expression_minimally_multi_mapping_rma_Data_obj)
           #Expression_rma_normalized_Data_obj = Expression_minimally_multi_mapping_rma_Data_obj
           
           ##Loess##
           #summarize with the multimappers but without offenders
           normalized_loess_Data_obj_precursor=normalize.loess(exprs(Data_obj))
           normalized_loess_Data_obj_multi =  Summarize_by_some_custom_ID(normalized_loess_Data_obj_precursor, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #try normalizing after you take the log
           normalized_loess_Data_obj_precursor_after_log = normalize.loess(log2(exprs(Data_obj)))
           normalized_loess_Data_obj_multi_after_log = Summarize_by_some_custom_ID(normalized_loess_Data_obj_precursor_after_log, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #Remove wildly variable probes
           normalized_loess_Data_obj_precursor_no_offenders = Rm_quantile_of_very_different_intensities(normalized_loess_Data_obj_precursor, 0.99999)
           probe_level_abundance_list[[length(probe_level_abundance_list)+1]] = normalized_loess_Data_obj_precursor_no_offenders
           names(probe_level_abundance_list)[length(probe_level_abundance_list)] = "Loess_probes_no_offenders"
           
           #summarize that
           normalized_loess_Data_obj_multi_no_offenders =  Summarize_by_some_custom_ID(normalized_loess_Data_obj_precursor_no_offenders, file2_probeset$probe_id, file2_probeset$gene_probeset_name)
           
           #summarize without the multimappers and without offenders
           idx_in_f2_with_no_multi = match(non_multiple_mapper_probes,file2_probeset$probe_id)
           idx1=which(!is.na(idx_in_f2_with_no_multi))
           idx_in_f2_with_no_multi=idx_in_f2_with_no_multi[idx1]
           gene_probeset_name_no_multi = file2_probeset$gene_probeset_name[idx_in_f2_with_no_multi]
           normalized_loess_Data_obj_no_multi =  Summarize_by_some_custom_ID(normalized_loess_Data_obj_precursor_no_offenders, non_multiple_mapper_probes, gene_probeset_name_no_multi)
           #num_probesets_using_ea_probe=ave(gene_probeset_name_no_multi,non_multiple_mapper_probes, FUN=function(x){length(unique(x[!is.na(x)]))})
           num_probes_using_ea_probeset=ave(non_multiple_mapper_probes,gene_probeset_name_no_multi, FUN=function(x){length(unique(x[!is.na(x)]))})
           sorted_match_counts<-table(num_probes_using_ea_probeset)[order(table(num_probes_using_ea_probeset),decreasing=T)]#table basically assigns them to bins with different counts of each occurrence
           
           
           #patch those missing from the multimapper removal with the those including multi-mappers but missing wildly variable probes using rbind
           missing_multi_only=normalized_loess_Data_obj_multi_no_offenders[!(normalized_loess_Data_obj_multi_no_offenders$probeset_id %chin% normalized_loess_Data_obj_no_multi$probeset_id),]
           rownames(missing_multi_only) = missing_multi_only$probeset_id
           Expression_minimally_multi_mapping_normalized_loess_Data_obj = rbind(normalized_loess_Data_obj_no_multi, missing_multi_only)
           
           #Give the Expression_minimally_multi_mapping_normalized_loess_Data_obj rows probeset_id names  and remove the column that currently contains that info.
           rnm = Expression_minimally_multi_mapping_normalized_loess_Data_obj[["probeset_id"]]
           Expression_minimally_multi_mapping_normalized_loess_Data_obj = Expression_minimally_multi_mapping_normalized_loess_Data_obj[,-ncol(Expression_minimally_multi_mapping_normalized_loess_Data_obj), with=F]
           rownames(Expression_minimally_multi_mapping_normalized_loess_Data_obj) = rnm
           
           #Give the normalized_loess_Data_obj_multi_no_offenders rows probeset_id names  and remove the column that currently contains that info.
           rnm2 = normalized_loess_Data_obj_multi_no_offenders[["probeset_id"]]
           normalized_loess_Data_obj_multi_no_offenders = normalized_loess_Data_obj_multi_no_offenders[,-ncol(normalized_loess_Data_obj_multi_no_offenders), with=F]
           rownames(normalized_loess_Data_obj_multi_no_offenders) = rnm2
           
           E_av_all_loess = rowMeans(Expression_minimally_multi_mapping_normalized_loess_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = as.data.frame(Expression_minimally_multi_mapping_normalized_loess_Data_obj)- E_av_all_loess
           rownames(E_lr_av_all_loess) = rownames(Expression_minimally_multi_mapping_normalized_loess_Data_obj)
           
           ##Combine##
           Expression_all_list = list(normalized_rma_Data_obj_no_offenders, Expression_minimally_multi_mapping_rma_Data_obj, normalized_loess_Data_obj_multi_no_offenders, Expression_minimally_multi_mapping_normalized_loess_Data_obj) #non_norm_Data_obj_with_multi,Expression_minimally_multi_mapping_no_norm_Data_obj, #Expression_rma_normalized_Data_obj
           names(Expression_all_list) = c("RMA_all_multi", "RMA_min._multi","Loess_all_multi", "Loess_min._multi") #"No norm all multi", "No norm min. multi",
           Expression_all_list = match_up_row_order(Expression_minimally_multi_mapping_no_norm_Data_obj, Expression_all_list, is_second_arg_a_list=T)
           Expression_lr_to_av_all_list = list(E_lr_av_all_rma, E_lr_av_all_loess) #E_lr_av_all_no_norm,  #list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("RMA_min._multi", "Loess_min._multi") #"No norm min multi",  #c("MAS 5", "RMA quantile", "LOESS")
           
           #library(hta20sttranscriptcluster.db)
           pckg_name="hta20sttranscriptcluster"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj,Expression_non_normalized_affy_obj, probe_level_abundance_list))
         },
         stop("Enter something that switches me!")
  )
}


Remove_mulimappers_and_return_probe_IDs = function(vec_of_probe_IDs, vec_of_probeset_IDs){
  #Takes as arguments a vector of probesetID and a vector of probes. Asks how many probesets are using each of the probes in that vector. Returns only those probes that are not used by more than one probeset.
  #Assumes that vec_of_probe_IDs and vec_of_probeset_IDs are ordered the same way and correspond.
  num_probesets_using_ea_probe=ave(vec_of_probeset_IDs,vec_of_probe_IDs, FUN=function(x){length(unique(x[!is.na(x)]))})
  sorted_match_counts<-table(num_probesets_using_ea_probe)[order(table(num_probesets_using_ea_probe),decreasing=T)]#table basically assigns them to bins with different counts of each occurrence
  non_multiple_mappers=vec_of_probe_IDs[!(num_probesets_using_ea_probe>1)]
  return(non_multiple_mappers)
}

Summarize_by_some_custom_ID = function(normalized_matrix_with_rownames, probe_ID_vec, custom_ID_vec){
  #WARNING: this function has not been convincingly troubleshooted and may contain at leas two potential issues: 1) the match statements may be reversed, although I don't think so and 2) hasn't been troubleshooted yet for what happens when nrow(normalized_matrix_with_rownames) is smaller than length(custom_ID_vec)
  #summarizes a normalized matrix (that is in linear space). Maps the probe IDs to the matrix. Assumes that the probe_ID_vec is ALREADY MAPPED CORRECTLY to the custom_ID_vec (i.e., they should be the same length and correspond to one another).
  #normalized_matrix_with_rownames is a matrix of normalized data on a linear scale.
  #probe_ID_vec is a vector of strings of IDs of the same type/ilk as rownames(normalized_matrix_with_rownames)
  #custom_ID_vec is a vector of strings of IDs to which those probes/IDs from probe_ID_vec are to be collapsed/summarized
  idx2=match(probe_ID_vec,rownames(normalized_matrix_with_rownames)) #position in y where x is
  idx1=which(!is.na(idx2))
  idx2=idx2[idx1]
  ready_for_summarization=matrix(NA, ncol=ncol(normalized_matrix_with_rownames), nrow=length(custom_ID_vec))  #nrow=min(nrow(normalized_matrix_with_rownames), length(custom_ID_vec))) #probably need to test this when ncol(matrix) is smaller than length(custom_ID_vec)
  ready_for_summarization[idx1,]=normalized_matrix_with_rownames[idx2,] 
  rownames(ready_for_summarization)=custom_ID_vec #[idx1] #should I index by idx1 here? #there will be NAs in here if normalized_matrix_with_rownames doesn't contain all probes in the array 
  colnames(ready_for_summarization) = colnames(normalized_matrix_with_rownames)
  sData = oligo::summarize(ready_for_summarization, method="medianpolish")
  out_dt = as.data.table(sData)
  out_dt$probeset_id = rownames(sData)
  rownames(out_dt) = rownames(sData)
  return(out_dt)
}

merge_annotation_dt_with_big_stats_mat = function(big_stats_mat, big_stats_mat_key, annotation_dt, annotation_dt_key){
  ##assumes that the annotation_dt has two columns: probeset and gene symbol and that big_stat_mat has a probeset ID column. The string representing the probeset colname in each case gets passed to the function
  big_stats_dt = as.data.table(big_stats_mat)
  setkeyv(big_stats_dt, c(big_stats_mat_key))
  setkeyv(annotation_dt, c(annotation_dt_key))
  return(annotation_dt[big_stats_dt])
}

merge_affy_pck_annotation_with_big_stats_mat = function (big_stats_mat, pckg_name){
  colnum_of_interest = which(colnames(big_stats_mat)=="Probeset_IDs")
  pckg_probesetIDs_df = as.data.frame(big_stats_mat[,colnum_of_interest])
  colnames(pckg_probesetIDs_df)=c(paste(pckg_name,"probesetIDs", sep="_"))
  
  x <- get(paste(pckg_name,"ENTREZID", sep="")) # Get the probe identifiers that are mapped to an ENTREZ Gene ID 
  mapped_probes <- mappedkeys(x) # Convert to a list 
  pckg_probesetID_to_pckg_entrez_ls <- as.list(x[mapped_probes])
  pckg_probesetIDs = names(pckg_probesetID_to_pckg_entrez_ls)
  pckg_entrez_IDs = unlist(pckg_probesetID_to_pckg_entrez_ls, use.names=F)
  pckg_probesetID_to_pckg_entrez_df=as.data.frame(cbind(pckg_probesetIDs, pckg_entrez_IDs))
  colnames(pckg_probesetID_to_pckg_entrez_df) = c(paste(pckg_name,"probesetIDs", sep="_"), paste(pckg_name,"entrez_IDs", sep="_"))
  
  pckg_merged = merge(pckg_probesetIDs_df, pckg_probesetID_to_pckg_entrez_df, all.x=TRUE, by.x = paste(pckg_name,"probesetIDs",sep="_"), by.y = paste(pckg_name,"probesetIDs",sep="_"))
  
  x <- get(paste(pckg_name,"SYMBOL", sep="")) # Get the probe identifiers that are mapped to an ENTREZ Gene ID 
  mapped_probes <- mappedkeys(x) # Convert to a list 
  pckg_probesetID_to_pckg_symbol_ls <- as.list(x[mapped_probes])
  pckg_probesetIDs = names(pckg_probesetID_to_pckg_symbol_ls)
  pckg_symbol_IDs = unlist(pckg_probesetID_to_pckg_symbol_ls, use.names=F)
  pckg_probesetID_to_pckg_symbol_df=as.data.frame(cbind(pckg_probesetIDs, pckg_symbol_IDs))
  colnames(pckg_probesetID_to_pckg_symbol_df) = c(paste(pckg_name,"probesetIDs", sep="_"), paste(pckg_name, "symbol_IDs", sep="_"))
  
  pckg_merged = merge(pckg_merged, pckg_probesetID_to_pckg_symbol_df, all.x=TRUE, by.x =paste(pckg_name,"probesetIDs",sep="_"), by.y = paste(pckg_name,"probesetIDs",sep="_"))
  big_stats_mat = rename_nonunique_cols(big_stats_mat)
  pckg_merged = merge(pckg_merged, big_stats_mat, all=T, by.x=paste(pckg_name,"probesetIDs",sep="_"), by.y="Probeset_IDs")
  return(pckg_merged)
}

map_mmu_entrez_to_mogene20st_probesetIDs = function(vec_of_mmu_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst){
  if (length(vec_of_mmu_entrez_IDs) != length(unique(vec_of_mmu_entrez_IDs))){
    print ("your mouse entrez id vector contains non-unique items")
  }
  print(vec_of_mmu_entrez_IDs[1])
  vec_of_mogene20st_probesetIDs=c()
  for (mmus_eID in vec_of_mmu_entrez_IDs){
    vec_of_mogene20st_probesetIDs=c(vec_of_mogene20st_probesetIDs, names(which(mogene20st_trans_clust_ENTREZID_lst==mmus_eID)))	
  }
  return (vec_of_mogene20st_probesetIDs)
}

map_hsa_entrez_to_mogene20st_probesetIDs = function(vec_of_hsa_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst){
  rm(vec_of_mogene20st_probesetIDs)
  ##translate hsa entrez IDs to mmu entrez IDs
  ls_of_mmu_entrez_IDs = inpIDMapper(vec_of_hsa_entrez_IDs, srcSpecies="HOMSA", destSpecies="MUSMU",srcIDType="EG",destIDType="EG", keepMultDestIDMatches = FALSE)
  ##this is a list where names() corresponds to human entrez IDs and the items in the list correspond to mouse entrez IDs. One problem will be if multiple mouse entrez IDs map to the same human entrez ID, in which case the length of that element in the list will be >1. I don't think that these are ones we want, but let's check whether we even have them##
  xs=lapply(ls_of_mmu_entrez_IDs, function(e) length(e)>1)
  xs=unlist(xs)
  print(paste("there are ",length(names(ls_of_mmu_entrez_IDs)[xs]), " mouse entrez IDs that have more than one human entrez ID mapped to them", sep=""))
  ##for now, subset the mouse entrez ids by those not included in this list##
  xsub=lapply(ls_of_mmu_entrez_IDs, function(e) e[1])
  xsub=unlist(xsub)
  vec_of_mmu_entrez_IDs=xsub[!xs]
  return (map_mmu_entrez_to_mogene20st_probesetIDs(vec_of_mmu_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst))
}


############################################
##Have been superseded by better functions##
############################################

get_lm_obj_stats = function (lm.obj, cntrl_str, other_treatment_str, dpath, filename_str){
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  p_vals_we_care_about = as.vector(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))])
  if(sum(is.na(p_vals_we_care_about))>0){
    rows_with_NAs = which(is.na(p_vals_we_care_about))
    print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
    #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
    #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
    #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
  }
  plot_p_val_hist(p_vals_we_care_about, cntrl_str, other_treatment_str, dpath, filename_str)
  q_vals = qvalue(p_vals_we_care_about)
  #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
  print(plot(q_vals, rng=c(0,0.8)))
  #dev.off()
  return(list(betas, p_vals, q_vals))
}

get_lm_obj_stats_two_factors = function (lm.obj, dpath){
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  p_vals_we_care_about = p_vals[,grep(paste("^","Tx.*","$",sep=""),colnames(p_vals))]
  class(p_vals_we_care_about) = "matrix"
  q_vals = list()
  for (i in 1:ncol(p_vals_we_care_about)){
    if(sum(is.na(p_vals_we_care_about[,i]))>0){
      rows_with_NAs = which(is.na(p_vals_we_care_about[,i]))
      print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
      #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
      #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
      #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
    }
    plot_p_val_hist_two_factor(p_vals_we_care_about[,i], dpath, colnames(p_vals_we_care_about)[i])
    q_vals[[i]] = qvalue(p_vals_we_care_about[,i])
    names(q_vals)[i] = colnames(p_vals_we_care_about)[i]
    #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
    print(plot(q_vals[[i]], rng=c(0,0.8)))
    #dev.off()
  }
  return(list(betas, p_vals, q_vals))
}


subset_and_lm_by_pair_with_percentP_covariate = function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix, percentP_vec){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("AA",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx+percentP_vec))
}


subset_and_lm_by_pair_with_second_treatment_covariate = function (cntrl_treatment_str_tr1, cntrl_treatment_str_tr2, tx, exprs_matrix){
  Tx_1 = tx$Treatment_text
  Tx_1[Tx_1==cntrl_treatment_str_tr1] = paste("AA",cntrl_treatment_str_tr1, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx_1 = as.factor(Tx_1)
  Tx_2 = tx$Treatment_text_2
  Tx_2[Tx_2==cntrl_treatment_str_tr2] = paste("AA",cntrl_treatment_str_tr2, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx_2 = as.factor(Tx_2)
  #print(Tx)
  return(lm(t(exprs_matrix)~Tx_1+Tx_2))
}

subset_and_lm_by_pair_with_MDS_covariates= function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix, MDS_x, MDS_y){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("aa",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx+MDS_x+MDS_y))
}

plot_p_val_hist_two_factor = function (pval_vec, dpath, filename_str){
  print(hist(pval_vec, breaks = round(length(pval_vec)/150), main = filename_str, xlab="P-value"))
}

subset_and_lm_by_pair= function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("aa",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx))
}


merge_affy_provided_annotation_with_big_stats_mat = function (big_stats_mat, annotation_csv_path){
  annotation_dt = fread(annotation_csv_path)
  big_stats_dt = as.data.table(big_stats_mat) #retains colnames, loses rownames
  setkey(annotation_dt, "Probe Set ID")
  setkey(big_stats_dt, Probeset_IDs)
  final_dt = annotation_dt[big_stats_dt]
  return(final_dt)
}


make_MA_plot_rows_cols = function(Expression_list, vec_of_treatments, cntrl_str){
  #Not really useful for knitr reports, but there's probably some salvageable code in here
  #A version of make_MA_plot that makes a matrix of plots, with nrow= #non-control treatments and ncol = #normalized matrices passed to the function
  ##Function to plot a panel of MA plots, with ncol = number of normalization algorithms and nrow=non-control treatments in array set.
  #Inputs are the list of expression matrices (each element represents a particular normalization algorithm), a vector with a treatment assigned to each array. E.g., if .CEL files were called 001.CEL, 002.CEL, 003.CEL, and 004.CEL, and 001 and 002 are treatment x and 003 and 004 are treatment y, vector would be c("x", "x", "y", "y"), and the string of the control treatment's name (e.g., "Mock" or "Control", or perhaps "X" in our example).
  #This could probably be split into two functions-one that generates log ratios to control for all non-control treatments, and one that plots MA plots.
  Treatment_wise_mean_expression_list = Generate_mean_by_treatment_matrix(Expression_list, vec_of_treatments)
  log_ratio_list = list()
  for(i in 1:length(Treatment_wise_mean_expression_list)){
    cntrl_col = which(colnames(Treatment_wise_mean_expression_list[[i]])==paste0(cntrl_str,"_Mn"))
    log_ratio_mat = NULL
    treatmnt_nm_vec = NULL
    for (j in 1:ncol(Treatment_wise_mean_expression_list[[i]])){
      if (j != cntrl_col){
        log_ratio_mat = cbind(log_ratio_mat, Treatment_wise_mean_expression_list[[i]][,j]-Treatment_wise_mean_expression_list[[i]][,cntrl_col])
        colnames(log_ratio_mat)[ncol(log_ratio_mat)] = paste0(colnames(Treatment_wise_mean_expression_list[[i]])[j], "_vs_", colnames(Treatment_wise_mean_expression_list[[i]])[cntrl_col])
        treatmnt_nm_vec = c(treatmnt_nm_vec,colnames(Treatment_wise_mean_expression_list[[i]])[j])
      }
    }
    rownames(log_ratio_mat) = rownames(Expression_list[[i]])
    log_ratio_list[[i]] = log_ratio_mat
    names(log_ratio_list)[i] = names(Expression_list)[i]
    rm(log_ratio_mat)
  }
  rs = ncol(Treatment_wise_mean_expression_list[[i]])-1 #number of non-control treatments
  cs = length(Treatment_wise_mean_expression_list) # number of normalization types
  par(mfcol = c(rs,cs), cex = 1.2, mar=c(4,4,0.5,0.5), oma=c(1,1,1,1), mgp=c(2,.6,0))
  for(e in 1:length(Treatment_wise_mean_expression_list)){
    for(i in 1:length(treatmnt_nm_vec)){
      xe = rowMeans(cbind(Treatment_wise_mean_expression_list[[e]][,paste0(cntrl_str,"_Mn")], Treatment_wise_mean_expression_list[[e]][,treatmnt_nm_vec[i]]))
      ye = log_ratio_list[[e]][,i]
      plot(x=xe, y=ye, xlab=paste0('Mn lg2 int, cntrl + ', treatmnt_nm_vec[i], ": ", names(Treatment_wise_mean_expression_list)[e]), ylab=colnames(log_ratio_list[[e]])[i], cex.lab=1.2, pch='.', col='blue')
      abline(h=0,col=grey(.5),lty="22")
    }
  }
  #dev.off()
}

make_heatmaps_exacloud = function(matrix_with_ID_symbol_and_betas_for_sig_IDs, beta_string_vec, ID_colname, Symbol_colname, gsub_str=''){
  #Assumes that matrix_with_ID_symbol_and_betas_for_sig_IDs already has rownames of ID (e.g., transcript cluster ID)
  subset_of_mat_of_interest = matrix_with_ID_symbol_and_betas_for_sig_IDs[,beta_string_vec]
  class(subset_of_mat_of_interest) = "numeric"
  c.lim = quantile(subset_of_mat_of_interest, probs=.95,na.rm=T)
  if(class(subset_of_mat_of_interest) != "matrix"){
    subset_of_mat_of_interest = t(as.matrix(subset_of_mat_of_interest))
  }
  if(ncol(subset_of_mat_of_interest)>2 & nrow(subset_of_mat_of_interest)>1){
    png(filename=paste0(params$dpath,"Heatmap_",Sys.Date(),".png"),width=5,height=5.4,units="in",res=144)
    print(tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_str, '',beta_string_vec),labCol=gsub(gsub_str, '',beta_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest))))
    #print(tmp)
    dev.off()
    indices_in_heat_map_order = rev(tmp$rowInd)
    ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
    indices_in_original_matrix = matrix_with_ID_symbol_and_betas_for_sig_IDs[,ID_colname] %in% ordered_probeset_ids
    ordered_gene_symbols = matrix_with_ID_symbol_and_betas_for_sig_IDs[indices_in_original_matrix,Symbol_colname]
    return(list(ordered_probeset_ids, ordered_gene_symbols))
  }else{
    return("")
  }
}


make_heatmaps = function(matrix_with_ID_symbol_and_betas_for_sig_IDs, beta_string_vec, ID_colname, Symbol_colname, gsub_str=''){
  #Assumes that matrix_with_ID_symbol_and_betas_for_sig_IDs already has rownames of ID (e.g., transcript cluster ID)
  subset_of_mat_of_interest = matrix_with_ID_symbol_and_betas_for_sig_IDs[,beta_string_vec]
  class(subset_of_mat_of_interest) = "numeric"
  c.lim = quantile(subset_of_mat_of_interest, probs=.95,na.rm=T)
  if(class(subset_of_mat_of_interest) != "matrix"){
    subset_of_mat_of_interest = t(as.matrix(subset_of_mat_of_interest))
  }
  if(ncol(subset_of_mat_of_interest)>2 & nrow(subset_of_mat_of_interest)>1){
    tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_str, '',beta_string_vec),labCol=gsub(gsub_str, '',beta_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
    print(tmp)
    indices_in_heat_map_order = rev(tmp$rowInd)
    ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
    indices_in_original_matrix = matrix_with_ID_symbol_and_betas_for_sig_IDs[,ID_colname] %in% ordered_probeset_ids
    ordered_gene_symbols = matrix_with_ID_symbol_and_betas_for_sig_IDs[indices_in_original_matrix,Symbol_colname]
    return(list(ordered_probeset_ids, ordered_gene_symbols))
  }else{
    return("")
  }
}

###########################
##Needs major improvement##
###########################

panel.smoothcust = function (x, y, col = par("col"), bg = NA, pch = ".", 
		cex = 1, col.smooth = "blue", span = 2/3, iter = 3, ...) 
{
	points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	ok <- is.finite(x) & is.finite(y)
	if (any(ok)) 
		#lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
		abline(0,1)
}


chart.Corr = function (R, histogram = FALSE, method = c("pearson", "kendall","spearman"), ...) 
{
	x = checkData(R, method = "matrix")
	if (missing(method)) 
		method = method[1]
	panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
			method, cex.cor, ...) {
		usr <- par("usr")
		on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- cor(x, y, use = use, method = method)
		txt <- format(c(r, 0.123456789), digits = digits)[1]
		txt <- paste(prefix, txt, sep = "")
		#print(txt)
		strwidth(txt)
		if (missing(cex.cor)) 
			cex <- 1.5 #/strwidth(txt)
		test <- cor.test(x, y, method = method)
		Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
				cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
						"**", "*", ".", " "))
		#text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
		text(0.5, 0.5, txt, cex = cex)
		text(0.7, 0.7, Signif, cex = cex-0.25, col = 1) #col = 2
	}
	f <- function(t) {
		dnorm(t, mean = mean(x), sd = sd.xts(x))
	}
	hist.panel = function(x, ...) {
		par(new = TRUE)
		hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
				main = "", breaks = "FD", pch=".")
		lines(density(x, na.rm = TRUE), col = "blue", lwd = 1)
		rug(x)
	}
	#print(x)
	#print(class(x))
	if (histogram) 
		pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
				diag.panel = hist.panel, method = method, ...) #(pch=".")
	else pairs(x, gap = 0, lower.panel = panel.smoothcust, upper.panel = panel.cor, 
				method = method, ...)
}

