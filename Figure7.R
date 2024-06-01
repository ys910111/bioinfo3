library(data.table)
library(stringr)
library(dplyr)
library(NMF)

samples_105 = readRDS('../list_suf_samples.rds')

# Figure 7A ####
cna = data.frame(fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct'), check.names = F, row.names = 1)
cna = cna[, samples_105]
mrna = data.frame(fread('../PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct'), check.names = F, row.names = 1)
mrna = mrna[, samples_105]
prot = data.frame(fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct'), check.names = F, row.names = 1)
prot = prot[, samples_105]
glyco = data.frame(fread('../PDAC_LinkedOmics_Data/N-glycoproteomics_Site_level_ratio_tumor.cct'), check.names = F)
rownames(glyco) = str_glue('{glyco$Modifications}_{glyco$Gene}')
glyco = glyco[,-(1:2)]
glyco = glyco[, samples_105]
phos = data.frame(fread('../PDAC_LinkedOmics_Data/phosphoproteomics_site_level_MD_abundance_tumor.cct'), check.names = F)
rownames(phos) = str_glue('{phos$Index}_{phos$Gene}')
phos = phos[,-(1:3)]
phos = phos[, samples_105]

fig7a_prepare_data = function(data_list){
  for (i in names(data_list)){
    data = data_list[[i]]
    # 0. Convert to ratio (if necessary)
    if (!(i %in% c('prot', 'glyco', 'phos'))){
      data = t(apply(data, 1, function(x){x/median(x)}))
    }
    # 1. Quantified in all tumors
    data_list[[i]] = data[complete.cases(data) & !is.infinite(rowSums(data)),]
  }
  data = do.call(rbind, data_list)
  # 2. SD < bottom 5th percentile
  stdev = apply(data, 1, sd)
  data = data[stdev >= quantile(stdev, 0.05),]
  # 3. Scale
  data = scale(data)
  # 4. Convert to non-negative matrix
  data_1 = data
  data_1[data_1 < 0] = 0
  data_2 = data
  data_2 = ifelse(data_2 > 0, 0, abs(data_2))
  data_final = rbind(data_1, data_2)
  data_final = data_final[apply(data_final, 1, function(x){sum(x == 0) != ncol(data_final)}),]
  
  return(data.frame(data_final, check.names = F))
}

nmf_data = fig7a_prepare_data(list('cna' = cna, 'mrna' = mrna, 'prot' = prot, 'glyco' = glyco, 'phos' = phos))
str_split_i(rownames(nmf_data), '\\.', 1) %>% table
apply(nmf_data, 1, function(x){sum(is.null(x)) == 105}) %>% which

# estim.r = nmf(nmf_data, 2:10, nrun = 1, seed = 1, maxIter = 50, .options = 'v')

nmf_res = nmf(nmf_data, 2, nrun = 500, seed = 1, .options='vP16', .pbackend = 'par')
saveRDS(nmf_res, '../etc/NMFresults_105samples_rank2.rds')
predict(nmf_res) %>% table
w = basis(nmf_res)
h = coef(nmf_res)

apply(h, 2, function(x){max(x/sum(x))})
