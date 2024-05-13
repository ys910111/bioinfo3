library(data.table)
library(dplyr)

samples = colnames(mut)[c(9, 10, 40)]
mut = fread('../PDAC_LinkedOmics_Data/Mutation_site_level.cgt')
mut_gene = fread('../PDAC_LinkedOmics_Data/Mutation_gene_level.cgt', header = T)

data.frame(mut, check.names = F)[, c('site', samples)] %>% View
data.frame(mut_gene, check.names = F)[, c('V1', samples)] %>% View

cnv = fread('../PDAC_LinkedOmics_Data/SCNA_log2_segment_level.cct')
cnv_gene = fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct')
View(data.frame(cnv_gene, check.names = F)[, c('V1', 'C3L-00017', samples[c(1,2)])])

rna = fread('../PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct')
   
data.frame(rna, check.names = F)[, c('V1', samples)] %>% View     

mirna = fread('../PDAC_LinkedOmics_Data/microRNA_TPM_log2_Tumor.cct')

data.frame(mirna, check.names = F)[, c('V1', samples)] %>% View     

methyl = fread('../PDAC_LinkedOmics_Data/methylation_betaValue_Tumor.cct')

data.frame(methyl, check.names = F)[, c('V1', samples)] %>% View   

prot = fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct') 

data.frame(prot, check.names = F)[, c('V1', samples)] %>% View 

phos = fread('../PDAC_LinkedOmics_Data/phosphoproteomics_site_level_MD_abundance_tumor.cct')
phos_gene = fread('../PDAC_LinkedOmics_Data/phosphoproteomics_gene_level_MD_abundance_tumor.cct')

data.frame(phos, check.names = F)[, c('Index', 'Gene', 'Peptide', samples)] %>% View 
data.frame(phos_gene, check.names = F)[, c('V1', samples)] %>% View 

gly = fread('../PDAC_LinkedOmics_Data/N-glycoproteomics_Site_level_ratio_tumor.cct')
gly_gene = fread('../PDAC_LinkedOmics_Data/N-glycoproteomics_peptide_level_ratio_tumor.cct')

data.frame(gly, check.names = F)[810:817,c('Modifications', 'Gene', samples[1])] %>% View 
data.frame(gly_gene, check.names = F)[30600:30607, c('Sequence','Gene', samples)] %>% View 

