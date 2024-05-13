library(maftools)
library(data.table)
library(OmnipathR)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pdftools)
library(parallel)
library(ltm)
library(ggplot2)
library(ggpubr)
library(EnrichmentBrowser)
library(stringr)

# Figure 2A #### 
genes = c('KRAS', 'TP53', 'CDKN2A', 'SMAD4', 'ARID1A', 'RNF43', 'GNAS', 'KMT2C', 'KMT2D', 'TGFBR2', 'RBM10')
samples = as.character(anno_df$sample[anno_df$`Tumor cellularity` == ">=15%, Sufficient Purity (n = 105)"])
# samples = c(as.character(anno_df[anno_df$`KRAS VAF` > 0.075,]$sample), 'C3N-00198', 'C3N-01380', 'C3N-01715', 'C3N-03426')
samples_orig = pdf_text('../purity_samples.pdf')


cnv_gene = fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct')
cnv_status = ifelse(cnv_gene[,-1] >= 0.4, 'CNV amplifications', ifelse(cnv_gene[,-1] <= -0.4, 'CNV deletions', NA))
cnv_status = data.frame(V1 = cnv_gene$V1, cnv_status, check.names = F)
cnv_status = melt(cnv_status, id.vars = 'V1')
cnv_status = cnv_status[!is.na(cnv_status$value),]
names(cnv_status) = c('Gene', 'Sample_name', 'CN')

mafandcn = read.maf(maf = '../PDAC_LinkedOmics_Data/Mutation_results.maf',
         cnTable = cnv_status,
         verbose = FALSE)

maf_sub = subsetMaf(mafandcn, genes = genes, query = "Tumor_Sample_Barcode %in% samples")
maf_sub@data$Variant_Classification %>% gsub('Missense_Mutation', 'Missense', .) %>% gsub('In_Frame.*', 'In-frame indel', .) %>% gsub('Frame_Shift.*', 'Frame shift', .) %>% gsub('Splice_Site', 'Splice site', .) %>% gsub('Nonsense_.*', 'Nonsense', .) -> maf_sub@data$Variant_Classification

hotspot = paste(c('G12', 'G13', 'Q61', 'K117', 'A146'), collapse = '|')

hotspot_df = subsetMaf(maf_sub, genes = 'KRAS', query = "grepl(hotspot, HGVSp_Short)")@data
hotspot_df$Hugo_Symbol = rep('KRAS hotspots', nrow(hotspot_df))
hotspot_df$Variant_Classification = sapply(hotspot_df$HGVSp_Short, function(x){sub('p.', '', x)})
hotspot_df$SYMBOL = rep('KRAS hotspots', nrow(hotspot_df))

maf_sub@data = rbind(maf_sub@data, hotspot_df)

mut_burden %>% filter(Tumor_Sample_Barcode %in% samples) %>% summarise(`Sample_Name`= `Tumor_Sample_Barcode`, `Log2 Mutation Burden`=log2(`Mutation burden`)) -> log2mutburden

col = c('#9aca4e', '#2cb67a', '#177196', '#343f7c', '#491653', '#ed1e22', '#3753a4',
        '#34549c', '#4eb78a', '#d3cd5d', '#c95836', '#d68750', '#50a050', '#CCCCCC', '#CCCCCC')

names(col) = c("Missense", "In-frame indel", "Frame shift", "Nonsense", "Splice site", "CNV amplifications", "CNV deletions",
               'G12D', 'G12V', 'G12R', 'G13D', 'Q61H', 'Q61R', 'Multi_Hit', 'Complex_Event')

genes_hotspot = c('KRAS hotspots', genes)
oncoplot(maf_sub, genes = genes_hotspot, colors =col, topBarData=log2mutburden, keepGeneOrder = T, drawRowBar=F, sampleOrder = samples, titleText = 'Sufficient tumor cellularity (n = 105)')

# Figure 2B #### 
rna = fread('../PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct')
prot = fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct')
rna = data.frame(rna, row.names = 1, check.names = F)
prot = data.frame(prot, row.names = 1, check.names = F)

int = import_omnipath_interactions(c('DEPOD', 'SIGNOR', 'Reactome_ProtMapper'))
int_sub = subset(int, grepl(paste0(genes, collapse = '|'), source_genesymbol) | grepl(paste0(genes, collapse = '|'), target_genesymbol))

## Interaction partner ####
corum = fread('../ref/humanComplexes.txt')
corum = subset(corum, grepl(paste0(genes, collapse = '|'), `subunits(Gene name)`))
int_list = list()
for (i in strsplit(corum$`subunits(Gene name)`, ';')){
  for (j in i){
    if (j %in% genes){
      if (j %in% names(int_list)){
        int_list[[j]] = c(int_list[[j]], i)
      } else {
        int_list[[j]] = i
      }
    }
  }
}

for (i in 1:nrow(int_sub)){
  source = strsplit(int_sub$source_genesymbol[i], '_')[[1]]
  target = strsplit(int_sub$target_genesymbol[i], '_')[[1]]
  for (s in source){
    if (s %in% names(int_list)){
      int_list[[s]] = c(int_list[[s]], source, target)
    }
  }
  for (t in target){
    if (t %in% names(int_list)){
      int_list[[t]] = c(int_list[[t]], target, source)
    }
  }
}

int_list = lapply(int_list, unique)


mut = data.frame(fread('../PDAC_LinkedOmics_Data/Mutation_gene_level.cgt', header = T), check.names = F, row.names = 1)
mut = mut[names(int_list),samples]
mut[mut != 'WT'] = "MUT"
apply(mut, 1, table)

gene_infig = c('ARID1A', 'ARID1B', 'CASP10', 'COL7A1', 'COPS6', 'EP300', 'FAS', 'FSTL3', 'JUN', 'MAPK3', 'MOXA1', 'POLD1', 'RBM10', 'RPS27L', 'SERPINE1', 'SMAD7', 'SMARCC2', 'STK17A', 'TOP2B')
rna_df = data.frame()
prot_df = data.frame()
for (gene in names(int_list)){
  mut = subset(maf_sub@data, Hugo_Symbol == gene)$Tumor_Sample_Barcode
  wt = setdiff(samples, mut)
  for (gene_ in int_list[[gene]]){
    if (gene_ %in% rownames(rna)){
      pval = wilcox.test(as.numeric(rna[gene_, mut]), as.numeric(rna[gene_, wt]))$p.value
      fc = 2 ** (median(as.numeric(rna[gene_, mut])) - median(as.numeric(rna[gene_, wt])))
      rna_df = rbind(rna_df, data.frame(gene, gene_, pval, fc))
    }
    if (gene_ %in% rownames(prot)){
      try({
        pval = wilcox.test(as.numeric(prot[gene_, mut]), as.numeric(prot[gene_, wt]))$p.value
        fc = median(as.numeric(prot[gene_, mut]))/median(as.numeric(prot[gene_, wt]))
        prot_df = rbind(prot_df, data.frame(gene, gene_, pval, fc))
      })
    }
  }
}

rna_df_ = data.frame()
prot_df_ = data.frame()
for (gene in names(int_list)){
  mut = subset(maf_sub@data, Hugo_Symbol == gene)$Tumor_Sample_Barcode
  wt = setdiff(samples, mut)
  for (gene_ in gene_infig){
    if (gene_ %in% rownames(rna)){
      pval = wilcox.test(as.numeric(rna[gene_, mut]), as.numeric(rna[gene_, wt]))$p.value
      fc = median(as.numeric(rna[gene_, mut])) - median(as.numeric(rna[gene_, wt]))
      rna_df_ = rbind(rna_df_, data.frame(gene, gene_, pval, fc))
    }
    if (gene_ %in% rownames(prot)){
      try({
        pval = wilcox.test(as.numeric(prot[gene_, mut]), as.numeric(prot[gene_, wt]))$p.value
        fc = median(as.numeric(prot[gene_, mut]))/median(as.numeric(prot[gene_, wt]))
        prot_df_ = rbind(prot_df_, data.frame(gene, gene_, pval, fc))
      })
    }
  }
}

rna_df_$fdr = p.adjust(rna_df_$pval, method = 'fdr')
prot_df_$fdr = p.adjust(prot_df_$pval, method = 'fdr')

subset(rna_df_, fdr < 0.05)

# Figure 2D ####
cnv_gene = fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct'))
cnv_gene = cnv_gene[,c('V1', as.character(samples))]

cnv_status = ifelse(cnv_gene[,-1] >= 0.4, 'amp', ifelse(cnv_gene[,-1] <= -0.4, 'del', NA))
genes = list('Amplification' = c('GATA6', 'AKT2', 'MYC', 'KRAS', 'ERBB2'), 'Deletion' = c('CDKN2A', 'SMAD4', 'ARID1A'))
df = data.frame()
for (i in seq_along(genes)){
  idx = match(genes[[i]], cnv_gene$V1)
  perc = apply(cnv_status[idx,] == tolower(substr(names(genes)[i], 1, 3)) & !is.na(cnv_status[idx,]), 1, sum)/ncol(cnv_status) * 100
  df = rbind(df, data.frame('gene' = genes[[i]], 'Percentage' = perc, 'CNV status' = names(genes)[i], check.names = F))
}

png('../Figures/Figure2D.png', width = 10, height = 4, units = 'in', res = 300)
ggplot(df, aes(x = gene, y = Percentage, fill = `CNV status`)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = rev(df$gene)) +
  coord_flip() +
  theme_classic() +
  labs(x = '')
dev.off()

# Figure 2E ####
cnv_seg = fread('../PDAC_LinkedOmics_Data/SCNA_log2_segment_level.cct')
cnv_seg_sub = subset(cnv_seg, Sample %in% samples)
fwrite(cnv_seg_sub, '../PDAC_LinkedOmics_Data/SCNA_log2_segment_level_105samples.txt', sep = '\t')

lesion = fread('../577862/all_lesions.conf_95.txt')
cytobands = subset(lesion, `q values` < 0.05)$Descriptor
score = fread('../577862/scores.gistic')
score_fdr0.05 = subset(score, `-log10(q-value)` > -log10(0.05))

fwrite(score_fdr0.05, '../577862/scores_fdr0.05.gistic')
gistic = readGistic(gisticAllLesionsFile = '../577862/all_lesions.conf_95.txt', gisticAmpGenesFile = '../577862/amp_genes.conf_95.txt', 
                    gisticDelGenesFile = '../577862/del_genes.conf_95.txt', gisticScoresFile = '../577862/scores_fdr0.05.gistic')

cyto = c('1p12', '9p21.3', '12p12.1', '18q11.2')
cytoband = cyto_peaks_scores[match(cyto, cyto_peaks_scores$Cytoband),]
cytoband$Cytoband = c('NOTCH2', 'CDKN2A', 'KRAS', 'GATA6')
cytoband$amp = sign(cytoband$amp) * (abs(cytoband$amp) + 0.02)
cyto_peaks_scores = rbind(cyto_peaks_scores, cytoband)

png('../Figures/Figure2E.png', width = 10, height = 5, units = 'in', res = 300)
gisticChromPlot(gistic, ref.build = 'hg38', fdrCutOff = 0.05, txtSize = 0.8) +
  text(x = cyto_peaks_scores$Start_Position_updated, y = (abs(cyto_peaks_scores$amp) + 0.01) * sign(cyto_peaks_scores$amp),
       labels = cyto_peaks_scores$Cytoband, font = 3, cex = 0.8)
dev.off() 


# Figure 2F #### 
cnv_gene = data.frame(fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct'), check.names = F, row.names = 1)
prot = data.frame(fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct'), check.names = F, row.names = 1)
rna = data.frame(fread('../PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct'), check.names = F, row.names = 1)

cnv_gene_105 = cnv_gene[, samples]
prot_105 = prot[, samples]
rna_105 = rna[, samples]

# 1st filtering: Quantifiable
genes = intersect(names(which(apply(cnv_gene_105, 1, function(x){sum(is.na(x)) != length(samples)}))), names(which(apply(prot_105, 1, function(x){sum(is.na(x)) != length(samples)}))))
genes = intersect(genes, names(which(apply(rna_105, 1, function(x){sum(is.na(x)) != length(samples)}))))

# 2nd filtering: Located in focal amplified regions (GISTIC Q value < 0.25)
# s2 = fread('../577862/amp_genes.conf_95.txt')
s2 = read_excel('../Table S2. Copy-number variation and methylation analysis, related to Figures 2 and S2.xlsx', 3)
amp_genes = unique(na.omit(unlist(c(s2[4:nrow(s2),-1]))))
amp_genes = amp_genes[!startsWith(amp_genes, '[')]
genes_1 = intersect(genes, amp_genes)

# 3rd filtering: Correlation with CN-RNA, CN-Protein (Spearman FDR < 0.05)
cor_df = data.frame()
cnv_genes = names(which(apply(cnv_gene_105, 1, function(x){sum(is.na(x)) != length(samples)})))
rna_genes = names(which(apply(rna_105, 1, function(x){sum(is.na(x)) != length(samples)})))
prot_genes = names(which(apply(prot_105, 1, function(x){sum(is.na(x)) != length(samples)})))

## Cis correlation ####

mc <- getOption("mc.cores",20)

# RNA
cal_rna_cor = function(cnv, rna){
  rna_test = cor.test(t(cnv_gene_105[cnv,]), t(rna_105[rna,]), method = 'spearman', exact = F)
  pval = rna_test$p.value
  rho = rna_test$estimate
  return(data.frame(t(c(cnv, rna, rho, pval))))
}
cnv_rna_genes = intersect(cnv_genes, rna_genes)
rna_cor = mcmapply(cal_rna_cor, cnv = cnv_rna_genes, rna = cnv_rna_genes, mc.cores = mc, SIMPLIFY = F)
rna_cor = rbindlist(rna_cor)
names(rna_cor) = c('CNV', 'RNA', 'Spearman', 'PValue')
rna_cor$PValue = as.numeric(rna_cor$PValue)

# Protein
cal_prot_cor = function(cnv, prot){
  tmp = prot_105[prot,]
  tmp[is.na(tmp)] = 0 # many NAs -> replaced NA to 0
  prot_test = cor.test(t(cnv_gene_105[cnv,]), t(tmp), method = 'spearman', exact = F)
  pval = prot_test$p.value
  rho = prot_test$estimate
  return(data.frame(t(c(cnv, prot, rho, pval))))
}
cnv_prot_genes = intersect(cnv_genes, prot_genes)
prot_cor = mcmapply(cal_prot_cor, cnv = cnv_prot_genes, prot = cnv_prot_genes, mc.cores = mc, SIMPLIFY = F)
prot_cor = rbindlist(prot_cor)
names(prot_cor) = c('CNV', 'Protein', 'Spearman', 'PValue')
prot_cor$PValue = as.numeric(prot_cor$PValue)

intersect(amp_genes, subset(rna_cor, PValue <= max(s2_rna_cor$PValue))$CNV) # n = 202
final_genes_mine = intersect(intersect(amp_genes, subset(rna_cor, PValue <= max(s2_rna_cor$PValue))$CNV), subset(prot_cor, PValue <= max(s2_prot_cor$PValue))$CNV)

## Cis + Trans correlation ####
rna_cor = data.frame()
prot_cor = data.frame()

rna_comb = data.frame()
for (i in cnv_genes){
  for (j in rna_genes){
    rna_cor = rbind(rna_cor, data.frame(t(mcmapply(cal_rna_cor, cnv = i, rna = rna_genes, mc.cores = mc)), check.names = F))
    prot_cor = rbind(prot_cor, data.frame(t(mcmapply(cal_prot_cor, cnv = i, prot = prot_genes, mc.cores = mc)), check.names = F))
}

names(cor_df) = c('CNV', 'RNA/Protein', 'RNA_pvalue', 'RNA_cor', 'Protein_pvalue', 'Protein_cor')
cor_df$RNA_FDR = p.adjust(cor_df$RNA_pvalue, method = 'BH')
cor_df$Protein_FDR = p.adjust(cor_df$Protein_pvalue, method = 'BH')

## Use Table S2 ####
s2_rna_cor = read_excel('../Table S2. Copy-number variation and methylation analysis, related to Figures 2 and S2.xlsx', 6)
s2_rna_cor = rbind(s2_rna_cor, read_excel('../Table S2. Copy-number variation and methylation analysis, related to Figures 2 and S2.xlsx', 7))
s2_prot_cor = read_excel('../Table S2. Copy-number variation and methylation analysis, related to Figures 2 and S2.xlsx', 8)

intersect(amp_genes, subset(s2_rna_cor, RNA == CNV & `Pvalue.adjusted` < 0.05)$CNV) # n = 95
final_genes = intersect(intersect(amp_genes, subset(s2_rna_cor, RNA == CNV & `Pvalue.adjusted` < 0.05)$CNV), subset(s2_prot_cor, Protein == CNV & `Pvalue.adjusted` < 0.05)$CNV) # n = 20

subset(s2_rna_cor, RNA == CNV & `Pvalue.adjusted` < 0.05)$CNV # n = 2914

## Compare mine and Table S2 ####
merged_rna = merge(rna_cor, subset(s2_rna_cor, RNA == CNV & `Pvalue.adjusted` < 0.05), by = 'RNA', all.x = T)
p = list()
p[[1]] = ggplot(merged_rna, aes(x = as.numeric(Spearman.x), y = as.numeric(Spearman.y))) +
  geom_point() +
  labs(x = 'RNA Spearman (Calculated)', y = 'RNA Spearman (Table S2)') +
  theme_bw()
p[[2]] = ggplot(merged_rna, aes(x = -log(as.numeric(PValue.x), 10), y = -log(as.numeric(PValue.x), 10))) +
  geom_point() +
  labs(x = 'RNA -log10(PValue) (Calculated)', y = 'RNA -log10(PValue) (Table S2)') +
  theme_bw()
merged_prot = merge(prot_cor, subset(s2_prot_cor, Protein == CNV & `Pvalue.adjusted` < 0.05), by = 'Protein', all.x = T)
p[[3]] = ggplot(merged_prot, aes(x = as.numeric(Spearman.x), y = as.numeric(Spearman.y))) +
  geom_point() +
  labs(x = 'Protein Spearman (Calculated)', y = 'Protein Spearman (Table S2)') +
  theme_bw()
p[[4]] = ggplot(merged_prot, aes(x = -log(as.numeric(PValue.x), 10), y = -log(as.numeric(PValue.x), 10))) +
  geom_point() +
  labs(x = 'Protein -log10(PValue) (Calculated)', y = 'Protein -log10(PValue) (Table S2)') +
  theme_bw()

ggplot(merged_prot, aes(x = Protein, y = PValue.x - PValue.y)) +
  geom_point() +
  labs(x = 'Protein -log10(PValue) (Calculated)', y = 'Protein -log10(PValue) (Table S2)') +
  theme_bw()

ggplot(merged_prot, aes(x = Protein, y = as.numeric(Spearman.x) - as.numeric(Spearman.y))) +
  geom_point() +
  labs(x = 'Protein -log10(PValue) (Calculated)', y = 'Protein -log10(PValue) (Table S2)') +
  theme_bw()

png('../Figures/Figure2F_Compare_Cal_S2.png', width = 10, height = 10, res = 300, units = 'in')
ggarrange(plotlist = p, nrow = 2, ncol = 2)
dev.off() 

library(ggvenn)
png('../Figures/Figure2F_Venn_Cal_S2.png', width = 8, height = 4, res = 300, units = 'in')
ggvenn(data = list('Table S2' = final_genes, 'Calculated' = final_genes_mine))
dev.off()

## GSA ####
final_genes_merged = union(final_genes, final_genes_mine)
final_genes_merged = union(final_genes, genes_2g)

ego = enrichGO(gene = final_genes_merged, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = 'BP', 
               pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize = 1200, universe = rownames(cnv_gene_105))

enriched_term = c('GO:0030029', 'GO:0007010', 'GO:0044087', 'GO:0032970', 'GO:0033043', 'GO:0032271', 'GO:0051493', 'GO:0051258')

subset(ego@result, p.adjust < 0.05)
ego@result[enriched_term[enriched_term %in% ego@result$ID],]

setdiff(s2_rna_cor$RNA %>% unique, rownames(cnv_gene_105))

get_geneset_Env <- function (term) {
  if (!exists(str_glue(".{term}_clusterProfiler_Env"), envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(str_glue(".{term}_clusterProfiler_Env"), new.env(), envir=envir)
  }
  get(str_glue(".{term}_clusterProfiler_Env"), envir = .GlobalEnv)
}

GO_db = get_geneset_Env('GO')
GO_db_list = split(GO_db$goAnno$SYMBOL, GO_db$goAnno$GOALL)
lapply(GO_db_list[enriched_term], function(x){genes_2g %in% x})
lapply(GO_db_list[enriched_term], length)

# Figure 2G ####
genes_2g = c('XPO5', 'USP14', 'KRAS', 'NOTCH2', 'CDK6', 'POLB')
setdiff(genes_2g, final_genes) # n = 4

cnv_df = data.frame(apply(cnv_gene_105[genes_2g,], 1, function(x){ifelse(x < 0.4, 'Control', 'CNV')}), check.names = F)
rna_df = t(rna_105[genes_2g,])
colnames(rna_df) = paste0(genes_2g, '_RNA')
prot_df = t(prot_105[genes_2g,])
colnames(prot_df) = paste0(genes_2g, '_Protein')
all_df = cbind(cnv_df, rna_df, prot_df)

p = list()
for (i in 1:length(genes_2g)){
  p[[i]] = ggplot(all_df, aes(x = factor(.data[[genes_2g[i]]], levels = c('Control', 'CNV')), y = scale(.data[[paste0(genes_2g[i], '_RNA')]]), fill = factor(.data[[genes_2g[i]]], levels = c('Control', 'CNV')))) +
    geom_violin(trim = F) +
    geom_point(size=1, position = position_jitterdodge(), color="black",alpha=1) +
    scale_fill_manual(values = list('Control' = '#4f749d', 'CNV' = '#d08a51')) +
    guides(fill=guide_legend(title="")) +
    labs(x = genes_2g[i]) +
    scale_y_continuous(name ='', breaks = seq(0,8, by=4), limits = c(-2, 8.1)) +
    stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 8) +
    stat_compare_means(method = "t.test", label.x = 1, label.y = 7.3) +
    theme_bw()
}
for (i in 1:length(genes_2g)){
  p[[i+length(genes_2g)]] = ggplot(all_df, aes(x = factor(.data[[genes_2g[i]]], levels = c('Control', 'CNV')), y = scale(.data[[paste0(genes_2g[i], '_Protein')]]), fill = factor(.data[[genes_2g[i]]], levels = c('Control', 'CNV')))) +
    geom_violin(trim = F) +
    geom_point(size=1, position = position_jitterdodge(), color="black",alpha=1) +
    scale_fill_manual(values = list('Control' = '#4f749d', 'CNV' = '#d08a51')) +
    guides(fill=guide_legend(title="")) +
    labs(x = genes_2g[i]) +
    scale_y_continuous(name ='', breaks = seq(-5,5, by=5), limits = c(-6, 10)) +
    stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 10) +
    stat_compare_means(method = "t.test", label.x = 1, label.y = 9.3) +
    theme_bw()
}


png('../Figures/Figure2G.png', width = 20, height = 8, res= 300, units = 'in')
figure = ggarrange(plotlist = p, nrow = 2, ncol = length(genes_2g), common.legend = TRUE, legend = "top")
annotate_figure(figure, left = textGrob("Expression Z-score", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
dev.off()


# Figure 2H ####
cnv = paste0(c('CDKN2A', 'SMAD4'), '_CNV')
genes = c('ARID1B', 'CDKN2A', 'FOXO3', 'HDAC1', 'JUN', 'NCOA3', 'SMAD4', 'SMAD7', 'WWTR1')
merged = Reduce(merge(by = 'row.names'), list(apply(cnv_gene_105[cnv,], 1, function(x){factor(ifelse(x < -0.4, 'CNV', 'Ctrl'))}), t(rna_105[genes,]), t(prot_105[genes,])))





      