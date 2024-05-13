library(ComplexHeatmap)
library(ggplot2)
library(data.table)
library(reshape2)
library(ggrepel)
library(dplyr)
library(maftools)
library(ComplexHeatmap)
library(readxl)
library(RColorBrewer)
library(colorRamp2)
library(ggnewscale)
library(ggpubr)
library(gridExtra)

meta = fread('../PDAC_LinkedOmics_Data/clinical_table_140.tsv', sep = '\t')

# Figure 1B ####

df = data.frame()
tmp = data.frame(melt(table(ifelse(meta$participant_country %in% c('Canada', 'China', 'Poland', 'Russia', 'United States'), meta$participant_country, 'Other'))))
tmp = cbind(tmp, 'type'=rep('Country of Origin', nrow(tmp)))
tmp = tmp[c(4,6,2,1,5,3),]
df = rbind(df, tmp)
tmp = melt(table(meta$tumor_stage_pathological, useNA = 'ifany'))
tmp = cbind(tmp, 'type'=rep('Cancer stage', nrow(tmp)))
df = rbind(df, tmp)
tmp = melt(table(ifelse(is.na(meta$tumor_site), NA, ifelse(meta$tumor_site %in% c('head', 'head and body'), 'Head', 'Tail')), useNA = 'ifany'))
tmp = cbind(tmp, 'type'=rep('Tumor site', nrow(tmp)))
df = rbind(df, tmp)
tmp = melt(table(meta$vital_status, useNA = 'ifany'))
tmp = cbind(tmp, 'type'=rep('Vital status', nrow(tmp)))
df = rbind(df, tmp)

# df2 = df %>% 
#   mutate(csum = rev(cumsum(rev(value))), 
#          pos = value/2 + lead(csum, 1),
#          pos = if_else(is.na(pos), value/2, pos))
# df2 = df2[order(df2$value, decreasing = T),]

# png('../Figures/Figure1B.png', width = 28, height = 7, units = 'in', res = 300)
# par(mfcol = c(1, 4))
# for (i in c('Country of Origin', 'Cancer stage', 'Tumor site', 'Vital status')){
#   df2 = subset(df, type == i)
#   p = pie(df2$value, labels = paste0(df2$Var1, ': ', round(df2$value / sum(df2$value) * 100, 1), "%"), main = i, clockwise = T, radius = 0.75, cex = 1.5, cex.main = 2, line =-1)
#   print(p)
# }
# dev.off()

tumor_site_table = table(meta$tumor_site, useNA='ifany')/nrow(meta)
tumor_site_table = tumor_site_table[order(tumor_site_table, decreasing = T)]
pie(tumor_site_table, labels = paste0(names(tumor_site_table), ': ', round(tumor_site_table * 100, 1), "%"), main = 'Tumor site', clockwise = T, radius = 0.75, cex = 1.5, cex.main = 2, line =-1)

# ggplot(df2, aes(x="", y=value, fill=Var1)) +
#   geom_bar(stat="identity", width=1, color="white", show.legend = F) +
#   coord_polar("y", start=0) +
#   theme_void() +
#   labs(title = 'Country of Origin') +
#   theme(plot.title = element_text(hjust = 0.5))
# 

df2 = data.frame(tumor_site_table, check.names = F)
df2$Freq = -df2$Freq
df2 = df2[c(1,5,2,4,3,6),]
df2$Var1 = factor(df2$Var1, levels = df2$Var1)
df2 = df2 %>% mutate(csum = rev(cumsum(rev(Freq))), 
                pos = Freq/2 + lead(csum, 1),
                pos = if_else(is.na(pos), Freq/2, pos))


ggplot(df2, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size = 1, stat = "identity") +
  coord_polar(theta = "y") +
  geom_text_repel(aes(y = pos,
                       label = paste0(Var1, ': ', round(Freq / sum(Freq) * 100, 1), "%")),
                   size = 7,
                   nudge_x = 1,
                   show.legend = FALSE) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 0),
        legend.position = "none", # Removes the legend
        panel.background = element_rect(fill = "white"))

# Figure 1C #### 
## Mutation status
mut = fread('../PDAC_LinkedOmics_Data/Mutation_gene_level.cgt', header = T)
mut_sub = data.frame(subset(mut, V1 %in% c('TP53', 'CDKN2A', 'SMAD4')), check.names = F, row.names = 1)
mut_sub = mut_sub[c('SMAD4', 'CDKN2A', 'TP53'),]
mut_sub = data.frame(apply(mut_sub, c(1,2), function(x){ifelse(grepl('splice|stop_gained|frameshift|missense|nonsense|inframe', x), x %>% gsub('splice.*|missense.*', 'Missense', .) %>% gsub('stop_gained.*|nonsense.*|frameshift.*', 'Nonsense/frameshift', .) %>% gsub('inframe.*', 'In-frame ins/del', .), NA)}),
                     check.names = F)

## KRAS VAF
maf = read.maf('../PDAC_LinkedOmics_Data/Mutation_results.maf')
maf_kras = subsetMaf(maf, genes = 'KRAS')
hotspot = paste(c('G12', 'G13', 'Q61', 'K117', 'A146'), collapse = '|')
maf_kras_hotspot = subsetMaf(maf_kras, query = "grepl(hotspot, HGVSp_Short)")
maf_kras_hotspot = maf_kras_hotspot@data[,c('Tumor_Sample_Barcode', 'HGVSp_Short', 't_alt_count', 't_depth')]
maf_kras_hotspot$HGVSp_Short = gsub('p.', '', maf_kras_hotspot$HGVSp_Short)
kras_vaf = data.frame(maf_kras_hotspot %>% summarise(`sample`= `Tumor_Sample_Barcode`, `KRAS VAF`=t_alt_count/t_depth, `Hotspot`=`HGVSp_Short`) %>% group_by(`sample`) %>% filter(`KRAS VAF` == max(`KRAS VAF`)), check.names = F)
# kras_vaf = merge(kras_vaf, maf_kras_hotspot, by = 'sample')
kras_vaf = rbind(kras_vaf, data.frame('sample'=setdiff(colnames(mut_sub), kras_vaf$sample), 'KRAS VAF' = rep(0, 5), 'Hotspot' = rep(NA, 5), check.names = F))

## Mutation Burden
mut_burden = data.frame(maf@data%>% group_by(`Tumor_Sample_Barcode`) %>% summarise(`Mutation burden` = n()), check.names = F)

## Tumor estimate, CNV index
molecular = read_excel('../Table S1.xlsx', 3)
molecular = data.frame(molecular, check.names = F)
molecular$epithelial_cancer_deconv = as.numeric(sapply(molecular$epithelial_cancer_deconv, function(x){gsub('NA', NA, x)}))

## Histological type
meta = data.frame(meta, check.names = F)

## CNV

cnv_seg = fread('../PDAC_LinkedOmics_Data/SCNA_log2_segment_level.cct')
cnv_gene = fread('../PDAC_LinkedOmics_Data/SCNA_log2_gene_level.cct')

cnv_amp = apply(cnv_gene[,-1], 2, function(x){sum(ifelse(na.omit(x)>0.2, 1, 0))})
cnv_del = apply(cnv_gene[,-1], 2, function(x){sum(ifelse(na.omit(x)< -0.2, -1, 0))})

cnv_seg %>% group_by(Sample) %>% summarise('CNV amplifications'=sum(ifelse(log2_copyRatio > 0.2, 1, 0)), 'CNV deletions'=sum(ifelse(log2_copyRatio < -0.2, -1, 0))) %>% View

cnv_df = data.frame('CNV amplifications'=cnv_amp, 'CNV deletions'=cnv_del, check.names = F)
cnv_df$sample = rownames(cnv_df)

cnv_index = subset(cnv_seg, !(Chromosome %in% c('chrX', 'chrY'))) %>% 
  group_by(Sample, Chromosome) %>% summarise(`seg_len` = (`End`-`Start`+1), `seg_cnv_index`=`seg_len`*abs(log2_copyRatio)/sum(seg_len)) %>%
  group_by(Sample) %>% summarise(`CNV index`=sum(`seg_cnv_index`))

cnv_df = merge(cnv_df, cnv_index, by.x = 'sample', by.y = 'Sample')

# apply(cnv_gene, )

# cnv_index = cnv[grepl('chr[0-9]+', cnv$Chromosome),] %>% group_by(Sample) %>% summarise(`CNV index`=sum(log2_copyRatio))
# data.frame(cnv_index[!duplicated(cnv_index$Sample),-2], row.names = 1)

## Tumor cellularity
anno_df = kras_vaf
anno_df = merge(anno_df, mut_burden, by.x = 'sample', by.y = 'Tumor_Sample_Barcode')
anno_df = merge(anno_df, molecular[,c('case_id', 'epithelial_cancer_deconv', 'neoplastic_cellularity_histology_estimate')], by.x = 'sample', by.y = 'case_id')
anno_df = merge(anno_df, meta[, c('case_id', 'histology_diagnosis')], by.x = 'sample', by.y = 'case_id')
anno_df = merge(anno_df, cnv_df, by = 'sample')
anno_df$`Tumor cellularity` = ifelse(anno_df$`KRAS VAF`> 0.075 | (anno_df$`KRAS VAF` == 0 & (anno_df$`CNV index` > 1 | anno_df$`Mutation burden` > 25)), '>=15%, Sufficient Purity (n = 105)', '<15%, Low purity (n = 35)')

names(anno_df)[5:7] = c('Tumor estimate: deconvolution', 'Tumor estimate: histology', 'Histological type')
anno_df = anno_df[order(anno_df$`Tumor estimate: histology`, decreasing = T),]
anno_df = anno_df[order(anno_df$`Tumor cellularity`),]

fwrite(anno_df, '../Figures/Figure1_anno_df.csv')

## Top annotation
top_ha = HeatmapAnnotation(`CNV amplifications`=anno_barplot(anno_df$`CNV amplifications`, height = unit(30, 'mm'), border = F, gp = gpar(fill = brewer.pal(3, 'Set3')[1])),
                           `CNV deletions`=anno_barplot(anno_df$`CNV deletions`, height = unit(30, 'mm'), border = F, gp = gpar(fill = brewer.pal(3, 'Set3')[2])),
                           `Mutation burden` = anno_barplot(anno_df$`Mutation burden`, height = unit(30, 'mm'), border = F, add_numbers=T, gp = gpar(fill = brewer.pal(3, 'Set3')[3])),
                           `CNV index` = anno_df$`CNV index`,
                           annotation_legend_param = list(`CNV index`=list(color_bar='discrete')))



## Bottom annotation
col_vaf = brewer.pal(9, 'Greens')[1:4]
names(col_vaf) = seq(0, 0.4, length.out = 4)
col_dec = brewer.pal(9, 'PuRd')[1:4]
names(col_dec) = seq(0.1, 0.6, length.out = 4)
col_his = brewer.pal(9, 'Blues')[1:4]
names(col_his) = seq(0.1, 0.7, length.out = 4)
col_histype = brewer.pal(3, 'Set3')[1:2]
names(col_histype) = na.omit(unique(meta$histology_diagnosis))
col_cell = brewer.pal(3, 'Set1')[1:2]
names(col_cell) = c('>=15%, Sufficient Purity (n = 105)', '<15%, Low purity (n = 35)')

bottom_ha = HeatmapAnnotation(df = anno_df[,c(2,5,6,7,11)],
                              # `CNV index`= ,
                              # `KRAS VAF` = kras_vaf$`KRAS VAF`, 
                              # `Tumor estimate: deconvolution` = molecular$epithelial_cancer_deconv, 
                              # `Tumor estimate: histology` = molecular$neoplastic_cellularity_histology_estimate,
                              # `Histological type` = meta$histology_diagnosis,
                              # `Tumor cellularity` = cellularity,
                              # col = list(`KRAS VAF` = col_vaf,mut
                              #            `Tumor estimate: deconvolution` = col_dec,
                              #            `Tumor estimate: histology` = col_his,
                              #            `Histological type` = col_histype,
                              #            `Tumor cellularity` = col_cell),
                              # annotation_legend_param = list(`KRAS VAF` = list(title='KRAS VAF', color_bar = "discrete", at = names(col_vaf), labels = c(0, '', '', 0.4)),
                              #                                 `Tumor estimate: deconvolution` = list(title='Tumor estimate: deconvolution', color_bar = "discrete", at = names(col_dec), labels = c(0.1, '', '', 0.6)),
                              #                                 `Tumor estimate: histology` = list(title='Tumor estimate: histology', color_bar = "discrete", at = names(col_his), labels = c(0.1, '', '', 0.7)),
                              #                                 `Histological type` = list(title='Histological type', color_bar = "discrete"),
                              #                                 `Tumor cellularity` = list(title='Tumor cellularity', color_bar = "discrete")),
                              annotation_legend_param = list(`KRAS VAF`=list(color_bar='discrete'),
                                                             `Tumor estimate: deconvolution`=list(color_bar='discrete'),
                                                            `Tumor estimate: histology`=list(color_bar='discrete'),
                                                            `Histological type`=list(color_bar='discrete'),
                                                            `Tumor cellularity`=list(color_bar='discrete')),
                              na_col = 'white')
mut_sub = mut_sub[, as.character(anno_df$sample)]

png('../Figures/Figure1C.png', width=25, height=8, units ='in', res=300)                                         
draw(Heatmap(mut_sub, name = 'Mutation status', width = ncol(mut_sub)*unit(4, "mm"), height = nrow(mut_sub)*unit(4, "mm"),
        top_annotation = top_ha, bottom_annotation = bottom_ha, heatmap_legend_param = list(color_bar = "discrete"), na_col = 'white', cluster_rows = F, cluster_columns = F), 
     annotation_legend_side = 'bottom', heatmap_legend_side = 'bottom', legend_grouping = 'original', merge_legend=T)
dev.off()

# Figure 1D ####
png('../Figures/Figure1D.png', width=15, height=5, units ='in', res=300)
ggplot(anno_df, aes(x = reorder(sample, -`KRAS VAF`), y=`KRAS VAF`, fill = Hotspot)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(limits=c(-0.01, 0.5), expand=expansion(mult=c(0,0))) +
  # scale_fill_manual(guide = guide_legend(order = 1)) +
  new_scale_fill() +
  geom_tile(inherit.aes = F, aes(x = 1:nrow(anno_df), y = -0.005, fill=anno_df$`Tumor cellularity`[order(anno_df$`KRAS VAF`, decreasing = T)], height = 0.008)) +
  labs(x = 'Sample IDs', y = 'KRAS VAF', fill = 'Tumor cellularity') +
  geom_hline(yintercept = 0.075, linetype = 'dashed') +
  # coord_cartesian(clip = "off", ylim = c(-0.05, 0.5)) +
  theme(
    axis.text.x = element_text(size=8, angle=90, vjust=0.4, hjust=0),
    plot.margin = margin(0.3, 0.3, 1, 0.3, "cm"),
    axis.title.x = element_text(vjust=-0.5),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "darkgrey", fill=NA, size=1))
dev.off()

# Figure S1C ####
png('../Figures/FigureS1C.png', width=10, height=5, units ='in', res=300)

a = ggplot(anno_df, aes(x = `KRAS VAF`, y = `Tumor estimate: histology`)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  stat_cor(method = 'spearman', aes(label = paste0('Spearman ', ..r.label..)), output.type = 'text') +
  labs(x = 'KRAS VAF', y = 'Tumor estimate: histology') +
  theme_bw()
b = ggplot(anno_df, aes(x = `KRAS VAF`, y = `Tumor estimate: deconvolution`)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  stat_cor(method = 'spearman', aes(label = paste0('Spearman ', ..r.label..)), output.type = 'text') +
  labs(x = 'KRAS VAF', y = 'Tumor estimate: deconvolution') +
  theme_bw()

grid.arrange(a,b, nrow=1, ncol=2)
dev.off()
