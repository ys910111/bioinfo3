library(data.table)
library(ggVolcano)

# Figure 3A ####

prot = data.frame(fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct'), check.names = F, row.names = 1)
nat = data.frame(fread('../PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_normal.cct'), check.names = F, row.names = 1)
samples_41 = intersect(samples_100, colnames(nat))

prot_41 = prot[,samples_41]
nat_41 = nat[,samples_41]

results = data.frame()
for (i in rownames(prot_41)){
  if (sum(!is.na(prot_41[i,])) > ncol(prot_41)/2 & sum(!is.na(nat_41[i,])) > ncol(prot_41)/2){
    P = wilcox.test(t(prot_41[i,]), t(nat_41[i,]), paired = T)$p.value
    nona = intersect(rownames(na.omit(t(prot_41[i,]))), rownames(na.omit(t(nat_41[i,]))))
    fc = median(t(prot_41[i,nona])) - median(t(nat_41[i,nona]))
    results = rbind(results, data.frame(i, fc, P))
  }
}
names(results) = c('Protein', 'fc', 'P')
results$FDR = p.adjust(results$P, 'BH')
results$sig = ifelse(results$FDR >= 0.01, 'FDR>=0.01', 
                     ifelse(results$fc < -1, '>2X Down',
                            ifelse(results$fc > 1, '>2X Up', ifelse(sign(results$fc) == -1, 'Down', 'Up'))))

ggplot(results, aes(x= fc, y = -log10(FDR), color = sig)) +
  geom_point() +
  scale_color_manual(values = c('FDR>=0.01'='grey', '>2X Down'='blue', '>2X Up'='red', 'Down'='#87cefa', 'Up'='#ffbfcb')) +
  xlim(-4, 4) + ylim(0, 11) +
  theme_bw() 

results$sig %>% table 


