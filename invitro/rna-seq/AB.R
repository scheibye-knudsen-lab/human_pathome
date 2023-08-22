# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport", ask=FALSE)
# BiocManager::install("DESeq2", ask=FALSE)
# BiocManager::install("EnhancedVolcano", ask=FALSE)
library(tximport)
library(DESeq2)
library(EnhancedVolcano)
# install.packages("readr")
library(readr)

# library(ggplot2)
# library(tximportData)
# library(PCAtools)
# library(stringr)
# # library(ggfortify)
# # install.packages("VennDiagram")
# library(VennDiagram)

quant = list.dirs('quants', recursive=F, full.names=F); quant
quant_files = quant
quant = paste('quants/', quant, '/quant.sf', sep=''); quant

quant_files[quant_files=='222C_DMSOIIIIII_quant'] = '222C_DMSO_IIIIII_quant'; quant_files
quant_files = gsub('_quant', '', quant_files); quant_files

cell_line = unlist(lapply(strsplit(quant_files, '_'), FUN=function(x){x[1]})); cell_line = factor(cell_line); cell_line
condition = unlist(lapply(strsplit(quant_files, '_'), FUN=function(x){x[2]})); condition = factor(condition); condition
replicate = unlist(lapply(strsplit(quant_files, '_'), FUN=function(x){x[3]})); replicate = factor(replicate); replicate
run = quant_files; run = factor(run); run

dose = rep(0, length(quant_files)); dose
dose[grepl('1', condition)] = 1; dose
dose[grepl('10', condition)] = 10; dose

drug = rep('TBD', length(quant_files)); drug
drug[grepl('A', condition)] = 'A'; drug
drug[grepl('B', condition)] = 'B'; drug
drug[grepl('DMSO', condition)] = 'DMSO'; drug
drug = factor(drug); drug

treatment = rep('IR', length(cell_line)); treatment
treatment[grepl('C', cell_line)] = 'NonIR'; treatment
treatment = factor(treatment); treatment

group = paste(treatment, condition, sep='_'); group = factor(group); group

data = data.frame(run=run, condition=condition, replicate=replicate, dose=dose, drug=drug, quant=quant, group=group, treatment=treatment); data


comparison_generator = function(condition_, group_, treatment_,
                                condition_case, group_case, treatment_case,
                                condition_ctrl, group_ctrl) {
  phenotype_ = rep('other', length(condition_)[1])
  phenotype_case = (condition_==condition_case)&(group_==group_case)&(treatment_==treatment_case)
  phenotype_case_name = paste(c(condition_case, group_case, treatment_case), collapse = '_')
  phenotype_[phenotype_case==TRUE] = phenotype_case_name
  phenotype_ctrl = (condition_==condition_ctrl)&(group_==group_ctrl)
  phenotype_ctrl_name = paste(c(condition_ctrl, group_ctrl), collapse = '_')
  phenotype_[phenotype_ctrl==TRUE] = phenotype_ctrl_name
  phenotype_name = paste(c(phenotype_case_name, 'vs', phenotype_ctrl_name), collapse = '_')
  phenotype_ = factor(phenotype_)
  #print(phenotype_name)
  return (phenotype_)
}

comparison_generator = function(data_, case_, control_) {
  comparison_ = rep('other', length(data$group))
  comparison_[data$group==case_]=case_
  comparison_[data$group==control_]=control_
  data[paste(c(case_, 'vs', control_), collapse = '_')]=factor(comparison_)
  return (data)
}

data = comparison_generator(data, 'IR_DMSO', 'NonIR_DMSO')
data = comparison_generator(data, 'IR_B10', 'IR_DMSO')
data = comparison_generator(data, 'IR_B1', 'IR_DMSO')
data = comparison_generator(data, 'IR_A10', 'IR_DMSO')
data = comparison_generator(data, 'IR_A1', 'IR_DMSO')
data = comparison_generator(data, 'NonIR_B1', 'NonIR_DMSO')
data = comparison_generator(data, 'NonIR_B10', 'NonIR_DMSO')
data = comparison_generator(data, 'NonIR_A1', 'NonIR_DMSO')
data = comparison_generator(data, 'NonIR_A10', 'NonIR_DMSO')

tx2gene = read.table('refgenie/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_txp2gene.tsv', sep = '\t', header=FALSE, nrows = 189154); names(tx2gene) = c("GENEID", "TXNAME"); head(tx2gene)
files = data$quant; names(files) = data$run; head(files)
txi = tximport(files, type="salmon", tx2gene=tx2gene)
dds = DESeqDataSetFromTximport(txi, colData = data, design = ~1); dds = DESeq(dds)
keep <- rowSums(counts(dds)) >= 10; dds <- dds[keep,]
vsd = vst(dds, blind=FALSE)
vsd = vst(dds, blind=FALSE); plotPCA(vsd, intgroup=c('treatment', 'condition'))

write.table(plotPCA(vsd, returnData=T), file='pca.tsv', sep = '\t')


# Results

comparison_names = names(data)[9:17]; comparison_names
for (comparison_name in comparison_names) {
  print(comparison_name)
  fmla = as.formula(paste('~',comparison_name)); fmla
  dds = DESeqDataSetFromTximport(txi, colData = data, design = fmla); dds = DESeq(dds)
  keep = rowSums(counts(dds)) >= 10; dds = dds[keep,]
  res = results(dds, contrast=c(comparison_name, unlist(strsplit(comparison_name, '_vs_'))[1], unlist(strsplit(comparison_name, '_vs_'))[2]))
  save(res, file = paste(c('results/', comparison_name, '.RData'), collapse = ''))
}

# Volcanos
x = comparison_names[1]
load(file = paste(c('results/', x, '.RData'), collapse = ''))
png(paste(c('results/', x, '_padj_R.png'), collapse = ''), width=1000, height=500)
plot(res$log2FoldChange, -log10(res$padj), xlab='', ylab='', pch='.', cex=3)
dev.off()

sapply(comparison_names, FUN=function(x){
  load(file = paste(c('results/', x, '.RData'), collapse = ''))
  png(paste(c('results/', x, '_padj_R.png'), collapse = ''), width=500, height=500)
  plot(res$log2FoldChange, -log10(res$padj), xlab='', ylab='', pch='.', cex=3, bty='l')
  dev.off()
})

sapply(comparison_names, FUN=function(x){
  load(file = paste(c('results/', x, '.RData'), collapse = ''))
  png(paste(c('results/', x, '_padj_EV.png'), collapse = ''), bg = "transparent", res=100)
  
  plot(EnhancedVolcano(res,
                     lab = rownames(res),
                     pCutoff = 0.05,
                     FCcutoff = 2,
                     title = x,
                     cutoffLineType = 'blank',
                     cutoffLineCol = 'black',
                     
                     selectLab = FALSE,
                     x = 'log2FoldChange',
                     y = 'padj',
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE))
  
  dev.off()
})

sapply(comparison_names, FUN=function(x){
  load(file = paste(c('results/', x, '.RData'), collapse = ''))
  png(paste(c('results/', x, '_pvalue.png'), collapse = ''))  
  plot(EnhancedVolcano(res,
                       lab = rownames(res),
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       title = x,
                       x = 'log2FoldChange',
                       y = 'pvalue'))
  
  dev.off()
  png(paste(c('results/', x, '_padj.png'), collapse = ''))  
  plot(EnhancedVolcano(res,
                       lab = rownames(res),
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       title = x,
                       x = 'log2FoldChange',
                       y = 'padj'))
  
  dev.off()
})

# GSEA files

expression_gct = 'results/expression.gct'
df = cbind.data.frame(data.frame(Description=rep('na',dim(dds)[1])), as.data.frame(counts(dds, normalized = TRUE)))
df = cbind.data.frame(data.frame(NAME = row.names(df)), df)
cat(paste('#1.2',paste(dim(dds),collapse='\t'),sep='\n'),file=paste(expression_gct,sep=''),sep="\n")
write_tsv(df,file=paste(expression_gct,sep=''),col_names=T,append=T)

phenotype_cls = function(phenotype_name, sample_phenotypes) {
  write_file(
    paste(paste(length(sample_phenotypes),length(unique(sample_phenotypes)),'1',sep=' '),
          paste('#',paste(unique(sample_phenotypes),collapse=' '),sep=' '),
          paste(sample_phenotypes,collapse=' '),sep='\n'),
    paste('results/',phenotype_name,'.cls',sep='')
  )
}

phenotype_cls('dose', data$dose)
phenotype_cls('group', data$group)
phenotype_cls('condition', data$condition)
phenotype_cls('treatment', data$treatment)
phenotype_cls('drug', data$drug)

gsea_cli = function(phenotype) {
  cat(
    paste(c('sh /Volumes/auditgroupdirs/SUND-ICMM-Vita/github/GSEA-cli/gsea-cli.sh GSEA -res /Users/lxc844/github/rna-seq/AB/results/expression.gct -cls /Users/lxc844/github/rna-seq/AB/results/group.cls#', phenotype,' -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label ', phenotype, ' -metric Signal2Noise -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Gene_Symbol_with_Remapping_MSigDB.v2022.1.Hs.chip -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /Users/lxc844/github/rna-seq/AB/results/gsea'), collapse = ''),
    file='results/gsea.sh', sep='\n', append = T)
}

gsea_cli('IR_DMSO_versus_NonIR_DMSO')

gsea_cli('IR_B10_versus_IR_DMSO')
gsea_cli('IR_B1_versus_IR_DMSO')
gsea_cli('IR_A10_versus_IR_DMSO')
gsea_cli('IR_A1_versus_IR_DMSO')

gsea_cli('NonIR_B1_versus_NonIR_DMSO')
gsea_cli('NonIR_B10_versus_NonIR_DMSO')
gsea_cli('NonIR_A1_versus_NonIR_DMSO')
gsea_cli('NonIR_A10_versus_NonIR_DMSO')


# comparisons = c(
#   'IR_B10_versus_IR_DMSO',
#   'IR_A10_versus_IR_DMSO',
#   'IR_DMSO_versus_NonIR_DMSO',
#   'NonIR_B10_versus_NonIR_DMSO',
#   'NonIR_A10_versus_NonIR_DMSO'
# )


# GSEA Analysis

gsea = list()
gsea_results = list.dirs('results/gsea', recursive=F, full.names=T); gsea_results
gsea_results = gsea_results[grepl('.Gsea.', gsea_results)]; gsea_results
for (gsea_result in gsea_results) {
  gsea_result_dir = gsea_result
  gsea_result_name = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(gsea_result)))
  gsea_files = list.files(gsea_result_dir, recursive=F, full.names=F)
  gsea_files = gsea_files[grepl('tsv', gsea_files)]
  gsea_files = gsea_files[grepl('gsea_report', gsea_files)]
  gsea_files_full_path = list.files(gsea_result_dir, recursive=F, full.names=T)
  gsea_files_full_path = gsea_files_full_path[grepl('tsv', gsea_files_full_path)]
  gsea_files_full_path = gsea_files_full_path[grepl('gsea_report', gsea_files_full_path)]
  gsea_files = gsea_files[grepl('tsv', gsea_files)]
  gsea_files = gsub('gsea_report_for_', '', gsea_files)
  gsea_files = unlist(lapply(strsplit(gsea_files, '_'), FUN=function(x){paste(x[1:(length(x)-1)], collapse='_')}))
  gsea[[gsea_result_name]] = list()
  gsea[[gsea_result_name]][[gsea_files[1]]] = read.table(gsea_files_full_path[1], sep='\t', h=T, check.names = F)
  gsea[[gsea_result_name]][[gsea_files[2]]] = read.table(gsea_files_full_path[2], sep='\t', h=T, check.names = F)  
}



# for (result_name1 in names(gsea)) {
#   print(result_name1)
#   for (result_name2 in names(gsea[[result_name1]])) {
#     print(result_name2)    
#     print(dim(gsea[[result_name1]][[result_name2]]))
#   }
# }





B10_DOWN = read.table('results/gsea/B10_DOWN.tsv', sep = '\t', h=T); dim(B10_DOWN)
B10_DOWN = B10_DOWN[c('NAME', 'NES', 'FDR.q.val')]
B10_DOWN = B10_DOWN[B10_DOWN$FDR.q.val<0.05,]
B10_DOWN = B10_DOWN[order(B10_DOWN$NES),]
B10_DOWN

B10_UP = read.table('results/gsea/B10_UP.tsv', sep = '\t', h=T); dim(B10_UP)
B10_UP = B10_UP[c('NAME', 'NES', 'FDR.q.val')]
B10_UP = B10_UP[B10_UP$FDR.q.val<0.05,]
B10_UP = B10_UP[order(B10_UP$NES),]
B10_UP

B10 = rbind(B10_UP, B10_DOWN)
B10 = B10 = B10[order(B10$NES, decreasing = T),]
write.table(B10, file = 'results/gsea/S1.tsv', sep='\t', quote=F, row.names = F, col.names = T)


NONIR_B10 = read.table('results/gsea/DOWN_IN_NonIR_B10_versus_NonIR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(NONIR_B10)
IR_B10 = read.table('results/gsea/DOWN_IN_IR_B10_versus_IR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(IR_B10)

# sum(NONIR_B10$`FDR q-val`<0.05)
# sum(IR_B10$`FDR q-val`<0.05)

sum(IR_B10$`FDR q-val`<0.05)
IR_B10$`FDR q-val`[IR_B10$NAME %in% NONIR_B10$NAME[NONIR_B10$`FDR q-val`<0.05]]=0.99
sum(IR_B10$`FDR q-val`<0.05)
write.table(IR_B10, file = paste(c('results/gsea/B10_DOWN.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)

NONIR_B10 = read.table('results/gsea/UP_IN_NonIR_B10_versus_NonIR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(NONIR_B10)
IR_B10 = read.table('results/gsea/UP_IN_IR_B10_versus_IR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(IR_B10)

# sum(NONIR_B10$`FDR q-val`<0.05)
# sum(IR_B10$`FDR q-val`<0.05)

sum(IR_B10$`FDR q-val`<0.05)
IR_B10$`FDR q-val`[IR_B10$NAME %in% NONIR_B10$NAME[NONIR_B10$`FDR q-val`<0.05]]=0.99
sum(IR_B10$`FDR q-val`<0.05)
write.table(IR_B10, file = paste(c('results/gsea/B10_UP.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)

do_plots2 = function(A, B, C, D, pCuttOff, gsea) {
  A_ = unlist(strsplit(A, '_versus_'))
  B_ = unlist(strsplit(B, '_versus_'))
  C_ = unlist(strsplit(C, '_versus_'))
  D_ = unlist(strsplit(D, '_versus_'))
  
  A1 = gsea[[A]][[A_[1]]]; A1 = A1[A1$`FDR q-val`<pCuttOff,]
  B1 = gsea[[B]][[B_[1]]]; B1 = B1[B1$`FDR q-val`<pCuttOff,]
  C1 = gsea[[C]][[C_[1]]]; C1 = C1[C1$`FDR q-val`<pCuttOff,]
  D1 = gsea[[D]][[D_[1]]]; DD1 = D1; D1 = D1[D1$`FDR q-val`<pCuttOff,]
  
  A2 = gsea[[A]][[A_[2]]]; A2 = A2[A2$`FDR q-val`<pCuttOff,]
  B2 = gsea[[B]][[B_[2]]]; B2 = B2[B2$`FDR q-val`<pCuttOff,]
  C2 = gsea[[C]][[C_[2]]]; C2 = C2[C2$`FDR q-val`<pCuttOff,]
  D2 = gsea[[D]][[D_[2]]]; DD2 = D2; D2 = D2[D2$`FDR q-val`<pCuttOff,]
  
  # D1ONLY = D1[((!(D1$NAME %in% A1$NAME))&(!(D1$NAME %in% B1$NAME))&(!(D1$NAME %in% C1$NAME))),]; print(dim(D1ONLY))
  # D2ONLY = D2[((!(D2$NAME %in% A2$NAME))&(!(D2$NAME %in% B2$NAME))&(!(D2$NAME %in% C2$NAME))),]; print(dim(D2ONLY))
  # DONLY = rbind.data.frame(D1ONLY, D2ONLY)
  # 
  # write.table(D1ONLY, file = paste(c('results/gsea/UP_',A,B,C,D,'.tsv'),collapse=''), sep='\t', quote=FALSE)
  # write.table(D2ONLY, file = paste(c('results/gsea/DOWN_',A,B,C,D,'.tsv'),collapse=''), sep='\t', quote=FALSE)
  # write.table(DONLY, file = paste(c('results/gsea/UP_AND_DOWN_',A,B,C,D,'.tsv'),collapse=''), sep='\t', quote=FALSE)
  
  DD1$`FDR q-val`[(DD1$NAME %in% A1$NAME)|(DD1$NAME %in% B1$NAME)|(DD1$NAME %in% C1$NAME)]=0.99; print(dim(DD1))
  DD2$`FDR q-val`[(DD2$NAME %in% A2$NAME)|(DD2$NAME %in% B2$NAME)|(DD2$NAME %in% C2$NAME)]=0.99; print(dim(DD2))
  # DD1$BETTERNAME = gsub('GOBP_', '', DD1$NAME) 
  # # DD1$BETTERNAME = gsub('_', ' ', DD1$BETTERNAME)   
  # DD2$BETTERNAME = gsub('GOBP_', '', DD2$NAME) 
  # # DD2$BETTERNAME = gsub('_', ' ', DD2$BETTERNAME)   
  
  write.table(DD1, file = paste(c('results/gsea/UP_IN_',D,'.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)
  write.table(DD2, file = paste(c('results/gsea/DOWN_IN_',D,'.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)
  
  
  # B3 = rbind(A1[!(A1$NAME %in% B1$NAME),],A2[!(A2$NAME %in% B2$NAME),])
  # B3 = rbind(B1[!(B1$NAME %in% A1$NAME),],B2[!(B2$NAME %in% A2$NAME),])
  
  png(paste(c('results/gsea/UP_',A,B,C,D,'_venn.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,4,4,4))
  x = list(A = A1$NAME, B = B1$NAME, C = C1$NAME, D = D1$NAME)
  plot(ggVennDiagram(x, set_size = 5, category.names = c(paste(c(A, '\nupregulated in ', A_[1]), collapse=''),
                                                         paste(c(B, '\nupregulated in ', B_[1]), collapse=''),
                                                         paste(c(C, '\nupregulated in ', C_[1]), collapse=''),
                                                         paste(c(D, '\nupregulated in ', D_[1]), collapse=''))))
  dev.off()

  png(paste(c('results/gsea/DOWN_',A,B,C,D,'_venn.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,4,4,4))
  x = list(A = A2$NAME, B = B2$NAME, C = C2$NAME, D = D2$NAME)
  plot(ggVennDiagram(x, set_size = 5, category.names = c(paste(c(A, '\ndownregulated in ', A_[2]), collapse=''),
                                                         paste(c(B, '\ndownregulated in ', B_[2]), collapse=''),
                                                         paste(c(C, '\ndownregulated in ', C_[2]), collapse=''),
                                                         paste(c(D, '\ndownregulated in ', D_[2]), collapse=''))))
  dev.off()
}

do_plots2(A = 'IR_A1_versus_IR_DMSO', B = 'IR_A10_versus_IR_DMSO', C = 'IR_B1_versus_IR_DMSO', D = 'IR_B10_versus_IR_DMSO', pCuttOff = 0.05, gsea=gsea)
do_plots2(A = 'NonIR_A1_versus_NonIR_DMSO', B = 'NonIR_A10_versus_NonIR_DMSO', C = 'NonIR_B1_versus_NonIR_DMSO', D = 'NonIR_B10_versus_NonIR_DMSO', pCuttOff = 0.05, gsea=gsea)


# do_plots2(A = 'IR_A1_versus_IR_DMSO', B = 'IR_B1_versus_IR_DMSO', C = 'IR_A10_versus_IR_DMSO', D = 'IR_B10_versus_IR_DMSO', pCuttOff = 0.05, gsea=gsea)
# do_plots2(A = 'NonIR_A1_versus_NonIR_DMSO', B = 'NonIR_B1_versus_NonIR_DMSO', C = 'NonIR_A10_versus_NonIR_DMSO', D = 'NonIR_B10_versus_NonIR_DMSO', pCuttOff = 0.05, gsea=gsea)

names(gsea)
do_plots = function(A, B, pCuttOff, gsea) {

  A_ = unlist(strsplit(A, '_versus_'))
  B_ = unlist(strsplit(B, '_versus_'))
  A1 = gsea[[A]][[A_[1]]]; A1 = A1[A1$FDR.q.val<pCuttOff,]
  B1 = gsea[[B]][[B_[1]]]; B1 = B1[B1$FDR.q.val<pCuttOff,]
  A2 = gsea[[A]][[A_[2]]]; A2 = A2[A2$FDR.q.val<pCuttOff,]
  B2 = gsea[[B]][[B_[2]]]; B2 = B2[B2$FDR.q.val<pCuttOff,]
  
  A3 = rbind(A1[!(A1$NAME %in% B1$NAME),],A2[!(A2$NAME %in% B2$NAME),])
  B3 = rbind(B1[!(B1$NAME %in% A1$NAME),],B2[!(B2$NAME %in% A2$NAME),])
  
  png(paste(c('results/gsea/UP_',A,B,'_venn.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,4,4,4))  
  x = list(A = A1$NAME, B = B1$NAME)
  plot(ggVennDiagram(x, set_size = 5, category.names = c(paste(c(A, '\nupregulated in ', A_[1]), collapse=''), paste(c(B, '\nupregulated in ', B_[1]), collapse=''))))
  dev.off()

  png(paste(c('results/gsea/DOWN_',A,B,'_venn.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,4,4,4))  
  x = list(A = A2$NAME, B = B2$NAME)
  plot(ggVennDiagram(x, set_size = 5, category.names = c(paste(c(A, '\ndownregulated in ', A_[2]), collapse=''), paste(c(B, '\ndownregulated in ', B_[2]), collapse=''))))
  dev.off()

  
  
  print(dim(A3))
  # png(paste(c('results/gsea/',A,'.barplot.png'),collapse=''), height = 2400, width = 2400, res=100)
  png(paste(c('results/gsea/',A,'_barplot.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,50,4,4))  
  names = A3$NAME; names = gsub('GOBP_', '', names); names = gsub('_', ' ', names)
  heights = A3$NES
  barplot(height = heights, names.arg = names, las=2, horiz = T, main = A, xlim=c(-3,3), xlab='NES')
  dev.off()
  
  print(dim(B3))
  # png(paste(c('results/gsea/',B,'.barplot.png'),collapse=''), height = 2400, width = 2400, res=100)
  png(paste(c('results/gsea/',B,'_barplot.png'),collapse=''), height = 2400, width = 2400, res=100)
  par(mar=c(4,50,4,4))  
  
  names = B3$NAME; names = gsub('GOBP_', '', names); names = gsub('_', ' ', names)
  heights = B3$NES
  barplot(height = heights, names.arg = names, las=2, horiz = T, main = B, xlim=c(-3,3), xlab='NES')
  dev.off()
}
# dev.off()  

do_plots(A = 'IR_A1_versus_IR_DMSO', B = 'IR_B1_versus_IR_DMSO', pCuttOff = 0.05, gsea=gsea)
do_plots(A = 'IR_A10_versus_IR_DMSO', B = 'IR_B10_versus_IR_DMSO', pCuttOff = 0.05, gsea=gsea)




for (condition in c('post', 'followup')) {
  for (group in c('COPD', 'Control')) {
    NR_ = paste(c(condition, group, 'NR', 'versus', 'baseline', group), collapse = '_')
    Placebo_ = paste(c(condition, group, 'Placebo', 'versus', 'baseline', group), collapse = '_')
    NR_UP_DOWN = print(unlist(strsplit(NR_, '_versus_')))
    Placebo_UP_DOWN = print(unlist(strsplit(Placebo_, '_versus_')))
    
    NR_label = paste(c(NR_, '\nupregulated in', NR_UP_DOWN[1]), collapse = ' ')
    Placebo_label = paste(c(Placebo_, '\nupregulated in', Placebo_UP_DOWN[1]), collapse = ' ')
    
    NR = gsea[[NR_]][[NR_UP_DOWN[1]]]
    NR = NR$NAME[NR$FDR.q.val<cuttoff]
    Placebo = gsea[[Placebo_]][[Placebo_UP_DOWN[1]]]
    Placebo = Placebo$NAME[Placebo$FDR.q.val<cuttoff]    
    x = list(
      NR = NR,
      Placebo = Placebo
    )
    png(paste(c('results/today/gsea/UP_', NR_, '_VS_', Placebo_, '.png'), collapse = ''), height = 1200, width = 1200)  
    print(ggVennDiagram(x, category.names = c(NR_label,Placebo_label), set_size = 6))
    dev.off()
    
    A=gsea[[NR_]][[NR_UP_DOWN[1]]]
    B=gsea[[Placebo_]][[Placebo_UP_DOWN[1]]]
    A=A[A$FDR.q.val<0.05,]; dim(A)
    B=B[B$FDR.q.val<0.05,]; dim(B)
    C = A[!(A$NAME %in% B$NAME),]
    write.table(C, file = paste(c('results/today/gsea/UP_', NR_, '.tsv'), collapse = ''), sep = '\t', quote = F)
    C = B[!(B$NAME %in% A$NAME),]
    write.table(C, file = paste(c('results/today/gsea/UP_', Placebo_, '.tsv'), collapse = ''), sep = '\t', quote = F)    
    
    
    NR_label = paste(c(NR_, '\ndownregulated in', NR_UP_DOWN[2]), collapse = ' ')
    Placebo_label = paste(c(Placebo_, '\ndownregulated in', Placebo_UP_DOWN[2]), collapse = ' ')
    NR = gsea[[NR_]][[NR_UP_DOWN[2]]]
    NR = NR$NAME[NR$FDR.q.val<cuttoff]
    Placebo = gsea[[Placebo_]][[Placebo_UP_DOWN[2]]]
    Placebo = Placebo$NAME[Placebo$FDR.q.val<cuttoff]    
    
    x = list(
      NR = NR,
      Placebo = Placebo
    )
    png(paste(c('results/today/gsea/DOWN_', NR_, '_VS_', Placebo_, '.png'), collapse = ''), height = 1200, width = 1200)  
    print(ggVennDiagram(x, category.names = c(NR_label,Placebo_label), set_size = 6))
    dev.off()    
    # gsea[[NR]][[NR_UP_DOWN[1]]]
    # gsea[[NR]][[NR_UP_DOWN[2]]]
    # gsea[[Placebo]][[Placebo_UP_DOWN[1]]]
    # gsea[[Placebo]][[Placebo_UP_DOWN[2]]]    
    
    # png(paste(c('results/today/gsea/',NR_,'.png'), collapse = ''), height = 2400, width = 1200, res=100)
    # par(mar=c(4,50,4,4))
    # barplot(height = c(C$NES,-C$NES), names.arg = c(C$NAME,C$NAME), las=2, horiz = T)
    # dev.off()
    
    A=gsea[[NR_]][[NR_UP_DOWN[2]]]
    B=gsea[[Placebo_]][[Placebo_UP_DOWN[2]]]
    A=A[A$FDR.q.val<0.05,]; dim(A)
    B=B[B$FDR.q.val<0.05,]; dim(B)
    C = A[!(A$NAME %in% B$NAME),]
    write.table(C, file = paste(c('results/today/gsea/DOWN_', NR_, '.tsv'), collapse = ''), sep = '\t', quote = F)
    C = B[!(B$NAME %in% A$NAME),]
    write.table(C, file = paste(c('results/today/gsea/DOWN_', Placebo_, '.tsv'), collapse = ''), sep = '\t', quote = F)    
    
    
  }
}








# grepl('1', condition)] = 1; dose

gsea_results = list.dirs('results/gsea', recursive=F, full.names=T); gsea_results

gsea = list()
for (comparison in comparisons) {
  gsea[comparison] = comparison
}
gsea

# read.table(gsea_results[grepl(comparisons[1], gsea_results)]

IR_DMSO_versus_NonIR_DMSO = read.table('results/gsea/IR_DMSO_versus_Control_DMSO.Gsea.1670756417660/gsea_report_for_IR_DMSO_1670756417660.tsv', sep='\t', h=T); dim(IR_DMSO_versus_Control_DMSO)
Control_B10_versus_Control_DMSO = read.table('/root/results/Control_B10_versus_Control_DMSO.Gsea.1670756752029/gsea_report_for_Control_B10_1670756752029.tsv', sep='\t', h=T); dim(Control_B10_versus_Control_DMSO)
IR_B10_versus_IR_DMSO = read.table('/root/results/IR_B10_versus_IR_DMSO.Gsea.1670755377035/gsea_report_for_IR_B10_1670755377035.tsv', sep='\t', h=T); dim(IR_B10_versus_IR_DMSO)
Control_A10_versus_Control_DMSO = read.table('/root/results/Control_A10_versus_Control_DMSO.Gsea.1670757043334/gsea_report_for_Control_A10_1670757043334.tsv', sep='\t', h=T); dim(Control_A10_versus_Control_DMSO)
IR_A10_versus_IR_DMSO = read.table('/root/results/IR_A10_versus_IR_DMSO.Gsea.1670755675579/gsea_report_for_IR_A10_1670755675579.tsv', sep='\t', h=T); dim(IR_A10_versus_IR_DMSO)

pCutoff = 0.25; FCcutoff = 2;

IR_DMSO_versus_Control_DMSO = IR_DMSO_versus_Control_DMSO[IR_DMSO_versus_Control_DMSO$FDR.q.val<pCutoff,]; dim(IR_DMSO_versus_Control_DMSO)
Control_B10_versus_Control_DMSO = Control_B10_versus_Control_DMSO[Control_B10_versus_Control_DMSO$FDR.q.val<pCutoff,]; dim(Control_B10_versus_Control_DMSO)
IR_B10_versus_IR_DMSO = IR_B10_versus_IR_DMSO[IR_B10_versus_IR_DMSO$FDR.q.val<pCutoff,]; dim(IR_B10_versus_IR_DMSO)
Control_A10_versus_Control_DMSO = Control_A10_versus_Control_DMSO[Control_A10_versus_Control_DMSO$FDR.q.val<pCutoff,]; dim(Control_A10_versus_Control_DMSO)
IR_A10_versus_IR_DMSO = IR_A10_versus_IR_DMSO[IR_A10_versus_IR_DMSO$FDR.q.val<pCutoff,]; dim(IR_A10_versus_IR_DMSO)


IR_B10_versus_IR_DMSO = IR_B10_versus_IR_DMSO[!(IR_B10_versus_IR_DMSO$NAME %in% Control_B10_versus_Control_DMSO$NAME),]; dim(IR_B10_versus_IR_DMSO)
IR_B10_versus_IR_DMSO = IR_B10_versus_IR_DMSO[!(IR_B10_versus_IR_DMSO$NAME %in% IR_DMSO_versus_Control_DMSO$NAME),]; dim(IR_B10_versus_IR_DMSO)
IR_B10_versus_IR_DMSO = IR_B10_versus_IR_DMSO[order(IR_B10_versus_IR_DMSO$FDR.q.val),]
write.table(IR_B10_versus_IR_DMSO, file = '/root/results/71.tsv', sep = '\t')
head(IR_B10_versus_IR_DMSO[c('NAME','FDR.q.val' )], 20)



V = list(
  IR_DMSO_versus_NonIR_DMSO = IR_DMSO_versus_Control_DMSO$NAME,
  NonIR_B10_versus_NonIR_DMSO = Control_B10_versus_Control_DMSO$NAME,
  IR_B10_versus_IR_DMSO = IR_B10_versus_IR_DMSO$NAME
)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(V)

V = list(
  IR_DMSO_versus_NonIR_DMSO = IR_DMSO_versus_Control_DMSO$NAME,
  NonIR_A10_versus_NonIR_DMSO = Control_A10_versus_Control_DMSO$NAME,
  IR_A10_versus_IR_DMSO = IR_A10_versus_IR_DMSO$NAME
)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(V)


########



########



plotPCA(vsd[,as.vector(data$drug!='A')], intgroup=c('treatment', 'condition'))
head(data)

pcaData <- plotPCA(vsd[,as.vector(data$drug=='B')], intgroup=c('condition', 'treatment'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()








