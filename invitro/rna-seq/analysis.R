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
  save(res, file = paste(c(comparison_name, '.RData'), collapse = ''))
}



sapply(comparison_names, FUN=function(x){
  load(file = paste(c(x, '.RData'), collapse = ''))
  png(paste(c(x, '_pvalue.png'), collapse = ''))  
  plot(EnhancedVolcano(res,
                       lab = rownames(res),
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       title = x,
                       x = 'log2FoldChange',
                       y = 'pvalue'))
  
  dev.off()
  png(paste(c(x, '_padj.png'), collapse = ''))  
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

expression_gct = 'expression.gct'
df = cbind.data.frame(data.frame(Description=rep('na',dim(dds)[1])), as.data.frame(counts(dds, normalized = TRUE)))
df = cbind.data.frame(data.frame(NAME = row.names(df)), df)
cat(paste('#1.2',paste(dim(dds),collapse='\t'),sep='\n'),file=paste(expression_gct,sep=''),sep="\n")
write_tsv(df,file=paste(expression_gct,sep=''),col_names=T,append=T)

phenotype_cls = function(phenotype_name, sample_phenotypes) {
  write_file(
    paste(paste(length(sample_phenotypes),length(unique(sample_phenotypes)),'1',sep=' '),
          paste('#',paste(unique(sample_phenotypes),collapse=' '),sep=' '),
          paste(sample_phenotypes,collapse=' '),sep='\n'),
    paste(phenotype_name,'.cls',sep='')
  )
}

phenotype_cls('dose', data$dose)
phenotype_cls('group', data$group)
phenotype_cls('condition', data$condition)
phenotype_cls('treatment', data$treatment)
phenotype_cls('drug', data$drug)


gsea_cli = function(phenotype) {
  cat(
    paste(c('sh gsea-cli.sh GSEA -res expression.gct -cls group.cls#', phenotype,' -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label ', phenotype, ' -metric Signal2Noise -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Gene_Symbol_with_Remapping_MSigDB.v2022.1.Hs.chip -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out results'), collapse = ''),
    file='gsea.sh', sep='\n', append = T)
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



# GSEA Analysis

gsea = list()
gsea_results = list.dirs('gsea', recursive=F, full.names=T); gsea_results
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

do_plots = function(A, B, C, D, pCuttOff, gsea) {
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
  
  DD1$`FDR q-val`[(DD1$NAME %in% A1$NAME)|(DD1$NAME %in% B1$NAME)|(DD1$NAME %in% C1$NAME)]=0.99; print(dim(DD1))
  DD2$`FDR q-val`[(DD2$NAME %in% A2$NAME)|(DD2$NAME %in% B2$NAME)|(DD2$NAME %in% C2$NAME)]=0.99; print(dim(DD2))

  write.table(DD1, file = paste(c('UP_IN_',D,'.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)
  write.table(DD2, file = paste(c('DOWN_IN_',D,'.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)
}

do_plots(A = 'IR_A1_versus_IR_DMSO', B = 'IR_A10_versus_IR_DMSO', C = 'IR_B1_versus_IR_DMSO', D = 'IR_B10_versus_IR_DMSO', pCuttOff = 0.05, gsea=gsea)
do_plots(A = 'NonIR_A1_versus_NonIR_DMSO', B = 'NonIR_A10_versus_NonIR_DMSO', C = 'NonIR_B1_versus_NonIR_DMSO', D = 'NonIR_B10_versus_NonIR_DMSO', pCuttOff = 0.05, gsea=gsea)


NONIR_B10 = read.table('DOWN_IN_NonIR_B10_versus_NonIR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(NONIR_B10)
IR_B10 = read.table('DOWN_IN_IR_B10_versus_IR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(IR_B10)

sum(IR_B10$`FDR q-val`<0.05)
IR_B10$`FDR q-val`[IR_B10$NAME %in% NONIR_B10$NAME[NONIR_B10$`FDR q-val`<0.05]]=0.99
sum(IR_B10$`FDR q-val`<0.05)
write.table(IR_B10, file = paste(c('B10_DOWN.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)

NONIR_B10 = read.table('UP_IN_NonIR_B10_versus_NonIR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(NONIR_B10)
IR_B10 = read.table('UP_IN_IR_B10_versus_IR_DMSO.tsv', sep='\t', h=T, check.names = F); dim(IR_B10)

sum(IR_B10$`FDR q-val`<0.05)
IR_B10$`FDR q-val`[IR_B10$NAME %in% NONIR_B10$NAME[NONIR_B10$`FDR q-val`<0.05]]=0.99
sum(IR_B10$`FDR q-val`<0.05)
write.table(IR_B10, file = paste(c('B10_UP.tsv'),collapse=''), sep='\t', quote=FALSE, row.names = F, col.names = T)



