library(pzfx) #install.packages('pzfx')
sapply(list.files(pattern='pzfx$'), function(file){
  sapply(pzfx_tables(file), function(table){
    write.table(read_pzfx(file, table=table), file = paste(c(file, table, 'tsv'), collapse = '.'), sep='\t')
  })
})