library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H')
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))

gcSample = split(all_gene_sets$gene_symbol,
                 all_gene_sets$gs_name)
names(gcSample)
file="Homo-H-examp.gmt"
gs=gcSample
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat( paste(c(i,'tmp',gs[[i]]),collapse='\t') )
    cat('\n')
  })
  sink()
}
write.gmt(gs,file)
