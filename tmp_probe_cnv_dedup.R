suppressPackageStartupMessages(library(TCGAbiolinks))
q <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number"
)
res <- q$results[[1]]
cat('rows', nrow(res), '\n')
cat('cols', paste(colnames(res), collapse='|'), '\n')
show_cols <- intersect(c('file_id','file_name','cases','sample.submitter_id','sample_type','workflow.type','experimental_strategy','platform','data_format'), colnames(res))
print(utils::head(res[, show_cols, drop=FALSE], 12))
cat('\nDUP_COUNTS\n')
for (k in c('cases','sample.submitter_id','workflow.type','platform','file_name')) {
  if (k %in% colnames(res)) {
    vals <- res[[k]]
    cat(k, 'unique=', length(unique(vals)), 'dups=', sum(duplicated(vals)), '\n')
  }
}
if ('cases' %in% colnames(res) && 'workflow.type' %in% colnames(res)) {
  cat('\nTOP_CASE_WORKFLOW_DUP\n')
  tab <- sort(table(paste(res$cases, res$workflow.type, sep=' || ')), decreasing=TRUE)
  print(utils::head(tab[tab > 1], 20))
}
if ('sample.submitter_id' %in% colnames(res) && 'workflow.type' %in% colnames(res)) {
  cat('\nTOP_SAMPLE_WORKFLOW_DUP\n')
  tab <- sort(table(paste(res[['sample.submitter_id']], res[['workflow.type']], sep=' || ')), decreasing=TRUE)
  print(utils::head(tab[tab > 1], 20))
}
