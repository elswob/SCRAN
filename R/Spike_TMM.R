#' TMM normalisation using spike in size factors
tmm_norm=function(geneData,spikeData,outDir){
  dir.create(outDir,showWarnings = F)
  sfSpike <<- estimateSizeFactorsForMatrix(spikeData)
  nCountsGenes <<- t( t(geneData) / sfSpike )
  write.table(x = nCountsGenes, file = paste0(outDir,"TMM_spike_normalised_counts.tsv"),sep="\t",quote = F, col.names = NA)
  return(nCountsGenes)
}