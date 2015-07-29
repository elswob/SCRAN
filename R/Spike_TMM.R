#' TMM normalisation using spike in size factors
tmm_norm=function(geneData,spikeData,outDir){
  print("TMM normalising...")
  dir.create(outDir,showWarnings = F)
  #add 1 to each count as estimateSizeFactorsForMatrix() fails if there are no genes with counts across all samples
  geneData=geneData+1
  sfGene<-estimateSizeFactorsForMatrix(geneData)
  sfSpike <- estimateSizeFactorsForMatrix(spikeData)
  normFactor<-sfGene/sfSpike
  #nCountsGenes <<- t( t(geneData) * sfSpike )
  nCountsGenes <<- t(t(geneData) * normFactor)
  #nCountsGenes <- t(t(geneData) / sfGene)
  #nCountsGenes <- nCountsGenes / sfSpike
  
  write.table(x = nCountsGenes, file = paste0(outDir,"TMM_spike_normalised_counts.tsv"),sep="\t",quote = F, col.names = NA)
  return(nCountsGenes)
}