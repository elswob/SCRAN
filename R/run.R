library(pheatmap)
library(ggbiplot)
library(reshape)
library(ggplot2)
library(edgeR)
library(biomaRt)
library(R.utils)
library(RCurl)
library(DESeq)
library(devtools)
library(genefilter)
library(EBImage)
library(statmod)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Rgraphviz)
library(SIBER)
library(edgeR)

scran_setup=function(){
  #filter
  print("Filter...")
  dFilt=filter_counts(d,sing_cols,cpmVal = 5)
  print(dim(dFilt))
  print(summary(rowSums(dFilt[,sing_cols])))
  
  #split
  print("Data split...")
  sep=sepCounts(dFilt,sing_cols, "ERCC")
  print(dim(sep$spikeData))
  print(dim(sep$geneData))
}

scran_qc=function(geneData, spikeData, geneCounts, spikeCounts, qcDir, sing_cols, dFilt, species){
  #plot raw counts
  print("Raw counts...")
  plot_raw_counts(geneCounts,spikeCounts,qcDir)
  
  #Read dist
  print("Read dist...")
  read_dist(geneData,sing_cols,qcDir,dFilt)
  
  #counts per gene
  print("Counts per gene...")
  cpg(geneData,sing_cols,qcDir)
  
  #gene counts
  print("Gene counts per sample...")
  gc_per_samp(geneData,qcDir)
  
  #biotypes
  print("Biotypes...")
  biotypes(geneData,sing_cols,species,qcDir)
  
  #ercc plots
  print("ERCC plots...")
  spike_in_check(spikeCounts, spikeData, sing_cols, qcDir)
}

scran_brennecke=function(){
  #brennecke analysis
  print("Brennecke analysis...")
  brDir=paste0(outDir,"/Brennecke/")
  brennecke(dCounts = d, species = "Mouse", outDir = brDir, spike_text = "ERCC", samp_pc = samp_pc)
}

scran_analysis=function(geneData, spikeData, sing_cols, top, outDir, species, counts){
  #PCA
  pca_heatmap(geneData,sing_cols,top,outDir)
  #spike-ins and HKGs
  spike_hkg(geneData = as.data.frame(geneData),spikeData = as.data.frame(spikeData), species = species, outDir = outDir)
  #siber
  siber_res=siber(geneData,sing_cols,outDir,counts)
}

#run_scran=function(counts,sing_cols,cpmVal=1,pc=5,spike_text="ERCC",species="Human",outDir){
scran_run=function(counts,sing_cols,outDir,cpmVal=1,samp_pc=5,top=50,spike_text="ERCC",species="Human"){
  ### Raw data ###
  #set out dirs
  dir.create(outDir,showWarnings = F)
  rawDir=paste0(outDir,"/Raw/")
  dir.create(rawDir,showWarnings = F)
  qcDir=paste0(rawDir,"/QC/")
  dir.create(qcDir,showWarnings = F)
  
  #convert sing_cols
  sing_cols=sing_col_convert(counts,sing_cols)
  
  #filte the counts
  print("Filter...")
  dFilt=filter_counts(counts,sing_cols,cpmVal,samp_pc)
  print(dim(dFilt))
  print(summary(rowSums(dFilt[,sing_cols])))
  
  #split into genes and spike-ins
  print("Data split...")
  sep=sepCounts(dFilt,sing_cols, spike_text)
  print(dim(sep$spikeData))
  print(dim(sep$geneData))
  
  #QC
  #scran_qc(geneData = sep$geneData, spikeData = sep$spikeData, geneCounts = sep$geneCounts, spikeCounts = sep$spikeCounts, qcDir = qcDir, sing_cols = sing_cols, dFilt = dFilt, species = species)
  #analysis
  aDir=paste0(rawDir,"/Analysis/")
  dir.create(aDir,showWarnings = F)
  #scran_analysis(geneData = sep$geneData, spikeData = sep$spikeData, species = species, counts = dFilt, sing_cols = sing_cols, top = top, outDir = aDir)
  
  ### Bennecke ###
  #remove symbol column
  d<-dFilt[,2:ncol(dFilt)]
  #brennecke(dCounts = d, species = species, outDir = outDir, spike_text = spike_text)
  
  ### TMM ###
  tDir=paste0(outDir,"/TMM_spike/")
  #t_norm=tmm_norm(geneData = sep$geneData,spikeData = sep$spikeData, outDir=tDir)
  #scran_analysis(geneData = t_norm, spikeData = sep$spikeData, sing_cols = sing_cols, top = top, outDir = tDir, species = species, counts = dFilt)
  
  ### BASiCS
  print("BASiCS normalisation...")
  baDir=paste0(outDir,"/BASiCS/")
  b_norm<-basics_norm(dFilt,sing_cols,baDir)
  scran_analysis(geneData = b_norm, spikeData = sep$spikeData, sing_cols = sing_cols, top=top, outDir = baDir, species = species, counts = dFilt)
}

test_run=function(){
  sing_cols=c(3:ncol(scran_test))
  scran_run(counts = scran_test, sing_cols = sing_cols, outDir = "~/Work/projects/single-cell_testing/scran/", species = "Mouse", samp_pc = 10)
}