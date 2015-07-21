#' Test function
run_test=function(){
  d=read.delim("test/genes.counts.filtered.test.txt",header=T,row.names = 1)
  sing_cols=c(3:ncol(d))
  sing_cols=sing_col_convert(d,sing_cols)
  print("Raw...")
  print(dim(d))
  
  #filter
  print("Filter...")
  dFilt=filter_counts(d,sing_cols)
  print(dim(dFilt))
  
  #split
  print("Data split...")
  sep=sepCounts(dFilt,sing_cols, "ERCC")
  print(dim(sep$spikeData))
  print(dim(sep$geneData))
  
  #set out dir
  outDir="test/raw/"
  dir.create(outDir,showWarnings = F)
  
  #plot raw counts
  print("Raw counts...")
  plot_raw_counts(sep$geneCounts,sep$spikeCounts,outDir)
  
  #Read dist
  print("Read dist...")
  read_dist(sep$geneData,sing_cols,outDir,dFilt)
  
  #counts per gene
  print("Counts per gene...")
  cpg(sep$geneData,sing_cols,outDir)
  
  #pca
  print("PCA...")
  pca_heatmap(sep$geneData,sing_cols,10,outDir)
  
  #unique genes
  print("Unique genes per cell...")
  u=uniq_genes(sep$geneData,sing_cols)
  print(u)
  
  #gene counts
  print("Gene counts per sample...")
  gc_per_samp(sep$geneData,outDir)
  
  #biotypes
  print("Biotypes...")
  biotypes(sep$geneData,sing_cols,"Mouse",outDir)
  
  #ercc plots
  print("ERCC plots...")
  spike_in_check(sep$spikeCounts, sep$spikeData, sing_cols, outDir)
  
  #spike in and hkg
  print("Spike-ins and HKGs...")
  spike_hkg(sep$geneData,sep$spikeData,"Mouse",sing_cols,outDir)
}