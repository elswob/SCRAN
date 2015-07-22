#' Test function
run_test=function(outDir){
  #load the test data
  print("Loading data...")
  #load("scran_test.Rdata")
  print(dim(scran_test))
  d=scran_test
  
  #set single cell columns and get column names
  sing_cols=c(3:ncol(d))
  sing_cols=sing_col_convert(d,sing_cols)
  
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
  
  #set out dirs
  dir.create(outDir,showWarnings = F)
  rawDir=paste0(outDir,"/Raw/")
  dir.create(rawDir,showWarnings = F)
  qcDir=paste0(rawDir,"/QC/")
  dir.create(qcDir,showWarnings = F)
  #plot raw counts
  print("Raw counts...")
  plot_raw_counts(sep$geneCounts,sep$spikeCounts,qcDir)
  
  #Read dist
  print("Read dist...")
  read_dist(sep$geneData,sing_cols,qcDir,dFilt)
  
  #counts per gene
  print("Counts per gene...")
  cpg(sep$geneData,sing_cols,qcDir)
  
  #pca
  print("PCA...")
  pca_heatmap(sep$geneData,sing_cols,10,qcDir)
  
  #unique genes
  print("Unique genes per cell...")
  #u=uniq_genes(sep$geneData,sing_cols)
  #print(u)
  
  #gene counts
  print("Gene counts per sample...")
  gc_per_samp(sep$geneData,qcDir)
  
  #biotypes
  print("Biotypes...")
  biotypes(sep$geneData,sing_cols,"Mouse",qcDir)
  
  #ercc plots
  print("ERCC plots...")
  spike_in_check(sep$spikeCounts, sep$spikeData, sing_cols, qcDir)
  
  #spike in and hkg
  print("Spike-ins and HKGs...")
  spike_hkg(sep$geneData,sep$spikeData,"Mouse",sing_cols,qcDir)
  
  #siber
  siber_res=siber(sep$geneData,sing_cols,rawDir,dFilt)
  
  #brennecke analysis
  print("Brennecke analysis...")
  #brennecke(dCounts = d, species = "Mouse", outDir = outDir, spike_text = "ERCC", pc = 10)
  
  #basics normalisation
  print("BASiCS normalisation...")
  baDir=paste0(outDir,"/BASiCS/")
  #b_norm<<-basics_norm(dFilt,sing_cols,baDir)
  #run QC and analysis steps on BASiCS
  #pca_heatmap(b_norm,sing_cols,top=50,outDir=baDir)
  #b_sep=sepCounts(b_norm,sing_cols,"ERCC")
  #spike_hkg(geneData = as.data.frame(b_sep$geneData),spikeData = as.data.frame(b_sep$spikeData), species = "Mouse", outDir = baDir)
  #siber_res=siber(b_sep$geneData,sing_cols,baDir,dFilt)
  
  #normalise using TMM and spike ins
  tDir=paste0(outDir,"/TMM_spike/")
  t_norm=tmm_norm(geneData = sep$geneData,spikeData = sep$spikeData, outDir=tDir)
  #run QC and analysis
  #pca_heatmap(t_norm,sing_cols,top=50,tDir)
  #spike_hkg(geneData = as.data.frame(t_norm), spikeData = sep$spikeData, species = "Mouse", sing_cols = sing_cols, outDir = tDir)
  #siber_res=siber(t_norm,sing_cols,tDir,dFilt)
}