#' Test function
run_test=function(outDir){
  #load the test data
  print("Loading data...")
  data(test_data)
  print(dim(d))
  
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
  
  #set out dir
  dir.create(outDir,showWarnings = F)
  
  rawDir=paste0(outDir,"/Raw/")
  dir.create(rawDir,showWarnings = F)
  #plot raw counts
  print("Raw counts...")
  plot_raw_counts(sep$geneCounts,sep$spikeCounts,rawDir)
  
  #Read dist
  print("Read dist...")
  read_dist(sep$geneData,sing_cols,rawDir,dFilt)
  
  #counts per gene
  print("Counts per gene...")
  cpg(sep$geneData,sing_cols,rawDir)
  
  #pca
  print("PCA...")
  pca_heatmap(sep$geneData,sing_cols,10,rawDir)
  
  #unique genes
  print("Unique genes per cell...")
  #u=uniq_genes(sep$geneData,sing_cols)
  #print(u)
  
  #gene counts
  print("Gene counts per sample...")
  gc_per_samp(sep$geneData,rawDir)
  
  #biotypes
  print("Biotypes...")
  biotypes(sep$geneData,sing_cols,"Mouse",rawDir)
  
  #ercc plots
  print("ERCC plots...")
  spike_in_check(sep$spikeCounts, sep$spikeData, sing_cols, rawDir)
  
  #spike in and hkg
  print("Spike-ins and HKGs...")
  spike_hkg(sep$geneData,sep$spikeData,"Mouse",sing_cols,rawDir)
  
  #brennecke analysis
  print("Brennecke analysis...")
  brDir=paste0(outDir,"/Brennecke/")
  brennecke(dCounts = d, species = "Mouse", outDir = brDir, spike_text = "ERCC", pc = 10)
  
  #basics normalisation
  print("BASiCS normalisation...")
  baDir=paste0(outDir,"/BASiCS/")
  b_norm<<-basics_norm(dFilt,sing_cols,baDir)
  
  #run QC steps on BASiCS
  pca_heatmap(b_norm,sing_cols,top=50,outDir=baDir)
  
  #run hkg check
  sep=sepCounts(b_norm,sing_cols,"ERCC")
  print(dim(sep$geneData))
  print(dim(sep$spikeData))
  spike_hkg(geneData = sep$geneData,spikeData = sep$spikeData, species = "Mouse", outDir = baDir)
}