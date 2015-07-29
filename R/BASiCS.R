#' BASiCS normalisation
#' @param dCounts The raw counts including ERCC spike in data
#' @param sing_cols The column names for the single cell data
basics_norm=function(dCounts,sing_cols,outDir){
  print("Running BASiCS...")
  #make directory
  dir.create(outDir,showWarnings = F)
  #generate Counts
  Counts=as.matrix(dCounts[,sing_cols])
  #convert counts data to integers
  Counts <- apply (Counts, c (1, 2), function (x) { (as.integer(x)) })
  
  #TRUE/FALSE for spike in data
  Tech=grepl("ERCC",rownames(dCounts))
  
  #Get spike in molecular data
  if(file.exists("ercc_data.txt")){
    print("Already got ERCC data...")
    e=read.delim("ercc_data.txt")
  }else{
    x=getURL("https://tools.lifetechnologies.com/content/sfs/manuals/cms_095046.txt")
    e=read.delim(text=x)
    write.table(e,file="ercc_data.txt",sep="\t",quote=F)
  }
  #only use the spike ins left after filtering counts
  ematch=e[,2] %in% rownames(dCounts)
  SpikeInput=e[ematch,][,4]
  
  #Data = newBASiCS_Data(Counts, Tech, SpikeInput)

  #filter_data_object}
  Filter = BASiCS_Filter(Counts, Tech, SpikeInput, 
                         MinTotalCountsPerCell = 2, MinTotalCountsPerGene = 2, 
                         MinCellsWithExpression = 2, MinAvCountsPerCellsWithExpression = 2)
  Data = newBASiCS_Data(Filter$Counts, Filter$Tech, Filter$SpikeInput)
  
  #fit the model
  MCMC_Output <- BASiCS_MCMC(Data, N = 400, Thin = 10, Burn = 200, StoreChains = T, StoreDir = outDir, RunName = "400")
  #MCMC_Output <- BASiCS_MCMC(Data, N = 40000, Thin = 10, Burn = 20000, StoreChains = T, StoreDir = outDir, RunName = "40k")
  
  #check the model
  plot(MCMC_Output, Param = "mu", Gene = 1)
  plot(MCMC_Output, Param = "delta", Gene = 1)
  plot(MCMC_Output, Param = "phi", Cell = 1)
  plot(MCMC_Output, Param = "s", Cell= 1)
  plot(MCMC_Output, Param = "nu", Cell = 1)
  plot(MCMC_Output, Param = "theta")
  
  MCMC_Summary <- Summary(MCMC_Output)
  #head(displaySummaryMu(MCMC_Summary))
  #head(displaySummaryDelta(MCMC_Summary))
  #head(displaySummaryPhi(MCMC_Summary))
  #head(displaySummaryS(MCMC_Summary))
  #head(displaySummaryNu(MCMC_Summary))
  #head(displaySummaryTheta(MCMC_Summary))
  
  #The following figures display posterior medians and the corresponding HPD 95% intervals for the normalising constants. The first figure, allows the assessment of cell-to-cell heterogeneity regarding to mRNA content. The second figure, displays possible evidence of cell-to-cell differeces in capture efficiency (and/or amplification biases).
  plot(MCMC_Summary, Param = "phi")
  plot(MCMC_Summary, Param = "s")
  
  #It is also possible to generate similar plots for other model parameters
  plot(MCMC_Summary, Param = "mu", main = "All genes")
  plot(MCMC_Summary, Param = "mu", Genes = 1:10, main = "First 10 genes")
  plot(MCMC_Summary, Param = "delta", main = "All genes")
  plot(MCMC_Summary, Param = "delta", Genes = c(2,5,10,50,100), main = "5 customized genes")
  plot(MCMC_Summary, Param = "nu", main = "All cells")
  plot(MCMC_Summary, Param = "nu", Cells = 1:5, main = "First 5 cells")
  plot(MCMC_Summary, Param = "theta")
  
  #To contrasts posterior medians of cell-specific parameters use
  plot(MCMC_Summary, Param = "phi", Param2 = "s")
  plot(MCMC_Summary, Param = "phi", Param2 = "nu")
  plot(MCMC_Summary, Param = "s", Param2 = "nu")
  
  #To display posterior medians of \(\delta_i\) (the parameters controlling the strength of the biological cell-to-cell expression heterogeneity of a gene \(i\) across the population of cells under study) againts overall gene-specific expression levels \(\mu_i\) use:
  plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x")
  
  #high_low
  VarDecomp = BASiCS_VarianceDecomp(MCMC_Output)
  print(head(VarDecomp))
  
  DetectHVG <<- BASiCS_DetectHVG(MCMC_Output, VarThreshold = 0.70, Plot = TRUE)
  #DetectHVG = DetectHVG[order(DetectHVG$GeneNames),]
  HVG_out<<-DetectHVG$Table
  #HVG_out$ensembl=rownames(dCounts)[HVG_out$GeneIndex]
  write.table(HVG_out, file=paste0(outDir,"basics_HVG.tsv"),sep="\t",quote=F,col.names=NA)
  
  DetectLVG <- BASiCS_DetectLVG(MCMC_Output, VarThreshold = 0.70, Plot = TRUE)
  #DetectLVG = DetectLVG[order(DetectLVG$GeneNames),]
  LVG_out=DetectLVG$Table
  #LVG_out$ensembl=rownames(dCounts)[LVG_out$GeneIndex]
  write.table(LVG_out, file=paste0(outDir,"basics_LVG.tsv"),sep="\t",quote=F,col.names=NA)
  
  DenoisedCounts <- BASiCS_DenoisedCounts(Data = Data, Chain = MCMC_Output)
  write.table(x = DenoisedCounts, file = paste0(outDir,"basics_normalised_counts.tsv"),sep="\t",quote = F, col.names = NA)
  return(DenoisedCounts)
}