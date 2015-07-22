#' SIBER bimodal analysis
siber=function(geneData,sing_cols,outDir,dCounts){
  print("Running SIBER...")
  outDir=paste0(outDir,"/SIBER/")
  dir.create(outDir,showWarnings = F)
  counts=geneData[,sing_cols]
  print(Sys.time())
  Dat=as.matrix(counts)
  TMM <- calcNormFactors(Dat, method='TMM')
  res=apply(Dat,1,function(x) SIBER(y=x, d=1/TMM, model='LN'))
  print(Sys.time())
  
  #reformat output
  res=as.data.frame(t(res))
  res_sort=res[order(-res$BI),]
  
  #top bimodal
  top_b=head(res_sort,n=50)
  
  #get original data
  r=rownames(geneData) %in% rownames(top_b)
  bi_geneData=geneData[r,]
  bi_geneData[bi_geneData==0]=1

  #plot clustering useing top bimodal genes
  pdf(paste0(outDir,"siber_bimodal_plot.pdf"),onefile=FALSE)
  pheatmap(log10(cpm(bi_geneData)),fontsize_row = 5,fontsize_col = 5)
  dev.off()
  
  #write table
  m=merge(top_b,dCounts[,1:2],by=0)
  m=m[order(-m$BI),]
  names(m)[1]="ensembl_id"
  write.table(m,file=paste0(outDir,"siber_top_bimodal_genes.tsv"),sep="\t",quote=F,row.names = F)
  write.table(res_sort,file=paste0(outDir,"siber_all_bimodal_data.tsv"),sep="\t",quote=F,row.names = F)
}