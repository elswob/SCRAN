#' Assign biotypes to each gene
#'
#' @param m A data frame with gene symbols as row names.
biotypes=function(m){
}

#' Filter the counts by CPM and percentage of samples
#'
#' @param dCounts The counts dataframe
#' @param sing_cols A vector of the single cell column names
#' @param cpmVal The minimum CPM value 
#' @param pc The percentage of cells to be expressed at the minumum CPM
#' @return The original and filtered counts
filter=function(dCounts,sing_cols,cpmVal,samp_pc){
  keep=rowSums(cpm(dCounts[,sing_cols])>cpmVal) >= (length(sing_cols)/100)*samp_pc
  summary(keep)
  dCounts=dCounts[keep,]
  dCounts_orig=dCounts
  return(list("dCounts" = dCounts, "dCounts_orig" = dCounts_orig))
}

#' Separate counts into genes and spike ins
sepCounts=function(dCounts_orig,sing_cols,spike_text){
  spikeData = dCounts_orig[grep(spike_text,rownames(dCounts_orig)),sing_cols]
  geneData = dCounts_orig[grep(spike_text,rownames(dCounts_orig),invert = T),sing_cols]
  geneCounts = colSums(geneData)
  spikeCounts=0
  if(nrow(spikeData)>0){
    spikeCounts = colSums(spikeData)
  }
  return(list("spikeData"=spikeData,"geneData"=geneData,"spikeCounts"=spikeCounts,"geneCounts"=geneCounts))
}

#' Plot the raw counts
plot_raw_counts=function(geneCounts,erccCounts,outDir){
  count_df=as.data.frame(geneCounts)
  if(sum(erccCounts)>0){
    count_df$erccCounts=erccCounts
    count_df$epc = (count_df$erccCounts/(count_df$geneCounts+count_df$erccCounts))*100
    count_df$gpc = 100-count_df$epc
  }
  summary(count_df)
  
  #plot
  count_df$sample=rownames(count_df)
  mv="geneCounts"
  if(sum(erccCounts)>0){
    mv=c("geneCounts","erccCounts")
  }
  m=melt(count_df,id.vars="sample",measure.vars=mv)
  pdf(paste0(outDir,"read_counts_per_sample.pdf"))
  g=ggplot(m, aes(x=sample,y=value,fill=variable)) + geom_bar(stat="identity") 
  g = g + xlab("Sample") + ylab("Read count") + theme(axis.text.x = element_text(size=5,angle = 45, hjust = 1))
  print(g)
  dev.off()
  print(g)
  write.table(m,file=paste0(outDir,"raw_counts_per_sample.tsv"),sep="\t",quote=F,row.names=F)
  
  if(sum(erccCounts)>0){  
    pdf(paste0(outDir,"read_pcs_per_sample.pdf"))
    m=melt(count_df,id.vars="sample",measure.vars=c("epc","gpc"))
    m=m[m$value>0 & m$value<101,]
    g=ggplot(m, aes(x=sample,y=value,fill=variable)) + geom_bar(stat="identity")
    g = g + xlab("Sample") + ylab("Percentage") + theme(axis.text.x = element_text(size=5, angle = 45, hjust = 1))
    print(g)
    dev.off()
    print(g)
    write.table(m,file=paste0(outDir,"raw_counts_pc_per_sample.tsv"),sep="\t",quote=F,row.names=F)
  }
}

#' Add fake bulk sample column created from mean of all single cells
fake_bulk=function(dCounts){
  Mean_SC=rowSums(dCounts[,sing_cols])/length(sing_cols)
  dCounts=cbind(dCounts,Mean_SC)
  #sing_cols = c(sing_cols,"Mean_SC")
  return(dCounts)
}

#' Read distributions
read_dist=function(dCounts,sing_cols,outDir){
  #how many genes identified at various CPM cutoffs
  l=list(c(0,1),c(1,5),c(5,10),c(10,50),c(50,100),c(100,500),c(500,1000),c(1000,5000),c(5000,10000),c(10000,10000000000))
  for(i in l){
    print(i[1])
    k=rowSums(dCounts[,sing_cols]>=i[1] & dCounts[,sing_cols]<i[2]) >= fNum
    print(table(k))
  }
  
  #how many reads mapped to genes ordered by most expressed genes
  r=rowSums(dCounts[,sing_cols])
  r_names=names(sort(r,decreasing = T))
  #print(r_names)
  mean_cols=mean(colSums(dCounts[,sing_cols]))
  total=sum(colSums(dCounts[,sing_cols]))
  sep=10
  counts=c()
  divs=c()
  m<<-0
  for (i in seq(1,length(r_names),by=sep)){
    if(i+sep<=length(r_names)){
      #remove one to stop overlap
      j=i+sep-1
      #print(paste0('i = ',i,' j = ',j))
      r_sub=r_names[i:j]
      #divide number of reads per gene (or genes) by total mean number of reads for each sample and convert to %
      #m=mean(rowSums(dCounts[r_sub,sing_cols]))
      m=m+sum(rowSums(dCounts[r_sub,sing_cols]))
      g=(m/total)*100
      #print(m)
      counts=c(counts,g)
      divs=c(divs,j)
    }
  }
  #print(divs)
  #print(counts)
  df=as.data.frame(divs)
  df$counts=counts
  df$divs=log10(df$divs)
  head(df)
  
  #plot
  #g=ggplot(df, aes(x=divs,y=counts)) + geom_bar(stat="identity")
  pdf(paste0(outDir,"cumulative_counts_per_gene.pdf"))
  g=ggplot(df,aes(x=divs,y=counts)) + geom_line() + xlab("Number of genes (log10)") + ylab("Cumulative percentage of counts")
  g=g + scale_y_continuous(breaks=seq(0, 100, 10))  + expand_limits(y=0) 
  g=g + scale_x_continuous(breaks=c(1,2,3,4,5))
  #g=g + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
  print(g)
  dev.off()
  
  #make file
  a=rowSums(dCounts[r_names,sing_cols])
  m=merge(a,dCounts_orig[,1,drop=FALSE],by=0)
  colnames(m)[1]="Ensembl"
  colnames(m)[2]="Count"
  colnames(m)[3]="Symbol"
  m=m[order(-m$Count),]
  m$PC=(m$Count/total)*100
  #reorder
  m=m[c("Ensembl","Symbol","Count","PC")]
  head(m)
  write.table(m,file=paste0(outDir,"sorted_counts_per_genes.tsv"), quote=F, row.names=F, sep="\t")
  
  #plot distribution of top 10
  t10=r_names[1:10]
  t10_counts=dCounts[t10,sing_cols]
  t10_pc=(t10_counts/colSums(dCounts[,sing_cols]))*100
  t10_data=merge(t10_pc,dCounts_orig[,1,drop=FALSE],by=0)
  colnames(t10_data)[ncol(t10_data)]="Symbol"
  tm=melt(t10_data,measure.vars=sing_cols,id.vars="Symbol")
  #tm$value=log10(tm$value)
  g=ggplot(tm, aes(x = reorder(Symbol, -value, FUN=median), y=value)) + geom_boxplot()
  g=g+ xlab("Top 10 expressed genes") + ylab("Percentage of counts")
  g=g+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(paste0(outDir,"top_10_expressed_gene_distribution.pdf"))
  print(g)
  dev.off()
  print(g)
}