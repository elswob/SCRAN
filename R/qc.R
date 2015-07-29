#' Convert single cell column numbers to names
sing_col_convert=function(dCounts,sing_cols){
  return(colnames(dCounts)[sing_cols])
}

#' Filter the counts by CPM and percentage of samples
#'
#' @param dCounts The counts dataframe
#' @param sing_cols A vector of the single cell column names
#' @param cpmVal The minimum CPM value 
#' @param pc The percentage of cells to be expressed at the minumum CPM
#' @return The original and filtered counts
filter_counts=function(dCounts,sing_cols,cpmVal,samp_pc){
  keep=rowSums(cpm(dCounts[,sing_cols])>cpmVal) >= (length(sing_cols)/100)*samp_pc
  summary(keep)
  return(dCounts[keep,])
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
plot_raw_counts=function(geneCounts,spikeCounts,outDir){
  count_df=as.data.frame(geneCounts)
  if(sum(spikeCounts)>0){
    count_df$spikeCounts=spikeCounts
    count_df$spc = (count_df$spikeCounts/(count_df$geneCounts+count_df$spikeCounts))*100
    count_df$gpc = 100-count_df$spc
  }
  summary(count_df)
  
  #plot
  count_df$sample=rownames(count_df)
  mv="geneCounts"
  if(sum(spikeCounts)>0){
    mv=c("geneCounts","spikeCounts")
  }
  m=melt(count_df,id.vars="sample",measure.vars=mv)
  pdf(paste0(outDir,"Read_counts_per_sample.pdf"))
  g=ggplot(m, aes(x=sample,y=value,fill=variable)) + geom_bar(stat="identity") 
  g = g + xlab("Sample") + ylab("Read count") + theme(axis.text.x = element_text(size=5,angle = 45, hjust = 1))
  print(g)
  dev.off()
  print(g)
  write.table(m,file=paste0(outDir,"Raw_counts_per_sample.tsv"),sep="\t",quote=F,row.names=F)
  
  if(sum(spikeCounts)>0){  
    pdf(paste0(outDir,"Read_pcs_per_sample.pdf"))
    m=melt(count_df,id.vars="sample",measure.vars=c("spc","gpc"))
    m=m[m$value>0 & m$value<101,]
    g=ggplot(m, aes(x=sample,y=value,fill=variable)) + geom_bar(stat="identity")
    g = g + xlab("Sample") + ylab("Percentage") + theme(axis.text.x = element_text(size=5, angle = 45, hjust = 1))
    print(g)
    dev.off()
    print(g)
    write.table(m,file=paste0(outDir,"Raw_counts_pc_per_sample.tsv"),sep="\t",quote=F,row.names=F)
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
read_dist=function(dCounts,sing_cols,outDir,dCounts_orig){
  #how many reads mapped to genes ordered by most expressed genes
  r=rowSums(dCounts[,sing_cols])
  r_names=names(sort(r,decreasing = T))
  #print(r_names)
  mean_cols=mean(colSums(dCounts[,sing_cols]))
  total=sum(colSums(dCounts[,sing_cols]))
  sep=10
  counts=c()
  divs=c()
  m<-0
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
  pdf(paste0(outDir,"Cumulative_counts_per_gene.pdf"))
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
  write.table(m,file=paste0(outDir,"Sorted_counts_per_genes.tsv"), quote=F, row.names=F, sep="\t")
  
  #plot distribution of top 10
  t10=r_names[1:10]
  t10_counts=dCounts[t10,sing_cols]
  t10_pc=(t10_counts/colSums(dCounts_orig[,sing_cols]))*100
  t10_data=merge(t10_pc,dCounts_orig[,1,drop=FALSE],by=0)
  colnames(t10_data)[ncol(t10_data)]="Symbol"
  tm=melt(t10_data,measure.vars=sing_cols,id.vars="Symbol")
  #tm$value=log10(tm$value)
  g=ggplot(tm, aes(x = reorder(Symbol, -value, FUN=median), y=value)) + geom_boxplot()
  g=g+ xlab("Top 10 expressed genes") + ylab("Percentage of counts")
  g=g+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(paste0(outDir,"Top_10_expressed_gene_distribution.pdf"))
  print(g)
  dev.off()
  print(g)
}

#' Counts per gene per cell
cpg=function(dCounts,sing_cols,outDir){
  if(length(sing_cols)<50){
    #plot each cells separately maybe?
    r=rowSums(dCounts[,sing_cols])
    r_names=names(sort(r,decreasing = T))
    #print(r_names)
    total=sum(colSums(dCounts[,sing_cols]))
    sep=10
    counts=c()
    divs=c()
    #cells=c()
    m<<-0
    for (i in seq(0,length(r_names),by=sep)){
      #print(i)
      if(i+sep<=length(r_names)){
        j=i+sep
        #add one to stop overlap at boundaries
        r_sub=r_names[(i+1):j]
        #divide number of reads per gene (or genes) by total mean number of reads for each sample and convert to %
        #m=mean(rowSums(dCounts[r_sub,sing_cols]))
        #m=m+sum(rowSums(dCounts[r_sub,sing_cols]))
        m=m+colSums(dCounts[r_sub,sing_cols])
        pcm=(m/colSums(dCounts[,sing_cols])*100)
        #print(m)
        counts=c(counts,pcm)
        divs=c(divs,rep(j,length(sing_cols)))
        #cells=c(cells,colnames(dCounts[sing_cols]))
      }
    }
    #print(divs)
    #print(counts)
    #print(cells)
    df=as.data.frame(divs)
    df$counts=counts
    df$cells=names(counts)
    df$divs=log10(df$divs)
    head(df)
    
    pdf(paste0(outDir,"Cumulative_counts_per_gene_per_cell.pdf"))
    g=ggplot(df,aes(x=divs,y=counts,color=cells)) + geom_point(size = 1) + xlab("Number of genes (log10)") + ylab("Cumulative percentage of counts")
    g=g + scale_y_continuous(breaks=seq(0, 100, 10))  + expand_limits(y=0) 
    g=g + scale_x_continuous(breaks=c(1,2,3,4,5))
    if (length(sing_cols)>20){
      g=g+guides(col=guide_legend(ncol=2))
    }
    print(g)
    dev.off()
  }else{
    print("Too many samples for plot!")
  }
}

#' PCA and heatmap
pca_heatmap=function(geneCounts,sing_cols,top,outDir,allCounts){
  #tpm transform
  tpmCounts=countToTpm(geneCounts,allCounts$Length)
  p.pca=prcomp(t(tpmCounts))
  #p.pca=prcomp(t(cpm(geneCounts)))
  pdf(paste0(outDir,"PCA.pdf"))
  g = ggbiplot(p.pca, obs.scale = 0, ellipse = TRUE, varname.size=0.001, labels=colnames(geneCounts))
  print(g)
  dev.off()
  pdf(paste0(outDir,"PCA_scree_plot.pdf"))
  screeplot(p.pca,type="lines",col=3)
  dev.off()
  
  #top=50
  load.rot=p.pca$rotation
  top_pca=names(load.rot[,1][order(abs(load.rot[,1]),decreasing=TRUE)][1:top])
  top_pca
  
  #redo PCA
  tpmCounts=tpmCounts[rownames(tpmCounts) %in% top_pca,]
  p.pca=prcomp(t(tpmCounts))
  #p.pca=prcomp(t(cpm(geneCounts)))
  pdf(paste0(outDir,"PCA_top_",top,".pdf"))
  g = ggbiplot(p.pca, obs.scale = 0, ellipse = TRUE, varname.size=0.001, labels=colnames(geneCounts))
  print(g)
  dev.off()
  pdf(paste0(outDir,"PCA_scree_plot_top_",top,".pdf"))
  screeplot(p.pca,type="lines",col=3)
  dev.off()
  
  #convert to log 
  tpmCounts[tpmCounts==0]=0.01
  g_log2=log2(tpmCounts[,sing_cols])
  #print(g_log10[0:5,0:5])
  #get relevant data
  dc=g_log2[rownames(g_log2) %in% top_pca,]
  
  #cluster and heatmap
  pdf(paste0(outDir,"Heatmap_PCA_top_",top,"_heatmap.pdf"),onefile=FALSE)
  if(nrow(dc)<=100){
    pheatmap(dc, show_rownames = T, fontsize_row=7, main="Unsupervised clustering of log2 TPM")
  }else{
    pheatmap(dc, show_rownames = F)
  }
  dev.off()
}

#' Convert counts to TPM
countToTpm <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

#' Number of unique genes expressed per sample
uniq_genes=function(geneCounts,sing_cols){
  geneCounts=geneCounts[,sing_cols]
  
  uc=function(rd){
    #find cases where number of cpms > limit for a gene is 1, i.e. unique for one sample 
    if(sum(rd)==1){
      #identify which sample
      which(rd>0)
    }
  }
  #run the function above on the data
  a=apply(geneCounts,1,uc) 
  #convert freq counts to df
  b=as.data.frame(table(unlist(a)))
  #add sample names to df
  b$Cell=colnames(geneCounts)[c(as.numeric(levels(b$Var1)))]
  b$Var1=NULL
  return(b)
  #mean(b$Freq)
}

#' Gene counts per sample
gc_per_samp=function(geneCounts,outDir){
  
  #calculate number of genes per sample greater than min CPM
  cc=colSums(geneCounts)
  
  #calculate mean of single-cells
  m=mean(cc)
  #cat("mean = ",m)
  
  #plot
  ccm<<-melt(cc)
  #print(head(ccm))
  pdf(paste0(outDir,"Expressed_genes_per_sample.pdf"))
  g=ggplot(ccm, aes(y=value,x=rownames(ccm))) + geom_bar(stat="identity")
  g=g+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Number of genes") + xlab("Samples") 
  g=g+geom_hline(yintercept=m, colour="red")
  g=g+scale_y_continuous(breaks=seq(0, 16000, 1000)) 
  print(g)
  dev.off()
}

#' Biotypes
biotypes=function(geneCounts,sing_cols,species,outDir){
  print(Sys.time())
  if(toupper(species)=='HUMAN'){
    print("Getting human biomart data...")
    ensembl <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  }else{
    print("Getting mouse biomart data...")
    ensembl <<- useEnsembl(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
  }
  btypes<-getBM(attributes=c('gene_biotype','ensembl_gene_id'), filters = "ensembl_gene_id", values=rownames(geneCounts), mart=ensembl)
  print(Sys.time())
  #print(nrow(btypes))
  #print(geneCounts[0:5,0:5])
  if(nrow(btypes)>0){
    dkb=merge(geneCounts,btypes[1:2],by.x=0,by.y="ensembl_gene_id")
    #print(dkb[0:5,0:5])
    print(dim(dkb))
    
    #relable less common biotypes as 'other' (default 1% across all samples)
    min_biotype=1
    a=as.data.frame(table(dkb$gene_biotype))
    a$pc=a$Freq/sum(a$Freq)*100
    #print(a)
    dkb=merge(dkb,a,by.x="gene_biotype",by.y="Var1")
    dkb[dkb$pc<min_biotype,]$gene_biotype="Other"
    dkb$gene_biotype=capitalize(gsub("_"," ",dkb$gene_biotype))
    print(table(dkb$gene_biotype))
    
    #melt the data and plot
    dkm<<-melt(dkb,measure.vars = sing_cols,id.vars="gene_biotype")
    dkm=dkm[dkm$value>0,]
    #print(head(dkm))
    pdf(paste0(outDir,"Biotypes_and_gene_counts.pdf"))
    g=ggplot(dkm, aes(x=variable,fill=gene_biotype)) + geom_bar(position="fill") + scale_y_continuous(labels = percent_format())
    g = g + xlab("Sample") + ylab("Percentage")+ scale_fill_discrete(name="Biotype")  + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5))
    print(g)
    dev.off()
    print(g)
    
    pdf(paste0(outDir,"Biotypes_percentage_per_cell.pdf"))
    g=ggplot(dkm, aes(x=variable,fill=gene_biotype)) + geom_bar() 
    g = g + xlab("Sample") + ylab("Count")+ scale_fill_discrete(name="Biotype")  + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5))
    g=g+geom_hline(yintercept=m, colour="red")
    g=g+scale_y_continuous(breaks=seq(0, 20000, 1000))
    print(g)
    dev.off()
    print(g)
  }else{
    print(paste0("No biotype data found, perhaps ",species," is the wrong species?"))
  }
}

#' Spike in plots
spike_in_check=function(spikeCounts,spikeData,sing_cols,outDir){
  if(file.exists("ercc_data.txt")){
    print("Already got ERCC data...")
    e=read.delim("ercc_data.txt")
  }else{
    x=getURL("https://tools.lifetechnologies.com/content/sfs/manuals/cms_095046.txt")
    e=read.delim(text=x)
    write.table(e,file="ercc_data.txt",sep="\t",quote=F)
  }
  m=merge(spikeData,e[,c(2,4)],by.x=0,by.y="ERCC.ID")
  #print(m[0:5,0:5])
  if(nrow(m)>0){
    mm=melt(m,measure.vars=sing_cols,id.vars=c(1,ncol(m)))
    colnames(mm)[2]="conc"
    #add value to zeroes to avoid log errors
    mm[mm==0]=0.01
    #maybe only include values above 1
    #print(head(mm))
    pdf(paste0(outDir,"ERCC_plot_log_log.pdf"))
    g=""
    if(length(sing_cols)<11){
      #print(length(sing_cols))
      g = ggplot(mm,aes(x=conc,y=value, color=variable)) + geom_point(shape=1) 
    }else{
      g = ggplot(mm,aes(x=conc,y=value, group=variable)) + geom_point(shape=1)
    }
    g = g + scale_y_log10() + scale_x_log10()
    g = g + geom_smooth(method=lm, se=F) 
    print(g)
    dev.off()
    print(g)
    
    g= ""
    pdf(paste0(outDir,"ERCC_plot_log_y.pdf"))
    if (length(sing_cols)<11){
      g = ggplot(mm,aes(x=conc,y=value, color=variable)) + geom_line()
    }else{
      g = ggplot(mm,aes(x=conc,y=value, group=variable)) + geom_line()
    }
    g = g + scale_y_log10()
    #g = g + geom_smooth(method=lm, se=F, )
    print(g)
    dev.off()
    print(g)
  }else{
    print("There were no ERCC spike ins found")
  }
}
#' Spike ins and house keeper genes
spike_hkg=function(geneData,spikeData,species,sing_cols,outDir){
  #gapdh
  gapdh="ENSG00000111640"
  actb="ENSG00000075624"
  b2m="ENSG00000166710"
  if(toupper(species)=='MOUSE'){
    gapdh="ENSMUSG00000057666"
    actb="ENSMUSG00000029580"
    b2m="ENSMUSG00000060802"
  }
  #geneData[,sing_cols]=cpm(geneData[,sing_cols])
  gapdhCounts=geneData[gapdh,sing_cols]
  actbCounts=geneData[actb,sing_cols]
  b2mCounts=geneData[b2m,sing_cols]
  
  #gapdh and actb
  #print(head(actbCounts))
  
  #merge
  m=melt(actbCounts,measure.vars=c(1:ncol(actbCounts)))
  rownames(m)=m[,1]
  m[,1]=NULL
  m=cbind(m,t(gapdhCounts),t(b2mCounts))
  colnames(m)[1:3]=c('ACTB','GAPDH','B2M')
  if(nrow(spikeData)>0){
    m$Spike=colMeans(spikeData)
  }
  m$samples=rownames(m)
  mm=melt(m,id.vars="samples")
  
  pdf(paste0(outDir,"HKG_spike_in_plot.pdf"))
  g=ggplot(mm,aes(x=samples,y=value,color=variable, group=variable))+geom_line()
  g=g+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) + ylab("Count") + xlab("Samples")
  g=g+geom_line(stat = "hline", yintercept = "mean", aes(colour = variable), linetype="dashed")
  print(g)
  dev.off()
  print(g)
}