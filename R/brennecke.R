#Modified from here - Brennecke et al, Accounting for technical noise in single-cell RNA-seq experiments, Nature Methods 2013
#http://www.nature.com/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf

#' Requires df of gene symbols, lengths and counts, rownames =  ensembl IDs and colnames = single cell IDs
brennecke=function(dCounts,species,outDir,spike_text,pc){
  #str(d)
  print(dim(d))
  #keep genes with counts in X% of cells and remove symbol column
  keep=rowSums(d[,3:ncol(d)]>0) >= (ncol(d)/100)*pc
  dCounts<-d[keep,2:ncol(d)]
  print(dim(dCounts))
  #print(head(dCounts))
  
  #add one to every count as estimateSizeFactorsForMatrix fails if no gene has counts across all samples
  #this might be a bit of a hack
  dCounts[,2:ncol(dCounts)]=dCounts[,2:ncol(dCounts)]+1
  
  #significane threshold
  sigVal=0.1
  
  outDir=paste0(outDir,"/Brennecke/")
  dir.create(outDir,showWarnings = F)
  
  #dCounts[ 1:10, 1:5 ]
  geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[ substr( rownames(dCounts), 1, 4 ) ] )
  print(table(geneTypes))
  spikeData <- dCounts[grep(spike_text,rownames(dCounts)),-1]
  geneData <- dCounts[grep(spike_text,rownames(dCounts),invert = T),-1]
  spikeLengths <- dCounts[grep(spike_text,rownames(dCounts)),1] 
  geneLengths <- dCounts[grep(spike_text,rownames(dCounts),invert = T), 1]
  #countsHsap <- dCounts[ which( geneTypes=="ENSG" ), -1 ]
  #countsERCC <- dCounts[ which( geneTypes=="ERCC" ), -1 ]
  #lengthsHsap <- dCounts[ which( geneTypes=="ENSG" ), 1 ]
  #lengthsERCC <- dCounts[ which( geneTypes=="ERCC" ), 1 ]
  print(head(geneData))
  
  #Calculate size factors
  sfHsap <- estimateSizeFactorsForMatrix( geneData )
  print(sfHsap)
  sfERCC <- estimateSizeFactorsForMatrix( spikeData )
  print(sfHsap)
  #print(rbind( sfHsap, sfERCC ))
  
  #test to show how that works
  mdat <- matrix(c(1,2,3, 1,2,3, 1,2,3 ), nrow = 3, ncol = 3, byrow = TRUE,
                 dimnames = list(c("row1", "row2", "row3"), c("C.1", "C.2", "C.3")))
  mdat
  e=estimateSizeFactorsForMatrix(mdat)
  e 
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) 
  }
  g=gm_mean(c(1,2,3))
  em=median(c(1/g,1/g,1/g))
  em
  
  #normalise by them
  nspikeData <- t( t(spikeData) / sfERCC )
  ngeneData <- t( t(geneData) / sfHsap )
  write.table(nspikeData,file=paste0(outDir,"/nspikeData.txt"),sep="\t",quote=F,col.names = NA)
  write.table(ngeneData,file=paste0(outDir,"/nCountsGenes.txt"),sep="\t",quote=F, col.names = NA)
  
  #We calculate the sample moments:
  meansERCC <- rowMeans( nspikeData )
  varsERCC <- rowVars( nspikeData )
  cv2ERCC <- varsERCC / meansERCC^2
  meansHsap <- rowMeans( ngeneData )
  varsHsap <- rowVars( ngeneData )
  cv2Hsap <- varsHsap / meansHsap^2
  
  #Normalize the mean counts by transcript length (i.e., "per kilobase", "PK"), too:
  meansERCCPK <- meansERCC / spikeLengths * 1e3
  #print(head(meansERCCPK))
  meansHsapPK <- meansHsap / geneLengths * 1e3
  
  #2.3 Plot of normalised counts
  colHsap <- "#00207040"
  colERCC <- "#70500040"
  
  #pairs( log10( .1 + rbind( ngeneData, nspikeData ) ), pch=19, cex=.2, col = c( rep( colHsap, nrow(ngeneData) ), rep( colERCC, nrow(nspikeData) ) ) )
  
  geneScatterplot <- function( x, y, xlab, ylab, col ) {
    plot( NULL, xlim=c( -.1, 6.2 ), ylim=c( -1, 6.2 ),
          xaxt="n", yaxt="n", xaxs="i", yaxs="i", asp=1,
          xlab=xlab, ylab=ylab )
    abline( a=-1, b=1, col = "lightgray", lwd=2 )
    abline( a=0, b=1, col = "lightgray", lwd=2 )
    abline( a=1, b=1, col = "lightgray", lwd=2 )
    abline( h=c(0,2,4,6), v=c(0,2,4,6), col = "lightgray", lwd=2 )
    points(
      ifelse( x > 0, log10(x), -.7 ),
      ifelse( y > 0, log10(y), -.7 ),
      pch=19, cex=.2, col = col )
    axis( 1, c( -.7, 0:6 ),
          c( "0", "1", "10", "100", expression(10^3), expression(10^4),
             expression(10^5), expression(10^6) ) )
    axis( 2, c( -.7, 0:6 ),
          c( "0", "1", "10", "100", expression(10^3), expression(10^4),
             expression(10^5), expression(10^6) ), las=2 )
    axis( 1, -.35, "//", tick=FALSE, line=-.7 )
    axis( 2, -.35, "\\\\", tick=FALSE, line=-.7 )
  }
  
  #geneScatterplot( geneData[,1], geneData[,3], "un-normalized read count, cell 1", "un-normalized read count, cell 3", colHsap )
  geneScatterplot( ngeneData[,1], ngeneData[,3], "normalized read count, cell 1", "normalized read count, cell 3", colHsap )
  geneScatterplot( nspikeData[,1], nspikeData[,3], "normalized read count, cell 1", "normalized read count, cell 3", colERCC )
  
  #4.2 Fit technical noise
  #We perform the fit as usual. However, as we have only rather few spikes, we have to be a bit more generous with the mean cut-off, now using the 80-percentile instead of the 95-percentile.
  minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
  useForFitA <- meansERCC >= minMeanForFitA
  print(paste0('minMeanForFitA=',minMeanForFitA))
  table( useForFitA )
  
  #Afterwards, we will compare with a fit using length-normalized counts. We prepare by finding the minimum for these, too:
  minMeanForFitB <- unname( quantile( meansERCCPK[ which( cv2ERCC > .3 ) ], .8 ) )
  useForFitB <- meansERCCPK >= minMeanForFitB
  print(paste0('minMeanForFitB=',minMeanForFitB))
  print(table( A=useForFitA, B=useForFitB ))
  
  #Note that the two lists overlap well.
  #We perform both fits.
  fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                      cv2ERCC[useForFitA] )
  fitB <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCCPK[useForFitB] ),
                      cv2ERCC[useForFitB] )
  
  #How much variance do the two fits explain?
  residualA <- var( log( fitted.values(fitA) ) - log( cv2ERCC[useForFitA] ) )
  totalA <- var( log( cv2ERCC[useForFitA] ) )
  residualB <- var( log( fitted.values(fitB) ) - log( cv2ERCC[useForFitB] ) )
  totalB <- var( log( cv2ERCC[useForFitB] ) )
  # explained variances of log CV^2 values
  c( A = 1 - residualA / totalA, B = 1 - residualB / totalB )
  
  #Fit B, which used the length-normalized counts, performed better.
  #As a second check, we plot both fits.
  par( mfrow=c(1,2) )
  plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA, main="A" )
  xg <- 10^seq( -3, 5, length.out=100 )
  lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
  segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
            meansERCC[useForFitA], fitA$fitted.values, col="gray" )
  
  plot( meansERCCPK, cv2ERCC, log="xy", col=1+useForFitB, main="B" )
  lines( xg, coefficients(fitB)["a0"] + coefficients(fitB)["a1tilde"]/xg )
  segments( meansERCCPK[useForFitB], cv2ERCC[useForFitB],
            meansERCCPK[useForFitB], fitB$fitted.values, col="gray" )
  
  #4.3 Test for High Variance
  #We start with the test using fit A:
  minBiolDisp <- .5^2
  xi <- mean( 1 / sfERCC )
  m <- ncol(geneData)
  psia1thetaA <- mean( 1 / sfERCC ) +
    ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfHsap )
  cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
  testDenomA <- ( meansHsap * psia1thetaA + meansHsap^2 * cv2thA ) / ( 1 + cv2thA/m )
  pA <- 1 - pchisq( varsHsap * (m-1) / testDenomA, m-1 )
  padjA <- p.adjust( pA, "BH" )
  table( padjA < sigVal )
  
  #Using fit B and the length-normalized counts, we get
  varsHsapPK <- rowVars( ngeneData / geneLengths * 1e3 )
  psia1thetaB <- mean( 1 / sfERCC ) +
    ( coefficients(fitB)["a1tilde"] - xi ) * mean( sfERCC / sfHsap )
  cv2thB <- coefficients(fitB)["a0"] + minBiolDisp + coefficients(fitB)["a0"] * minBiolDisp
  testDenomB <- ( meansHsapPK * psia1thetaB + meansHsapPK^2 * cv2thB ) / ( 1 + cv2thB/m )
  pB <- 1 - pchisq( varsHsapPK * (m-1) / testDenomB, m-1 )
  padjB <- p.adjust( pB, "BH" )
  table( B = padjB < sigVal )
  
  #4.4 A diagnostic plot
  #We have a closer look at the reliability of the technical noise predictions from fit A. The variance from technical noise,
  #predicted for a biological gene is given by Omega( 0, mu ), where mu is the normalized mean count for the gene, and
  #Omega is the function defined in the Online Methods. Dividing by mu2 to get CV2 values (and ignoring the negligible term
  # a0/m), we get
  
  predictedNoiseCV2 <- psia1thetaA / meansHsap + coefficients(fitA)["a0"]
  
  #We plot the ratio of observed total CV2 to predicted technical CV2 against transcript length, using only genes with a mean
  #count above the cut-off also used for the fit. This is Supplementary Figure 7.
  
  useInPlot <- meansHsap>minMeanForFitA
  plot( geneLengths[useInPlot], ( cv2Hsap / predictedNoiseCV2 )[useInPlot], log="xy",
        pch=20, cex=.2, col = "#705000A0", xlab = "length", ylab="total variance / predicted noise" )
  points( spikeLengths[useForFitA], cv2ERCC[useForFitA] / fitted.values(fitA), pch=20, cex=1)
  
  #4.5 Plot of results
  #plot fitA results
  #need to get a0 from 2.4
  #To get the actual noise coefficients, we need to subtract Xi (see Supplementary Note 6 for the difference between a1tilde and a1).
  xi <- mean( 1 / sfHsap )
  a0 <- unname( fitA$coefficients["a0"] )
  a1 <- unname( fitA$coefficients["a1tilde"] - xi )
  c( a0, a1 )
  
  pdf(paste0(outDir,"var_vs_counts_normalised.pdf"))
  #To produce Figure 3, which depicts the results of fit and test A, the following code was used
  plot( NULL, xaxt="n", yaxt="n", log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
        xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                         expression(10^4), expression(10^5) ) )
  axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  # Plot the plant genes, use a different color if they are highly variable
  points( meansHsap, cv2Hsap, pch=20, cex=.2,
          col = ifelse( padjA < sigVal, "#C0007090", "#70500040" ) )
  # Add the technical noise fit, as before
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, coefficients(fitA)["a1tilde"] / xg + a0, col="#FF000080", lwd=3 )
  # Add a curve showing the expectation for the chosen biological CV^2 thershold
  lines( xg, psia1thetaA/xg + coefficients(fitA)["a0"] + minBiolDisp,
         lty="dashed", col="#C0007090", lwd=3 )
  # Add the normalised ERCC points
  points( meansERCC, cv2ERCC, pch=20, cex=1, col="#0060B8A0" )
  dev.off()
  
  #plot fitB results
  a0 <- unname( fitB$coefficients["a0"] )
  a1 <- unname( fitB$coefficients["a1tilde"] - xi )
  c( a0, a1 )
  pdf(paste0(outDir,"var_vs_counts_length_normalised.pdf"))
  #To produce Figure 3, which depicts the results of fit and test A, the following code was used
  plot( NULL, xaxt="n", yaxt="n", log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
        xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                         expression(10^4), expression(10^5) ) )
  axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  # Plot the plant genes, use a different color if they are highly variable
  points( meansHsap, cv2Hsap, pch=20, cex=.2,
          col = ifelse( padjB < sigVal, "#C0007090", "#70500040" ) )
  # Add the technical noise fit, as before
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, coefficients(fitB)["a1tilde"] / xg + a0, col="#FF000080", lwd=3 )
  # Add a curve showing the expectation for the chosen biological CV^2 thershold
  lines( xg, psia1thetaB/xg + coefficients(fitB)["a0"] + minBiolDisp,
         lty="dashed", col="#C0007090", lwd=3 )
  # Add the normalised ERCC points
  points( meansERCCPK, cv2ERCC, pch=20, cex=1, col="#0060B8A0" )
  dev.off()
  
  #2.5.2 Table of higly variable genes
  log2RelExprAt <- log2( ngeneData / meansHsap )
  print(dim(log2RelExprAt))
  
  sig <- padjB < sigVal
  print(table(sig))
  
  print(log2RelExprAt[sig,])
  
  print(colnames(log2RelExprAt))
  #strongest = factor( colnames( log2RelExprAt )[ apply( log2RelExprAt[ sig, ], 1, which.max ) ] )
  a=apply( log2RelExprAt[ sig, ], 1, which.max )
  print(paste0('a = ',a))
  
  highVarTable <- data.frame(
    row.names = NULL,
    geneID = rownames(geneData)[ sig ],
    #geneSymbol = geneSymbols[ sig ],
    meanNormCount = meansHsap[ sig ],
    strongest = factor( colnames( log2RelExprAt )[
      apply( log2RelExprAt[ sig, ], 1, which.max ) ] ),
    log2RelExprAt[ sig, ],
    check.names=FALSE )
  
  d$ens=rownames(d)
  d_symbols=d[,c(1,ncol(d))]
  m=merge(highVarTable,d_symbols,by.x="geneID",by.y="ens")
  
  highVarTable=m
  colnames(highVarTable)[ncol(highVarTable)]="geneSymbol"
  
  head( highVarTable ) 
  write.table( highVarTable, file=paste0(outDir,"/length_adjusted_highly_variant_genes.tsv"), row.names=FALSE, sep="\t", quote=F )
  
  #2.6 GO analysis
  minCountForEnrichment <- 300
  
  topGOAnalysis <- function( geneIDs, inUniverse, inSelection )
    sapply( c( "MF", "BP", "CC" ), function( ont ) {
      alg <<- factor( as.integer( inSelection[inUniverse] ) )
      names(alg) <<- geneIDs[inUniverse]
      if(toupper(species) == 'HUMAN'){
        tgd <<- new( "topGOdata", ontology=ont, allGenes = alg, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl" )
      }else{
        tgd <<- new( "topGOdata", ontology=ont, allGenes = alg, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl" )
      }
      resultTopGO <<- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
      #plot go structures
      pdf(paste0(outDir,"sig_go_tree_",ont,".pdf"))
      showSigOfNodes(tgd, score(resultTopGO), firstSigNodes = 5, useInfo = 'all')
      dev.off()
      GenTable( tgd, resultTopGO, topNodes=20 )
    },
    simplify=FALSE )
  
  goResults <-
    topGOAnalysis(
      rownames(geneData),
      meansHsap >= minCountForEnrichment & !is.na(padjB),
      padjB < sigVal )
  
  #write tables
  write.table(goResults$MF, file=paste0(outDir,"topGO_MF_table.tsv"), sep="\t", quote=F, col.names=T, row.names = F)
  write.table(goResults$BP, file=paste0(outDir,"topGO_BP_table.tsv"), sep="\t", quote=F, col.names=T, row.names = F)
  write.table(goResults$CC, file=paste0(outDir,"topGO_CC_table.tsv"), sep="\t", quote=F, col.names=T, row.names = F)
  
  #2.7 Heatmap
  pdf(paste0(outDir,"high_variant_clustering.pdf"),onefile = F)
  relSig=log2RelExprAt[sig,]
  #relSig[relSig=='-Inf']=0
  relSig[ relSig < -4 ] <- -4
  pheatmap(relSig,show_rownames=F)
  dev.off()
}