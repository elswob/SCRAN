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
options( max.print=300, width=100 )

run_scran=function(counts,sing_cols,cpmVal=1,pc=5,spike_text="ERCC",species="Human"){
  
}