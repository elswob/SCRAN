# SCRAN
Single Cell Rna-seq ANalysis

### Installation

Install libraries

```
install.packages("ggplot2")
install.packages("pheatmap")
install.packges("devtools")
install.packages("R.utils")
install.packages("RCurl")

install_github("vqv/ggbiplot")

source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("biomaRt")
biocList("DESeq")
```

```
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
library( genefilter )
library( EBImage )
library( statmod )
library( topGO )
library( org.Hs.eg.db )
library( org.Mm.eg.db )
library(Rgraphviz)
options( max.print=300, width=100 )
```

Install and load SCRAN
```
install_github("elswob/SCRAN")
libary(SCRAN)
```