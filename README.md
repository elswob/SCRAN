# SCRAN
Single Cell Rna-seq ANalysis

### Installation

Install libraries

```
#CRAN
install.packages("ggplot2")
install.packages("pheatmap")
install.packges("devtools")
install.packages("R.utils")
install.packages("RCurl")
install.packages("statmod")

#GitHub
install_github("vqv/ggbiplot")

#Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("biomaRt")
biocLite("DESeq")
biocLite("genefilter")
biocLite("EBImage")
biocLite("topGO")
biocLite("org.Hs.eg.db")
biocLite("org.Mm.eg.db")
biocLite("Rgraphviz")
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
library(genefilter)
library(EBImage)
library(statmod)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Rgraphviz)
```

Install and load BASiCS normalisation
```
install_github('catavallejos/BASiCS')
library(BASiCS)
```

Install and load SCRAN
```
install_github("elswob/SCRAN")
libary(SCRAN)
```