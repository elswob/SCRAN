# SCRAN
Single Cell Rna-seq ANalysis

Essentially a collection of single-cell related analysis tools. Including:

- Basic QC
- Brennecke analysis
- SIBER bimodal analysis
- BASiCS normalisation.

Currently works with human and mouse data only.

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
biocLite("edgeR")

#other
source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
oompainstall(groupName="siber")
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
library(edgeR)
library(SIBER)
```

Install and load BASiCS normalisation
```
install_github('catavallejos/BASiCS')
library(BASiCS)
```

Install and load SCRAN:

From Github:
```
install_github("elswob/SCRAN")
libary(SCRAN)
```

From Stash:

1. Download using the 'Download' link
2. Unzip
3. Install

```
setwd("~/Downloads/scran_master")
install(".")
```

### Run

Requires counts data with symbol and length columns labelled 'Symbol' and 'Length' loaded with rownames as ensembl IDs and cell names as headers. An example data set is included:

```
scran_test[0:5,0:5]
                   Symbol  Length SRR1033853 SRR1033854 SRR1033855
ENSMUSG00000000001  Gnai3 3262.00        386         24          1
ENSMUSG00000000028  Cdc45 1574.00          0          0          0
ENSMUSG00000000031    H19 1268.60          0          0          0
ENSMUSG00000000037  Scml2 3297.14          0          0          0
ENSMUSG00000000056   Narf 1785.00          0         96          0
```

To test with the preloaded data just call the test_run() function with a directory to place the output, e.g.
```
test_run("/path/to/somewhere")
```

On other data
```
a=read.delim("file.tsv",header=T,row.names=1)
sing_cols=c(3:ncol(a))
scran_run(counts=a, sing_cols=sing_cols, outDir="~/", species="mouse")
```

The demo data is the first 20 cells from the Treutlin et al data set http://www.nature.com/nature/journal/v509/n7500/abs/nature13173.html