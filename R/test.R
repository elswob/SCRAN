#d=read.delim("test/genes.counts.filtered.txt",header=T,row.names = 1)
sing_cols=c(3:ncol(d))
sing_cols=sing_col_convert(d,sing_cols)
print(dim(d))

#filter
dFilt=filter_counts(d,sing_cols)
print(dim(dFilt))

#split
sep=sepCounts(dFilt,sing_cols, "ERCC")
print(dim(sep$spikeData))
print(dim(sep$geneData))

#set out dir
outDir="test/"

#plot raw counts
plot_raw_counts(sep$geneCounts,sep$spikeCounts,outDir)

print("Read dist")
read_dist(sep$geneData,sing_cols,outDir,dFilt)