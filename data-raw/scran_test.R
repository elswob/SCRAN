scran_test=read.delim('data-raw/genes.counts.filtered.test.txt', header=T, row.names=1)
save(scran_test,file='data/scran_test.rdata', compress='xz')
