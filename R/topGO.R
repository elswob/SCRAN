#' TopGO Analysis
topGOAnalysis <- function( geneIDs, inUniverse, inSelection ,species, outDir){
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
}