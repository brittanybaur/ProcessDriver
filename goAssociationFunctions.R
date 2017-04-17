# GO-term association step 
# Associate cis genes with GO terms through the trans genes with the topGO package
# Cluster with the GOSim package
# TopGO and GOSim available through bioconductor
# cisExp - matrix of cis gene expression (rows genes, cols samples)
# transExp - matrix of trans gene expression (rows genes, cols samples)
# The matrices cisExp and transExp must have the entrez IDs for the corresponding genes as the row names
# There are two functions in this file:
# goAssociation - which takes in the cis and trans gene expression and returns a table [goTable]
# of GO terms that are associated with cis genes through the trans genes 
# and goClustering which takes in the goTable and returns a clustering of the unique terms [clusterGO] 

library('topGO')
library('GOSim')

#Example usage
#goTable = goAssociation(cisExp, transExp)
#clusterGO = goClustering(goTable)

goAssociation <- function(cisExp, transExp) {
# go Table stores the significant terms associated with each cis gene
goTable = data.frame(cisGene = character(), pathway = character(), 
annotation = character(),numberGenes = double(), pval=double())

#correlation matrix
X = cor(t(transExp),t(cisExp))

#for every cis gene perform KS test
for (y in 1:ncol(X)) {

ex1 = sort(abs(X[,y]), decreasing=TRUE)
c1 = colnames(X)[y]


GOdata <- new("topGOdata", ontology = "BP", allGenes = ex1, geneSel = function(p) p > 0.3, description = "Test", 
annot = annFUN.org, mapping="org.Hs.eg.db", ID="entrez")


#kolmogorov-smirnov
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder = "decreasing")


allRes <- GenTable(GOdata, KS = resultKS,orderBy = "KS", ranksOf = "classic", topNodes = 10)
idx2 = which(as.numeric(allRes[,"KS"]) < .05)
goID <- allRes[idx2, "GO.ID"]
pv = allRes[idx2, "KS"]
term = allRes[idx2, "Term"]
ta = allRes[idx2, "Annotated"]

# if there are any significant terms associated with the cis gene, add to the table
if (length(idx2) > 0 ){ 
for (r in 1:length(idx2)) {
goTable <- rbind(goTable, data.frame(cisGene=c1,  pathway = goID[r], annotation=term[r],
numberGenes=as.numeric(ta[r]), pval=as.numeric(pv[r]) ))
} #add a row for each significant term of the cis gene
} #if significant terms

} #KS test for every cis gene


return(goTable)
}


 
goClustering <- function(goTable) {
clust=getTermSim(unique(as.character(goTable$pathway)), method = "relevance", verbose = FALSE)
h=hclust(as.dist(1-clust),"average")
clusterGO <- cutree(h, h=max(h$height/2))
return(clusterGO)
}
