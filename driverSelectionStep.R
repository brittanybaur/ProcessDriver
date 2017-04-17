library('PMA')
library('glmnet')
library('topGO')
library('org.Hs.eg.db')
source('MultiTaskFunctions.R')

# This script can be run to create the moduleTable, which is the final
# result showing the modules and the ranking of the cis genes as drivers.
# This script assumes you have already run the 2 functions in 
# goAssociationFunctions.R and have the goTable and the clusterGO objects
# descripted in that file. 
# Also requires the matrices of cis gene expression and trans gene expression
# described in goAssociationFunctions.R

moduleTable = data.frame(goTerms = character(), goAnnotation=character(), canonicalCorrelation=numeric(),
cisGenesSelected = character(), resamplingScore=character(), lambdaPathScore = character())

#set the lowest number of trans genes to be remaining in the module
numTrans = 3

# for each GO term cluster in the aberration

for (i in 1:max(clusterGO)){ 

idx = which(clusterGO == i)

#go terms associated with the cluster
goTerms  = names(clusterGO)[idx]

#get the cis genes in the cluster
idx = which(goTable$pathway %in% goTerms)
goAnnot = unique(goTable$annotation[idx])
cisGenesInModule = goTable$cisGene[idx]

#get the trans genes in the cluster. We can use topGO to find the
#genes that were annotated with these terms
xx2 = numeric(length=nrow(transExp))
names(xx2) = rownames(transExp)
GOdata <- new("topGOdata", ontology = "BP", allGenes = xx2, geneSel = function(p) p > 0.3, description = "Test", 
annot = annFUN.org, mapping="org.Hs.eg.db", ID="entrez")
transGenesInModule = genesInTerm(GOdata, as.character(goTerms))
transGenesInModule = unique(unlist(transGenesInModule))

idx = which(rownames(transExp) %in% transGenesInModule)
transExpModule = t(transExp[idx,])
idx = which(rownames(cisExp) %in% cisGenesInModule)
cisExpModule = t(cisExp[idx,])

## Perform Sparse CCA ##
perm.out <- CCA.permute(transExpModule, cisExpModule)
c = CCA(transExpModule, cisExpModule, penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)

# if canonical correlation is greater than 0.7, remove cis and trans genes with zero coefficients
# and perform multi-task ranking
if (c$cor[1] > .7){

idx = which(c$u != 0)
nonZeroTargetExpression = transExpModule[,idx]
idx = which(c$v != 0)
nonZeroRegExpression = cisExpModule[,idx]

if (is.null(ncol(nonZeroRegExpression))==FALSE & is.null(ncol(nonZeroTargetExpression)) == FALSE) {
if (ncol(nonZeroTargetExpression) >= numTrans) {
results.multi.resamp = multiTaskResampling(nonZeroRegExpression, nonZeroTargetExpression)
inclusionProbabilities = results.multi.resamp
rankings = multiTaskRank(nonZeroRegExpression, nonZeroTargetExpression)

## Add GO term module with rankings to the module table

goTerms2 = paste(goTerms, sep=",", collapse=",")
goAnnot2 = paste(goAnnot, sep=",", collapse=",")
selectedCisGenes = paste(as.character(colnames(nonZeroRegExpression)), sep=",", collapse=",")
lambdaPathScore = paste(rankings, sep=",", collapse=",")
resamplingRank = paste(inclusionProbabilities, sep=",", collapse=",")


#this is where the results would be added to file
moduleTable <- rbind(moduleTable, data.frame(goTerms=goTerms2, goAnnotation=goAnnot2, 
canonicalCorrelation = c$cor[1], cisGenesSelected=selectedCisGenes, resamplingScore = resamplingRank, 
lambdaPathScore=lambdaPathScore))

} #if number of remaining targets greater than or equal to numTrans
}
} # if canoncor > 0.7

} # for each go term module

#write out the module table
#write.table(moduleTable, file="moduleTable.txt", sep="\t")

