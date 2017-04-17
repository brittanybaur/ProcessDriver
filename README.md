ProcessDriver is a tool to detect copy number cancer drivers and associated biological processes
Baur and Bozdag, "ProcessDriver: A computational pipeline to identify copy number drivers and associated disrupted biological processes in cancer", Accepted to Genomics.

GO Association function are in the goAssocationFunctions.R file. You will need a matrix of cis gene expression and a matrix of trans gene expression. This will associate the cis genes with biological processes through the trans genes and cluster similar go terms. It depends on topGO and GOSim, which you can obtain in bioconductor. 

The driverSelectionStep.R file has code that depends on MultiTaskFuncions.R. You will need the matrices from the GO assocation step, as well as the "go table" produced from the goAssocation() function and the clustering of goTerms produced from the goClustering() function. This step runs sparse CCA (the version from Witten et al.'s PMA package) and if the canonical correlation is greater than 0.7 will rank the remaining cis genes as drivers based on two different methods utilizing multi-task LASSO. 

Additional comments and descriptions are in each of the files. 
