###******* multi-task resampling ***** ###
#Calculate "inclusion probabilities" ie the number of times that the 
#cis gene was selected in the multi-task LASSO model when resampling
#the population at 90% 50 timees
#regExpression - non-zero (in sparse CCA) cis genes in GO term module
#targExpression - non-zero (in sparse CCA) trans in the GO term module

multiTaskResampling <- function(regExpression, targExpression) {

inclusionProbabilities = numeric(length = ncol(regExpression))
nSamps = length(targExpression[,1])

for (k in 1:50){
samples = sample(1:nSamps, ceiling(.9*nSamps))
predictors = regExpression[samples,]
responses = targExpression[samples,]
#group lasso to select the cis drivers
r.cv = cv.glmnet(predictors, responses, family = 'mgaussian', nfolds = 10, alpha = 1)
lambdaMin = coef(r.cv, s = "lambda.1se")
idx = which(lambdaMin[[1]][1:ncol(regExpression)+1,1] != 0)
if (length(idx) > 0) {
inclusionProbabilities[idx] = inclusionProbabilities[idx]+1
} else if (length(idx) == 0) {
lambdaMin = coef(r.cv, s = "lambda.min")
idx = which(lambdaMin[[1]][1:ncol(regExpression)+1,1] != 0)
inclusionProbabilities[idx] = inclusionProbabilities[idx]+1
} #elseif cis gene selection 
} # 50 resamples 
return(inclusionProbabilities)
}

###****Multi-task ranking*****###
#Rank cis genes by the order that they appear in the lambda path
multiTaskRank = function(regExpression, targExpression) {

r.cv = cv.glmnet(regExpression, targExpression, family = 'mgaussian', nfolds = 10, alpha = 1)
rankings1 = numeric(length=ncol(regExpression))
cc = 1
s = vector(length=length(r.cv[[8]]$beta[[1]][1,]))
#looping over columns/lambda for a gene
for (x in 1:length(r.cv[[8]]$beta[[1]][1,])){

#if there is non-zero coefficients
s[x] = sum(r.cv[[8]]$beta[[1]][,x] != 0)
if (s[x] > 0){

#find which regulators are not zero
idx = which(r.cv[[8]]$beta[[1]][,x] != 0)
for (k in 1:length(idx)) {
ix = idx[k]

#if there is not a number assigned, assign it the current ranking
if ((rankings1[ix])==0){
rankings1[ix] = cc
} #
} #
if (s[x] > s[x-1]){
cc = cc + 1
}
} #if
} #lambdaFor
return(rankings1)
}
