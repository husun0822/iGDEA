To run a negative binomial model for gene differential expression analysis, data inputs are:

1. gene expression count matrix: n by m matrix. n is the number of genes and m is the number of samples.
note: i. the gene expression count matrix should contain non-negative integers only.
     ii. our current implementation does not allow for a gene to have zero counts in all of the sample, so please check 
         if a certain gene has zero counts in the whole row.
         
2. design matrix: m by k matrix. m is the number of samples and k is the number of covariates.
note: i. the design matrix's first column must be all 1s, which serves as an intercept of the GLM model.
     ii. other covariates should describe the treatment condition of each sample, our current version does not allow for "expanded"
         design matrix, which means that a binary variable can only be represented in 0-1 form in a single covariate, but not
         in two separate dummy variables. The "expanded" matrix is supported by DESeq2, but not here.
    iii. the column name of the design matrix must be: c("(Intercept)","Treatment1",...)
    
    
A. Example vignette:
# suppose we have a "countdata" object for gene count data and a "modelMatrix" object for design matrix information
source("main.R")
mobject = DES_deseq(countdata,modelMatrix)

    
B. Example call of the NB-GLM model from simulated data:

library(DESeq2)
source("main.R")

goodsample = F
while(goodsample==F){
  simdata = makeExampleDESeqDataSet(betaSD = 1)
  countdata = DESeq2::counts(simdata)
  # our function cannot handle the case of allzero gene
  allzero = apply(countdata,1,function(x) ifelse(sum(x==0)==length(x),0,1))
  if(sum(allzero==0)==0){
    goodsample = T
  }
}

designmatrix = ifelse(as.matrix(simdata$condition)=="A",0,1)
newdesign = cbind(rep(1,length(designmatrix)),designmatrix)
colnames(newdesign) = c("(Intercept)","Treatment")

mobject = DEA_deseq(countdata,newdesign)
res = data.frame(dispGeneEst = mobject$dispGeneEst, dispFit = mobject$dispFit, dispMAP = mobject$dispersions,
                 pval = mobject$Wald_Pvalue[,2])
