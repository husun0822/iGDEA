library(DESeq2)
goodsample = F
while(goodsample==F){
  simdata = makeExampleDESeqDataSet(betaSD = 1)
  countdata = counts(simdata)
  # our function cannot handle the case of allzero gene
  allzero = apply(countdata,1,function(x) ifelse(sum(x==0)==length(x),0,1))
  if(sum(allzero==0)==0){
    goodsample = T
  }
}

designmatrix = ifelse(as.matrix(simdata$condition)=="A",0,1)

## try DESeq2
object = DESeq(simdata, betaPrior = F, useT = T)
deres = as.data.frame(mcols(object))

newcountdata = countdata


## try our package
source("main.R")


newdesign = cbind(rep(1,length(designmatrix)),designmatrix)
colnames(newdesign) = c("(Intercept)","Treatment")
mobject = DEA_deseq(newcountdata,newdesign)
res = data.frame(dispGeneEst = mobject$dispGeneEst, dispFit = mobject$dispFit, dispMAP = mobject$dispersions,
                 pval = mobject$Wald_Pvalue[,2])

# compare the p-value between our negative binomial implementation and the DESeq2 package

library(ggplot2)
plotdata = data.frame(DESeq2 = deres$WaldPvalue_condition_B_vs_A, NBGLM = res$pval)
p = ggplot(data = plotdata, aes(x = NBGLM, y = DESeq2)) + 
  geom_point()+
  theme_bw()

p

# compare the running time

generate_data = function(samplesize = 100){
  goodsample = F
  while(goodsample==F){
    simdata = makeExampleDESeqDataSet(n = 500,m = samplesize, betaSD = 1)
    countdata = DESeq2::counts(simdata)
    allzero = apply(countdata,1,function(x) ifelse(sum(x==0)==length(x),0,1))
    if(sum(allzero==0)==0){
      goodsample = T
    }
  }
  designmatrix = ifelse(as.matrix(simdata$condition)=="A",0,1)
  
  return(list(simdata = simdata, countdata = countdata, designmatrix = designmatrix))
}

deseq_time = NULL
nbglm_time = NULL

size = round(seq.int(from = 100, to = 1000, length.out =50))

for (s in size){
  idx = which(size==s)
  sprintf("Progress: %d/%d",idx,length(size))
  L = generate_data(s)
  simdata = L$simdata
  countdata = L$countdata
  designmatrix = L$designmatrix
  
  start_time = Sys.time()
  object = DESeq(simdata, betaPrior = F, useT = T, quiet = T)
  end_time = Sys.time()
  
  deseq_time = c(deseq_time,as.numeric(end_time - start_time))
  
  newdesign = cbind(rep(1,length(designmatrix)),designmatrix)
  colnames(newdesign) = c("(Intercept)","Treatment")
  
  start_time = Sys.time()
  mobject = DEA_deseq(countdata,newdesign)
  end_time = Sys.time()
  
  nbglm_time = c(nbglm_time,as.numeric(end_time - start_time))
}


# plot the running time
timedata = data.frame(gene_size = size, DESeq2 = deseq_time[1:50], NBGLM = nbglm_time[1:50])
ptime = ggplot(data = timedata)+
  geom_line(aes(x = gene_size, y = DESeq2,color="DESeq2"))+
  geom_line(aes(x = gene_size, y = NBGLM, color = "NBGLM"))+
  theme_bw()+
  ylab("running time (s)")

ptime

library(gridExtra)
resultplot = grid.arrange(p,ptime,nrow=2)
resultplot
