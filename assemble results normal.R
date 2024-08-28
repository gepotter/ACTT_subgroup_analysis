# Gail Potter
# Aug 2024
# SEAMOS + Qmax power type 1 error
# assemble results from cluster

library(xtable)

setwd("G:/")
rm(list=ls())


ss1=rep(30,4); ss2=c(48,24,24,24)

get.results = function(means, ss, result='error')  {
  if (result=='error'){
  load(paste('Type1_error___means', 
             paste(means, collapse='_'), 
             'ss_',paste(ss, collapse='_'), sep='','.Rda')) 
    } else
 if (result=='power')
   load(paste('Power___means', 
              paste(means, collapse='_'), 
              'ss_',paste(ss, collapse='_'), sep='','.Rda')) else
                error("Result must be error or power")

  seamos.reject=NULL
  Qmax.reject=NULL
  for (s in 1:length(reject.vectors)) {
    Qmax.reject = c(Qmax.reject, reject.vectors[[s]]$Qmax)
    seamos.reject = c(seamos.reject, reject.vectors[[s]]$SEAMOS)
  }
  length(seamos.reject)
  length(Qmax.reject)
  
  
  return(c('SEAMOS rejections'=mean(seamos.reject, na.rm=TRUE),
        'Qmax rejections'=mean(Qmax.reject, na.rm=TRUE)))
  
}


get.error.and.power = function(means, ss){
  power.results =
    get.results(means=means, ss=ss, result='power')
  error.results =
    get.results(means=means, ss=ss, result='error')
  list(results= c("SEAMOS pow" = power.results[1], "Qmax pow"=power.results[2],
    "SEAMOS err" = error.results[1], "Qmax err"=error.results[2]))
}

means1=c(0,0,0,0)
means2=c(1000,0,0,0)
r1=get.error.and.power(means=means1,ss=ss1)$results
r1
r2=get.error.and.power(means=means1,ss=ss2)$results
r2
r3=get.error.and.power(means=means2,ss=ss1)$results
r3
r4=get.error.and.power(means=means2,ss=ss2)$results
r4
sss = rbind(toString(ss1),toString(ss2),toString(ss1),toString(ss2))
means = rbind(cbind(toString(means1),toString(means1+c(1,0,0,0))),
              cbind(toString(means1),toString(means1+c(1,0,0,0))),
              cbind(toString(means2),toString(means2+c(1,0,0,0))),
              cbind(toString(means2),toString(means2+c(1,0,0,0))))
tab=cbind(sss, means, round(rbind(r1,r2,r3,r4),3))              
colnames(tab)=c("Sample sizes", "Control Means","Treated means" ,
                "Seamos power", "Qmax power", "Seamos err", "Qmax err")
tab

xtable(tab)
