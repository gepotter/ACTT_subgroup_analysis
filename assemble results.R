# Gail Potter
# Aug 2024
# SEAMOS + Qmax power type 1 error
# assemble results from cluster

library(xtable)

setwd("G:/")
rm(list=ls())


event.rate.control1 = rep(.4,4)
event.rate.treated1 = c( .05, .4,.4, .4)

event.rate.control2 = c(.4, .05,.2, .6)
event.rate.treated2 = c(.05,.05, .2, .6)

ss1=rep(30,4)
ss2=c(48,24,24,24)

get.power = function(treatment.rates, control.rates, ss){
  load(paste('Power___rates', 
             paste(control.rates, collapse='_'), 
             'ss_',paste(ss, collapse='_'), sep='','.Rda'))
  
  seamos.reject=NULL
  Qmax.reject=NULL
  for (s in 1:length(H0.reject.vectors)){
    # maybe should be reject.vectors
    Qmax.reject = c(Qmax.reject, H0.reject.vectors[[s]]$Qmax.reject)
    seamos.reject = c(seamos.reject, H0.reject.vectors[[s]]$seamos.reject)
  }
  length(seamos.reject)
  length(Qmax.reject)
  
  
  return(c('SEAMOS rejections'=mean(seamos.reject, na.rm=TRUE),
           'Qmax rejections'=mean(Qmax.reject, na.rm=TRUE),
           'SEAMOS NA rate '=mean(is.na(seamos.reject)),
           'Qmax NA rate '=mean(is.na(Qmax.reject)))_
         
         
}

get.results = function(control.rates, treated.rates, ss, result='error')  {
  if (result=='error'){
  load(paste('Type1_error___rates', 
             paste(control.rates, collapse='_'), 
             'ss_',paste(ss, collapse='_'), sep='','.Rda')) 
    reject.vectors=H0.reject.vectors 
    } else
 if (result=='power')
   load(paste('Power___rates', 
              paste(control.rates, collapse='_'), 
              'ss_',paste(ss, collapse='_'), sep='','.Rda')) else
                error("Result must be error or power")

  seamos.reject=NULL
  Qmax.reject=NULL
  for (s in 1:length(H0.reject.vectors)) {
    Qmax.reject = c(Qmax.reject, reject.vectors[[s]]$Qmax.reject)
    seamos.reject = c(seamos.reject, reject.vectors[[s]]$seamos.reject)
  }
  length(seamos.reject)
  length(Qmax.reject)
  
  
  return(c('SEAMOS rejections'=mean(seamos.reject, na.rm=TRUE),
        'Qmax rejections'=mean(Qmax.reject, na.rm=TRUE),
        'SEAMOS NA rate '=mean(is.na(seamos.reject)),
        'Qmax NA rate '=mean(is.na(Qmax.reject))))
  
}


get.error.and.power = function(control.rates, treated.rates, ss){
  power.results =
    get.results(control.rates=control.rates, treated.rates = treated.rates, 
                ss=ss, result='power')
  error.results =
    get.results(control.rates=control.rates, treated.rates = control.rates, 
                ss=ss, result='error')
  list(results= c("SEAMOS pow" = power.results[1], "Qmax pow"=power.results[2],
    "SEAMOS err" = error.results[1], "Qmax err"=error.results[2]),
    pct.na = c(power.results[3:4], error.results[3:4]))
}

r1=get.error.and.power(control.rates=event.rate.control1, 
                       treated.rates=event.rate.treated1,ss=ss1)$results
r1
r2=get.error.and.power(control.rates=event.rate.control1, 
                       treated.rates=event.rate.treated1,ss=ss2)$results
r2
r3=get.error.and.power(control.rates=event.rate.control2, 
                       treated.rates=event.rate.treated2,ss=ss1)$results
r3
r4=get.error.and.power(control.rates=event.rate.control2, 
                       treated.rates=event.rate.treated2,ss=ss2)$results
r4
sss = rbind(toString(ss1),toString(ss2),toString(ss1),toString(ss2))
rates = rbind(toString(event.rate.control1),toString(event.rate.treated1),
              toString(event.rate.control2),toString(event.rate.treated2))
              
tab = cbind(sss, rates, round(rbind(r1,r2,r3,r4),3))
colnames(tab)=c("Sample sizes", "event rates", "Seamos power", "Qmax power", "Seamos err", "Qmax err")
tab

xtable(tab)
na.rates = rbind(get.error.and.power(control.rates=event.rate.control1, 
                 treated.rates=event.rate.treated1,ss=ss1)$pct.na,
                 get.error.and.power(control.rates=event.rate.control1, 
                       treated.rates=event.rate.treated1,ss=ss2)$pct.na,
                 get.error.and.power(control.rates=event.rate.control2, 
                       treated.rates=event.rate.treated2,ss=ss1)$pct.na,
                 get.error.and.power(control.rates=event.rate.control2, 
                       treated.rates=event.rate.treated2,ss=ss2)$pct.na)
na.rates


# tab=cbind (sample.sizes, event.rates, 
#            round(rbind(c(r11[1:2],r12[1:2]), 
#                        c(r21[1:2],r22[1:2]), 
#                        c(r31[1:2], r32[1:2]), 
#                        c(r41[1:2],r42[1:2])),3))
# colnames(tab)=c("Subgroup sizes" ,"Event rates", "SEAMOS power", "Qmax power", "SEAMOS error", "Qmax error")
# tab
# 
# 
# data.frame(tab) %>%
#   rtf_body() %>%
#   rtf_encode() %>%
#   write_rtf(paste(outdir,'non_normal_tab_propdiff.rtf',sep='/'))
# 
# save.image(paste(outdir,'/seamos_nonnormal_propdiff.Rdata', sep=''))
