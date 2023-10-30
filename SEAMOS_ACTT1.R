# Gail Potter
# MAY 2021
# ACTT1 analyze subgroup effect with SEAMOS
# R version 4.0.4


library(haven)
library(dplyr)
library(survival)
library(survminer)
library(MASS)

rm(list=ls())
indata1="N:/COVID-19/ACTT/ACTT-1/Data/"
outdir='C:/Users/potterge/Box/Mike Proschan/'
actt1=read_sas(paste(indata1,"adtte.sas7bdat",sep=""))
mort1 = actt1 %>% filter(PARAMCD=='D29DTHEF' & ITTFL=='Y')


# Data for patients who died were censored at day 29.
mort1$CNSR[mort1$AVAL>28]=1
mort1$AVAL[mort1$AVAL>28]=28

mod <- survfit(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort1)

mort=mort1 %>% filter(!is.na(BCSOSN))

mort$os = factor(mort$BCSOSN)

delta.i  = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==4)))$coef[1,1],
             summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5)))$coef[1,1],
             summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6)))$coef[1,1],
             summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7)))$coef[1,1])

delta.ses = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==4)))$coef[1,3],
              summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5)))$coef[1,3],
              summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6)))$coef[1,3],
              summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7)))$coef[1,3])


## For SEAMOS, the overall effect, delta.bar, is the estimate from a model ignoring subgroups:
mod.observed.ignore.subgroup = coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort)

delta = mod.observed.ignore.subgroup$coef
delta

standardized.effects.observed = (delta.i - delta)/delta.ses
min(standardized.effects.observed)
max(standardized.effects.observed)

## Get permutation distribution
nsim=5000
low=NULL
high=NULL

set.seed(10)
for (s in 1:nsim){
  mort$permuted.os = mort$os[sample(1:nrow(mort))]
  
  delta.ses = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==4)))$coef[1,3],
                summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==5)))$coef[1,3],
                summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==6)))$coef[1,3],
                summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==7)))$coef[1,3])
  
  
  delta.i  = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==4)))$coef[1,1],
               summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==5)))$coef[1,1],
               summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==6)))$coef[1,1],
               summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(permuted.os==7)))$coef[1,1])
  
  std.effects = (delta.i - delta)/delta.ses  
  ## delta estimate is same as that for the observed data since it ignores subgroup
  
  low[s]=min(std.effects)
  high[s]=max(std.effects)
}

standardized.effects.observed
mean(low<=min(standardized.effects.observed))
mean(high >= max(standardized.effects.observed))


ci90 = c(quantile(low, .05), quantile(high, .95))
ci95 = c(quantile(low, .025), quantile(high, .975))
round(ci90,2)
round(ci95,2)
round(standardized.effects.observed,2)
names(standardized.effects.observed)=paste("OS", 4:7)

png(paste(outdir,'actt1_seamos.png',sep=''),res=300,
    width=1700,height=1400)
plot(standardized.effects.observed, c(2,1,3,4), ylim=c(0,6), xlim=c(-3.25,3.25), 
     pch=16,ylab='',xlab='Standardized difference from overall effect',axes=FALSE)
abline(v=0, col='gray', lty=2) ;
text(standardized.effects.observed, c(2,1,3,4),
     labels=names(standardized.effects.observed), pos=2)
axis(side=1, at=seq(-3,3,), pos=0)
segments(x0=ci95, x1=ci95, y0=0, y1=6, col='gold')
segments(x0=ci90, x1=ci90, y0=0, y1=5.25, col='blue')
text(0,5,'90% interval', col='blue')
text(0,5.75,'95% interval', col='gold')

arrows(x0=c(-1, 1), x1=ci90, y0=5, y1=5, col='blue', length=.05)
arrows(x0=c(-1, 1), x1=ci95, y0=5.75, y1=5.75, col='gold', length=.05)

dev.off()
