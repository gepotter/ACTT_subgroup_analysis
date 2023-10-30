# Gail Potter
# MAY 2021
# ACTT1, analyze subgroup effect with Q and Qmax
# INPUT: adtte.sas7bdat
# R version 4.0.4

library(haven)
library(dplyr)
library(survival)
library(survminer)
library(MASS)
library(r2rtf)

rm(list=ls())
indata1="N:/COVID-19/ACTT/ACTT-1/Data/"
outdir = "C:/Users/potterge/Box/Mike Proschan"

actt1=read_sas(paste(indata1,"adtte.sas7bdat",sep=""))
mort1 = actt1 %>% filter(PARAMCD=='D29DTHEF' & ITTFL=='Y')

# Data for patients who died were censored at day 29.
mort1$CNSR[mort1$AVAL>28]=1
mort1$AVAL[mort1$AVAL>28]=28

mod <- survfit(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort1)


# REPLICATE PUBLISHED ANALYSES
# hazard ratios were calculated from the stratified Cox model
# The P value for this ratio was calculated with the stratified 
# log-rank test (overall model stratified by actual disease severity).

# stratified log rank test
p=pchisq(survdiff(Surv(AVAL, (1-CNSR)) ~ 
      TRTP+strata(STRATAV), data = mort1)$chisq,
      df=1,lower.tail=FALSE)
p

mort=mort1 %>% filter(!is.na(BCSOSN))
summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort))
summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==4)))
summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5)))
summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6)))
summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7)))

## OS 5, 6, 7 only

delta.i = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5)))$coef[1,1],
              summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6)))$coef[1,1],
              summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7)))$coef[1,1])

delta.ses = c(summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5)))$coef[1,3],
        summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6)))$coef[1,3],
        summary(coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7)))$coef[1,3])


w.i  = 1/(delta.ses^2)
w = sum(w.i)

delta = (1/w)*sum(w.i*delta.i)
delta
delta.i - delta


(q.observed =  sum(w.i*(delta.i-delta)^2))
(qmin.observed = min((delta.i-delta)/sqrt((w-w.i)/(w.i*w))))


## Compare Qmin.observed to a normal distribution with Bonferroni correction
## Critical value
critval = qnorm(0.025/3, lower.tail=TRUE)
pnorm(critval, lower.tail=TRUE)*3
pnorm(qmin.observed, lower.tail=TRUE)*3

pchisq(q.observed, lower.tail=FALSE, df=2)

# test statistic for OS 5
os5stat = (delta.i-delta)/sqrt((w-w.i)/(w.i*w))
pnorm(os5stat)
