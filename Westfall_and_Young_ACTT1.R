# Gail Potter
# FEB 2022
# ACTT1: calculate WFY p values for OS subgroup effects.
# R version 4.2.1

library(haven)
library(dplyr)
library(survival)
library(survminer)
library(MASS)
library(ggplot2)
library(ggsci)
library(patchwork)
library(NRejections)

rm(list=ls())
indata1="N:/COVID-19/ACTT/ACTT-1/Data/"
outdir="C:/Users/potterge/Box/Mike Proschan"
actt1=read_sas(paste(indata1,"adtte.sas7bdat",sep=""))
mort1 = actt1 %>% filter(PARAMCD=='D29DTHEF' & ITTFL=='Y')



# Data for patients who died were censored at day 29.
mort1$CNSR[mort1$AVAL>28]=1
mort1$AVAL[mort1$AVAL>28]=28
mort=mort1 %>% filter(!is.na(BCSOSN))


mod4=coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==4))
mod5=coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==5))
mod6=coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==6))
mod7=coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort %>% filter(BCSOSN==7))

## Two sided p-values for subgroup-specific effects on the ACTT data
two.sided.pvals = c(summary(mod4)$coef[5], summary(mod5)$coef[5], 
                    summary(mod6)$coef[5], summary(mod7)$coef[5])
round(two.sided.pvals, 3)

## One-sided p-values for subgroup-specific effects on the ACTT data
## We are only interested in negative coefficients (reduction in mortality)
pvals = c(pnorm(summary(mod4)$coef[4], lower.tail=TRUE),
          pnorm(summary(mod5)$coef[4], lower.tail=TRUE),
          pnorm(summary(mod6)$coef[4], lower.tail=TRUE),
          pnorm(summary(mod7)$coef[4], lower.tail=TRUE))
pvals
round(pvals, 3)

ordered.indices = order(pvals)


## Get permutation distribution.
nsim=10000
pvalsp=matrix(nrow=nsim, ncol=4)
colnames(pvalsp)=4:7
set.seed(99)
for (s in 1:nsim){
  mort$permuted.TRTP=NA
  mort$permuted.TRTP[mort$BCSOSN==4] = sample(mort$TRTP[mort$BCSOSN==4])
  mort$permuted.TRTP[mort$BCSOSN==5] = sample(mort$TRTP[mort$BCSOSN==5])
  mort$permuted.TRTP[mort$BCSOSN==6] = sample(mort$TRTP[mort$BCSOSN==6])
  mort$permuted.TRTP[mort$BCSOSN==7] = sample(mort$TRTP[mort$BCSOSN==7])
  mod4p=coxph(Surv(AVAL, (1-CNSR)) ~ permuted.TRTP, data = mort %>% filter(BCSOSN==4))
  mod5p=coxph(Surv(AVAL, (1-CNSR)) ~ permuted.TRTP, data = mort %>% filter(BCSOSN==5))
  mod6p=coxph(Surv(AVAL, (1-CNSR)) ~ permuted.TRTP, data = mort %>% filter(BCSOSN==6))
  mod7p=coxph(Surv(AVAL, (1-CNSR)) ~ permuted.TRTP, data = mort %>% filter(BCSOSN==7))
  
  # p values from permutation distribution
  pvalsp[s,] = c(pnorm(summary(mod4p)$coef[4], lower.tail=TRUE),
                 pnorm(summary(mod5p)$coef[4], lower.tail=TRUE),
                 pnorm(summary(mod6p)$coef[4], lower.tail=TRUE),
                 pnorm(summary(mod7p)$coef[4], lower.tail=TRUE))
}

# Retrieve the name of the subgroup with the minimum value
get.min.groupname = function(x, groupnames) groupnames[x==min(x)]
  
minp=apply(pvalsp, 1, min)  ## min p-value from each simulation
groups1= apply(pvalsp, 1, get.min.groupname,groupnames=colnames(pvalsp))
# Name of group with min p value from each simulation

## Remove group with smallest p-value from permutation set
pvals.minus1 = pvalsp[, -which(rank(pvals)==1)]
minp.minus1=apply(pvals.minus1, 1, min)
groups.minus1= apply(pvals.minus1, 1, 
                     get.min.groupname,groupnames=colnames(pvals.minus1))

pvals.minus2 = pvalsp[, -which(rank(pvals) %in% (1:2))]
minp.minus2=apply(pvals.minus2, 1, min)
groups.minus2= apply(pvals.minus2, 1, 
                     get.min.groupname,groupnames=colnames(pvals.minus2))

minp.minus3 = pvalsp[, -which(rank(pvals) %in% (1:3))]
## Only 1 group left.
groups.minus3=   rep(colnames(pvalsp)[-which(rank(pvals) %in% (1:3))], nsim) 


unadj.p = sort(round(pvals, 3))

areas = round(c(mean(minp <= pvals[ordered.indices[1]]),
          mean(minp.minus1 <= pvals[ordered.indices[2]]),
          mean(minp.minus2 <= pvals[ordered.indices[3]]),
          mean(minp.minus3 <= pvals[ordered.indices[4]])),3)

adj.p = sort(round(adj_Wstep(p=pvals, p.bt = t(pvalsp)),3))

dat =   data.frame(pval=c(minp, minp.minus1, minp.minus2, minp.minus3), 
                   OS=factor(c(groups1, groups.minus1, groups.minus2, groups.minus3),
                             levels=4:7, labels=4:7),
                   perm=factor(rep(1:4, each=nsim),
                               labels=c("Minimum", "Minimum excluding OS 5",
                                        "Minimum excluding OS 4,5",
                                        "Minimum excluding OS 4,5,6")))

data.segm = 
  data.frame(x=sort(pvals), xend=sort(pvals), y=0, yend=1800, 
             perm=unique(dat$perm))



textdat = data.frame(unadj.p = unadj.p, areas=areas,
                     Area=paste("Area =", areas),
                     op = c('=','=','=','<'),
                     p.adj.text = paste('p[(', c(1:3,3), ')]^{adj}'), 
                     adj.p=adj.p, 
                     unadj.p.labels=paste('p[(', 1:4, ')]==',unadj.p, sep=''),
                     perm = factor(1:4, labels=c("Minimum", "Minimum excluding OS 5",
                                            "Minimum excluding OS 4,5",
                                            "Minimum excluding OS 4,5,6")),
                     y=rep(1900, 4), x=c(.1, .12, .18, .26),
                     y.area = c(1700, 1500, 900, 900),
                     # exception = paste('p[(4)]^{adj}==p[(3)]^{adj}'),
                     # exceptionx = c(-2,-2,-2,.26),
                     exception = c('','','',paste('p[(4)]^{adj}==p[(3)]^{adj}')),
                     exceptionx =.26,
                     exceptiony=650)

textdat$adj.p2 = paste('=',textdat$adj.p)
textdat$adj.p2[1:3]=''

setwd(outdir)

png(file='wfpvals.png', res=300, width=2400,height=1900)
ggplot()+ 
  geom_histogram(aes(x=pval, fill=OS), data=dat, bins=30)+
  theme_classic() + scale_fill_aaas()+xlab("P-value") + ylab("")+#xlim(0,1.04)+
  facet_wrap(~perm)+
  theme(legend.position="bottom")+ 
  geom_segment(data=data.segm, aes(x=x,xend=xend,y=y,yend=yend),
               color='orange')+
  geom_text(data=textdat, color='orange',
            aes(x=unadj.p+.075, y=y, label=unadj.p.labels), parse=TRUE)+
  geom_text(data=textdat, aes(x=x, y=y.area, label=Area))+
  geom_text(data=textdat, aes(x=x+.175, y=y.area, label=op))+
geom_text(data=textdat, aes(x=x+.245, y=y.area, label=p.adj.text),parse=TRUE)+
  geom_text(data=textdat, aes(x=exceptionx, y=exceptiony, label=exception), parse=TRUE)+
  #xlim(-.05,1)+
  geom_text(data=textdat, aes(x=exceptionx+.2, y=exceptiony, label=adj.p2))
dev.off()

## Histograms of minimum p-values from the permutation distribution
## shaded by baseline ordinal scale


## For group whose p-value is rank j, take the max p-value 
## over i=1,..,j:

adjusted.pvalues1=NULL
for (i in 1:length(adj.pval.0))  adjusted.pvalues1[i]=max(adj.pval.0[1:i])

adjusted.pvalues=adjusted.pvalues1[ordered.indices]
adjusted.pvalues
adjusted.pvalues1

# table of unadjusted, Bonferroni-adjusted, and WFY adjusted p values
Bonferroni.adjusted = pvals*3
Bonferroni.adjusted[Bonferroni.adjusted>1]=1
tab = cbind(pvals, Bonferroni.adjusted, adjusted.pvalues)
colnames(tab)=c("Unadjusted", "Bonferroni", "WFY")
rownames(tab)=4:7


overall.mod = coxph(Surv(AVAL, (1-CNSR)) ~ TRTP, data = mort)
overall.p = pnorm(summary(overall.mod)$coef[4], lower.tail=TRUE)
round(overall.p,3)
round(tab,3)

getwd()
save.image("wfy.rdata")
