# Programmer: Gail Potter
# Date: AUG 2024
# SEAMOS and QMAX power on biowulf



library(dplyr)
library(ggplot2)
library(haven)
library(logistf)
library(stringr)
library(getopt)
library(parallel)

rm(list=ls())

event.rate.control1 = rep(.4,4); event.rate.treated1 = c( .05, .4,.4, .4)
event.rate.control2 = c(.4, .05,.2, .6); event.rate.treated2 = c(.05,.05, .2, .6)
ss1=rep(30,4); ss2=c(48,24,24,24)


ss=ss2
event.rate.control=event.rate.control2
event.rate.treated=event.rate.treated2



cmd_options= commandArgs(TRUE);
seed=getopt(matrix(c('seed','s', 1, "integer"), byrow=TRUE, ncol=4))$seed


source('powerfun_q_seamos_diff_prop_cluster.R')

n.trials = 200 ## number of trials per task.  
## Want 10000, so do 50 tasks of 200 trials
nsim=500  ## number of permutations



detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(2)
  } 
  return(ncores) 
}

mc=detectBatchCPUs()
cl=makeCluster(mc)
clusterSetRNGStream(cl, seed)
clusterExport(cl, varlist=list('powerfun', 'ss', 'event.rate.control', 
                               'event.rate.treated', 'n.trials', 'nsim'))

reject.vectors = clusterEvalQ(cl, {
  ## returns a list: (1) vector of H0 rejections for n.trials by SEAMOS,
  ##                 (2) vector of H0 rejections for n.trials by Qmax
  ## compile the list and calculate % of NA
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(haven)
  library(logistf)
  library(stringr)

  powerfun(
  outcome.rates.control=event.rate.control, 
  outcome.rates.treated=event.rate.treated, n.trials=n.trials, nsim=nsim, 
  n.per.group.control=ss,  n.per.group.treated=ss, alpha=0.025,
  plots=FALSE, outdir="G:/") 
  })
stopCluster(cl)

save.image("aug16.Rdata")

if (sum(event.rate.control==event.rate.treated) ==0)
  outcome='Type1_error' else outcome='Power'

filename =  paste(outcome,'___rates', 
            paste(event.rate.control, collapse='_'), 
            'ss_',paste(ss, collapse='_'), sep='','.Rda')

save(reject.vectors, file=filename)