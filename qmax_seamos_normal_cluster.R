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

means1=c(0,0,0,0)
means2=c(1000,0,0,0)

ss1=rep(30,4); ss2=c(48,24,24,24)


ss=ss1
means=means2
delta1=0



cmd_options= commandArgs(TRUE);
seed=getopt(matrix(c('seed','s', 1, "integer"), byrow=TRUE, ncol=4))$seed


source('powerfun_normal_cluster.R')

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
clusterExport(cl, varlist=list('powerfun', 'ss', 'means', 'delta1', 'n.trials', 'nsim'))

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

  powerfun(n.per.arm.per.subgroup=ss, delta1=delta1, control.means=means, 
           sigma=1, alpha=.025, n.trials=n.trials, nsim=nsim, 
           plots=FALSE, outdir="G:/")
  })
stopCluster(cl)

save.image("aug16.Rdata")

if (delta1==0)  outcome='Type1_error' else outcome='Power'

filename =  paste(outcome,'___means', 
            paste(means, collapse='_'), 
            'ss_',paste(ss, collapse='_'), sep='','.Rda')

save(reject.vectors, file=filename)