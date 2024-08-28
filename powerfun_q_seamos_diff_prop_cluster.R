# function to calculate power for binary data with treatment
# effect measured as difference in proportions

powerfun = function(outcome.rates.control, outcome.rates.treated, n.trials, nsim, 
                    n.per.group.control, n.per.group.treated, alpha=0.025,
                    plots=FALSE, outdir) {
  
  # outcome.rates.control = subgroup-specific outcome rates in controls
  # outcome.rates.treated = subgroup-specific outcome rates in treated
  # n.trials = number of simulated trials for power approximation
  # nsim = number of resamples for SEAMOS
  # n.per.group.treated = vector of numbers in each subgroup in treatment group
  # n.per.group.controls = vector of numbers in each subgroup in control group

  if (sum(outcome.rates.treated!=outcome.rates.control)==0) calculation='Type 1 error' else calculation ='power'
  plotname = paste(calculation, "SS", toString(n.per.group.control),"rates",toString(outcome.rates.control), ' ')
  
  rates = as.vector(rbind(outcome.rates.control, outcome.rates.treated))
  ns =    as.vector(rbind(n.per.group.control, n.per.group.treated))
  
  indiv.probs = rep(rates, ns)
  
  os=rep(4:7, n.per.group.control + n.per.group.treated)
  
  treated=rep(rep(0:1,4), ns)
  
  n.total = length(os)
  
  reject=NULL
  reject.Qmax=NULL; 
  emptymat = matrix(NA, nrow=n.trials, ncol=4)
  pc=emptymat; pt=emptymat; effects=emptymat; 
  semat = emptymat;Qstatmat=emptymat;rejectmat=emptymat;
  seamos.test.stats = emptymat;
  num.discarded=NULL; Qmin=NULL
  
  for (i in 1:n.trials){
    x=rbinom(n=n.total, prob=indiv.probs, size=1)
    
    p.control = tapply(x[treated==0], os[treated==0], mean)
    p.treated = tapply(x[treated==1], os[treated==1], mean)
    
    pc[i,]=p.control
    pt[i,]=p.treated
    
    delta.i = p.treated - p.control
    effects[i,]=delta.i
    delta.ses = sqrt(
      p.treated*(1-p.treated)/as.vector(table(os[treated==1]))+
      p.control*(1-p.control)/as.vector(table(os[treated==0]))
    )
    v.i = delta.ses^2
    w.i = 1/v.i
    w = sum(w.i)
    delta = (1/w)*sum(w.i * delta.i)  
    ## overall effect is for Qmax.
    
    Qstats = (delta.i - delta)/sqrt((w-w.i)/(w.i*w))
    rejectmat[i,] = pnorm(Qstats) <= alpha/(length(outcome.rates.control))

    Qmin[i]=min(Qstats)
    ## 19JUL
    Qstatmat[i,]=Qstats
    ##
    
    ## Power is to detect a difference from overall effect in first subgroup only (Qstats[1])
    if (calculation=='power')     
      ### test subgroup 1.  ACTT tests subgroup 2
      reject.Qmax[i] = pnorm(Qstats[1]) <= alpha/4 else
      if (calculation=='Type 1 error')     reject.Qmax[i] = pnorm(min(Qstats)) <= alpha/4
    
    
    ## SEAMOS TEST STATISTIC
    overall.effect  = mean(x[treated==1])-mean(x[treated==0])
    
    observed.test.stat = (delta.i - overall.effect)/sqrt(v.i)
    
    seamos.test.stats[i,] = observed.test.stat
    
    ##    shuffle subgroup labels
    min.test.stat=NULL

    if (i==1){
      permuted.test.stats=matrix(NA, nrow=nsim, ncol=4)      
      shuffle.effect.mat=matrix(NA, ncol=4, nrow=nsim)
      shuffle.pc=matrix(NA, ncol=4, nrow=nsim)
      shuffle.pt=matrix(NA, ncol=4, nrow=nsim)
    }
    for (s in 1:nsim){
      sim.os = sample(os)
      
      ## need at some of each treatment group in each OS
      while(min(table(sim.os,treated))==0)      sim.os = sample(os)
      
      p.control = tapply(x[treated==0], sim.os[treated==0], mean)
      p.treated = tapply(x[treated==1], sim.os[treated==1], mean)
      
      if (i==1){
        shuffle.pc[s,]=p.control
        shuffle.pt[s,]=p.treated
      }

      sim.subgroup.effects = p.treated - p.control
      ses = sqrt(
        p.treated*(1-p.treated)/as.vector(table(sim.os[treated==1]))+
          p.control*(1-p.control)/as.vector(table(sim.os[treated==0]))
      )

      min.test.stat[s] = 
        min((sim.subgroup.effects - overall.effect)/ses)
      if (i==1) {
        shuffle.effect.mat[s,]=sim.subgroup.effects
        permuted.test.stats[s,] = (sim.subgroup.effects - overall.effect)/ses
      }
    }
    ref.bound=quantile(min.test.stat, alpha, na.rm=TRUE)
    num.discarded[i] = sum(is.na(min.test.stat))
    ## number of permutations discarded because effect was inestimable
    
    if (calculation=='power')   reject[i] = observed.test.stat[1] < ref.bound else
      if (calculation=='Type 1 error')     reject[i] = min(observed.test.stat) < ref.bound
  }
  
  apply(pc, 2, mean)
  apply(pt, 2, mean)

  apply(pc, 2, sd)
  apply(pt, 2, sd)
  apply(shuffle.pc, 2, sd)
  apply(shuffle.pt, 2, sd)
  apply(seamos.test.stats, 2, sd)
  apply(effects, 2, sd)
  
  power = list("seamos.reject" = reject, 
            "Qmax.reject" = reject.Qmax) 
            # "% NA SEAMOS" = mean(is.na(reject)),
            # "% NA Qmax" = mean(is.na(reject.Qmax)))
  power
  
  if (plots) {
    dat=data.frame(outcome=c(pc, pt, shuffle.pc, shuffle.pt),
                   treated=rep(c("control","treated", "control", "treated"), 
                               4*c(n.trials, n.trials, nsim, nsim)),
                   os = factor(c(rep(rep(paste("OS", 4:7), each=n.trials), 2), 
                               rep(rep(paste("OS", 4:7), each=nsim),2))),
                   source=rep(c("Actual", "reshuffled"), 8*c(n.trials, nsim))) %>%
      mutate(treated.source = paste(treated, source))
    
                     
    png(paste(outdir,'/', plotname,'_diffprop_1.png',sep=''), 
        res=300, width=2200, height=2200)
    print(dat %>% ggplot()+geom_histogram(aes(x=outcome))+
      facet_grid(treated.source~os) +
      geom_vline(data=dat %>% group_by(treated.source, os) %>%
                   summarise(mean=mean(outcome)),aes(xintercept=mean),
                 colour='red') + theme_minimal() + xlab("Prob(event)"))
      dev.off()
      
      dat2=data.frame(effects = c(effects, apply(effects, 1, min),
                                  shuffle.effect.mat, 
                                  apply(shuffle.effect.mat, 1, min)),
                      os = c(rep(c(paste("OS", 4:7), "Min"), each=n.trials), 
                             rep(c(paste("OS", 4:7), "Min"), each=nsim)),
                      source=rep(c("Actual", "Reshuffled"), c(n.trials*5, nsim*5)))
                      
      dat3=data.frame(seamos.test.stats = c(seamos.test.stats, apply(seamos.test.stats, 1, min),
                                            permuted.test.stats, apply(permuted.test.stats, 1, min)),
                      os = c(rep(c(paste("OS", 4:7), "Min"), each=n.trials), 
                             rep(c(paste("OS", 4:7), "Min"), each=nsim)),
                      source=rep(c("Actual", "Reshuffled"), c(n.trials*5, nsim*5)))
      percentiles = dat3 %>% group_by(source, os) %>%
        summarise(ref.025=quantile(seamos.test.stats, .025, na.rm=TRUE))
      
      percentiles$ref.025[percentiles$os!='Min']=NA
      
      
      p1 =dat2 %>% ggplot()+geom_histogram(aes(x=effects))+
        facet_grid(source~os) +
        geom_vline(data=dat2 %>% group_by(source, os) %>%
                     summarise(mean=mean(effects)),aes(xintercept=mean),
                   colour='red') + theme_minimal() + xlab("Difference in proportions")
      
      p2= dat3 %>% ggplot()+geom_histogram(aes(x=seamos.test.stats))+
        facet_grid(source~os) +
        geom_vline(data=percentiles,aes(xintercept=ref.025),
                   colour='red') + theme_minimal() + xlab("SEAMOS test statistic")
      

      png(paste(outdir,'/', plotname,'_diffprop_2.png',sep=''), res=300, width=2200, height=2200)
      print(p1/p2)
      dev.off()
      
      
    } 
  return(power)
  
}
