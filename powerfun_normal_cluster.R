## Gail Potter
## June  2024
## function to calculate power for normal outcom data

powerfun = function(n.per.arm.per.subgroup, delta1, control.means, 
                    sigma=1, alpha=.025, n.trials=200, nsim=200, 
                    plots=FALSE, outdir){

  ## n.per.arm.per.subgroup is a vector- one element for each subgroup
  
  k=length(n.per.arm.per.subgroup)
    print (paste("number of subgroups=",k))
    
    if (delta1==0) calculation='Type 1 error' else calculation ='power'
    plotname = paste(calculation, "SS", toString(n.per.arm.per.subgroup),"means",toString(control.means), ' ')
    
    n.total=2*sum(n.per.arm.per.subgroup)
    n.per.arm = n.total/2
    
    subgroups = rep (rep(1:k, n.per.arm.per.subgroup), 2) 
    subgroup.sizes = as.vector(table(subgroups))
    
    individual.control.means = rep(control.means, n.per.arm.per.subgroup)
    individual.treated.means = individual.control.means + 
      rep(c(delta1,rep(0, k-1)), n.per.arm.per.subgroup)
    individual.means = c(individual.control.means, individual.treated.means)
    treat = rep(c("control", "treated"), each=n.per.arm)
    print(table(individual.means, paste(subgroups,treat)))
    
    d = delta1/sqrt(2/n.per.arm.per.subgroup[1])
    
    reject=NULL
    reject.Qmax=NULL
    effects=matrix(NA, ncol=k, nrow=n.trials)
    Qstats=matrix(NA, ncol=k, nrow=n.trials)
    rejectmat=matrix(NA, ncol=k, nrow=n.trials)
    SEAMOS.stats=matrix(NA, ncol=k, nrow=n.trials)

    for (i in 1:n.trials){
      x=rnorm(n.total, mean=c(individual.control.means,individual.treated.means) , sd=sigma)
      sample.treated.means = tapply(x[treat=='treated'], subgroups[treat=='treated'], mean)
      sample.control.means = tapply(x[treat=='control'], subgroups[treat=='control'], mean)
      
      subgroup.effects = sample.treated.means - sample.control.means
      effects[i,]=subgroup.effects
      overall.effect = mean(x[treat=='treated']) - mean(x[treat=='control'])
      delta.i = subgroup.effects
      v.i = 2*sigma^2 / n.per.arm.per.subgroup
      w.i = 1/v.i
      w = sum(w.i)
      delta = (1/w)*sum(w.i * delta.i)  
    
      Qstats[i,] = (delta.i - delta)/sqrt((w-w.i)/(w.i*w))

      ## Power is to detect a difference from overall effect in first subgroup only (Qstats[1])
      if (calculation=='power')     
        reject.Qmax[i] = pnorm(Qstats[i,1], lower.tail=FALSE) <= alpha/k else
        if (calculation=='Type 1 error')     
          reject.Qmax[i] = pnorm(max(Qstats[i,]), lower.tail=FALSE) <= alpha/k else
          Error("Calculation must be power or Type 1 error")

      
      rejectmat[i,] = pnorm(Qstats[i,], lower.tail=FALSE) <= alpha/k
      
      ## SEAMOS test statistic
      
      observed.test.stat = (subgroup.effects - overall.effect)/sqrt(v.i)
      SEAMOS.stats[i,]=observed.test.stat
      
      ## Shuffle subgroup labels
      max.test.stat=NULL
      
      if (i==1) {
        reshuffled.seamos.stats = matrix(NA, nrow=nsim, ncol=k)
        shuffled.effects=matrix(NA, nrow=nsim, ncol=k)
      }
      
      for (s in 1:nsim){
        sim.subgroups = sample(subgroups)
        sim.sample.treated.means = tapply(x[treat=='treated'], sim.subgroups[treat=='treated'], mean)
        sim.sample.control.means = tapply(x[treat=='control'], sim.subgroups[treat=='control'], mean)
        sim.subgroup.effects = sim.sample.treated.means - sim.sample.control.means
        ns =  table(treat,sim.subgroups)
        treated.vars = tapply(x[treat=='treated'], sim.subgroups[treat=='treated'], var)
        control.vars = tapply(x[treat=='control'], sim.subgroups[treat=='control'], var)
        sim.v.i = control.vars / ns[1,] + treated.vars / ns[2,]
        
        if (i==1) {
          shuffled.effects[s,] = sim.subgroup.effects
          reshuffled.seamos.stats[s,]=(sim.subgroup.effects - overall.effect)/sqrt(sim.v.i)
        }
        max.test.stat[s] = max((sim.subgroup.effects - overall.effect)/sqrt(sim.v.i))
      }
      ref.bound=quantile(max.test.stat, .975, na.rm=TRUE)

      if (calculation=='power')   reject[i] = observed.test.stat[1] > ref.bound else
        if (calculation=='Type 1 error')     reject[i] = max(observed.test.stat) > ref.bound
      
    }

    power = c("SEAMOS" = mean(reject), "Qmax" = mean(reject.Qmax),
              "Qmax.analytical" = pnorm(d*sqrt((w-w.i[1])/w)-qnorm(1-alpha/k)))
    print(power)
    
    if (plots){
      pa=data.frame(x=x,treat=factor(treat),subgroups=factor(subgroups),
                    sim.subgroups=factor(sim.subgroups)) %>%
        ggplot()+geom_histogram(aes(x=x)) + facet_grid(treat~subgroups)+theme_classic()+
        ggtitle("Distribution of responses in 1 simulated trial")
      
      pb=data.frame(x=x,treat=factor(treat),subgroups=factor(subgroups),
                    sim.subgroups=factor(sim.subgroups)) %>%
        ggplot()+geom_histogram(aes(x=x)) + facet_grid(treat~sim.subgroups)+theme_classic()+
        ggtitle("Distribution of responses with reshuffle subgroups")
      
      png(paste(outdir,'/', plotname,' 1.png',sep=''), 
          res=300, width=2200, height=2200)
      print(  pa/pb)
      dev.off()
      
      dat2=data.frame(effects = c(effects, apply(effects, 1, max),
                                  shuffled.effects, 
                                  apply(shuffled.effects, 1, max)),
                      subgroup= c(rep(c(paste("Subgroup", 1:k), "Max"), each=n.trials), 
                                  rep(c(paste("Subgroup", 1:k), "Max"), each=nsim)),
                      source=rep(c("Actual", "Reshuffled"), c(n.trials*(k+1), nsim*(k+1))))
      
      dat3=data.frame(seamos.test.stats = c(SEAMOS.stats, apply(SEAMOS.stats, 1, min),
                                            reshuffled.seamos.stats, apply(reshuffled.seamos.stats, 1, min)),
                      subgroup = c(rep(c(paste("Subgroup", 1:k), "Max"), each=n.trials), 
                                   rep(c(paste("Subgroup", 1:k), "Max"), each=nsim)),
                      source=rep(c("Actual", "Reshuffled"), c(n.trials*(k+1), nsim*(k+1))))
      
      percentiles = dat3 %>% group_by(source, subgroup) %>%
        summarise(ref=quantile(seamos.test.stats, .975, na.rm=TRUE))
      
      percentiles$ref[percentiles$subgroup!='Max']=NA
      
    
      p1 =dat2 %>% ggplot()+geom_histogram(aes(x=effects))+
        facet_grid(source~subgroup) +
        geom_vline(data=dat2 %>% group_by(source, subgroup) %>%
                     summarise(mean=mean(effects)),aes(xintercept=mean),
                   colour='red') + theme_minimal() + xlab("Difference in means")
      
      p2= dat3 %>% ggplot()+geom_histogram(aes(x=seamos.test.stats))+
        facet_grid(source~subgroup) +
        geom_vline(data=percentiles,aes(xintercept=ref),
                   colour='red') + theme_minimal() + xlab("SEAMOS test statistic")
      png(paste(outdir,'/', plotname,' 2.png',sep=''), 
          res=300, width=2200, height=2200)
      print(p1/p2)
      dev.off()
    }
    
    reject.vectors = list("SEAMOS" = reject, "Qmax" = reject.Qmax,
              "Qmax.analytical" = pnorm(d*sqrt((w-w.i[1])/w)-qnorm(1-alpha/k)))
    return(reject.vectors)
}