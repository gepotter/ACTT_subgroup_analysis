FisherLSD<-function(x){
# function needed to compute
# power for QLSD test, the analog of Fisher's LSD
# for interactions
chisqcrit<-qchisq(1-2*alpha,k-1) 
f1<-dnorm(x,ddd*sqrt((w-w1)/w),1)
f2<-1-pchisq(chisqcrit-x^2,k-2)
return(f1*f2)
}

setwd('C:/Users/potterge/Box/Mike Proschan')

alpha<-0.025 # one-tailed alpha
dd<-seq(0,6,0.2) # standardized difference
dat=data.frame(d=NA, pow=NA, method=NA, k=NA)

# The scenario we consider is when the treatment 
# effects in categories 2,...,k are 0 and the 
# treatment effect in category 1 is 
# dd*sd(hat delta 1)


for (k in c(2:10)){
  
  # number of categories for the factor variable
  w1<-1 # weight for the stratum-specific estimator
  w<-k # weight for the overall estimator
  
  # Asymptotic power for the Qmax statistic under 
  # above scenario is:
  powQmax<-pnorm(dd*sqrt((w-w1)/w)-qnorm(1-alpha/k))
  
  # Asymptotic power for the Q statistic under above 
  # scenario is (no longer relevant):
  noncenpar<-dd^2*(w-w1)/w # noncentrality parameter 
  # for the chi-squared variable
  powQ<-1-pchisq(qchisq(1-2*alpha,k-1),k-1,
                 ncp=noncenpar)
  powQLSD=NULL
  for(ddd in dd){
    powQLSD<-c(powQLSD,integrate(FisherLSD,lower=
                                   qnorm(1-alpha),upper=Inf)$value)
  }
  
  newdat=data.frame(d=rep(dd,3), 
                    pow=c(powQmax, powQLSD, powQ), 
                    method=rep(c("Qmax","QLSD", "Q"), 
                               each=length(dd)), k=k)
  dat=rbind(dat, newdat)
  
}


dat2=dat %>% filter(!is.na(d)) 
dat2$ngp = factor(paste(dat2$k, "groups"),
                  levels=paste(dat2$k, "groups"),
                  labels=paste(dat2$k, "groups"))
  
png("power2thru10_plus_q.png", width=1800, height=1800, res=300)
ggplot(data=dat2)+
  geom_line(aes(x=d, y=pow, col=factor(method), linetype=factor(method)), size=.8)+
  facet_wrap(~ngp, nrow=3)+
  theme_classic()+ylab("Power")+
  xlab("Standardized treatment effect in Group 1")+
  scale_color_manual(values=c("gray", 'navyblue', 'red'),
                     labels=expression(Q,Q[LSD], Qmax))+
  scale_linetype(labels=expression(Q,Q[LSD], Qmax))+
  theme(legend.title = element_blank(), 
        legend.position='bottom') 

dev.off()
