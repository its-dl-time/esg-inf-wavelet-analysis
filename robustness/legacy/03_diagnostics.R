#!/usr/bin/env Rscript
# Diagnostics do not determine which wavelet findings are retained.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(urca))
p=readRDS('results/prepared.rds')
longest <- function(d,col) {
 ok=is.finite(d$epu_mom)&is.finite(d[[col]]);rr=rle(ok);en=cumsum(rr$lengths)
 can=which(rr$values);k=can[which.max(rr$lengths[can])]
 d[seq(en[k]-rr$lengths[k]+1,en[k]),]
}
desc=list();tests=list()
for(s in c(p$retained,'SCMP China','Mainland China')) {
 col=paste0('inf_',p$chosen_sa);d=longest(p$series[[s]],col)
 for(v in c('epu_mom',col)) {
  z=d[[v]];mu=mean(z);ss=sd(z);zz=(z-mu)/sqrt(mean((z-mu)^2))
  desc[[length(desc)+1]]=data.frame(series=s,variable=v,first=d$month[1],last=tail(d$month,1),
   n=length(z),mean=mu,sd=ss,min=min(z),max=max(z),skewness=mean(zz^3),kurtosis=mean(zz^4))
  ad=ur.df(z,type='drift',lags=12,selectlags='AIC');crit=as.numeric(ad@cval['tau2','5pct']);st=as.numeric(ad@teststat[1,'tau2'])
  kp=ur.kpss(z,type='mu',lags='long');kst=as.numeric(kp@teststat);kcrit=as.numeric(kp@cval[1,'5pct'])
  za=ur.za(z,model='intercept',lag=ad@lags);zst=as.numeric(za@teststat);zcrit=as.numeric(za@cval[2])
  tests[[length(tests)+1]]=data.frame(series=s,variable=v,n=length(z),adf_lags=ad@lags,
   adf_stat=st,adf_critical_5pct=crit,adf_reject_unitroot=st<crit,
   kpss_stat=kst,kpss_critical_5pct=kcrit,kpss_reject_stationarity=kst>kcrit,
   ZA_stat=zst,ZA_critical_5pct=zcrit,ZA_reject_unitroot=zst<zcrit,
   ZA_break_month=d$month[za@bpoint])
 }
}
write.csv(do.call(rbind,desc),'results/descriptive_statistics.csv',row.names=FALSE)
write.csv(do.call(rbind,tests),'results/stationarity_tests.csv',row.names=FALSE)
cat('Descriptives and ADF/KPSS/ZA completed for',length(p$retained)+2,'series.\n')
