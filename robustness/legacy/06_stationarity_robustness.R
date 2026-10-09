#!/usr/bin/env Rscript
# Sensitivity to low-frequency drift when KPSS rejects level stationarity.
options(stringsAsFactors=FALSE)
src=readLines('02_run_wavelet.R');cut=grep('^# Validate the statistical wrapper',src)[1]
eval(parse(text=src[seq_len(cut-1)]))
tests=read.csv('results/stationarity_tests.csv');col=paste0('inf_',prep$chosen_sa)
affected=tests$series[tests$variable==col&tests$kpss_reject_stationarity&tests$series%in%prep$retained]
run_one <- function(s) {
 id=paste0(s,'__INF_drift_ns3');path=paste0('results/models/',slug(id),'.rds')
 d=longest_window(prep$series[[s]],c('epu_mom',col));tt=seq_len(nrow(d))
 y=residuals(lm(d[[col]]~splines::ns(tt,df=3)));d[[col]]=as.numeric(y)
 j=list(id=id,series=s,scenario='INF_drift_ns3',role='stationarity_robustness',sa=prep$chosen_sa,data=d)
 if(file.exists(path))r=readRDS(path)$model else {
  r=mc_wavelet(d$epu_mom,y,seed=BASE_SEED+9000+match(s,affected));r$dates=d$month
  atomic_rds(list(model=r,job=j),path)
 }
 cat('DONE drift robustness:',s,'\n');summarize(r,j)
}
out=if(.Platform$OS.type=='windows')lapply(affected,run_one) else parallel::mclapply(affected,run_one,mc.cores=2)
z=do.call(rbind,out);z$band_mc_p_Holm=p.adjust(z$band_mc_p,'holm')
write.csv(z,'results/stationarity_robustness.csv',row.names=FALSE)
