#!/usr/bin/env Rscript
# R wavelet analysis. All transforms and smoothing come from WaveletComp/biwavelet.
# Statistical wrapper adds COI-restricted band tests and finite-simulation p-values.
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages(library(WaveletComp))
if(!requireNamespace('biwavelet',quietly=TRUE)) stop('Install biwavelet from CRAN.')
dir.create('results/models',showWarnings=FALSE,recursive=TRUE)
dir.create('results/jobs',showWarnings=FALSE,recursive=TRUE)
dir.create('results/detail',showWarnings=FALSE,recursive=TRUE)
B <- as.integer(Sys.getenv('WAVELET_NSIM','1000'))
CORES <- as.integer(Sys.getenv('WAVELET_CORES','6'))
BASE_SEED <- 24414L
atomic_rds <- function(object,path) {
 tmp=paste0(path,'.tmp_',Sys.getpid());saveRDS(object,tmp)
 invisible(readRDS(tmp))
 if(!file.rename(tmp,path))stop('Could not commit checkpoint: ',path)
}
prep <- readRDS('results/prepared.rds')
slug <- function(x) gsub('[^A-Za-z0-9_]+','_',x)
bands <- function(p) ifelse(p<8,'2-<8',ifelse(p<32,'8-<32','32+'))
band_levels <- c('2-<8','8-<32','32+')
safe_mean <- function(x) if(length(x)) mean(x,na.rm=TRUE) else NA_real_
coi_mask <- function(r,conservative=FALSE) {
 if(conservative) {
  dist=pmin(seq_len(r$nc)-1,r$nc-seq_len(r$nc))
  curve=(2*pi/6)/sqrt(2)*dist
 } else curve=2^r$coi.2[seq_len(r$nc)+2]
 outer(r$Period,curve,'<=')
}
transform <- function(x,y,dj=1/20,win=5,upper=min(128,floor(length(x)/3))) {
 analyze.coherency(data.frame(EPU=x,INF=y),my.pair=c('EPU','INF'),loess.span=0,
  dt=1,dj=dj,lowerPeriod=2,upperPeriod=upper,window.type.t=1,window.type.s=1,
  window.size.t=win,window.size.s=.25,make.pval=FALSE,verbose=FALSE)
}
bandstat <- function(mat,mask,p) {
 vapply(band_levels,function(b){ix=which(bands(p)==b);safe_mean(mat[ix,,drop=FALSE][mask[ix,,drop=FALSE]])},numeric(1))
}
mc_wavelet <- function(x,y,method='AR',seed=BASE_SEED,dj=1/20,win=5,n.sim=B) {
 set.seed(seed);r=transform(x,y,dj,win);mask=coi_mask(r);cons=coi_mask(r,TRUE)
 count=list(co=matrix(0,r$nr,r$nc),x=matrix(0,r$nr,r$nc),y=matrix(0,r$nr,r$nc))
 observed=bandstat(r$Coherence,mask,r$Period);observed_cons=bandstat(r$Coherence,cons,r$Period)
 null=matrix(NA_real_,n.sim,3);null_cons=null
 par=list(AR=list(p=1))
 for(i in seq_len(n.sim)) {
  xs=WaveletComp::SurrogateData(x,method=method,params=par)
  ys=WaveletComp::SurrogateData(y,method=method,params=par)
  s=WaveletComp::WaveletCoherency(xs,ys,dt=1,dj=dj,lowerPeriod=2,upperPeriod=min(128,floor(length(x)/3)),
   window.type.t=1,window.type.s=1,window.size.t=win,window.size.s=.25)
  count$co=count$co+(s$Coherence>=r$Coherence)
  count$x=count$x+(s$Power.x>=r$Power.x);count$y=count$y+(s$Power.y>=r$Power.y)
  null[i,]=bandstat(s$Coherence,mask,r$Period);null_cons[i,]=bandstat(s$Coherence,cons,r$Period)
 }
 r$Coherence.pval=(count$co+1)/(n.sim+1)
 r$Power.x.pval=(count$x+1)/(n.sim+1);r$Power.y.pval=(count$y+1)/(n.sim+1)
 r$band_p=vapply(1:3,function(k)(1+sum(null[,k]>=observed[k]))/(n.sim+1),numeric(1))
 r$band_p_conservative=vapply(1:3,function(k)(1+sum(null_cons[,k]>=observed_cons[k]))/(n.sim+1),numeric(1))
 r$null_band_means=null;r$null_band_means_conservative=null_cons
 r$n.sim=n.sim;r$null_method=method;r$seed=seed
 # BY accounts for arbitrary dependence among time-period cells; BH is also reported.
 for(typ in c('co','x','y')) {
  p=switch(typ,co=r$Coherence.pval,x=r$Power.x.pval,y=r$Power.y.pval)
  for(adj in c('BH','BY')) {
   q=matrix(NA_real_,r$nr,r$nc);q[mask]=p.adjust(p[mask],method=adj)
   r[[paste0('q_',typ,'_',adj)]]=q
  }
 }
 r
}
longest_window <- function(d,columns=c('epu_mom','inf_stl_log')) {
 ok=Reduce('&',lapply(d[columns],is.finite))
 rr=rle(ok);en=cumsum(rr$lengths);candidates=which(rr$values)
 k=candidates[which.max(rr$lengths[candidates])];ids=seq(en[k]-rr$lengths[k]+1,en[k])
 d[ids,]
}
events_at_scale <- function(sig) {
 out=matrix(FALSE,nrow(sig),ncol(sig))
 for(k in seq_len(nrow(sig))) {
  rr=rle(sig[k,]);en=cumsum(rr$lengths)
  for(i in which(rr$values & rr$lengths>=6)) out[k,seq(en[i]-rr$lengths[i]+1,en[i])]=TRUE
 }
 out
}
summarize <- function(r,job,threshold=.70,alpha=.05,conservative=FALSE) {
 mask=coi_mask(r,conservative);sig=mask & r$Coherence.pval<=alpha & r$Coherence>=threshold
 events=events_at_scale(sig);out=list()
 for(b in band_levels) {
  ix=which(bands(r$Period)==b);m=mask[ix,,drop=FALSE];c=r$Coherence[ix,,drop=FALSE]
  s=sig[ix,,drop=FALSE];ev=events[ix,,drop=FALSE];angle=r$sAngle[ix,,drop=FALSE]
  ns=sum(s);nv=sum(m);ang=angle[s];delay=(angle*r$Period[ix]/(2*pi))[s]
  circ=if(ns)Arg(mean(exp(1i*ang))) else NA_real_
  p=if(conservative)r$band_p_conservative[match(b,band_levels)] else r$band_p[match(b,band_levels)]
  out[[b]]=data.frame(series=job$series,scenario=job$scenario,role=job$role,sa=job$sa,
   first=job$data$month[1],last=tail(job$data$month,1),n_months=nrow(job$data),band=b,
   period_max=max(r$Period),valid_cells=nv,mean_coherence=safe_mean(c[m]),
   point_sig_cells=ns,point_sig_pct=100*ns/nv,event_cells=sum(ev),event_pct=100*sum(ev)/nv,
   event_months=sum(colSums(ev)>0),band_mc_p=p,
   epulead_pct=if(ns)100*mean(ang>0) else NA_real_,inflead_pct=if(ns)100*mean(ang<0) else NA_real_,
   inphase_pct=if(ns)100*mean(cos(ang)>=0) else NA_real_,circular_phase=circ,
   phase_concentration=if(ns)Mod(mean(exp(1i*ang))) else NA_real_,median_phase_delay_months=if(ns)median(delay) else NA_real_,
   coherence_BH_sig_cells=sum(r$q_co_BH[ix,,drop=FALSE]<=.05,na.rm=TRUE),
   coherence_BY_sig_cells=sum(r$q_co_BY[ix,,drop=FALSE]<=.05,na.rm=TRUE),
   EPU_power_sig_cells=sum(m&r$Power.x.pval[ix,,drop=FALSE]<=alpha),
   INF_power_sig_cells=sum(m&r$Power.y.pval[ix,,drop=FALSE]<=alpha),
   conservative_coi=conservative,coherence_threshold=threshold,alpha=alpha)
 }
 do.call(rbind,out)
}
export_detail <- function(r,job) {
 m=coi_mask(r);grid=expand.grid(period_idx=seq_len(r$nr),time_idx=seq_len(r$nc))
 grid$month=job$data$month[grid$time_idx];grid$period_months=r$Period[grid$period_idx]
 grid$band=bands(grid$period_months);grid$in_coi=as.vector(m);grid$coherence=as.vector(r$Coherence)
 grid$coherence_p=as.vector(r$Coherence.pval);grid$coherence_q_BY=as.vector(r$q_co_BY)
 grid$phase_x_minus_y=as.vector(r$sAngle);grid$epu_power=as.vector(r$Power.x);grid$inf_power=as.vector(r$Power.y)
 grid$epu_power_p=as.vector(r$Power.x.pval);grid$inf_power_p=as.vector(r$Power.y.pval)
 grid$significant=as.vector(m & r$Coherence.pval<=.05&r$Coherence>=.7)
 con=gzfile(paste0('results/detail/',slug(job$id),'.csv.gz'),'wt');write.csv(grid,con,row.names=FALSE);close(con)
 # Events at each scale: sums are cells, not independent tests or unique calendar months.
 events=list()
 for(k in seq_len(r$nr)) {
  s=m[k,]&r$Coherence.pval[k,]<=.05&r$Coherence[k,]>=.7;rr=rle(s);en=cumsum(rr$lengths)
  for(i in which(rr$values & rr$lengths>=6)) {
   ids=seq(en[i]-rr$lengths[i]+1,en[i]);a=r$sAngle[k,ids]
   events[[length(events)+1]]=data.frame(series=job$series,scenario=job$scenario,period_months=r$Period[k],
    band=bands(r$Period[k]),start=job$data$month[ids[1]],end=job$data$month[tail(ids,1)],
    duration_months=length(ids),mean_coherence=mean(r$Coherence[k,ids]),mean_p=mean(r$Coherence.pval[k,ids]),
    circular_phase=Arg(mean(exp(1i*a))),phase_concentration=Mod(mean(exp(1i*a))))
  }
 }
 if(length(events))write.csv(do.call(rbind,events),paste0('results/detail/',slug(job$id),'_events.csv'),row.names=FALSE)
}
# Validate the statistical wrapper against the package's own Monte Carlo engine.
if(!file.exists('results/mc_validation.csv')) {
 set.seed(77);x=rnorm(90);y=rnorm(90)
 set.seed(123);native=WaveletComp::wc(x,y,upperPeriod=30,dj=1/20,make.pval=TRUE,method='AR',params=list(AR=list(p=1)),n.sim=20)
 # mc_wavelet uses the same upperPeriod=min(128,N/3) and exact package surrogates.
 ours=mc_wavelet(x,y,seed=123,n.sim=20)
 delta=max(abs(native$Coherence.pval*20-(ours$Coherence.pval*21-1)))
 stopifnot(delta<1e-10)
 write.csv(data.frame(check='MC counts equal native WaveletComp',max_count_error=delta,B=20), 'results/mc_validation.csv',row.names=FALSE)
 # Verify phase sign on a known lead: x=sin(t), y=sin(t-3) at a 24-month period.
 tt=1:240;r=transform(sin(2*pi*tt/24),sin(2*pi*(tt-3)/24));k=which.min(abs(r$Period-24))
 a=Arg(mean(exp(1i*r$sAngle[k,61:180])));stopifnot(a>0)
 write.csv(data.frame(true_EPU_lead_months=3,period=r$Period[k],estimated_phase=a,
 estimated_delay=a*r$Period[k]/(2*pi)), 'results/phase_validation.csv',row.names=FALSE)
}
jobs=list();attempts=list();windows=list()
add <- function(s,scenario,role,sa,d,engine='WaveletComp',method='AR',dj=1/20,win=5) {
 id=paste(s,scenario,sep='__')
 jobs[[length(jobs)+1]] <<- list(id=id,series=s,scenario=scenario,role=role,sa=sa,data=d,
  engine=engine,method=method,dj=dj,win=win,seed=BASE_SEED+length(jobs)*101L)
}
for(s in c(prep$retained,'SCMP China','Mainland China')) {
 d=prep$series[[s]];col=paste0('inf_',prep$chosen_sa)
 ok=is.finite(d$epu_mom)&is.finite(d[[col]])
 d=d[seq(which(ok)[1],tail(which(ok),1)),]
 error=tryCatch({z=transform(d$epu_mom,d[[col]]);if(any(!is.finite(z$Coherence)))stop('Non-finite wavelet output');NA_character_},error=function(e)conditionMessage(e))
 attempts[[s]]=data.frame(series=s,full_first=d$month[1],full_last=tail(d$month,1),
  full_n=nrow(d),missing_pairs=sum(!ok[match(d$month,prep$series[[s]]$month)]),full_success=is.na(error),error=error)
 role=if(s%in%c('SCMP China','Mainland China'))'proxy_sensitivity' else if(s=='US'&&!is.na(error))'provisional_US' else 'main'
 if(!is.na(error))d=longest_window(d,c('epu_mom',col))
 windows[[s]]=d
 scen=if(s=='US'&&!is.na(error))'observed_contiguous' else if(s=='SCMP China'&&!is.na(error))'proxy_contiguous' else 'main'
 add(s,scen,role,prep$chosen_sa,d)
 if(s%in%prep$retained) {
  for(sa in setdiff(c('classical','stl_log','stl_level','nsa'),prep$chosen_sa)) {
   alt=longest_window(d,c('epu_mom',paste0('inf_',sa)))
   add(s,paste0('SA_',sa),'robustness',sa,alt)
  }
  add(s,'Fourier_null','robustness',prep$chosen_sa,d,method='Fourier.rand')
  add(s,'smooth_time_9','parameter_sensitivity',prep$chosen_sa,d,engine='deterministic',win=9)
  add(s,'grid_12voices','parameter_sensitivity',prep$chosen_sa,d,engine='deterministic',dj=1/12)
  add(s,'biwavelet','library_sensitivity',prep$chosen_sa,d,engine='biwavelet')
 }
}
# China proxy comparison uses the exact same continuous time window for all three indexes.
cn=c('China','SCMP China','Mainland China');common=Reduce(intersect,lapply(windows[cn],function(d)d$month))
for(s in cn) {
 d=windows[[s]];d=d[match(common,d$month),];stopifnot(all(diff(as.integer(d$date))>=28 & diff(as.integer(d$date))<=31))
 add(s,'China_common_window','proxy_matched',prep$chosen_sa,d)
}
write.csv(do.call(rbind,attempts),'results/full_sample_attempts.csv',row.names=FALSE,na='')
job_manifest=do.call(rbind,lapply(jobs,function(j)data.frame(id=j$id,series=j$series,scenario=j$scenario,role=j$role,
 sa=j$sa,engine=j$engine,null=j$method,dj=j$dj,window_t=j$win,seed=j$seed,n=nrow(j$data),
 first=j$data$month[1],last=tail(j$data$month,1),n_sim=if(j$engine=='deterministic')0 else B)))
write.csv(job_manifest,'results/job_manifest.csv',row.names=FALSE)
run_job <- function(j) {
 path=paste0('results/jobs/',slug(j$id),'.rds')
 modelpath=paste0('results/models/',slug(j$id),'.rds')
 cached=if(file.exists(modelpath))tryCatch(readRDS(modelpath),error=function(e)NULL) else NULL
 same=FALSE
 if(!is.null(cached)) {
  oj=cached$job
  same=identical(oj$data$month,j$data$month)&&
   isTRUE(all.equal(oj$data$epu_mom,j$data$epu_mom,tolerance=0))&&
   isTRUE(all.equal(oj$data[[paste0('inf_',oj$sa)]],j$data[[paste0('inf_',j$sa)]],tolerance=0))&&
   identical(oj$engine,j$engine)&&identical(oj$method,j$method)&&identical(oj$dj,j$dj)&&identical(oj$win,j$win)
  if(j$engine=='WaveletComp')same=same&&identical(cached$model$n.sim,B)
 }
 if(file.exists(path)&&same) { old=readRDS(path);if(old$status=='ok')return(old) }
 t0=proc.time()[3];cat('START',j$id,'n=',nrow(j$data),'\n');flush.console()
 ans=tryCatch({
  x=j$data$epu_mom;y=j$data[[paste0('inf_',j$sa)]];stopifnot(all(is.finite(x)),all(is.finite(y)))
  if(j$engine=='deterministic') {
   r=transform(x,y,j$dj,j$win);m=coi_mask(r)
   out=data.frame(series=j$series,scenario=j$scenario,band=band_levels,
    mean_coherence=bandstat(r$Coherence,m,r$Period),n_months=length(x))
   atomic_rds(list(model=r,job=j),paste0('results/models/',slug(j$id),'.rds'))
   list(status='ok',id=j$id,parameters=out)
  } else if(j$engine=='biwavelet') {
   set.seed(j$seed);ff=4*pi/(6+sqrt(38))
   r=biwavelet::wtc(cbind(seq_along(x),as.numeric(scale(x))),cbind(seq_along(y),as.numeric(scale(y))),
    dj=1/20,s0=2/ff,max.scale=min(128,floor(length(x)/3))/ff,mother='morlet',nrands=B,quiet=TRUE)
   m=outer(r$period,as.vector(r$coi),'<=');sig=m&r$rsq>=r$signif&r$rsq>=.70
   out=do.call(rbind,lapply(band_levels,function(b){k=which(bands(r$period)==b);mm=m[k,,drop=FALSE]
    data.frame(series=j$series,scenario=j$scenario,band=b,n_months=length(x),valid_cells=sum(mm),
     mean_coherence=safe_mean(r$rsq[k,,drop=FALSE][mm]),point_sig_pct=100*sum(sig[k,,drop=FALSE])/sum(mm))}))
   atomic_rds(list(model=r,job=j),paste0('results/models/',slug(j$id),'.rds'))
   list(status='ok',id=j$id,library=out)
  } else {
   if(same)r=cached$model else r=mc_wavelet(x,y,j$method,j$seed,j$dj,j$win)
   r$dates=j$data$month
   stopifnot(all(r$Coherence>=-1e-8),all(r$Coherence<=1+1e-8),length(r$dates)==r$nc)
   atomic_rds(list(model=r,job=j),paste0('results/models/',slug(j$id),'.rds'))
   main=summarize(r,j);cons=summarize(r,j,conservative=TRUE)
   threshold=NULL
   if(j$role%in%c('main','provisional_US')) {
    threshold=do.call(rbind,lapply(c(.6,.7,.8),function(th)do.call(rbind,lapply(c(.01,.05,.10),function(al)summarize(r,j,th,al)))))
   }
   if(j$role%in%c('main','provisional_US','proxy_matched','proxy_sensitivity'))export_detail(r,j)
   list(status='ok',id=j$id,summary=main,conservative=cons,threshold=threshold)
  }
 },error=function(e)list(status='error',id=j$id,error=conditionMessage(e)))
 ans$elapsed=proc.time()[3]-t0;atomic_rds(ans,path)
 cat('DONE',j$id,ans$status,round(ans$elapsed,1),'sec\n');flush.console();ans
}
# Main and Chinese proxy results are checkpointed first. Independent jobs use R's process pool.
priority=c('main','provisional_US','proxy_sensitivity','proxy_matched','robustness','library_sensitivity','parameter_sensitivity')
jobs=jobs[order(match(vapply(jobs,`[[`,character(1),'role'),priority))]
all=if(.Platform$OS.type=='windows')lapply(jobs,run_job) else
 parallel::mclapply(jobs,run_job,mc.cores=CORES,mc.preschedule=FALSE,mc.set.seed=FALSE)
combine <- function(field) {
 vals=lapply(all,function(a)a[[field]]);vals=Filter(Negate(is.null),vals)
 if(length(vals))do.call(rbind,vals) else data.frame()
}
summary=combine('summary');conservative=combine('conservative')
# Holm family: 18 country x 3 band tests, including the explicitly provisional US contiguous window.
main=summary$role%in%c('main','provisional_US')
summary$band_mc_p_Holm=NA_real_;summary$band_mc_p_Holm[main]=p.adjust(summary$band_mc_p[main],method='holm')
for(sc in unique(summary$scenario[!main])) {
 k=which(!main&summary$scenario==sc);summary$band_mc_p_Holm[k]=p.adjust(summary$band_mc_p[k],method='holm')
}
conservative$band_mc_p_Holm=NA_real_;k=conservative$role%in%c('main','provisional_US')
conservative$band_mc_p_Holm[k]=p.adjust(conservative$band_mc_p[k],method='holm')
write.csv(summary,'results/all_band_summary.csv',row.names=FALSE)
write.csv(summary[main,],'results/main_band_summary.csv',row.names=FALSE)
write.csv(conservative,'results/conservative_coi_summary.csv',row.names=FALSE)
write.csv(combine('threshold'),'results/threshold_sensitivity.csv',row.names=FALSE)
write.csv(combine('library'),'results/biwavelet_summary.csv',row.names=FALSE)
write.csv(combine('parameters'),'results/parameter_sensitivity.csv',row.names=FALSE)
status=do.call(rbind,lapply(all,function(a)data.frame(id=a$id,status=a$status,elapsed_sec=a$elapsed,error=if(is.null(a$error))'' else a$error)))
write.csv(status,'results/run_status.csv',row.names=FALSE)
writeLines(capture.output(sessionInfo()),'results/sessionInfo.txt')
if(any(status$status!='ok'))stop('Some jobs failed. See run_status.csv and checkpoints.')
cat('ALL',nrow(status),'JOBS COMPLETED; MC draws per inference job:',B,'\n')
