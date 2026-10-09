#!/usr/bin/env Rscript
# Run from the directory containing this script. All preprocessing is in R.
options(stringsAsFactors=FALSE, warn=1)
dir.create('results',showWarnings=FALSE)
P <- read.csv('data/panel.csv',na.strings=c('', 'NA'))
C <- read.csv('data/cpi_raw.csv',na.strings=c('', 'NA'))
iso <- c(Brazil='BR',Canada='CA',Chile='CL',China='CN',France='FR',Germany='DE',
 Greece='GR',India='IN',Ireland='IE',Italy='IT',Japan='JP',Korea='KR',Mexico='MX',
 Russia='RU',Spain='ES',Singapore='SG',UK='GB',US='US',Australia='AU',Pakistan='PK',
 'SCMP China'='CN','Mainland China'='CN')
monthid <- function(x) as.integer(substr(x,1,4))*12L+as.integer(substr(x,6,7))
mom <- function(x) { x[!is.finite(x)|x<=0]=NA_real_;c(NA_real_,100*diff(log(x))) }
bridge_short <- function(y,maxgap=2) {
 rr <- rle(is.na(y)); ends <- cumsum(rr$lengths); starts <- ends-rr$lengths+1
 for(k in which(rr$values & rr$lengths<=maxgap)) {
  a=starts[k]; b=ends[k]
  if(a>1 && b<length(y) && is.finite(y[a-1]) && is.finite(y[b+1]))
   y[a:b] <- seq(y[a-1],y[b+1],length.out=b-a+3)[2:(b-a+2)]
 }
 y
}
stl_sa <- function(x,dates,logscale=TRUE) {
 work <- bridge_short(if(logscale)log(x) else x)
 ans <- rep(NA_real_,length(x)); rr=rle(is.finite(work)); en=cumsum(rr$lengths)
 for(k in which(rr$values & rr$lengths>=36)) {
  ids=seq(en[k]-rr$lengths[k]+1,en[k]); first=dates[ids[1]]
  fit=stats::stl(ts(work[ids],frequency=12,start=c(as.integer(substr(first,1,4)),as.integer(substr(first,6,7)))),
   s.window=13,s.degree=1,t.degree=1,l.degree=1,s.jump=1,t.jump=1,l.jump=1,robust=TRUE)
  z=work[ids]-fit$time.series[,'seasonal'];ans[ids]=if(logscale)exp(z) else z
 }
 ans[is.na(x)] <- NA_real_; ans
}
seasonality_diag <- function(y,dates) {
 ok=is.finite(y); y=y[ok]; dates=dates[ok]; n=length(y)
 m=factor(substr(dates,6,7)); tt=monthid(dates)-monthid(dates[1])
 # Joint calendar-month test after a flexible trend, Newey-West covariance (12 lags).
 fit=lm(y~splines::ns(tt,df=6)+m); X=model.matrix(fit); u=residuals(fit)
 A=X*as.numeric(u); meat=crossprod(A)
 for(l in 1:min(12,n-1)) {
  g=crossprod(A[(l+1):n,,drop=FALSE],A[1:(n-l),,drop=FALSE]);
  meat=meat+(1-l/13)*(g+t(g))
 }
 bread=solve(crossprod(X)); V=bread%*%meat%*%bread*n/(n-ncol(X))
 ix=grep('^m',colnames(X)); b=coef(fit)[ix]
 W=tryCatch(as.numeric(t(b)%*%solve(V[ix,ix,drop=FALSE],b)),error=function(e) NA_real_)
 # Partial R2 quantifies residual seasonality; lag-12 ACF is descriptive only.
 f0=lm(y~splines::ns(tt,df=6)); pr2=1-sum(residuals(fit)^2)/sum(residuals(f0)^2)
 c(month_wald=W,month_p=pchisq(W,length(ix),lower.tail=FALSE),month_partial_r2=pr2,
  acf12=acf(y,lag.max=12,plot=FALSE)$acf[13])
}
cpi <- list(); diagnostics <- list(); revisions <- list()
for(code in unique(C$iso2)) {
 d=C[C$iso2==code,]; d=d[order(d$month),]
 dates=format(seq(as.Date(paste0(d$month[1],'-01')),as.Date(paste0(tail(d$month,1),'-01')),by='month'),'%Y-%m')
 d=d[match(dates,d$month),]; d$month=dates;d$iso2=code
 d$cpi_stl_log=stl_sa(d$cpi_nsa,dates,TRUE)
 d$cpi_stl_level=stl_sa(d$cpi_nsa,dates,FALSE)
 for(method in c('classical','stl_log','stl_level','nsa')) {
  y=mom(d[[paste0('cpi_',method)]]);d[[paste0('inf_',method)]]=y
 }
 cpi[[code]]=d
}
N <- list(); coverage <- list(); stationarity <- list()
countries=setdiff(unique(P$series),c('GEPU_current','GEPU_ppp','SCMP China','Mainland China'))
for(s in c(countries,'SCMP China','Mainland China')) {
 a=P[P$series==s,];a=a[order(a$month),];code=unname(iso[s]);d=cpi[[code]]
 a$date=as.Date(paste0(a$month,'-01'));a$epu_mom=mom(a$epu_level)
 # Any skipped EPU month is a real gap, never an adjacent observation.
 a$epu_mom[c(FALSE,diff(monthid(a$month))!=1)]=NA_real_
 j=match(a$month,d$month)
 for(m in c('classical','stl_log','stl_level','nsa')) {
  a[[paste0('cpi_',m)]]=d[[paste0('cpi_',m)]][j]
  a[[paste0('inf_',m)]]=d[[paste0('inf_',m)]][j]
 }
 a$epu_mom[!is.finite(a$epu_mom)]=NA_real_
 N[[s]]=a
 n=sum(is.finite(a$epu_mom)&is.finite(a$inf_classical))
 coverage[[s]]=data.frame(series=s,n_pairs=n,first=a$month[1],last=tail(a$month,1),iso2=code)
 if(s %in% countries && !s %in% c('Australia','Pakistan')) {
  for(m in c('classical','stl_log','stl_level','nsa')) {
   z=a[[paste0('inf_',m)]]; v=seasonality_diag(z,a$month)
   diagnostics[[paste(s,m)]]=data.frame(series=s,method=m,n=sum(is.finite(z)),t(v))
  }
  # End-point revision: compare the last 24 overlapping growth rates when 12 latest months are withheld.
  dd=d[seq_len(nrow(d)-12),];
  for(m in c('stl_log','stl_level')) {
   z=mom(stl_sa(dd$cpi_nsa,dd$month,m=='stl_log'))
   full=d[[paste0('inf_',m)]][seq_len(nrow(dd))];k=tail(which(is.finite(z)&is.finite(full)),24)
   revisions[[paste(s,m)]]=data.frame(series=s,method=m,revision_rmse_pp=sqrt(mean((z[k]-full[k])^2)),n=length(k))
  }
 }
}
coverage=do.call(rbind,coverage);rownames(coverage)=NULL
ranked=coverage[coverage$series%in%countries,];ranked=ranked[order(ranked$n_pairs,ranked$series),]
excluded=head(ranked$series,2);stopifnot(identical(excluded,c('Australia','Pakistan')))
retained=setdiff(countries,excluded)
diagnostics=do.call(rbind,diagnostics);rownames(diagnostics)=NULL
# Decide before wavelet results: residual calendar seasonality, not statistical significance of coherence.
scores=do.call(rbind,lapply(c('classical','stl_log','stl_level'),function(m){
 z=diagnostics[diagnostics$method==m,]
 data.frame(method=m,rejecting_countries=sum(z$month_p<.05,na.rm=TRUE),
 median_partial_r2=median(z$month_partial_r2),median_abs_acf12=median(abs(z$acf12)))
}))
scores=scores[order(scores$rejecting_countries,scores$median_partial_r2),]
chosen=scores$method[1]
write.csv(coverage,'results/coverage.csv',row.names=FALSE)
write.csv(diagnostics,'results/seasonality_diagnostics.csv',row.names=FALSE)
write.csv(scores,'results/seasonality_scores.csv',row.names=FALSE)
write.csv(do.call(rbind,revisions),'results/endpoint_revisions.csv',row.names=FALSE)
processed=do.call(rbind,N);rownames(processed)=NULL
write.csv(processed,'results/processed_panel.csv',row.names=FALSE,na='')
saveRDS(list(series=N,cpi=cpi,retained=retained,excluded=excluded,chosen_sa=chosen,scores=scores), 'results/prepared.rds')
cat('EXCLUDED:',paste(excluded,collapse=', '),'\nCHOSEN SA:',chosen,'\n');print(scores)
