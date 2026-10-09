#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
src=readLines('02_run_wavelet.R');cut=grep('^# Validate the statistical wrapper',src)[1]
eval(parse(text=src[seq_len(cut-1)]))
regrid <- function(mat,period,to_period,columns) {
 vapply(columns,function(t)approx(log2(period),mat[,t],xout=log2(to_period),rule=1,ties='ordered')$y,numeric(length(to_period)))
}
main=list()
for(f in list.files('results/models',pattern='\\.rds$',full.names=TRUE)) {
 o=readRDS(f)
 if(o$job$role%in%c('main','provisional_US'))main[[o$job$series]]=o
}
out=list()
for(f in list.files('results/models',pattern='\\.rds$',full.names=TRUE)) {
 o=readRDS(f);j=o$job
 if(!j$role%in%c('robustness','library_sensitivity','parameter_sensitivity','stationarity_robustness'))next
 if(!j$series%in%names(main))next
 a=main[[j$series]];r=a$model;s=o$model
 dates=intersect(a$job$data$month,j$data$month);ia=match(dates,a$job$data$month);ib=match(dates,j$data$month)
 isbi=j$role=='library_sensitivity'
 pp=if(isbi)s$period else s$Period
 mat=if(isbi)s$rsq else s$Coherence
 aa=r$Coherence[,ia,drop=FALSE];bb=regrid(mat,pp,r$Period,ib)
 if(isbi) {
  am=coi_mask(r,TRUE)[,ia,drop=FALSE]
  bm=outer(r$Period,as.vector(s$coi)[ib],'<=')
  critical=regrid(s$signif,pp,r$Period,ib)
  bs=bb>=critical & bb>=.7
 } else {
  am=coi_mask(r)[,ia,drop=FALSE]
  curve=2^s$coi.2[seq_len(s$nc)+2];bm=outer(r$Period,curve[ib],'<=')
  if(is.null(s$Coherence.pval))bs=matrix(NA,nrow(bb),ncol(bb)) else
   bs=regrid(s$Coherence.pval,pp,r$Period,ib)<=.05 & bb>=.7
 }
 mm=am&bm&is.finite(aa)&is.finite(bb)
 asig=r$Coherence.pval[,ia,drop=FALSE]<=.05&aa>=.7
 for(b in band_levels) {
  k=which(bands(r$Period)==b);m=mm[k,,drop=FALSE]
  av=aa[k,,drop=FALSE][m];bv=bb[k,,drop=FALSE][m]
  sa=asig[k,,drop=FALSE][m];sb=bs[k,,drop=FALSE][m]
  valid_sig=all(!is.na(sb));union=sum(sa|sb,na.rm=TRUE)
  out[[length(out)+1]]=data.frame(series=j$series,scenario=j$scenario,band=b,
   first=dates[1],last=tail(dates,1),n_common_months=length(dates),shared_COI_cells=sum(m),
   coherence_map_correlation=cor(av,bv),mean_abs_coherence_difference=mean(abs(av-bv)),
   main_mean_coherence=mean(av),alternative_mean_coherence=mean(bv),
   main_point_sig_pct=100*mean(sa),alternative_point_sig_pct=if(valid_sig)100*mean(sb) else NA_real_,
   significant_region_Jaccard=if(valid_sig&&union>0)sum(sa&sb)/union else NA_real_,
   mask=if(isbi)'shared conservative/native biwavelet COI' else 'intersection of package COIs')
 }
}
z=do.call(rbind,out);write.csv(z,'results/robustness_comparison.csv',row.names=FALSE)
# Endpoint revisions are an additional diagnostic; choice of SA was made before wavelet analysis.
cat('Matched robustness comparisons:',nrow(z),'rows.\n')
