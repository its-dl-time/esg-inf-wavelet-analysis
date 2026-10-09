#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
# Load shared functions without executing the job scheduler.
src=readLines('02_run_wavelet.R');cut=grep('^# Validate the statistical wrapper',src)[1]
eval(parse(text=src[seq_len(cut-1)]))
dir.create('figures',showWarnings=FALSE)
files=list.files('results/models',pattern='\\.rds$',full.names=TRUE)
files=files[file.exists(sub('/models/','/jobs/',files,fixed=TRUE))]
objects=lapply(files,readRDS)
main=Filter(function(o)o$job$role%in%c('main','provisional_US'),objects)
main=main[order(vapply(main,function(o)o$job$series,character(1)))]
if(length(main)!=18)stop('Expected 18 main/provisional models before making final figures.')
pal=viridisLite::viridis(100)
draw_map <- function(r,j,type='coherence',compact=FALSE,significantOnly=FALSE,significanceMask=NULL,conservativeLine=TRUE) {
 pp=log2(r$Period);xx=seq_len(r$nc);mask=coi_mask(r)
 mat=if(type=='coherence')r$Coherence else if(type=='EPU')log2(r$Power.x+1e-8) else log2(r$Power.y+1e-8)
 if(type=='coherence')zlim=c(0,1) else zlim=range(mat[mask],finite=TRUE)
 sig=if(is.null(significanceMask))r$Coherence.pval<=.05 else significanceMask
 if(significantOnly)mat[!mask|!sig|r$Coherence<.70]=NA_real_
 image(xx,pp,t(mat),col=pal,zlim=zlim,axes=FALSE,xlab='',ylab='',useRaster=TRUE)
 at=unique(round(seq(1,r$nc,length.out=if(compact)4 else 7)))
 axis(1,at=at,labels=substr(j$data$month[at],1,4),cex.axis=if(compact).7 else .9)
 yt=c(2,4,8,16,32,64,128);yt=yt[yt<=max(r$Period)+1e-7]
 axis(2,at=log2(yt),labels=yt,las=1,cex.axis=if(compact).65 else .9)
 if(!compact){mtext('Period (months)',side=2,line=2.5);mtext('Calendar year',side=1,line=2.5)}
 # Gray shading uses the exact package COI; dashed line is a stricter e-folding boundary.
 curve=2^r$coi.2[seq_len(r$nc)+2];top=max(pp);bottom=min(pp)
 cclip=pmin(top,pmax(bottom,log2(curve)))
 polygon(c(xx,rev(xx)),c(cclip,rep(top,r$nc)),col=adjustcolor('grey85',alpha.f=.65),border=NA)
 lines(xx,cclip,col='grey25',lwd=1)
 if(conservativeLine) {
  conservative=(2*pi/6)/sqrt(2)*pmin(xx-1,r$nc-xx)
  lines(xx,pmin(top,pmax(bottom,log2(pmax(conservative,1e-8)))),col='white',lty=2,lwd=1)
 }
 pv=if(type=='coherence')r$Coherence.pval else if(type=='EPU')r$Power.x.pval else r$Power.y.pval
 flag=if(type=='coherence'&&!is.null(significanceMask))significanceMask&mask else (pv<=.05)&mask
 if(any(flag)&&any(!flag))contour(xx,pp,t(flag*1),levels=.5,drawlabels=FALSE,add=TRUE,col='black',lwd=.65)
 if(type=='coherence'&&!compact) {
  tstep=max(1,floor(r$nc/27));sstep=max(1,floor(r$nr/15))
  for(k in seq(1,r$nr,by=sstep))for(t in seq(1,r$nc,by=tstep))if(flag[k,t]&&r$Coherence[k,t]>=.7) {
   a=r$sAngle[k,t];dx=.01*r$nc*cos(a);dy=.06*sin(a)
   arrows(t-dx/2,pp[k]-dy/2,t+dx/2,pp[k]+dy/2,length=.035,col='black',lwd=.75)
  }
 }
 abline(h=log2(c(8,32)),col=adjustcolor('white',.65),lty=3)
 if(!compact) {
  labs=if(type=='coherence')c('R² = 0','R² = 0.5','R² = 1') else paste('log2 power',round(c(zlim[1],mean(zlim),zlim[2]),2))
  legend('topright',legend=labs,fill=pal[c(1,50,100)],bg='white',cex=.65,bty='o',horiz=TRUE)
 }
 box();title(main=paste(j$series,if(j$role=='provisional_US')'*' else '',if(compact)'' else paste('|',type)),cex.main=if(compact).95 else 1.2)
}
# Publication vector PDF plus matching PNG for each country and the China matched-window sensitivity.
plot_models=c(main,Filter(function(o)o$job$role=='proxy_matched',objects))
for(o in plot_models)for(type in c('coherence','EPU','INF')) {
 f=paste0('figures/',slug(o$job$id),'_',type)
 for(ext in c('pdf','png')) {
  if(ext=='pdf')pdf(paste0(f,'.pdf'),width=11,height=5.5) else png(paste0(f,'.png'),width=1980,height=990,res=180)
  par(mar=c(4.2,4.2,3,1));draw_map(o$model,o$job,type)
  mtext('Black contour: pointwise AR(1) p <= 0.05; gray: package edge region; dashed white: conservative COI.',side=1,line=3.3,cex=.7)
  dev.off()
 }
}
for(ext in c('pdf','png')) {
 if(ext=='pdf')pdf('figures/all_country_coherence.pdf',width=14,height=17) else png('figures/all_country_coherence.png',width=2100,height=2550,res=150)
 par(mfrow=c(6,3),mar=c(2.5,2.3,2,0.5),oma=c(2,1,3,1))
 for(o in main)draw_map(o$model,o$job,compact=TRUE)
 mtext('EPU-inflation wavelet coherence | 18 countries',outer=TRUE,side=3,line=1,font=2,cex=1.3)
 mtext('Colors: 0 (purple) to 1 (yellow). * US ends September 2025; full sample awaits missing CPI.',outer=TRUE,side=1,line=.5,cex=.9)
 dev.off()
}
png('figures/all_country_significant_coherence.png',width=2100,height=2550,res=150)
par(mfrow=c(6,3),mar=c(2.5,2.3,2,0.5),oma=c(2,1,3,1))
for(o in main)draw_map(o$model,o$job,compact=TRUE,significantOnly=TRUE)
mtext('Local coherence: p <= 0.05 and R2 >= 0.70 | 18 countries',outer=TRUE,side=3,line=1,font=2,cex=1.2)
mtext('Pointwise significance before multiplicity correction. * US sample ends September 2025.',outer=TRUE,side=1,line=.5,cex=.9)
dev.off()
biobjs=Filter(function(o)o$job$role=='library_sensitivity',objects)
biobjs=biobjs[order(vapply(biobjs,function(o)o$job$series,character(1)))]
as_bi_plot <- function(o) {
 s=o$model
 list(Coherence=s$rsq,Period=s$period,nc=length(s$t),nr=length(s$period),
  coi.2=c(0,0,log2(as.vector(s$coi)),0,0),sAngle=s$phase)
}
for(o in biobjs) {
 r=as_bi_plot(o);j=o$job;if(j$series=='US')j$role='provisional_US'
 for(ext in c('pdf','png')) {
  f=paste0('figures/',slug(j$id),'_coherence.',ext)
  if(ext=='pdf')pdf(f,width=11,height=5.5) else png(f,width=1980,height=990,res=180)
  par(mar=c(4.2,4.2,3,1))
  draw_map(r,j,significanceMask=o$model$rsq>=o$model$signif,conservativeLine=FALSE)
  mtext('biwavelet; black contour: native AR(1) 95% threshold; gray: native COI edge region.',side=1,line=3.3,cex=.75)
  dev.off()
 }
}
png('figures/all_country_biwavelet_coherence.png',width=2100,height=2550,res=150)
par(mfrow=c(6,3),mar=c(2.5,2.3,2,0.5),oma=c(2,1,3,1))
for(o in biobjs){j=o$job;if(j$series=='US')j$role='provisional_US';draw_map(as_bi_plot(o),j,compact=TRUE,
 significanceMask=o$model$rsq>=o$model$signif,conservativeLine=FALSE)}
mtext('biwavelet coherence | independently calibrated library sensitivity',outer=TRUE,side=3,line=1,font=2,cex=1.2)
mtext('Colors: 0 (purple) to 1 (yellow). Black: AR(1) 95% threshold. * US ends September 2025.',outer=TRUE,side=1,line=.5,cex=.85)
dev.off()
if(Sys.getenv('WAVELET_PREVIEW','0')=='1')quit(save='no',status=0)
sum=read.csv('results/main_band_summary.csv');all=read.csv('results/all_band_summary.csv')
tableVI=list();tableVII=sum[,c('series','first','last','n_months','band','mean_coherence','valid_cells',
 'point_sig_cells','point_sig_pct','event_cells','event_months','epulead_pct','inflead_pct','inphase_pct','phase_concentration','band_mc_p','band_mc_p_Holm')]
for(o in main) {
 r=o$model;j=o$job;m=coi_mask(r);total=sum(m)
 counts=c(WTC=sum(m&r$Coherence.pval<=.05&r$Coherence>=.7),
  CWT_EPU=sum(m&r$Power.x.pval<=.05),CWT_INF=sum(m&r$Power.y.pval<=.05))
 for(k in names(counts))tableVI[[length(tableVI)+1]]=data.frame(series=j$series,method=k,
  valid_COI_cells=total,significant_cells=unname(counts[k]),percentage=100*counts[k]/total)
}
write.csv(do.call(rbind,tableVI),'results/Table_VI_power_coherence.csv',row.names=FALSE)
write.csv(tableVII,'results/Table_VII_period_direction.csv',row.names=FALSE)
# Side-by-side country summary heatmaps, with explicit Monte Carlo multiplicity correction.
for(ext in c('pdf','png')) {
 if(ext=='pdf')pdf('figures/band_summary.pdf',width=12,height=9) else png('figures/band_summary.png',width=1800,height=1350,res=150)
 layout(matrix(1:2,nrow=1));par(mar=c(4,7,4,1))
 ss=sort(unique(sum$series));ss=rev(ss)
 for(metric in c('point_sig_pct','mean_coherence')) {
  z=matrix(NA_real_,18,3)
  for(i in 1:18)for(k in 1:3)z[i,k]=sum[sum$series==ss[i]&sum$band==band_levels[k],metric]
  image(1:3,1:18,t(z),col=pal,zlim=if(metric=='mean_coherence')c(0,1) else c(0,max(z)),axes=FALSE,xlab='',ylab='')
  axis(1,at=1:3,labels=c('2-<8','8-<32','32+'),cex.axis=.85);axis(2,at=1:18,labels=ss,las=1,cex.axis=.85)
  for(i in 1:18)for(k in 1:3)text(k,i,if(metric=='mean_coherence')sprintf('%.2f',z[i,k]) else sprintf('%.1f%%',z[i,k]),cex=.8,
   col=if((z[i,k]-(if(metric=='mean_coherence')0 else 0))/(if(metric=='mean_coherence')1 else max(z))>.55)'black' else 'white')
  title(main=if(metric=='mean_coherence')'Mean coherence inside package COI' else 'Pointwise significant, coherence >= 0.70')
  mtext('Period band (months)',side=1,line=2.6);box()
 }
 dev.off()
}
# Seasonal-adjustment diagnostic figure.
sd=read.csv('results/seasonality_diagnostics.csv');sc=read.csv('results/seasonality_scores.csv')
png('figures/seasonality_comparison.png',width=1500,height=1000,res=150)
par(mfrow=c(1,2),mar=c(6,4.5,3,1))
cols=c(classical='#AA5A46',stl_level='#7496AC',stl_log='#224E69')
methods=c('classical','stl_level','stl_log')
barplot(sc$rejecting_countries[match(methods,sc$method)],names.arg=c('Classical','STL level','STL log'),col=cols[methods],
 ylim=c(0,18),ylab='Countries with residual month effect (5%)',main='Residual seasonality: 18 countries')
boxplot(month_partial_r2~factor(method,levels=methods),data=sd[sd$method!='nsa',],names=c('Classical','STL level','STL log'),
 col=cols[methods],ylab='Partial R² of calendar-month indicators',main='Magnitude after flexible trend control')
dev.off()
# China common-window proxy comparison.
cc=all[all$role=='proxy_matched',]
png('figures/china_proxy_comparison.png',width=1500,height=950,res=150)
par(mar=c(5,5,3,2));mat=matrix(NA_real_,3,3,dimnames=list(c('China','SCMP China','Mainland China'),band_levels))
for(s in rownames(mat))for(b in colnames(mat))mat[s,b]=cc$point_sig_pct[cc$series==s&cc$band==b]
bp=barplot(mat,beside=TRUE,col=c('#15334A','#6289A2','#B56A51'),ylab='Significant, high-coherence cells (%)',
 xlab='Period band (months)',main=paste('China proxy sensitivity | matched window',cc$first[1],'to',cc$last[1]),ylim=c(0,max(mat)*1.25+.1))
legend('topright',legend=rownames(mat),fill=c('#15334A','#6289A2','#B56A51'),bty='n');dev.off()
# Save structured report inputs for the workbook and reproducibility documentation.
tables=list(main=sum,all_bands=all,tableVI=do.call(rbind,tableVI),tableVII=tableVII,
 coverage=read.csv('results/coverage.csv'),seasonality=sd,SA_scores=sc,
 revisions=read.csv('results/endpoint_revisions.csv'),descriptives=read.csv('results/descriptive_statistics.csv'),
 stationarity=read.csv('results/stationarity_tests.csv'),attempts=read.csv('results/full_sample_attempts.csv'),
 conservative=read.csv('results/conservative_coi_summary.csv'),thresholds=read.csv('results/threshold_sensitivity.csv'),
 biwavelet=read.csv('results/biwavelet_summary.csv'),parameters=read.csv('results/parameter_sensitivity.csv'),
 china_matched=cc,processed=read.csv('results/processed_panel.csv'),jobs=read.csv('results/run_status.csv'),
 robustness=read.csv('results/robustness_comparison.csv'),drift=read.csv('results/stationarity_robustness.csv'))
jsonlite::write_json(tables,'results/report_tables.json',dataframe='rows',auto_unbox=TRUE,na='null',digits=15)
cat('Figures and tables complete. Holm significant main bands:',sum(sum$band_mc_p_Holm<=.05),'of',nrow(sum),'\n')
