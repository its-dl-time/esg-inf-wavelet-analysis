#!/usr/bin/env Rscript
# Run from the repository root. Inference is over prespecified country x band maps.
options(stringsAsFactors=FALSE, warn=1)
if (!requireNamespace('WaveletComp', quietly=TRUE)) stop('Install WaveletComp first (see README).')
source('R/wtc_functions.R')

env_int <- function(name, default) {
 z <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
 if (is.na(z) || z < 1) stop(name, ' must be a positive integer')
 z
}
B_CAL <- env_int('WTC_CAL', 250L)
B_TEST <- env_int('WTC_TEST', 999L)
NULL_METHOD <- Sys.getenv('WTC_NULL', 'tvar1')
if (!NULL_METHOD %in% c('tvar1','ar1','fourier')) stop('Invalid WTC_NULL')
WIN_T <- env_int('WTC_TIME', 9L)
WIN_S <- as.numeric(Sys.getenv('WTC_SCALE', '1'))
if (!is.finite(WIN_S) || WIN_S <= 0) stop('Invalid WTC_SCALE')
OUT <- Sys.getenv('WTC_OUT', file.path('results','core'))
dir.create(file.path(OUT,'checkpoints'), recursive=TRUE, showWarnings=FALSE)
dat <- read.csv('data/processed/core_18.csv', na.strings=c('','NA'))
countries <- sort(unique(dat$series))
stopifnot(length(countries)==18L)
requested <- Sys.getenv('WTC_COUNTRIES', '')
if (nzchar(requested)) {
 countries <- trimws(strsplit(requested, ',', fixed=TRUE)[[1]])
 if (any(!countries %in% unique(dat$series)) || anyDuplicated(countries)) stop('Invalid WTC_COUNTRIES')
}
band_names <- c('2-<8','8-<32','32+')
cfg <- list(B_CAL=B_CAL,B_TEST=B_TEST,null=NULL_METHOD,time=WIN_T,scale=WIN_S,
  sa='stl_log',dj=1/20,lowerPeriod=2,band_names=band_names,threshold=.70,
  episode_min_months=6L,seed=24414L)
saveRDS(cfg, file.path(OUT,'configuration.rds'))

run_one <- function(s, index) {
 d <- dat[dat$series==s,]
 d <- d[order(d$month),]
 ok <- is.finite(d$epu_mom) & is.finite(d$inf_stl_log)
 rr <- rle(ok); ends <- cumsum(rr$lengths); k <- which(rr$values)[which.max(rr$lengths[rr$values])]
 if (!length(k)) stop('No complete window for ',s)
 ids <- seq.int(ends[k]-rr$lengths[k]+1L, ends[k]); d <- d[ids,]
 if (any(diff(as.integer(substr(d$month,1,4))*12L+as.integer(substr(d$month,6,7)))!=1L)) stop('Non-contiguous months: ',s)
 x <- d$epu_mom; y <- d$inf_stl_log
 upper <- min(128, floor(length(x)/3))
 params <- list(t=WIN_T, s=WIN_S, dj=1/20, upper=upper)
 set.seed(cfg$seed + index*10007L)
 obs <- wt_map(x,y,params)
 mask <- wt_coi(obs)
 band <- cut(obs$Period,c(-Inf,8,32,Inf),right=FALSE,labels=band_names)
 valid_band <- vapply(band_names,function(b)any(mask[band==b,,drop=FALSE]),logical(1))
 if (any(!valid_band)) stop('Empty COI band: ',s)
 source_x <- make_surrogate(x,NULL_METHOD); source_y <- make_surrogate(y,NULL_METHOD)
 # Independent calibration and testing Monte Carlo batches. No observed map selects thresholds.
 cal <- array(NA_real_, c(obs$nr,obs$nc,B_CAL))
 for (j in seq_len(B_CAL)) {
  z <- wt_pair(source_x(), source_y(), params)
  cal[,,j] <- z$Coherence
 }
 threshold <- matrix(apply(cal,c(1,2),stats::quantile,probs=.95,na.rm=TRUE),obs$nr,obs$nc)
 threshold <- pmax(.70,threshold)
 obs_stat <- wt_statistics(obs$Coherence,mask,band,threshold,6L)
 null_stat <- array(NA_real_, c(B_TEST,3L,3L),dimnames=list(NULL,band_names,c('band_mean','cluster_mass','episode_mean')))
 for (j in seq_len(B_TEST)) {
  z <- wt_pair(source_x(), source_y(), params)
  null_stat[j,,] <- wt_statistics(z$Coherence,mask,band,threshold,6L)
 }
 p <- matrix(NA_real_,3,3,dimnames=dimnames(obs_stat))
 for (a in seq_len(3L)) for (b in seq_len(3L))
  p[a,b] <- (1+sum(null_stat[,a,b]>=obs_stat[a,b]))/(B_TEST+1)
 rows <- expand.grid(band=band_names,statistic=colnames(obs_stat),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
 rows$series <- s; rows$first <- d$month[1]; rows$last <- tail(d$month,1)
 rows$n_months <- nrow(d); rows$observed <- as.vector(obs_stat); rows$p <- as.vector(p)
 rows$null_method <- NULL_METHOD; rows$window_t <- WIN_T; rows$window_s <- WIN_S
 rows$valid_cells <- vapply(band_names,function(b)sum(mask[band==b,,drop=FALSE]),numeric(1))[match(rows$band,band_names)]
 list(summary=rows, month=d$month, period=obs$Period, coherence=obs$Coherence,
  phase=obs$sAngle, coi=mask, threshold=threshold, cfg=cfg)
}

out <- list()
for (i in seq_along(countries)) {
 s <- countries[i]; path <- file.path(OUT,'checkpoints',paste0(s,'.rds'))
 old <- if (file.exists(path)) tryCatch(readRDS(path),error=function(e)NULL) else NULL
 # Input hash is essential: never silently reuse a checkpoint after data change.
 input <- dat[dat$series==s,c('month','epu_mom','inf_stl_log')]
 signature <- serialize(list(cfg=cfg,input=input),NULL,version=3)
 if (!is.null(old) && identical(old$signature,signature)) {
  result <- old$result
 } else {
  cat('Running',s,'(',i,'/',length(countries),')\n');flush.console()
  result <- run_one(s, match(s,sort(unique(dat$series))))
  temp <- paste0(path,'.tmp_',Sys.getpid())
  saveRDS(list(signature=signature,result=result),temp,compress='xz')
  if (!file.rename(temp,path)) stop('Cannot save checkpoint: ',path)
 }
 out[[s]] <- result$summary
}
summary <- do.call(rbind,out); rownames(summary) <- NULL
summary$q_BH <- NA_real_; summary$q_BY <- NA_real_; summary$p_Holm <- NA_real_
for (metric in unique(summary$statistic)) {
 k <- which(summary$statistic==metric)
 # 18 x 3 = 54 tests for a full run. Subset runs are marked exploratory in README.
 summary$q_BH[k] <- p.adjust(summary$p[k],method='BH')
 summary$q_BY[k] <- p.adjust(summary$p[k],method='BY')
 summary$p_Holm[k] <- p.adjust(summary$p[k],method='holm')
}
write.csv(summary,file.path(OUT,'summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
cat('Completed',length(countries),'countries. Output:',OUT,'\n')
