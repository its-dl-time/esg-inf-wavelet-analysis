#!/usr/bin/env Rscript
# Descriptive conditional WTC, not calibrated conditional significance.
options(stringsAsFactors=FALSE)
if (!requireNamespace('biwavelet',quietly=TRUE)) stop('Install biwavelet')
root <- 'extensions/uct'
u <- read.csv(file.path(root,'data/UCT.csv'),check.names=FALSE)
o <- read.csv(file.path(root,'data/MCOILBRENTEU.csv'),check.names=FALSE)
panel <- read.csv('data/processed/core_18.csv')

# Inspect the upstream header. The explicit mapping is intentional: silent selection
# of an arbitrary numeric column would change the control variable.
uct_date <- Sys.getenv('UCT_DATE_COLUMN','')
uct_value <- Sys.getenv('UCT_VALUE_COLUMN','')
if (!nzchar(uct_date) || !nzchar(uct_value))
 stop('Set UCT_DATE_COLUMN and UCT_VALUE_COLUMN after inspecting UCT.csv header. Columns: ',paste(names(u),collapse=', '))
if (!all(c(uct_date,uct_value)%in%names(u))) stop('UCT column mapping not found')
uct <- data.frame(month=substr(as.character(u[[uct_date]]),1,7),uct=as.numeric(u[[uct_value]]))
if (!all(c('DATE','MCOILBRENTEU')%in%names(o))) stop('Unexpected FRED header')
oil <- data.frame(month=substr(o$DATE,1,7),brent=as.numeric(o$MCOILBRENTEU))
if (anyDuplicated(uct$month) || anyDuplicated(oil$month)) stop('Duplicate control months')
transform_uct <- Sys.getenv('UCT_TRANSFORM','diff')
if (!transform_uct%in%c('diff','level')) stop('UCT_TRANSFORM must be diff or level')
uct <- uct[order(uct$month),];oil <- oil[order(oil$month),]
if(transform_uct=='diff') uct$uct <- c(NA_real_,diff(uct$uct))
oil$brent <- c(NA_real_,100*diff(log(oil$brent)))

ff <- 4*pi/(6+sqrt(38))
pair <- function(a,b) {
 n <- length(a); biwavelet::wtc(cbind(seq_len(n),as.numeric(scale(a))),
  cbind(seq_len(n),as.numeric(scale(b))),dj=1/20,s0=2/ff,
  max.scale=min(128,floor(n/3))/ff,mother='morlet',nrands=0,quiet=TRUE)
}
rho <- function(w) sqrt(pmax(0,pmin(1,w$rsq)))*exp(1i*w$phase)
partial <- function(ab,ac,bc) {
 den <- sqrt((1-Mod(ac)^2)*(1-Mod(bc)^2))
 z <- (ab-ac*Conj(bc))/den
 z[!is.finite(den)|den<.05] <- NA_complex_
 Mod(z)^2
}
multiple <- function(R,variables) {
 dims <- dim(R[[1]]);out <- matrix(NA_real_,dims)
 for(k in seq_len(prod(dims))) {
  C <- diag(1+0i,length(variables)+1L)
  labs <- c('INF',variables)
  for(i in seq_along(labs)) for(j in seq_along(labs)) if(i!=j) {
   key <- paste(sort(c(labs[i],labs[j])),collapse=':')
   z <- R[[key]][k]
   C[i,j] <- if(labs[i]<labs[j]) z else Conj(z)
  }
  Q <- C[-1,-1,drop=FALSE];v <- C[1,-1,drop=FALSE]
  if(any(!is.finite(C)) || min(Re(eigen(Q,only.values=TRUE)$values))<.05) next
  a <- Re(v%*%solve(Q,Conj(t(v))))
  if(is.finite(a) && a>=-1e-5 && a<=1+1e-5)out[k] <- max(0,min(1,a))
 }
 out
}
dir.create(file.path(root,'results'),showWarnings=FALSE)
for(country in c('US','China')) {
 d <- panel[panel$series==country,c('month','epu_mom','inf_stl_log')]
 d <- merge(merge(d,uct,by='month'),oil,by='month')
 d <- d[order(d$month),];d <- d[complete.cases(d),]
 if (any(diff(as.integer(substr(d$month,1,4))*12L+as.integer(substr(d$month,6,7)))!=1L))
  stop('Matched data have calendar gaps: ',country)
 stopifnot(nrow(d)>100)
 v <- list(INF=d$inf_stl_log,EPU=d$epu_mom,UCT=d$uct,Oil=d$brent)
 var_names <- names(v); R <- list();W <- list()
 for(i in seq_len(length(v)-1L)) for(j in (i+1L):length(v)) {
  a <- var_names[i];b <- var_names[j]; key <- paste(sort(c(a,b)),collapse=':')
  w <- pair(v[[a]],v[[b]]);z <- rho(w)
  # Pairwise phase is first series minus second; orient to lexical key.
  R[[key]] <- if(a<b) z else Conj(z); W[[key]] <- w
 }
 ref <- W[['EPU:INF']];period <- ref$period
 coi <- Reduce('&',lapply(W,function(w)outer(w$period,as.vector(w$coi),'<=') ))
 epu_inf <- Mod(R[['EPU:INF']])^2
 p_uct <- partial(R[['EPU:INF']],R[['EPU:UCT']],R[['INF:UCT']])
 p_oil <- partial(R[['EPU:INF']],R[['EPU:Oil']],R[['INF:Oil']])
 models <- list(unconditional=epu_inf,partial_UCT=p_uct,partial_oil=p_oil,
  multiple_EPU_UCT=multiple(R,c('EPU','UCT')),
  multiple_EPU_oil=multiple(R,c('EPU','Oil')),
  multiple_EPU_UCT_oil=multiple(R,c('EPU','UCT','Oil')))
 rows <- list()
 for (label in names(models)) for (band in c('2-<8','8-<32','32+')) {
  ix <- if(band=='2-<8') period<8 else if(band=='8-<32')period>=8&period<32 else period>=32
  z <- models[[label]][ix,,drop=FALSE];valid <- coi[ix,,drop=FALSE]&is.finite(z)
  rows[[length(rows)+1L]] <- data.frame(country=country,model=label,band=band,
   first=d$month[1],last=tail(d$month,1),n=nrow(d),valid_cells=sum(valid),
   mean_coherence=if(any(valid))mean(z[valid]) else NA_real_,
   uct_transform=transform_uct)
 }
 write.csv(do.call(rbind,rows),file.path(root,'results',paste0(country,'_',transform_uct,'.csv')),row.names=FALSE)
}
