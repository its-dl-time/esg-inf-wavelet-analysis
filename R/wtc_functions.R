# Morlet WTC and matched-surrogate statistics. No package-internal objects are modified.
wt_map <- function(x,y,p) {
 WaveletComp::analyze.coherency(data.frame(EPU=x,INF=y),my.pair=c('EPU','INF'),
  loess.span=0,dt=1,dj=p$dj,lowerPeriod=2,upperPeriod=p$upper,
  window.type.t=1,window.type.s=1,window.size.t=p$t,window.size.s=p$s,
  make.pval=FALSE,verbose=FALSE)
}
wt_pair <- function(x,y,p) {
 WaveletComp::WaveletCoherency(x,y,dt=1,dj=p$dj,lowerPeriod=2,upperPeriod=p$upper,
  window.type.t=1,window.type.s=1,window.size.t=p$t,window.size.s=p$s)
}
wt_coi <- function(r) outer(r$Period,2^r$coi.2[seq_len(r$nc)+2L],'<=' )

make_surrogate <- function(x,method) {
 if (method=='ar1') return(function() WaveletComp::SurrogateData(x,method='AR',params=list(AR=list(p=1))))
 if (method=='fourier') return(function() WaveletComp::SurrogateData(x,method='Fourier.rand'))
 # Heteroskedastic AR(1): fixed phi and centered, month-specific innovation variance.
 # This is a local-volatility null, not an exact wavelet-phase nonstationary surrogate.
 mu <- mean(x); v <- x-mu
 phi <- sum(v[-length(v)]*v[-1])/sum(v[-length(v)]^2)
 phi <- max(-.98,min(.98,phi))
 eps <- v[-1]-phi*v[-length(v)]
 w <- 25L; local_var <- vapply(seq_along(eps),function(i) {
  ids <- max(1L,i-w%/%2L):min(length(eps),i+w%/%2L)
  mean(eps[ids]^2)
 },numeric(1))
 sigma <- sqrt(pmax(local_var,1e-12)); n <- length(x)
 function() {
  z <- numeric(n); z[1] <- rnorm(1,sd=stats::sd(x))
  for (i in 2:n) z[i] <- phi*z[i-1]+sigma[i-1]*rnorm(1)
  z+mu
 }
}

cluster_max <- function(sig,excess) {
 nr <- nrow(sig); nc <- ncol(sig); seen <- matrix(FALSE,nr,nc); best <- 0
 # Eight-connected components in the period-time lattice, restricted to one band.
 for (start in which(sig)) {
  if (seen[start]) next
  queue <- integer(length(sig)); queue[1] <- start; head <- 1L; tail <- 1L
  seen[start] <- TRUE; mass <- 0
  while (head<=tail) {
   j <- queue[head];head <- head+1L; mass <- mass+excess[j]
   row <- (j-1L)%%nr+1L; col <- (j-1L)%/%nr+1L
   for (a in max(1L,row-1L):min(nr,row+1L))
    for (b in max(1L,col-1L):min(nc,col+1L)) {
     k <- a+(b-1L)*nr
     if (sig[k] && !seen[k]) {seen[k] <- TRUE;tail <- tail+1L;queue[tail] <- k}
    }
  }
  best <- max(best,mass)
 }
 best
}

episode_max <- function(sig,co,minimum=6L) {
 best <- 0
 for (i in seq_len(nrow(sig))) {
  r <- rle(sig[i,]);ends <- cumsum(r$lengths)
  for (j in which(r$values & r$lengths>=minimum)) {
   ids <- seq.int(ends[j]-r$lengths[j]+1L,ends[j])
   best <- max(best,mean(co[i,ids]))
  }
 }
 best
}

wt_statistics <- function(co,mask,band,threshold,minimum=6L) {
 labels <- c('2-<8','8-<32','32+')
 out <- matrix(0,3,3,dimnames=list(labels,c('band_mean','cluster_mass','episode_mean')))
 for (i in seq_along(labels)) {
  ix <- which(band==labels[i]); m <- mask[ix,,drop=FALSE]
  c <- co[ix,,drop=FALSE]; th <- threshold[ix,,drop=FALSE]
  out[i,1] <- mean(c[m]); sig <- m & is.finite(c) & c>=th
  out[i,2] <- cluster_max(sig,pmax(0,c-th))
  out[i,3] <- episode_max(sig,c,minimum)
 }
 out
}
