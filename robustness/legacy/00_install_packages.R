# Install standard CRAN libraries. Use an R release supported by CRAN.
packages=c('WaveletComp','biwavelet','urca','jsonlite','viridisLite')
missing=packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
for(p in packages)cat(p,as.character(packageVersion(p)),'\n')
