#JTVD - 2020

#-------------------------------------------------------------------------------
# libraries and set-up
#-------------------------------------------------------------------------------

#libraries
suppressMessages(library(data.table))
suppressMessages(library(sparr))
suppressMessages(library(raster))

#parameters
options(warn=-1)
year <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

#-------------------------------------------------------------------------------
# read from stdin
#-------------------------------------------------------------------------------

#surname data
pop.sur  <- na.omit(fread('file:///dev/stdin',col.names=c('surname','id','x','y')))
agg.sur <- pop.sur[,.(n=.N),by=.(id,x,y)]
name.sur <- pop.sur$surname[1]

#-------------------------------------------------------------------------------
# read from disk
#-------------------------------------------------------------------------------

#load xy grid (xy), isonymy classes (iso), population density surface (pop.rdf)
load('data/xy.RData')
load('data/iso.RData')
load(paste0('data/ipop',year,'.RData'))

#-------------------------------------------------------------------------------
# calculation surname KDE
#-------------------------------------------------------------------------------

#prepare surname
gb.win  <- as.mask(owin(c(-19693.0,676307.0),c(-15372.0,1238628.0),unitname='meters'),eps=1000)
ppp.sur <- ppp(x=agg.sur$x,y=agg.sur$y,window=gb.win,checkdup=FALSE)
num.sur <- nrow(pop.sur)

#bandwidth selection
iso.sur <- iso[surname==name.sur]
if (nrow(iso.sur) > 0) {
  #widespread, large
  if(iso.sur$iso <= 9 & num.sur > 10000) {
    bwd.sur <- 18000
  #generally small
  } else if(num.sur < 500) {
    bwd.sur <- 8000
  #concentrated, medium
  } else {bwd.sur <- 12000}
#other
} else {bwd.sur <- 10000}

#kde
sur.dns <- bivariate.density(ppp.sur,h0=bwd.sur,weights=agg.sur$n,intensity=TRUE,xy=list(x=xy$x,y=xy$y),adapt=FALSE)

#-------------------------------------------------------------------------------
# weigh surname by total population
#-------------------------------------------------------------------------------

#weigh surname by population
sur.rdf <- as.data.table(sur.dns$z)
sur.wgh <- cbind(sur.rdf,pop_value=pop.rdf[,value])
sur.wgh[,w := value]
sur.wgh[,w := ifelse(value!=0 & pop_value!=0,value/-log(pop_value),
              ifelse(value!=0 & pop_value==0,value/min(-log(pop.rdf$value)),value))]

#-------------------------------------------------------------------------------
# normalise values 0-1 with output normalisation factor
#-------------------------------------------------------------------------------

#output normalisation factor using median value (40,000 sample)
if (num.sur <= 1000) {
  #max sample value freq <= 1,000, population correction
  corfac <- 4.2e-08 * (num.sur/1000)
} else if(num.sur > 1000 & num.sur < 10000) {
  #max sample value freq <= 10,000, population correction
  corfac <- 2.2e-07 * (num.sur/10000)
} else if (num.sur >= 10000 & num.sur < 40000) {
  #max sample value freq <= 40,000, population correction
  corfac <- 4.7e-07 * (num.sur/40000)
} else if (num.sur >= 40000 & num.sur < 100000) {
  #max sample value freq <= 100,000, population correction
  corfac <- 1.1e-06 * (num.sur/100000)
} else {
  #max sample value freq <= 100,000, no population correction
  corfac <- 1.1e-06
}

#normalise
sur.nor <- sur.wgh[,.(x,y,w)]
sur.nor[,n := w/corfac]

#-------------------------------------------------------------------------------
# reclassify 
#-------------------------------------------------------------------------------

#kde contours
sur.nor[,r := .I]
sur.nor[,c := ifelse(n <= 0.10,0,
              ifelse(n >  0.10 & n <= 0.20,1,
              ifelse(n >  0.20 & n <= 0.30,2,
              ifelse(n >  0.30,3,0))))]

#-------------------------------------------------------------------------------
# prepare output, write to stdout
#-------------------------------------------------------------------------------

#filter
sur.kde <- sur.nor[c != 0]

#collapse data
out.idx <- paste0('',sur.kde$r,collapse=',')
out.kde <- paste0('',sur.kde$c,collapse=',')

#stdout
cat(paste0(name.sur,';',year,';',num.sur,';',bwd.sur,';',out.idx,';',out.kde))