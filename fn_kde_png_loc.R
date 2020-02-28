#JTVD - 2020

#-------------------------------------------------------------------------------
# libraries and set-up
#-------------------------------------------------------------------------------

#libraries
suppressMessages(library(data.table))
suppressMessages(library(sparr))
suppressMessages(library(raster))
suppressMessages(library(leaflet))

#parameters
year <- 2016

#grid extent
x.min <- -19693.0
x.max <- 676307.0
y.min <- -15372.0
y.max <-  1238628.0

#-------------------------------------------------------------------------------
# read from stdin
#-------------------------------------------------------------------------------

#surname data
pop.sur <- na.omit(fread('data_jones',col.names=c('surname','id','x','y')))
agg.sur <- pop.sur[,.(n=.N),by=.(id,x,y)]
name.sur <- pop.sur$surname[1]

#-------------------------------------------------------------------------------
# read from disk
#-------------------------------------------------------------------------------

#load xy grid (xy), gb outline (gb), isonymy classes (iso), population density surface (pop.rdf)
load('rdata/xy.RData')
load('rdata/gb.RData')
load('rdata/iso.RData')
load(paste0('rdata/ipop',year,'.RData'))

#-------------------------------------------------------------------------------
# calculation individual KDE's
#-------------------------------------------------------------------------------

#prepare surname
gb.win <- as.mask(owin(c(x.min,x.max),c(y.min,y.max),unitname='meters'),eps=1000)
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
#to raster
wgh.rdf <- sur.wgh[,.(x,y,w)]
wgh.rst <- rasterFromXYZ(wgh.rdf,crs=CRS('+init=epsg:27700'))

#-------------------------------------------------------------------------------
# normalise value
#-------------------------------------------------------------------------------

#clip to gb coastline
clp.rst <- mask(wgh.rst,gb)
prj.rst <- projectRaster(clp.rst,crs=CRS('+init=epsg:3857'))

#normalise values 0-1, with population adjustment
nval <- as.data.table(values(clp.rst)/(corfac))

#map colours to values
nval[,c := ifelse(is.na(V1),'transparent',
           ifelse(V1 > 0.10 & V1 <= 0.20,'#9ecae1',
           ifelse(V1 > 0.20 & V1 <= 0.30,'#6baed6',
           ifelse(V1 > 0.30,'#4292c6','transparent'))))]

#lookup colors for scaled values, convert to raw and write to file
rgb_data <- col2rgb(nval$c,alpha=TRUE)
raw_data <- as.raw(rgb_data)
dim(raw_data) <- c(4,ncol(prj.rst),nrow(prj.rst))

#output parameters
kde_png <- 'web/kde_name.png'

#output
png::writePNG(raw_data,kde_png,dpi=c(10))

#meta data to stdout
cat(paste0(name.sur,',',year,',',num.sur,',',bwd.sur,'\n'))
