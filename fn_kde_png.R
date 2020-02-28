#JTVD - 2020

#-------------------------------------------------------------------------------
# libraries and set-up
#-------------------------------------------------------------------------------

#libraries
suppressMessages(library(data.table))
suppressMessages(library(sparr))
suppressMessages(library(raster))
suppressMessages(library(leaflet))
suppressMessages(library(RColorBrewer))

#parameters
args <- commandArgs(trailingOnly=TRUE)
year <- as.numeric(args[1])

#-------------------------------------------------------------------------------
# read from stdin
#-------------------------------------------------------------------------------

#surname data
pop.sur  <- na.omit(fread('file:///dev/stdin',col.names=c('surname','id','x','y')))
agg.sur  <- pop.sur[,.(n=.N),by=.(id,x,y)]
name.sur <- pop.sur$surname[1]
num.sur  <- nrow(pop.sur)

#-------------------------------------------------------------------------------
# read from disk
#-------------------------------------------------------------------------------

#load xy grid (xy), gb outline (gb), population density surface (pop.rdf)
load('rdata/xy.RData')
load('rdata/gb.RData')
load(paste0('rdata/pop',year,'.RData'))

#-------------------------------------------------------------------------------
# calculation individual KDE's
#-------------------------------------------------------------------------------

#prepare surname
gb.win  <- as.mask(owin(c(-19693.0,676307.0),c(-15372.0,1238628.0),unitname='meters'),eps=1000)
ppp.sur <- ppp(x=agg.sur$x,y=agg.sur$y,window=gb.win,checkdup=FALSE)

#bandwidth calculation
if (num.sur > 5000) {
  smp.sur <- pop.sur[sample(.N,size=5000,replace=FALSE)]
  ppp.smp <- ppp(x=smp.sur$x,y=smp.sur$y,window=gb.win,checkdup=FALSE)
  bwd.sur <- round(LIK.density(ppp.smp,hlim=c(10000,60000),verbose=FALSE))
} else {
  bwd.sur <- round(LIK.density(ppp.sur,hlim=c(10000,60000),verbose=FALSE))
}

#kde
sur.dns <- bivariate.density(ppp.sur,h0=bwd.sur,weights=agg.sur$n,intensity=FALSE,xy=list(x=xy$x,y=xy$y),adapt=FALSE)

#weigh surname by population
sur.rdf <- as.data.table(sur.dns$z)
sur.wgh <- cbind(sur.rdf,pop_value=pop.rdf[,value])
sur.wgh[,w := ifelse(value!=0 & pop_value!=0,value/-log(pop_value),
              ifelse(value!=0 & pop_value==0,value/min(-log(pop.rdf$value)),
                     value))]

#to raster
wgh.rdf <- sur.wgh[,.(x,y,w)]
wgh.rst <- rasterFromXYZ(wgh.rdf,crs=CRS('+init=epsg:27700'))

#clip to gb coastline, project
clp.rst <- mask(wgh.rst,gb)
prj.rst <- projectRaster(clp.rst,crs=CRS('+init=epsg:3857'))

#map values onto colours
colors <- leaflet::colorNumeric(brewer.pal(n=11,name='Spectral')[c(1:6,10,11)],-1:1001,na.color='transparent',reverse=TRUE)
cols   <- c(colors(-1:1001),colors(NA))
cols[c(1:155)] <- 'transparent'

#scale values
cval <- scale(values(prj.rst))   
minmax <- quantile(cval,probs=c(0.01,0.99),na.rm=TRUE)

#map values to colours
vals <- round((((cval-minmax[1])/(minmax[2]-minmax[1]))*1000))
vals[vals < 0] <- 0
vals[vals > 1000] <- 1000
vals[is.na(vals)] <- 1002

#lookup colors for scaled values, convert to raw, and write to file
valcolors <- cols[vals+2]
rgb_data <- col2rgb(valcolors,alpha=TRUE)
raw_data <- as.raw(rgb_data)
dim(raw_data) <- c(4,ncol(prj.rst),nrow(prj.rst))

#output parameters
folder <- paste0('kde/',substr(name.sur,1,1))
subfolder <- paste0('kde/',substr(name.sur,1,1),'/',name.sur)
kde_png <- paste0(subfolder,'/kde',year,'.png')

#output
dir.create(folder,showWarnings=FALSE)
dir.create(subfolder,showWarnings=FALSE)
png::writePNG(raw_data,kde_png)

#meta data to stdout
cat(paste0(name.sur,',',year,',',num.sur,',',bwd.sur,'\n'))