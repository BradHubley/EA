## Spatial Group
##
##  Purpose of the spatial working group: To incorporate environmental/ecosystem data into a single species stock assessment using a spatial approach. 
##  Most environmental data such as depth, temperature, substrate and species distribution is spatial in nature. Survey data for single species stock 
##  assessments are also collected across space and time but space is rarely considered in most assessments. Levaging spatial environmental data would 
##  allow us to better estimate species distribution and abundance patterns, understand how fishing patterns affect different components of the stock 
##  and produce better science advice. We already do this informally in the CSAS process but in this workshop we can begin to explore how to do it 
##  explicitly by incorporating avaialable enviromental and ecosystem data sources.
##
## Spatial Approaches						Feed into Management Advice (What is spatial Management?)
##
## 	- Species distribution models			->	quantify available habitat
##
##  - Habitat suitability models			->	quantify available habitat, explicitly spatial
##
##  - Spatial Abundance models				->	quantify abundance and distribution, construct an abundance index that can feed into a fisheries model
##
##  - Spatial-temporal population dynamics	->	explicit incorporation of spatial abundance and distribution patterns into a fisheries model
##
##
##
##

# set wd
setwd("/home/hubleyb/bio/EA/")

# basic R packages
library(dismo)
library(mgcv)
library(fields)
library(RandomFields)
library(Matrix)
library(lattice)
library(spatstat)
library(lubridate)
library(devtools)

# gitHub R packages
library(TMB) # https://github.com/kaskr/adcomp, don't worry about this one right now
library(bio.utilities) # https://github.com/Beothuk/bio.utilities
#to install: install_github("Beothuk/bio.utilities")


## Abundance Data (Survey)

lobsterRV = read.csv("data/lobsterRV.csv")
scallopRV = read.csv("data/scallopRV.csv")
silverhakeRV = read.csv("data/silverhakeRV.csv")
halibutRV = read.csv("data/halibutRV.csv")

## Environmental Data

load("data/baseline.rdata") # baseLine
load("data/isobaths.rdata") # isoBaths
load("data/coastline.rdata") # coastLine
load("data/predspace.rdata") # predSpace
load("data/predspacetime.rdata") # predSpaceTime


load_all() # loads functions


# interpolate abundance for visuals
year=2014

	# Lobster
	xyz = subset(lobsterRV,year(sdate)==year,c('plon','plat','totwgt'))


		datarange = c(min(xyz$totwgt[xyz$totwgt>0])*0.5,quantile(xyz$totwgt,0.99))
		datascale = seq(datarange[1],datarange[2],l=50)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
	# add zeros
    #xyz =na.omit( zeroInflate(xyz,corners=corners,type=2,type.scaler=0.5,eff=datarange[1]) )

	EAmap( xyz,  pts=xyz[,c("plon","plat")], fn="LobsterRV", annot=year, loc="output", datascale=datascale , corners=corners, interpolation=T,log.variable=T,display=T,add.zeros=T)

	# Scallop
	xyz = subset(scallopRV,year(sdate)==year,c('plon','plat','totwgt'))

		datarange = c(min(xyz$totwgt[xyz$totwgt>0])*0.5,quantile(xyz$totwgt,0.99))
		datascale = seq(datarange[1],datarange[2],l=50)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
	# add zeros
    #xyz =na.omit( zeroInflate(xyz,corners=corners,type=2,type.scaler=0.5,eff=datarange[1]) )

	EAmap( xyz,  pts=xyz[,c("plon","plat")], fn="ScallopRV", annot=year, loc="output", datascale=datascale , corners=corners, interpolation=T,log.variable=T,display=T,add.zeros=T)

	# silverhake
	xyz = subset(silverhakeRV,year(sdate)==year,c('plon','plat','totwgt'))

		datarange = c(min(xyz$totwgt[xyz$totwgt>0])*0.5,quantile(xyz$totwgt,0.99))
		datascale = seq(datarange[1],datarange[2],l=50)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
	# add zeros
    #xyz =na.omit( zeroInflate(xyz,corners=corners,type=2,type.scaler=0.5,eff=datarange[1]) )

	EAmap( xyz,  pts=xyz[,c("plon","plat")], fn="silverhakeRV", annot=year, loc="output", datascale=datascale , corners=corners, interpolation=T,log.variable=T,display=T,add.zeros=T)

	# Halibut
	xyz = subset(halibutRV,year(sdate)==year,c('plon','plat','totwgt'))

		datarange = c(min(xyz$totwgt[xyz$totwgt>0])*0.5,quantile(xyz$totwgt,0.99))
		datascale = seq(datarange[1],datarange[2],l=50)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
	# add zeros
    #xyz =na.omit( zeroInflate(xyz,corners=corners,type=2,type.scaler=0.5,eff=datarange[1]) )

	EAmap( xyz,  pts=xyz[,c("plon","plat")], fn="HalibutRV", annot=year, loc="output", datascale=datascale , corners=corners, interpolation=T,log.variable=T,display=T,add.zeros=T)


## Environmental Data

# Bathymetry 

	# Depth
	xyz = predSpace[c('plon','plat','z')]
	datascale = seq(10,1000,l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="Depth", loc="output", datascale=datascale , corners=corners, log.variable=T, display=T,rev=T)


	# Slope
	xyz = predSpace[c('plon','plat','dZ')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99))
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="Slope", loc="output", datascale=datascale , corners=corners, display=T,log.variable=T)

	# Curvature
	xyz = predSpace[c('plon','plat','ddZ')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99))
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="Curvature", loc="output", datascale=datascale , corners=corners, display=T,log.variable=T)

  

# Substrate

	# Grain Size
	xyz = predSpace[c('plon','plat','log.substrate.grainsize')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="Grain.Size", loc="output", datascale=datascale , corners=corners, display=T)


# Temperature

	# Climatology mean
	xyz = predSpace[c('plon','plat','tmean.climatology')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="tmean.climatology", loc="output", datascale=datascale , corners=corners, display=T)

	# Climatology sd
	xyz = predSpace[c('plon','plat','tsd.climatology')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="tsd.climatology", loc="output", datascale=datascale , corners=corners, display=T)

	# Climatology min
	xyz = predSpace[c('plon','plat','tmin.climatology')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="tmin.climatology", loc="output", datascale=datascale , corners=corners, display=T)

	# Climatology max
	xyz = predSpace[c('plon','plat','tmax.climatology')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="tmax.climatology", loc="output", datascale=datascale , corners=corners, display=T)

	# Climatology amplitude
	xyz = predSpace[c('plon','plat','amplitude.climatology')]

	datarange = quantile(xyz[,3],probs=c(0.01,0.99),na.rm=T)
	datascale = seq(datarange[1],datarange[2],l=50)
	corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

	EAmap( xyz, fn="amplitude.climatology", loc="output", datascale=datascale , corners=corners, display=T)


## 	Species distribution / Habitat suitability models

#dismo

modelformula = formula( paste( 
      varname, ' ~ s(t, k=3, bs="ts") + s(tmean.climatology, k=3, bs="ts") + s(tsd.climatology, k=3, bs="ts")  ', 
      ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
      ' + s( log(mr), k=3, bs="ts") + s( Npred, k=3, bs="ts") + s( smr, k=3, bs="ts")  ',
      ' + s(log.substrate.grainsize, k=3, bs="ts") + s(ca1, k=3, bs="ts") + s(ca2, k=3, bs="ts")   ' ))  # no space 



##  with space




##  Spatial Abundance models




##  Spatial-temporal population dynamics