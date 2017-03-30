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
library(sp)
library(rgdal)
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
library(bio.base) # https://github.com/Beothuk/bio.base
#to install: install_github("Beothuk/bio.base")
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

Years = unique(year(predSpaceTime$timestamp))

load_all() # loads functions


## interpolate abundance for visuals
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
	#Climate Envelop Model Booth 2014 Diversity and Distributions 20:1-9

	# Description

	# The Bioclim algorithm has been extensively used for species distribution modeling. Bioclim is the
	# classic ’climate-envelope-model’. Although it generally does not perform as good as some other
	# modeling methods (Elith et al. 2006) and is unsuited for predicting climate change effects (Hijmans
	# and Graham, 2006). It is still used, however, among other reasons because the algorithm is easy to
	# understand and thus useful in teaching species distribution modeling.

	# The BIOCLIM algorithm computes the similarity of a location by comparing the values of environmental
	# variables at any location to a percentile distribution of the values at known locations of
	# occurrence (’training sites’). The closer to the 50th percentile (the median), the more suitable the
	# location is. The tails of the distribution are not distinguished, that is, 10 percentile is treated as
	# equivalent to 90 percentile.

	# In this R implementation, percentile scores are between 0 and 1, but predicted values larger than 0.5
	# are subtracted from 1. Then, the minimum percentile score across all the environmental variables
	# is computed (i.e. this is like Liebig’s law of the minimum, except that high values can also be
	# limiting factors). The final value is subtracted from 1 and multiplied with 2 so that the results are
	# between 0 and 1. The reason for this transformation is that the results become more like that of
	# other distribution modeling methods and are thus easier to interpret. The value 1 will rarely be
	# observed as it would require a location that has the median value of the training data for all the
	# variables considered. The value 0 is very common as it is assigned to all cells with a value of an
	# environmental variable that is outside the percentile distribution (the range of the training data) for
	# at least one of the variables.
	
	# In the predict function, you can choose to ignore one of the tails of the distribution (e.g, to make
	# low rainfall a limiting factor, but not high rainfall),

		dat = lobsterRV
		dat = scallopRV
		dat = silverhakeRV
		dat = halibutRV
		dat$Y = ifelse(dat$totwgt>0,1,0)
		dat$dyear = decimal_date(dat$sdate)
		dat$t = dat$bottom_temperature

		dat.p = subset(dat,Y==1)

		# just temp, depth, slope, curvature
		bc = bioclim(dat.p[,c('t','z','dZ','ddZ')])
		pI = predSpace[,c('tmean.climatology','z','dZ','ddZ')]
		names(pI)[1] = 't'
		bcp = predict(bc,pI) 

		xyz = cbind(baseLine,z=bcp)

		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))

		EAmap( xyz, fn="lobster.biocl.pred", loc="output",corners=corners, display=T)


		# by year
		bc = bioclim(dat.p[,c('dyear','t','z','dZ','ddZ')])

		for(i in 1:length(Years)){
			pdat = subset(predSpaceTime,year(timestamp)==Years[i])
			pdat$dyear = decimal_date(pdat$timestamp)
			pI = cbind(pdat[,c('dyear','t')],predSpace[,c('z','dZ','ddZ')])
			bcp = predict(bc,pI) 
			xyz = cbind(baseLine,z=bcp)
			corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
			EAmap( xyz, fn=paste("lobster.biocl.pred",Years[i],sep='.'), annot=Years[i],loc="output",corners=corners, display=T)
		}



	# GAM (generalized additive models)
	
	# Habitat models	
		
		dat = lobsterRV
		dat = scallopRV
		dat = silverhakeRV
		dat = halibutRV
		dat$Y = ifelse(dat$totwgt>0,1,0)
		dat$dyear = decimal_date(dat$sdate)
		dat$t = dat$bottom_temperature
	
	##  with space

	    Mf = formula( Y ~ s(t) +  s(z) + s(dZ) + s(ddZ) + s(plon, plat)  )

	    Md = dat[,c('plon','plat','t','z','dZ','ddZ','Y')]

		Mo = gam( Mf, data=Md, family=binomial())
		summary(Mo)


		pI = predSpace[,c('plon','plat','tmean.climatology','z','dZ','ddZ')]
		names(pI)[3] = 't'
		bcp = predict(Mo,pI,type='response') 
		xyz = cbind(baseLine,z=bcp)
		corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
		EAmap( xyz, fn="lobster.gambi.pred", loc="output",corners=corners, display=T)

	####################################
	# 
	#	to incorporate into stock assessment: Use habitat layer to create strata


	##  with space and time 

	    Mf = formula( Y ~ s(dyear) + s(t) +  s(z) + s(dZ) + s(ddZ) + s(plon, plat)  )

	    Md = dat[,c('plon','plat','dyear','t','z','dZ','ddZ','Y')]

		Mo = gam( Mf, data=Md, family=binomial())
		summary(Mo)


		for(i in 1:length(Years)){
			pdat = subset(predSpaceTime,year(timestamp)==Years[i])
			pdat$dyear = decimal_date(pdat$timestamp)
			pI = cbind(pdat[,c('plon','plat','dyear','t')],predSpace[,c('z','dZ','ddZ')])
			bcp = predict(Mo,pI,type='response') 
			xyz = cbind(baseLine,z=bcp)
			corners = data.frame(lon=c(-67.54,-56.5),lat=c(41,47.2))
			EAmap( xyz, fn=paste("lobster.gambi.pred",Years[i],sep='.'), annot=Years[i],loc="output",corners=corners, display=T)
		}








##  Spatial Abundance models


	modelformula = formula( paste( 
	      varname, ' ~ s(t, k=3, bs="ts") + s(tmean.climatology, k=3, bs="ts") + s(tsd.climatology, k=3, bs="ts")  ', 
	      ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
	      ' + s( log(mr), k=3, bs="ts") + s( Npred, k=3, bs="ts") + s( smr, k=3, bs="ts")  ',
	      ' + s(log.substrate.grainsize, k=3, bs="ts") + s(ca1, k=3, bs="ts") + s(ca2, k=3, bs="ts")   ' ))  # no space 



##  Spatial-temporal population dynamics