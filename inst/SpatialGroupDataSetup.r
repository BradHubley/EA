#SpatialGroupDataSetup.r

# set wd
setwd("/home/hubleyb/bio/EA/")

#https://github.com/jae0
#https://github.com/Beothuk


# load R packages
require(bio.groundfish)
require(bio.snowcrab)
require(bio.indicators)
require(bio.spacetime)

## Survey Data

	## snowcrab
	#p = bio.snowcrab::load.environment( year.assessment=2016 )
	#snowcrabENS = snowcrab.db( DS ="set.complete", p=p )
	#write.csv(snowcrabENS,"data/snowcrabENS.csv")

	### lobster
	#lobster34 = read.csv("data/lobster34.csv")
	### set space-time coordinates
	#lobster34 = lonlat2planar(lobster34,"utm20",input_names=c("SET_LONG", "SET_LAT"))
	#lobster34$timestamp = as.Date(lobster34[,"SET_DATE"])

	### scallop
	#scallopBBn = read.csv("data/BBnSurvey9114.csv")
	### set space-time coordinates
	#scallopBBn = lonlat2planar(scallopBBn,"utm20")
	#scallopBBn$timestamp = as.Date(paste(scallopBBn$year,"06-01",sep='-'))

	## rv survey
	p = bio.groundfish::load.groundfish.environment( "bio.survey")

	p$strat=440:495
	p$series =c('summer')# p$series =c('4vswcod');p$series =c('georges')
	p$years.to.estimate = c(1999:2014)
	p$vessel.correction = T
	p$vessel.correction.fixed = 1.2
	p$length.based = T
	p$size.class= c(1,1000)
	p$by.length= T
	p$by.sex = F
	p$strata.efficiencies=F
	p$bootstrapped.ci=F
	p$strata.files.return=T
	p$functional.groups = F
	p$clusters = c( rep( "localhost", 7) )
	
	p$species = 2550 #lobster
	p = make.list(list(v=p$species, yrs=p$years.to.estimate),Y=p)
	p$runs = p$runs[order(p$runs$v),]
	aout= groundfish.analysis(DS='stratified.estimates.redo',p=p)
	st.out = groundfish.analysis(DS='species.set.data',p=p)
    al = lapply(st.out,"[[",2)
    lobsterRV = do.call('rbind',al)
    lobsterRV$slong = lobsterRV$slong * -1
	lobsterRV = lonlat2planar(lobsterRV,"utm20", input_names=c("slong", "slat"))

	p$species = 14 #silver hake
	p = make.list(list(v=p$species, yrs=p$years.to.estimate),Y=p)
	p$runs = p$runs[order(p$runs$v),]
	aout= groundfish.analysis(DS='stratified.estimates.redo',p=p)
	st.out = groundfish.analysis(DS='species.set.data',p=p)
    al = lapply(st.out,"[[",2)
    silverhakeRV = do.call('rbind',al)
    silverhakeRV$slong = silverhakeRV$slong * -1
	silverhakeRV = lonlat2planar(silverhakeRV,"utm20", input_names=c("slong", "slat"))

	p$species =  4321 #scallop
	p = make.list(list(v=p$species, yrs=p$years.to.estimate),Y=p)
	p$runs = p$runs[order(p$runs$v),]
	aout= groundfish.analysis(DS='stratified.estimates.redo',p=p)
	st.out = groundfish.analysis(DS='species.set.data',p=p)
    al = lapply(st.out,"[[",2)
    scallopRV = do.call('rbind',al)
    scallopRV$slong = scallopRV$slong * -1
	scallopRV = lonlat2planar(scallopRV,"utm20", input_names=c("slong", "slat"))

	p$species =  30 #halibut
	p = make.list(list(v=p$species, yrs=p$years.to.estimate),Y=p)
	p$runs = p$runs[order(p$runs$v),]
	aout= groundfish.analysis(DS='stratified.estimates.redo',p=p)
	st.out = groundfish.analysis(DS='species.set.data',p=p)
    al = lapply(st.out,"[[",2)
    halibutRV = do.call('rbind',al)
    halibutRV$slong = halibutRV$slong * -1
	halibutRV = lonlat2planar(halibutRV,"utm20", input_names=c("slong", "slat"))





## Environmental Data

	p = spatial_parameters( type = "SSE" )  
	p$yrs = 1970:2016
    p$nw = 10 # number of intervals in time within a year
    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks

	# baseline grid use as a predictive surface
	baseLine = bathymetry.db(p=p, DS="baseline")
	save(baseLine,file= "data/baseline.rdata")

	# depth contours for plotting
	isoBaths = isobath.db( p=p, depths=seq( 50, 1000, 50), crs=p$internal.crs )
	save(isoBaths,file= "data/isobaths.rdata")

    # coastline data for plotting
    coastLine = coastline.db(p=p, crs=p$internal.crs)
	save(coastLine,file= "data/coastline.rdata")

    # spatial only (no time component) varialbes for predicting
    predSpace = indicators.lookup( p=p, DS="spatial", locs=baseLine)
	save(predSpace,file= "data/predspace.rdata")

	times = as.Date(paste(1999:2014,"06-01",sep='-'))
	preddata = merge(baseLine,data.frame(timestamp=times),all=T)

    predSpaceTime = indicators.lookup( p=p, DS="spatial.annual", locs=baseLine, timestamp=preddata$timestamp)
 	predSpaceTime = cbind(preddata,predSpaceTime)
 	save(predSpaceTime,file= "data/predspacetime.rdata")
   



## Joining Fishery and Environmental Data
setwd("/home/hubleyb/bio/EA/")

	load("data/baseline.rdata") # baseLine
	load("data/isobaths.rdata") # isoBaths
	load("data/coastline.rdata") # coastLine
	load("data/predspace.rdata") # predSpace
	load("data/predspacetime.rdata") # predSpaceTime
	

	## Lobster


    # identify locations of data relative to baseline for envionmental data
    locsmap = match( 
        lbm::array_map( "xy->1", lobsterRV[,c("plon","plat")], gridparams=p$gridparams ), 
        lbm::array_map( "xy->1", baseLine, gridparams=p$gridparams ) )
      
    # assign environmental data based on space-time coordinates
    btemp = indicators.lookup( p=p, DS="spatial.annual.seasonal", locsmap=locsmap, timestamp=as.Date(lobsterRV$sdate), varnames="tmean" )
    spatial.annual = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=as.Date(lobsterRV$sdate))
    spatial = indicators.lookup( p=p, DS="spatial", locsmap=locsmap)

    # join with dataset and save a csv
    lobsterRV = cbind( lobsterRV[,-which(names(lobsterRV)=='z')],  btemp,  spatial[,-which(names(spatial)%in%c('plon','plat'))],spatial.annual)
	write.csv(lobsterRV,"data/lobsterRV.csv")


	## scallop

    # identify locations of data relative to baseline for envionmental data
    locsmap = match( 
        lbm::array_map( "xy->1", scallopRV[,c("plon","plat")], gridparams=p$gridparams ), 
        lbm::array_map( "xy->1", baseLine, gridparams=p$gridparams ) )
      
    # assign environmental data based on space-time coordinates
    btemp = indicators.lookup( p=p, DS="spatial.annual.seasonal", locsmap=locsmap, timestamp=as.Date(scallopRV$sdate), varnames="tmean" )
    spatial.annual = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=as.Date(scallopRV$sdate))
    spatial = indicators.lookup( p=p, DS="spatial", locsmap=locsmap)
      

    # join with dataset and save a csv
    scallopRV = cbind( scallopRV[,-which(names(scallopRV)=='z')],  btemp,  spatial[,-which(names(spatial)%in%c('plon','plat'))], spatial.annual)
	write.csv(scallopRV,"data/scallopRV.csv")


	## silverhake

    # identify locations of data relative to baseline for envionmental data
    locsmap = match( 
        lbm::array_map( "xy->1", silverhakeRV[,c("plon","plat")], gridparams=p$gridparams ), 
        lbm::array_map( "xy->1", baseLine, gridparams=p$gridparams ) )
      
    # assign environmental data based on space-time coordinates
    btemp = indicators.lookup( p=p, DS="spatial.annual.seasonal", locsmap=locsmap, timestamp=as.Date(silverhakeRV$sdate), varnames="tmean" )
    spatial.annual = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=as.Date(silverhakeRV$sdate))
    spatial = indicators.lookup( p=p, DS="spatial", locsmap=locsmap)
      

    # join with dataset and save a csv
    silverhakeRV = cbind( silverhakeRV[,-which(names(silverhakeRV)=='z')],  btemp,  spatial[,-which(names(spatial)%in%c('plon','plat'))], spatial.annual)
	write.csv(silverhakeRV,"data/silverhakeRV.csv")


	## halibut

    # identify locations of data relative to baseline for envionmental data
    locsmap = match( 
        lbm::array_map( "xy->1", halibutRV[,c("plon","plat")], gridparams=p$gridparams ), 
        lbm::array_map( "xy->1", baseLine, gridparams=p$gridparams ) )
      
    # assign environmental data based on space-time coordinates
    btemp = indicators.lookup( p=p, DS="spatial.annual.seasonal", locsmap=locsmap, timestamp=as.Date(halibutRV$sdate), varnames="tmean" )
    spatial.annual = indicators.lookup( p=p, DS="spatial.annual", locsmap=locsmap, timestamp=as.Date(halibutRV$sdate))
    spatial = indicators.lookup( p=p, DS="spatial", locsmap=locsmap)
      

    # join with dataset and save a csv
    halibutRV = cbind( halibutRV[,-which(names(halibutRV)=='z')],  btemp,  spatial[,-which(names(spatial)%in%c('plon','plat'))], spatial.annual)
	write.csv(halibutRV,"data/halibutRV.csv")




## Simulated Abundance Data

	Dim = c("n_x"=100, "n_y"=100)
	loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
	Scale = 10
	Sigma2 = (0.5) ^ 2

	# Simulate spatial process
	RMmodel = RMexp(var=Sigma2, scale=Scale)
	epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)
	image( z=epsilon_xy )

	beta0 = 3
	prob_missing = 0.2

	# SImulate counts
	c_xy = array(NA, dim=dim(epsilon_xy))
	for(x in 1:nrow(c_xy)){
	for(y in 1:ncol(c_xy)){
	  c_xy[x,y] = rpois(1, exp(beta0 + epsilon_xy[x,y]) )
	  #if( rbinom(n=1, size=1, prob=prob_missing)==1) c_xy[x,y] = NA
	}}
	true_abundance =  exp(beta0 + epsilon_xy) 
