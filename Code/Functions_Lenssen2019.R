###############################################################################
###############################################################################
# Functions from the Lenssen et al. (2019) analysis. As many of these are used
# as possible in the analysis, but some are superceeded by new functions.

# NOTE: Load this function script prior to the Functions.R script!

# GISTEMP Uncertainty Ensemble
# Version 1.0.0 (August 21, 2024)
# Nathan Lenssen (lenssen@mines.edu)
# https://data.giss.nasa.gov/gistemp/
###############################################################################
###############################################################################

###############################################################################
# Functions used in ProcessMERRA.R
###############################################################################

zoneAssign <- function(lon,lat){
	latBands <- c(-90, -64.2, -44.4, -23.6, 0, 23.6, 44.4, 64.2, 90)

	latInt <- which(diff(findInterval(latBands,lat))==1)

	if(length(latInt)>0){
		latBand <- latBands[latInt:(latInt+1)]
	} else{
		latBand <- latBands[1:2]
		latInt <- 1
	}
	

	if(latInt %in% c(1,8)) lonBands <- seq(-180,180,length=5)
	if(latInt %in% c(2,7)) lonBands <- seq(-180,180,length=9)
	if(latInt %in% c(3,6)) lonBands <- seq(-180,180,length=13)
	if(latInt %in% c(4,5)) lonBands <- seq(-180,180,length=17)
	
	lonInt <- which(diff(findInterval(lonBands,lon))==1)

	if(length(lonInt)>0){
		lonBand <- lonBands[lonInt:(lonInt+1)]
	} else{
		lonBand <- lonBands[1:2]
		lonInt <- 1
	}
	# now create a crude system to convert the bands into the numbering
	# system used in the 87 paper

	latVals <- c(1,5,13,25,41,57,69,77)

	zone <- latVals[latInt] + (lonInt-1)
	
	return(zone)
}


###############################################################################
# Functions used in ProcessStationData.R 
###############################################################################
parseDataString <- function(dataString){
	sapply(seq(from=1, to=nchar(dataString), by=8), function(i) as.numeric(substr(dataString, i, i+4)))
}


###############################################################################
# Functions used in CreateMask.R 
###############################################################################
coverageCheck <- function(station,decade,id){
	yearsPresent <- station[,1]

	yearInd <- rep(NA, 10)

	# loop over the years
	for(i in 2:11){

		# Initialize a temporary meterological year (Dec-Nov) with NA
		metYearTemp <- rep(NA,12)

		# Populate the (year-1) december if year in data
		if(decade[i-1] %in% yearsPresent){
				metYearTemp[1] <- station[station[,1] == decade[i-1],13]
			} 

		# Populate the remainder of the year	
		if(decade[i] %in% yearsPresent){
				metYearTemp[-1] <- station[station[,1] == decade[i],2:12]
			} 
		
		monthIndicator <- !is.na(metYearTemp)

		seasonInd <- rep(NA,4)

		for(s in 1:4){
			seasonInd[s] <- sum(monthIndicator[(1+(s-1)*3 ): (3+(s-1)*3)]) > 1
		}

		yearInd[i-1] <- sum(seasonInd) >2
	}

	out <- data.frame(as.character(substr(id,1,11)),sum(yearInd)>4)
	names(out) <- c("id", "coverage")

	return(out)
}



# This function builds a strict mask using only locations that have
# a station inside their grid box
merraGridMask <- function(coveredStations,merraLon,merraLat){
	mask <- matrix(0,nrow=length(merraLon),ncol=length(merraLat))

	for(i in 1:nrow(coveredStations)){
		singleStation <- coveredStations[i,]
		lonInd <- which.min(abs(singleStation$lon - merraLon))
		latInd <- which.min(abs(singleStation$lat - merraLat))

		mask[lonInd,latInd] <- 1
	}

	return(mask)
}


###############################################################################
# Functions used in seaIceMask.R 
###############################################################################

landMaskNewGrid <- function(reanalysis,lonRef,latRef,seasonalLandMaskRef,saving=TRUE){
	# load in the new grid
	load(file=sprintf('Data/%s/anomalyData_%s.Rda',reanalysis,reanalysis))
	lon <- anomalyData$lon
	lat <- anomalyData$lat
	rm(anomalyData)

	# make the new mask
	newMask <- array(NA, dim=c(length(lon),length(lat),12))

	for(sInd in 1:12){
		for(i in 1:length(lon)){
			for(j in 1:length(lat)){
				lonInd <- which.min(abs(lon[i] - lonRef))	
				latInd <- which.min(abs(lat[j] - latRef))

				newMask[i,j,sInd] <- seasonalLandMaskRef[lonInd,latInd,sInd]
			}
		}
	}

	# rename, arrange and save to disk
	seasonalLandMask <- newMask
	maximalLandMask <- ifelse(apply(seasonalLandMask,c(1,2),function(x) any(!is.na(x))),1,NA)

	if(saving) save(seasonalLandMask,maximalLandMask,file=sprintf('Data/%s/landMasks_%s.Rda',reanalysis,reanalysis))

	return(list(seasonalLandMask=seasonalLandMask,maximalLandMask=maximalLandMask))
}




###############################################################################
# Functions used in interpolation
###############################################################################
interpolateFieldPointwise <- function(decadalMask,radius,nCores,merraData,landMask){
	# unpack the list for more readable code
	lon <- merraData$lon
	nlon <- length(lon)
	lat <- merraData$lat
	nlat <- length(lat)
	nt <- dim(merraData$anomalyField)[3]

	# make the anomaly data a matrix with non-land removed
	tempAnomaly <- apply(merraData$anomalyField,3L,c)

	goodInds    <- which(!is.na(landMask) | decadalMask==1)
	goodIndsMat <- which(!is.na(landMask) | decadalMask==1,arr.ind=T)

	anomalyMat <- tempAnomaly[goodInds,]

	# get the matrix locations of the tables
	stationIndsMat <- which(decadalMask==1,arr.ind=T)

	partialMask <- c(decadalMask)[goodInds]
	stationInds <- which(partialMask==1)

	# using parallel do with foreach
	fnExport  <- c('singlePoint','rdist.earth','triangleRBFVec')
	varExport <- c('radius','lon','lat','goodIndsMat','stationInds','stationIndsMat','anomalyMat')

	if(nCores > 1){
		cl <- makeCluster(nCores)
		registerDoParallel(cl)
		interpolatedMatPar <- foreach(i = 1:nrow(goodIndsMat),
									 .export=c(fnExport,varExport)) %dopar% singlePoint(i)
		stopCluster(cl)
		gc()
	} else{
		interpolatedMatPar <- foreach(i = 1:nrow(goodIndsMat),
									 .export=c(fnExport,varExport)) %do% singlePoint(i)
	}


	interpolatedField <- array(NA,dim=dim(merraData$anomalyField))

	for(i in 1:nrow(goodIndsMat)){
		interpolatedField[goodIndsMat[i,1],goodIndsMat[i,2],] <- interpolatedMatPar[[i]]
	}


	return(interpolatedField)
	
}



interpolateFieldPointwiseLowMem <- function(radius,lon,lat,goodIndsMat,stationInds,stationIndsMat,anomalyMat){

	# using parallel do with foreach
	fnExport  <- c('singlePoint','rdist.earth','triangleRBFVec')
	varExport <- c('radius','lon','lat','goodIndsMat','stationInds','stationIndsMat','anomalyMat')

	if(nCores > 1){
		cl <- makeCluster(nCores)
		registerDoParallel(cl)
		interpolatedMatPar <- foreach(i = 1:nrow(goodIndsMat),
									 .export=c(fnExport,varExport)) %dopar% singlePoint(i)
		stopCluster(cl)
		gc()
	} else{
		interpolatedMatPar <- foreach(i = 1:nrow(goodIndsMat),
									 .export=c(fnExport,varExport)) %do% singlePoint(i)
	}


	interpolatedField <- array(NA,dim=c(length(lon),length(lat),ncol(anomalyMat)))

	for(i in 1:nrow(goodIndsMat)){
		interpolatedField[goodIndsMat[i,1],goodIndsMat[i,2],] <- interpolatedMatPar[[i]]
	}


	return(interpolatedField)
	
}


# helper functions
singlePoint <- function(i){
	dists <- rdist.earth(matrix(c(lon[goodIndsMat[i,1]],lat[goodIndsMat[i,2]]),nrow=1,ncol=2),
						cbind(lon[stationIndsMat[,1]],lat[stationIndsMat[,2]]),miles=FALSE)[1,]

	distInds <- which(dists < radius)
	
	# return NA series if no stations in range
	if(length(distInds)==0) return(rep(NA,ncol(anomalyMat)))

	# If stations, calcuate the weighted average 
	weights <- triangleRBFVec(dists[distInds],radius)
	tempStationInds <- stationInds[distInds]

	subAnomalyMat <- matrix(anomalyMat[tempStationInds,],ncol=ncol(anomalyMat))
	return(apply(subAnomalyMat,2,weighted.mean,w=weights,na.rm=T))
}


triangleRBFVec <- function(dist,radius){
	mat <- cbind(rep(0,length(dist)),1 - dist/radius)
	return(apply(mat,1,max))
}


# new function for the ensemble analysis (parallize the 14 at once and run each serially)
interpolateFieldSerial <- function(radius,lon,lat,goodIndsMat,stationInds,stationIndsMat,anomalyMat){
	require(fields)

	# functions needed
	triangleRBFVec <- function(dist,radius){
		mat <- cbind(rep(0,length(dist)),1 - dist/radius)
		return(apply(mat,1,max))
	}
	
	singlePoint <- function(i){
		dists <- rdist.earth(matrix(c(lon[goodIndsMat[i,1]],lat[goodIndsMat[i,2]]),nrow=1,ncol=2),
							cbind(lon[stationIndsMat[,1]],lat[stationIndsMat[,2]]),miles=FALSE)[1,]

		distInds <- which(dists < radius)
		
		# return NA series if no stations in range
		if(length(distInds)==0) return(rep(NA,ncol(anomalyMat)))

		# If stations, calcuate the weighted average 
		weights <- triangleRBFVec(dists[distInds],radius)
		tempStationInds <- stationInds[distInds]

		subAnomalyMat <- matrix(anomalyMat[tempStationInds,],ncol=ncol(anomalyMat))
		return(apply(subAnomalyMat,2,weighted.mean,w=weights,na.rm=T))
	}




	# loop over the single points
	interpolatedMatPar <- list()

	for(i in 1:nrow(goodIndsMat)){
		interpolatedMatPar[[i]] <- singlePoint(i)
	}


	interpolatedField <- array(NA,dim=c(length(lon),length(lat),ncol(anomalyMat)))



	for(i in 1:nrow(goodIndsMat)){
		interpolatedField[goodIndsMat[i,1],goodIndsMat[i,2],] <- interpolatedMatPar[[i]]
	}


	return(interpolatedField)
	
}


###############################################################################
# Functions used in Mean Calculations
###############################################################################
globalMean <- function(anomalyField,mask=NULL,lat,nCores=2,zoneMask,simpleMean=FALSE, ALn, ALs){

	require(foreach)
	require(doParallel)

	nt <- dim(anomalyField)[3]

	meanFull  <- rep(NA, nt)
	meanNorth <- rep(NA, nt)
	meanSouth <- rep(NA, nt)
	meanBands <- matrix(NA,nrow=nt,ncol=8)

	fnExport <- c('gistempGlobalMeanUpdated','meanSlice','explicitMean')
	varExport <- c('mask','anomalyField','lat','zoneMask')

	if(nCores > 1){
		cl <- makeCluster(nCores)
		registerDoParallel(cl)
		outList <- foreach(time=1:nt, .export=c(fnExport)) %dopar% meanSlice(time,anomalyField,mask,lat,zoneMask,simpleMean, ALn, ALs)

		stopCluster(cl)
		gc()
	} else{
		outList <- foreach(time=1:nt, .export=c(fnExport)) %do% meanSlice(time,anomalyField,mask,lat,zoneMask,simpleMean, ALn, ALs)
		gc()
	}



	for(time in 1:nt){
		meanFull[time]  <- outList[[time]]$global
		meanNorth[time] <- outList[[time]]$nh
		meanSouth[time] <- outList[[time]]$sh
		meanBands[time,] <- outList[[time]]$bands
	}

	outList <- list(global=meanFull, nh=meanNorth, sh=meanSouth, bands=meanBands)

	return(outList)
}

meanSlice <- function(time,anomalyField,mask,lat,zoneMask,simpleMean=FALSE, ALn, ALs){
	if(is.null(mask)){
		trueField <- anomalyField[,,time]
	} else{
		trueField <- anomalyField[,,time] * mask
	}
	# take the mean
	if(simpleMean){
		outObj <- explicitMean(trueField,lat)
	} else{
		outObj <- gistempGlobalMeanUpdated(trueField,lat,zoneMask, ALn, ALs)
	}

	return(outObj)
}

gistempGlobalMeanUpdated <- function(fieldSlice,lat,zoneMask, ALn, ALs){
	nx <- dim(fieldSlice)[1]
	ny <- dim(fieldSlice)[2]

	latVals <- c(1,5,13,25,41,57,69,77,81)
	
	# create a weight mask
	weights    <- cos(lat * (pi/180))
	weightMask <- matrix(rep(weights,nx),nx,ny,byrow=T)

	bandWeight <- rep(NA, length(latVals)-1)
	bandMeans <- rep(NA, length(latVals)-1)
	zoneMeans <- rep(NA,80)

	for(band in 1:(length(latVals)-1)){
		boxInds <- latVals[band]:(latVals[band+1]-1)
		boxMeans <- rep(NA, length(boxInds))
		boxWeights <- rep(NA, length(boxInds))

		for(i in 1:length(boxInds)){
			#First mask the data for only the specific box
			vec <- rep(NA, nx*ny)
			vec[ which( zoneMask==boxInds[i] & !is.na(c(fieldSlice)) ) ] <- 1
			subSlice   <- fieldSlice * matrix(vec, nrow=nx, ncol=ny)

			subWeights <- weightMask * matrix(vec, nrow=nx, ncol=ny)
			# weight by the appropriate cos(lat)

			summedData <- sum(subSlice*weightMask,na.rm=TRUE)
			summedWeights <- sum(subWeights,na.rm=TRUE)

			boxMeans[i] <- summedData/summedWeights


			# see the prop of gridboxes that have data to properly
			# weight the band mean
			boxWeights[i] <- sum(!is.na(subSlice))/sum(zoneMask==boxInds[i],na.rm=TRUE)
		}
		zoneMeans[latVals[band]:(latVals[band]+length(boxMeans)-1)] <- boxMeans
		
		bandWeight[band] <- sum(boxWeights,na.rm=TRUE)
		bandMeans[band] <- sum(boxMeans*boxWeights,na.rm=TRUE)/bandWeight[band]
	}

	# take the data-avaiablility weighted means of the boxes within a band
	weightedBandMeans <- bandMeans * bandWeight

	nHemisMean <- sum(weightedBandMeans[5:8],na.rm=T)/sum(bandWeight[5:8]) 
	sHemisMean <- sum(weightedBandMeans[1:4],na.rm=T)/sum(bandWeight[1:4])

	# area weighted mean of hemispheres (following GISTEMP algorithm)
	globalMean <- sum(ALn * nHemisMean + ALs * sHemisMean)/(ALn + ALs) 

	outList <- list(global=globalMean,nh=nHemisMean,sh=sHemisMean,bands=bandMeans)
	return(outList)
}

mergeFields <- function(landField,merraData){
	lon <- merraData$lon
	lat <- merraData$lat
	anomalyField <- merraData$anomalyField
	landMask <- merraData$landMask

	# make anti-land Mask 
	oceanMask <- ifelse(is.na(landMask),1,0)

	fullField <- array(NA, dim=dim(landField))

	for(i in 1:dim(landField)[3]){
		landSlice <- landField[,,i]

		landSlice[is.na(landSlice) & oceanMask==1] <- 0
		fullField[,,i] <- landSlice + anomalyField[,,i]*oceanMask
	}

	return(fullField)
}

explicitMean <- function(fieldSlice,lat){
	naMask <- ifelse(is.na(fieldSlice),NA,1)
	weightMat <- matrix(cos(lat* pi/180), nrow=dim(fieldSlice)[1], ncol=length(lat),byrow=T)

	globalMean <- sum(fieldSlice * weightMat,na.rm=TRUE)/sum(naMask * weightMat,na.rm=TRUE)

	northCols  <- which(lat>0)
	southCols  <- which(lat<0)
	
	nHemisMean <- sum(fieldSlice[,northCols] * weightMat[,northCols],na.rm=TRUE)/
						sum(naMask[,northCols] * weightMat[,northCols],na.rm=TRUE)
	sHemisMean <- sum(fieldSlice[,southCols] * weightMat[,southCols],na.rm=TRUE)/
						sum(naMask[,southCols] * weightMat[,southCols],na.rm=TRUE)

	bandMeans <- rep(NA, 8)
	latBands <- c(-90, -64.2, -44.4, -23.6, 0, 23.6, 44.4, 64.2, 90)

	for(i in 1:length(bandMeans)){
		bandCols <- which(lat >= latBands[i] & lat <= latBands[i+1])
		bandMeans[i] <- sum(fieldSlice[,bandCols] * weightMat[,bandCols],na.rm=TRUE)/
						 sum(naMask[,bandCols] * weightMat[,bandCols],na.rm=TRUE)
	}

	outList <- list(global=globalMean,nh=nHemisMean,sh=sHemisMean,bands=bandMeans)
}

###############################################################################
# Functions used in Plotting
###############################################################################
totalLandUncertainty <- function(resultsList,tDec,normalVal=1.96){
	# global mean calculations
	trueGlobal <- monthToYearMeans(resultsList$true$global)

	maskGlobal <- matrix(NA, nrow=length(tDec), ncol=length(trueGlobal))

	for(i in 1:nrow(maskGlobal)) maskGlobal[i,] <- monthToYearMeans(resultsList$mask[[i]]$global)

	# get the beta and sigma estimates
	sigma <- rep(NA, length(tDec))
	beta  <- rep(NA, length(tDec))
	beta2 <- rep(NA, length(tDec))

	se     <- rep(NA, length(tDec))
	se2    <- rep(NA, length(tDec))

	variance  <- rep(NA, length(tDec))
	variance2 <- rep(NA, length(tDec))

	diffVar  <- rep(NA, length(tDec))

	for(i in 1:length(tDec)){
		obj <- differenceSeriesLM(trueGlobal,maskGlobal[i,],plotting=F)
		sigma[i] <- obj$lmSD
		beta[i]  <- obj$slope
		beta2[i] <- obj$slope2

		se[i]    <- obj$se
		se2[i]   <- obj$se2

		variance[i]  <- obj$se^2
		variance2[i] <- obj$se2^2

		diffVar[i] <- var(trueGlobal - maskGlobal[i,])
	}


	globalLandCI <- normalVal * sigma

	return(list(ci=globalLandCI,sigma=sigma,beta=beta,se=se,variance=variance,diffVar=diffVar,
				beta2 = beta2, se2 = se2, variance2 = variance2))
}

# NEED TO REWRITE OUR DECADE TO YEAR FUNCTION TO HANDLE GENERAL CASES
decadeToYearNA <- function(series,tYear,tDec){
	yearMeans <- rep(NA, length(tYear))
	
	# handle the first decade
	tempInds <- which(tYear < tDec[2])
	lenFirst <- length(tempInds)
	yearMeans[1:(lenFirst-1)] <- series[1]
	

	# handle the middle decades
	for(i in 1:(length(tDec)-2)){
		yearMeans[((lenFirst+1)+(i-1)*10):((lenFirst)+(i)*10-1)] <- series[i+1]
	}

	# handle the last decade
	yearMeans[((lenFirst+1)+(length(tDec)-2)*10):length(tYear)] <- series[length(tDec)]
	return(yearMeans)
}

decadeToYear <- function(series,tYear,tDec){
	
	yearMeans <- rep(NA, length(tYear))
	
	# handle the first decade
	tempInds <- which(tYear < tDec[2])
	yearMeans[tempInds] <- series[1]
	lenFirst <- length(tempInds)

	# handle the middle decades
	for(i in 1:(length(tDec)-2)){
		yearMeans[((lenFirst+1)+(i-1)*10):((lenFirst)+(i)*10)] <- series[i+1]
	}

	# handle the last decade
	yearMeans[((lenFirst+1)+(length(tDec)-2)*10):length(tYear)] <- series[length(tDec)]
	return(yearMeans)
}

differenceSeriesLM <- function(trueYear, maskYear,plotting=TRUE,...){
	tYear <- 1980:2016
	diffSD <- sd(trueYear-maskYear)
	fit <- lm(trueYear~maskYear)
	fit2 <- lm(maskYear~trueYear)

	lmSD <- summary(fit)$sigma
	lmSD2 <- summary(fit)$sigma

	if(plotting){
		set.panel(1,3)
		yr <- range(c(trueYear,maskYear))
		plot(tYear,trueYear,type='l',col='black',ylim=yr,...)
		points(tYear,maskYear,type='l',col='red')

		plot(tYear,trueYear-maskYear,type='l')
		legend("bottomright",paste("Difference SD:",round(diffSD,3)))

		plot(trueYear~maskYear)
		abline(fit)
		legend("bottomright",paste("Regression SD:",round(lmSD,3)))
	}

	return(list(lmSD=lmSD,diffSD=diffSD,slope=fit$coefficients[2], se=summary(fit)$coefficients[2,2],
				lmSD2 = lmSD2, slope2 = fit2$coefficients[2], se2 = summary(fit2)$coefficients[2,2]))
}

monthToYearMeans <- function(monthMeans){
	yearMeans <- rep(NA, length(monthMeans)/12)

	# take yearly means
	for(year in 1:length(yearMeans)){
		yearInds <- (1+(year-1)*12):(year*12)
		yearMeans[year] <- mean(monthMeans[yearInds],na.rm=TRUE)
	}

	return(yearMeans)
}
###############################################################################
# From the ERSST visFunctions.R (used to extract the ERSST results from 
# the 1000 ensemble runs)
###############################################################################

extratTSRdat <- function(ensemble,ddir,tYear){
	annMat <- matrix(NA,nrow=length(tYear),ncol=length(ensemble))

	for(i in 1:length(ensemble)){
		load(sprintf("%s/%s",ddir,ensemble[i]))
		annMat[,i] <- monthToYearMeans(oceanMean$global) 
	}

	return(annMat)
}

extratTSRdatNH <- function(ensemble,ddir,tYear){
	annMat <- matrix(NA,nrow=length(tYear),ncol=length(ensemble))

	for(i in 1:length(ensemble)){
		load(sprintf("%s/%s",ddir,ensemble[i]))
		annMat[,i] <- monthToYearMeans(oceanMean$nh) 
	}

	return(annMat)
}

extratTSRdatSH <- function(ensemble,ddir,tYear){
	annMat <- matrix(NA,nrow=length(tYear),ncol=length(ensemble))

	for(i in 1:length(ensemble)){
		load(sprintf("%s/%s",ddir,ensemble[i]))
		annMat[,i] <- monthToYearMeans(oceanMean$sh) 
	}

	return(annMat)
}


###############################################################################
# Function used to perform calcualations in section 6
###############################################################################

resultsMC <- function(meanSeries,ciSeries,tYear,alpha=0.5,M=1000){
	randWarmest <- rep(NA,M)
	arWarmest   <- rep(NA,M)

	# MC loop
	for(i in 1:M){
		randDraw <- rnorm(length(tYear),mean=meanSeries,
										sd=ciSeries/1.96)
		
		arDraw   <- meanSeries + arima.sim(model=list(ar=alpha),
			n=length(tYear), mean = 0, sd=ciSeries/1.96)

		# get which year is warmest year
		randWarmest[i] <- tYear[which.max(randDraw)]
		arWarmest[i]   <- tYear[which.max(arDraw)]
	}

	# return the warmest years
	return(list(ind = table(randWarmest)/M, ar = table(arWarmest)/M)
	)
}


###############################################################################
# Function used to read the ERSST binary files
###############################################################################

readGridded <- function(ddir,fileName){
	# make the 2x2 grid
	nx <- 180
	ny <- 89
	nt <- 12*(2014-1854+1)


	shift <- 0 # Needed to get the images to align with boundries. Should check on this
	rawLon <- c(seq(0+shift,180,by=2),seq(-178,-2+shift,by=2))
	perm   <- rank(rawLon)
	lon <- seq(-178,180,by=2)
	lat <- seq(-88,88,by=2)

	# The full ssta array for a given sample
	datArray <- array(NA, dim=c(nx,ny,nt))
	# the data is presented in 4 byte floats (big endian) with start and end caps (ints)
	# open the binary file stream in read binary mode (rb)
	fileStream <- file(sprintf("%s/%s",ddir, fileName),open="rb")

	for(i in 1:nt){
		startCap <- readBin(fileStream,integer(),size=4,n=1,endian='big')
		datRaw   <- readBin(fileStream,numeric(),size=4,n=(nx*ny),endian='big')
		endCap   <- readBin(fileStream,integer(),size=4,n=1,endian='big')
		datMat   <- matrix(NA, nx,ny)
		datMat[perm,] <- matrix(datRaw,nx,ny)
		datMat[datMat < -900] <- NA

		datArray[,,i] <- datMat
	}

	close(fileStream)

	return(list(lon=lon,lat=lat,sst=datArray))
}


