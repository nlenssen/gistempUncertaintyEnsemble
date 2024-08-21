###############################################################################
###############################################################################
# New or modified functions for the GISTEMP ensemble analysis.
# NOTE: load this after the Lenssen2019 functions.

# GISTEMP Uncertainty Ensemble
# Version 1.0.0 (August 21, 2024)
# Nathan Lenssen (lenssen@mines.edu)
# https://data.giss.nasa.gov/gistemp/
###############################################################################
###############################################################################

###############################################################################
# Functions used in GHCN parsing
###############################################################################

parseDataString <- function(dataString){
	sapply(seq(from=1, to=nchar(dataString), by=9), 
				function(i) as.numeric(substr(dataString, i, i+4)))
}

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
gridMask <- function(coveredStations,lon,lat){
	mask <- matrix(0,nrow=length(lon),ncol=length(lat))

	for(i in 1:nrow(coveredStations)){
		singleStation <- coveredStations[i,]
		lonInd <- which.min(abs(singleStation$lon - lon))
		latInd <- which.min(abs(singleStation$lat - lat))

		mask[lonInd,latInd] <- 1
	}

	return(mask)
}


###############################################################################
# Helper functions for the full ensemble generation 
###############################################################################

maskEnsembleData <- function(arr, timeMap, mask, tDec){
	outArray <- array(NA, dim=dim(arr))

	for(dec in 1:length(tDec)){
		timeInds <- which(timeMap[,1] >= tDec[dec] & timeMap[,1] < tDec[dec] + 10)
		
		# hard coded to account for the code through 2020
		if(dec == 14){
			timeInds <- which(timeMap[,1] >= tDec[dec] & timeMap[,1] < tDec[dec] + 11)
		}

		for(i in 1:length(timeInds)){
			ind <- timeInds[i]
			outArray[,,ind] <- arr[,,ind] * mask[,,timeMap[ind,2],dec]
		}
	}

	return(outArray)
}


###############################################################################
# Functions used for various plots
###############################################################################
addVgram <- function(vlist, linecolor){
	xVals <- vlist[[1]]$centers

	vgramMean <- matrix(NA, length(xVals),length(vlist))

	for(i in 1:length(vlist)){
		vgramMean[,i] <- vlist[[i]]$stats[2,]
	}

	vgramQuants <- apply(vgramMean,1,quantile,prob=c(0.1, 0.25, 0.5, 0.75, 0.9))

	
	points(xVals, vgramQuants[3,], type='l', col=linecolor, lwd=3)

	points(xVals, vgramQuants[2,], type='l', lty=2, col=linecolor, lwd=1.5)
	points(xVals, vgramQuants[4,], type='l', lty=2, col=linecolor, lwd=1.5)

	# points(xVals, vgramQuants[1,], type='l', lty=3, col=linecolor, lwd=1.5)
	# points(xVals, vgramQuants[5,], type='l', lty=3, col=linecolor, lwd=1.5)
}


divMap <- function(x,y,z,zForce=NULL, zCap=NULL, brewerPal = 'Spectral', rev=TRUE,n=256,...){
	require(fields)
	require(RColorBrewer)
	
	if(!is.null(zCap)){
		z[which(z>zCap)] <- zCap
		z[which(z < -zCap)] <- -zCap
	}
	

	zMax <- max(abs(z),zForce,zCap,na.rm=T)
	zr <- c(-zMax,zMax)
	
	pal <- designer.colors(n,brewer.pal(11,brewerPal))

	if(rev) pal <- rev(pal)

	worldMap(x,y,z,zlim=zr,col=pal,...)
}

worldMap <- function(x,y,z,...){
	image.plot(x,y,matrix(z,length(x),length(y)),
		xlab='',ylab='',...)
	world(add=T)
}

plotSeries <- function(countryInd,colInd,plotdir='Figures',yr=NULL,oldSeries=TRUE,verbose=FALSE){
	meanSeries <- ensembleMeanFull[,colInd,countryInd]
	ciMat <- empiricalCIFull[,,colInd,countryInd]

	title <- paste('Annual Met. Mean Temperature Series for',nameVec[countryInd])

	if(is.null(yr)) yr <- range(ciMat,na.rm=T)

	pdf(sprintf('%s/annualSeries_%s.pdf',plotdir,countryCodes[countryInd]),8,4)
	plot(tYear,meanSeries,type='l',ylim=yr,lwd=2,main=title,
		xlab='Year', ylab='Temperature Anomaly (°C)')
	grid(lwd=1.5)
	polygon(c(tYear,rev(tYear)),c(ciMat[1,],rev(ciMat[2,])),
		col=adjustcolor('black',alpha=0.3),
		border=NA)
	if(oldSeries) points(tYear,legacyDatArr[,colInd,countryInd],type='l',col='red',lwd=1.5)
	points(tYear,meanSeries,type='l',lwd=2)
	dev.off()

	if(verbose){
		cat(paste(title,'\n'))
		print(cbind(tYear,apply(ciMat,2,diff)))
	}
}

worldPlot <- function(plottingSeries,zr,barLabel,main='',pal='RdBu',dir=1,cex=15,lsize=10){
  meanTemp2017 <- data.frame(gaulFao=countryCodes, temp=plottingSeries)

  # hardcode fill in a few regions that are difficult to merge
  meanTemp2017 <- rbind(meanTemp2017, c(900, meanTemp2017$temp[meanTemp2017$gaulFao==226]))
  meanTemp2017 <- rbind(meanTemp2017, c(901, meanTemp2017$temp[meanTemp2017$gaulFao==2648]))

  # first merge the tables based on gaul code
  firstJoin <- left_join(meanTemp2017,gaulTab,by='gaulFao')


  # add a few rows manually
  # then join the tables based on ISO as this is carried by the sf dataframe
  names(firstJoin)[3] <- mergeCol
  finalJoin <- left_join(world, firstJoin, by = mergeCol)

  val <- finalJoin$temp

  ggplot() +
    geom_sf(data = finalJoin,aes(fill = val), colour = "transparent") +
    scale_fill_distiller(palette = pal,limits = zr,direction=dir) +

    theme(legend.position = "bottom", 
          legend.key.width = unit(1.5,"cm"),
          legend.title=element_text(size=lsize),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,-5,-5,-5),
          text=element_text(size = cex)) +
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank()) +

    guides(fill = guide_colorbar(title.position = "top")) +

    labs(fill = barLabel,tag = expression(bold(""))) +

    coord_sf(datum = NA) + 

    ggtitle(main)+
  	
  	ylim(-55,90)
}


###############################################################################
# Functions used in GISTEMP Mean Calculations
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
# General Helper Functions
###############################################################################

monthToYearMeans <- function(monthMeans){
	yearMeans <- rep(NA, length(monthMeans)/12)

	# take yearly means
	for(year in 1:length(yearMeans)){
		yearInds <- (1+(year-1)*12):(year*12)
		yearMeans[year] <- mean(monthMeans[yearInds],na.rm=TRUE)
	}

	return(yearMeans)
}

monthToSeasonalMeansRough <- function(monthMeans){
	seasonalMeans <- rep(NA, length(monthMeans)/12)

	# take yearly means
	for(s in 1:length(seasonalMeans)){
		sInds <- (1+(s-1)*3):(s*3)
		seasonalMeans[s] <- mean(monthMeans[sInds],na.rm=TRUE)
	}

	return(seasonalMeans)
}

annualMeanArray <- function(arr){
	nx <- dim(arr)[1]
	ny <- dim(arr)[2]

	outArr <- array(NA, dim=c(dim(arr)[1],dim(arr)[2],dim(arr)[3]/12))

	for(i in 1:nx){
		for(j in 1:ny){
			outArr[i,j,] <- monthToYearMeans(arr[i,j,])
		}
	}

	return(outArr)
}

read22File <- function(loc){
	require(ncdf4)
	handle <- nc_open(loc)
	rawLon <- ncvar_get(handle,'lon')
	lat <- ncvar_get(handle,'lat')

	nlon <- length(rawLon)
	lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

	rawAnom <- ncvar_get(handle,'tempAnom')

	anomArray <- rawAnom[c((nlon/2+1):nlon,1:(nlon/2)),,]
	nc_close(handle)
	return(list(lon=lon,lat=lat,anomArray=anomArray))
}


