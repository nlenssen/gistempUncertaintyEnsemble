

###############################################################################
# load and preprocess the 2x2 ERA analysis
###############################################################################

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/GistempProduction/coverageInfo.Rda',scratchDir))

# unpack the ERA list
lon <- anomalyData$lon
lat <- anomalyData$lat
nlon <- length(lon)
nlat <- length(lat)
timeMapERA <- anomalyData$timeMap
anomalyField <- anomalyData$anomalyField
zoneMask <- anomalyData$zoneMask

rm(anomalyData)

# time stuff
nt <- dim(anomalyField)[3]

# land-mask the annualField (could eventually make this a monthly mask)
landMask <- landMaskList$maximalMask

for(i in 1:12){
	subInds  <- which(timeMapERA[,2]==i)
	for(k in subInds){
		anomalyField[,,k] <- anomalyField[,,k] * landMask # drop in monthly mask here indexed by i
	}
}
 



###############################################################################
# Generate the difference series for each decade of interest
###############################################################################

# Get the interpolated field paths
files <- system(sprintf('ls %s/Intermediate/ERA5/InterpolatedFieldsNetcdf/interpField*_2x2.nc',
	scratchDir),intern=T)


fullIndList <- make.surface.grid(list(lon=lon,lat=lat))
landIndList <- fullIndList[which(landMask==1),]
# save the point estimate of the uncertainty variance
uncEstimate <- array(NA, dim=c(length(lon),length(lat),length(decInds)))

diffMatList <- list()
dataIndsList <- list()

for(d in 1:length(decInds)){
	# Load and process the interpolated field
	interpolatedField <- read22File(files[decInds[d]])$anomArray

	# mask by coverage in production gistemp
	for(k in 1:nt){
		interpolatedField[,,k] <- interpolatedField[,,k] * decadalCoverageMask[,,d]
	}

	# place to store all of the difference series
	diffSeriesArr <- array(NA, dim=c(length(lon),length(lat),nt))


	for(i in 1:nrow(landIndList)){
		lonInd <- which(lon==landIndList[i,1])
		latInd <- which(lat==landIndList[i,2])

		trueVec <- anomalyField[lonInd,latInd,]
		maskVec <- interpolatedField[lonInd,latInd,]
		if(all(is.na(trueVec)))	next()
		if(all(is.na(maskVec))){
			if(allLocationsLandUncertainty){
				uncEstimate[lonInd,latInd,d] <- sd(trueVec,na.rm=T)
			} else{
				uncEstimate[lonInd,latInd,d] <- NA
			}
			next()
		}
	
		diffVec <- trueVec - maskVec
		diffSeriesArr[lonInd,latInd,] <- diffVec
		uncEstimate[lonInd,latInd,d] <- sd(diffVec,na.rm=T)

	}

	# store locations where we have data
	dataInds     <- which(!is.na(uncEstimate[,,d]),arr.ind=T)
	numLocations <- nrow(dataInds)

	diffMat <- matrix(NA,nrow=numLocations,ncol=nt)

	for(k in 1:nrow(diffMat)){
		tempInd <- dataInds[k,]
		diffMat[k,] <- diffSeriesArr[tempInd[1],tempInd[2],]
	}

	diffMatList[[d]] <- diffMat 
	dataIndsList[[d]] <- dataInds
}



# save the important output
if(allLocationsLandUncertainty){
	ofname <- sprintf('%s/Intermediate/CovarianceMats/spatialAnalysisFullData_Empirical_allLoc.Rda',scratchDir)
} else{
	ofname <- sprintf('%s/Intermediate/CovarianceMats/spatialAnalysisFullData_Empirical.Rda',scratchDir)
}
save(lon, lat, dataIndsList, diffMatList, uncEstimate,
	file = ofname)
