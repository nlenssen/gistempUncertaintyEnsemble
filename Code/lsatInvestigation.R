nCores <- 10


###############################################################################
# Make some mask stuff to make each decades coverage consistent for sampling
# uncertainty reasons
###############################################################################

# Make a decadal mask of land locations with coverage (based on GHCN coverage)
decadalLandCoverageArray <- array(NA, dim=c(nlon,nlat,nDec))

for(i in 1:nDec){
	dataMatDec <- dataIndsList[[i]]
	for(j in 1:nrow(dataMatDec)){
		decadalLandCoverageArray[dataMatDec[j,1], dataMatDec[j,2], i] <- 1
	}
}

# Make a monthly mask of ocean locations (based on sea ice)
monthlyOceanMask <- array(NA, dim=c(nlon,nlat,12))

for(i in 1:12){
	monthlyOceanMask[,,i] <- ifelse(is.na(landMaskList$monthlyMask[,,i]), 1, NA)
}

# combine the two to have a decadal, monthly mask of locations with either
# land data (add sampling uncertainty) or ocean data (samping uncertainty
# already contained in ERSST ensemble output)
possibleDataMask <- array(NA, dim=c(nlon,nlat,12,nDec))

for(i in 1:12){
	for(d in 1:nDec){
		possibleDataMask[,,i,d] <- ifelse(!is.na(monthlyOceanMask[,,i]) |
										  !is.na(decadalLandCoverageArray[,,d]), 1, NA)
	}
}

###############################################################################
# Get an estimate of the global homogenization LSAT uncertainty
###############################################################################
filePaths <- system(sprintf('ls %s/Intermediate/GistempLandEnsemble/F*/*.nc',scratchDir),intern=T)[lsatInds] 
nens <- length(filePaths)


# get the land mask info
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

ALn <- landMaskList$landAreaList$ALn
ALs <- landMaskList$landAreaList$ALs

# get the zone mask
load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
zoneMask <- anomalyData$zoneMask
rm(anomalyData)


# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

nlon <- length(lon)
nlat <- length(lat)
nt         <- handle$dim$time$len
nc_close(handle)



###############################################################################
# Quantify the global homogenization uncertainty
###############################################################################

# take the global mean of each (long running!)

homogGlobalMean <- array(NA, dim=c(nt, nens))
homogNHemMean   <- array(NA, dim=c(nt, nens))
homogSHemMean   <- array(NA, dim=c(nt, nens))
homogBandMean   <- array(NA, dim=c(nt, 8, nens))	


pb   <- txtProgressBar(0, nens, style=3)

for(i in 1:nens){
	setTxtProgressBar(pb, i)

	# load step
	handle <- nc_open(filePaths[i])
	anomArray <- ncvar_get(handle, 'tempanomaly')
	nc_close(handle)

	anomArrayMasked <- maskEnsembleData(anomArray, timeMap, possibleDataMask, tDec)

	# mean step
	meanTemp <- globalMean(anomArrayMasked,mask=landMaskList$maximalMask,lat,
		nCores=nCores, zoneMask,simpleMean=FALSE, ALn, ALs)

	# package results step
	homogGlobalMean[,i] <- meanTemp$global
	homogNHemMean[,i]   <- meanTemp$nh
	homogSHemMean[,i]   <- meanTemp$sh
	homogBandMean[,,i]  <- meanTemp$bands

	# run a garbage collect in case
	rm(anomArray)
	gc()
}


save(homogGlobalMean,homogNHemMean,homogSHemMean,homogBandMean,
	file=sprintf('%s/Output/LsatAnalysis/homogMeans.Rda',scratchDir))


###############################################################################
# Quantify the global mean 
###############################################################################
nSampling <- 200

# load the sampling variance estiamtes
load(sprintf('%s/Intermediate/CovarianceMats/spatialAnalysisFullData_Empirical.Rda',scratchDir))

# some controls for the analysis
nDraws <- nSampling*12*10
timeInds <- 1:(12*10)


###############################################################################
# Re-generate the sampling ensemble and calculate means
###############################################################################


samplingGlobalMean <- array(NA, dim=c(nDec, 12*10, nSampling))
samplingNHemMean   <- array(NA, dim=c(nDec, 12*10, nSampling))
samplingSHemMean   <- array(NA, dim=c(nDec, 12*10, nSampling))
samplingBandMean   <- array(NA, dim=c(nDec, 12*10, 8, nSampling))
  
# build a full grid list
fullGridInds <- expand.grid(1:nlon, 1:nlat)

for(d in 1:nDec){

	# pull the correct dataInds
	dataInds <- dataIndsList[[d]]

	# get the right positions of a vector for the grid list
	gridListInds <- c()
	for(i in 1:nrow(dataInds)){
		gridListInds[i] <- which(dataInds[i,1] == fullGridInds[,1] & dataInds[i,2] == fullGridInds[,2])
	}

	# get the timepoints from the ERA analysis to draw from
	diffMatInds <- sample(1:ncol(diffMatList[[d]]), nDraws, replace = TRUE)

	# construct a sample of the error matricies to use in the ensemble
	rawSample <- diffMatList[[d]][,diffMatInds]

	# empty single time slice to help rearrange data
	
	sampleArray <- array(NA, dim=c(nlon,nlat,12*10,nSampling))

	for(j in 1:nSampling){
		
		i <- 1

		while(i <= 12*10){
			# generate a random chunk size
			chunkSize <- sample(1:18,1)

			# empty single time slice to help rearrange data
			sampleSlice <- matrix(NA, nlon,nlat)

			# rearrange and save into the out array
			sampleSlice[gridListInds] <- rawSample[,(((i-1)*nSampling) + j)]

			tempYearInds <- timeInds[i:(i+chunkSize-1)]


			if((i+chunkSize) >= length(timeInds)){
				tempYearInds <- timeInds[i:length(timeInds)]
				i <- length(timeInds) + 1
			} else{
				i <- i + chunkSize
			}

			# loop over the months that this particular sampling draw applies to
			# and add in the sampling unceratinty
			for(k in 1:length(tempYearInds)){
				sampleArray[,,tempYearInds[k],j] <- sampleSlice
			}	

		}
	}


	# mask the data prior to taking the mean
	sampleArrayMasked <- array(NA, dim=dim(sampleArray))
	monthVec <- rep(1:12, 10)

	for(i in 1:nSampling){
		for(j in 1:dim(sampleArray)[3])
			sampleArrayMasked[,,j,i] <- sampleArray[,,j,i] * possibleDataMask[,,monthVec[j],d]
	}

	# take the mean
	for(i in 1:nSampling){
		samplingMean <- globalMean(sampleArrayMasked[,,,i],mask=landMaskList$maximalMask,lat,
			nCores=nCores, zoneMask,simpleMean=FALSE, ALn, ALs)

		# packages results
		samplingGlobalMean[d,,i] <- samplingMean$global
		samplingNHemMean[d,,i]   <- samplingMean$nh
		samplingSHemMean[d,,i]   <- samplingMean$sh
		samplingBandMean[d,,,i]  <- samplingMean$bands
	}	
	# run a garbage collect in case

	gc()
}

save(samplingGlobalMean,samplingNHemMean,samplingSHemMean,samplingBandMean,
	file=sprintf('%s/Output/LsatAnalysis/samplingMeans.Rda',scratchDir))


