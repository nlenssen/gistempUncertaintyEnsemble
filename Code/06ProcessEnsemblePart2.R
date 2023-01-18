# A flag to restart the analysis from the middle if debugging
restart <- FALSE

###############################################################################
# get the important metadata before running anymore
###############################################################################
# load in the land and zone mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
zoneMask <- anomalyData$zoneMask
rm(anomalyData)

# get the gistemp metadata
filePaths <- system(sprintf('ls %s/*.nc',ensembleOutDir),intern=T)
nens <- length(filePaths)

# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

anomArray <- ncvar_get(handle, 'tempAnom', start=c(1,1,1,1), count =c(-1, -1, -1, 1))

nlon <- length(lon)
nlat <- length(lat)

nt         <- handle$dim$time$len

# NOTE: nsamples set in the namelist

nc_close(handle)

# make a coverage mask
countMat <- apply(anomArray, c(1,2), function(x) sum(!is.na(x)))

# make a time map
timeMap <- cbind(rep(startYear:endYear,each=12),1:12,NA)
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12
timeMap <- timeMap[timeMap[,3] <= endYear + (endMonth-1)/12,]

if(nt!=nrow(timeMap)) warning('Potential Time Dimension Mismatch')

# get the Aln and Als values in eary to reach variables
ALn <- landMaskList$landAreaList$ALn
ALs <- landMaskList$landAreaList$ALs


###############################################################################
# calculate the global and band means
###############################################################################

if(restart){
	load('Data/ProcessedEnsemble/meanSeries.Rda')
} else{
	ensembleGlobalMean <- array(NA, dim=c(nt, nens, nSamples))
	ensembleNHemMean   <- array(NA, dim=c(nt, nens, nSamples))
	ensembleSHemMean   <- array(NA, dim=c(nt, nens, nSamples))
	ensembleBandMean   <- array(NA, dim=c(nt, 8, nens, nSamples))	
}




for(j in 1:nSamples){
	print(sprintf('Starting Major Loop %s of %s', j, nSamples))
	pb   <- txtProgressBar(0, 100, style=3)

	for(i in 1:nens){
		setTxtProgressBar(pb, i)
		# load step
		handle <- nc_open(filePaths[i])
		anomArray <- ncvar_get(handle, 'tempAnom', start=c(1,1,1,j), count =c(-1, -1, -1, 1))
		nc_close(handle)

		# mean step
		meanTemp <- globalMean(anomArray,mask=NULL,lat,nCores=nCores,zoneMask,simpleMean=FALSE, ALn, ALs)

		# package results step
		ensembleGlobalMean[,i,j] <- meanTemp$global
		ensembleNHemMean[,i,j]   <- meanTemp$nh
		ensembleSHemMean[,i,j]   <- meanTemp$sh
		ensembleBandMean[,,i,j]   <- meanTemp$bands

		# run a garbage collect in case
		rm(anomArray)
		gc()
	}
}

save(ensembleGlobalMean, ensembleNHemMean , ensembleSHemMean, ensembleBandMean,
	file= sprintf('%s/Output/EnsembleStats/meanSeries.Rda',scratchDir))
