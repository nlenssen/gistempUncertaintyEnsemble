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



indexList <- expand.grid(list(i=1:nens, j=1:nSamples))

chunkTimeMeans <- function(k){
	require(ncdf4)


	i <- indexList[k,1]
	j <- indexList[k,2]

	handle <- nc_open(filePaths[i])
	anomArray <- ncvar_get(handle, 'tempAnom', start=c(1,1,1,j), count =c(-1, -1, -1, 1))
	nc_close(handle)

	# mean step
	meanTemp <- globalMean(anomArray, mask=NULL, lat, nCores=1,
		zoneMask, simpleMean=FALSE, ALn, ALs)

	rm(anomalyData)
	gc()

	return(meanTemp)
}


# run the analysis
cl <- makeCluster(nCores_Step6b)
registerDoParallel(cl)

outList <- foreach(k=1:nrow(indexList)) %dopar% chunkTimeMeans(k)

stopCluster(cl)



# Process the output
ensembleGlobalMean <- array(NA, dim=c(nt, nens, nSamples))
ensembleNHemMean   <- array(NA, dim=c(nt, nens, nSamples))
ensembleSHemMean   <- array(NA, dim=c(nt, nens, nSamples))

ensembleBandMean   <- array(NA, dim=c(nt, 8, nens, nSamples))	


for(k in 1:length(outList)){

		# pull a single results list
		tempList <- outList[[k]]

		# get the array inds
		i <- indexList[k,1]
		j <- indexList[k,2]

		# package results step
		ensembleGlobalMean[,i,j] <- tempList$global
		ensembleNHemMean[,i,j]   <- tempList$nh
		ensembleSHemMean[,i,j]   <- tempList$sh
		ensembleBandMean[,,i,j]   <- tempList$bands

}

save(ensembleGlobalMean, ensembleNHemMean , ensembleSHemMean, ensembleBandMean,
	file= sprintf('%s/Output/EnsembleStats/meanSeries_NEW.Rda',scratchDir))

