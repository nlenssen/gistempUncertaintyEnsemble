# load in the land and zone mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
zoneMask <- anomalyData$zoneMask
rm(anomalyData)

###############################################################################
# get the important metadata before running anymore
###############################################################################

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
# nSamples   <- handle$dim$record$len
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
# pull a region instead of just a single point to speed things up
###############################################################################

# How we are chunking up the lon/lat grid
lonChunkSize <- 9
latChunkSize <- 9

lonStarts <- seq(1,nlon,by=lonChunkSize)
latStarts <- seq(1,nlat,by=latChunkSize)

quantProbs <- c(0.025,0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)


chunkStatistics <- function(i){
	require(ncdf4)

	anomMat <- array(NA, dim=c(lonChunkSize, latChunkSize, nt, nens*nSamples))

	for(j in 1:length(latStarts)){

		# load the anomMat with all data in the lon/lat range
		for(k in 1:nens){
			lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
			latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)
			ensembleInds <- (1+(k-1)*nSamples):(k*nSamples)
			

			handle <- nc_open(filePaths[k])
			anomMat[,,,ensembleInds] <- ncvar_get(handle, 'tempAnom',
													start=c(lonStarts[i],latStarts[j],1,1),
													count =c(lonChunkSize, latChunkSize, -1, nSamples))
			nc_close(handle)	
		}

		# any analysis here on the full lon/lat chunk
		monthlyEnsembleSize<- apply(anomMat,c(1,2,3),function(x) sum(!is.na(x)))
		monthlyEnsembleMean <- apply(anomMat,c(1,2,3),mean, na.rm=T)
		monthlyEnsembleSd  <- apply(anomMat,c(1,2,3),sd, na.rm=T)

		monthlyEnsembleQuantiles <- apply(anomMat,c(1,2,3),quantile,
			probs= quantProbs, na.rm=T)

	}

	return(list(monthlyEnsembleSize, monthlyEnsembleMean, monthlyEnsembleSd, monthlyEnsembleQuantiles))
}


# run the analysis
cl <- makeCluster(nCores_Step6a)
registerDoParallel(cl)

outList <- foreach(i=1:length(lonStarts)) %dopar% chunkStatistics(i)

stopCluster(cl)

# process the output

monthlyEnsembleSize <- array(NA, dim=c(nlon, nlat, nt))
monthlyEnsembleMean <- array(NA, dim=c(nlon, nlat, nt))
monthlyEnsembleSd   <- array(NA, dim=c(nlon, nlat, nt))

monthlyEnsembleQuantiles <- array(NA, dim=c(nlon, nlat, nt, length(quantProbs)))


for(i in 1:length(lonStarts)){
	for(j in 1:length(latStarts)){

		lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
		latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)

			
		tempList <- outList[[i]]


		# any analysis here on the full lon/lat chunk
		monthlyEnsembleSize[lonInds,latInds,] <- tempList[[1]]
		monthlyEnsembleMean[lonInds,latInds,] <- tempList[[2]]
		monthlyEnsembleSd[lonInds,latInds,]   <- tempList[[3]]

		quantileTemp <- tempList[[4]]

		for(q in 1:length(quantProbs)){
			monthlyEnsembleQuantiles[lonInds,latInds,,q] <- quantileTemp[q,,,]
		}

	}
}



# save the final
save(monthlyEnsembleSize, monthlyEnsembleMean, monthlyEnsembleSd, monthlyEnsembleQuantiles,
	file= sprintf('%s/griddedSummaryStatistics.Rda',ensembleOutDir))








###############################################################################
# check SST uncetainty in a few points
###############################################################################

# lonInd <- which(lon==1)
# latInd <- which(lat==89)


# # get the matrix of time x ensemble member at a single location
# anomMat <- matrix(nrow=nt,ncol=nens*nSamples)
# for(k in 1:nens){
# 	colInds <- (1+(k-1)*nSamples):(k*nSamples)

# 	handle <- nc_open(filePaths[k])
# 	anomMat[,colInds] <- ncvar_get(handle, 'tempAnom', start=c(lonInd,latInd,1,1), count =c(1, 1, -1, -1))
# 	nc_close(handle)	
# }

# # safety check to make sure we have data in this cell
# if(all(is.na(anomMat))) next()

# # calculate statistics
# monthlyEnsembleSize[i,j,] <- apply(anomMat,1,function(x) sum(!is.na(x)))
# monthlyEnsembleMean[i,j,] <- apply(anomMat,1,mean, na.rm=T)
# monthlyEnsembleSd[i,j,]   <- apply(anomMat,1,sd, na.rm=T)

# monthlyEnsembleQuantiles[i,j,,] <- t(apply(anomMat,1,quantile,
# 	probs= quantProbs, na.rm=T))

