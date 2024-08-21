# pull step03 analysis for the sampling unc
load(sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis.Rda',scratchDir))

# pull the raw ERA data and key grid info
load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/GistempProduction/coverageInfo.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/possibleDataMask.Rda',scratchDir))

nCores <- 10
nlon <- length(lon)
nlat <- length(lat)

nt <- nrow(timeMap)

###############################################################################
###############################################################################
# Quantify the field LSAT homogenization uncertainty
###############################################################################
###############################################################################


###############################################################################
# get the important metadata before running 
###############################################################################

filePaths <- system(sprintf('ls %s/Intermediate/GistempLandEnsemble/F*/*.nc',scratchDir),intern=T)[lsatInds] 
nens <- length(filePaths)

# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

anomArray <- ncvar_get(handle, 'tempanomaly', start=c(1,1,1), count =c(-1, -1, -1))

nlon <- length(lon)
nlat <- length(lat)

nt         <- handle$dim$time$len
# nSamples   <- handle$dim$record$len
nc_close(handle)

# make a time map
timeMap <- cbind(rep(startYear:endYear,each=12),1:12,NA)
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12
timeMap <- timeMap[timeMap[,3] <= endYear + (endMonth-1)/12,]

if(nt!=nrow(timeMap)) warning('Potential Time Dimension Mismatch')

# get the Aln and Als values in eary to reach variables
ALn <- landMaskList$landAreaList$ALn
ALs <- landMaskList$landAreaList$ALs


###############################################################################
# run the homog analysis
###############################################################################

# How we are chunking up the lon/lat grid
lonChunkSize <- 9
latChunkSize <- 9

lonStarts <- seq(1,nlon,by=lonChunkSize)
latStarts <- seq(1,nlat,by=latChunkSize)

quantProbs <- c(0.025,0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)


chunkStatistics <- function(i){
	require(ncdf4)

	anomMat <- array(NA, dim=c(lonChunkSize, latChunkSize, nt, nens))

	monthlyEnsembleSize      <- list()
	monthlyEnsembleMean      <- list()
	monthlyEnsembleSd        <- list()
	monthlyEnsembleQuantiles <- list()

	for(j in 1:length(latStarts)){

		# load the anomMat with all data in the lon/lat range
		for(k in 1:nens){
			lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
			latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)			

			handle <- nc_open(filePaths[k])
			anomMat[,,,k] <- ncvar_get(handle, 'tempanomaly',
													start=c(lonStarts[i],latStarts[j],1),
													count =c(lonChunkSize, latChunkSize, -1))
			nc_close(handle)	
		}

		# any analysis here on the full lon/lat chunk
		monthlyEnsembleSize[[j]] <- apply(anomMat,c(1,2,3),function(x) sum(!is.na(x)))
		monthlyEnsembleMean[[j]] <- apply(anomMat,c(1,2,3),mean, na.rm=T)
		monthlyEnsembleSd[[j]]   <- apply(anomMat,c(1,2,3),sd, na.rm=T)

		monthlyEnsembleQuantiles[[j]] <- apply(anomMat,c(1,2,3),quantile,
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

homogEnsembleSize <- array(NA, dim=c(nlon, nlat, nt))
homogEnsembleMean <- array(NA, dim=c(nlon, nlat, nt))
homogEnsembleSd   <- array(NA, dim=c(nlon, nlat, nt))

homogEnsembleQuantiles <- array(NA, dim=c(nlon, nlat, nt, length(quantProbs)))


for(i in 1:length(lonStarts)){
	for(j in 1:length(latStarts)){

		lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
		latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)



		# any analysis here on the full lon/lat chunk
		homogEnsembleSize[lonInds,latInds,] <- outList[[i]][[1]][[j]]
		homogEnsembleMean[lonInds,latInds,] <- outList[[i]][[2]][[j]]
		homogEnsembleSd[lonInds,latInds,]   <- outList[[i]][[3]][[j]]

		quantileTemp <- outList[[i]][[4]][[j]]

		for(q in 1:length(quantProbs)){
			homogEnsembleQuantiles[lonInds,latInds,,q] <- quantileTemp[q,,,]
		}

	}
}

# save the final
save(homogEnsembleSize, homogEnsembleMean, homogEnsembleSd, homogEnsembleQuantiles,
	file= sprintf('%s/Output/LsatAnalysis/homogGriddedLsatAnalysis.Rda',scratchDir))


###############################################################################
###############################################################################
# Quantify the field LSAT Sampling uncertainty
###############################################################################
###############################################################################


###############################################################################
# get the important metadata before running anymore
###############################################################################

filePaths <- system(sprintf('ls %s/Intermediate/SamplingUncertainty/SamplingEnsemble/*.nc',scratchDir),intern=T)
nens <- length(filePaths)

# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

anomArray <- ncvar_get(handle, 'tempAnom', start=c(1,1,1,1), count =c(-1, -1, -1, 1))

nlon <- length(lon)
nlat <- length(lat)

nt         <- handle$dim$time$len
nSamples   <- handle$dim$record$len
nc_close(handle)

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

	monthlyEnsembleSize      <- list()
	monthlyEnsembleMean      <- list()
	monthlyEnsembleSd        <- list()
	monthlyEnsembleQuantiles <- list()

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
		monthlyEnsembleSize[[j]] <- apply(anomMat,c(1,2,3),function(x) sum(!is.na(x)))
		monthlyEnsembleMean[[j]] <- apply(anomMat,c(1,2,3),mean, na.rm=T)
		monthlyEnsembleSd[[j]]   <- apply(anomMat,c(1,2,3),sd, na.rm=T)

		monthlyEnsembleQuantiles[[j]] <- apply(anomMat,c(1,2,3),quantile,
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

samplingEnsembleSize <- array(NA, dim=c(nlon, nlat, nt))
samplingEnsembleMean <- array(NA, dim=c(nlon, nlat, nt))
samplingEnsembleSd   <- array(NA, dim=c(nlon, nlat, nt))

samplingEnsembleQuantiles <- array(NA, dim=c(nlon, nlat, nt, length(quantProbs)))


for(i in 1:length(lonStarts)){
	for(j in 1:length(latStarts)){

		lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
		latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)



		# any analysis here on the full lon/lat chunk
		samplingEnsembleSize[lonInds,latInds,] <- outList[[i]][[1]][[j]]
		samplingEnsembleMean[lonInds,latInds,] <- outList[[i]][[2]][[j]]
		samplingEnsembleSd[lonInds,latInds,]   <- outList[[i]][[3]][[j]]

		quantileTemp <- outList[[i]][[4]][[j]]

		for(q in 1:length(quantProbs)){
			samplingEnsembleQuantiles[lonInds,latInds,,q] <- quantileTemp[q,,,]
		}

	}
}



# save the final
save(samplingEnsembleSize, samplingEnsembleMean, samplingEnsembleSd, samplingEnsembleQuantiles,
	file= sprintf('%s/Output/LsatAnalysis/samplingGriddedLsatAnalysis.Rda',scratchDir))


