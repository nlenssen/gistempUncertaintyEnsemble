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
# Quantify the global/hemi LSAT homogenization uncertainty
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
nt   <- handle$dim$time$len
nc_close(handle)


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
# 02) Quantify the global/hemi LSAT sampling uncertainty
###############################################################################

filePaths <- system(sprintf('ls %s/Intermediate/SamplingUncertainty/SamplingEnsemble/',scratchDir),intern=T)
nens <- length(filePaths)


# get some grid data
handle <- nc_open(sprintf('%s/Intermediate/SamplingUncertainty/SamplingEnsemble/%s',scratchDir,filePaths[1]))

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

nlon <- length(lon)
nlat <- length(lat)
nt   <- handle$dim$time$len
nmem <- handle$dim$record$len
nc_close(handle)


# take the global mean of each (long running!)

samplingGlobalMean <- array(NA, dim=c(nt, nens*nmem))
samplingNHemMean   <- array(NA, dim=c(nt, nens*nmem))
samplingSHemMean   <- array(NA, dim=c(nt, nens*nmem))
samplingBandMean   <- array(NA, dim=c(nt, 8, nens*nmem))	


pb   <- txtProgressBar(0, nens, style=3)

for(i in 1:nens){
	setTxtProgressBar(pb, i)

	# load step
	handle <- nc_open(sprintf('%s/Intermediate/SamplingUncertainty/SamplingEnsemble/%s',scratchDir,filePaths[i]))
	anomArray <- ncvar_get(handle, 'tempAnom')
	nc_close(handle)

	for(j in 1:nmem){

		# mean step
		meanTemp <- globalMean(anomArray[,,,j],mask=landMaskList$maximalMask,lat,
			nCores=nCores, zoneMask,simpleMean=FALSE, ALn, ALs)

		# package results step
		indOut <- (i-1)*nmem+j

		samplingGlobalMean[,indOut] <- meanTemp$global
		samplingNHemMean[,indOut]   <- meanTemp$nh
		samplingSHemMean[,indOut]   <- meanTemp$sh
		samplingBandMean[,,indOut]  <- meanTemp$bands
	}
	# run a garbage collect in case
	rm(anomArray)
	gc()
}


save(samplingGlobalMean,samplingNHemMean,samplingSHemMean,samplingBandMean,
	file=sprintf('%s/Output/LsatAnalysis/samplingMeans.Rda',scratchDir))



