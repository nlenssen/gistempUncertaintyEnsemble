# read the land mask that includes sea ice first
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

lon <- landMaskList$lon
lat <- landMaskList$lat

nlon <- length(lon)
nlat <- length(lat)

landMaskWorking <- landMaskList$minimalMask
# get the total number of years in the ensemble
nt <- length(startYear:endYear)*12

# read in the full gistemp data and mask
handle <- nc_open(sprintf('%s/Raw/GistempProduction/gistemp1200_GHCNv4_ERSSTv5.nc',scratchDir))
anomalyArray <- ncvar_get(handle, 'tempanomaly', start=c(1,1,1), count = c(-1,-1,nt))
nc_close(handle)

maskedArray <- array(NA, dim=c(nlon,nlat,nt))
for(i in 1:nt){
	maskedArray[,,i] <- anomalyArray[,,i] * landMaskWorking
}

# loop through each decade and get the approx grid-cell coverage
decadalCoverageMask <- array(NA, dim=c(nlon,nlat,nDec))

for(i in 1:nDec){
	timeInds <- which(timeMap[,1] >= tDec[i] & timeMap[,1] < tDec[i+1])

	if(i ==nDec){
		timeInds <- min(which(timeMap[,1] == tDec[i])):nrow(timeMap)
	}

	tempCounts <- apply(maskedArray[,,timeInds], c(1,2), function(x) sum(!is.na(x)) )

	decadalCoverageMask[,,i] <- ifelse(tempCounts >= length(timeInds)/2 , 1, NA)
}

###############################################################################
# Get the official land mask
###############################################################################
gistempLandMaskRaw <- read.table(sprintf('%s/Raw/GistempProduction/landmask.2degx2deg.txt',scratchDir),
							skip=1, header=T)

prodLandMask <- matrix(NA, nlon, nlat)

for(i in 1:nrow(gistempLandMaskRaw)){
	prodLandMask[gistempLandMaskRaw[i,1],gistempLandMaskRaw[i,2]] <- gistempLandMaskRaw[i,5]
}

# par(mfrow = c(1,3))
# image(landMaskWorking)
# image(decadalCoverageMask[,,1])
# image.plot(prodLandMask)


###############################################################################
# save output
###############################################################################

save(landMaskWorking,decadalCoverageMask,
	file=sprintf('%s/Intermediate/GistempProduction/coverageInfo.Rda',scratchDir))

