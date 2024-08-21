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

# save the reconstruction differences as an array
differenceArray <- array(NA, dim=c(length(lon),length(lat),nt,length(decInds)))

# save the monthly uncertainty cycle for each decade
uncEstimate <- array(NA, dim=c(length(lon),length(lat),12,length(decInds)))

# loop over decade's analysis
for(d in 1:nDec){
	# Load and process the interpolated field
	interpolatedField <- read22File(files[d])$anomArray

	# mask by coverage in production gistemp
	for(i in 1:nt){
		interpolatedField[,,i] <- interpolatedField[,,i] * decadalCoverageMask[,,d]
	}

	# Loop over each location and calculate the method's error series
	for(i in 1:nlon){
		for(j in 1:nlat){

			trueVec <- anomalyField[i,j,]
			maskVec <- interpolatedField[i,j,]

			monthseq <- rep(1:12, length=length(trueVec))

			if(all(is.na(trueVec)))	next()
			if(all(is.na(maskVec))){
				if(allLocationsLandUncertainty){
					for(k in 1:12){
						uncEstimate[i,j,k,d] <- sd(trueVec[monthseq==k],na.rm=T)
					}
				} else{
					uncEstimate[i,j,,d] <- NA
				}

				next()
			}

			differenceSeries <- trueVec - maskVec

			differenceArray[i,j,,d] <- differenceSeries

			for(k in 1:12){
				monthInds <- which(timeMapERA[,2] == k)
				uncEstimate[i,j,k,d] <- sd(differenceSeries[monthInds], na.rm=T)
			}
		}
	}

}
# save the important output
if(allLocationsLandUncertainty){
	ofname <- sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis_allLoc.Rda',scratchDir)
} else{
	ofname <- sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis.Rda',scratchDir)
}

save(lon, lat, uncEstimate, differenceArray, timeMapERA, file = ofname)