# seed for reproducibility
set.seed(12345)

###############################################################################
# Load in some data
###############################################################################
# load in the land mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

# load the sampling variance estimate
if(allLocationsLandUncertainty){
	ofname <- sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis_allLoc.Rda',scratchDir)
} else{
	ofname <- sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis.Rda',scratchDir)
}

load(ofname)

###############################################################################
# get the GHCN-ERSST ensemble paths
###############################################################################
filePaths <- system(sprintf('ls %s/Intermediate/GistempLandEnsemble/F*/*.nc',scratchDir),intern=T)[lsatInds] 
nens <- length(filePaths)

###############################################################################
# open the first one and get some important metadata
###############################################################################
handle <- nc_open(filePaths[1])

lon <- ncvar_get(handle, 'lon')
lat <- ncvar_get(handle, 'lat')
anom <- ncvar_get(handle, 'tempanomaly')

nc_close(handle)

nlon <- length(lon)
nlat <- length(lat)
nt   <- dim(anom)[3]

nYear <- ceiling(nt/12)

# Build a time map to be able to pick out needed time range
timeMap <- cbind(rep(startYear:(startYear+nYear-1),each=12),1:12,NA)
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12
timeMap <- timeMap[timeMap[,3] <= endYear + (endMonth-1)/12,]

# TEST 1: timeMap is same length as data
dim(anom)[3]==nrow(timeMap)

# build a full grid list
fullGridInds <- expand.grid(1:nlon, 1:nlat)


# figure out the time stuff for the shortened time domain
timeInds <- which(timeMap[,1] >= ensembleStartYear)
timeMapSub <- timeMap[timeInds,]
nt <- nrow(timeMapSub)


###############################################################################
# Make some mask stuff to make each decades coverage consistent for sampling
# uncertainty reasons (handle in a monthly manner)
###############################################################################

# Make a decadal mask of land locations with coverage (based on GHCN coverage)
decadalLandCoverageArray <- array(NA, dim=c(nlon,nlat,12,nDec))

for(i in 1:12){
	for(j in 1:nDec){
		decadalLandCoverageArray[,,i,j] <- ifelse(!is.na(uncEstimate[,,i,j]), 1, NA)
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
										  !is.na(decadalLandCoverageArray[,,i,d]), 1, NA)
	}
}


###############################################################################
# Master function for generating ensemble
###############################################################################

singleEnsemble <- function(en){
	require(ncdf4)

	cat(paste('Ensemble',en,'of',length(filePaths),'\n'))

	# the single GHCN ensemble output to run and then save
	# Doing inside the loop to make sure that we are starting empty each time
	samplingEnsemble <- array(NA, dim=c(nlon,nlat,nt,nSampling))

	# load and time-trim the homog ensemble member
	handle <- nc_open(filePaths[en])
	anom <- ncvar_get(handle, 'tempanomaly')
	anomSub <- anom[,,timeInds]
	rm(anom)
	nc_close(handle)

	# mask the ensemble member with the possible data mask
	anomSubMasked <- maskEnsembleData(anomSub, timeMapSub, possibleDataMask, tDec)

	# loop over each decade
	for(d in 1:length(decInds)){

		# get the inds of the decade we are working on while also
		# accounting for the last decade not being complete here
		timeIndsSub <- ((d-1)*120 + 1):min((d*120),nt)

		if(d==14){
			timeIndsSub <- ((d-1)*120 + 1):min((d*120)+12,nt)
		}

		# rearrange and add to the homog uncertainty array
		for(j in 1:nSampling){
			
			i <- 1

			while(i <= length(timeIndsSub)){
				# generate a random chunk size
				chunkSize <- sample(1:12,1)

				# get the time inds within the decade
				tempYearInds <- timeIndsSub[i:(i+chunkSize-1)]

				# get the start month and a sample
				currentMonth <- timeMap[tempYearInds[1],2]

				# get the sample ind in a way that the full time series doen't run off the 
				# end of the difference series analysis
				sampleInd <- 1e10
				while((sampleInd + chunkSize - 1) > nrow(timeMapERA)){
					sampleInd <- sample(which(timeMapERA[,2]==currentMonth), 1, replace = TRUE)
				}

				# get the inds into the sampling ensemble to insert the draw
				# check for the end of the decade and increment i
				if((i+chunkSize) >= length(timeIndsSub)){
					tempYearInds <- timeIndsSub[i:length(timeIndsSub)]
					i <- length(timeIndsSub) + 1
				} else{
					i <- i + chunkSize
				}

				# now get the time inds from the difference array
				diffArrayInds <- sampleInd:(sampleInd + length(tempYearInds) - 1)


				# rearrange and save into the out array
				samplingUncertaintyArray <- differenceArray[,,diffArrayInds,d]

				# make it a 3-array if only length of one
				if(length(tempYearInds) == 1){
					samplingUncertaintyArray <- array(samplingUncertaintyArray, dim=c(nlon,nlat,1))
				}

				# make the NA values 0 to not overwrite the sampling vals
				samplingUncertaintyArray[is.na(samplingUncertaintyArray)] <- 0

				for(k in 1:length(tempYearInds)){
					samplingEnsemble[,,tempYearInds[k],j] <- anomSubMasked[,,tempYearInds[k]] + samplingUncertaintyArray[,,k]
				}
			}
		}
	}

	##############################
	# Save step
	##############################
	memberInds <- 1:nSampling + (en-1)*nSampling
	time <- timeMapSub[,3]

	# define dimensions
	londim      <- ncdim_def("lon","degrees_east",as.integer(lon)) 
	latdim      <- ncdim_def("lat","degrees_north",as.integer(lat)) 

	dayInds <- c(15.5, 45.0, 74.5, 105.0, 135.5, 166.0,
		196.5, 227.5, 258.0, 288.5, 319.0, 349.5)
	timeVecOut <- dayInds + rep(seq(0, nYear-1), each=12)*365

	timeDim     <- ncdim_def("time",'days since 1880-01-01 00:00:00',as.integer(timeVecOut))

	# define variables
	fillvalue <- -9999
	longname <- "blended LSAT and SST temperature anomalies (1951-1980 climatology)"
	temp_def <- ncvar_def("tas",'Kelvin',list(londim,latdim,timeDim),
							fillvalue,longname, prec="single")


	# loop over members
	for(m in 1:length(memberInds)){
		# create netCDF file and put arrays
		ncfname <- sprintf('%s/ensembleChunk_%04d.nc',ensembleOutDir,memberInds[m])

		# remove exisitng netcdf 
		if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

		ncout <- nc_create(ncfname,list(temp_def),force_v4=TRUE)

		# put variables
		ncvar_put(ncout,temp_def,samplingEnsemble[,,,m])

		# put additional attributes into dimension and data variables
		ncatt_put(ncout,"lon","axis","X")
		ncatt_put(ncout,"lat","axis","Y")
		ncatt_put(ncout,"time","axis","T")


		# add global attributes
		ncatt_put(ncout,0,"title",'GISTEMP v4 Uncertainty Ensemble')
		ncatt_put(ncout,0,"institution",'NASA GISTEMP')
		history <- paste("N. Lenssen", date(), sep=", ")
		ncatt_put(ncout,0,"history",history)

		# CRITICAL: close the netcdf so it is readable
		nc_close(ncout)
	}
}



cl <- makeCluster(nCores_Step5)
registerDoParallel(cl)

foreach(en=1:length(filePaths)) %dopar% singleEnsemble(en)

stopCluster(cl)