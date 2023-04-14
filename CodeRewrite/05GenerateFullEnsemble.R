# to pull the ensemble data off of discover, use:
# rsync -R nlenssen@discover.nccs.nasa.gov:/discover/nobackup/mjhendri/gistemp_ensemble/data/run_ensemble_icemask_100/FLs\*/\*.nc .
# mv ./discover/nobackup/mjhendri/gistemp_ensemble/data/run_ensemble_icemask_100/* .
# rm -r ./discover/

# if the csv mean series files are needed, use:
# rsync -R nlenssen@discover.nccs.nasa.gov:/discover/nobackup/mjhendri/gistemp_ensemble/data/run_ensemble_icemask_100/FLs\*/\*.csv .
# mv ./discover/nobackup/mjhendri/gistemp_ensemble/data/run_ensemble_icemask_100/* .
# rm -r ./discover/


# remove once running for reals
source('Namelists/awsNew.Rnl')


###############################################################################
# Load in some data
###############################################################################
# load in the land mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

# load the sampling variance estiamtes
load(sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis.Rda',scratchDir))

###############################################################################
# get the ensembles
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
				chunkSize <- sample(1:18,1)

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
	timeDim     <- ncdim_def("time",'Year',as.integer(time))
	recordDim   <- ncdim_def("record",'ensemble_member',as.integer(memberInds),unlim=TRUE)

	# define variables
	fillvalue <- -9999
	dlname <- "LSAT Anomaly Uncertainty Ensemble"
	temp_def <- ncvar_def("tempAnom",'Kelvin',list(londim,latdim,timeDim,recordDim),
							fillvalue,dlname, prec="single")

	# create netCDF file and put arrays
	ncfname <- sprintf('%s/ensembleChunk_%03d.nc',ensembleOutDir,en)

	# remove exisitng netcdf 
	if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

	ncout <- nc_create(ncfname,list(temp_def),force_v4=TRUE)

	# put variables
	ncvar_put(ncout,temp_def,samplingEnsemble)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"lon","axis","X")
	ncatt_put(ncout,"lat","axis","Y")
	ncatt_put(ncout,"time","axis","T")


	# add global attributes
	ncatt_put(ncout,0,"title",'UncertaintyEnsemble')
	ncatt_put(ncout,0,"institution",'NASA GISTEMP')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# CRITICAL: close the netcdf so it is readable
	nc_close(ncout)

	# save(samplingEnsemble,file=sprintf('Data/Ensemble/ensembleChunk_%03d.Rda',en))
}

# for(en in 1:length(filePaths)){
# 	singleEnsemble(en)
# }


cl <- makeCluster(nCores_Step5)
registerDoParallel(cl)

foreach(en=1:length(filePaths)) %dopar% singleEnsemble(en)

stopCluster(cl)



# concat all the netcdf files and remove the chunks (not sure if this is how I want
# to proceed right now)

# currentDir <- getwd()

# setwd(sprintf('cd %s',ensembleOutDir))
# system('ncrcat ensembleChunk_*.nc -O ensembleFull.nc')
# system('rm ensembleChunk_*')
# setwd(sprintf('cd %s',currentDir))

