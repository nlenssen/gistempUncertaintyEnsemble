# pull step03 analysis for the sampling unc
load(sprintf('%s/Intermediate/SamplingUncertainty/samplingUncertaintyAnalysis.Rda',scratchDir))

# pull the raw ERA data and key grid info
load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/GistempProduction/coverageInfo.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/possibleDataMask.Rda',scratchDir))

nlon <- length(lon)
nlat <- length(lat)

nt <- nrow(timeMap)


# overwrite nSampling if we want to get a better sample size of this analysis
# as we don't have the extra 100 replicates from the LSAT analysis
nSampling <- 200

blocksize <- 10

for(bk in 1:(nSampling/blocksize)){

	samplingEnsembleSim <- array(NA, dim=c(nlon,nlat,nt,blocksize))

	# build a full grid list
	for(d in 1:length(decInds)){

		# get the inds of the decade we are working on while also
		# accounting for the last decade not being complete here
		timeIndsSub <- ((d-1)*120 + 1):min((d*120),nt)

		if(d==14){
			timeIndsSub <- ((d-1)*120 + 1):min((d*120)+12,nt)
		}

		# rearrange and add to the homog uncertainty array
		for(j in 1:blocksize){
			
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

				# get the full "true" anomaly fields from ERA5
				fullAnomData <- anomalyData$anomalyField[,,diffArrayInds]

				# make it a 3-array if only length of one
				if(length(tempYearInds) == 1){
					samplingUncertaintyArray <- array(samplingUncertaintyArray, dim=c(nlon,nlat,1))
					fullAnomData             <- array(fullAnomData, dim=c(nlon,nlat,1))
				}

				# merge the "true" data with the estimated data for this time point and 
				# save to the full output array
				# Note: currently doing in a crude way where every non-estimated point,
				# regardless of surface type, is filled in with ERA5 data as there
				# is a mask step in the mean calculation analysis
				infillMask <- is.na(samplingUncertaintyArray[,,1])
				for(k in 1:length(tempYearInds)){
					tempMat <- samplingUncertaintyArray[,,k]
					tempInfill <- fullAnomData[,,k]
					tempMat[infillMask] <- tempInfill[infillMask]

					samplingEnsembleSim[,,tempYearInds[k],j] <- tempMat
				}
			}
		}
	}

	##############################
	# Save step
	##############################
	memberInds <- 1:blocksize + (bk-1)*blocksize
	time <- timeMap[,3]

	# define dimensions
	londim      <- ncdim_def("lon","degrees_east",as.integer(lon)) 
	latdim      <- ncdim_def("lat","degrees_north",as.integer(lat)) 
	timeDim     <- ncdim_def("time",'Year',as.integer(time))
	recordDim   <- ncdim_def("record",'ensemble_member',as.integer(memberInds),unlim=TRUE)

	# define variables
	fillvalue <- -9999
	dlname <- "LSAT *SAMPLING* Uncertainty Ensemble"
	temp_def <- ncvar_def("tempAnom",'Kelvin',list(londim,latdim,timeDim,recordDim),
							fillvalue,dlname, prec="single")

	# create netCDF file and put arrays
	ncfname <- sprintf('%s/Intermediate/SamplingUncertainty/SamplingEnsemble/member_%02d.nc',scratchDir,bk)

	# remove exisitng netcdf 
	if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

	ncout <- nc_create(ncfname,list(temp_def),force_v4=TRUE)

	# put variables
	ncvar_put(ncout,temp_def,samplingEnsembleSim)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"lon","axis","X")
	ncatt_put(ncout,"lat","axis","Y")
	ncatt_put(ncout,"time","axis","T")

	# add global attributes
	ncatt_put(ncout,0,"title",'*SAMPLING* Uncertainty Ensemble')
	ncatt_put(ncout,0,"institution",'NASA GISTEMP')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# CRITICAL: close the netcdf so it is readable
	nc_close(ncout)

}
