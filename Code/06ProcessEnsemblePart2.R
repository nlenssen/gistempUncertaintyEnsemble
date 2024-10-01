###############################################################################
###############################################################################
# Calculate the key time-seres statistics from the ensemble

# GISTEMP Uncertainty Ensemble
# Version 1.0.0 (August 21, 2024)
# Nathan Lenssen (lenssen@mines.edu)
# https://data.giss.nasa.gov/gistemp/
###############################################################################
###############################################################################


# calculate land-only and ocean-only means as well as the combined means
landMeans  <- TRUE
oceanMeans <- TRUE

###############################################################################
# get the important metadata before running anymore
###############################################################################
# load in the land and zone mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
zoneMask <- anomalyData$zoneMask
rm(anomalyData)

# get the gistemp metadata
filePaths <- system(sprintf('ls %s/ensembleChunk*',ensembleOutDir),intern=T)
nens <- length(filePaths)

# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

anomArray <- ncvar_get(handle, 'tas')

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

chunkTimeMeans <- function(k, mask=NULL){
	require(ncdf4)

	handle <- nc_open(filePaths[k])
	anomArray <- ncvar_get(handle, 'tas')
	nc_close(handle)

	# mean step (weights are 1/2,1/2 as we have full coverage here (SST+LSAT))
	meanTemp <- globalMean(anomArray, mask=mask, lat, nCores=1,
		zoneMask, simpleMean=FALSE, 1/2, 1/2)

	rm(anomalyData)
	gc()

	return(meanTemp)
}


# run the analysis (Total Global)
cl <- makeCluster(nCores_Step6b)
registerDoParallel(cl)

outList <- foreach(k=1:nens) %dopar% chunkTimeMeans(k)

stopCluster(cl)


# Process the output
ensembleGlobalMean <- array(NA, dim=c(nt, nens))
ensembleNHemMean   <- array(NA, dim=c(nt, nens))
ensembleSHemMean   <- array(NA, dim=c(nt, nens))

ensembleBandMean   <- array(NA, dim=c(nt, 8, nens))	


for(k in 1:length(outList)){

		# pull a single results list
		tempList <- outList[[k]]

		# package results step
		ensembleGlobalMean[,k] <- tempList$global
		ensembleNHemMean[,k]   <- tempList$nh
		ensembleSHemMean[,k]   <- tempList$sh
		ensembleBandMean[,,k]   <- tempList$bands

}

# save as Rda
save(ensembleGlobalMean, ensembleNHemMean , ensembleSHemMean, ensembleBandMean,
	file= sprintf('%s/meanSeries_NEW.Rda',ensembleOutDir))

##############################
# Save as netcdf
##############################

fillvalue <- -9999

# define dimensions
ensDim		<- ncdim_def("ens",'Ensemble Member',as.integer(1:nens))
bandDim		<- ncdim_def("band",'Zonal Band',as.integer(1:8))	

# Define the time dimension ~properly~ now to account for leap years
dateSeq <- seq(as.Date('1880-01-15'), as.Date(sprintf('%s-12-15',endYear)), by='month')
timeVecOut <- dateSeq - as.Date('1880-01-01')

timeDim     <- ncdim_def("time",'days since 1880-01-01 00:00:00',as.integer(timeVecOut))
		  
# define loop stuff
nameVec <- sprintf("Ensemble SST+LSAT Mean Series (%s)",c('Global Mean', 'NH Mean', 'SH Mean', 'Zonal Band Mean'))

ofNameVec <- c('Global', 'NH', 'SH', 'Band')


objectVector <- c('ensembleGlobalMean', 'ensembleNHemMean',
				  'ensembleSHemMean', 'ensembleBandMean')

for(i in 1:length(nameVec)){
	dlname <- nameVec[i]

	if(i < 4){
		field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,ensDim),
					fillvalue,dlname, prec="single")		
	} else{
		field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,bandDim, ensDim),
					fillvalue,dlname, prec="single")	
	}


	# create netCDF file and put arrays
	ncfname <- sprintf('%s/ensembleCombinedSeries_%s.nc',ensembleOutDir,ofNameVec[i])

	# remove exisitng netcdf 
	if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

	ncout <- nc_create(ncfname, list(field_def),force_v4=TRUE)

	# put variables
	ncvar_put(ncout,field_def,get(objectVector[i]))

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"time","axis","T")

	# add global attributes
	ncatt_put(ncout,0,"title",'Uncertainty Ensemble')
	ncatt_put(ncout,0,"institution",'NASA GISTEMP')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# CRITICAL: close the netcdf so it is readable
	nc_close(ncout)
}


###############################################################################
# calculate the LSAT ONLY global and band means
###############################################################################

if(landMeans){
	# create the land mask (currently not accounting for monthly sea ice)
	landMask <- landMaskList$maximalMask


	# run the analysis (Land Only)
	cl <- makeCluster(nCores_Step6b)
	registerDoParallel(cl)

	outList <- foreach(k=1:nens) %dopar% chunkTimeMeans(k, mask=landMask)

	stopCluster(cl)


	# Process the output
	ensembleGlobalMean <- array(NA, dim=c(nt, nens))
	ensembleNHemMean   <- array(NA, dim=c(nt, nens))
	ensembleSHemMean   <- array(NA, dim=c(nt, nens))

	ensembleBandMean   <- array(NA, dim=c(nt, 8, nens))	


	for(k in 1:length(outList)){

			# pull a single results list
			tempList <- outList[[k]]

			# package results step
			ensembleGlobalMean[,k] <- tempList$global
			ensembleNHemMean[,k]   <- tempList$nh
			ensembleSHemMean[,k]   <- tempList$sh
			ensembleBandMean[,,k]   <- tempList$bands

	}

	save(ensembleGlobalMean, ensembleNHemMean , ensembleSHemMean, ensembleBandMean,
		file= sprintf('%s/meanSeries_NEW_Land.Rda',ensembleOutDir))

	##############################
	# Save as netcdf
	##############################

	fillvalue <- -9999

	# define dimensions
	ensDim		<- ncdim_def("ens",'Ensemble Member',as.integer(1:nens))
	bandDim		<- ncdim_def("band",'Zonal Band',as.integer(1:8))

	# Define the time dimension ~properly~ now to account for leap years
	dateSeq <- seq(as.Date('1880-01-15'), as.Date(sprintf('%s-12-15',endYear)), by='month')
	timeVecOut <- dateSeq - as.Date('1880-01-01')

	timeDim     <- ncdim_def("time",'days since 1880-01-01 00:00:00',as.integer(timeVecOut))

	# define loop stuff
	nameVec <- sprintf("Ensemble LSAT Mean Series (%s)",c('Global Mean', 'NH Mean', 'SH Mean', 'Zonal Band Mean'))

	ofNameVec <- c('Global', 'NH', 'SH', 'Band')


	objectVector <- c('ensembleGlobalMean', 'ensembleNHemMean',
					  'ensembleSHemMean', 'ensembleBandMean')

	for(i in 1:length(nameVec)){
		dlname <- nameVec[i]

		if(i < 4){
			field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,ensDim),
						fillvalue,dlname, prec="single")		
		} else{
			field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,bandDim, ensDim),
						fillvalue,dlname, prec="single")	
		}


		# create netCDF file and put arrays
		ncfname <- sprintf('%s/ensembleLSATSeries_%s.nc',ensembleOutDir,ofNameVec[i])

		# remove exisitng netcdf 
		if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

		ncout <- nc_create(ncfname, list(field_def),force_v4=TRUE)

		# put variables
		ncvar_put(ncout,field_def,get(objectVector[i]))

		# put additional attributes into dimension and data variables
		ncatt_put(ncout,"time","axis","T")

		# add global attributes
		ncatt_put(ncout,0,"title",'Uncertainty Ensemble')
		ncatt_put(ncout,0,"institution",'NASA GISTEMP')
		history <- paste("N. Lenssen", date(), sep=", ")
		ncatt_put(ncout,0,"history",history)

		# CRITICAL: close the netcdf so it is readable
		nc_close(ncout)
	}
}



###############################################################################
# calculate the SST ONLY global and band means
###############################################################################

if(oceanMeans){
	# create the land mask (currently not accounting for monthly sea ice)
	oceanMask <- ifelse(is.na(landMaskList$maximalMask),1,NA)


	# run the analysis (Land Only)
	cl <- makeCluster(nCores_Step6b)
	registerDoParallel(cl)

	outList <- foreach(k=1:nens) %dopar% chunkTimeMeans(k, mask=oceanMask)

	stopCluster(cl)


	# Process the output
	ensembleGlobalMean <- array(NA, dim=c(nt, nens))
	ensembleNHemMean   <- array(NA, dim=c(nt, nens))
	ensembleSHemMean   <- array(NA, dim=c(nt, nens))

	ensembleBandMean   <- array(NA, dim=c(nt, 8, nens))	


	for(k in 1:length(outList)){

			# pull a single results list
			tempList <- outList[[k]]

			# package results step
			ensembleGlobalMean[,k] <- tempList$global
			ensembleNHemMean[,k]   <- tempList$nh
			ensembleSHemMean[,k]   <- tempList$sh
			ensembleBandMean[,,k]   <- tempList$bands

	}

	save(ensembleGlobalMean, ensembleNHemMean , ensembleSHemMean, ensembleBandMean,
		file= sprintf('%s/meanSeries_NEW_Ocean.Rda',ensembleOutDir))

	##############################
	# Save as netcdf
	##############################

	fillvalue <- -9999

	# define dimensions
	ensDim		<- ncdim_def("ens",'Ensemble Member',as.integer(1:nens))
	bandDim		<- ncdim_def("band",'Zonal Band',as.integer(1:8))	

	# Define the time dimension ~properly~ now to account for leap years
	dateSeq <- seq(as.Date('1880-01-15'), as.Date(sprintf('%s-12-15',endYear)), by='month')
	timeVecOut <- dateSeq - as.Date('1880-01-01')

	timeDim     <- ncdim_def("time",'days since 1880-01-01 00:00:00',as.integer(timeVecOut))
	
	# define loop stuff
	nameVec <- sprintf("Ensemble SST Mean Series (%s)",c('Global Mean', 'NH Mean', 'SH Mean', 'Zonal Band Mean'))

	ofNameVec <- c('Global', 'NH', 'SH', 'Band')


	objectVector <- c('ensembleGlobalMean', 'ensembleNHemMean',
					  'ensembleSHemMean', 'ensembleBandMean')

	for(i in 1:length(nameVec)){
		dlname <- nameVec[i]

		if(i < 4){
			field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,ensDim),
						fillvalue,dlname, prec="single")		
		} else{
			field_def <- ncvar_def('tas', 'Kelvin',list(timeDim,bandDim, ensDim),
						fillvalue,dlname, prec="single")	
		}


		# create netCDF file and put arrays
		ncfname <- sprintf('%s/ensembleSSTSeries_%s.nc',ensembleOutDir,ofNameVec[i])

		# remove exisitng netcdf 
		if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

		ncout <- nc_create(ncfname, list(field_def),force_v4=TRUE)

		# put variables
		ncvar_put(ncout,field_def,get(objectVector[i]))

		# put additional attributes into dimension and data variables
		ncatt_put(ncout,"time","axis","T")

		# add global attributes
		ncatt_put(ncout,0,"title",'Uncertainty Ensemble')
		ncatt_put(ncout,0,"institution",'NASA GISTEMP')
		history <- paste("N. Lenssen", date(), sep=", ")
		ncatt_put(ncout,0,"history",history)

		# CRITICAL: close the netcdf so it is readable
		nc_close(ncout)
	}
}

