# load in the land and zone mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
zoneMask <- anomalyData$zoneMask
rm(anomalyData)

###############################################################################
# get the important metadata before running anymore
###############################################################################

filePaths <- system(sprintf('ls %s/ensembleChunk*.nc',ensembleOutDir),intern=T)
nens <- length(filePaths)

# get some grid data
handle <- nc_open(filePaths[1])

lon    <- ncvar_get(handle, 'lon')
lat    <- ncvar_get(handle, 'lat')

anomArray <- ncvar_get(handle, 'tas')

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
			anomMat[,,,k] <- ncvar_get(handle, 'tas',
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
monthlyEnsembleSize <- array(NA, dim=c(nlon, nlat, nt))
monthlyEnsembleMean <- array(NA, dim=c(nlon, nlat, nt))
monthlyEnsembleSd   <- array(NA, dim=c(nlon, nlat, nt))

monthlyEnsembleQuantiles <- array(NA, dim=c(nlon, nlat, nt, length(quantProbs)))


for(i in 1:length(lonStarts)){
	for(j in 1:length(latStarts)){

		lonInds      <- lonStarts[i]:(lonStarts[i] + lonChunkSize - 1)
		latInds      <- latStarts[j]:(latStarts[j] + latChunkSize - 1)



		# any analysis here on the full lon/lat chunk
		monthlyEnsembleSize[lonInds,latInds,] <- outList[[i]][[1]][[j]]
		monthlyEnsembleMean[lonInds,latInds,] <- outList[[i]][[2]][[j]]
		monthlyEnsembleSd[lonInds,latInds,]   <- outList[[i]][[3]][[j]]

		quantileTemp <- outList[[i]][[4]][[j]]

		for(q in 1:length(quantProbs)){
			monthlyEnsembleQuantiles[lonInds,latInds,,q] <- quantileTemp[q,,,]
		}

	}
}


##############################
# save the final grided analysis as R files
##############################


save(monthlyEnsembleSize, monthlyEnsembleMean, monthlyEnsembleSd, monthlyEnsembleQuantiles,
	file= sprintf('%s/griddedSummaryStatistics.Rda',ensembleOutDir))




##############################
# save each final grided analysis as netcdf
##############################

fillvalue <- -9999

# define dimensions
londim      <- ncdim_def("lon","degrees_east",as.integer(lon)) 
latdim      <- ncdim_def("lat","degrees_north",as.integer(lat)) 
timeDim     <- ncdim_def("time",'Months since 12/1879',as.integer(1:nt))
quantDim    <- ncdim_def("quantile",'0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975',as.integer(1:length(quantProbs)))



# define loop stuff
nameVec <- sprintf("Ensemble Field %s",c('Sample Size', 'Mean', 'SD', 'Quantiles'))
varNameVec <- c('sampleSize', 'mean', 'sd', 'quants')
unitVec    <- c('# of obs', rep('Kelvin',3))

objectVector <- c('monthlyEnsembleSize', 'monthlyEnsembleMean',
				  'monthlyEnsembleSd', 'monthlyEnsembleQuantiles')

for(i in 1:length(nameVec)){
	dlname <- nameVec[i]

	if(i < 4){
		field_def <- ncvar_def(varNameVec[i], unitVec[i],list(londim,latdim,timeDim),
					fillvalue,dlname, prec="single")		
	} else{
		field_def <- ncvar_def(varNameVec[i], unitVec[i],list(londim,latdim,timeDim, quantDim),
			fillvalue,dlname, prec="single")
	}


	# create netCDF file and put arrays
	ncfname <- sprintf('%s/griddedSummary_%s.nc',ensembleOutDir,varNameVec[i])

	# remove exisitng netcdf 
	if(file.exists(ncfname) & overwriteFiles) file.remove(ncfname)

	ncout <- nc_create(ncfname, list(field_def),force_v4=TRUE)

	# put variables
	ncvar_put(ncout,field_def,get(objectVector[i]))

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"lon","axis","X")
	ncatt_put(ncout,"lat","axis","Y")
	ncatt_put(ncout,"time","axis","T")

	# add global attributes
	ncatt_put(ncout,0,"title",'Uncertainty Ensemble')
	ncatt_put(ncout,0,"institution",'NASA GISTEMP')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# CRITICAL: close the netcdf so it is readable
	nc_close(ncout)
}