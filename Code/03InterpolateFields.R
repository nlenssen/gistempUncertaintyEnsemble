# load in the necessary data
# anomaly and land mask data from this workflow
load(sprintf('%s/Intermediate/LandMasks/landMask_merraGrid.Rda',scratchDir))

# decadal masks from the station ensemble
load(sprintf('%s/Intermediate/GHCNv4/decadalMasks_Ensemble.Rda',scratchDir))

# shortcut variable for the era5 intermediate dir to make code more readable
eraDir <- sprintf('%s/Intermediate/ERA5/',scratchDir)

###############################################################################
###############################################################################
# Run the Interpolation Analysis
###############################################################################
###############################################################################

# unpack the list for more readable code
load(sprintf('%s/anomalyData_ERA_merraGrid.Rda',eraDir))
lon <- anomalyData$lon
nlon <- length(lon)
lat <- anomalyData$lat
nlat <- length(lat)
nt <- dim(anomalyData$anomalyField)[3]

timeMap <- anomalyData$timeMap

rm(anomalyData)

nCores <- nCores_Step3

# perform the analysis decade by decade
for(i in 1:length(tDec)){

	print(paste('Decade',i,'of',length(tDec)))

	# load the data to extract to minimize hit on memory and make the dat into a matrix
	load(sprintf('%s/anomalyData_ERA_merraGrid.Rda',eraDir))
	tempAnomaly <- apply(anomalyData$anomalyField,3L,c)
	rm(anomalyData)

	# trim the data to the land surface
	goodInds    <- which(!is.na(landMaskList$maximalMask) | decadalMasks[,,i]==1)
	goodIndsMat <- which(!is.na(landMaskList$maximalMask) | decadalMasks[,,i]==1,arr.ind=T)

	anomalyMat <- tempAnomaly[goodInds,]


	# get the matrix locations of the tables
	stationIndsMat <- which(decadalMasks[,,i]==1,arr.ind=T)
	partialMask <- c(decadalMasks[,,i])[goodInds]
	stationInds <- which(partialMask==1)

	rm(tempAnomaly)
	gc()

	# interpolate the field to the maximal land mask
	interpolatedField <- interpolateFieldPointwiseLowMem(
							radius=radius,
							lon=lon,
							lat=lat,
							goodIndsMat=goodIndsMat,
							stationInds=stationInds,
							stationIndsMat=stationIndsMat,
							anomalyMat=anomalyMat)


	# save the output each loop to allow restarts
	ofname <- sprintf('%s/InterpolatedFields/interpolatedField%02d_%s.Rda',eraDir,i,radius)

	save(interpolatedField,file=ofname)
	
	# clean up to prevent memory issues
	rm(interpolatedField)
	gc()
}


###############################################################################
###############################################################################
# Write fields to netCDF and convert to 2x2 grid
###############################################################################
###############################################################################


files <- system(sprintf('ls %s/InterpolatedFields',eraDir),intern=T)

# load in the anomaly data
load(sprintf('%s/anomalyData_ERA_merraGrid.Rda',eraDir))

# loop to create for each interpolated field
for(dec in 1:length(files)){
	load(sprintf('%s/interpolatedFields/%s',eraDir, files[dec]))

	# unpack
	lon <- anomalyData$lon
	lat <- anomalyData$lat

	time <- anomalyData$timeMap[,3]

	# define dimensions
	londim   <- ncdim_def("lon","degrees_east",as.double(lon)) 
	latdim   <- ncdim_def("lat","degrees_north",as.double(lat)) 
	timeDim   <- ncdim_def("time",'months from 01/01/1950',as.double(0:(length(time)-1)),unlim=TRUE)

	# define variables
	fillvalue <- 1e32
	dlname <- "ERA5 SAT Anomaly"
	temp_def <- ncvar_def("tempAnom",'Kelvin',list(londim,latdim,timeDim),
							fillvalue,dlname,prec="single")

	# create netCDF file and put arrays
	ncfname <- sprintf('%s/InterpolatedFieldsNetcdf/interpolatedField%02d.nc',eraDir,dec)
	ncout <- nc_create(ncfname,list(temp_def),force_v4=T)

	# put variables
	ncvar_put(ncout,temp_def,interpolatedField)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"lon","axis","X")
	ncatt_put(ncout,"lat","axis","Y")
	ncatt_put(ncout,"time","axis","T")

	# add global attributes
	ncatt_put(ncout,0,"title",'Interpolated Field')
	ncatt_put(ncout,0,"institution",'NASA GISTEMP')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# close the file, writing data to disk
	nc_close(ncout)
}

# then use nco to convert each of them to a different grid size
files <- system(sprintf('ls %s/InterpolatedFieldsNetcdf/interpolated*.nc',eraDir),intern=TRUE)

for(i in 1:length(files)){
	system(sprintf('cdo remapcon,r180x90 %s %s/InterpolatedFieldsNetcdf/interpField%02d_2x2.nc',files[i],eraDir,i))
}



