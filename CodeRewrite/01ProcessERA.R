###############################################################################
###############################################################################
# (I) Convert Full ERA Fields
###############################################################################
###############################################################################

# CDO command used to map the ERA to 2x2 (perform in gistempRaw directory)
# cdo remapcon,r180x90 era5_t2m_1979_present.nc era5_t2m_1979_present_2x2.nc

# cdo remapcon,r576x361 era5_t2m_1979_present.nc era5_t2m_1979_present_merraGrid.nc
# cdo remapcon,r576x361 era5_t2m_1950_1978.nc era5_t2m_1950_1978_merraGrid.nc
# cdo remapcon,r576x361 era5_seaIce_1979_present.nc era5_seaIce_1979_present_merraGrid.nc



################################
# STEP (1) Load all the 2x2 data
################################

# get the early record first
handle <- nc_open(sprintf('%s/Raw/ERA5/era5_t2m_1950_1978_2x2.nc',scratchDir))

rawLon <- ncvar_get(handle,'lon')
nlon <- length(rawLon)
lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

rawLat <- ncvar_get(handle,'lat')
nlat <- length(rawLat)
lat <- rawLat
rawTempEarly <- ncvar_get(handle,'t2m')

nc_close(handle)

tempArrayEarly <- rawTempEarly[c((nlon/2+1):nlon,1:(nlon/2)),,]
rm(rawTempEarly)

# then get the modern record
handle <- nc_open(sprintf('%s/Raw/ERA5/era5_t2m_1979_present_2x2.nc',scratchDir))

rawTempLate <- ncvar_get(handle,'t2m',start=c(1,1,1,1),count=c(-1,-1,1,-1))

nc_close(handle)

tempArrayLate <- rawTempLate[c((nlon/2+1):nlon,1:(nlon/2)),,]
rm(rawTempLate)

# join the two together
tempArrayFull <- abind(tempArrayEarly, tempArrayLate, along=3)


# make a time map
ntFull <- dim(tempArrayEarly)[3] + dim(tempArrayLate)[3]

lastYear <- startYearERA + ceiling(ntFull/12) - 1

timeMapFull <- cbind(rep(startYearERA:lastYear,each=12),1:12,NA)[1:ntFull,]
timeMapFull[,3] <- timeMapFull[,1] + (timeMapFull[,2]-1)/12

# package the data together and clip to the desired end time point
inds <- which(timeMapFull[,1] <= endYearERA)

tempArray <- tempArrayFull[,,inds]
timeMap   <- timeMapFull[inds,]

########################################
# STEP (3) Compute Anomalies on the grid
########################################
anomalyField <- array(NA, dim=dim(tempArray))

monthNorm <- array(NA,dim=c(dim(tempArray)[1:2],12))

for(month in 1:12){

	inds <- which(timeMap[,2]==month)
	climInds <- which(timeMap[,2]==month & timeMap[,1] %in% climPeriod)

	monthNorm[,,month] <- apply(tempArray[,,climInds],c(1,2),mean)

	for(j in 1:length(inds)){	
		anomalyField[,,inds[j]] <- tempArray[,,inds[j]] - monthNorm[,,month]
	}
}

###########################
# STEP (4) Create Zone Mask
###########################

zoneMask <- matrix(NA, nrow=length(lon), ncol=length(lat))

for(i in 1:length(lon)){
	for(j in 1:length(lat)){
		zoneMask[i,j] <- zoneAssign(lon[i], lat[j])
	}
}

anomalyData <- list(anomalyField = anomalyField,
				 lon = lon,
				 lat = lat,
				 zoneMask = zoneMask,
				 timeMap = timeMap)


save(anomalyData,
	file=sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))

rm(anomalyData,anomalyField)

######################################
# STEP (5) Load all the merraGrid data
######################################

# get the early record first
handle <- nc_open(sprintf('%s/Raw/ERA5/era5_t2m_1950_1978_merraGrid.nc',scratchDir))

rawLon <- ncvar_get(handle,'lon')
nlon <- length(rawLon)
lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

rawLat <- ncvar_get(handle,'lat')
nlat <- length(rawLat)
lat <- rawLat
rawTempEarly <- ncvar_get(handle,'t2m')

nc_close(handle)

tempArrayEarly <- rawTempEarly[c((nlon/2+1):nlon,1:(nlon/2)),,]
rm(rawTempEarly)

# then get the modern record
handle <- nc_open(sprintf('%s/Raw/ERA5/era5_t2m_1979_present_merraGrid.nc',scratchDir))

rawTempLate <- ncvar_get(handle,'t2m',start=c(1,1,1,1),count=c(-1,-1,1,-1))

nc_close(handle)

tempArrayLate <- rawTempLate[c((nlon/2+1):nlon,1:(nlon/2)),,]
rm(rawTempLate)

# join the two together
tempArrayFull <- abind(tempArrayEarly, tempArrayLate, along=3)

# use the same inds as the 2x2 case
inds <- which(timeMapFull[,1] <= endYearERA)
tempArray <- tempArrayFull[,,inds]


####################################################
# STEP (6) Compute Anomalies on the grid (merraGrid)
####################################################

anomalyField <- array(NA, dim=dim(tempArray))

monthNorm <- array(NA,dim=c(dim(tempArray)[1:2],12))

for(month in 1:12){

	inds <- which(timeMap[,2]==month)
	climInds <- which(timeMap[,2]==month & timeMap[,1] %in% climPeriod)

	monthNorm[,,month] <- apply(tempArray[,,climInds],c(1,2),mean)

	for(j in 1:length(inds)){	
		anomalyField[,,inds[j]] <- tempArray[,,inds[j]] - monthNorm[,,month]
	}
}

#######################################
# STEP (7) Create Zone Mask (merraGrid)
#######################################

zoneMask <- matrix(NA, nrow=length(lon), ncol=length(lat))

for(i in 1:length(lon)){
	for(j in 1:length(lat)){
		zoneMask[i,j] <- zoneAssign(lon[i], lat[j])
	}
}

anomalyData <- list(anomalyField = anomalyField,
				 lon = lon,
				 lat = lat,
				 zoneMask = zoneMask,
				 timeMap = timeMap)


save(anomalyData,
	file=sprintf('%s/Intermediate/ERA5/anomalyData_ERA_merraGrid.Rda',scratchDir))

rm(anomalyData,anomalyField)

###############################################################################
###############################################################################
# (II) Make Land-Only Masks
###############################################################################
###############################################################################

# open the MERRA constant field to pull in grid info and make the land mask
handle <- nc_open(
	sprintf("%s/Raw/MERRA/MERRA2_101.const_2d_asm_Nx.00000000.nc4",scratchDir))

lonRef  <- ncvar_get(handle,'lon')
nlon <- length(lon)
latRef  <- ncvar_get(handle,'lat')
nlat <- length(lat)

landField  <- ncvar_get(handle,"FRLAND")
iceField   <- ncvar_get(handle,"FRLANDICE")
lakeField  <- ncvar_get(handle,"FRLAKE")

landMaskMerra <- ifelse(landField > propCutoff | iceField > propCutoff |
					    lakeField > 0, 1, NA)

nc_close(handle)


##################################################
# get the maximal monthly sea ice extent from ERA5
##################################################

# make a time map
timeMapLate <- timeMap[timeMap[,1]>1978,]

# load in the data
handle <- nc_open(sprintf('%s/Raw/ERA5/era5_seaIce_1979_present_merraGrid.nc',scratchDir))

rawLon <- ncvar_get(handle,'lon')
nlon <- length(rawLon)
lon <- c(rawLon[(nlon/2+1):nlon]-360,rawLon[1:(nlon/2)])

rawLat <- ncvar_get(handle,'lat')
nlat <- length(rawLat)
lat <- rawLat
seaIceLateRaw <- ncvar_get(handle,'siconc',
	start=c(1,1,1,1),count=c(-1,-1,1,nrow(timeMapLate)))

nc_close(handle)

seaIceLate <- ifelse(seaIceLateRaw[c((nlon/2+1):nlon,1:(nlon/2)),,] > propCutoff, 1, NA)
rm(seaIceLateRaw)

monthlySeaIceMask <- array(NA, dim=c(nlon,nlat,12))

# get the extent of sea ice that occurs in the top 25% of years
for(i in 1:12){
	subArray <- seaIceLate[,,timeMapLate[,2]==i]
	monthlySeaIceMask[,,i] <- ifelse(apply(subArray,c(1,2),sum,na.rm=T) > (dim(subArray)[3]/4),1,NA)
}

maximalSeaIce <- ifelse(apply(monthlySeaIceMask,c(1,2),sum,na.rm=T) > 0,1,NA)
minimalSeaIce <- ifelse(apply(monthlySeaIceMask,c(1,2),sum,na.rm=T) == 12,1,NA)

# make a 12-month full mask (land + sea ice)
monthlyMask <- monthlySeaIceMask

for(i in 1:12){
 	monthlyMask[,,i] <- ifelse(!is.na(monthlyMask[,,i]) | !is.na(landMaskMerra),1,NA)
}

# make a 1-month maximal full mask
maximalMask <- ifelse(!is.na(maximalSeaIce) | !is.na(landMaskMerra),1,NA)
minimalMask <- ifelse(!is.na(minimalSeaIce) | !is.na(landMaskMerra),1,NA)

#######################
# Land Area Calculation
#######################
cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))
totalArea <- sum(cosMat)

nLats <- which(lat>=0)
sLats <- which(lat<=0)

# calculate the land area for the monthly masks
ALMonthly  <- rep(NA, 12)
ALnMonthly <- rep(NA, 12)
ALsMonthly <- rep(NA, 12)

for(i in 1:12){
	ALMonthly[i]  <- sum(monthlyMask[,,i] * cosMat,na.rm=T)/totalArea
	ALnMonthly[i] <- sum(monthlyMask[,nLats,i] * cosMat[,nLats],na.rm=T)/sum(cosMat[,nLats])
	ALsMonthly[i] <- sum(monthlyMask[,sLats,i] * cosMat[,sLats],na.rm=T)/sum(cosMat[,sLats])
}

ASMonthly  <- 1 - ALMonthly
ASnMonthly <- 1 - ALnMonthly
ASsMonthly <- 1 - ALsMonthly


# take the annual values as the mean of the monthly values for when needed in calculations

AL  <- mean(ALMonthly)
ALn <- mean(ALnMonthly)
ALs <- mean(ALsMonthly)

AS  <- 1 - AL
ASn <- 1 - ALn
ASs <- 1 - ALs

landAreaList <- list(AL=AL,ALn=ALn,ALs=ALs,AS=AS,ASn=ASn,ASs=ASs,
	ALMonthly=ALMonthly,ALnMonthly=ALnMonthly,ALsMonthly=ALsMonthly,
	ASMonthly=ASMonthly,ASnMonthly=ASnMonthly,ASsMonthly=ASsMonthly)

#########################
# Save Mask and Land Area
#########################

landMaskList <- list(lon=lon,
					 lat=lat,
					 monthlyMask=monthlyMask,
					 maximalMask=maximalMask,
					 minimalMask=minimalMask,
					 landAreaList=landAreaList)

save(landMaskList,
	file=sprintf('%s/Intermediate/LandMasks/landMask_merraGrid.Rda',scratchDir))

###########################
# transform to the 2x2 grid
###########################

# NEED TO HAVE PULLED THE LAND ENSEMBLE FROM DISCOVER FIRST!
handle <- nc_open(
	sprintf('%s/Intermediate/GistempLandEnsemble/FLs_700/gistemp1200.nc',scratchDir))

lonRef <- ncvar_get(handle, 'lon')
latRef <- ncvar_get(handle, 'lat')

nc_close(handle)

# make the new mask (landMaskNewGrid() from production code)
monthlySeaIceMask22 <- array(NA, dim=c(length(lonRef),length(latRef),12))


# create the 2x2 sea ice mask
for(i in 1:length(lonRef)){
	for(j in 1:length(latRef)){
		lonInd <- which.min(abs(lonRef[i] - lon))	
		latInd <- which.min(abs(latRef[j] - lat))

		for(s in 1:12){
			monthlySeaIceMask22[i,j,s] <- monthlySeaIceMask[lonInd,latInd,s]
		}
	}
}

maximalSeaIce22 <- ifelse(apply(monthlySeaIceMask22,c(1,2),sum,na.rm=T) > 0,1,NA)
minimalSeaIce22 <- ifelse(apply(monthlySeaIceMask22,c(1,2),sum,na.rm=T) == 12,1,NA)

###################################################################
# Intersect with the production GISTEMP landmask to stay consistent
###################################################################
gistempLandMaskRaw <- read.table(sprintf('%s/Raw/GistempProduction/landmask.2degx2deg.txt',scratchDir),
							skip=1, header=T)

prodLandMask <- matrix(NA, length(lonRef), length(latRef))

for(i in 1:nrow(gistempLandMaskRaw)){
	prodLandMask[gistempLandMaskRaw[i,1],gistempLandMaskRaw[i,2]] <- gistempLandMaskRaw[i,5]
}

prodLandMaskFinal <- ifelse(prodLandMask > 0, 1, NA)


monthlyMask22 <- array(NA, dim=c(length(lonRef),length(latRef),12))
maximalMask22 <- array(NA, dim=c(length(lonRef),length(latRef)))
minimalMask22 <- array(NA, dim=c(length(lonRef),length(latRef)))

for(i in 1:length(lonRef)){
	for(j in 1:length(latRef)){
		isLand <- prodLandMaskFinal[i,j] == 1
	
		for(s in 1:12){
			monthlyMask22[i,j,s] <- ifelse(isLand | monthlySeaIceMask22[i,j,s] == 1,1,NA)
		}		
	}
}


# make a 1-month maximal full mask
maximalMask22 <- ifelse(!is.na(maximalSeaIce22) | !is.na(prodLandMaskFinal),1,NA)
minimalMask22 <- ifelse(!is.na(minimalSeaIce22) | !is.na(prodLandMaskFinal),1,NA)

#######################
# Land Area Calculation
#######################
cosMat <- cos(matrix(latRef,nrow=length(lonRef),ncol=length(latRef),byrow=T)*(pi/180))
totalArea <- sum(cosMat)

nLats <- which(latRef>=0)
sLats <- which(latRef<=0)

# calculate the land area for the monthly masks
ALMonthly  <- rep(NA, 12)
ALnMonthly <- rep(NA, 12)
ALsMonthly <- rep(NA, 12)

for(i in 1:12){
	ALMonthly[i]  <- sum(monthlyMask22[,,i] * cosMat,na.rm=T)/totalArea
	ALnMonthly[i] <- sum(monthlyMask22[,nLats,i] * cosMat[,nLats],na.rm=T)/sum(cosMat[,nLats])
	ALsMonthly[i] <- sum(monthlyMask22[,sLats,i] * cosMat[,sLats],na.rm=T)/sum(cosMat[,sLats])
}

ASMonthly  <- 1 - ALMonthly
ASnMonthly <- 1 - ALnMonthly
ASsMonthly <- 1 - ALsMonthly


# take the annual values as the mean of the monthly values for when needed in calculations

AL  <- mean(ALMonthly)
ALn <- mean(ALnMonthly)
ALs <- mean(ALsMonthly)

AS  <- 1 - AL
ASn <- 1 - ALn
ASs <- 1 - ALs

landAreaList22 <- list(AL=AL,ALn=ALn,ALs=ALs,AS=AS,ASn=ASn,ASs=ASs,
	ALMonthly=ALMonthly,ALnMonthly=ALnMonthly,ALsMonthly=ALsMonthly,
	ASMonthly=ASMonthly,ASnMonthly=ASnMonthly,ASsMonthly=ASsMonthly)

##################
# package and save
##################
landMaskList <- list(lon=lonRef,
					 lat=latRef,
					 monthlyMask=monthlyMask22,
					 maximalMask=maximalMask22,
					 minimalMask=minimalMask22,
					 landAreaList=landAreaList22,
					 monthlySeaIce = monthlySeaIceMask22)

save(landMaskList, file=sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))

