# The GHCN ensemble number used


########################################################################
# Import the GHCN station metadata into an easy to work with format
########################################################################
# Read in the (currently) relevant data from the station data

rawInfo <- readLines(sprintf('%s/Raw/GHCNv4/Operational/v4.inv',scratchDir))

# Vectorized to run over all rows of the station file

id  <- substr(rawInfo,1,11) # as string this time...
lat <- as.numeric(substr(rawInfo,13,20))
lon <- as.numeric(substr(rawInfo,22,30))


# Place in a data frame to mirror the raw, human-readable organization
stationID <- data.frame(id,lon,lat, stringsAsFactors=FALSE)

########################################################################
# Import the GHCN station data into an easy to work with format
########################################################################
# use `cat ./* > allStation.dat` to build the station file
raw <- readLines(sprintf('%s/Raw/GHCNv4/Operational/ghcnm.tavg.qcf.dat',scratchDir))

dataID <- rep(NA,length(raw))
stationData <- matrix(NA,nrow=length(raw),ncol=13)
colnames(stationData) <- c('year', 1:12)

# the function for the ghcn ensemble
parseDataString <- function(dataString){
        sapply(seq(from=1, to=nchar(dataString), by=9),
                                function(i) as.numeric(substr(dataString, i, i+4)))
}
# the function for the ghcn operational
parseDataString2 <- function(dataString){
        sapply(seq(from=1, to=nchar(dataString), by=8),
                                function(i) as.numeric(substr(dataString, i, i+4)))
}


# helper function used to parse the data portion of the input string
# 11 or 12 digit station?
pb   <- txtProgressBar(1, length(raw), style=3)

for(i in 1:length(raw)){
	setTxtProgressBar(pb, i)

	temp <- raw[i]
	station <- as.character(substr(temp,1,11))
	year    <- as.numeric(substr(temp,12,15))

	data <- parseDataString2(substr(temp,20,nchar(temp)))

	dataID[i] <- station
	stationData[i,] <- c(year,data)
}

save(stationID, stationData, dataID,  
	file = sprintf('%s/Intermediate/GHCNv4/stationData.Rda',scratchDir))

########################################################################
# Make coverage masks
########################################################################

# Load reanalysis grid metadata
load(sprintf('%s/Intermediate/LandMasks/landMask_merraGrid.Rda',scratchDir))
lon <- landMaskList$lon
lat <- landMaskList$lat
rm(landMaskList)

# Create some objects for the loop
startYears <- seq(1880,2010,by=10)

decadalMasks <- array(NA, dim=c(length(lon),
								length(lat),
								length(startYears)))
coverageList <- list()

# loop over decades to determine coverage and map to the reanalysis grid
pb   <- txtProgressBar(1, length(startYears), style=3)
for(i in 1:length(startYears)){
	setTxtProgressBar(pb, i)

	decade <- (startYears[i]-1):(startYears[i]+9)

	# start by subsetting the data that appears in this data
	partialInds <- which(stationData[,1] %in% decade)

	partial <- stationData[partialInds,]

	# replace all -9999 with NA
	partial[partial==-9999] <- NA

	# pull the unique station ids and figure out what years they 
	ids <- unique(dataID[partialInds])


	# Compute which stations have coverage in the decade
	coverageMatFull <- data.frame()

	for(j in 1:length(ids)){
		station <- subset(partial, dataID[partialInds] == ids[j])
		coverageMatFull <- rbind(coverageMatFull,
			coverageCheck(station,decade,ids[j]))
	}

	# Now take the stations with coverage and map to lon/lat
	coveredStations <- subset(stationID, stationID$id %in% 
						 coverageMatFull[coverageMatFull[,2]==1,][,1])

	coverageList[[i]] <- coveredStations

	# Now we need to fill out the merra grid based on these stations
	decadalMasks[,,i] <- gridMask(coveredStations,lon,lat)
}

save(decadalMasks,
	file = sprintf('%s/Intermediate/GHCNv4/decadalMasks_Ensemble.Rda',scratchDir))
