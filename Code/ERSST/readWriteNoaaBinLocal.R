# GISTEMP Uncertainty Analysis
# Version 2.0.0 (Summer 2021)
# Nathan Lenssen (n.lenssen@columbia.edu)

library(doParallel)
library(foreach)

nCores <- 2

cdir <- 'Code/ERSST'

# location of data
ddir <- '/Volumes/Archive/ERSSTv5/operationalEnsemble'

# location to write the sbbx files
ofdir <- 'SBBX_Output_2020_IceMask'

# location to write out
earlyDir <- 'ensemble.1854-2017'
lateDir  <- 'ensemble.2001.present'

earlyFilePrefix <- 'sst2d.ano.1854.2017.ensemble'
lateFilePrefix  <- 'sst2d.ano.2001.last.ensemble'

source('Code/ERSST/Functions.R')

# load in the climatology correcction calculated in correctErsstClimatology.R
correctClimatology <- TRUE
climPeriod <- 'Official'

if(correctClimatology){
	load(sprintf('Data/sstClimCorrectionArray_%s.Rda',climPeriod))
}

# Take every other ensemble member from the early ensemble
# to get an even parameter sweep of the larger ensemble
earlyIndVec <- 1:500*2
lateIndVec  <- 1:500

# The lengths of the raw data
startYear1 <- 1854
endYear1   <- 2009

startYear2 <- 2001
endYear2   <- 2020

# the year parameters to write out
startYear  <- 1854
switchYear <- 2010
endYear    <- 2020
endMonth   <- 12

##############################
# Deal with some sea ice stuff
##############################
earlyRecord <- sprintf('%s/%s.%04d.dat',earlyDir,earlyFilePrefix,earlyIndVec[1])

# load in the data
dat1 <- readGridded(ddir, earlyRecord,
			startYear=1854,endYear=2009)
lon <- dat1$lon
lat <- dat1$lat

load('Data/landMask_merraGrid.Rda')

merraLon <- landMaskList$lon
merraLat <- landMaskList$lat

monthlyMask22 <- array(NA, dim=c(length(lon),length(lat),12))
maximalMask22 <- array(NA, dim=c(length(lon),length(lat)))
minimalMask22 <- array(NA, dim=c(length(lon),length(lat)))

for(i in 1:length(lon)){
	for(j in 1:length(lat)){
		lonInd <- which.min(abs(lon[i] - merraLon))	
		latInd <- which.min(abs(lat[j] - merraLat))

		for(s in 1:12){
			monthlyMask22[i,j,s] <- landMaskList$monthlyMask[lonInd,latInd,s]
		}
		maximalMask22[i,j] <- landMaskList$maximalMask[lonInd,latInd]
		minimalMask22[i,j] <- landMaskList$minimalMask[lonInd,latInd]
	}
}



processERSSTMember <- function(ensembleNumber){
# get the file paths as the data is already downloaded
earlyRecord <- sprintf('%s/%s.%04d.dat',earlyDir,earlyFilePrefix,earlyIndVec[ensembleNumber])
lateRecord  <- sprintf('%s/%s.%04d.dat',lateDir,lateFilePrefix,lateIndVec[ensembleNumber])

# load in the data
dat1 <- readGridded(ddir, earlyRecord,
			startYear=1854,endYear=2009)

timeMap1 <- cbind(rep(startYear1:endYear1,each=12),1:12,NA)
timeMap1[,3] <- timeMap1[,1] + (timeMap1[,2]-1)/12


dat2 <- readGridded(ddir, lateRecord,
			startYear=2001,endYear=2020)

timeMap2 <- cbind(rep(startYear2:endYear2,each=12),1:12,NA)
timeMap2[,3] <- timeMap2[,1] + (timeMap2[,2]-1)/12


lon <- dat1$lon
lat <- dat1$lat

#timeMap
inds1 <- which(timeMap1[,1] >= startYear & timeMap1[,1] < switchYear)
inds2 <- which(timeMap2[,1] >= switchYear & timeMap2[,1] <= endYear)

sstOut <- array(NA, dim=c(length(lon),length(lat),length(inds1)+length(inds2)))

sstOut[,,1:length(inds1)] <- dat1$sst[,,inds1]
sstOut[,,(length(inds1)+1):dim(sstOut)[3]] <- dat2$sst[,,inds2]

timeMap <- rbind(timeMap1[inds1,],timeMap2[inds2,])

# trims to mid-year if not through the end of a year
if(endMonth != 12){
	lastTimePoint <- which(timeMap[,1] == endYear & timeMap[,2] == endMonth)
	sstOut <- sstOut[,,1:lastTimePoint]
	timeMap <- timeMap[1:lastTimePoint,]
}


# correct if needed (ADDING the correction factor as it is defined)
if(correctClimatology){
	for(i in 1:dim(sstOut)[3]){
		monthInd <- timeMap[i,2]
		sstOut[,,i] <- sstOut[,,i] + climCorrectionArray[,,monthInd]
	}
	
}

# remove all of the land and sea ice regions from the sst
for(i in 1:nrow(timeMap)){
	monthInd <- timeMap[i,2]
	sstOut[,,i] <- sstOut[,,i] * ifelse(is.na(monthlyMask22[,,monthInd]),1,NA)
}


outList <- list(lon=dat1$lon,lat=dat1$lat,sst=sstOut)

# Pull the SBBX info
fileName <- "SBBX.SAMPLE"
templateList <- readFile(fileName,ddir)
allInfo <- pullInfo(templateList)

# merge grids and write
writeList <- mergeGrids(outList,allInfo)
writeListFix <- fixInfo(writeList)

# make the title
newTitle <- sprintf("Monthly Sea Surface Temperature anom (C) ERSSTv5 01/1880 - %02d/%04d              ",endMonth,endYear)
headInfo <- allInfo[1,]
headInfo[c(1,4)] <- writeListFix[[2]]$info[1]
headInfo[5]      <- writeListFix[[2]]$info[1] + 8

writeListFix[[1]] <- list(title= newTitle, info=headInfo)


# write the binary file (to the python working directory)
ofname <- sprintf("%s/SBBX.ensemble.%04d",ofdir,ensembleNumber)
writeFile(ddir,ofname,writeListFix)

}


# foreach call

cl <- makeCluster(nCores)
registerDoParallel(cl)

foreach(i=359:500) %dopar% processERSSTMember(i)

stopCluster(cl)



# rsync -avzh /Volumes/Archive/ERSSTv5/operationalEnsemble/SBBX_Output_2020_IceMask/ nlenssen@discover.nccs.nasa.gov:/discover/nobackup/nlenssen/ErsstEnsembleSBBX_2020_IceMask













