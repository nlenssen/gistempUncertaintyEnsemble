# GISTEMP Uncertainty Analysis
# Version 2.0.0 (Summer 2021)
# Nathan Lenssen (n.lenssen@columbia.edu)

library(doParallel)
library(foreach)

nCores <- 4

cdir <- 'Code/ERSST'

# location for the ERSST Data
ddir <- '/Users/lenssen/Documents/ERSSTv5_Ensemble'

earlyDir <- '/pub/data/cmb/ersst/v5/ensemble.1854-2017'
lateDir  <- '/pub/data/cmb/ersst/v5/ensemble.2001.present'

earlyFilePrefix <- 'sst2d.ano.1854.2017.ensemble'
lateFilePrefix  <- 'sst2d.ano.2001.last.ensemble'

source('Code/ERSST/Functions.R')

startYear1 <- 1854
endYear1   <- 2009

startYear2 <- 2001
endYear2   <- 2020

startYear  <- 1854
switchYear <- 2010
endYear    <- 2019

midclip   <- TRUE
clipYear  <- 2019
clipMonth <- 9

processERSSTMember <- function(ensembleNumber){
# download the data needed
earlyRecord <- downloadEnsemble(ensembleNumber, earlyDir, earlyFilePrefix, cdir, ddir)
lateRecord  <- downloadEnsemble(ensembleNumber, lateDir,  lateFilePrefix,  cdir, ddir)

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
# optional clip to mid year
if(midclip){
	lastTimePoint <- which(timeMap[,1] == clipYear & timeMap[,2] == clipMonth)
	sstOut <- sstOut[,,1:lastTimePoint]
	timeMap <- timeMap[1:lastTimePoint,]
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
newTitle <- sprintf("Monthly Sea Surface Temperature anom (C) ERSSTv5 01/1880 - %02d/%04d              ",clipMonth,clipYear,ensembleNumber)
headInfo <- allInfo[1,]
headInfo[c(1,4)] <- writeListFix[[2]]$info[1]
headInfo[5]      <- writeListFix[[2]]$info[1] + 8

writeListFix[[1]] <- list(title= newTitle, info=headInfo)


# write the binary file (to the python working directory)
ofname <- sprintf("SBBX_Output/SBBX.ensemble.%04d",ensembleNumber)
writeFile(ddir,ofname,writeListFix)

# remove the source data from the directory to save disk space
system(sprintf('rm %s/%s',ddir,earlyRecord))
system(sprintf('rm %s/%s',ddir,lateRecord))
}


# foreach call

cl <- makeCluster(nCores)
registerDoParallel(cl)

foreach(i=81) %dopar% processERSSTMember(i)

stopCluster(cl)

