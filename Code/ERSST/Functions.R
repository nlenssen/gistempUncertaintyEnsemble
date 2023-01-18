###############################################################################
# Function used to read and write  the ERSST binary files
###############################################################################

readGridded <- function(ddir,fileName,startYear,endYear){
	# make the 2x2 grid
	nx <- 180
	ny <- 89
	nt <- 12*(endYear-startYear+1)


	shift <- 0 # Needed to get the images to align with boundries. Should check on this
	rawLon <- c(seq(0+shift,180,by=2),seq(-178,-2+shift,by=2))
	perm   <- rank(rawLon)
	lon <- seq(-178,180,by=2)
	lat <- seq(-88,88,by=2)

	# The full ssta array for a given sample
	datArray <- array(NA, dim=c(nx,ny,nt))
	# the data is presented in 4 byte floats (big endian) with start and end caps (ints)
	# open the binary file stream in read binary mode (rb)
	fileStream <- file(sprintf("%s/%s",ddir, fileName),open="rb")

	for(i in 1:nt){
		startCap <- readBin(fileStream,integer(),size=4,n=1,endian='big')
		datRaw   <- readBin(fileStream,numeric(),size=4,n=(nx*ny),endian='big')
		endCap   <- readBin(fileStream,integer(),size=4,n=1,endian='big')

		if(startCap != endCap) stop('Cap Mismatch')
		
		datMat   <- matrix(NA, nx,ny)
		datMat[perm,] <- matrix(datRaw,nx,ny)
		datMat[datMat < -900] <- NA

		datArray[,,i] <- datMat
	}

	close(fileStream)

	return(list(lon=lon,lat=lat,sst=datArray))
}


writeGridded <- function(dataArray, ddir,fileName){
	# make the 2x2 grid
	nx <- 180
	ny <- 89
	nt <- 12*(2014-1854+1)

	writeSize <- as.integer(nx*ny*4)

	targetLon <- c(seq(0,180,by=2),seq(-178,-2,by=2))
	dataLon <- seq(-178,180,by=2)

	ind <- which(targetLon==180)

	perm <- c((ind-1):nx,1:(ind-2))


	# the data is presented in 4 byte floats (big endian) with start and end caps (ints)
	# open the binary file stream in read binary mode (rb)
	fileStream <- file(sprintf("%s/%s",ddir, fileName),open="wb")

	for(i in 1:nt){
		writeMat <- c(dataArray[perm,,i])
		writeMat[is.na(writeMat)] <- -999.9

		writeBin(writeSize,fileStream,size=4,endian='big')
		writeBin(writeMat,fileStream,size=4,endian='big')
		writeBin(writeSize,fileStream,size=4,endian='big')
	}

	close(fileStream)

}



################################################################################
# Functions used to merge the two grids
################################################################################

mergeGrids <- function(ensembleDat,allInfo){
	glon <- ensembleDat$lon
	glat <- ensembleDat$lat
	grid <- ensembleDat$sst


	writeList <- list()
	for(i in 1:8000){
	# pull the lat/lon boundries of the subbox
	infoRow <- allInfo[i+1,]
	sblon <- infoRow[4:5]/100
	sblat <- infoRow[2:3]/100

	# find the corresponding gridbox(s) of the gridded file
	# MAY NEED TO REDO THIS!!! think about more
	lonRange <- c(floor(sblon[1]/2),ceiling((sblon[2])/2))*2
	latRange <- c(floor(sblat[1]/2),ceiling(sblat[2]/2))*2

	lonInds <- which(glon %in% lonRange)
	latInds <- which(glat %in% latRange)

	# check for one one range (means poles or lon discontinuity)
	if(length(lonInds)<2){
		lonInds <- sort( c(lonInds, ifelse(lonRange[1]>0,length(glon),1)) )
	}
	if(length(latInds)==1){
		latInds <- sort( c(latInds,ifelse(latRange[1]>0,length(glat),1)) )
	}


	lonIndsRange <- lonInds[1]:lonInds[2]
	latIndsRange <- latInds[1]:latInds[2]

	subGrid <- grid[lonIndsRange,latIndsRange,]

	# Now we need to area weight the boxes to combine into the sbbx time series

	tsweight <- matrix(NA,nrow(subGrid),ncol(subGrid))

	for(j in 1:nrow(subGrid)){
		for(k in 1:ncol(subGrid)){
			tempGriddedLoc <- c(lonIndsRange[j],latIndsRange[k])

			tempglon <- glon[tempGriddedLoc[1]]
			tempglat <- glat[tempGriddedLoc[2]]

			tempglonRange <- c(tempglon-1,tempglon+1)
			tempglatRange <- c(tempglat-1,tempglat+1)

			fullgArea <- gbArea(tempglonRange,tempglatRange)

			# calculate the intersection of the two gridboxes
			overlap <- gbIntersection(tempglonRange,tempglatRange,sblon,sblat)
			overlapArea <- gbArea(overlap$lon,overlap$lat)

			# calculate the weigting (set to zero if negative)
			prop <- overlapArea/fullgArea
			tsweight[j,k] <- ifelse(prop>0, prop, 0)
		}
	}


	# now return the subbox time series
	sbTS <- timeSeriesCalc(subGrid,tsweight)

	# set the missing values to 9999
	sbTS[is.na(sbTS)] <- 9999
	sbTS[is.nan(sbTS)] <- 9999


	#trim the early years off
	sbTSTrim <- sbTS[(12*(1880-1854)+1):length(sbTS)]


	# write totally missing data as a single 9999
	if(all(sbTSTrim == 9999)) sbTSTrim <- 9999

	
	# update the info file
	infoRowFinal <- infoRow
	infoRowFinal[1] <- length(sbTSTrim)
	infoRowFinal[7] <- length(sbTSTrim) - sum(sbTSTrim == 9999)

	writeList[[i+1]] <- list(data=sbTSTrim,info=infoRowFinal)
	} #end big loop

return(writeList)
}



compareTwo <- function(list1,list2,index){
	y1 <- list1[[index]]$data[1:(length(list1[[index]]$data)-24)]
	y2 <- list2[[index]]$data[(12*(1880-1854)+1):length(list2[[index]]$data)]

	y1[y1==9999] <- NA
	y2[y2==9999] <- NA

	yr <- range(y1,y2,na.rm=TRUE)
	par(mfrow=c(2,1))
	plot(y1,type='l',col='black',lwd=2)
	points(y2,type='l',col='blue',lwd=2)

	plot(y1-y2,type='l',col='red')
	abline(h=0)
}


# functions used in analysis
timeSeriesCalc <- function(mat,weights){
	tsOut <- rep(NA, dim(mat)[3])
	for(t in 1:dim(mat)[3]){
		tsOut[t] <- weighted.mean(mat[,,t],weights,na.rm=TRUE)
	}
	return(tsOut)
}


gbArea <- function(lon,lat){
	degToRad <- pi/180
	return((lon[2]*degToRad-lon[1]*degToRad) * (sin(lat[2]*degToRad) - sin(lat[1]*degToRad)))
}


gbIntersection <- function(lon1,lat1,lon2,lat2){
	lonOut <- c(max(c(lon1[1],lon2[1])), min(c(lon1[2],lon2[2])))
	latOut <- c(max(c(lat1[1],lat2[1])), min(c(lat1[2],lat2[2])))

	# crude fixing the intervals...
	if(lonOut[2] < lonOut[1]) lonOut[2] <- lonOut[1]
	if(latOut[2] < latOut[1]) latOut[2] <- latOut[1]

	return(list(lon=lonOut,lat=latOut))
}

################################################################################
# Functions used to read the binary files
################################################################################

readFile <- function(fileName, ddir){
	handle <- file(sprintf("%s/%s",ddir,fileName),open='rb')

	outList <- list()

	# read the first segment
	titleLength <- readBin(handle,integer(),n=1,endian='big')
	info <- readBin(handle,integer(),n=8,endian='big')
	title <- readChar(handle,nchars=80)

	outList[[1]] <- list(title=title,info=info)
	titleEnd <- readBin(handle,integer(),n=1,endian='big')

	# handle the rest of the boxes
	for(i in 1:8000){
		readLength <- readBin(handle,integer(),n=1,endian='big')
		outList[[i+1]] <- readBox(handle,readLength)
	}

	close(handle)

	return(outList)
}

readBox <- function(file,readLength){
	boxInfo <- readBin(file,integer(),n=8,endian='big')
	data    <- readBin(file,double(),size=4,n=readLength/4-8,endian='big')
	endCap <- readBin(file,integer(),n=1,endian="big")
	
	try(if(endCap != readLength) stop("Wrong read Lenth somewhere!!!"))

	return(list(data=data,info=boxInfo))
}


pullInfo <- function(megaList){
	numBoxes <- length(megaList)

	infoMat <- matrix(NA, nrow=numBoxes,ncol=8)

	for(i in 1:numBoxes){
		infoMat[i,] <- megaList[[i]]$info
	}
	return(infoMat)
}


################################################################################
# Functions used to write the binary files
################################################################################
writeFile <- function(ofdir, fileName, writeList){
	# open the conneciton with the file
	handle <- file(sprintf("%s/%s",ofdir,fileName),open='wb')

	# write the head
	writeHead(handle,writeList[[1]]$info,writeList[[1]]$title)

	# write the data
	for(i in 1:8000){
		slice <- writeList[[i+1]]
		writeBox(handle,slice$info,slice$data)	
	}
	

	# close the connection with the file
	close(handle)
}

writeHead <- function(handle,info, title){
	# calculate the length of the write
	byteSize <- length(info)*4 + nchar(title)

	# Write the headcap and the info 
	writeBin(as.integer(c(byteSize,info)),handle,endian='big')
	
	# write the title (eos: End of String=NULL for no end cap)
	writeChar(title,handle,eos=NULL)

	# write the end cap
	writeBin(as.integer(byteSize),handle,endian='big')
}

writeBox <- function(handle, info, ts){
	byteSize <- length(info) * 4 + length(ts) * 4

	# write the headcap and the info (integers)
	writeBin(as.integer(c(byteSize,info)),size=4,handle,endian='big')

	# write the data (Re*4)
	writeBin(ts,handle,size=4,endian='big')

	# write the endcap
	writeBin(as.integer(byteSize),size=4,handle,endian='big')

}

# an attempt to shift the first index of all of the data info vectors to the
# previous record. (May fix the issues I'm having with the python code)
fixInfo <- function(writeList){

	for(i in 1:(length(writeList)-1)){
		writeList[[i]]$info[1] <- writeList[[i+1]]$info[1]
	}

	writeList[[length(writeList)]]$info[1] <- 0

	return(writeList)
}

################################################################################
# Functions used to download the ensemble data from the NOAA server
################################################################################

downloadEnsemble <- function(ensembleNumber, remotedir, filePrefix,
							cdir="/Users/nlenssen/Dropbox/NASA/ERSST/Code",
							ddir="/Users/nlenssen/Dropbox/NASA/ERSST/Data/gridded"){
	
	info <- sprintf("%s.%04d.dat", filePrefix, ensembleNumber)
	outName <- sprintf("%s/%s.%04d.dat",ddir, filePrefix, ensembleNumber)

	# download the file from noaa and save to ddir
	system(sprintf("%s/pullData.sh %s %s %s", cdir, remotedir, info, outName),intern=TRUE)

	# return the file name
	return(sprintf("%s.%04d.dat", filePrefix, ensembleNumber))
}


################################################################################
# Function to foreach the merge step
################################################################################
indexedMerge <- function(ensembleNumber, allInfo, delete=FALSE,
					ddir="/Users/nlenssen/Documents/ensemble",
					ofdir="/Users/nlenssen/Documents/ensemble/sbbx"){
	cat(paste("Index:",ensembleNumber,"\n"))
	tick <- proc.time()
	
	# Make the name and try to unzip the file
	ensembleName <- sprintf("sst2d.ano.1854.2014.ensemble.%04d.dat",ensembleNumber)
	try(system(sprintf("gunzip %s/%s.gz",ddir,ensembleName)),silent=TRUE)
	
	# Read the gridded file
	ensembleDat <- readGridded(ddir,ensembleName)

	# Option to delete the ensemble file
	if(delete) system(sprintf("rm %s/%s",ddir,ensembleName))

	# Merge, and fix the sbbx list object
	writeList <- mergeGrids(ensembleDat,allInfo)
	writeListFix <- fixInfo(writeList)

	# make the title
	newTitle <- sprintf("Monthly Sea Surface Temperature anom (C) ERSSTv4 01/1880 - 12/2016 ensemble%04d ",ensembleNumber)
	headInfo <- allInfo[1,]
	headInfo[c(1,4)] <- writeListFix[[2]]$info[1]
	headInfo[5]      <- writeListFix[[2]]$info[1] + 8

	writeListFix[[1]] <- list(title= newTitle, info=headInfo)

	# write the binary file (to the python working directory)
	ofname <- sprintf("SBBX.ensemble.%04d",ensembleNumber)
	writeFile(ofdir,ofname,writeListFix)

	tock <- proc.time()
	gc()

	return(tock[3]-tick[3])
}


globalMixed <- function(ensembleNumber,
				ddir      = "/Users/nlenssen/Documents/ensemble/sbbx",
				pyDir     = "/Users/nlenssen/Documents/gistemp1.0",
				cdir      = "/Users/nlenssen/Dropbox/NASA/ERSST/Code",
				resultDir = "/Users/nlenssen/Documents/resultsSST"){

	cat(paste("Working Ensemble:" , ensembleNumber, "\n" ))
	tick <- proc.time()
	###############################################################################
	# copy sbbx file to python wd
	###############################################################################

	system(sprintf("cp %s/SBBX.ensemble.%04d %s/tmp/input",ddir,ensembleNumber,pyDir))

	###############################################################################
	# Run the ocean portion of GISTEMP analysis
	###############################################################################

	# run steps 4 and 5 of gistemp [60 sec]
	setwd(pyDir) # Need to be in the root directory of the project
	system(sprintf("python3 tool/run.py -s 4-5 -p ocean_source=ensemble.%04d",ensembleNumber))


	# make a new directory for the output
	system(sprintf("mkdir %s/results%04d",resultDir,ensembleNumber))

	# move all of the mixed, ensemble results
	system(sprintf("mv %s/tmp/result/mixed*ENSEMBLE.%04d*.txt %s/results%04d",pyDir,ensembleNumber,resultDir,ensembleNumber))

	###############################################################################
	# Cleanup
	###############################################################################

	# remove the SBBX file from the python file tree
	system(sprintf("rm %s/tmp/input/SBBX.ensemble.%04d",pyDir,ensembleNumber))

	# cleanup the /temp/result directory (save SBBX1880.Ts.GHCN.CL.PA.1200.npz)
	system(sprintf("rm %s/tmp/result/%s",pyDir,"SBBX.SST.npz"))
	system(sprintf("rm %s/tmp/result/%s",pyDir,"land*"))
	system(sprintf("rm %s/tmp/result/%s",pyDir,"GHCNv3BoxesLand.1200.txt"))


	system(sprintf("rm %s/tmp/result/%s",pyDir,"mixed*npz"))

	tock <- proc.time()

	return(tock[3]-tick[3])
}
