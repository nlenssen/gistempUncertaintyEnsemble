
# variogram stuff
maxDistance <- 850
nBins       <- 10

# matern stuff
smoothness <- 1
earlyRange <- 250
lateRange  <- 200

###############################################################################
# load and preprocess the 2x2 ERA analysis
###############################################################################

load(sprintf('%s/Intermediate/ERA5/anomalyData_ERA_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))
load(sprintf('%s/Intermediate/GistempProduction/coverageInfo.Rda',scratchDir))

# unpack the ERA list
lon <- anomalyData$lon
lat <- anomalyData$lat
nlon <- length(lon)
nlat <- length(lat)
timeMapERA <- anomalyData$timeMap
anomalyField <- anomalyData$anomalyField
zoneMask <- anomalyData$zoneMask

rm(anomalyData)

# time stuff
nt <- dim(anomalyField)[3]

# land-mask the annualField (could eventually make this a monthly mask)
landMask <- landMaskList$maximalMask

for(i in 1:12){
	subInds  <- which(timeMapERA[,2]==i)
	for(k in subInds){
		anomalyField[,,k] <- anomalyField[,,k] * landMask # drop in monthly mask here indexed by i
	}
}
 



###############################################################################
# Generate the difference series for each decade of interest
###############################################################################

# Get the interpolated field paths
files <- system(sprintf('ls %s/Intermediate/ERA5/InterpolatedFieldsNetcdf/interpField*_2x2.nc',
	scratchDir),intern=T)


fullIndList <- make.surface.grid(list(lon=lon,lat=lat))
landIndList <- fullIndList[which(landMask==1),]
# save the point estimate of the uncertainty variance
uncEstimate <- array(NA, dim=c(length(lon),length(lat),length(decInds)))

diffMatList <- list()
dataIndsList <- list()

for(d in 1:length(decInds)){
	# Load and process the interpolated field
	interpolatedField <- read22File(files[decInds[d]])$anomArray

	# mask by coverage in production gistemp
	for(k in 1:nt){
		interpolatedField[,,k] <- interpolatedField[,,k] * decadalCoverageMask[,,d]
	}

	# place to store all of the difference series
	diffSeriesArr <- array(NA, dim=c(length(lon),length(lat),nt))


	for(i in 1:nrow(landIndList)){
		lonInd <- which(lon==landIndList[i,1])
		latInd <- which(lat==landIndList[i,2])

		trueVec <- anomalyField[lonInd,latInd,]
		maskVec <- interpolatedField[lonInd,latInd,]
		if(all(is.na(trueVec)))	next()
		if(all(is.na(maskVec))){
			uncEstimate[lonInd,latInd,d] <- NA
			next()
		}
	
		diffVec <- trueVec - maskVec
		diffSeriesArr[lonInd,latInd,] <- diffVec
		uncEstimate[lonInd,latInd,d] <- sd(diffVec,na.rm=T)

	}

	# store locations where we have data
	dataInds     <- which(!is.na(uncEstimate[,,d]),arr.ind=T)
	numLocations <- nrow(dataInds)

	diffMat <- matrix(NA,nrow=numLocations,ncol=nt)

	for(k in 1:nrow(diffMat)){
		tempInd <- dataInds[k,]
		diffMat[k,] <- diffSeriesArr[tempInd[1],tempInd[2],]
	}

	diffMatList[[d]] <- diffMat 
	dataIndsList[[d]] <- dataInds
}



###############################################################################
# Generate the Matern Covariance Mat explicity
###############################################################################

sdVecList <- list()

for(d in 1:length(decInds)){
	dataInds <- dataIndsList[[d]]

	lonLatGrid <- cbind(lon[dataInds[,1]],lat[dataInds[,2]])

	# create the matern matrix with unit variances (diagonal) 
	# Note: currently fitting in a very ad hoc method to get initial results
	# 250,1 are the variogram-fitted values for the error fields
	if(tDec[d] < 1960){
		range <- earlyRange
	} else{
		range <- lateRange
	}
	
	covMat <- stationary.cov(lonLatGrid, Covariance = "Matern",
									Distance = "rdist.earth", range=range,smoothness=smoothness)

	# construct the matrix with 0 nugget as this seems to be a reasonable approx
	# from the variogram testing
	sdVec <- c(uncEstimate[,,d])
	sdVec <- sdVec[!is.na(sdVec)]

	sdVecList[[d]] <- sdVec

	# scale the Matern covariance matrix by the empirical gridbox variance
	covMatFinal <- diag(sdVec) %*% covMat %*% t(diag(sdVec))
	save(covMatFinal, file=sprintf('%s/Intermediate/CovarianceMats/covMatFinal_%02d.Rda',scratchDir,d))
}

# save the important output
save(lon,lat, dataIndsList, diffMatList, uncEstimate, sdVecList,
	file=sprintf('%s/Intermediate/CovarianceMats/spatialAnalysisFullData.Rda',scratchDir))


###############################################################################
# (Optional) Calculate some variograms to get a sense of what is going on
###############################################################################

if(variogramTesting){

	varInds <- 1:14
	nPoints <- 60
	timePoints <- sample(1:ncol(diffMat),nPoints)

	makeVarGrams <- function(dec){
		require(Rfast)
		require(fields)
		# pick some time points
		

		# load matching decades and get the necessary objects
		load(sprintf('Data/CovarianceMats/covMatFinal_%02d.Rda',dec))
		workingDiffMat <- diffMatList[[dec]] 
		dataInds <- dataIndsList[[dec]]
		lonLatGrid <- cbind(lon[dataInds[,1]],lat[dataInds[,2]])

		# the variograms for the various difference fields and the simulated
		# matern covariance
		variogramSingleTime <- list()
		simGram <- list()

		covariogramSingleTime <- list()
		simCovGram <- list()

		# simulate draws from the multivariate normal as given by the covariance
		# matrix made in the previous step (USE Rfast TO MAKE SUPER QUICK)
		test <- rmvnorm(n=nPoints, mu=rep(0,nrow(covMatFinal)),covMatFinal)


		for(i in 1:length(timePoints)){
			variogramSingleTime[[i]] <- vgram(lonLatGrid,workingDiffMat[,timePoints[i]], N= nBins,
													lon.lat=TRUE, dmax= maxDistance)


			simGram[[i]]             <- vgram(lonLatGrid,test[i,], N= nBins,
													lon.lat=TRUE, dmax= maxDistance)

			# covariogramSingleTime[[i]] <- vgram(lonLatGrid,workingDiffMat[,timePoints[i]], N= nBins,
			# 										lon.lat=TRUE,dmax= maxDistance,
			# 										type='correlogram')


			# simCovGram[[i]]             <- vgram(lonLatGrid,test[i,], N= nBins,
			# 										lon.lat=TRUE,dmax= maxDistance,
			# 										type='correlogram')
		}
		

		save(variogramSingleTime,simGram,covariogramSingleTime,simCovGram,
			file=sprintf('%s/Intermediate/CovarianceMats/varCorGramData_%02d.Rda',scratchDir,dec))
	}


	for(i in varInds){
		makeVarGrams(varInds[i])
		gc()
	}


	# system('say -v Tessa "all done!"')

}


# Variogram Plotting stuff
if(plotting & variogramTesting){
	plotInds <- 1:14

	for(dec in plotInds){
		load(sprintf('%s/Intermediate/CovarianceMats/varCorGramData_%02d.Rda',scratchDir,dec))

		tDecSub <- seq(1880,2010, by=10)

		pdf(sprintf('%s/varCorGram_%02d.pdf',plotdir,dec),10,6)
		par(mfrow=c(1,1))
		# variogram
		plot(NULL, xlim=c(0,maxDistance), ylim=c(-0.2,4),
			xlab='Distance (km)', ylab='Variance',
			main = sprintf('Variogram for %ss LSAT Error Field',tDecSub[dec]))

		addVgram(variogramSingleTime, 'black')
		addVgram(simGram, 'red')	
		legend('bottomright', c('Empirical', 'Gaussian Process Approximation'),
			lwd=2, col=c('black', 'red'))

		# covariogram
		# plot(NULL, xlim=c(0,maxD), ylim=c(-0.2,4),
		# 	xlab='Distance (km)', ylab='Covariance',
		# 	main = sprintf('Correlogram for %ss LSAT Error Field',tDecSub[dec]))

		# addVgram(covariogramSingleTime, 'black')
		# addVgram(simCovGram, 'red')
		# legend('bottomright', c('Empirical', 'Gaussian Process Approximation'),
		# 	lwd=2, col=c('black', 'red'))
		
		dev.off()	

	}




	# make a final plot for the paper
	plotInds <- c(3,8,12)
	cexScale <- 1.25

	pdf(sprintf('%s/varCorGram_Paper.pdf',plotdir),11,4)
		par(mfrow=c(1,3))
	
	tDecSub <- seq(1880,2010, by=10)

	for(dec in plotInds){
		load(sprintf('%s/Intermediate/CovarianceMats/varCorGramData_%02d.Rda',scratchDir,dec))


		# variogram
		plot(NULL, xlim=c(0,maxDistance), ylim=c(-0.2,4),
			xlab='Distance (km)', ylab='Variance (°C^2)',
			main = sprintf('Variogram for %ss LSAT Error Field',tDecSub[dec]), 
 		    cex.lab=cexScale, cex.axis=cexScale, cex.main=cexScale)


		addVgram(variogramSingleTime, 'black')
		addVgram(simGram, 'red')	
		if(dec==8) legend(0,4, c('Empirical Variogram', 'Simulated Mátern Covariogram'),
			lwd=2, col=c('black', 'red'), bty='n')	
	}

	dev.off()	
}
















