source('Namelists/fullEnsemble.Rnl')

# Overwrite the variogram settings from the namelist and set some scripts params
variogramTesting <- FALSE
maxDistance      <- 850
nSims            <- 60
nBins            <- 15


###############################################################################
# Functions used
###############################################################################

# requires:
# - lon,lat
# - dataIndsList
# - sdVecList

vgramTest <- function(d, covFunction, dmax=550, nSims = 40, ...){
	dataInds <- dataIndsList[[d]]

	lonLatGrid <- cbind(lon[dataInds[,1]],lat[dataInds[,2]])

	# create the matern matrix with unit variances (diagonal) 
	# Note: currently fitting in a very ad hoc method to get initial results
	# 250,1 are the variogram-fitted values for the error fields
	covMat <- stationary.cov(lonLatGrid, Covariance = covFunction,
								Distance = "rdist.earth", ...)

	# construct the matrix with 0 nugget as this seems to be a reasonable approx
	# from the variogram testing
	sdVec <- c(uncEstimate[,,d])
	sdVec <- sdVec[!is.na(sdVec)]

	sdVecList[[d]] <- sdVec

	# scale the Matern covariance matrix by the empirical gridbox variance
	covMatFinal <- diag(sdVec) %*% covMat %*% t(diag(sdVec))

	# simulate from the cov mat to see how the variogram lines up
	test <- rmvnorm(n=nSims, mu=rep(0,nrow(covMatFinal)),covMatFinal)

	simGramTest <- list()

	for(i in 1:nSims){
		simGramTest[[i]] <- vgram(lonLatGrid,test[i,], N=15, lon.lat=TRUE,dmax=dmax)	
	}

	return(simGramTest)
}




###############################################################################
# populate the workspace with objects from the spatial uncertainty analysis
###############################################################################
source('Code/04SpatialUncertaintyCoarse.R')

# Overwrite the variogram settings from the namelist and set some scripts params
variogramTesting <- FALSE
maxDistance      <- 850
nSims            <- 60
nBins            <- 15



###############################################################################
# run a sweep over parameters generating simulated variograms
###############################################################################

rangeSweep <- seq(100, 300, by = 50)
smoothnessSweep <- seq(1, 2, by = 0.5)
trialMat <- as.matrix(expand.grid(list(range=rangeSweep, smoothness=smoothnessSweep)))

for(dec in 1:14){
	outList <- list()

	for(i in 1:nrow(trialMat)){
			 tempGram <- vgramTest(d=dec, dmax=maxDistance, nSims = nSims,
				covFunction= 'Matern',range=trialMat[i,1], smoothness=trialMat[i,2])

			 # remove the full info from the vgram object
			 tempGram$d <- NULL
			 tempGram$vgram <- NULL

			 outList[[i]] <- tempGram
	}

	save(outList, file=sprintf('Data/CovSweepTests/outlist%02d.Rda',dec))
}



###############################################################################
# Get the empirical variograms from the observed difference matricies
###############################################################################

timePoints <- sample(1:ncol(diffMat), nSims)

for(dec in 1:14){

	workingDiffMat <- diffMatList[[dec]] 
	dataInds <- dataIndsList[[dec]]
	lonLatGrid <- cbind(lon[dataInds[,1]],lat[dataInds[,2]])

	# the variograms for the various difference fields and the simulated
	# matern covariance
	variogramSingleTime <- list()


	for(i in 1:nSims){
		variogramSingleTime[[i]] <- vgram(lonLatGrid,workingDiffMat[,timePoints[i]], N= nBins,
													lon.lat=TRUE, dmax= maxDistance)

	}

	# clear the really mem/disk intesive parts of the vgram object
	variogramSingleTime$d <- NULL
	variogramSingleTime$vgram <- NULL

	save(variogramSingleTime, file=sprintf('Data/CovSweepTests/obsGram%02d.Rda',dec))

}

###############################################################################
# Plotting stuff
###############################################################################


for(dec in 1:14){
	load(sprintf('Data/CovSweepTests/outlist%02d.Rda',dec))
	load(sprintf('Data/CovSweepTests/obsGram%02d.Rda',dec))

	pdf(sprintf('Figures/CovSweepTests/testPlot%02d.pdf',dec),6,6)

	plot(NULL, xlim=c(0,600), ylim=c(-0.2,4),
		xlab='Distance (km)', ylab='Variance',
		main = sprintf('Variogram for %ss LSAT Error Field',tDec[dec]))

	pal <- designer.colors(length(outList), brewer.pal(9, 'YlOrBr'))

	for(i in 1:length(outList)){
		addVgram(outList[[i]], pal[i])
	}

	# the ones from the actual data so far (need to run the full other code)
	addVgram(variogramSingleTime, 'black')

	dev.off()
}


###############################################################################
# Try minimizing MSE of each of them to see which is the best
###############################################################################
resultsList <- list()

for(dec in 1:14){

	load(sprintf('Data/CovSweepTests/outlist%02d.Rda',dec))
	load(sprintf('Data/CovSweepTests/obsGram%02d.Rda',dec))


	# get the empirical varigram info
	xValsEmp <- variogramSingleTime[[1]]$centers

	vgramEmpMean <- matrix(NA, length(xValsEmp),length(variogramSingleTime))

	for(i in 1:length(variogramSingleTime)){
		vgramEmpMean[,i] <- variogramSingleTime[[i]]$stats[2,]
	}

	vgramEmpQuants <- apply(vgramEmpMean,1,quantile,prob=c(0.1, 0.25, 0.5, 0.75, 0.9))
	medianEmpEst <- vgramEmpQuants[3,]

	# loop over all of the trial matern patterns
	medianSweepEstMat <- matrix(NA, length(xValsEmp), nrow(trialMat))

	for(i in 1:nrow(trialMat)){

		# get the sweep variogram info
		vlist <- outList[[i]]
		xVals <- vlist[[1]]$centers

		# run some tests to make sure the empirical and simulated vgrams are compatible
		if(!all(xVals == xValsEmp)){
			print('Variogram Bin Mismatch')
			break()
		}

		if(length(variogramSingleTime) != length(vlist)){
			print('Number of time points mismatch')
			break()
		}

		# calculate the statistics of the sweep variogram
		vgramSweepMean <- matrix(NA, length(xVals),length(vlist))

		for(j in 1:length(vlist)){
			vgramSweepMean[,j] <- vlist[[j]]$stats[2,]
		}

		vgramSweepQuants <- apply(vgramSweepMean,1,quantile,prob=c(0.1, 0.25, 0.5, 0.75, 0.9))
		medianSweepEstMat[,i] <- vgramSweepQuants[3,]
	}

	# calculate the rmse for each of the cases (unweighted for now, could weight by obs)
	rmseVec <- rep(NA, ncol(medianSweepEstMat))

	for(i in 1:ncol(medianSweepEstMat)){
		rmseVec[i] <- sqrt(sum((medianEmpEst - medianSweepEstMat[,i])^2))
	}

	resultsList[[dec]] <- list(empirical = medianEmpEst, sweep = medianSweepEstMat,
							  rmse = rmseVec, sweepParams = trialMat)
}

save(resultsList, file='Data/CovSweepTests/resultsList.Rda')

# plot them in 2d space
par(mfrow=c(4,4))

for(dec in 1:14){
	image.plot(rangeSweep, smoothnessSweep,
		matrix(resultsList[[dec]]$rmse,length(rangeSweep), length(smoothnessSweep)),
		zlim=c(0,2), main=tDec[dec])

	minInd <- which.min(resultsList[[dec]]$rmse)
	points(trialMat[minInd,1],trialMat[minInd,2], col='red', cex=2)
}








