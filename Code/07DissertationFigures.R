###############################################################################
# Get the filepaths and grid
###############################################################################
filePaths <- system(sprintf('ls %s/*.nc',ensembleOutDir),intern=T)

handle <- nc_open(filePaths[1])

lon <- ncvar_get(handle, 'lon')
lat <- ncvar_get(handle, 'lat')

nc_close(handle)


###############################################################################
# Plot the mean global series
###############################################################################

# get this one
prodMeanTable <- read.csv('Data/ZonAnn.Ts+dSST.csv')
prodVars <- names(prodMeanTable)

# load the uncertainty ensemble mean series
load(sprintf('%s/Output/EnsembleStats/meanSeries.Rda',scratchDir))

# the dim of the uncertainty ensemble mean files
ndim1 <- dim(ensembleGlobalMean)[2]
ndim2 <- dim(ensembleGlobalMean)[3]

global <- apply(ensembleGlobalMean, 1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
nh     <- apply(ensembleNHemMean, 1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
sh     <- apply(ensembleSHemMean, 1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
bands  <- apply(ensembleBandMean, c(1,2), quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)

globalAnnual <- apply(apply(ensembleGlobalMean, c(2,3), monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
nhAnnual     <- apply(apply(ensembleNHemMean, c(2,3), monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
shAnnual     <- apply(apply(ensembleSHemMean, c(2,3), monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)

bandAnnualTemp <- array(NA, dim=c(nYear, 8, ndim1, ndim2))
for(i in 1:ndim1){
	for(j in 1:ndim2){
		for(b in 1:8){
			bandAnnualTemp[,b,i,j] <- monthToYearMeans(ensembleBandMean[,b,i,j])
		}
	}
}

bandsAnnual <- apply(bandAnnualTemp,c(1,2), quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)

xlabel <- 'Year'
ylabel <- 'Annual Mean Temperature Anomaly (°C)'
yr <- c(-0.55, 1.35)


# global and hemi lighter weight
pdf(sprintf('%s/Gistemp_GlobalHemi_Annual.pdf',plotdir),6,9)
set.panel(3,1)
#global
plot(Glob~Year, data=prodMeanTable,type='l',
	xlab=xlabel, ylab= ylabel, main='Global Mean',ylim=yr)

points(tYear, globalAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalAnnual[1,], rev(globalAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(Glob~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

#NH
plot(NHem~Year, data=prodMeanTable,type='l',
	xlab=xlabel, ylab= ylabel, main='Northern Hemisphere Mean',ylim=yr)

points(tYear, nhAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(nhAnnual[1,], rev(nhAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(NHem~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

#SH
plot(SHem~Year, data=prodMeanTable,type='l',
	xlab=xlabel, ylab= ylabel, main='Southern Hemisphere Mean',ylim=yr)

points(tYear, shAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(shAnnual[1,], rev(shAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(SHem~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

dev.off()


# bands


colInds <- c(8,15,9,14,10,13,11,12)
bandsInds <- c(8,1,7,2,6,3,5,4)

mainVec <- paste(rep(c('NH', 'SH'),4), 
				 rep(c('Polar Mean (64-90)', 'Mid-Lat Mean (44-64)',
				       'Sub-Tropic Mean (24-44)', 'Tropic Mean (0-24)'), each=2))
xlabel <- 'Year'
ylabel <- 'Annual Mean Temperature Anomaly (°C)'

antarcicaInds <- which(prodMeanTable$Year > 1959)
antarcicaIndsEnsemble <- which(tYear > 1959)
yr <- c(-0.75,1.75)

pdf(sprintf('%s/Gistemp_Bands_Annual.pdf',plotdir),10,12)
set.panel(4,2)

i <- 1
plot(prodMeanTable[,1],prodMeanTable[,colInds[i]],
	type='l',xlim=range(prodMeanTable$Year),
	xlab=xlabel, ylab= ylabel, main=mainVec[i])

points(tYear, bandsAnnual[2,,bandsInds[i]],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(bandsAnnual[1,,bandsInds[i]], rev(bandsAnnual[3,,bandsInds[i]])),
	col=adjustcolor('red',alpha=0.3), border=NA)

points(prodMeanTable[,1],prodMeanTable[,colInds[i]],type='l',lwd=1.5)
abline(h=0,lty=3)
legend(1880,2,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')


i <- 2
plot(prodMeanTable[antarcicaInds,1],prodMeanTable[antarcicaInds,colInds[i]],
	type='l',xlim=range(prodMeanTable$Year),ylim=yr,
	xlab=xlabel, ylab= ylabel, main=mainVec[i])

points(tYear[antarcicaIndsEnsemble], bandsAnnual[2,antarcicaIndsEnsemble,bandsInds[i]],col='red',lwd=1.5,type='l')
polygon(c(tYear[antarcicaIndsEnsemble], rev(tYear[antarcicaIndsEnsemble])),
	c(bandsAnnual[1,antarcicaIndsEnsemble,bandsInds[i]], rev(bandsAnnual[3,antarcicaIndsEnsemble,bandsInds[i]])),
	col=adjustcolor('red',alpha=0.3), border=NA)

points(prodMeanTable[antarcicaInds,1],prodMeanTable[antarcicaInds,colInds[i]],type='l',lwd=1.5)
abline(h=0,lty=3)
legend(1880,1,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')


for(i in 3:length(colInds)){
	plot(prodMeanTable[,1],prodMeanTable[,colInds[i]],
		type='l',xlim=range(prodMeanTable$Year),ylim=yr,
	xlab=xlabel, ylab= ylabel, main=mainVec[i])

	points(tYear, bandsAnnual[2,,bandsInds[i]],col='red',lwd=1.5,type='l')
	polygon(c(tYear, rev(tYear)), c(bandsAnnual[1,,bandsInds[i]], rev(bandsAnnual[3,,bandsInds[i]])),
		col=adjustcolor('red',alpha=0.3), border=NA)


	points(prodMeanTable[,1],prodMeanTable[,colInds[i]],type='l',lwd=1.5)
	abline(h=0,lty=3)
	legend(1880,1,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
}


dev.off()


###############################################################################
# Compare the unc Estimte of the ensemble and Lenssen et al 2019
###############################################################################

prodGlobalCI <- read.csv('Data/totalCI_ERA.csv')

best <- read.table("Data/Best_uncertainty.txt",skip=48,header=FALSE)
bestInds <- which(best[,1] > 1879 & best[,1] < 2021)
had  <- read.csv("Data/HadCRUT.5.0.1.0.analysis.summary_series.global.annual.csv")
hadInds <- which(had[,1] > 1879 & had[,1] < 2021)

ensembleGlobalCI <- (globalAnnual[3,] - globalAnnual[1,])/2

pdf(sprintf('%s/uncComp.pdf',plotdir),10,5)
plot(prodGlobalCI$year, prodGlobalCI$ci95, type='l', lwd=2, ylim=c(0,0.22), col='red',
	ylab="95% Total Uncertainty (ºC)",xlab="Year",main="Comparison of Uncertainty Estimates")
points(best[bestInds,1],best[bestInds,3], type='l',col='blue',lwd=1, lty=2)
points(had[hadInds,1],(had[hadInds,4]-had[hadInds,3])/2, type='l', col='blue',lwd=1)
points(tYear, ensembleGlobalCI, type='l', lwd=3)

legend(1970,0.2, c('GISTEMP Ensemble', 'GISTEMP (Lenssen et al. 2019)', 'HadCRUT5 Infilled', 'Berkeley Earth'),
	lwd=c(3,2,1,1), col=c('black', 'red', 'blue', 'blue'), lty=c(1,1,1,2), bty='n')
dev.off()
###############################################################################
# Compare the field statistics
###############################################################################
load(sprintf('%s/Output/EnsembleStats/griddedSummaryStatistics.Rda',scratchDir))


yearVec <- c(1910, 1940, 1970)
timeInds <- c(which(timeMap[,1] == 1910 & timeMap[,2] == 1),
			  which(timeMap[,1] == 1940 & timeMap[,2] == 1),
			  which(timeMap[,1] == 1970 & timeMap[,2] == 1),
			  which(timeMap[,1] == 2010 & timeMap[,2] == 1))

pdf(sprintf('%s/fieldSdVis.pdf',plotdir),10,10)
mainVec <- sprintf('Temperature Uncertainty SD (January %s)',yearVec)


zCap <- 3.5
zr <- range(0,min(max(monthlyEnsembleSd[,,timeInds],na.rm=T),zCap))

pal <- designer.colors(256,brewer.pal(9,'YlOrBr'))
par(mfrow=c(3,2))
for(i in 1:3){

	plotMat <- monthlyEnsembleSd[,,timeInds[i]]
	plotMat[plotMat > zCap] <- zCap

	worldMap(lon,lat,plotMat,zlim=zr,col=pal,
		main=mainVec[i])
	hist(plotMat,xlim=zr, breaks=20,main='',
		xlab='Monthly Temperature SD')
}
dev.off()


###############################################################################
# Compare the sampling and homogenization uncertainty global series
###############################################################################

# load in the analyses
load(sprintf('%s/Output/LsatAnalysis/homogMeans.Rda',scratchDir))
load(sprintf('%s/Output/LsatAnalysis/samplingMeans.Rda',scratchDir))

# calculate the global mean uncertainty for the homg ensemble
homogAnnual <- apply(homogGlobalMean, 2, monthToYearMeans)
homogCI     <- apply(homogGlobalMean,1,sd) * 1.96

# calculate the global mean uncertainty for the sampling draws

samplingAnnual <- apply(samplingGlobalMean, c(1,3), monthToYearMeans)

samplingCI <- apply(samplingAnnual, 2, sd)* 1.96

plot
###############################################################################
# Compare the sampling and homogenization uncertainty fields (TODO?)
###############################################################################