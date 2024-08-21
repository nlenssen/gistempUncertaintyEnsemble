###############################################################################
###############################################################################
# Generate all figures in the Lenssen et al. 2024 paper

# GISTEMP Uncertainty Ensemble
# Version 1.0.0 (August 21, 2024)
# Nathan Lenssen (lenssen@mines.edu)
# https://data.giss.nasa.gov/gistemp/
###############################################################################
###############################################################################


# Source the local namelist for now
source('Namelists/awsNew_LocalTest.Rnl')


###############################################################################
# Get the filepaths and grid
###############################################################################
filePaths <- system(sprintf('ls %s/*.nc',ensembleOutDir),intern=T)

handle <- nc_open(filePaths[1])

lon <- ncvar_get(handle, 'lon')
lat <- ncvar_get(handle, 'lat')

nc_close(handle)


# and land mask
load(sprintf('%s/Intermediate/LandMasks/landMask_2x2.Rda',scratchDir))
landMask <- landMaskList$maximalMask


###############################################################################
# Plot the mean global series
###############################################################################

# get the true dat and trim to the same length as the record
prodMeanTable <- read.csv('Data/ZonAnn.Ts+dSST.csv')
prodVars <- names(prodMeanTable)

prodMeanTable <- prodMeanTable[which(prodMeanTable$Year <= endYear),]

# load the uncertainty ensemble mean series
load(sprintf('%s_AllLand/meanSeries_NEW.Rda',ensembleOutDir))

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

globalJan	 <- apply(ensembleGlobalMean[timeMap[,2]==1,,],1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
globalJul    <- apply(ensembleGlobalMean[timeMap[,2]==7,,],1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)


# get the calc for each month
globalMonthly <- array(NA, dim=c(3,nYear, 12))
globalMonthlyCI <- array(NA, dim=c(nYear,12))

for(i in 1:12){
	globalMonthly[,,i] <- apply(ensembleGlobalMean[timeMap[,2]==i,,],1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
	globalMonthlyCI[,i] <- (globalMonthly[3,,i]-globalMonthly[1,,i])/2
}



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
	xlab=xlabel, ylab= ylabel, main='Global Mean Anomaly (1950-1979 Climatology)',ylim=yr)

points(tYear, globalAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalAnnual[1,], rev(globalAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(Glob~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

legend('topleft',c('(a)'),bty='n')


#NH
plot(NHem~Year, data=prodMeanTable,type='l',
	xlab=xlabel, ylab= ylabel, main='Northern Hemisphere Mean Anomaly (1950-1979 Climatology)',ylim=yr)

points(tYear, nhAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(nhAnnual[1,], rev(nhAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(NHem~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

legend('topleft',c('(b)'),bty='n')

#SH
plot(SHem~Year, data=prodMeanTable,type='l',
	xlab=xlabel, ylab= ylabel, main='Southern Hemisphere Mean Anomaly (1950-1979 Climatology)',ylim=yr)

points(tYear, shAnnual[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(shAnnual[1,], rev(shAnnual[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(SHem~Year, data=prodMeanTable,type='l',lwd=1.5)
legend(1900,1.5,c('Production GISTEMP', 'GISTEMP Ensemble'), lwd=c(2,1), col=c('black', 'red'),bty='n')
abline(h=0,lty=3)

legend('topleft',c('(c)'),bty='n')

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

legend('topleft',c('(a)'),bty='n')
text(1950, 3, "1950-1979 Climatology")

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

text(1950, 1.6, "1950-1979 Climatology")


legend('topleft',c('(b)'),bty='n')

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

	legend('topleft',c(sprintf('(%s)',letters[i])),bty='n')

	text(1950, 1.6, "1950-1979 Climatology")
}


dev.off()


# quick calc for text
tInd <- c(121,141)

(bandsAnnual[3,tInd,] - bandsAnnual[1,tInd,])/2

###############################################################################
# Get monthly stuff
###############################################################################
# use era5 clim to calculate the base vals
load(sprintf('%s/Intermediate/ERA5/climData_ERA_2x2.Rda',scratchDir))
globalMeanTemps <- rep(NA,12)
cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

for(i in 1:12){
	globalMeanTemps[i] <- weighted.mean(eraClimList$clim[,,i],w=cosMat) - 273.15
}

globalSeasonalCycle <- globalMeanTemps - mean(globalMeanTemps)


# monthly

xlabel <- 'Year'
ylabel <- 'Monthly Mean Temperature Anomaly (°C)'
yr <- c(-0.65, 1.35)
yr2 <- pm(3)

pdf(sprintf('%s/Gistemp_Global_JanJul.pdf',plotdir),12,5)

set.panel(1,2)
# first anomalies
plot(NULL, xlim=range(timeMap[,1]),
	xlab=xlabel, ylab= ylabel, main='Monthly Global Mean Anomaly (1950-1979 Climatology)',ylim=yr)

points(tYear, globalJan[2,],col='blue',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalJan[1,], rev(globalJan[3,])),col=adjustcolor('blue',alpha=0.3), border=NA)

points(tYear, globalJul[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalJul[1,], rev(globalJul[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

legend(1900,1.5,c('January', 'July'), lwd=c(2,2), col=c('blue', 'red'),bty='n')

legend(1867, 1.5,c('(a)'),bty='n')


# the estimated raw values (using ERA5 clim) anomalies
plot(NULL, xlim=range(timeMap[,1]),
	xlab=xlabel, ylab= 'Monthly Anomaly Relative to Seasonal Cycle (°C)', main='Global Mean Anomaly + ERA5 Seasonal Cycle',ylim=yr2)

points(tYear, globalJan[2,] + globalSeasonalCycle[1],col='blue',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalJan[1,]+ globalSeasonalCycle[1], rev(globalJan[3,]+ globalSeasonalCycle[1])),col=adjustcolor('blue',alpha=0.3), border=NA)

points(tYear, globalJul[2,]+ globalSeasonalCycle[7],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalJul[1,]+ globalSeasonalCycle[7], rev(globalJul[3,]+ globalSeasonalCycle[7])),col=adjustcolor('red',alpha=0.3), border=NA)

legend(1900,3.45,c('January', 'July'), lwd=c(2,2), col=c('blue', 'red'),bty='n')


legend(1867, 3.45,c('(b)'),bty='n')


dev.off()



monthlyUncertaintyAverage <- rep(NA, 12)
for(i in 1:12){
	monthlyUncertaintyAverage[i] <- mean(globalMonthlyCI[which(tYear %in% 1960:2020),i])
}

pdf(sprintf('%s/Gistemp_Global_Monthly_Unc.pdf',plotdir),6,8)

set.panel(2,1)
# first anomalies
plot(NULL, xlim=range(timeMap[,1]),
	xlab='Year', ylab= 'Monthly 2-sigma Uncertainty', main='Monthly Uncertaity',ylim=range(globalMonthlyCI))


for(i in 1:12){
	points(tYear, rollmean(globalMonthlyCI[,i],10,na.pad=T),
		col=rainbow(12)[i],lwd=1.5,type='l')
}


legend(1950, 0.18, month.name[1:6], col=rainbow(12)[1:6],lwd=2, bty='n')
legend(1980, 0.18, month.name[7:12], col=rainbow(12)[7:12],lwd=2, bty='n')


plot(1:12,monthlyUncertaintyAverage, type='b', ylim=c(0,0.1),
	xlab='Month', ylab= 'Monthly 2-sigma Uncertainty', main='Monthly Mean (1960-2020) 2-sigma Uncertainty')


dev.off()

###############################################################################
# Compare the unc Estimte of the monthly estimates
###############################################################################
pal <- brewer.pal(8,'Dark2')

best <- read.table("Data/Best_uncertainty_Monthly.txt",skip=86,header=FALSE)
had  <- read.csv("Data/HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")
noaa<- read.csv("Data/globalMonthlyCI_NOAA_1850_2015.csv")


# For January
monthInd <- 1

bestInds <- which(best[,1] > 1879 & best[,1] < 2021 & best[,2] == monthInd)
hadInds <- which(had[,1] > 1880 & had[,1] < 2021 & as.numeric(substr(had[,1],6,7)) == monthInd)
noaaInds <- which(noaa[,1] > 1879 & noaa[,1] < 2016 & noaa[,2] == monthInd)



pdf(sprintf('%s/uncComp_JanJuly.pdf',plotdir),12,5)
set.panel(1,2)
plot(NULL, ylim=c(0,0.25), xlim=c(1880,2020),
	ylab="95% Global January Uncertainty (ºC)",xlab="Year",
	main="Comparison of January Global Uncertainty")

points(best[bestInds,1],best[bestInds,4], type='l',col=pal[2],lwd=2, lty=1)
points(best[bestInds,1],(had[hadInds,4]-had[hadInds,3])/2, type='l', col=pal[3],lwd=2)
points(noaa[noaaInds,1],(noaa[noaaInds,5]-noaa[noaaInds,3])/2, type='l', col=pal[4],lwd=2)

points(tYear, (globalJan[3,] - globalJan[1,])/2, type='l', lwd=3)

legend(1950,0.25, c('GISTEMP Ensemble', 'Berkeley Earth', 'HadCRUT5 Infilled', 'NOAA GlobalTemp'),
	lwd=c(3,2,2,2), col=c('black', pal[2:4]), lty=1, bty='n')


legend(1867, 0.265,c('(a)'),bty='n')

# For July
monthInd <- 7

bestInds <- which(best[,1] > 1879 & best[,1] < 2021 & best[,2] == monthInd)
hadInds <- which(had[,1] > 1880 & had[,1] < 2021 & as.numeric(substr(had[,1],6,7)) == monthInd)
noaaInds <- which(noaa[,1] > 1879 & noaa[,1] < 2016 & noaa[,2] == monthInd)




plot(NULL, ylim=c(0,0.25), xlim=c(1880,2020),
	ylab="95% Global July Uncertainty (ºC)",xlab="Year",
	main="Comparison of July Global Uncertainty")

points(best[bestInds,1],best[bestInds,4], type='l',col=pal[2],lwd=2, lty=1)
points(best[bestInds,1],(had[hadInds,4]-had[hadInds,3])/2, type='l', col=pal[3],lwd=2)
points(noaa[noaaInds,1],(noaa[noaaInds,5]-noaa[noaaInds,3])/2, type='l', col=pal[4],lwd=2)

points(tYear, (globalJul[3,] - globalJul[1,])/2, type='l', lwd=3)


legend(1867, 0.265,c('(b)'),bty='n')


dev.off()



###############################################################################
# Compare the unc Estimte of the ensemble and Lenssen et al 2019
###############################################################################

pal <- brewer.pal(8,'Dark2')

prodGlobalCI <- read.csv('Data/totalCI_ERA.csv')

best <- read.table("Data/Best_uncertainty.txt",skip=48,header=FALSE)
bestInds <- which(best[,1] > 1879 & best[,1] < 2021)
had  <- read.csv("Data/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")
hadInds <- which(had[,1] > 1879 & had[,1] < 2021)
noaa<- read.csv("Data/globalAnnualCI_NOAA_1850_2015.csv")
noaaInds <- which(had[,1] > 1879 & had[,1] < 2021)

dkrz<- read.csv("Data/globalAnnualCI_DKRZ_1850_2022.csv")
dkrzInds <- which(had[,1] > 1879 & had[,1] < 2021)

ensembleGlobalCI <- (globalAnnual[3,] - globalAnnual[1,])/2

pdf(sprintf('%s/uncComp.pdf',plotdir),10,5)
plot(prodGlobalCI$year, prodGlobalCI$ci95, type='l', lwd=2, ylim=c(0,0.22), col=pal[1],
	ylab="95% Total Uncertainty (ºC)",xlab="Year",main="Comparison of Global Annual Uncertainty Estimates")
points(best[bestInds,1],best[bestInds,3], type='l',col=pal[2],lwd=2, lty=1)
points(had[hadInds,1],(had[hadInds,4]-had[hadInds,3])/2, type='l', col=pal[3],lwd=2)
points(noaa[noaaInds,1],(noaa[noaaInds,4]-noaa[noaaInds,2])/2, type='l', col=pal[4],lwd=2)
points(dkrz[dkrzInds,1],(dkrz[dkrzInds,4]-dkrz[dkrzInds,2])/2, type='l', col=pal[5],lwd=2)

points(tYear, ensembleGlobalCI, type='l', lwd=3)

legend(1970,0.2, c('GISTEMP Ensemble', 'GISTEMP (Lenssen et al. 2019)', 'Berkeley Earth', 'HadCRUT5 Infilled', 'NOAA GlobalTemp', 'DKRZ CNN'),
	lwd=c(3,2,2,2,2,2), col=c('black', pal[1:5]), lty=c(1,1,1,1,1,1), bty='n')
dev.off()
###############################################################################
# Compare the field statistics
###############################################################################
load(sprintf('%s/griddedSummaryStatistics.Rda',ensembleOutDir))


yearVec <- c(1910, 1940, 1970)
timeInds <- c(which(timeMap[,1] == 1910 & timeMap[,2] == 1),
			  which(timeMap[,1] == 1940 & timeMap[,2] == 1),
			  which(timeMap[,1] == 1970 & timeMap[,2] == 1),
			  which(timeMap[,1] == 2010 & timeMap[,2] == 1))

cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

pdf(sprintf('%s/fieldSdVis.pdf',plotdir),10,10)
mainVec <- sprintf('Temperature Uncertainty SD (January %s)',yearVec)

zCap <- 3.5
zr <- range(0,min(max(monthlyEnsembleSd[,,timeInds],na.rm=T),zCap))


zBreaks <- seq(0,zCap, by=0.5)
pal <- brewer.pal(length(zBreaks)-1, 'YlOrBr')
par(mfrow=c(3,2))
for(i in 1:3){

	plotMat <- monthlyEnsembleSd[,,timeInds[i]]
	plotMat[plotMat > zCap] <- zCap

	robinsonPlot(plotMat,range=c(0,zCap), main=mainVec[i], col=pal)

	text(-1.35e7,8.5e6,sprintf('(%s)',letters[(i-1)*2 + 1]))

	goodInds <- which(!is.na(plotMat))

	wtd.hist(plotMat[goodInds],weight = cosMat[goodInds], ylim=c(0,5500), xlim=zr, breaks=15,main='',
		xlab='Monthly Temperature SD',ylab='Area-weighted Frequency', col = "lightgray")

	legend(-0.3,5700,sprintf('(%s)',letters[(i)*2]),bty='n')
}

dev.off()


###############################################################################
# Compare the sampling and homogenization uncertainty global series
###############################################################################


load(sprintf('%s/Output/LsatAnalysis/samplingMeans.Rda',scratchDir))

# get the annual mean empirical CIs
globalAnnualSampling <- apply(apply(samplingGlobalMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
nhAnnualSampling     <- apply(apply(samplingNHemMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
shAnnualSampling     <- apply(apply(samplingSHemMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)


load(sprintf('%s/Output/LsatAnalysis/homogMeans.Rda',scratchDir))

# get the annual mean empirical CIs
globalAnnualHomog <- apply(apply(homogGlobalMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
nhAnnualHomog     <- apply(apply(homogNHemMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)
shAnnualHomog     <- apply(apply(homogSHemMean, 2, monthToYearMeans),1, quantile, probs=c(0.025, 0.5, 0.975),na.rm=T)



xlabel <- 'Year'
ylabel <- 'Annual Mean Temperature Anomaly (°C)'
yr <- c(-.25, .25)
yr2 <- c(0,0.75)

pdf(sprintf('%s/Gistemp_GlobalHemi_Annual_Sampling_Unc.pdf',plotdir),6,9)
set.panel(2,1)
# plot the sampling unc series for all three
plot(tYear, shAnnualSampling[2,],col='blue',lwd=1.5,type='l',ylim=yr)
polygon(c(tYear, rev(tYear)), c(shAnnualSampling[1,], rev(shAnnualSampling[3,])),col=adjustcolor('blue',alpha=0.3), border=NA)

points(tYear, nhAnnualSampling[2,],col='red',lwd=1.5,type='l',ylim=yr)
polygon(c(tYear, rev(tYear)), c(nhAnnualSampling[1,], rev(nhAnnualSampling[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(tYear, globalAnnualSampling[2,],col='black',lwd=1.5,type='l',ylim=yr)
polygon(c(tYear, rev(tYear)), c(globalAnnualSampling[1,], rev(globalAnnualSampling[3,])),col=adjustcolor('black',alpha=0.3), border=NA)

# plot the 95% magnitude for them
plot(tYear, (globalAnnualSampling[3,] - globalAnnualSampling[1,])/2, type='l', lwd=2, col='black',ylim=yr2)
points(tYear, (nhAnnualSampling[3,] - nhAnnualSampling[1,])/2, type='l', lwd=2, col='red')
points(tYear, (shAnnualSampling[3,] - shAnnualSampling[1,])/2, type='l', lwd=2, col='blue')
dev.off()


xlabel <- 'Year'
ylabel <- 'Annual Mean Temperature Anomaly (°C)'
yr3 <- c(-1, 2)
yr4 <- c(0,0.25)

pdf(sprintf('%s/Gistemp_GlobalHemi_Annual_Homog_Unc.pdf',plotdir),6,9)
set.panel(2,1)
# plot the sampling unc series for all three
plot(tYear, shAnnualHomog[2,],col='blue',lwd=1.5,type='l',ylim=yr3)
polygon(c(tYear, rev(tYear)), c(shAnnualHomog[1,], rev(shAnnualHomog[3,])),col=adjustcolor('blue',alpha=0.3), border=NA)

points(tYear, nhAnnualHomog[2,],col='red',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(nhAnnualHomog[1,], rev(nhAnnualHomog[3,])),col=adjustcolor('red',alpha=0.3), border=NA)

points(tYear, globalAnnualHomog[2,],col='black',lwd=1.5,type='l')
polygon(c(tYear, rev(tYear)), c(globalAnnualHomog[1,], rev(globalAnnualHomog[3,])),col=adjustcolor('black',alpha=0.3), border=NA)

# plot the 95% magnitude for them
plot(tYear, (globalAnnualHomog[3,]- globalAnnualHomog[1,])/2, type='l', lwd=2, col='black',ylim=yr4)
points(tYear, (nhAnnualHomog[3,]- nhAnnualHomog[1,])/2, type='l', lwd=2, col='red')
points(tYear, (shAnnualHomog[3,]- shAnnualHomog[1,])/2, type='l', lwd=2, col='blue')
dev.off()



load('Data/fig3bList.Rda')
samplingOldTrim <- fig3bList$sampling[-c(1:30,168)]
pdf(sprintf('%s/Gistemp_GlobalHemi_Annual_Comp.pdf',plotdir),10,5)
par(mfrow=c(1,2))
plot(fig3bList$t, samplingOldTrim, ylim=yr2, type='l', lwd=1, lty=1,xlim=c(1880,2020),
	xlab='Year', ylab='95% LSAT Confidence Interval (°C)', main = 'LSAT Uncertainty Decomposition',col='blue')
grid()
abline(h=0)
points(fig3bList$t,fig3bList$total,type='l',col='black',lwd=1,lty=1)
points(fig3bList$t,fig3bList$homog,type='l',col='red',lwd=1,lty=1)


samplingNew <- (globalAnnualSampling[3,] - globalAnnualSampling[1,])/2
homogNew    <- (globalAnnualHomog[3,]- globalAnnualHomog[1,])/2
totalNew    <- 1.96*sqrt((samplingNew/1.96)^2 + (homogNew/1.96)^2)

points(tYear, samplingNew, type='l', lwd=2, col='blue')
points(tYear, homogNew, type='l', lwd=2, col='red')
points(tYear, totalNew, type='l', lwd=2, col='black')

legend('topright', 
	c('Total LSAT Uncertainty', 'Homogenization Uncertainty', 'Sampling Uncertainty'),
	col=c('black', 'red', 'blue'), lwd=2,lty=c(1,1,1),bty='n')

text(1880,0.75,sprintf('(a)'))



plot(NULL, xlim=c(1880,2020),ylim=pm(0.25),
	xlab='Year', ylab='Difference of 95% CI (°C)', main = 'Difference from Lenssen et al. (2019)')
grid()
abline(h=0)
points(tYear,totalNew - c(fig3bList$total,rep(NA,4)),type='l',col='black',lwd=1.5,lty=1)
points(tYear,homogNew - c(fig3bList$homog,rep(NA,4)),type='l',col='red',lwd=1.5,lty=1)
points(tYear,samplingNew - c(samplingOldTrim,rep(NA,4)),type='l',col='blue',lwd=1.5,lty=1)


legend('topright', 
	c('Total LSAT Uncertainty', 'Homogenization Uncertainty', 'Sampling Uncertainty'),
	col=c('black', 'red', 'blue'), lwd=2,lty=c(1,1,1),bty='n')

text(1880,0.25,sprintf('(b)'))

dev.off()



###############################################################################
# Compare the sampling and homogenization uncertainty fields
###############################################################################
load(sprintf('%s/Output/LsatAnalysis/samplingGriddedLsatAnalysis.Rda',scratchDir))
load(sprintf('%s/Output/LsatAnalysis/homogGriddedLsatAnalysis.Rda',scratchDir))

# check at a single time point
pal <- designer.colors(128,brewer.pal(11,'Spectral'))

pdf(sprintf('%s/LSAT_ratio_map.pdf',plotdir),14,5)
par(mfrow=c(1,2))
timeInd <- 1681

plotMat <- log((samplingEnsembleSd[,,timeInd])^2/(homogEnsembleSd[,,timeInd])^2)*landMask
zMax <- max(abs(plotMat),na.rm=T)
zr <- c(-zMax,zMax)
robinsonPlot(plotMat,range=zr, main='Log ratio of Sampling Unc/Homog Unc (January 2020)', col=pal)

text(-1.35e7,8.5e6,sprintf('(a)'))



timeInd <- 1687

plotMat <- log((samplingEnsembleSd[,,timeInd])^2/(homogEnsembleSd[,,timeInd])^2)*landMask
zMax <- max(abs(plotMat),na.rm=T)
zr <- c(-zMax,zMax)
robinsonPlot(plotMat,range=zr, main='Log ratio of Sampling Unc/Homog Unc (July 2020)', col=pal)


text(-1.35e7,8.5e6,sprintf('(b)'))


dev.off()

