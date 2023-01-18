nCoreVec <- c(30, 20, 10, 5)

timingArray <- array(NA, dim=c(length(nCoreVec), 5))

preserveVarVec <- c('step', 'nCoreVec', 'timingArray', 'preserveVarVec')


for(step in 1:length(nCoreVec)){
	print(paste('Loop', step, 'of', length(nCoreVec)))
	

	source('Namelists/awsFullEnsemble.Rnl')
	nCores_Step3 <- nCoreVec[step]
	timingArray[step,] <- system.time(source('Code/03InterpolateFields.R'))
	
	objs <- ls()
	objs <- objs[-which(objs %in% preserveVarVec)]
	rm(list=objs)

}

#save(nCoreVec, timingArray, file='Data/timingResultsTimeSeries.Rda')


