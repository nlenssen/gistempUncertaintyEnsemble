rnCoreVec <- rev(c(1,2,5,10,15,20))

timingArray <- array(NA, dim=c(length(nCoreVec), 5, 3))

preserveVarVec <- c('step', 'nCoreVec', 'nCores', 'timingArray', 'preserveVarVec')


for(step in 1:length(nCoreVec)){
	print(paste('Loop', step, 'of', length(nCoreVec)))
	nCores <- nCoreVec[step]

	source('Namelists/awsFullEnsemble.Rnl')
	timingArray[step,,1] <- system.time(source('Code/05GenerateFullEnsemble_Empirical.R'))

	objs <- ls()
	objs <- objs[-which(objs %in% preserveVarVec)]
	rm(list=objs)

	source('Namelists/awsFullEnsemble.Rnl')
	timingArray[step,,2] <- system.time(source('Code/06ProcessEnsemble.R'))
	
	objs <- ls()
	objs <- objs[-which(objs %in% preserveVarVec)]
	rm(list=objs)


	source('Namelists/awsFullEnsemble.Rnl')
	timingArray[step,,3] <- system.time(source('Code/06ProcessEnsemblePart2.R'))
	
	objs <- ls()
	objs <- objs[-which(objs %in% preserveVarVec)]
	rm(list=objs)

}


save(nCoreVec, timingArray, file='Data/timingResults.Rda')


