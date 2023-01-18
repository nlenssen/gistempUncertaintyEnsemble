nCoreVec <- rev(c(1,2,5,10,15,20))

timingArray <- array(NA, dim=c(length(nCoreVec), 5, 3))

preserveVarVec <- c('step', 'nCoreVec', 'nCores', 'timingArray', 'preserveVarVec')


source('Namelists/awsSuperEnsemble.Rnl')
nCores <- 12
source('Code/05GenerateFullEnsemble_Empirical.R')
gc()


source('Namelists/awsSuperEnsemble.Rnl')
nCores <- 12
source('Code/06ProcessEnsemble.R')
gc()

source('Namelists/awsSuperEnsemble.Rnl')
nCores <- 12
source('Code/06ProcessEnsemblePart2.R')
gc()