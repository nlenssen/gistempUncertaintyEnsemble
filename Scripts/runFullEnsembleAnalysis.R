###############################################################################
###############################################################################
# A script to run the entire GISTEMP Ensemble generation and analysis

# GISTEMP Uncertainty Ensemble
# Version 1.0.0 (August 21, 2024)
# Nathan Lenssen (lenssen@mines.edu)
# https://data.giss.nasa.gov/gistemp/
###############################################################################
###############################################################################

timingVec <- rep(NA, 11)

# objects to keep in the R instance between steps
preserveVarVec <- c('preserveVarVec', 'timingVec')

# Step 1
source('Namelists/awsNamelist.Rnl')
timingVec[1] <- system.time(source('Code/01ProcessERA.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# # Step 2
source('Namelists/awsNamelist.Rnl')
timingVec[2] <- system.time(source('Code/02ProcessGHCN.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# # Step 3
source('Namelists/awsNamelist.Rnl')
timingVec[3] <- system.time(source('Code/03InterpolateFields.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 4
source('Namelists/awsNamelist.Rnl')
print('Starting Step 4')
timingVec[4] <- system.time(source('Code/04SamplingUncertainty.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 5
source('Namelists/awsNamelist.Rnl')
print('Starting Step 5')
timingVec[5] <- system.time(source('Code/05GenerateFullEnsemble.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 6: Calculate key statistics of the ensemble

# Step 6a
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6a')
timingVec[6] <- system.time(source('Code/06ProcessEnsemble.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 6b
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6b')
timingVec[7] <- system.time(source('Code/06ProcessEnsemblePart2.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 7: Decompose the LSAT uncertainty into homogenization and sampling
# uncertainties by drawing simulations from the sampling unc. distribution

# Step 7a
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6a')
timingVec[8] <- system.time(source('Code/07lsatDataObjects.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 7b
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6a')
timingVec[9] <- system.time(source('Code/07lsatUncertainty.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 7c
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6a')
timingVec[10] <- system.time(source('Code/07lsatUncertaintyPart2.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 8
source('Namelists/awsNamelist.Rnl')
print('Starting Step 6a')
timingVec[10] <- system.time(source('Code/08PaperFigures.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)


print(cbind(c(1:6, 6.5, 7, 7.3, 7.6, 8),timingVec/60))

# save the timing vector to a file to review
save(timingVec, file='./finalTimingResults.Rda')


