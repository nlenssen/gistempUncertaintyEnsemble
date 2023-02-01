# Final timing test on AWS run on Nov 29th

timingVec <- rep(NA, 7)

preserveVarVec <- c('preserveVarVec', 'timingVec')


# Step 1
# source('Namelists/awsFullEnsemble.Rnl')
# timingVec[1] <- system.time(source('Code/01ProcessERA.R'))[3]

# objs <- ls()
# objs <- objs[-which(objs %in% preserveVarVec)]
# rm(list=objs)

# # Step 2
# source('Namelists/awsFullEnsemble.Rnl')
# timingVec[2] <- system.time(source('Code/02ProcessGHCN.R'))[3]

# objs <- ls()
# objs <- objs[-which(objs %in% preserveVarVec)]
# rm(list=objs)

# # Step 3
# source('Namelists/awsFullEnsemble.Rnl')
# timingVec[3] <- system.time(source('Code/03InterpolateFields_NEW.R'))[3]

# objs <- ls()
# objs <- objs[-which(objs %in% preserveVarVec)]
# rm(list=objs)

# Step 4
source('Namelists/awsFullEnsemble.Rnl')
timingVec[4] <- system.time(source('Code/04SpatialUncertaintyCoarse_Empirical.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 5
source('Namelists/awsFullEnsemble.Rnl')
timingVec[5] <- system.time(source('Code/05GenerateFullEnsemble_Empirical.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 6a
source('Namelists/awsFullEnsemble.Rnl')
timingVec[6] <- system.time(source('Code/06ProcessEnsemble.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)


# Step 6b
source('Namelists/awsFullEnsemble.Rnl')
timingVec[7] <- system.time(source('Code/06ProcessEnsemblePart2_NEW.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)




save(timingVec, file='Data/finalTimingResults.Rda')


