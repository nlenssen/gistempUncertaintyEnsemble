# Final timing test on AWS run on Nov 29th

timingVec <- rep(NA, 7)

preserveVarVec <- c('preserveVarVec', 'timingVec')


# Step 1
source('Namelists/awsNew.Rnl')
timingVec[1] <- system.time(source('CodeRewrite/01ProcessERA.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 2
source('Namelists/awsNew.Rnl')
timingVec[2] <- system.time(source('CodeRewrite/02ProcessGHCN.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 3
source('Namelists/awsNew.Rnl')
timingVec[3] <- system.time(source('CodeRewrite/03InterpolateFields.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 4
source('Namelists/awsNew.Rnl')
timingVec[4] <- system.time(source('CodeRewrite/04SamplingUncertainty.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 5
source('Namelists/awsNew.Rnl')
timingVec[5] <- system.time(source('CodeRewrite/05GenerateFullEnsemble.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)

# Step 6a
source('Namelists/awsNew.Rnl')
timingVec[6] <- system.time(source('CodeRewrite/06ProcessEnsemble.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)


# Step 6b
source('Namelists/awsNew.Rnl')
timingVec[7] <- system.time(source('CodeRewrite/06ProcessEnsemblePart2.R'))[3]

objs <- ls()
objs <- objs[-which(objs %in% preserveVarVec)]
rm(list=objs)



print(cbind(1:7,timingVec/60))

save(timingVec, file='Data/finalTimingResults_Rewrite.Rda')


