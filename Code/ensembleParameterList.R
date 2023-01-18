sstInds  <- 1:500
lsatInds <- 700:800

ensembleSize <- 101

ensembleTab <- matrix(NA, nrow=ensembleSize,ncol=2)
colnames(ensembleTab) <- c('lsat_ind','sst_ind')

ensembleTab[,2] <- sample(sstInds, ensembleSize)

ensembleTab[,1] <- rep(sample(lsatInds),5)[1:ensembleSize]

ensembleTab[,1] <- lsatInds

# check if any two rows are the same
badRows <- c()
for(i in 1:(ensembleSize/2)){
	ind1 <- 2*(i-1)+1
	ind2 <- ind1+1
	if(ensembleTab[ind1,1] == ensembleTab[ind2,1]){
		print('Found One!')
		badRows <- c(badRows,ind1)
		ensembleTab[c(ind1,ind2),1] <- sample(lsatInds,2)
	}
}


write.csv(ensembleTab,file=sprintf('Data/ensembleParameterList_%d.csv',ensembleSize))


# ensembleTab <- cbind(1:101,1)
# colnames(ensembleTab) <- c('lsat_ind','sst_ind')
# write.csv(ensembleTab,file=sprintf('Data/ensembleParameterList_LSAT_ENSEMBLE.csv'))