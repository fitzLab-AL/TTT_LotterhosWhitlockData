
# Load libraries ---------------------------------------------------------------
#library(devtools)
#install_github("whitlock/OutFLANK")

library(LEA)
library(adegenet)
library(gdm)
library(gradientForest)
library(foreach)
library(doParallel)
library(OutFLANK)
library(pbapply)
library(gdata)
library(data.table)
################################################################################


# FUNCTIONS & SOURCED SCRIPTS --------------------------------------------------
#gfOutObj <- setClass("gfOutObj", slots = c(alFreq="data.frame", imp="list"))

source("/Users/mfitzpatrick/code/plantGenome/FstByRowtoGDMmatrices.R")

#####
# function to calculate pop-level allele counts / frequencies
# alleleTab = input allele data
# popSize = number of individuals sampled
# numPops = number of locations sampled = nrow(geo)
alleleDat <- function(alleleTab, popSize, numPops){
  mat <- matrix(NA, nrow=numPops, ncol=ncol(alleleTab))
  ind <- data.frame(i1=seq(1, popSize*numPops, by=popSize),
                    i2=seq(0, popSize*numPops, by=popSize)[-1])
  for(i in 1:numPops){
    mat[i,] <- apply(alleleTab[ind[i,1]:ind[i,2],], 2, sum)
  }
  colnames(mat) <- colnames(alleleTab)
  return(list(ifelse(1-(mat/popSize)<0.5, 1-(mat/popSize), mat/popSize), mat))}
#####

#####
# function to create a matrix for each locus of the counts of each 
# allele (columns) in each population (rows)
buildMats <- function(counts, popSize){
  bMat <- cbind(counts, popSize-counts)
  return(bMat)}
#####

#####
# function to calculate population pairwise Fst for a single locus
pwFst <- function(tab, numPops, ind1, ind2){
  mat <- matrix(0, numPops, numPops)
  for(i in 1:length(ind1)){
    mat[ind1[i], ind2[i]] <- WC_FST_FiniteSample_Haploids_2AllelesB_MCW(tab[c(ind1[i], ind2[i]),])[3]
    mat[ind1[i], ind2[i]] <- ifelse(mat[ind1[i], ind2[i]]<0, 0, mat[ind1[i], ind2[i]])
  }
  return(upperTriangle(mat, diag=F))}
#####

#####
# Function to add genetic distance to site-pair, remove NAs, and scale if desired.
# Scaling can improve model fitting in some instances by increasing the range 
# of Fst values.
finalPrepFst <-  function(x, sitePair, scale=F){
  fst <- sitePair 
  fst$distance <- x
  fst <- na.omit(fst)
  if(scale==T){
    return(scaleDist(fst)) # function sourced from above
  } else {return(fst)}}
#####

#####
# build output GF data frames
gfR2tab <- function(gfMods.list, alFreqs){
  i=1
  while(is.null(gfMods.list[[i]])){i=i+1}
  tab <- do.call(rbind, gfMods.list)
  vrNm <- rep(row.names(tab)[1:nrow(gfMods.list[[i]])], 
              nrow(tab)/nrow(gfMods.list[[i]]))
  tab <- data.frame(variable=vrNm, tab)
  tab <- dcast(tab, SNP~variable, value.var="imps")
  envR2 <- rowSums(tab[,-1])
  R2Tab <- data.frame(tab, envR2=envR2)
  
  # get name of SNP if it has a positive R2
  posR2 <- unlist(lapply(gfMods.list, function(x){
    return(as.character(unique(x[,2])))}))
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  negR2 <- !(colnames(alFreqs) %in% posR2)
  negR2 <- colnames(alFreqs)[negR2] 
  
  noGF <- data.frame(matrix(0, nrow=length(negR2), ncol=ncol(R2Tab)))
  colnames(noGF) <- colnames(R2Tab)           
  noGF$SNP <- negR2
  
  R2Tab <- rbind(R2Tab, noGF)
  snpID <- sapply(strsplit(as.character(R2Tab$SNP), "V"),function(x){as.numeric(x[2])})
  R2Tab$SNP <- snpID
  return(R2Tab[order(R2Tab$SNP),])}
#####
################################################################################


# Load and prep data (env, allelic) ---------------------------------------
# "background" environment
bgEnv <- read.table(paste(getwd(),"/results_AdaptreeEnviFor_R90.txt", sep="")) 

# sims
simFiles <- list.files(path=paste(getwd(), "/simfiles", sep=""), full.names=T)
simIDs <- unique(sapply(strsplit(sapply(strsplit(simFiles, "_NumPops"), function(x){
  x[1]}), "/simfiles/", fixed=T), function(x){
    x[2]}))

# x-y and environment
sim <- simFiles[grep(simIDs[9], simFiles)]
envSelect <- read.table(sim[grep("env", sim)])
names(envSelect) <- "envSelect"

# allele data
# columns = loci
# rows = total # of individuals (#populations x #inds sampled)
allelic <- fread(sim[grep("lfmm", sim)], header=F, data.table=F)

# data stats - used for indexing, etc
popSize <- nrow(allelic)/nrow(bgEnv)
numPops <- nrow(bgEnv)

# build env table 
envPop <- data.frame(popID=bgEnv[,1], x=bgEnv$X_Pops, y=bgEnv$Y_Pops, 
                     envSelect=unique(envSelect), bgEnv[,8:27])

# Calculate minor allele frequencies
alFreq.x <- alleleDat(allelic, popSize, numPops)
alFreq <- alFreq.x[[1]] # minor allele frequencies
alCount <- alFreq.x[[2]] # allele counts
rm(alFreq.x)

# returns a list n loci long, each element is a 2-column matrix with 
# major & minor allele counts
locusMats <- lapply(1:ncol(alCount), function(x, tab){
  buildMats(tab[,x], popSize)
}, tab=alCount)

# Minor allele frequencies using GF
# gfAllele.freq <- gradientForest(data.frame(envPop, alFreq), 
#                                 predictor.vars=colnames(envPop)[-c(1:3)],
#                                 response.vars=colnames(alFreq), ntree=500, 
#                                 trace=T)$result

# fit gf model to each SNP individually
cl <- makeCluster(12)
registerDoParallel(cl)

gfAllele.freq <- foreach(k=1:ncol(alFreq), .verbose=F, .packages=c("gradientForest")) %dopar%{
  locus <- data.frame(alFreq[,k])
  names(locus) <- colnames(alFreq)[k]
  gfLocus <- gradientForest(data.frame(envPop, locus), predictor.vars=colnames(envPop)[-c(1:3)], 
                           response.vars=colnames(alFreq)[k], 
                           corr.threshold=0.5, ntree=500, trace=F)
  if(!is.null(gfLocus)){
    imps <- importance(gfLocus)
    imps <- imps[order(names(imps))]
    data.frame(imps, SNP = colnames(alFreq)[k])
  }
}

stopCluster(cl)
#ttt <- gfOutObj(alFreq = data.frame(alFreq), imp = gfAllele.freq)

# find outliers
# table of importance values for each SNP (rows) and each var (columns)
gfAllele.R2 <- gfR2tab(gfAllele.freq, alFreq)
barplot(sort(apply(gfAllele.R2[,-c(1,23)],2, mean), decreasing=F), horiz=T)

write.csv(gfAllele.R2, paste("gfResults_", simIDs[9], ".csv", sep=""), row.names=F)


# Now, calculate pairwise Fst values between all populations
# Inefficient to do all possible pairs, so reduce to only unique
allPairs <- expand.grid(1:numPops, 1:numPops)
for(i in 1:nrow(allPairs)){
  if(allPairs[i,1] > allPairs[i,2]){allPairs[i,] <- allPairs[i,c(2,1)]}
}

# find unique combinations of population indices
uniqPairs <- unique(allPairs)
ind1 <- uniqPairs[,1]
ind2 <- uniqPairs[,2]

# calculate population pairwise Fst for all loci
Fst <- mclapply(locusMats, function(tab){pwFst(tab, numPops, ind1, ind2)},
                mc.cores=12, mc.cleanup = T)

# build site-pair tables needed for GDM
# easier to start with fake data then add Fst values
fstMat <- matrix(0, numPops, numPops)
fstMat <- data.frame(popID=envPop$popID, fstMat)
sitePair <- formatsitepair(bioData=fstMat, bioFormat=3, siteColumn="popID", 
                           predData=envPop, XColumn="x", 
                           YColumn="y")

# create complete site-pair tables for all Fst matrices
fstGDM <- pblapply(Fst, finalPrepFst, sitePair, scale=T)


# Fst using GDM
gdmMods <- mclapply(fstGDM, gdm, geo=F, mc.cores=12, mc.cleanup = T)





# x-y and environment
env <- read.table(sim[grep("env", sim)])

sims <- c("1R_R90_1351142986_950_9_NumPops=90_NumInd=20", 
          "2R_R90_1351142986_950_10_NumPops=90_NumInd=20",
          "IBD_R90_1351142986_950_11_NumPops=90_NumInd=20", 
          "IM_R90_1351142986_950_12_NumPops=90_NumInd=20")

for(k in 1:length(sims)){
  simsFiles <- list.files(pattern=sims[k])
  
  # x-y and environment
  env <- read.table(simsFiles[grep("env", simsFiles)])
  names(env) <- "env"
  geo <- read.table("/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/SchemeRandom1.txt")
  geo <- geo[which(geo$R90==TRUE),]
  
  
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- read.table(simsFiles[grep("lfmm", simsFiles)])
  
  # data stats - used for indexing, etc
  popSize <- nrow(allelic)/nrow(geo)
  numPops <- nrow(geo)
  
  # build env table w/ x-y coords
  xC <- geo$X_Pops[sort(rep(1:numPops, popSize))]
  yC <- geo$Y_Pops[sort(rep(1:numPops, popSize))]
  env.ind <- data.frame(x=xC, y=yC, env=env)
  envPop <- unique(env.ind)
  envPop <- data.frame(popID=geo$PopID, envPop)
  
  # function to calculate pop-level allele counts / frequencies
  # alleleTab = input allele data
  # popSize = number of individuals sampled
  # numPops = number of locations sampled = nrow(geo)
  alleleDat <- function(alleleTab, popSize, numPops){
    mat <- matrix(NA, nrow=numPops, ncol=ncol(alleleTab))
    ind <- data.frame(i1=seq(1, popSize*numPops, by=popSize),
                      i2=seq(0, popSize*numPops, by=popSize)[-1])
    for(i in 1:numPops){
      mat[i,] <- apply(alleleTab[ind[i,1]:ind[i,2],], 2, sum)
    }
    colnames(mat) <- colnames(alleleTab)
    return(list(ifelse(1-(mat/popSize)<0.5, 1-(mat/popSize), mat/popSize), mat))
  }
  
  # Minor allele frequencies
  alFreq.x <- alleleDat(allelic, popSize, numPops)
  alFreq <- alFreq.x[[1]]
  
  # Allele counts
  alCount <- alFreq.x[[2]]
  rm(alFreq.x)
  
  # For each locus, create a matrix of the counts of each allele (columns) in 
  # each population (rows)
  buildMats <- function(counts, popSize){
    bMat <- cbind(counts, popSize-counts)
    return(bMat)
  }
  
  # returns a list n loci long, each element is a 2-column matrix with 
  # major & minor allele counts
  locusMats <- pblapply(1:ncol(alCount), function(x, tab){
    buildMats(tab[,x], popSize)
  }, tab=alCount)
  
  # Now, calculate pairwise Fst values between all populations
  # Inefficient to do all possible pairs, so reduce to only unique
  allPairs <- expand.grid(1:numPops, 1:numPops)
  for(i in 1:nrow(allPairs)){
    if(allPairs[i,1] > allPairs[i,2]){allPairs[i,] <- allPairs[i,c(2,1)]}
  }
  
  # find unique combinations of population indices
  uniqPairs <- unique(allPairs)
  ind1 <- uniqPairs[,1]
  ind2 <- uniqPairs[,2]
  
  # function to calculate population pairwise Fst for a single locus
  pwFst <- function(tab, numPops, ind1, ind2){
    mat <- matrix(0, numPops, numPops)
    for(i in 1:length(ind1)){
      mat[ind1[i], ind2[i]] <- WC_FST_FiniteSample_Haploids_2AllelesB_MCW(tab[c(ind1[i], ind2[i]),])[3]
      mat[ind1[i], ind2[i]] <- ifelse(mat[ind1[i], ind2[i]]<0, 0, mat[ind1[i], ind2[i]])
    }
    return(upperTriangle(mat, diag=F))
  }
  
  # calculate population pairwise Fst for all loci
  Fst <- pblapply(locusMats, function(tab){pwFst(tab, numPops, ind1, ind2)})
  
  # build site-pair tables needed for GDM
  # easier to start with fake data then add Fst values
  fstMat <- matrix(0, numPops, numPops) 
  sitePair <- formatsitepair(bioData=fstMat, bioFormat=3, siteColumn="popID", 
                             predData=envPop, XColumn="x", 
                             YColumn="y")
  
  # loop to add genetic distance to site-pair, remove NAs, and scale if desired.
  # Scaling can improve model fitting in some instances by increasing the range 
  # of Fst values.
  finalPrepFst <-  function(x, sitePair, scale=F){
    fst <- sitePair 
    fst$distance <- x
    fst <- na.omit(fst)
    if(scale==T){
      return(scaleDist(fst)) # function sourced from above
    } else {return(fst)}
  }
  
  # create complete site-pair tables for all Fst matrices
  fstGDM <- pblapply(Fst, finalPrepFst, sitePair, scale=T)
  
  save.image(file=paste("dataPrepped4models_", sims[k], ".Rdata", sep=""))
}
################################################################################


# GF & GDM modeling ------------------------------------------------------------
sims <- c("1R_R90_1351142986_950_9_NumPops=90_NumInd=20", 
          "2R_R90_1351142986_950_10_NumPops=90_NumInd=20",
          "IBD_R90_1351142986_950_11_NumPops=90_NumInd=20", 
          "IM_R90_1351142986_950_12_NumPops=90_NumInd=20")

for(k in 1:length(sims)){
  prepDat <- list.files(pattern=sims[k])
  load(prepDat[grep("dataPrepped", prepDat)])
  
  # Allele presence-absence using GF
  # Need to change 0/1 to factor for classification (otherwise defaults to regression) 
  allelic <- apply(allelic, 2, factor)
  
  # fit individual GF models by locus.
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  gfAllele.pa <- foreach(k=1:ncol(allelic), .verbose=T, .packages=c("gradientForest")) %dopar%{
    locus <- data.frame(allelic[,k])
    names(locus) <- colnames(allelic)[k]
    gradientForest(data.frame(env.ind, locus), predictor.vars=c("x", "y", "env"),
                   response.vars=colnames(allelic)[k], ntree=500, trace=T)$result}
  
  stopCluster(cl)
  
  # Minor allele frequencies using GF
  gfAllele.freq <- gradientForest(cbind(unique(envPop), alFreq), predictor.vars=c("x", "y", "env"),
                                  response.vars=colnames(alFreq), ntree=500, trace=T)$result
  
  # Fst using GDM
  gdmMods <- pblapply(fstGDM, gdm, geo=T)
  
  save.image(file=paste("modelOutput_", sims[k], ".Rdata", sep=""))
}
################################################################################


# Prep & plot model output -----------------------------------------------------
sims <- c("1R_R90_1351142986_950_9_NumPops=90_NumInd=20", 
          "2R_R90_1351142986_950_10_NumPops=90_NumInd=20",
          "IBD_R90_1351142986_950_11_NumPops=90_NumInd=20", 
          "IM_R90_1351142986_950_12_NumPops=90_NumInd=20")

errorRates <- list()

# loop through each simulation and plot gf and gdm models
for(k in 1:length(sims)){
  modOut <- list.files(pattern=sims[k])
  load(modOut[grep("modelOutput", modOut)])
  
  ########## GF on allele P/A ##########
  # Allele presence-absence using GF
  allelePA <- unlist(gfAllele.pa)
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  negR2loci <- !(colnames(allelic) %in% names(allelePA))
  negR2loci <- colnames(allelic)[negR2loci] 
  noGF <- rep(0, length(negR2loci))
  names(noGF) <- negR2loci
  
  # R2 for all loci with GF model + 0's for loci with no GF model, sorted by name
  allelePA <- c(allelePA, noGF)
  
  # index to keep track of alleles under selection >9900
  ind <- NULL
  for(w in 1:length(names(allelePA))){
    ind[w] <- as.numeric(strsplit(names(allelePA)[w], "V")[[1]][2])
  }
  
  # data frame 
  allelePA <- data.frame(ind, R2=allelePA)
  allelePA <- allelePA[order(allelePA$ind),]
  
  # alpha to calculate false postives / false negatives
  alpha <- seq(0, 0.1, 0.001)
  quants <- quantile(allelePA$R2, probs=1-alpha)
  
  # false neg/pos rates at alpha = 0.01
  fnrPlot <- sum(allelePA$R2[allelePA$ind>9900] < quantile(allelePA$R2, probs=0.99))/sum(allelePA$ind>9900)
  fprPlot <- sum(allelePA$R2[allelePA$ind<=9900] >= quantile(allelePA$R2, probs=0.99))/sum(allelePA$ind<=9900)
  
  fnr <- NULL
  fpr <- NULL
  for(w in 1:length(quants)){
    fnr[w] <- sum(allelePA$R2[allelePA$ind>9900] < quants[w])/sum(allelePA$ind>9900)
    fpr[w] <- sum(allelePA$R2[allelePA$ind<=9900] >= quants[w])/sum(allelePA$ind<=9900)
  }
  
  errorRates[[k]] <- list()
  errorRates[[k]][[1]] <- data.frame(fnr, fpr)

  # assign colors for plotting
  bg <- (allelePA$R2 >= quantile(allelePA$R2, probs=0.99))+1
  bg <- ifelse(bg==1, rgb(0,0,0,0.15), rgb(1,0,0,0.5))
  bg[allelePA$ind>9900] <- rgb(0,0,1,0.5)
  
  # plot & save 
  filename <- paste("gf.allele.pa_", sims[k], ".tif", sep="")
  tiff(file=filename, width=18, height=12, units="in", compression="lzw", res=150,
       type="cairo")
  
  par(mai=c(1.25,1.25,1.25,1.25))
  
  col <- rgb(0,0,0,0.5)
  
  plot(allelePA$R2, bg=bg, pch=21, xlab="locus", ylab=expression(R^2), cex.lab=3.5, 
       cex.axis=2, col=col, ylim=c(0, 1), main="GF on allele p/a", 
       col.axis="grey50", cex.main=2)
  
  text(5000, 1, sims[k], cex=1.5)
  text(2000, 0.8, paste("false negative rate =", round(fnrPlot, 3), sep=" "), cex=3) 
  text(2000, 0.74, paste("false positive rate =", round(fprPlot, 4), sep=" "), cex=3)  
  dev.off()
  ####################
  
  ########## GF on allele freqs ##########
  gfAlfreq <- gfAllele.freq
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  negR2loci <- !(colnames(alFreq) %in% names(gfAlfreq))
  negR2loci <- colnames(alFreq)[negR2loci] 
  noGF <- rep(0, length(negR2loci))
  names(noGF) <- negR2loci
  
  # R2 for all loci with GF model + 0's for loci with no GF model, sorted by name
  gfAlfreq <- c(gfAlfreq, noGF)
  
  # index to keep track of alleles under selection >9900
  ind <- NULL
  for(w in 1:length(names(gfAlfreq))){
    ind[w] <- as.numeric(strsplit(names(gfAlfreq)[w], "V")[[1]][2])
  }
  
  # data frame
  gfAlfreq <- data.frame(ind, R2=gfAlfreq)
  gfAlfreq <- gfAlfreq[order(gfAlfreq$ind),]
  
  # alpha to calculate false postives / false negatives
  alpha <- seq(0, 0.1, 0.001)
  quants <- quantile(gfAlfreq$R2, probs=1-alpha)
  
  # false neg/pos rates at alpha = 0.01
  fnrPlot <- sum(gfAlfreq$R2[gfAlfreq$ind>9900] < quantile(gfAlfreq$R2, probs=0.99))/sum(gfAlfreq$ind>9900)
  fprPlot <- sum(gfAlfreq$R2[gfAlfreq$ind<=9900] >= quantile(gfAlfreq$R2, probs=0.99))/sum(gfAlfreq$ind<=9900)
  
  fnr <- NULL
  fpr <- NULL
  for(w in 1:length(quants)){
    fnr[w] <- sum(gfAlfreq$R2[gfAlfreq$ind>9900] < quants[w])/sum(gfAlfreq$ind>9900)
    fpr[w] <- sum(gfAlfreq$R2[gfAlfreq$ind<=9900] >= quants[w])/sum(gfAlfreq$ind<=9900)
  }
  
  errorRates[[k]][[2]] <- data.frame(fnr, fpr)
  
  # colors for plotting
  bg <- (gfAlfreq$R2 >= quantile(gfAlfreq$R2, probs=0.99))+1
  bg <- ifelse(bg==1, rgb(0,0,0,0.15), rgb(1,0,0,0.5))
  bg[gfAlfreq$ind>9900] <- rgb(0,0,1,0.5)
  
  # plot & save 
  filename <- paste("gf.allele.freq_", sims[k], ".tif", sep="")
  tiff(file=filename, width=18, height=12, units="in", compression="lzw", res=150,
       type="cairo")
  
  par(mai=c(1.25,1.25,1.25,1.25))
  
  col <- rgb(0,0,0,0.5)
  
  plot(gfAlfreq$R2, bg=bg, pch=21, xlab="locus", ylab=expression(R^2), cex.lab=3.5, 
       cex.axis=2, col=col, ylim=c(0, 1), main="GF on minor allele frequencies", 
       col.axis="grey50", cex.main=2)
  
  text(5000, 1, sims[k], cex=1.5)
  text(2000, 0.8, paste("false negative rate =", round(fnrPlot, 3), sep=" "), cex=3) 
  text(2000, 0.74, paste("false positive rate =", round(fprPlot, 4), sep=" "), cex=3)  
  dev.off()
  ####################
  
  
  ########## GDM on Fst ##########
  gdmRes <- NULL
  for(i in 1:length(gdmMods)){
    gdmRes[i] <- sum(gdmMods[[i]]$coefficients)
  }
  
  gdmRes <- ifelse(gdmRes<0,0,gdmRes)
  gdmRes <- data.frame(ind=1:length(gdmRes), sumCoeff=gdmRes)
  
  # alpha to calculate false postives / false negatives
  alpha <- seq(0, 0.1, 0.001)
  quants <- quantile(gdmRes$sumCoeff, probs=1-alpha)
  
  # false neg/pos rates at alpha = 0.01
  fnrPlot <- sum(gdmRes$sumCoeff[gdmRes$ind>9900] < quantile(gdmRes$sumCoeff, probs=0.99))/sum(gdmRes$ind>9900)
  fprPlot <- sum(gdmRes$sumCoeff[gdmRes$ind<=9900] >= quantile(gdmRes$sumCoeff, probs=0.99))/sum(gdmRes$ind<=9900)
  
  fnr <- NULL
  fpr <- NULL
  for(w in 1:length(quants)){
    fnr[w] <- sum(gdmRes$sumCoeff[gdmRes$ind>9900] < quants[w])/sum(gdmRes$ind>9900)
    fpr[w] <- sum(gdmRes$sumCoeff[gdmRes$ind<=9900] >= quants[w])/sum(gdmRes$ind<=9900)
  }
  
  errorRates[[k]][[3]] <- data.frame(fnr, fpr)
  
  # colors for plotting
  bg <- (gdmRes$sumCoeff >= quantile(gdmRes$sumCoeff, probs=0.99))+1
  bg <- ifelse(bg==1, rgb(0,0,0,0.15), rgb(1,0,0,0.5))
  bg[gdmRes$ind>9900] <- rgb(0,0,1,0.5)
  
  # plot & save 
  filename <- paste("gdm.fst_", sims[k], ".tif", sep="")
  tiff(file=filename, width=18, height=12, units="in", compression="lzw", res=150,
       type="cairo")
  
  par(mai=c(1.25,1.25,1.25,1.25))
  
  col <- rgb(0,0,0,0.5)
  
  plot(gdmRes$sumCoeff, bg=bg, pch=21, xlab="locus", ylab="Sum of coeff.", cex.lab=3.5, 
       cex.axis=2, col=col, ylim=c(0, 1), main="GDM on Fst", 
       col.axis="grey50", cex.main=2)
  
  text(5000, 1, sims[k], cex=1.5)
  text(2000, 0.8, paste("false negative rate =", round(fnrPlot, 3), sep=" "), cex=3) 
  text(2000, 0.74, paste("false positive rate =", round(fprPlot, 4), sep=" "), cex=3)  
  dev.off()
  ####################
}


filename <- "power.tif"
tiff(file=filename, width=12, height=12, units="in", compression="lzw", res=150,
     type="cairo")

par(mai=c(1.25,1.25,1.25,1.25))

col <- rgb(0,0,0,0.5)
plot(0,0, type="n", xlim=c(0,0.10), ylim=c(0,1), xlab="Alpha", ylab="Power", cex.lab=3.5, 
     cex.axis=2, main="Test", col=col, col.axis="grey50", cex.main=2)
for(i in 1:length(sims)){
  for(j in 1:3){
  lines(alpha, 1-errorRates[[i]][[j]]$fnr)
  }
}
dev.off()


filename <- "fpr.tif"
tiff(file=filename, width=12, height=12, units="in", compression="lzw", res=150,
     type="cairo")

par(mai=c(1.25,1.25,1.25,1.25))

col <- rgb(0,0,0,0.5)
plot(0,0, type="n", xlim=c(0,0.10), ylim=c(0,1), xlab="Alpha", ylab="False positive rate", 
     cex.lab=3.5, cex.axis=2, main="Test", col.axis="grey50", cex.main=2)
for(i in 1:length(sims)){
  for(j in 1:3){
    lines(alpha, errorRates[[i]][[j]]$fpr)
  }
}
dev.off()






# Sims for NSF proposal ---------------------------------------------------
sims <- list.files(path="/Users/mfitzpatrick/Rdatafiles/adaptiveSims/nemooutputCoreSet/IM_2patch_base_stnvar_withfreq/lfmm_files/",
                    pattern="10000", recursive=F)
inds <- list.files(path="/Users/mfitzpatrick/Rdatafiles/adaptiveSims/nemooutputCoreSet/IM_2patch_base_stnvar_withfreq/ind_info/",
                    pattern="10000", recursive=F)
sims <- sims[grep("1e-05", sims)]
sims <- sims[-grep("RData", sims)]
inds <- inds[grep("1e-05", inds)]

unlist(strsplit(sims, ".lfmm"))
unlist(strsplit(inds, ".indinfo"))

keep <- which(unlist(strsplit(inds, ".indinfo")) %in% unlist(strsplit(sims, ".lfmm")))


files <- list.files(path="/Users/mfitzpatrick/Rdatafiles/adaptiveSims/nemooutputCoreSet/IM_2patch_base_stnvar_withfreq/",
                    pattern="10000", recursive=T, full.names=T)
files <- files[grep("1e-05", files)]

#keep <- which(unlist(strsplit(inds, "ind_info")) %in% unlist(strsplit(sims, "lfmm_files")))

sims <- files[grep("lfmm", files)]
sims <- sims[-grep("RData", sims)]
inds <- files[grep("indinfo", files)][keep]

simsDat <- list()
for(i in 1:length(sims)){
  simsDat[[i]] <- data.frame(sims=sims[i], inds=inds[i])
}
  
gfSims <- function(simDat, nTree){
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- read.table(as.character(simDat$sims))
  indInfo <- read.table(as.character(simDat$inds), header=T)
  
  # removed fixed alleles
  notFixed <- which(unlist(lapply(apply(allelic, 2, unique), length))>1)
  allelic  <- allelic[,notFixed]
  
  # Allele presence-absence using GF
  # Need to change 0/1 to factor for classification (otherwise GF defaults to regression) 
  allelic <- apply(allelic, 2, factor)
  
  # predictors - need to add junk pred as GF crashes on one predictor...
  #preds <- data.frame(P1=indInfo$P1, pop=indInfo$pop)#junk=rep(0, length(indInfo$P1)))
  preds <- data.frame(P1=indInfo$P1, junk=rep(0, length(indInfo$P1)))
  
  gfMod <- gradientForest(data.frame(preds, allelic), 
                          predictor.vars=c("P1", "junk"),
                        response.vars=colnames(allelic), 
                        ntree=nTree, trace=T)
  
  gf.cumImp <- cumimp(gfMod, "P1", type="Species")
  
  out <- structure(list(R2=gfMod$result, Imp=gf.cumImp, vars=importance(gfMod)))
  
  outName <- paste(as.character(simDat$sims), ".RData", sep="")
  
  save(out, file=outName)
  
  #return(out)
}


gfOut <- mclapply(simsDat, gfSims, nTree=250, mc.cores=22)


nTree <- 250
for(k in 1:length(simsDat)){
  
  #outName <- paste(as.character(simsDat[[k]]$sims), ".RData", sep="")
  #outName <- strsplit(outName, "IM_2patch_base_stnvar_withfreq//")[[1]][2]
  #Rdatas <- list.files(path="/Users/mfitzpatrick/RdataFiles/adaptiveSims/nemooutputCoreSet/IM_2patch_base_stnvar_withfreq",
  #           pattern="RData", recursive=T, full.names=F)
  #if!(outName %in% Rdatas){
    
    
  
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- read.table(as.character(simsDat[[k]]$sims))
  indInfo <- read.table(as.character(simsDat[[k]]$inds), header=T)
  
  # removed fixed alleles
  notFixed <- which(unlist(lapply(apply(allelic, 2, unique), length))>1)
  allelic  <- allelic[,notFixed]
  
  # Allele presence-absence using GF
  # Need to change 0/1 to factor for classification (otherwise GF defaults to regression) 
  allelic <- apply(allelic, 2, factor)
  
  # predictors - need to add junk pred as GF crashes on one predictor...
  #preds <- data.frame(P1=indInfo$P1, pop=indInfo$pop)#junk=rep(0, length(indInfo$P1)))
  preds <- data.frame(P1=indInfo$P1, junk=rep(0, length(indInfo$P1)))
  
  gfMod <- gradientForest(data.frame(preds, allelic), 
                          predictor.vars=c("P1", "junk"),
                          response.vars=colnames(allelic), 
                          ntree=nTree, trace=T)
  
  gf.cumImp <- cumimp(gfMod, "P1", type="Species")
  
  out <- structure(list(R2=gfMod$result, Imp=gf.cumImp, vars=importance(gfMod)))
  outName <- paste(as.character(simsDat[[k]]$sims), ".RData", sep="")
  
  save(out, file=outName)
}


rDats <- list.files(path="/Users/mfitzpatrick/Rdatafiles/adaptiveSims/nemooutputCoreSet/IM_2patch_base_stnvar_withfreq/gfResults/RData",
                    pattern=".RData", recursive=T, full.names=T)

maketheOuts <- function(rDat){
  load(rDat)
  vvv <- paste("V", 1:9240, sep="")
  r2Out <- data.frame(loci=vvv, R2=0)
  
  #####################
  checkYNA <- function(lst, num){
    ###########
    #lst <- out[[2]]
    #num <- 1
    ###########
    
    ##returns TRUE if has NA, returns FALSE if does not have NAs
    if(sum(is.na(lst[[num]]$y))>0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  #####################
  
  #####################
  impRangeX <- function(lst, num){
    xR <- range(lst[[num]]$x)
    #yR <- range(lst[[num]]$y)
    return(xR)
  }
  #####################
  
  #####################
  impRangeY <- function(lst, num){
    #xR <- range(lst[[num]]$x)
    yR <- range(lst[[num]]$y)
    return(yR)
  }
  #####################
  
  thing <- lapply(1:length(out[[2]]), checkYNA, lst=out[[2]])
  
  R2 <- out[[1]][!unlist(thing)]
  cumImp <- out[[2]][!unlist(thing)]
  
  r2Out$R2[match(names(R2), r2Out$loci)] <- as.numeric(R2)
  
  namStub <- strsplit(rDat, ".lfmm.RData")[[1]][1]
  
  write.csv(r2Out, paste(namStub, "_R2.csv", sep=""), row.names=F)
  
  xRange <- lapply(1:length(cumImp), impRangeX, lst=cumImp)
  xmax <- max(unlist(xRange))
  xmin <- min(unlist(xRange))
  
  yRange <- lapply(1:length(cumImp), impRangeY, lst=cumImp)
  ymax <- max(unlist(yRange))
  ymin <- min(unlist(yRange))
  
  ylab <- "Cumulative importance"
  filename <- paste(namStub, ".pdf", sep="")
  pdf(file=filename, width=12, height=8)
  
  par(mai=c(1.25,1.25,1.25,1.25))
  plot(cumImp[[1]]$x, cumImp[[1]]$y, type="l", 
       xlab="", ylab=ylab, lwd=1, lty="solid", cex.lab=3.5, cex.axis=2, 
       col.axis="grey50", col=rgb(0,0,0,0.25), xlim=c(xmin,xmax),
       ylim=c(ymin,ymax))
  #mtext(xlabs[i], side=1, cex=3, line=4)
  for(j in 2:length(cumImp)){
    lines(cumImp[[j]]$x, cumImp[[j]]$y, lwd=1, lty="solid", 
          col=rgb(0,0,0,0.25))
  }
  dev.off()
}


sapply(rDats, maketheOuts)
