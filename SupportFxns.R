#######################################################
# SAMPLE CODE FOR DYNAMIC SURVIVAL PREDICTION
#######################################################
####################################################
## ----------------------------------------------------------------------------
##  Code to accompany "Dynamic Survival Prediction Combining Landmarking 
##                     with a Machine Learning Ensemble: 
##                     Methodology and Empirical Comparison"
##
##  Authors: Kamaryn Tanner, Linda Sharples, Rhian Daniel, Ruth Keogh
##  Contact for code: kamaryn.tanner1@lshtm.ac.uk
## ----------------------------------------------------------------------------
##  Date        :       15-May-2020
##  R Version   :       3.5.3                                                           
##  
##  This file contains the helper functions called from the main file           
##  
## ------------------------------------------------------------------
#####################################################################


##
## --- Function to prepare the pbc dataset in the format we need for dynamic survival prediction
##
## pbcDat will be a data.frame in counting process format
## start, stop represent the start and stop times that apply to logBili and the event indicator
## fold designates the fold number for each id in the cross-validation
##
prepPBC <- function(){
## We use the sequential data from the PBC dataset {survival} 
## to illustrate the methods.
##  Longitudinal predictor:   bili
##  Time-dependent predictor: spiders
##  Time-fixed predictors:    age, sex 
##  
  pbcDat <- pbcseq[c("id", "status", "futime", "day", "age", "sex", "bili", "trt", "spiders")]
  ## Create composite outcome
  pbcDat$statusComp <- ifelse(pbcDat$status==0, 0, 1)
  ## Work with log(serum bilirubin)
  pbcDat$logBili <- log(pbcDat$bili)
  ##females=1
  pbcDat$sex <- if_else(pbcDat$sex=="f", 1, 0)
  ## Convert timescale to years
  pbcDat$measTime <- pbcDat$day/365.25
  pbcDat$survTime <- pbcDat$futime/365.25
  ## Format with start/stop times 
  ##  (this code from JMbayes package Examples; find by typing ?jointModelBayes)
  ##
  pbcDat$start <- pbcDat$measTime
  splitID <- split(pbcDat[c("start", "survTime")], pbcDat$id)
  pbcDat$stop <- unlist(lapply(splitID,
                             function (d) c(d$start[-1], d$survTime[1]) ))
  pbcDat$event <- with(pbcDat, ave(statusComp, id,
                                 FUN = function (x) c(rep(0, length(x)-1), x[1])))
  pbcDat <- pbcDat[!is.na(pbcDat$spiders), ]
  pbcDat <- subset(pbcDat, select= -c(futime, day, status, bili))

  ## Assign data (by id) to one of n folds for cross-validation
  numFolds <- 10
  ids <- unique(pbcDat$id)
  set.seed(1212)
  #Randomly shuffle the data
  ids<-ids[sample(length(ids))]
  #Create 10 equal size folds
  folds <- cut(seq(1,length(ids)),breaks=numFolds,labels=FALSE)
  idFold <- data.frame(id=ids, fold=folds)
  # Join this to dat.train
  pbcDat <- join(pbcDat, idFold, by="id", type="left")
  return(pbcDat)
}


##
## --- Function to fit a joint model and obtain dynamic survival predictions
##

runModelJoint <- function(testDat, trainDat, w, tseq ){
  ##NOTE: Many interesting examples of joint models fit to PBC data can be found in
  ## "Joint Models for Longitudinal and Time-to-Event Data With Applications in R"
  ## by Dmitris Rizopoulos, CRC Press, 2012.
  
  # longitudinal sub-model
  lmeFit <- lme(logBili ~ trt + measTime,
                random = ~measTime | id, data = pbcDat)
  # survival sub-model
  coxFit <- coxph(Surv(start, stop, event) ~ trt + spiders + cluster(id),
                  data = pbcDat, x = TRUE, model = TRUE)
  # joint model
  jointFit <- jointModel(lmeFit, coxFit, timeVar = "measTime", method="spline-PH-aGH")
  
  ##Get predictions at same grid of landmark times used in other methods
  maxTime <- 15
  testIds <- unique(testDat$id)
  ##set up storage for predicted survival probabilities
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds
  
  for( i in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    ## limit to ids that were at risk and only predict using data before landmark time
    sData <- testDat[testDat$survTime>tseq[i],]
    sData <- sData[sData$measTime <= tseq[i],]
    sIds <- as.character(unique(sData$id))
    predProbs <- matrix(nrow=3, ncol=length(sIds))
    colnames(predProbs) <- sIds
    #This loop calculates the conditional survival probabilities by id
    for( l in 1:length(sIds)){
      ND <- sData[sData$id %in% sIds[l],]
      ND$LM <- tseq[i]
      p1 <- JM::survfitJM(jointFit, newdata=ND, survTimes=c(tseq[i], (tseq[i]+w)), last.time = "LM")
      ##capture mean survival probability (column 2)
      sp <- p1[[1]][[1]][,2]
      if(p1[[1]][[1]][1,1] != tseq[[i]]){ #then we need to add a 1 to the front of sp
        sp <- c(1.0, sp)
      }
      testSurvPreds[i, which(testIds==sIds[l])] <- t(sp)[2]
    }
  }
  ##Predicted survival probabilities on the test set for this CV fold
  return(testSurvPreds)  
}



##
## --- Function to fit Cox landmarking and obtain dynamic survival predictions
##

runModelCox_LM <- function(testDat, trainDat, w, tseq){
  ## get stacked data for landmarking
  lmData <- genStackedData(dat=subset(trainDat, select=-c(stop, event, fold)), w, tseq)
  
  ##run Cox regression stratified on landmark time
  LMcox <- coxph(Surv(lmTime, survTime, event)~ trt + spiders + logBili +
                   strata(lmTime) + cluster(id), data=lmData, method="breslow")
  
  ##Get predictions along landmark time grid
  testIds <- unique(testDat$id)
  ##set up storage for predicted survival probabilities
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds
  
  for( i in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    sData <- testDat[testDat$survTime>tseq[i],]
    sData <- sData[sData$measTime <= tseq[i],]
    sData$survTime[sData$survTime > (tseq[i]+w)] <- (tseq[i]+w)
    sData <- ddply(sData,.(id), tail,1)  # Keep only the most recent record for each id
    newData <- subset(sData, select= -c(id, measTime, survTime, statusComp, event, fold))
    
    ## Predicted curves from a coxph model have one row for each stratum in the Cox model fit 
    sfhor <- survfit(LMcox,newdata=newData)
    ##The data I want will be between indices:
    if(i == 1){
      index1 <- 1  
    }else {
      index1 <- cumsum(sfhor$strata)[i-1] + 1
    }
    index2 <- index1 + sfhor$strata[i] - 1
    
    nhor <- index2
    ids <- unique(sData$id)
    if( length(ids)>1){
      for(j in 1:length(ids)){#loop over individuals
        testSurvPreds[i, which(testIds==ids[j])] <- sfhor$surv[nhor,j]
      }
    }else {
      testSurvPreds[i,ids[1]] <- sfhor$surv[nhor]
    }
  }
  ##Predicted survival probabilities on the test set for this CV fold
  return(testSurvPreds)    
}


##
## --- Function to fit Super Learner landmarking and obtain dynamic survival predictions
##

runModelSL_LM <- function(testDat, trainDat, w, tseq, SL.lib){
  
  ## get stacked data for landmarking with SL
  sd <- genStackedSLData(dat=subset(trainDat, select=-c(stop, event, fold)), w, tseq)
  stackedData.disc <- sd
#  newXstackedData.disc <- sd[[2]]
  
  SLmodel <- stackedSL(SL.lib, stackedData= stackedData.disc, #newXdata = newXstackedData.disc, 
                       verbose=verbose, numFolds=numSLFolds, method="method.NNLS")
  
  
  
  ##Get predictions along landmark time grid
  testIds <- unique(testDat$id)
  ##set up storage for predicted survival probabilities
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds
  colN <- paste0("LM", tseq)
  colN2 <- paste0("LM", tseq, "sq")
  
  for( i in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    
    thor <- tseq[i]+w
    thisLM <- paste0("LM", tseq[i])
    thisLM2 <- paste0("LM", tseq[i], "sq")
    sData <- testDat[testDat$survTime>tseq[i],]
    sData <- sData[sData$measTime <= tseq[i],] # delete all records where measurement collected after tseq[i]
    sData <- ddply(sData,.(id), tail,1) # Keep only the most recent record for each id
    ids2 <- as.character(unique(sData$id))
    delta.bound <- sort(unique(stackedData.disc[thisLM][stackedData.disc[thisLM]!=0] ))
    nhor <- length(delta.bound)
    ## Must put the data for prediction into the same format as we used to fit the model
    ND <- createPredictData(dataX= subset(sData[sData$id %in% ids2,], 
                                          select= -c(id, survTime, event, measTime)),time=delta.bound)
    
    ## We need to know delta.bound but don't want it in the data
    ND2 <- ND[1:(ncol(ND)-1)]
    ## add the landmark time dummy variables
    dum <- matrix(0, nrow=nrow(ND2), ncol=length(tseq))
    colnames(dum) <- colN
    rownames(ND2) <- NULL
    ND2 <- cbind(ND2, dum)
    # Add the landmark time as an interaction with time
    ND2[,thisLM ]<- ND$time   
    sfhor <- predict(SLmodel, newdata = ND2 ) 
    ## sfhor contains pred and library.predict -- these are conditional hazards
    ## there will be one row per row in dhor.disc
    #loop over individuals to calc survival probs
    for(j in 1:length(ids2)){
      beg <- ((j-1)*nhor) +1
      en <- j*nhor
      survP <- cumprod((1-sfhor$pred[beg:en]))
      testSurvPreds[i, ids2[j]] <- survP[nhor]
    }  
  }
  return(testSurvPreds)
}


##
## --- Function to create stacked dataset for landmarking
##

genStackedData <- function(dat, w, tseq){
  
  stackedData <- data.frame()
  dat$event <- dat$statusComp
  for(i in 1:length(tseq)){
    dat_i <- dat[dat$survTime>tseq[i],]   #keep only those still at risk at s  
    dat_i$event[dat_i$survTime>(tseq[i]+w)] <- 0   #ignore events after s+w
    dat_i$survTime[dat_i$survTime>(tseq[i]+w)] <- (tseq[i]+w)  #administratively censor at s+w
    dat_i <- dat_i[dat_i$measTime <= tseq[i],] # delete all records where measurement collected after s
    dat_i <- ddply(dat_i,.(id), tail,1)  # Keep only the most recent record for each id
    dat_i$lmTime <- tseq[i]                   # Add a column for the landmark time
    stackedData <- rbind(stackedData, dat_i) 
  }
  return(stackedData)
}   ##closes function



##
## --- Function to create stacked dataset for Super Learner landmarking
##     The difference is that this dataset will be discretised
##

genStackedSLData <- function(dat, w, tseq){
  n.delta <- 5  # number discrete intervals to use
  stackedData <- data.frame()
  colN <- paste0("LM", tseq)
  dat$event <- dat$statusComp
  for (i in 1:length(tseq)) {
    dat_i <- dat[dat$survTime>tseq[i],]   #keep only those still at risk at s
    dat_i$event[dat_i$survTime>(tseq[i]+w)] <- 0   #ignore events after s+w
    dat_i$survTime[dat_i$survTime>(tseq[i]+w)] <- (tseq[i]+w)  #administratively censor at s+w
    dat_i <- dat_i[dat_i$measTime <= tseq[i],] # delete all records where measurement collected after s
    dat_i <- ddply(dat_i,.(id), tail,1) #get the last record only
  
    # Discretize it based on quantiles of the event times
    delta.upper <- createDiscreteIntervals(time= dat_i$survTime, event=dat_i$event, 
                                         n.delta=n.delta)
    dat_i.newX <- createDiscrete2(time= dat_i$survTime, event=dat_i$event, 
                                 dataX = subset(dat_i, select=-c(survTime, event, statusComp, measTime, start)), 
                                 delta.upper = delta.upper, s=tseq[i])
    
    dat_i.disc <- dat_i.newX[!is.na(dat_i.newX$N.delta), ] 
    
    ## Add dummies for landmark time and interact with delta.lower
    dum <- matrix(0, nrow=nrow(dat_i.disc), ncol=length(tseq))
    colnames(dum) <- colN

    dat_i.disc <- cbind(dat_i.disc, as.data.frame(dum))
    thisLM <- paste0("LM", tseq[i])
    dat_i.disc[,thisLM ]<- dat_i.disc$delta.lower         # Add the landmark time as an interaction with time

    #rbind is known to be slow...
    stackedData <- rbind(stackedData, dat_i.disc)
  } 

return(stackedData)

}   


##
## --- Function to determine the boundaries of each discrete interval
##

createDiscreteIntervals <- function(time, event, n.delta) {
  # time is the vector of survival times
  # n.delta is the number of discrete intervals
  n <- length(time)
  n.delta <- min(n.delta, length(unique(time[event==1])))
  probs.delta <- seq(from=0, to=1, length.out=(n.delta+1))
  delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
  return(delta.upper[-1])	
}

## This function was adapted from code from the Github account of Eric Polley
## https://github.com/ecpolley/SuperLearner_Old/blob/master/R/createDiscrete.R
## Also see Chapter 16 "Super Learning for Right-Censored Data" in the book
## "Targeted Learning: Causal Inference for Observational and Experimental Data"
## by van der Laan and Rose, Springer, 2011.

createDiscrete2 <- function(time, event, dataX, delta.upper, s) {
  
  n <- length(time)
  delta.lower <- c(s, delta.upper[-length(delta.upper)])
  n.delta <- length(delta.upper)
  ID <- dataX$id
  dat <- cbind(ID, time, event, dataX)
  
  ## Note: Can use an apply statement here instead of pre-allocating and looping but
  ## apply converts data.frame to a matrix and uses least restrictive class which may be character
  long.dat <- matrix(NA, n.delta*n, ncol(dat))
  long.dat <- as.data.frame(long.dat)
  colnames(long.dat) <- colnames(dat)
  for( k in 1:nrow(dat)){ ##loop over ids
    b <- (n.delta*(k-1))+1
    e <- n.delta*k
    long.dat[b:e,] <- dat[k,]
  }
  
  N.delta <- rep(NA, nrow(long.dat))
  long.dat <- cbind(long.dat, delta.lower, delta.upper, N.delta)
  # Include censored people in the interval in which they were censored  
  long.dat$N.delta <- ifelse(long.dat$time > long.dat$delta.upper, 0, 
                              ifelse(long.dat$event==1, ifelse(long.dat$time <= long.dat$delta.lower, NA, 1), 
                                     ifelse(long.dat$time>long.dat$delta.lower, 0, NA)))
  
  m <- delta.upper[n.delta]
  long.dat$N.delta <- ifelse(long.dat$time == m & long.dat$delta.upper == m, 
                              ifelse(is.na(long.dat$N.delta), 0, long.dat$N.delta), long.dat$N.delta)
  
  return(long.dat)	
}


##
## --- This function is similar to the above but formats data we will make predictions on
##

createPredictData <- function(dataX, time){
  # Need to create a prediction set with the relevant data for this window
  n <- nrow(dataX)
  n.delta <- length(time)
  
  long.DATA <- apply(dataX, 2, function(x) rep(x, times=rep(n.delta, times=n)))
  long.DATA <- cbind(long.DATA, time)
  long.DATA <- as.data.frame(long.DATA)
  return(long.DATA)	
} 


##
## --- Function to apply SuperLearner to our prepared dicretised stacked dataset
##

stackedSL <- function(SL.lib, stackedData, newXdata=NULL, verbose=FALSE, 
                      numFolds=10, method="method.NNLS"){  
  
  
  ## fit the SuperLearner to the single stacked dataset
  set.seed(12)
  df.X <- subset(stackedData, select= -c(ID, time, event, id, delta.lower, delta.upper, N.delta ))

  cat("\nBegin SuperLearner ")
  cvControl = list(V=numFolds)

  SL.hor<- SuperLearner( Y = stackedData$N.delta, X = df.X,  
                           SL.library=SL.lib, 
                           method=method, verbose=verbose,
                           id=stackedData$id, family=binomial(), cvControl=cvControl)

  return(SL.hor)
}
