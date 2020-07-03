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
##  Required files:     SupportFxns.R which contains all functions called in this main file           
##  
## ------------------------------------------------------------------
#####################################################################


## ---- Load R libraries ---- ##
library(survival) ##PBC data, Cox prortional hazards model
library(JM)  ##joint model  (could also us JMbayes library)
library(plyr)
library(dplyr)
library(SuperLearner)

## Remember to set your working directory to the location of this file
source("SupportFxns.R") #Contains all the helper functions called in this file

## ---- Prepare dataset ---- ##

## Use the pbc dataset from the survival package and format for this analysis
pbcDat <- prepPBC()
head(pbcDat)

## ---- Set up storage for results ---- ##

methodJoint <- list()
methodCox_LM <- list()
methodSL_LM <- list()

## ---- Parameters for analysis ---- ##

#prediction horizon (years)
w <- 5

#grid of landmark times used for methods and evaluation
tseq <- seq(0,10,by=0.5)

## ---- SuperLearner Settings ---- ##

#Set verbose=TRUE to provide more output while algorithm runs
verbose = FALSE 
#Set number of cross-validation folds within Super Learner to 5 for faster run time
numSLFolds = 5 

#SL.libQuick is a small library of candidate algorithms chosen to have fast run time
# Note: must have gam and ranger packages installed 
SL.libQuick = list("SL.mean", "SL.glm", "SL.gam", "SL.ranger")



## ---- External cross-validation loop ---- ##

ECVfold <- 1  ##"external cross validation fold"
for( ECVfold in 1:10){
  cat("\nFold = ", ECVfold)
  
  # ---- Step 1: Segement data by fold to create test/training sets
  
  testIndexes <- which(pbcDat$fold==ECVfold)
  testDat <- pbcDat[testIndexes, ]
  trainDat <- pbcDat[-testIndexes, ]
  
  # ---- Step 2: Obtain predicted survival probabilities from each method ---- ##
  
  ## For each method, we compute predicted dynamic survival probabilities using:
    ##  Longitudinal predictor:   bili (endogenous)
    ##  Time-dependent predictor: spiders (assumed to be exogenous for this example)
    ##  Time-fixed predictors:    age, sex 
  
  ## Joint Model
  methodJoint[[ECVfold]] <- runModelJoint(testDat, trainDat, w, tseq )
  
  ## Cox landmarking
  methodCox_LM[[ECVfold]] <- runModelCox_LM(testDat, trainDat, w, tseq)
  
  ## Super Learner landmarking
  methodSL_LM[[ECVfold]] <- runModelSL_LM(testDat, trainDat, w, tseq, SL.lib=SL.libQuick)
  
  ##Predicted survival probabilities for each fold are stored in a list by method
  ##and can be compared using Brier score, C-index, etc.
}






  
  



  