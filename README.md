# DynamicPrediction
R Code to accompany the manuscript "Dynamic Survival Prediction Combining Landmarking with a Machine Learning Ensemble: Methodology and Empirical Comparison" by
Kamaryn T. Tanner, Linda D. Sharples, Rhian M. Daniel and Ruth H. Keogh. (2021) Journal of the Royal Statistical Society Series A, 184 (1) p3-30. https://doi.org/10.1111/rssa.12611 

For questions and comments regarding the code, please contact Kamaryn Tanner (kamaryn.tanner1@lshtm.ac.uk). As code is updated, it will be posted to https://github.com/KamTan/DynamicPrediction. 

This code was written in R v3.5.3 running on a Windows 10 PC.

We ilustrate the use of this code with the Mayo Clinic Primary Biliary Cirrhosis dataset publicly available via the R package survival. Although we are not permitted to make the UK CF Registry data used in the manuscript public, interested researchers may apply for this data by following the instructions provided here:
https://www.cysticfibrosis.org.uk/the-work-we-do/uk-cf-registry/apply-for-data-from-the-uk-cfregistry

The code is organised in two files:

 * In Main.R, we load the necessary libraries, set parameters used in the analysis, and run through a 10-fold external validation loop calculating predicted dynamic survival probabilities for each of the 10 test folds using the three methods compared in the paper: joint modeling, Cox landmarking and Super Learner landmarking. 
 
 * SupportFxns.R contains all of the functions called in Main.R to implement the three methods.
