# apply LBSPR to different simulated data sets
# need to calculate the 'true' SBPR for each
# scenarios: constant recruitment, single month
# scenarios: sporadic recruitment, single month
# scenarios: constant recruitment, double month
# scenarios: sporadic recruitment, double month

# Pull data from mantis_om, populate below
# compare the actual Bzero at a given time to the estimated for different scenarios
# and 

devtools::install_github("AdrianHordyk/LBSPR")
library(LBSPR)

MyPars <- new("LB_pars")

MyPars@Species <- "MySpecies"
MyPars@Linf <- 100 
MyPars@L50 <- 66 
MyPars@L95 <- 70
MyPars@MK <- 1.5 
MyPars@L_units <- "mm"
datdir <- DataDir()
Len1 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LFreq_MultiYr.csv"), 
            dataType="freq")
plotSize(Len1)

myFit1 <- LBSPRfit(MyPars, Len1)
myFit1@Ests

plotMat(myFit1)
plotEsts(myFit1)



