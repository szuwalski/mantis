# mantis

This repo contains an operating model to simulate a size-structured population on a monthly time step coupled with an size structured assessment method that can be applied to data drawn from the operating model. To implement this analysis, the user must specify the desired characteristics of the population and fishery in 'mantis_OM.R' and run the script. This will write .PIN and .DAT files to '/admb', after which the model 'mantis.exe' can be run at the command line. The model fits can then be examined by running the file 'mantis_EM.R'.

# things to be done

Add fleets, distinguish between sexes, specify size bins to distribute recruitment, implement a stock recruit relationship in OM.
