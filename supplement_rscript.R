##################################################################
##          MODEL FREE ESTIMATES FOR SURVIVAL ANALYSIS          ##
##################################################################
# 01/04/2021
# This rscript is a supplement to the paper:

# Model-free estimates that complement information obtained from the Hazard Ratio. 

# Salil V Deo, Vaishali Deo, Varun Sundaram
#
# Louis Stokes Cleveland VA Medical Center,Cleveland Ohio
# svd14 at case.edu

# This is an rscript that is a supplement file to the paper on model-free estimates  #
# for reporting results of survival analysis. This file will provide 
# readers an introduction to the package and functions that can be used
# to obtain model free estimates for their data.

# Packages used:
# The packages that are beneficial for these calculations are:
# survRM2, surv2sampleComp, and ComparisonSurv.
# All these packages can be installed from CRAN.

install.packages("survRM2", dependencies = T)

install.packages("surv2sampleComp", dependencies = T)

install.packages("ComparisonSurv", dependencies = T)

# apart from these, we will load some useful packages 

library(pacman)

p_load(tidyverse, survival, survminer, survRM2, 
       surv2sampleComp, ComparisonSurv, rms, Hmisc)


# we will use the lung package in the survival library to 
# work on the functions provided in these packages.

df <- lung 

# see the data.

glimpse(df)

# for these functions, we need to make some changes
# status - needs to be coded as 0/1

df$died = df$status - 1 # this converts it to 0/1 with 1 = died.

df$female = df$sex - 1 # here now , 0 = males, 1 = females.

# surv2sample function in the surv2sampleComp can be used to 
# obtain ratio/difference in survival percentile estimates,
# ratio/difference in median survival time, and 
# it can also be used to obtain ratio/difference in t-year survival rates.

# using the function requires input of the following values:

# time - time till event or censoring.
# status - 0/1 censored/observed event.
# arm - two arms to be compared, they have to fit as 0/1.
# tau - this is the maximum time at which we would like to calculate the RMST.
# However, for the RMST, we will use the survRM2 package and not this function.
# quantprob - if we want to obtain survival probabilities at specific quantile values.
# 
# There are some other options like 'tau_start'. But, we will not discuss
# them here. They are not needed for routine right censored data. 
# The rest of the parameters can be kept at their default values.
# The conf.int paramter can be changed to provide results at different
#  confidence intervals. However, default = 0.95.

# We would recommend that the result for surv2sample be saved as an object.

# As we have also provided a supplemental file to provide examples in STATA, we will
# replicate the same analysis that is done there.

# thus after converting the variables female and status to 0/1 format, we will covert
# time to years.

df$surv_years <- (1 + df$time)/365.24

res <- surv2sample(time = df$surv_years,
                   status = df$died,
                   arm = df$female,
                   tau = 1,
                   quanprobs = 0.5, # to calculate the median survival time
                   SEED = 1974,
                   timepoints = c(1))

# as this depends upon bootstrapped resampling, please provide a seed.
# please note that if the status and arm are not coded as 0/1
# the function will throw an error. The error is not helpful at
# understanding the mistake. So, we would recommend to take note of this 
# important step prior to fitting the function.
# After the function is fit, 

# the output object is very large and contains many parts.
# please see all the parts at once.

str(res)

# This function provides values for each parameter according to group and then we can obtain the contrast
# between groups. Hence, unlike STATA output, we need only run this function to obtain
# all results.

# to see the results for group 0 - males.

res$group0

#                                Est. Lower 95% Upper 95%         SE
# RMST                      0.6632342 0.6078744 0.7162000 0.02796118
# Loss time                 0.3367658 0.2838000 0.3921256 0.02796118
# Prob at 1                 0.3360878 0.2575620 0.4206296 0.04272376
# Quantile at 50 %          0.7419779 0.5777023 0.8405432 0.07662344
# Ave of t-year event rates 0.3360878 0.2575620 0.4206296 0.04272376
# Ave percentiles           0.7419779 0.5777023 0.8405432 0.07662344


# Here -- the results we are interested in are (1) Prob at 1 - this is the 
# survival probability at 1 year, hence in males 0.33 % survive at 1 year. 
# (2) quantile at 50% - 0.74 , hence 50% are surviving at 0.74 years from the start.

res$group1 # this will similarly give the same values for the group = 1 (female)

# as can be seen from results, they are almost exactly similar to those
# provided by STATA.

# This function also provides the RMST at the time point specified by 
# input to tau. If we compare the RMST reported here, again, it is almost
# exactly the same as that reported by stmrst command in STATA. The small
# variation in the final decimal places is due to computation and resampling.

res$contrast.diff10
#         
#                                               Est.   Lower 95%   Upper 95%        p-val
# RMST Group1-Group0                       0.1528465  0.07255657  0.23313634 0.0001906082
# Loss time Group1-Group0                 -0.1528465 -0.23313634 -0.07255657 0.0001906082
# Prob at 1 Group1-Group0                  0.1903752  0.04464380  0.33610659 0.0104556172
# Quantile at 50 % Group1-Group0           0.4271164  0.08983419  0.76439864 0.0130649109
# Ave of t-year event rates Group1-Group0  0.1903752  0.04464380  0.33610659 0.0104556172
# Ave percentiles Group1-Group0            0.4271164  0.08983419  0.76439864 0.0130649109


# The result for difference between female and males with confidence intervals and p-value.

# This package surv2sampleComp provides all the parameters that we have discussed in the paper.
# The survRM2 is useful for computing the rmst and performing covariate adjusted rmst calculations.

# References:
        
# 1. https://cran.r-project.org/web/packages/survRM2/vignettes/survRM2-vignette3-2.html - 
# An excellent vignette prepared by developers of survRM2 to explain their main function rmst2.

