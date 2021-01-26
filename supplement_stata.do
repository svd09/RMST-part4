/* Supplement to the paper - 
Model-free estimates that complement information obtained from the Hazard Ratio. 

Salil V Deo MD1,2, Vaishali S Deo MBBS DNB MPH1, Varun Sundaram MD1,2

1 – Louis Stokes VA Medical Center, NorthEast Ohio Veteran Affairs Healthcare system
2 – School of Medicine, Case Western Reserve University, Cleveland Ohio 

This supplemental file provides an example of how to measure the model free 
parameters described in the paper using STATA. 
The dataset used for the analysis is the lung cancer dataset provided in the survival
package in R. This dataset contains mortality data regarding males and females with advanced lung cancer that underwent surgery.
This data is available for download as a .csv file and can be imported into STATA. */

import delimited "H:\\ICTVS_review_biostatistics\\surv4\\data.csv" // place your filepath here.

// now to stset the data prior to doing further analyses.
// am going to convert the time to years to make the analysis easier to understand.
// also going to add 1 to ensure that we do not have surv_years = 0.


gen surv_years = (time + 1)/365.25

// for calculations, it is easier to have status and sex, our main variables coded as 0/1.
// to convert

gen status2 = status - 1 // died = 1

gen female = sex - 1 // female = 1

//stset the data.

stset surv_years, fail(status2) // once the data is stset, we will graph the KM curve
// to see it before further calculations.

sts graph, by(female) surv // The graph shows that females have better survival than males.
// As the patients have advanced lung cancer, survival is poor and almost everyone is dead // @ 2 years.

// 1. Median survival time -- The median survival time is that time at which 50% 
// of the cohort is surviving.

stci, by(female)

/*
         failure _d:  status2
   analysis time _t:  surv_years

             | Number of 
female       |  subjects         50%      Std. Err.    [95% Conf. Interval]
-------------+-------------------------------------------------------------
           0 |       138    .7419575      .0733424      .577686     .84052
           1 |        90    1.169062      .1210294      .947296    1.43737
-------------+-------------------------------------------------------------
       Total |       228    .8514716      .0596099      .780287    .991102

As seen from the results, the median survival time for men = 0.74(0.57 - 0.84) years 
and for women 1.16 (0.94 - 1.43) years. It is clear that women have a significantly longer median survival time than men as the Confidence intervals do not overlap. However, we can 
perform a formal hypothesis test using laplace approximation methods. */

laplacereg surv_years female, failure(status2) q(50) reps(500)

/*

Laplace regression can be used to compare median survival time (or any quantile survival times). The command requires the time & group followed by the failure variable if censored observations. Additionally, the q(50) denotes median survival time. Default setting is 500 reps.


Laplace regression                               No. of subjects  =        228
                                                 No. of failures  =        165
------------------------------------------------------------------------------
             |              Bootstrap
  surv_years |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
q50          |
      female |   .3267923   .1290863     2.53   0.011     .0737879    .5797967
       _cons |   .7392203   .0753544     9.81   0.000     .5915284    .8869121
------------------------------------------------------------------------------

As seen, median survival in females is 0.32 (0.07 - 0.57) years more than males (p = 0.011) */

// 2. Quantile survival: This same function can be used to obtain quantile survival times //and compare them too. 


// 3. Obtaining the difference in survival estimates between the two groups at specific 't'
// time periods.

stsurvdiff female, gen(diff) 

/* stsurvdiff can generate the difference in survival at each time point from the stset data.
Using gen(varname), we can create 3 columns with the difference difference_UL and difference_LL. Confidence interval can be specified. */

// 4. RMST - RMST can be obtained using STATA in multiple ways. Here, we will present the 
// the command strmst2 which will provide rmst using the non-parametric method. It also allows for adjustment of covariates in the same manner that the survRM2 package in R provides. In STATA there is another command available, strmst which uses the flexible PH parametric model introduced by Royston and Parmar. */

// once the data is stset.

strmst2 female, tau(1)
/*
Number of observations for analysis = 228
 
The truncation time: tau = 1 was specified.

Restricted Mean Survival Time (RMST) by arm
-----------------------------------------------------------
   Group |  Estimate    Std. Err.      [95% Conf. Interval]
---------+-------------------------------------------------
   arm 1 |     0.816       0.029        0.758        0.874
   arm 0 |     0.663       0.028        0.608        0.719
-----------------------------------------------------------

Between-group contrast (arm 1 versus arm 0) 
------------------------------------------------------------------------
           Contrast  |  Estimate       [95% Conf. Interval]     P>|z|
---------------------+--------------------------------------------------
RMST (arm 1 - arm 0) |     0.153        0.073        0.233      0.000
RMST (arm 1 / arm 0) |     1.230        1.103        1.373      0.000
------------------------------------------------------------------------

As seen here, for females rmst = 0.81 (0.75 - 0.87). That means, that at the end of 1 year, on average women survive for 0.81 years. For men, rmst = 0.66 (0.6 - 0.71) i.e. at the end of 1 year, men on average survive for 0.66 years. The difference in rmst between groups is 0.15(0.07 - 0.23)  with p < 0.001. Thus, we can conclude that, at the end of 1 year, women will live longer than men. On average, they will live 0.15 years i.e. 54 days longer. 

Readers are welcome to use the strmst command on their own.

References:

1. https://www.stata-journal.com/article.html?article=st0294 - A command for Laplace regression in STATA. - Bottai, Orsini.

2. https://journals.sagepub.com/doi/pdf/10.1177/1536867X1501500409 - Estimating the treatment effect in a clinical
trial using difference in restricted mean survival time. - Royston.

3. https://journals.sagepub.com/doi/pdf/10.1177/1536867X1601600310 - strmst2 and strmst2pw: New commands to
compare survival curves using the restricted mean survival time. Cronin, Tian, Uno.
*/



