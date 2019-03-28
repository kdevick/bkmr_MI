****************************************************
BKMR plot functions for multiply imputed data
****************************************************

Code developed by:   Katrina Devick (kdevick@hsph.harvard.edu)
Last updated:   28 March 2019


******************
Getting Started
******************

The files example_code_for_MI_BKMR_fits.R and bkmr_MI_source.R that can be downloaded from the repository combine information from multiple Bayesian kernel machine regression (BKMR) models fit using the bkmr R package (Bobb et al. 2015, Valeri et al. 2017, Bobb et al. 2018, Anglen Bauer et al. 2019).  The file bkmr_MI_source.R contains functions to be used with MI BMKR fits to create a data frame for plotting with ggplot. The file example_code_for_MI_BKMR.R shows an example of how the functions are used with simulated data (stored in .RData file BKMRfit_MI_500.RData).


*******************
Important Notes
*******************
  
1) This code assumes you have K BKMR fits and that each of these fits were ran for the same number if MCMC iterations 
2) All effects are calculated as the average change in the outcome for a change in the exposure elements from a particular quantile to another quantile calculated across ALL imputed datasets. If there are no missing values in the Z matrix (in the mixture exposure), this is same contrast considered when only using the observed Z.  
3) All functions have the option to choose between an "approx" or "exact" method. The "exact" method combines the posterior samples from all MI fits and uses this posterior chain of length #iterations times #MI datasets for inference. The "approx" method uses the bkmr approx estimates and std errors from each MI fit and calculates an overall estimate and sd using Rubin's 1987 method.  
4) When using the "exact" method, know the functions take a while to run, so make sure you save the data frames to be used for plotting.


**************
References 
**************
Anglen Bauer J, Devick KL, Bobb JF, Coull BA, Zoni S, Fedrighi C, Benedetti C, Guazzetti S, White R, Bellinger D, Yang Q, Webster T, Wright RO, Smith D, Lucchini R, Claus Henn. Associations from a mixture of manganese, lead, copper and chromium and adolescent neurobehavior. 

Bobb JF, Claus Henn B, Valeri L, Coull BA. 2018. Statistical software for analyzing the health effects of multiple concurrent exposures via Bayesian kernel machine regression. Environ Health 17:67; doi:10.1186/s12940-018-0413-y.

Bobb JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, et al. 2015. Bayesian kernel machine regression for estimating the health effects of multi-pollutant mixtures. Biostatistics 16:493–508; doi:10.1093/biostatistics/kxu058.

Rubin DB. 1987. Multiple imputation for nonresponse in surveys. Wiley.

Valeri L, Mazumdar M, Bobb J, Claus Henn B, Sharif O, Al. E. 2017. The joint effect of prenatal exposure to metal mixtures on neurodevelopmental outcomes at 24 months: evidence from rural Bangladesh. Env Heal Perspect 125; doi:DOI: 10.1289/EHP614.