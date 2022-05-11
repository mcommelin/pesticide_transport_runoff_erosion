# Analysis on pesticide transport during rainfall runoff events
This repository contains the code and secondary data used for analysis and calculations on runoff, erosion
and pesticide transport in particulate and dissolved phase. The analysis is done on an observational dataset
in which rainfall-runoff events where monitored in a small agricultural catchment in Limburg, The Netherlands.

**Author:**  
M.C. Commelin  
Soil Physics and Land Management group, Wageningen University, 6708 WG Wageningen, The Netherlands   
Contact: meindert.commelin@wur.nl

**Related publication:**  
Pesticides are substantially transported in particulate phase, driven by land use, 
rainfall event and pesticide characteristics – a runoff and erosion study in a small 
agricultural catchment. M.C. Commelin et al. (2022).  
DOI: 10.3389/fenvs.2022.830589 

**Related dataset:**  
Runoff, erosion and pesticide transport in a small catchment in Limburg - The Netherlands  
DOI: 10.4121/19690684  

## Performing all calculations and analysis
The main calculations are performed in 'Paper_overview.R', this code follows the order of the related paper as much as possible.
Analysis of LC-MS/MS data to derive pesticide concentrations is done in 'pesticide_calculations.R'. The final
pesticide cocnetration are stored in 'sources/lc_all_data.csv' and already provided in the repository.

The secondary data used for calculations in this research are stored in 'ext_data/'.
Functions in R and intermediate data files are provided in 'sources/'

To run the code and analysis for this research the following software and data is required:
 - R, download from: https://cran.r-project.org/
 - RStudio, download from: https://www.rstudio.com/
 - the primary dataset
 
This repository should be downloaded, and the primary dataset added in a subfolder 'data/'.
Then open the R project file: 'paper_analysis.Rproj'. From there the different R codes can be executed, after installation of required packages with: 'install.packages("package name")'.

## Session info  
The calculations for the original manuscript where done using the following software and package versions in R.

**R version 4.1.2 (2021-11-01)**

**Versions of loaded packages:** 
_rstatix(v.0.7.0)_, _ggpubr(v.0.4.0)_, _pander(v.0.6.4)_, _broom(v.0.7.12)_, _ggforce(v.0.3.3)_, _extrafont(v.0.17)_, _hms(v.1.1.1)_, _ggrepel(v.0.9.1)_, _sf(v.1.0-7)_, _scales(v.1.1.1)_, _zoo(v.1.8-9)_, _cowplot(v.1.1.1)_, _lubridate(v.1.8.0)_, _aomisc(v.0.647)_, _multcompView(v.0.1-8)_, _car(v.3.0-12)_, _carData(v.3.0-5)_, _plyr(v.1.8.6)_, _drc(v.3.0-1)_, _MASS(v.7.3-55)_, _forcats(v.0.5.1)_, _stringr(v.1.4.0)_, _dplyr(v.1.0.8)_, _purrr(v.0.3.4)_, _readr(v.2.1.2)_, _tidyr(v.1.2.0)_, _tibble(v.3.1.6)_, _ggplot2(v.3.3.5)_ and _tidyverse(v.1.3.1)_
