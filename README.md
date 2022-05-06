# Runoff, erosion and related pesticide transport
This repository contains the code and secondary data used for analysis and calculations in:

Pesticides are substantially transported in particulate phase, driven by land use, rainfall event and pesticide characteristics – a runoff and erosion study in a small agricultural catchment. M.C. Commelin et al. (....).

The primary data can be obtained from the 4TU data repository at: [fill link and DOI]

The main calculations are performed in 'Paper_overview.R', this code follows the order of the paper as much as possible.
Calculations of LC-MS/MS data to pesticide concentrations is done in 'pesticide_calculations.R'

The secondary data is stored in 'ext_data/' and functions in R and intermediate data files are provide in 'sources/'

# Performing all calculations and analysis
To run the code and analysis for this research the following software and data is required:
 - R, download from: https://cran.r-project.org/
 - RStudio, download from: https://www.rstudio.com/
 - the primary dataset
 
This repository should be downloaded, and the primary dataset added in a subfolder 'data/'. Then open the file 'paper_analysis.Rproj'. From the different R codes can be executed, after installation of required packages with: 'install.packages("package name")'.


# Session info
The calculations for the original manuscript where done using the following software and package versions in R.
**R version 4.1.2 (2021-11-01)**

**Platform:** x86_64-w64-mingw32/x64 (64-bit) 

**locale:**
_LC_COLLATE=English_Netherlands.1252_, _LC_CTYPE=English_Netherlands.1252_, _LC_MONETARY=English_Netherlands.1252_, _LC_NUMERIC=C_ and _LC_TIME=English_Netherlands.1252_

**attached base packages:** 
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_rstatix(v.0.7.0)_, _ggpubr(v.0.4.0)_, _pander(v.0.6.4)_, _broom(v.0.7.12)_, _ggforce(v.0.3.3)_, _extrafont(v.0.17)_, _hms(v.1.1.1)_, _ggrepel(v.0.9.1)_, _sf(v.1.0-7)_, _scales(v.1.1.1)_, _zoo(v.1.8-9)_, _cowplot(v.1.1.1)_, _lubridate(v.1.8.0)_, _aomisc(v.0.647)_, _multcompView(v.0.1-8)_, _car(v.3.0-12)_, _carData(v.3.0-5)_, _plyr(v.1.8.6)_, _drc(v.3.0-1)_, _MASS(v.7.3-55)_, _forcats(v.0.5.1)_, _stringr(v.1.4.0)_, _dplyr(v.1.0.8)_, _purrr(v.0.3.4)_, _readr(v.2.1.2)_, _tidyr(v.1.2.0)_, _tibble(v.3.1.6)_, _ggplot2(v.3.3.5)_ and _tidyverse(v.1.3.1)_

**loaded via a namespace (and not attached):** 
_fs(v.1.5.2)_, _bit64(v.4.0.5)_, _httr(v.1.4.2)_, _tools(v.4.1.2)_, _backports(v.1.4.1)_, _utf8(v.1.2.2)_, _R6(v.2.5.1)_, _KernSmooth(v.2.23-20)_, _DBI(v.1.1.2)_, _colorspace(v.2.0-3)_, _withr(v.2.5.0)_, _tidyselect(v.1.1.2)_, _bit(v.4.0.4)_, _compiler(v.4.1.2)_, _extrafontdb(v.1.0)_, _cli(v.3.2.0)_, _rvest(v.1.0.2)_, _xml2(v.1.3.3)_, _sandwich(v.3.0-1)_, _classInt(v.0.4-3)_, _mvtnorm(v.1.1-3)_, _proxy(v.0.4-26)_, _digest(v.0.6.29)_, _pkgconfig(v.2.0.3)_, _plotrix(v.3.8-2)_, _dbplyr(v.2.1.1)_, _rlang(v.1.0.2)_, _readxl(v.1.3.1)_, _rstudioapi(v.0.13)_, _generics(v.0.1.2)_, _farver(v.2.1.0)_, _jsonlite(v.1.8.0)_, _gtools(v.3.9.2)_, _vroom(v.1.5.7)_, _magrittr(v.2.0.2)_, _Matrix(v.1.4-0)_, _Rcpp(v.1.0.8)_, _munsell(v.0.5.0)_, _fansi(v.1.0.2)_, _abind(v.1.4-5)_, _lifecycle(v.1.0.1)_, _stringi(v.1.7.6)_, _multcomp(v.1.4-18)_, _grid(v.4.1.2)_, _parallel(v.4.1.2)_, _crayon(v.1.5.0)_, _lattice(v.0.20-45)_, _haven(v.2.4.3)_, _splines(v.4.1.2)_, _pillar(v.1.7.0)_, _ggsignif(v.0.6.3)_, _codetools(v.0.2-18)_, _reprex(v.2.0.1)_, _glue(v.1.6.2)_, _modelr(v.0.1.8)_, _vctrs(v.0.3.8)_, _tzdb(v.0.2.0)_, _tweenr(v.1.0.2)_, _Rttf2pt1(v.1.3.10)_, _cellranger(v.1.1.0)_, _gtable(v.0.3.0)_, _polyclip(v.1.10-0)_, _assertthat(v.0.2.1)_, _e1071(v.1.7-9)_, _class(v.7.3-20)_, _survival(v.3.2-13)_, _units(v.0.8-0)_, _TH.data(v.1.1-0)_ and _ellipsis(v.0.3.2)_