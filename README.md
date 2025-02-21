# EPICC immune analysis scripts
Scripts performing immune-related analysis within the Evolutionary Predictions on Colorectal Cancer (EPICC) cohort.  
Our manuscript describing and analysing the results is available on bioRxiv: https://www.biorxiv.org/content/10.1101/2024.02.12.579956v1. 

Files named X.proc.* contain processing scripts that create high-level summary tables from raw/low-level data.
Files named X.fig.* contain scripts that plot and interpret the high-level summary tables.

All data used in these scripts are available on Mendeley: https://doi.org/10.17632/cjfmmc95dm.2

## Requirements and running instructions
All package dependencies are listed and evoked in 0.basics.R. For installation of these packages, refer to their respective user manuals. Installation time for all packages together should take 5-20 minutes on a standard laptop.  
Scripts were written and tested under R version 4.4.2 (Mac, Apple silicon), but should be compatible with all versions after 4.0.

To run all scripts, it is essential to set up the path to the full dataset by modifying the following line in 0.basics.R.
```
setwd('~/EPICC_immune_analysis/EPICC_example/') # Change this folder depending on the location of the downloaded files
```
Once the main folder (where data from Mendely was saved to) is specified, all subfolder structure should be automatically correct. Additionally, the line
```
source('0.basics.R')
```
may have to be changed depending on the folder structure used.

A subset of the data, sufficient to run 0.basics.R and 1.fig.overview.R, is included in the repository in subfolder EPICC_example. To test if the data folders are appropriately identified in R and packages are evoked correctly, first run 0.basics.R and then run 1.fig.overview.R. The code is expected to run within seconds, and produces unprocessed sub-panels of Figure 1b&c of the manuscript (see [here](https://www.biorxiv.org/content/10.1101/2024.02.12.579956v1.full.pdf)).



