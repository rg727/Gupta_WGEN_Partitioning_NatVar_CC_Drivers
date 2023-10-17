# Gupta_WGEN_Partitioning_NatVar_CC_Drivers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7693324.svg)](https://doi.org/10.5281/zenodo.7693324)

This repository contains example code and data for the following paper:

Gupta, R.S., Steinschneider,S. & Reed, P.M. (2023). Understanding Contributions of Paleo-Informed Natural Variability and Climate Changes on Hydroclimate Extremes in the Central Valley region of California (Accepted in Earth's Future).

## Contents

Below we provide steps and code (and example data where possible) that are applied in this paper. Please refer to the Zenodo zipped archive for example data. 

### Step 1: Fit the Non-Homogeneous Hidden Markov Model (NHMM) 

* Folder: `./NHMM` <br />
* Input Data: `hgt500.synoptic.region` (too big to share, ask if needed) <br />
* Script: `fit_nhmm.R`<br />
* Description: The 5-state NHMM is fit with the first 9 principal components of 500 hpa GPH from the period of 1950-2013 and the covariates are the 4 principal components of annual weather regime occurrence (reconstruction from Gupta et al. (2022) that overlap with the same time period (1950-2013) which we generate within the script). The NHMM creates a Viterbi sequence of states using the covariates. You can plot the frequency of each WR. Once the NHMM is fit, you can retrieve the symbolic transition matrix formulas. We set up code to populate that matrix for each day of the 600 year period, using the daily covariates defined from the reconstruction. Finally, we simulate daily Monte Carlo chains of states through the whole reconstruction period. <br />
* Output Product: a matrix of 50 ensembles of 600-year (daily resolution) Monte Carlo chains of weather regimes <br />


### Step 2: Run the baseline weather generator
* Folder: `./WGEN` <br />
* Input Data: `TLG_Prcp.txt`, `TLG_Temp.txt` (regional, historic, precipitation and temperature for the Tuolumne Basin), the MC chains from step 1 (too big to share, ask if needed) <br />
* Scripts: `config.wgen.simulations.R`, `run_generator.sh` (bash script for a cluster), helper functions <br />
* Description: We have the weather generator set up to run across the five sub-basins in the study, but we’re just sharing an example script for running the generator for the Tuolumne basin. Helper functions are included. The most up to date (and generalized) generator can be found here: https://github.com/nassernajibi/WGEN-v2.0 . As the chains of WRs are read in, local weather from the historical period is bootstrapped for the cold season (November – April). The whole warm season for each year (May-October) is bootstrapped randomly from the historical data (since WRs aren’t active during this time). <br />
* Output Product: .rds files of precipitation, tmin, and tmax for each grid cell, for the whole 600-year period (too large to share, ask if needed)<br />
### Step 3: Apply climate change perturbations  
* Folder: `./CC_Changes`<br />
* Input Data: the temperature and precipitation .rds files from Step 2 (too big to share)<br />
* Scripts: `perturb.R`, helper functions, bash scripts for cluster submission<br />
* Description: This script takes the baseline weather and applies specified climate changes to it (either temperature shifts or precipitation scaling). <br />
* Output Product: Creates parquet files for each grid cell, for the 50 ensemble members for each climate scenario (combinations of temperature trend/cc scaling). <br /> 
### Step 4: Run the hydrologic model 
* Folder: `./Hydrologic_Modeling`<br />
* Input Data: the parquet files of temperature and precipitation from Step 3<br />
* Scripts: `run_sacsma.m`, helper functions and data, bash scripts for cluster submission<br />
* Description: This script takes the temperature and precipitation for each grid cell and uses a calibrated version of SAC-SMA and then routes the flow to a single gauged location in the basin. <br />
* Output Product: Creates streamflow at an outlet gauge (in this example, Don Pedro in the Tuolumne Basin) for the 50 ensemble members.   
### Step 5: Generate Metrics 
* Folder: `./Metrics`<br />
* Input Data: streamflow text files for each basin (providing Tuolumne), either baseline or with climate changes applied<br />
* Scripts:<br />
  * create_flood_metrics.R, plot_flood_metrics.py<br />
    * Description: this script takes a streamflow ensemble and fits a GEV distribution to 3-day peak annual flow across 30 or 100-year moving windows. It then calculates the 10-year and 100-year flood for each moving window.<br />
  * create_drought_metrics.R, plot_drought_metrics.py<br />
    * Description: this script takes monthly streamflow and transforms it into an SSI index. Then we apply methods of characterizing drought occurrence, severity, and duration across 30 or 100-year moving windows  <br />
  * create_joint_flood_metrics.R, plot_joint_metrics.py<br />
    * Description: This script fits a copula to 3-day peak annual flows across a combination of basins. Then we calculate the likelihood of exceeding the historical 100-year flood.      
### Step 6: Run Variance Decomposition 
* Folder: `./Variance_Decomposition`<br />
* Input Data: all 50 ensemble members of the 25 climate scenarios (too large to share)<br />
* Script: `run_anova.R`<br />
* Description: in order to perform a variance decomposition, you need the metric data for all the scenarios. The script we are attaching shows how to set up the ANOVA decomposition (for any metric, but using the 10-year flood as an example) and plot the final results (if you generate your own full suite of scenarios).  <br />


