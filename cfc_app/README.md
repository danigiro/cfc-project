# Forecasting Australian daily electricity generation

This README provides step-by-step instructions for reproducing the forecasting experiments 
detailed in the paper. The experiments utilize the Australian daily electricity generation 
dataset (available at `data/energy`) and involve the implementation of various forecasting 
approaches.


## Instructions for reproducing the experiment

To replicate the analysis, please execute the following scripts in the specified order or 
execute the bash file `energy.sh`.

 1. `Rscript data/energy/extract_data.R`: extract the data from the csv file. 
 
 2. `Rscript R/00_base_forecast.R energy tbats`: compute the tbats base forecasts\
    `Rscript R/00_base_forecast.R energy stlf`: compute the stlf base forecasts\
    `Rscript R/00_base_forecast.R energy arima`: compute the arima base forecasts\
    `Rscript R/00_base_forecast.R energy snaive`: compute the snaive base forecasts

 3. `Rscript R/01_seq_reco.R energy sa tbats stlf arima`: sequential approaches using ew combination\
    `Rscript R/01_seq_reco.R energy var tbats stlf arima`: sequential approaches using ow(var) combination\
    `Rscript R/01_seq_reco.R energy cov0 tbats stlf arima`: sequential approaches using ow(cov) combination\
    `Rscript R/01_opt_reco.R energy tbats stlf arima`: optimal coherent forecast combination
  
  4. `Rscript R/02_score.R energy`: evaluation scores
  
  5. `Rscript R/03_tables.R energy`: replicate AvgRelMSE and AvgRelMAE tables\
     `Rscript R/03_mcb.R energy`: replicate Multiple Comparison with the Best\
     `Rscript R/03_dm.R energy`: replicate Diebold-Mariano test\
     `Rscript R/03_mcs.R energy`: replicate Model Confidence Sets tables