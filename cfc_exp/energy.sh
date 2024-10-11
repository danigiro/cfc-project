#!/bin/bash
#Rscript data/energy/extract_data.R
#Rscript R/00_base_forecast.R energy tbats
#Rscript R/00_base_forecast.R energy stlf
#Rscript R/00_base_forecast.R energy arima
#Rscript R/00_base_forecast.R energy snaive
#Rscript R/01_seq_reco.R energy sa tbats stlf arima
#Rscript R/01_seq_reco.R energy var tbats stlf arima
#Rscript R/01_seq_reco.R energy cov0 tbats stlf arima
#Rscript R/02_opt_reco.R energy tbats stlf arima
#Rscript R/03_score.R energy
Rscript R/04_tables.R energy
Rscript R/04_mcb.R energy
Rscript R/04_dm.R energy


read -rsp $'Press enter to continue...\n'
