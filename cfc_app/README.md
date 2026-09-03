# Forecasting Australian daily electricity generation

This README provides step-by-step instructions for reproducing the forecasting experiments 
detailed in the paper. The experiments utilize the Australian daily electricity generation 
dataset (available at `data/energy`) and involve the implementation of various forecasting 
approaches.


## Data

The experiment uses daily electricity generation data for the Australian
National Electricity Market (NEM), broken down by generation source, covering
366 days from 11 June 2019 to 10 June 2020.

The data were obtained from the replication material of Panagiotelis,
Gamakumara, Athanasopoulos and Hyndman,
[*Probabilistic Forecast Reconciliation*](https://github.com/PuwasalaG/Probabilistic-Forecast-Reconciliation).
They originate from the Australian Energy Market Operator (AEMO) and were
retrieved through [OpenElectricity](https://explore.openelectricity.org.au),
formerly OpenNEM.

`data/energy/daily.csv` is the raw file as obtained from that source: a CSV in
UTF-8 with byte-order mark, 42 columns, of which this experiment uses the 15
generation series in GWh and discards the emission-volume ones. It is included
in this package, so **the raw data do not have to be downloaded**.

### Pre-processing

`data/energy/extract_data.R` (step 1) turns the raw file into the three
files the experiment reads, following the same construction as the source
replication material:

-   the last observation is dropped;
-   the negative values of `Biomass` and `Distillate`, which are physically
    non-negative, are set to zero;
-   `Pumps` and `Battery (Charging)` are sign-reversed, as they enter the
    aggregation with a negative sign;
-   the 15 bottom-level series are combined into 8 upper-level ones through the
    aggregation matrix of Figure 2 of the manuscript.

The outputs are `data.csv` (the 23 series, upper and bottom), `agg_mat.csv` (the
aggregation matrix) and `inds.rds` (the dates). All three are included, and
`extract_data.R` regenerates them identically from `daily.csv`.

### Terms of use

The data are redistributed here as part of the replication material of a
published article, for the sole purpose of reproducing its results. They
originate from AEMO, whose market data are publicly available. **The source
repository declares no licence**, so the terms under which these data may be
further redistributed are not formally established; anyone intending to reuse
them beyond the reproduction of this paper should refer to the
[AEMO](https://aemo.com.au) and [OpenElectricity](https://explore.openelectricity.org.au)
terms of use, and cite the original replication material.

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

## Tables and figures of the paper

| Paper | Element | Script (step) | File |
|-------|---------|---------------|------|
| Manuscript | Figure 2 (right panel) | `R/00_plot_ts.R` | `img/ts_energy/agg_mat_small.pdf` |
| Manuscript | Table 7 | `R/03_tables.R` (5) | `tables/p_energy.tex` |
| Manuscript | Figure 3 | `R/03_dm.R` (5) | `img/p_energy_dm.pdf` |
| Manuscript | Table 8 | `R/03_mcs.R` (5) | `tables/p_energy_mcs.tex` |
| Online Appendix | Figure H.1 | `R/00_plot_ts.R` | `img/ts_energy/plot_all_ts_energy.pdf` |
| Online Appendix | Table H.1 | `R/03_tables.R` (5) | `tables/oa_energy_mae.tex` |
| Online Appendix | Table H.2 | `R/03_tables.R` (5) | `tables/oa_energy_mse.tex` |
| Online Appendix | Figure H.2 | `R/03_mcb.R` (5) | `img/p_energy_mcb.pdf` |
| Online Appendix | Figure H.3 | `R/03_dm.R` (5) | `img/oa_energy_dm.pdf` |
| Online Appendix | Table H.3 | `R/03_mcs.R` (5) | `tables/oa_energy_mae_mcs.tex` |
| Online Appendix | Table H.4 | `R/03_mcs.R` (5) | `tables/oa_energy_mse_mcs.tex` |
| Online Appendix | Figure H.4 | `R/03_mcb.R` (5) | `img/oa_energy_mae_mcb.pdf` |
| Online Appendix | Figure H.5 | `R/03_mcb.R` (5) | `img/oa_energy_mse_mcb.pdf` |

The left panel of Figure 2 (the generation hierarchy) is drawn in the manuscript
itself and is not produced by any script.

The scripts also write some files that do not appear in the paper: the
`tables/*_nored.tex` tables, identical to their counterparts above but without
the red highlighting of the entries worse than the benchmark, and
`img/ts_energy/plot_uts_ts_energy.pdf` and `plot_bts_ts_energy.pdf`, which show
the upper and the bottom time series separately.

## Intermediate files

Re-running the whole pipeline is not necessary to reproduce the results of the
paper. The following intermediate files are available on request from the
authors:

| Directory | Produced by | Content | Size |
|-----------|-------------|---------|------|
| `data/energy/` | step 1 (`extract_data.R`) | pre-processed data (`data.csv`, `agg_mat.csv`, `inds.rds`) | 184 KB |
| `fc/energy/base/` | step 2 (`00_base_forecast.R`) | base forecasts of the four models, one `.RData` per rolling window | 295 MB |
| `fc/energy/{sa,var,cov0,opt,time}/` | step 3 (`01_seq_reco.R`, `01_opt_reco.R`) | reconciled and combined forecasts, and the timings of each approach | 23 MB |
| `fc/energy/energy_*.rds` | step 4 (`02_score.R`) | scores and the inputs of the MCB, DM and MCS tests | 37 MB |

Depending on which files you have, you can enter the pipeline at any point:

-   with `fc/energy/base/`: start from step 3;
-   with the reconciled forecasts: start from step 4;
-   with `fc/energy/energy_*.rds`: start from step 5, which reproduces every
    table and figure of the paper in less than a minute.

## Runtime

Measured wall-clock time of the original run (Intel Core i7-10700 @ 2.90GHz,
64 GB RAM, Windows 10 build 19045, R 4.4.0):

| Step | Runtime |
|------|---------|
| Step 2, base forecasts — `tbats` | 86 min |
| Step 2, base forecasts — `arima` | 8 min |
| Step 2, base forecasts — `stlf` | 4 min |
| Step 2, base forecasts — `snaive` | < 1 min |
| Steps 3–4, reconciliation, combination and scores | ~ 2 min |
| Step 5, tables and figures | < 1 min |

The four base forecasting models are independent and can be run concurrently, in
which case step 2 takes as long as `tbats` alone. The timings of the
reconciliation and combination step are recorded by the scripts themselves in
`fc/energy/time/`.

**The tables and figures of the paper can be reproduced with step 5 alone**,
without recomputing the forecasts, using the intermediate files described above.

## Hardware used

Intel Core i7-10700 @ 2.90GHz (8 cores / 16 threads), 64 GB RAM, Windows 10 x64
(build 19045), R 4.4.0. See the *Computing environment* section of the README
at the root of the repository for the full details and the package versions.
