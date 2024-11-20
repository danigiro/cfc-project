# Simulation experiment

This README provides step-by-step instructions for reproducing the simulation experiments 
detailed in the paper.

To replicate the analysis, please execute the bash file `run_bal.sh` (balance framework) and
`run_unb.sh` (unbalance framework). Then, to replicate tables and figures:

 1. `Rscript R/extract_score.R`
 
 2. `Rscript R/tables.R`