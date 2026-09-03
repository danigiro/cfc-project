#!/bin/bash

settlist="sett1 sett2 sett3 sett4 sett5 sett6"
plist="4 10 20"
nlist="50 100 200"
structure="hier"
cor_hat="cor diag"

for cor in $cor_hat
  do
  for s in $structure
  do
    for n in $nlist
    do
      for p in $plist
      do
        for sett in $settlist
        do
          echo "Running $sett $p $n $s $cor"
          Rscript R/sim_unb_parallel.R $sett $p $n 100 $s $cor 
        done
      done
    done
  done
done

