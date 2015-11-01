#!/bin/bash
start=`date`

#run=100k_fixed
#run=100k_exp_4
#run=100k_moroi
#run=tune357_10k
#run=pythia6_tune350_noBW_100TeV_100k
#run=pythia6_tune350_100TeV_10k
#run=pythia6_tune350_1000TeV_100k
#run=pythia6_def_100TeV_10k
#run=pythia8_def_BW_mDM.50_ww
#run=pythia8_noBW_1000TeV_100k
#run=def_100TeV_100k
#run=def_100TeV_1m
run=test

#./run_mode-parallel.sh $run 1 1
./run_dm_ann_parallel.sh $run tautau 0 1

echo $start
date