#!/bin/bash
selfdir=$(cd $(dirname $0);pwd)
start=`date`

run=test

###  simulate tau tau mode with run name "$run". Send mail when jobs are finished. 
$selfdir/run_general_parallel.sh $run 1 tautau 0 1

echo $start
date