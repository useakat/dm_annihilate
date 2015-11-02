#!/bin/bash

run=$1
imkdir=$2
mail=$3

job_system=bsub
que=l
job=annDM$RANDOM
submit_mode=h

 ./submit_jobs.sh $job_system $que 1 $job "./run_parallel.sh $run ww 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 2 $job "./run_parallel.sh $run zz 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 3 $job "./run_parallel.sh $run hh 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 4 $job "./run_parallel.sh $run tautau 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 5 $job "./run_parallel.sh $run ddbar 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 6 $job "./run_parallel.sh $run uubar 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 7 $job "./run_parallel.sh $run ssbar 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 8 $job "./run_parallel.sh $run ccbar 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 9 $job "./run_parallel.sh $run bbbar 1 0" $submit_mode
 ./submit_jobs.sh $job_system $que 10 $job "./run_parallel.sh $run check 1 0" $submit_mode

./monitor

makedir.sh rslt_$run $imkdir
i=1
while [ $i -le 10 ];do
    mv par_$i/rslt_$run/* rslt_${run}/.
    i=`expr $i + 1`
done
rm -rf par_*

if [ $mail -eq 1 ];then
    bsub -q e -J annDM -u takaesu@post.kek.jp nulljob.sh >/dev/null 2>&1
fi