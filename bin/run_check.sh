#!/bin/bash
run=$1
decay=$2
imkdir=$3
mail=$4
###### MODIFY HERE: running parameters #################
output=hadron_dist_$decay.dat
jobname=dm_ann_check
submit_mode=0 # 0:serial submittion 1:parallel submission
job_system=bsub
que=l

min=10
max=1000
imin=1
ndiv=30
ext=check

nevents=1000
logflag=1

imax=`expr $ndiv + 1`
######################################################
start=`date`
echo $start

i=$imin
while [ $i -le $imax ];do
    job=${jobname}_$i
    if [ $i -eq 1 ];then
	x=$min
    elif [ $i -le $ndiv ];then
	if [ $logflag -eq 1 ];then
	    x=`echo "scale=5; e( (l($min)/l(10) +(l($max)/l(10) -l($min)/l(10))/$ndiv*($i-1))*l(10) )" | bc -l`
	else
	    x=`echo "scale=5; $min +($max -$min)/$ndiv*($i-1)" | bc -l`
	fi
    else
	x=$max
    fi
    ./submit_job_dm_ann.sh $job_system $que $i $job "./run_dm_ann_general.sh run_$i $x $nevents $decay > allprocess.log" $submit_mode
    i=`expr $i + 1`
done
n=$i

if [ $submit_mode -eq 1 ];then
    ./monitor
fi

rsltdir=rslt_$run
makedir.sh $rsltdir $imkdir
rm -rf $rsltdir/$output

x=$min
i=$imin
while [ $i -lt $n ];do
### MODIFY HERE for preparing result files ###############
    if [ $i -eq $imin ];then
	cat par_$i/data/run_$i/Nch.dat > $rsltdir/Nch_$ext.dat
    else
	cat par_$i/data/run_$i/Nch.dat >> $rsltdir/Nch_$ext.dat
    fi
###################################################
    i=`expr $i + 1`
done

### MODIFY HERE for saving files relatee to this run
cp -rf par_$imin/run_dm_ann_general.sh $rsltdir/.
cp -rf submit_job_dm_ann.sh $rsltdir/.
cp -rf run_dm_ann_parallel_moroi.sh $rsltdir/.
###################################################

echo "finished!"
echo $start
echo `date`

rm -rf par_*

if [ $mail -eq 1 ];then
    bsub -q e -J annDM -u takaesu@post.kek.jp nulljob.sh >/dev/null 2>&1
fi