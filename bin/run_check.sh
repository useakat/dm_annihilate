#!/bin/bash
selfdir=$(cd $(dirname $0);pwd)

run=$1
channel=$2
imkdir=$3
mail=$4
###### MODIFY HERE: running parameters #################
output=hadron_dist_$channel.dat
jobname=dm_check
submit_mode=0 # 0:serial submittion 1:parallel submission
#cluster=kekcc
cluster=icrr  # name of computer cluster: kekcc/icrr
que=l

min=10
max=1000
imin=1
ndiv=30
ext=check

nevents=1000
logflag=1

imax=`expr $ndiv + 1`

# working space for jobs on a remote server
if [ $cluster == "icrr" ];then
    work_dir=/disk/th/work/takaesu/$run
    mkdir $work_dir
elif [ $cluster == "kekcc" ];then
    work_dir=./
fi
######################################################
start=`date`
echo $start

bin_dir=$selfdir

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
    $bin_dir/submit_jobs.sh $cluster $que $i $job "./run_general.sh run_$i $x $nevents $channel" $submit_mode $work_dir
    i=`expr $i + 1`
done
n=$i

if [ $submit_mode -eq 1 ];then
    $bin_dir/monitor
fi
if [ $cluster == "icrr" ];then
    mv $work_dir/* .
    rm -rf $work_dir
fi

rsltdir=$bin_dir/../results/rslt_$run
$bin_dir/makedir.sh $rsltdir $imkdir
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

$bin_dir/mail_notify $mail $cluster $jobname