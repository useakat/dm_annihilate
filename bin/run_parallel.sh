#!/bin/bash
run=$1
decay=$2
imkdir=$3
mail=$4
###### MODIFY HERE: running parameters #################
output=hadron_dist_$decay.dat
jobname=dm_ann
submit_mode=1 # 0:serial submittion 1:parallel submission
#cluster=kekcc # name of computer cluster: kekcc/icrr
cluster=icrr  # name of computer cluster: kekcc/icrr
que=l

#Emax=1000000
Emax=1000000  # default value: 1000000
if [ $decay == "ww" ];then
    min=100  # default 100 GeV
    max=$Emax
    imin=1
    ndiv=20  # default 5points/1order
#    ndiv=3
    ext=ww
elif [ $decay == "zz" ];then
    min=100
    max=$Emax
    imin=1
    ndiv=20
#    ndiv=15
    ext=zz
elif [ $decay == "hh" ];then
    min=100
    max=$Emax
    imin=2
    ndiv=20
#    ndiv=15
    ext=hh
elif [ $decay == "tautau" ];then
    min=1
    max=$Emax
    imin=3
#    ndiv=25
    ndiv=30
    ext=tautau
elif [ $decay == "uubar" ];then
    min=1
    max=$Emax
    imin=2
    ndiv=30
#    ndiv=25
    ext=uubar
elif [ $decay == "ddbar" ];then
    min=1
    max=$Emax
    imin=2
    ndiv=30
#    ndiv=25
    ext=ddbar
elif [ $decay == "ssbar" ];then
    min=1
    max=$Emax
    imin=2
    ndiv=30
#    ndiv=25
    ext=ssbar
elif [ $decay == "ccbar" ];then
    min=1
    max=$Emax
    imin=3
    ndiv=30
#    ndiv=25
    ext=ccbar
elif [ $decay == "bbbar" ];then
    min=1
    max=$Emax
    imin=5
    ndiv=30
#    ndiv=25
    ext=bbbar
elif [ $decay == "ttbar" ];then
    min=100
    max=$Emax
    imin=3
    ndiv=20
    ext=ttbar
elif [ $decay == "check" ];then
    min=5
    max=100000
    imin=1
    ndiv=20
    ext=check
fi
logflag=1
#nevents=1000000
#nevents=100000
#nevents=1000

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

    ix=`echo $x |cut -d. -f1`
    if [ $ix -le 10 ];then
	nevents=1000000
    elif [ $ix -le 100 ];then
	nevents=1000000
    elif [ $ix -le 1000 ];then
	nevents=1000000
    elif [ $ix -le 10000 ];then
	nevents=1000000
    elif [ $ix -le 100000 ];then
	nevents=100000
    elif [ $ix -le 1000000 ];then
	nevents=10000
    else 
	nevents=1000
    fi

    ./submit_job_dm_ann.sh $cluster $que $i $job "./run_dm_ann_general.sh run_$i $x $nevents $decay" $submit_mode $work_dir
    i=`expr $i + 1`
done
n=$i
if [ $decay == "check" ];then
    job=${jobname}_$i
    x=91.2
    ./submit_job_dm_ann.sh $cluster $que $i $job "./run_dm_ann_general.sh run_$i $x $nevents $decay" $submit_mode $work_dir
    i=`expr $i + 1`
fi

if [ $submit_mode -eq 1 ];then
    ./monitor $work_dir
fi
if [ $cluster == "icrr" ];then
    mv $work_dir/* .
    rm -rf $work_dir
fi

rsltdir=rslt_$run
./makedir.sh $rsltdir $imkdir
rm -rf $rsltdir/$output

echo "1, 141, 300, 0, 30, 1" > $rsltdir/$output
echo "%%%%%" >> $rsltdir/$output
x=$min
i=$imin
while [ $i -lt $n ];do
### MODIFY HERE for preparing result files ###############
    mass=`cat par_$i/data/run_$i/mass.dat`
    Evis_tot=`cat par_$i/data/run_$i/Evis_tot.dat`
    echo $mass 0 0 $Evis_tot >> $rsltdir/$output
    cat par_$i/data/run_$i/Edist.dat >> $rsltdir/$output
    echo "%%%%%" >> $rsltdir/$output
    if [ $i -eq $imin ];then
# prepare outputfiles and write the first line
	cat par_$i/data/run_$i/np_sptrm.dat > $rsltdir/np_sptrm_$ext.dat
	cat par_$i/data/run_$i/nini.dat > $rsltdir/nini_$ext.dat
	cat par_$i/data/run_$i/Evis.dat > $rsltdir/Evis_$ext.dat
#	cat par_$i/data/run_$i/Nch.dat > $rsltdir/Nch_$ext.dat
	cat par_$i/data/run_$i/np_sptrm_norm.dat > $rsltdir/Ekin_z_$ext.dat
    else
# Add data line to the outputfiles
	cat par_$i/data/run_$i/np_sptrm.dat >> $rsltdir/np_sptrm_$ext.dat
	cat par_$i/data/run_$i/nini.dat >> $rsltdir/nini_$ext.dat
	cat par_$i/data/run_$i/Evis.dat >> $rsltdir/Evis_$ext.dat
#	cat par_$i/data/run_$i/Nch.dat >> $rsltdir/Nch_$ext.dat
	cat par_$i/data/run_$i/np_sptrm_norm.dat >> $rsltdir/Ekin_z_$ext.dat
    fi
    # copy log file
    cp -rf par_$i/allprocess.log $rsltdir/allprocess_$i.log
###################################################
    i=`expr $i + 1`
done

### MODIFY HERE for saving files relatee to this run
#cp -rf par_1/param_card.dat $rsltdir/.
#cp -rf par_1/run_card.dat $rsltdir/.
cp -rf par_$imin/run_dm_ann_general.sh $rsltdir/.
cp -rf submit_job_dm_ann.sh $rsltdir/.
cp -rf run_dm_ann_parallel.sh $rsltdir/.
###################################################

echo "finished!"
echo $start
echo `date`

rm -rf par_*

./mail_notify $mail $cluster $jobname