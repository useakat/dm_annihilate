#!/bin/bash
run=$1
decay=$2
imkdir=$3
mail=$4
###### MODIFY HERE: running parameters #################
#output=width_mass.dat
jobname=dm_ann
submit_mode=0 # 0:serial submittion 1:parallel submission
job_system=bsub
que=l

if [ $decay == "ww" ];then
    min=100
    max=1000000
    imin=1
    ndiv=20
    ext=ww
elif [ $decay == "zz" ];then
    min=100
    max=1000000
    imin=1
    ndiv=20
    ext=zz
elif [ $decay == "hh" ];then
    min=100
    max=1000000
    imin=2
    ndiv=20
    ext=hh
elif [ $decay == "tautau" ];then
    min=10
    max=1000000
    imin=3
    ndiv=25
    ext=tautau
elif [ $decay == "uubar" ];then
    min=1
    max=1000000
    imin=2
    ndiv=30
    ext=uubar
elif [ $decay == "ddbar" ];then
    min=1
    max=1000000
    imin=2
    ndiv=30
    ext=ddbar
elif [ $decay == "bbbar" ];then
    min=1
    max=1000000
    imin=5
    ndiv=30
    ext=bbbar
fi
logflag=1

imax=`expr $ndiv + 1`
#mg5dir=grv_decay
######################################################
start=`date`
echo $start

i=$imin
while [ $i -le $imax ];do
    if [ $logflag -eq 1 ];then
	x=`echo "scale=5; e( (l($min)/l(10) +(l($max)/l(10) -l($min)/l(10))/$ndiv*($i-1))*l(10) )" | bc -l`
    else
	x=`echo "scale=5; $min +($max -$min)/$ndiv*($i-1)" | bc -l`
    fi

    job=${jobname}_$i
    if [ $i -le 15 ];then
	nevents=1000
#	nevents=100000
    else
	nevents=1000
#	nevents=100000
    fi
    ./submit_job_dm_ann.sh $job_system $que $i $job "./run_dm_ann_general.sh run_$i $x $nevents $decay > allprocess.log" $submit_mode

    i=`expr $i + 1`
done
n=$i

if [ $submit_mode -eq 1 ];then
    ./monitor
    # i=$imin
    # while [ $i -lt $n ];do
    # 	cd par_$i
    # 	if [ -e done.bjob$i ];then
    # 	    a=3
    # 	else
    # 	    ./bjob$i
    # 	fi
    # 	cd ..
    # 	i=`expr $i + 1`
    # done
fi

rsltdir=rslt_$run
makedir.sh $rsltdir $imkdir

#touch $rsltdir/np_sptrm_$ext.dat
#touch $rsltdir/nini_$ext.dat
#touch $rsltdir/Evis_$ext.dat

#echo "1, 141, 300, 0, 30, 1" >> $rsltdir/$output
#echo "%%%%%" >> $rsltdir/$output
x=$min
i=$imin
while [ $i -lt $n ];do
### MODIFY HERE for preparing result files ###############
#    grvinfo=`cat par_$i/grvinfo.dat`
#    Evis_tot=`cat par_$i/data/run_$i/Evis_tot.dat`
#    cat par_$i/grvinfo.dat >> $rsltdir/$output
#    echo $grvinfo $Evis_tot >> $rsltdir/$output
#    echo $grvinfo >> $rsltdir/$output
    if [ $i -eq $imin ];then
	cat par_$i/data/run_$i/np_sptrm.dat > $rsltdir/np_sptrm_$ext.dat
	cat par_$i/data/run_$i/nini.dat > $rsltdir/nini_$ext.dat
	cat par_$i/data/run_$i/Evis.dat > $rsltdir/Evis_$ext.dat
    else
	cat par_$i/data/run_$i/np_sptrm.dat >> $rsltdir/np_sptrm_$ext.dat
	cat par_$i/data/run_$i/nini.dat >> $rsltdir/nini_$ext.dat
	cat par_$i/data/run_$i/Evis.dat >> $rsltdir/Evis_$ext.dat
    fi
#    echo "%%%%%" >> $rsltdir/$outpu
#    echo " " >> $rsltdir/np_sptrm_$ext.dat
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

if [ $mail -eq 1 ];then
    bsub -q e -J cbscan -u takaesu@post.kek.jp nulljob.sh >/dev/null 2>&1
fi