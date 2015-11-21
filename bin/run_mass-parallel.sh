#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: run_mass_parallel.sh [run name] [run mode] [channel] [imkdir] [mail]"
    echo ""
    echo "run name: run name"
    echo "run mode: 1:DM decay 2:DM annihilation"
    echo "run mode: ww/zz/tautau/uubar/ddbar/ccbar/ssbar/bbbar/ttbar"
    echo "  imkdir: 0:do not delete an existing dir with the same name when making a dir  1:delete the existing dir"
    echo "    mail: 0:no mail notification  1:mail notification when all jobs are finished"
    echo ""
    exit
fi
selfdir=$(cd $(dirname $0);pwd)
start=`date`
echo $start

run=$1       
run_mode=$2  
channel=$3   
imkdir=$4    
mail=$5      
###### MODIFY HERE: running parameters #################
output=hadron_dist_$channel.dat
submit_mode=1
#cluster=kekcc
cluster=icrr
que=l

# working space for jobs on a remote server
if [ $cluster == "icrr" ];then
    work_dir=/disk/th/work/takaesu/$run
    mkdir $work_dir
elif [ $cluster == "kekcc" ];then
    work_dir=./
fi

#masses=("0" "10" "20" "30" "10000" "100000" "1000000")  # Dark matter masses
masses=("0" "30" "100" "1000" "10000" "100000" "1000000")  # Dark matter masses
#masses=("0" "60" "200" "2000" "20000" "200000" "2000000")  # For check with annihilation case
nevents=("0" "1750000000" "1630000000" "1460000000" "1270000000" "1140000000" "1000000000")  # Event numbers to be generated
#nevents=("0" "10000000" "10000000" "10000000" "10000000" "10000000" "10000000")  # generated events in tautau
#nevents=("0" "10000000" "10000000" "10000000" "10000000" "10000000" "10000000")  # generated events in ww
#nevents=("0" "40000" "30000" "10000" "3000" "500" "100") # optimized ratio for bbbar
#nevents=("0" "1000000" "500000" "140000" "30000" "6000" "1000") # optimized ratio for qqbar
imin=1
imax=6
################# Main Code #####################################
bin_dir=$selfdir
if [ $run_mode -eq 1 ];then
    job=decay_general_mass$RANDOM
elif [ $run_mode -eq 2 ];then
    job=ann_general_mass$RANDOM
else 
    "ERROR: Invalid run_mode is entered! run_mode = 1 or 2."
fi

i=1
while [ $i -le $imax ];do 
    $bin_dir/submit_jobs.sh $cluster $que $i $job "$bin_dir/run_general.sh run_$i $run_mode ${masses[$i]} ${nevents[$i]} $channel" $submit_mode $work_dir
    i=`expr $i + 1`
done


if [ $submit_mode -eq 1 ];then
    $bin_dir/monitor $work_dir
fi
if [ $cluster == "icrr" ];then
    mv $work_dir/* .
    rm -rf $work_dir
fi

rsltdir=$bin_dir/../results/rslt_$run
$bin_dir/makedir.sh $rsltdir $imkdir
rm -rf $rsltdir/$output

echo "1, 141, 300, 0, 30, 1" > $rsltdir/$output
echo "%%%%%" >> $rsltdir/$output
i=$imin
while [ $i -le $imax ];do
### MODIFY HERE for preparing result files ###############
    mass=${masses[$i]}
    ext=${channel}_${mass}
    Evis_tot=`cat par_$i/data/run_$i/Evis_tot.dat`
    echo $mass 0 0 $Evis_tot >> $rsltdir/$output
    cat par_$i/data/run_$i/Edist.dat >> $rsltdir/$output
    echo "%%%%%" >> $rsltdir/$output
# Prepare outputfiles and write the first line
    cat par_$i/data/run_$i/np_sptrm.dat > $rsltdir/np_sptrm_$ext.dat
    cat par_$i/data/run_$i/nini.dat > $rsltdir/nini_$ext.dat
    cat par_$i/data/run_$i/Evis.dat > $rsltdir/Evis_$ext.dat
    cat par_$i/data/run_$i/np_sptrm_norm.dat > $rsltdir/Ekin_z_$ext.dat
    cp -rf par_$i/allprocess.log $rsltdir/allprocess_$i.log
###################################################
    i=`expr $i + 1`
done
### MODIFY HERE for saving files relatee to this run
cp -rf $bin_dir/run_general.sh $rsltdir/.
cp -rf $bin_dir/submit_jobs.sh $rsltdir/.
cp -rf $bin_dir/run_mass-parallel.sh $rsltdir/.
###################################################
echo "finished!"
rm -rf par_*

$bin_dir/mail_notify $mail $cluster $job

echo $start
echo `date`