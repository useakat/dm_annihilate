#!/bin/bash
selfdir=$(cd $(dirname $0);pwd)

run=$1
run_mode=$2  # 1: DM decay  2: DM annihilation 
mass=$3  # in GeV
nevents=$4
channel=$5 

run_name=${run}

echo "start $run_name"
echo ""

cd pythia
if [ $run_mode -eq 1 ];then
    echo "Generic DM decay and hadronization with pythia..."
    ECM=`echo "scale=5; $mass" | bc`
elif [ $run_mode -eq 2 ];then
    echo "Generic DM annihilation and hadronization with pythia..."
    ECM=`echo "scale=5; 2*$mass" | bc`
fi
echo ""   
if [ $channel == "zz" ];then
    id1=23
    id2=$id1 
elif [ $channel == "ww" ];then
    id1=24
    id2=-24
elif [ $channel == "hh" ];then
    id1=25
    id2=$id1
elif [ $channel == "tautau" ];then
    id1=15
    id2=-15
elif [ $channel == "uubar" ];then
    id1=2
    id2=-2
elif [ $channel == "ddbar" ];then
    id1=1
    id2=-1
elif [ $channel == "ssbar" ];then
    id1=3
    id2=-3
elif [ $channel == "ccbar" ];then
    id1=4
    id2=-4
elif [ $channel == "bbbar" ];then
    id1=5
    id2=-5
elif [ $channel == "ttbar" ];then
    id1=6
    id2=-6
elif [ $channel == "check" ];then
    id1=0
    id2=$id1
elif [ $channel == "lhe" ];then
    id1=-1
    id2=-1
fi

sed -e "s/999999:addChannel .*/999999:addChannel = 1 1.00 101 $id1 $id2/" \
    -e "s/999999:all .*/999999:all = GeneralResonance void 1 0 0 $ECM 1. 0. 0. 0./" generic_resonance.cmnd > tmp.cmnd
mv tmp.cmnd generic_resonance.cmnd

    rm -rf *.dat 
    make hadron_dist_moroi
#    gfortran hadron_dist_pythia6.f pythia-6.4.28.f -o hadron_dist_pythia6
#    ./hadron_dist $mass $nevents
    ./hadron_dist_moroi $run_mode $mass $nevents $id1
#    ./hadron_dist_pythia6 $mass $nevents $id1    
cd ..

$selfdir/makedir.sh data 0
cd data
if [ -e $run_name ];then
    a=3
else
    mkdir $run_name
fi
mv ../pythia/*.dat $run_name/.
mv ../pythia/generic_resonance.cmnd $run_name/.
cd ..

echo "$run_name finished!"
echo ""