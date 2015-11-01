#!/bin/bash
selfdir=$(cd $(dirname $0);pwd)

run=$1
mass=$2  # in GeV
nevents=$3
mg5dir=$4
chi0mass=$5

run_name=${run}

echo "start $run_name"
echo ""

cd $mg5dir
#sed -e "s/NNEVENTS/$nevents/" Cards/run_card_temp.dat > Cards/run_card.dat
#sed -e "s/MMGRV/$mass/" Cards/param_card_HSUSY_noHdecay.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_HSUSY_nodecay.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_HSUSY.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_Hsfermion_noHdecay.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_Hsfermion.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_case1.dat > Cards/param_card.dat
#    sed -e "s/MMGRV/$mass/" Cards/param_card_splitsusy.dat > Cards/param_card.dat
cd ..

nopsflag=0
if [ $mass -gt $chi0mass ];then
    echo "generating events..."
    cd $mg5dir
    echo ""
    line=""
    i=1
    while [ -z "$line" ];do 
	echo "event generation $i"
#    ./bin/generate_events $run_name -f | tee ./mg5.log
	touch mg5.log
	./bin/generate_events $run_name -f > mg5.log
	line=`grep "Width :" mg5.log`
	rm -rf RunWeb
	sleep 10
	i=`expr $i + 1`
#    if [ $i -eq 4 ];then
#	line=3.24
#	echo "forced to finish event generation"
#	echo ""
#    fi
    done
    cd ..
    width=${line#*:}
    width=${width%+-*}
    exp=${width#*e}
    rr=${width%e*}
    exp=`expr -25 - $exp` 
    taur=`echo "scale=5; 6.667/$rr" | bc`
    tau=${taur}e$exp
    echo $mass $width $tau > grvinfo.dat
else
    echo "m3/2 is smaller than m_chi0. skip event generation..."
    nopsflag=1
    echo $mass 0.0 1e99 > grvinfo.dat
fi

cd pythia
#rm -rf decayLHE_Evis decayLHE decayLHE_nops
if [ $nopsflag -eq 0 ];then
    echo "hadronizing with pythia..."
    echo ""
    cp -rf ../$mg5dir/Events/$run_name/unweighted_events.lhe.gz .
    gunzip -f unweighted_events.lhe.gz
    ./modify_lhe.sh $mass $mg5dir
    make decayLHE_Evis
    ./decayLHE_Evis $nevents
    make decayLHE
    ./decayLHE $nevents
else
    echo "skip pythia ..."
    echo "0.0" > Evis_tot.dat
    make decayLHE_nops
    ./decayLHE_nops
fi
cd ..

cd data
if [ -e $run_name ];then
    a=3
else
    mkdir $run_name
fi
mv ../pythia/*.dat $run_name/.
echo "$run_name finished!"
echo ""