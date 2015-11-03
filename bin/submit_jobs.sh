#!/bin/bash
selfdir=$(cd $(dirname $0);pwd)
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: submit_job.sh [job_system] [que] [it] [jobname] [command]"
    echo ""
    exit
fi

cluster=$1
que=$2
i=$3
jobname=$4
command="$5"
submit_mode=$6 # 0:serial 1:parallel
work_dir=$7

njob=bjob

dir=par_$i
cd $work_dir
mkdir $dir
cd $dir

echo "#!/bin/bash" > $njob$i
echo "" >> $njob$i
if [ $cluster == "icrr" ];then
    echo '#------ pjsub option --------#' >> $njob$i
    echo '#PJM -L "rscunit=common"' >> $njob$i
    echo '#PJM -L "rscgrp=B"' >> $njob$i
#    echo '#PJM -L "rscunit=group"' >> $njob$i
#    echo '#PJM -L "rscgrp=th"' >> $njob$i
    echo '#PJM -L "vnode=1"' >> $njob$i
    echo '#PJM -L "vnode-core=1"' >> $njob$i
#    echo '#PJM -L "vnode-mem=3Gi"' >> $njob$i
#    echo '#PJM -L "elapse=00:15:00"' >> $njob$i # A:<3h B:<24h C:<1week th:no limit
fi
echo '#------- Program execution -------#' >> $njob$i
echo "date >allprocess.log" >> $njob$i
echo "rm -rf wait.${njob}$i" >> $njob$i
echo "touch run.${njob}$i" >> $njob$i
echo "cp -rf $selfdir/../pythia ." >> $njob$i
#echo "cp -rf $selfdir ." >> $njob$i
echo "mkdir data" >> $njob$i
echo "$command >>allprocess.log 2>&1" >> $njob$i
echo "rm -rf run.${njob}$i" >> $njob$i
echo "touch done.${njob}$i" >> $njob$i
echo "rm -rf pythia" >> $njob$i

chmod +x $njob$i
touch wait.$njob$i

if [ $submit_mode -eq 0 ];then
    echo "job$i launched"
    ./$njob$i 1>/dev/null
    echo "job$i finished"
    echo
else
    if [ $cluster == "kekcc" ];then
	bsub -q $que -J $jobname ./$njob$i
    elif [ $cluster == "icrr" ];then
	pjsub -N $jobname -j $njob$i
    fi
    echo ""
fi
