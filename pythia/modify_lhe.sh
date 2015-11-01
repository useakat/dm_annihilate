#!/bin/bash 

mass=$1
mg5dir=$2

file=unweighted_events.lhe
outfile=unweighted_events_mod.lhe
pid=1000049

hmass=`echo "scale=5; $mass/2." | bc`
nhmass=`echo "scale=5; -1*$mass/2." | bc`

echo '<LesHouchesEvents version="1.0">' > $outfile

headerflag=0
initflag=0
writeflag=0
eventflag=0
eventstartflag=0

while read line;
do
    a=($line)

    if [ "${a[0]}" = "</slha>" ];then    
    	slhaflag=0
    	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "</header>" ];then    
    	headerflag=0
    	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "</init>" ];then    
    	initflag=0
    	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "</event>" ];then
    	eventflag=0
    	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "</LesHouchesEvents>" ];then
    	writeflag=0
    	echo ${a[@]} >> $outfile
    fi
    
    if [ $writeflag -eq 1 ];then
     	if [ $initflag -eq 1 ];then
     	    if [ "${a[0]}" = "$pid" ];then
     		a[0]=11
     		a[1]=-11
     		a[2]=$hmass
     		a[3]=$hmass
     	    fi
     	    echo ${a[@]} >> $outfile
    	elif [ $eventstartflag -eq 1 ];then
    	    a[0]=`expr ${a[0]} + 2`
    	    echo ${a[@]} >> $outfile
    	    echo "11 -1 0 0 0 0 0.00000000000E+00 0.00000000000E+00" $hmass $hmass "0.00000000000E+00 0. -1." >> $outfile
    	    echo "-11 -1 0 0 0 0 0.00000000000E+00 0.00000000000E+00" $nhmass $hmass "0.00000000000E+00 0. 1." >> $outfile
    	    eventstartflag=0
    	elif [ $eventflag -eq 1 ];then	    
    	    if [ "${a[1]}" = "-1" ]; then
    		a[1]=2
    		a[2]=1
    		a[3]=2
    	    elif [ "${a[2]}" = "1" ]; then
    		if [ "${a[3]}" = "0" ]; then
    		    a[2]=3
    		    a[3]=3
    		fi
    	    fi
    	    echo ${a[@]} >> $outfile
     	fi
    fi

    if [ "${a[0]}" = "<header>" ];then
	headerflag=1
	writeflag=1
	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "<slha>" ];then
	slhaflag=1
	echo ${a[@]} >> $outfile
	cat ../$mg5dir/Cards/param_card.dat >> $outfile
    elif [ "${a[0]}" = "<init>" ];then
	initflag=1
	echo ${a[@]} >> $outfile
    elif [ "${a[0]}" = "<event>" ];then
	eventflag=1
	eventstartflag=1
	echo ${a[@]} >> $outfile
    fi
done < $file