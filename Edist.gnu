if (exist("ii")==0 || ii<0) ii=1
###################### Options ###########################################
set logscale x
set logscale y
set format x '10^{%L}'
#set format x '%L'
set format y '10^{%L}'
#set xtics (2,3,4,5,6)
#set ytics (1,10,1E2,1E3,1E4,1E5,1E6,1E7,1E8,1E9,1E10)
#set tics scale 2
#set key at 1.0E3,1.0E7 samplen 2
#set xrange [1:7]
#set yrange [1E-5:2E8]
####################### Definitions ######################################
#file1 = 'rslt_10k/np_sptrm_ww.dat'
file1 = 'rslt_100k_2/np_sptrm_ww.dat'
c1 = 'red'
c2 = 'blue'
c3 = '#006400' # dark green
c4 = 'purple'
c5 = '#ff33ff'
c6 = '#cc6600' # dark orange
##########################################################################
set terminal postscript eps enhanced 'Times-Roman' color 20
set grid
set key spacing 1.5 samplen 2
#set multiplot

set output 'plots/Edist_'.ii.'.eps' 
start = 1 +802*(ii-1) 
mass = int(exp((log10(100) +(log10(1000000) -log10(100))/20*(ii-1))*log(10)))
end = start +800
set xrange [1E-2:1000000]
set yrange [1E-5:0.1]
#set yrange [0:10]
set title '{/=28 Edist}'
set label 'm3/2 = '.mass.' GeV' at graph 0.05, graph 0.92
set xlabel '{/=24 Ekin [GeV]}'
set ylabel '{/=24 Number of particles / event}' offset 0,0
#################### plot ##########################################
plot \
file1 every :::ii::ii u 2:3 title "neutron" w l lt 1 lw 3 lc rgb c1, \
file1 every :::ii::ii u 2:4 title "proton" w l lt 1 lw 3 lc rgb c2, \
file1 every :::ii::ii u 2:5 title "pi+" w l lt 1 lw 3 lc rgb c3, \
file1 every :::ii::ii u 2:6 title "pi-" w l lt 1 lw 3 lc rgb c4, \
file1 every :::ii::ii u 2:7 title "K+" w l lt 1 lw 3 lc rgb c5, \
file1 every :::ii::ii u 2:8 title "K-" w l lt 1 lw 3 lc rgb c6, \
file1 every :::ii::ii u 2:9 title "KLong" w l lt 2 lw 3 lc rgb c1, \
file1 every :::ii::ii u 2:10 title "anti-neutron" w l lt 2 lw 3 lc rgb c2, \
file1 every :::ii::ii u 2:11 title "anti-proton" w l lt 2 lw 3 lc rgb c3
###########################################################################

#set nomultiplot
unset label

if (ii<20) pause 0.1; ii=ii+1; \
reread
ii=-1

reset

