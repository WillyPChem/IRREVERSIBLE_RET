#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 20 


set output 'RET_EFFICIENCY.eps'
set xrange [0:1]
set key bottom right
set ylabel "Energy Transferred (eV)"
set xlabel '{/Symbol g}_A / {/Symbol g}_D'
plot 'SCAN_1D_VARY_GAMMA_ACCEPTOR_COUPLED.txt' u ($4/$3):($5*20) w l lw 3 title '20x Energy Transferred from Au', \
'SCAN_1D_Ag_VARY_GAMMA_ACCEPTOR_COUPLED.txt' u ($4/$3):5 w l lw 3 title 'Energy Transferred from Ag'


set output 'RET_EFFICIENCY_SCALED.eps'
set xrange [0:4e-5]
set ylabel "Energy Transferred (eV)"
set xlabel 'Scaled {/Symbol g}_A / {/Symbol g}_D'
plot 'SCAN_1D_VARY_GAMMA_ACCEPTOR_COUPLED.txt' u ($4/($3*42*42*$2*$2)):($5*20) w l lw 3 title '20x Energy Transferred from Au', \
'SCAN_1D_Ag_VARY_GAMMA_ACCEPTOR_COUPLED.txt' u ($4/($3*72*72*$2*$2)):5 w l lw 3 title 'Energy Transferred from Ag'

