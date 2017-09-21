#!/usr/bin/gnuplot
set terminal postscript enhanced color 'Helvetica' 22
set output 'Pop_Superposition.eps'
set ylabel 'Sink Population'
set xlabel 'Time (atomic units)'
plot 'PopulationMG_50_halfAg.dat' u 1:3 w l lw 3 title '{/Symbol g}_A = 0.5{/Symbol g}_{Ag}', \
'PopulationMG_50_1Ag.dat' u 1:3 w l lw 3 title '{/Symbol g}_A = {/Symbol g}_{Ag}', \
'PopulationMG_50_2Ag.dat' u 1:3 w l lw 3 title '{/Symbol g}_A = 2{/Symbol g}_{Ag}', \
'PopulationMG_50_3Ag.dat' u 1:3 w l lw 3 title '{/Symbol g}_A = 3{/Symbol g}_{Ag}'

set xrange [0:2000]
set output 'DP_Superposition.eps'
set ylabel 'Dipole Moment (a.u.)'
plot 'DipoleMomentMG_50_halfAg.dat' u 1:2 w l lw 3 title '{/Symbol g}_A = 0.5{/Symbol g}_{Ag}', \
'DipoleMomentMG_50_1Ag.dat' u 1:2 w l lw 3 title '{/Symbol g}_A = {/Symbol g}_{Ag}', \
'DipoleMomentMG_50_2Ag.dat' u 1:2 w l lw 3 title '{/Symbol g}_A = 2{/Symbol g}_{Ag}', \
'DipoleMomentMG_50_3Ag.dat' u 1:2 w l lw 3 title '{/Symbol g}_A = 3{/Symbol g}_{Ag}'

