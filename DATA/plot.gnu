#!/usr/bin/gnuplot
set terminal postscript enhanced color 'Helvetica' 22
set output 'Scat_MG_Compare.eps'
set ylabel 'Cross Section (m^2)'
set xlabel 'Energy (eV)'
plot 'AbsorptionSpectrum_Au_50.dat' u 1:3 w l lw 2 title 'Rad Decay', \
'AbsorptionSpectrum_3Au_50.dat' u 1:3 w l lw 2 dt 2 title 'Rad+NonRad Decay'

set output 'Abs_MG_Compare.eps'
set ylabel 'Cross Section (m^2)'
set xlabel 'Energy (eV)'
plot 'AbsorptionSpectrum_Au_50.dat' u 1:5 w l lw 2 title 'Rad Decay', \
'AbsorptionSpectrum_3Au_50.dat' u 1:5 w l lw 2 dt 2 title 'Rad+NonRad Decay'


set output 'Scat_Au_Compare.eps'
set ylabel 'Cross Section (m^2)'
set xlabel 'Energy (eV)'
plot 'AbsorptionSpectrum_Au_50.dat' u 1:2 w l lw 2 title 'Rad Decay', \
'AbsorptionSpectrum_3Au_50.dat' u 1:2 w l lw 2 dt 2 title 'Rad+NonRad Decay'

set output 'Abs_Au_Compare.eps'
set ylabel 'Cross Section (m^2)'
set xlabel 'Energy (eV)'
plot 'AbsorptionSpectrum_Au_50.dat' u 1:4 w l lw 2 title 'Rad Decay', \
'AbsorptionSpectrum_3Au_50.dat' u 1:4 w l lw 2 dt 2 title 'Rad+NonRad Decay'


