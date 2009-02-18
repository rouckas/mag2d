set term x11
#set term postscript color
#set output "pot.ps"
set pm3d
unset surface
set multiplot
set yrange [0:0.015]
set xrange [0:0.065]
set ytics 0.005
set view map
set size ratio 0.5
splot "output/potential.dat" u 2:1:3
unset pm3d
set surface
splot "output/traj1.dat" u 2:1:(0) w lines
splot "output/traj2.dat" u 2:1:(0) w lines
unset multiplot

pause -1
