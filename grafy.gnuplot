set pm3d
unset surface
set multiplot
set xrange [0:0.015]
set yrange [0:0.04]
set xtics 0.005
set view map
set size ratio 2
splot "output/out.dat"
unset pm3d
set surface
splot "output/traj.dat" u 1:2:(0) w lines
unset multiplot

pause -1

