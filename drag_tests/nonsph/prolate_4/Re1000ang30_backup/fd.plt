# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set   autoscale                        # scale axes automatically
# unset log                              # remove any log-scaling
# unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Force"
set xlabel "Timestep"
set ylabel "Fd"
# set key 0.01,100
# set label "Yield Point" at 0.003,260
# set arrow from 0.0028,250 to 0.003,280
# set xr [0.0:0.022]
# set yr [0:325]
set terminal png
set output 'xyz.png'
plot "./Production/md-cfg_out_p********-0430869971.asc" using 1:($4/(0.5*pi*23.81**2*0.09**2)) title 'Column' with lines
