
set multiplot layout 3,1
set grid

set xlabel "t"
set ylabel "Lyapnov"
#set yr[0:1.5]
#set xr[780:800]
plot "lyapnov.dat" using 1:2 with line lc rgb "#ff0000"
#plot "./time_series/Lyapnov_time_series_P.dat" using 1:2 w dots lc rgb "#ff0000"


set xlabel "t"
set ylabel "abs_dV"
#set yr[0:1.5]
#set xr[780:800]
plot "a.dat" using 1:3 with line lc rgb "#ff0000"


set xlabel "t"
set ylabel "X"
#set yr[-20:20]
#set xr[780:800]
#plot "a.dat" using 1:4 with line lc rgb "#ff0000"

unset multiplot

