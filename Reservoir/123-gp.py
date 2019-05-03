#set terminal x11
set grid
set tics font "Arial,10"
set title font "Arial,15"
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set y2label font "Arial,15"
set zlabel font "Arial,15"
set title "NU=0.00"
set ylabel "lyapnov"
set xlabel "G"
set zlabel "z"
#set xrange[100:500]
set yrange[-2.0:2.0]
#set y2tics
##! 右側のY軸の範囲を設定する
#set y2range [0:1.0]
#set y2label "ERROR"
#set yrange [-3.70:-3.67]


#!============================================================
#!============================================================
#plot "./data_out/lyapnov_end_2.dat" using 1:3 axis x1y1  with linespoints title "LYAPNOV(NU=0.00)" lc 2
plot "./data_out/lyapnov_end.dat" using 1:3 axis x1y1  with linespoints title "LYAPNOV(NU=0.00)" lc 2 lw 2
#replot "./data_out/lyapnov_end.dat" using 1:4 axis x1y2  with linespoints title "ERROR" lc 15
#replot "./data_out/lyapnov_end_2.dat" using 1:5 axis x1y1  with linespoints title "LYAPNOV(NU=1.00)" lc 2
#replot "./data_out/lyapnov_end_2.dat" using 1:6 axis x1y2  with linespoints title "ERROR(NU=1.00)" lc 15
#plot "./data_renban/lyapnov.100000" using 1:2 axis x1y1  with line title "LYAPNOV G=1.00"
#replot "./data_renban/lyapnov.100000" using 1:3 axis x1y2  with line title "LYAPNOV G=1.00"
#plot "./data_renban2/rc_out.040000" using 1 axis x1y1  with line title "LYAPNOV G=1.20"
#replot "./data_renban/lyapnov.200000" using 1:2 axis x1y1  with line title "LYAPNOV G=2.00"
#replot "./data_renban/lyapnov.100000" using 1:3 axis x1y2  with line title "ERROR"
#replot "./data_renban2/rc_out.0480" using 3 with line



#!============================================================
#!============================================================
#plot "./data_renban2/rc_out.140200" using 1 with line
#replot "./data_renban2/rc_out.140200" using 2 with line
#replot "./data_renban2/rc_out.140200" using 3 with line
#replot "./data_renban2/rc_out.140200" using 4 with line
#plot "./data_renban2/rc_out.080100" using 4 with line
#replot "./data_renban2/rc_out.100000" using 5 with line
#
#!============================================================
#!============================================================
#splot "./data_renban2/rc_out.005" using 1:2:3 with line
#plot "./data_renban2/rc_out.0060" using 1 with line
#replot "./data_renban2/rc_out.0100" using 2 with line
#replot "./data_renban2/rc_out.0100" using 3 with line
