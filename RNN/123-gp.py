#set terminal x11
#set grid
set tics font "Arial,10"
set title font "Arial,15"
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set y2label font "Arial,15"
set zlabel font "Arial,15"
#set ylabel "y"
set xlabel "epoch"
set title "RNN"
set ylabel "LYAPNOV"
set y2label "ERROR"
set y2tics
set logscal y2
set ytics nomirror
set xrange[0:100]
set yrange[:0.3]
rc_err=0.07140703008014705
rc_lyapnov=-0.656803
#splot "./data/output_Runge_Lorenz.dat" using 2:3:4 with line
#plot "./data/output_Runge_Lorenz.dat" using 1:2 with line
#replot "./data/output_Runge_Lorenz.dat" using 1:3 with line


#plot "./data_out/output_traning1.dat" using 1:2 with line
#plot "./data_out/rc_out.dat" using 1:3 with line
#replot "./data_out/rc_out.dat" using 1:4 with line
#plot "tmp2.dat" using 1 with line
#plot "./data_out/lyapnov.dat" using 1:2 with line

#plot "tmp2.dat" using 1 with line
#replot "tmp2.dat" using 2 with line
#replot "tmp2.dat" using 3 with line


#! 右側のY軸を使用可能にする
#set y2tics
#! 右側のY軸の範囲を設定する
#set y2range [0:2.2]
#set yrange [-3.70:-3.67]
plot "./data_out/lyapnov_end.dat" using 1:2 axis x1y1  with linespoints pt 7 lc 2 lw 2 title "RNN-LYAPNOV"
replot rc_lyapnov axis x1y1  with lines dt (5,5) lc 2 lw 2 title "RC -LYAPNOV"
replot "./data_out/lyapnov_end.dat" using 1:3 axis x1y2  with linespoints pt 7 lc 15 lw 2 title "RNN-ERROR"
replot rc_err axis x1y2  with lines dt (5,5) lc 15 lw 2 title "RC -ERROR"

#! データをプロットする。
#! 「axis x1y1」＝左側のY軸を使用して描画
#! 「axis x1y2」＝右側のY軸を使用して描画
#plot "./data_image/node005/lyapnov_end.dat" using 1:2 axis x1y1  with linespoints pt 7 lc 2 lw 2 title "RNN-LYAPNOV"
#replot rc_lyapnov axis x1y1  with lines dt (5,5) lc 2 lw 2 title "RC -LYAPNOV"
#replot "./data_image/node005/lyapnov_end.dat" using 1:3 axis x1y2  with linespoints pt 7 lc 15 lw 2 title "RNN-ERROR"
#replot rc_err axis x1y2  with lines dt (5,5) lc 15 lw 2 title "RC -ERROR"
#!write(41,*) iepoch,lyap,ERR