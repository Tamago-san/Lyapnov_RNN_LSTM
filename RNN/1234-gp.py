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
#set ylabel "LYAPNOV"
set ylabel "ERROR"
#set y2tics
set logscal y
set ytics nomirror
set xrange[0:200]
#set yrange[:1]
rc_err=0.024269096781920487


#! データをプロットする。
#! 「axis x1y1」＝左側のY軸を使用して描画
#! 「axis x1y2」＝右側のY軸を使用して描画
plot "./data_image/node005/lyapnov_end.dat" using 1:3 axis x1y1  with l  lw 2 title "node-005"
replot "./data_image/node010/lyapnov_end.dat" using 1:3 axis x1y2  with l  lw 2 title "node-010"
replot "./data_image/node015/lyapnov_end.dat" using 1:3 axis x1y2  with l  lw 2 title "node-015"
replot "./data_image/node020/lyapnov_end.dat" using 1:3 axis x1y2  with l  lw 2 title "node-020"
replot "./data_image/node050/lyapnov_end.dat" using 1:3 axis x1y2  with l  lw 2 title "node-050"
replot "./data_image/node100/lyapnov_end.dat" using 1:3 axis x1y2  with l  lw 2 title "node-100"


replot rc_err axis x1y2  with lines dt (5,5) lc 15 lw 2 title "RC(node100)-ERROR"
#!write(41,*) iepoch,lyap,ERR