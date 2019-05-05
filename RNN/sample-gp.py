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
#set logscal y2
set ytics nomirror
set xrange[0:100]
#set xrange[20:35]
set yrange[:0.3]
rc_err=0.07140703008014705
rc_lyapnov=-0.656803

#! 右側のY軸を使用可能にする
#set y2tics
#! 右側のY軸の範囲を設定する
#set y2range [0:2.2]
#set yrange [-3.70:-3.67]

#plot "./data_image/node010/sample_num/lyapnov_end_01.dat" using 1:2 axis x1y1  with linespoints pt 7 lc 2 lw 2 title "RNN-LYAPNOV"
#replot rc_lyapnov axis x1y1  with lines dt (5,5) lc 2 lw 2 title "RC -LYAPNOV"
#replot "./data_image/node010/sample_num/lyapnov_end_01.dat" using 1:3 axis x1y2  with linespoints pt 7 lc 15 lw 2 title "RNN-ERROR"
#replot rc_err axis x1y2  with lines dt (5,5) lc 15 lw 2 title "RC -ERROR"

plot "./data_image/node010/sample_num/lyapnov_end_sample.0001" using 1:3 axis x1y2 with linespoints pt 7 lc 15 lw 2 title "1"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0002" using 1:3 axis x1y2 with linespoints pt 7 lc 1 lw 2 title "0002"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0003" using 1:3 axis x1y2 with linespoints pt 7 lc 2 lw 2 title "0003"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0004" using 1:3 axis x1y2 with linespoints pt 7 lc 3 lw 2 title "0004"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0005" using 1:3 axis x1y2 with linespoints pt 7 lc 4 lw 2 title "0005"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0006" using 1:3 axis x1y2 with linespoints pt 7 lc 5 lw 2 title "0006"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0007" using 1:3 axis x1y2 with linespoints pt 7 lc 6 lw 2 title "0007"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0008" using 1:3 axis x1y2 with linespoints pt 7 lc 7 lw 2 title "0008"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0009" using 1:3 axis x1y2 with linespoints pt 7 lc 8 lw 2 title "0009"
replot "./data_image/node010/sample_num/lyapnov_end_sample.0011" using 1:3 axis x1y2 with linespoints pt 7 lc 9 lw 2 title "10"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0011" using 1:3 axis x1y2 with linespoints pt 7 lc 10 lw 2 title "0011"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0012" using 1:3 axis x1y2 with linespoints pt 7 lc 1 lw 2 title "0012"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0013" using 1:3 axis x1y2 with linespoints pt 7 lc 2 lw 2 title "0013"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0014" using 1:3 axis x1y2 with linespoints pt 7 lc 3 lw 2 title "0014"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0015" using 1:3 axis x1y2 with linespoints pt 7 lc 4 lw 2 title "0015"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0016" using 1:3 axis x1y2 with linespoints pt 7 lc 5 lw 2 title "0016"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0017" using 1:3 axis x1y2 with linespoints pt 7 lc 6 lw 2 title "0017"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0018" using 1:3 axis x1y2 with linespoints pt 7 lc 7 lw 2 title "0018"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0019" using 1:3 axis x1y2 with linespoints pt 7 lc 8 lw 2 title "0019"
replot "./data_image/node010/sample_num/lyapnov_end_sample.0019" using 1:3 axis x1y2 with linespoints pt 7 lc 6 lw 2 title "200"
replot "./data_image/node010/sample_num/lyapnov_end_sample.0030" using 1:3 axis x1y2 with linespoints pt 7 lc 2 lw 2 title "300"
#replot "./data_image/node010/sample_num/lyapnov_end_sample.0040" using 1:3 axis x1y2 with linespoints pt 7 lc 3 lw 2 title "400"
replot "./data_image/node010/sample_num/lyapnov_end_sample.0050" using 1:3 axis x1y2 with linespoints pt 7 lc 4 lw 2 title "500"

#replot rc_err axis x1y2  with lines dt (5,5) lc 15 lw 2 title "RC -ERROR"