#unset key
set grid
#set terminal windows enhanced
set output
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#カラーの用意   ex) ls 1 なら赤                               !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set style line 1 linecolor rgbcolor "#FF0000" # red
set style line 2 linecolor rgbcolor "#0000FF" # blue
set style line 3 linecolor rgbcolor "#808000" # olive
set style line 4 linecolor rgbcolor "#800080" # purple
set style line 5 linecolor rgbcolor "#008000" # green
set style line 6 linecolor rgbcolor "#00ffff" # aqua
set style line 7 linecolor rgbcolor "#000080" # navy
set style line 8 linecolor rgbcolor "#800000" # maroon

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#範囲、ラベルの指定　                                                 !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#set xrange[85000:90000]
#set xrange[2000000:3000000]
#set xlabel "パラメータ : b"
#set xlabel "step(dt=10^(-4), skip=10time)"
#set ylabel "LYAPNOV"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#使用ファイル                                                 !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot "lyapnov.dat" using 1:2 with line
#plot "lyapnov_end.dat" using 1:2 with lp ls 1 lw 2 pt 5 ps 2 title "b-lyapnov"
#plot "lyapnov_notstadard.dat" using 1:2 with line
#plot "Ly_cal.dat" using 1:6 with line
#write(30,*) step,lyap,abs_dr,r(NOW,2),r(NOW,3),dr(1)

#plot "Ly_skip.dat" using 3:4 with line
#write(30,*) step,abs_dr,r(NOW,1),r(NOW,2),dr(1)

#plot "time_series.dat" using 1:2 with line

#pause -1

