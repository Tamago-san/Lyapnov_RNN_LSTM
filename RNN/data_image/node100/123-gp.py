#set terminal x11
#set grid
set tics font "Arial,10"
set title font "Arial,15"
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set y2label font "Arial,15"
set zlabel font "Arial,15"
#set ylabel "y"
set xlabel "G"
set title "NU = 0.5 "
set ylabel "LYAPNOV"
set y2label "ERROR"
set y2tics
set xrang[0.8:1.6]
set logscal y2
set ytics nomirror
plot 0 lc "#000" lt 0
replot "nu05-ly-err.tsv" using 1:2 axis x1y1  with linespoints pt 7 lc 2 lw 2 title "LYAPNOV"
replot "nu05-ly-err.tsv" using 1:3 axis x1y2  with linespoints pt 7 lc 15 lw 2 title "ERROR"
