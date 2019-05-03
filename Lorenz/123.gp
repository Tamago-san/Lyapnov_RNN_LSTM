#set terminal x11
set grid
set tics font "Arial,10"
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set zlabel font "Arial,15"
set ylabel "y"
set xlabel "x"
set zlabel "z"
#splot "output_Runge_Lorenz.dat" using 2:3:4 with line
#splot "output_001_40_x7y1z6.dat" using 1:2:3 with line
plot "output_Runge_Lorenz.dat" using 2:3 with line
