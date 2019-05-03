set pm3d
set pm3d map
set ticslevel 0
set cbrange[:1.0]
set yrange[0.1:]
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set zlabel font "Arial,15"
set ylabel "NU"
set xlabel "G"
#set logscal cb
#set logscale cb
#set zlabel "z"
#set pm3d interpolate 5, 5
#set palette defined (-1 "blue", 0 "white", 1 "red")
#set palette defined (0 "blue", 1 "white")
#set palette rgbformulae 22, 13, -31
#splot "./data_out/lyapnov_end_ly.dat" with pm3d
splot "./data_out/lyapnov_end_err.dat" with pm3d
#splot "./data_out_g_nu/Book1.dat" with pm3d
#splot "./data_out/error.tsv" with pm3d
#splot "./data_out/ly_g-00-34_nu-00-30.tsv" with pm3d
#splot "./data_out/err_g-00-34_nu-00-30.tsv" with pm3d