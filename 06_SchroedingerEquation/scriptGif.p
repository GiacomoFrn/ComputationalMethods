set terminal gif animate delay 1
set xrange [A_min_x:A_max_x]
set yrange [A_min_y-0.01:A_max_y+0.01]
set xlabel "position [0:L]"
set ylabel "|psi|^2"
set terminal gif animate delay 2
set output "psi2.gif"
do for [i=0:int(A_blocks-1)]{ plot "phi2.txt" index i w l lw 2.5 }