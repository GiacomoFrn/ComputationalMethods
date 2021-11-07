unset colorbox
set palette grey
stats "spinMat.txt" name "A"
set xrange [0:50]
set yrange [0:50]
set terminal gif animate delay 2
set output "IsingLattice2D.gif"
do for [i=0:int(A_blocks-1)]{ plot "spinMat.txt" index i matrix with image }