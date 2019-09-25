set title "Convergence process"
set xlabel "Iterations"
set ylabel "P at selected points"

set style line 1 lt 1 lw 2 lc rgb "#FF0000"
set style line 2 lt 1 lw 2 lc rgb "#00C000"
set style line 3 lt 1 lw 2 lc rgb "#0000FF"
set style line 4 lt 1 lw 2 lc rgb "#FF9C20"
set style line 5 lt 1 lw 2 lc rgb "#9400D3"

set mytics 10
set grid
set grid mytics


plot  "pressProbes/0/p" using 1:2 with line ls 1 title 'x=0.5'

pause mouse
reread

