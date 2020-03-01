set xrange[-0.1:1.1]
set yrange[-0.1:1.1]
set size ratio 1
plot "output.dat" u 1:2:($3/20):($4/20):(sqrt($3*$3+$4*$4)) w vectors lc palette ti ""