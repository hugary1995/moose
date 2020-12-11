set datafile separator ','
set yrange [0:]

plot 'min_ev.csv' using (-log($1)):2 with linespoints pt 7

pause -1
