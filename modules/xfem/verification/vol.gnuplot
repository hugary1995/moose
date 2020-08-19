set datafile separator ','
set tics font 'Verdana,24'
set key off

plot 'volume_change_out.csv' using 1:2 with linespoints lw 3 pt 7 ps 1
