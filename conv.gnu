#
#   GNUPLOT v3.6 beta multiplot script file
#
plot 'conv.dat' using 0:1 with line title 'h',\
     'conv.dat' using 0:2 with line title 'hu',\
     'conv.dat' using 0:3 with line title 'hv'
#==========================================
pause -1 "Hit return to continue"

