#!/bin/csh -f
#
# analyse phimaps
ionintm193 <<***
0.,$1,0.,$1
5
dapin
dapin
1
0
dapout
dapout
0
***
