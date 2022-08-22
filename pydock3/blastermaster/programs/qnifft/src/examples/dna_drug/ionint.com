#!/bin/csh -f
#
# analyse phimaps
ionint <<***
$1
5
dapin
dapin
1
0
dapout
dapout
0
***
