#!/bin/csh -f
#
# analyse phimaps
ionint <<***
$1
5
sphin
sphin
1
0
sphout
sphout
0
***
