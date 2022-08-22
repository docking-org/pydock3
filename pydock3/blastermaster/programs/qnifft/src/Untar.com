#!/bin/csh -f
# extract qnifft
#
mkdir delphi
cd delphi
tar xvf ../qnifft.tar
