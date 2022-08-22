#!/bin/csh
#
grep -i 'net assigned' $1
grep '0.5\*sum(' $1
grep 'rho.phi                    i:' $1
grep 'dPI                        i:' $1
