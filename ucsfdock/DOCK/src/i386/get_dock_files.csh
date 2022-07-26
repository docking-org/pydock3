#/bin/csh

grep '$base' $1 | grep -v then | sed 's#\s*$base/##' | awk '{print $1}'

