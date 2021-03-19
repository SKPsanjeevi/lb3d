#!/bin/bash

for f in *.f; do TF=`echo $f | sed -e 's/f$/F90/'`; echo "$f => $TF" ; cat $f | sed -e 's/^c/!/' | sed -e '1~2 N;s/\(.*\)\n\([ \t]*\)\$/\1 \&\n\2\&/g' | sed -e '2~2 N;s/\(.*\)\n\([ \t]*\)\$/\1 \&\n\2\&/g' > $TF; done
