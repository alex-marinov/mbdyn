#!/bin/sh

set -x
aclocal 
libtoolize --force --copy
autoheader
for i in `find . -name 'Makefile.am'` ; do
	j=`echo $i | sed "s/\.am//"`
	automake --foreign --add-missing --copy $j
done
autoconf

