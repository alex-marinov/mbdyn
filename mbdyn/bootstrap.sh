#!/bin/sh

set -x
aclocal 
libtoolize --force --copy
autoheader
for i in `find . -name 'Makefile.am'` ; do
	j=`echo $i | sed "s/\.am//"`
	automake --foreign --add-missing --copy $j
done
if test -d contrib ; then
	for i in `find contrib -name 'bootstrap.sh'` ; do
		dir=`echo $i | sed "s/\(.*\)\/bootstrap\.sh/\1/"`
		olddir=`pwd`
		echo "=> entering subdir '$dir'..."
		cd $dir
		sh bootstrap.sh
		echo "<= reverting to dir '$olddir'..."
		cd $olddir
	done
fi
autoconf

