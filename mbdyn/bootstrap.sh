#!/bin/sh
# $Header$

set -x
aclocal 
libtoolize --automake --force --copy
automake --foreign --add-missing --copy --force-missing
autoheader
autoconf

if test -d contrib ; then
	for i in `find contrib -name 'bootstrap.sh'` ; do
		dir=`echo $i | sed "s/\(.*\)\/bootstrap\.sh/\1/"`
		olddir=`pwd`
		echo "=> entering subdir '$dir'..."
		( cd $dir && sh bootstrap.sh )
		echo "<= leaving subdir '$dir'; reverting to dir '$olddir'..."
	done
fi

