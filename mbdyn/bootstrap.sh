#!/bin/sh

set -x
aclocal 
libtoolize --force --copy
autoheader
for i in `find . -name 'Makefile.am'` ; do
	j=`echo $i | sed "s/\.am$//" | grep -v "^\./contrib"`
	if test "${j:+set}" = "set" ; then
		automake --foreign --add-missing --copy $j
	fi
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

#
# Copyright update
#
# for i in `grep -Erl 'Copyright.*-200[^3]' .`; do sed 's;\(Copyright.*\)200[0-2];\12004;' $i >x; mv -f x $i; done
#
