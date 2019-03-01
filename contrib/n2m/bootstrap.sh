#!/bin/sh

set -x
aclocal 
libtoolize --automake --force --copy
autoheader
automake --foreign --add-missing --copy --force-missing
autoconf

