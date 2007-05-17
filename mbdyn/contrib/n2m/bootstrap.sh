#!/bin/sh

set -x
aclocal 
libtoolize --automake --force --copy
automake --foreign --add-missing --copy --force-missing
autoheader
autoconf

