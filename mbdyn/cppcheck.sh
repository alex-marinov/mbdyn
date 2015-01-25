#!/bin/sh

#CPPCHECK=/home/masarati/tmp/cppcheck-1.46.1/cppcheck
CPPCHECK=cppcheck
SRCDIR=/home/masarati/Lavoro/mbdyn/mbdyn-src
BLDDIR=/home/masarati/Lavoro/mbdyn/mbdyn-devel

$CPPCHECK --enable=all --force \
	-I $BLDDIR/include/ \
	-I $SRCDIR/include/ \
	-I $SRCDIR/libraries/libann/ \
	-I $SRCDIR/libraries/libcolamd/ \
	-I $SRCDIR/libraries/libmbc \
	-I $SRCDIR/libraries/libmbmath \
	-I $SRCDIR/libraries/libmbutil/ \
	-I $SRCDIR/libraries/libmbwrap/ \
	-I $SRCDIR/libraries/libnaive/ \
	-I $SRCDIR/libraries/libobjs/ \
	-I $SRCDIR/libraries/liby12/ \
	-I $SRCDIR/mbdyn/ \
	-I $SRCDIR/mbdyn/aero/ \
	-I $SRCDIR/mbdyn/base/ \
	-I $SRCDIR/mbdyn/elec/ \
	-I $SRCDIR/mbdyn/hydr/ \
	-I $SRCDIR/mbdyn/struct/ \
	-I $SRCDIR/mbdyn/thermo/ \
	-I $SRCDIR/modules/ \
	$SRCDIR/libraries/ $SRCDIR/mbdyn/ $SRCDIR/utils/ $SRCDIR/modules/ \
	> x.log 2> x.err

