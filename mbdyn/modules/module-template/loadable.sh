#!/bin/sh
SRC=../../src
UTILS=../../utils
SOURCE=module-template.so
LOADABLE=$UTILS/loadable
DEMANGLED=$SRC/demangled.h
MANGLEDTMP=mangled.tmp
MANGLED=mangled.h

if !(test -f $SOURCE) ; then
	echo "build "$SOURCE" first";
	exit 1;
fi
if !(test -f $LOADABLE) ; then
	echo "build "$LOADABLE" first";
	exit 1;
fi

nm $SOURCE | grep " T " | awk '{print $3}' > $MANGLEDTMP
$LOADABLE $DEMANGLED $MANGLEDTMP $MANGLED
mv $MANGLED $SRC
rm $MANGLEDTMP
