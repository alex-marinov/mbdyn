#!/bin/sh

if test x"$MBDYN_PATH" = x ; then
	SCRIPT=/home/masarati/Lavoro/mbdyn/mbdyn-mfree/contrib/OpenDX-post/GeoOut.awk
else
	SCRIPT="$MBDYN_PATH/var/GeoOut.awk"
fi

ARGS=""

while test $# -gt 0 ; do
        case "$1" in
		-incr)
                        ARGS="$ARGS -vINCR=$2"
                        shift; shift ;;
		-start)
                        ARGS="$ARGS -vSTART=$2"
                        shift; shift ;;
                -skip)
                        ARGS="$ARGS -vSKIP=$2"
                        shift; shift ;;
                -stop)
                        ARGS="$ARGS -vSTOP=$2"
                        shift; shift ;;
		*)
			FILE=$1
			shift ;;
	esac
done

if ! [ -f "${FILE}.log" ] ; then
	echo "missing log file ${FILE}.log"
	exit 1
fi

if ! [ -f "${FILE}.mov" ] ; then
	echo "missing node data file ${FILE}.mov"
	exit 1
fi

echo "======" | cat "${FILE}.log" - "${FILE}.mov" | awk $ARGS -f "$SCRIPT"
