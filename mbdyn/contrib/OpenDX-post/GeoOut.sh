#!/bin/sh

FILE=${1:?"usage: GeoOut.sh <file>"}
SCRIPT=/home/masarati/Lavoro/mbdyn/mbdyn-mfree/contrib/OpenDX-post/GeoOut.awk

if ! [ -f "${FILE}.log" ] ; then
	echo "missing log file ${FILE}.log"
	exit 1
fi

if ! [ -f "${FILE}.mov" ] ; then
	echo "missing node data file ${FILE}.mov"
	exit 1
fi

echo "======" | cat "${FILE}.log" - "${FILE}.mov" | awk -f "$SCRIPT"
