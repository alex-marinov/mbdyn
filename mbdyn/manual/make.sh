#!/bin/sh

VERSIONFILE=../build/version

DATE=`date "+%a %b %e, %Y"`
if [ -r $VERSIONFILE ] ; then
	VERSION="`cat $VERSIONFILE`"
	VV="-$VERSION"
fi
LOG="make.log"
IN="input"
OUT="mbdyn-input$VV"

rm -rf "$IN" "$OUT" "$LOG"

echo "making dvi..."
which latex > /dev/null 2>&1
RC=$?
if [ "$RC" != 0 ] ; then
	echo "unable to find latex"
	exit 1
fi
latex "$IN" >> "$LOG" 2>&1
which bibtex > /dev/null 2>&1
RC=$?
if [ "$RC" = 0 ] ; then
	bibtex "$IN" >> "$LOG" 2>&1
else
	echo "warning: no bibtex available..."
fi
latex "$IN" >> "$LOG" 2>&1
latex "$IN" >> "$LOG" 2>&1

echo "making ps..."
which dvips > /dev/null 2>&1
RC=$?
if [ "$RC" != 0 ] ; then
	echo "unable to find dvips..."
else
	dvips -o "$OUT.ps" "$IN" >> "$LOG" 2>&1
	
	echo "making pdf..."
	which ps2pdf > /dev/null 2>&1
	RC=$?
	CONVERT=""
	if [ "$RC" = 0 ] ; then
		CONVERT="ps2pdf"
	else
		which convert > /dev/null 2>&1
		RC=$?
		if [ "$RC" = 0 ] ; then
			CONVERT="convert"
		else
			echo "unable to find either ps2pdf or convert"
		fi
	fi
	if [ "$CONVERT" != "" ] ; then
		"$CONVERT" "$OUT.ps" "$OUT.pdf" >> "$LOG" 2>&1
	fi
fi

echo "making html..."
which latex2html > /dev/null 2>&1
RC=$?
if [ "$RC" != 0 ] ; then
	echo "unable to find latex2html"
else
	latex2html \
		-split 4 \
		-toc_depth 3 \
		-local_icons \
		-address "<a href=\"http://www.aero.polimi.it/~mbdyn\">MBDyn</a>: MultiBody Dynamics Software<br>Document version: $VERSION <br>Last update: $DATE <br>Maintained by <a href=\"mailto:mbdyn@aero.polimi.it\">mbdyn@aero.polimi.it</a>" \
		"$IN" >> "$LOG" 2>&1
	mv -f "$IN" "$OUT"
fi

