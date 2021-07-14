#!/bin/sh
# $Header$

VERSIONFILE=../build/version

DATE=`date "+%a %b %e, %Y"`
if [ -r $VERSIONFILE ] ; then
	VERSION="`cat $VERSIONFILE`"
	VV="-$VERSION"
fi
LOG="make.log"

MAKEPS=no
MAKEPDF=yes
MAKEHTML=no
MAKESRC=no

TGT="input install tecman"

while test x$1 != x ; do
	case $1 in
	--ps=no) MAKEPS=no ; shift ;;
	--ps=yes|--ps) MAKEPS=yes ; shift ;;
	--pdf=no) MAKEPDF=no ; shift ;;
	--pdf=yes|--pdf) MAKEPDF=yes ; shift ;;
	--html=no) MAKEHTML=no ; shift ;;
	--html=yes|--html) MAKEHTML=yes ; shift ;;
	--src=no) MAKESRC=no ; shift ;;
	--src=yes|--src) MAKESRC=yes ; shift ;;
	--all) MAKEPDF=yes ; MAKEHTML=yes; MAKESRC=yes ; shift ;;
	--tgt) shift; TGT=$1; shift ;;
	*) echo "unknown arg $1" ; exit 1 ;;
	esac
done

rm -f "$LOG"

for IN in $TGT; do
	if ! test -d $IN ; then
		echo "invalid tgt=$IN"
		exit 1
	fi

	cd $IN

	OUT="mbdyn-${IN}$VV"

	rm -rf "$IN" "$OUT"

	if test $MAKEPDF = yes ; then
		echo "making $OUT.dvi..."
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

		echo "making $OUT.ps(.gz)..."
		which dvips > /dev/null 2>&1
		RC=$?
		if [ "$RC" != 0 ] ; then
			echo "unable to find dvips..."
		else
			dvips -o "$OUT.ps" "$IN" >> "$LOG" 2>&1
	
			echo "making $OUT.pdf..."
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
			if test "x$MAKEPS" = "xyes" ; then
				gzip "$OUT.ps"
			else
				rm -f "$OUT.ps"
			fi
		fi
	fi

	if test $MAKEHTML = yes -a $IN != "tecman" ; then
		echo "making $OUT.html..."
		which latex2html > /dev/null 2>&1
		RC=$?
		if [ "$RC" != 0 ] ; then
			echo "unable to find latex2html"
		else
			latex2html \
				-no_math -html_version 3.2,math \
				-split 4 \
				-toc_depth 3 \
				-local_icons \
				-address "<a href=\"http://www.aero.polimi.it/~mbdyn\">MBDyn</a>: MultiBody Dynamics Software<br>Document version: $VERSION <br>Last update: $DATE <br>Maintained by <a href=\"mailto:mbdyn@aero.polimi.it\">mbdyn@aero.polimi.it</a>" \
				"$IN" >> "$LOG" 2>&1
			mv -f "$IN" "$OUT"
		fi
	fi

	cd -
done

which doxygen > /dev/null 2>&1
RC1=$?
which dot > /dev/null 2>&1
RC2=$?
if [ "$RC1" != 0 ] ; then
	echo "unable to find doxygen"
elif [ "$RC2" != 0 ] ; then
	echo "unable to find dot"
else
	if test $MAKESRC = yes ; then
		echo "Generating html source documentation...."
		doxygen
	fi
fi
