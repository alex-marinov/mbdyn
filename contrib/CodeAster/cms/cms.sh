#!/bin/sh

if test x"$1" = x ; then
	echo "need file number"
	exit 1
fi
COMMNUMBER="$1"
FNAME=`printf "cms%05d" $COMMNUMBER`

if test x"$2" != x ; then
	MAILNUMBER="$2"
	MAILFNAME=`printf "cms%05d.mail" $MAILNUMBER`
elif test -f "${FNAME}.mail" ; then
	MAILFNAME="${FNAME}.mail"
else
	MAILFNAME="cms.mail"
fi

set +xv

### modify as appropriate
ASTERVERSION="9.3"
ASTERBINPATH="/home/masarati/Lavoro/mbdyn/ASTER/aster-9.3.0-2/outils"
#ASTERVERSION="10.2"
#ASTERBINPATH="/opt/aster/bin"
MBDYNSRCDIR="/home/masarati/Lavoro/mbdyn/mbdyn-src"

ASTERCMD="$ASTERBINPATH/as_run"
WORKDIR="$MBDYNSRCDIR/contrib/CodeAster/cms"
EXPORTINFNAME="cms.export.in"
EXPORTFNAME="${FNAME}.export"
COMMFNAME="${FNAME}.comm"

cat "$EXPORTINFNAME" | sed \
	-e "s;@ASTERVERSION@;$ASTERVERSION;g" \
	-e "s;@WORKDIR@;$WORKDIR;g" \
	-e "s;@COMMFNAME@;$COMMFNAME;g" \
	-e "s;@MAILFNAME@;$MAILFNAME;g" \
	> "$EXPORTFNAME"

$ASTERCMD $EXPORTFNAME

rm -f "$EXPORTFNAME"

