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

ASTERCMD="/home/masarati/Lavoro/mbdyn/ASTER/aster-9.3.0/outils/as_run"
WORKDIR="/home/masarati/Lavoro/mbdyn/ASTER/2008-07-17-RomanelliSerioli/Code_Aster/TUTORIAL/cms"
EXPORTINFNAME="cms.export.in"
EXPORTFNAME="${FNAME}.export"
COMMFNAME="${FNAME}.comm"

cat "$EXPORTINFNAME" | sed \
	-e "s;@WORKDIR@;$WORKDIR;g" \
	-e "s;@COMMFNAME@;$COMMFNAME;g" \
	-e "s;@MAILFNAME@;$MAILFNAME;g" \
	> "$EXPORTFNAME"

$ASTERCMD $EXPORTFNAME

rm -f "$EXPORTFNAME"

