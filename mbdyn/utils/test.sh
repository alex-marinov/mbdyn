#!/bin/sh

FILEEXT="mov jnt act"

# TIME=`awk -v STEP=$2 '/Step/ {j++; if (j==STEP) { print $5; exit }}' $1.out`
STEP=`awk -v TIME=$2 '/Step/ {j++; if ($5>=TIME) { print $2; exit }}' $1.out`
echo "Time "$2" corresponds to step "$STEP
for i in $FILEEXT ; do
    echo "    doing "$1"."$i" ..." ;
    awk -v K=$STEP 'BEGIN {getline;i=$1;j=0;} {if ($1==i) {j++;} if ( j==K) { print $0; } if (j==K+1) { exit;}}' $1.$i >$1.$i.log ;
done
echo "done"
