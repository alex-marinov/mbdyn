#!/bin/sh
# $Header$

RCSITE=etc/rcsite
echo "Sourcing $RCSITE:"
cat $RCSITE
source $RCSITE
./configure 2>&1 | tee make.log
make 2>&1 | tee -a make.log

