#!/bin/sh

# set as appropriate
NVCC=nvcc
SRCDIR=../..
BLDDIR=../..

#SRCDIR=/home/masarati/Lavoro/mbdyn/mbdyn-src
#BLDDIR=/home/masarati/Lavoro/mbdyn/mbdyn-devel

$NVCC \
-I $SRCDIR/include/ \
-I $SRCDIR/mbdyn/ \
-I $SRCDIR/libraries/libmbutil/ \
-I $SRCDIR/libraries/libmbmath/ \
-I $BLDDIR/include/ \
--shared \
-arch=sm_20 \
--compiler-options '-fPIC' \
-o $BLDDIR/modules/module-cudatest/cudatest.o \
-c cudatest.cu

cp -pf $BLDDIR/modules/module-cudatest/cudatest.o \
	$BLDDIR/modules/module-cudatest/.libs/cudatest.o

