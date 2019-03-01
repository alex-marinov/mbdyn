#!/bin/sh
# $Header$

for i in `find . -name Makefile`; do 
	sed "s;^subdir = \([^.]\);subdir = ./\1;" $i > x
	mv -f x $i
done
