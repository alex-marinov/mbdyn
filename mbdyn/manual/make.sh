#!/bin/sh

#	-address 'Maintained by <a href="mailto:mbdyn@aero.polimi.it">mbdyn@aero.polimi.it</a>'

latex2html \
	-split 4 \
	-toc_depth 3 \
	-local_icons \
	input

