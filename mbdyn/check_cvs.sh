#!/bin/sh

cvs status -v 2>/dev/null | \
awk '
	/Repository revision/ {
		c = $3;
		f = $4;
	} 
	/MBDYN_REL_ENG/ {
		if (match($3, "([.0-9]+)")) {
			r = substr($3, RSTART, RLENGTH);
			if (c != r) {
				printf("%s: %s -> %s\n", f, c, r)
			}
		}
	}' | \
	sed "s;.*mbdyn-1\.0/;;" | \
	sed "s;,v;;"

