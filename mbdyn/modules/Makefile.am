## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

.PHONY: all

all:
	for m in module-* ; do \
		cd $m ; \
		$(MAKE) ; \
	done

