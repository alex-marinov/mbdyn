## $Header$

### set FG version
FGVERSION = "2.4.0"
FGSRCPATH = /home/masarati/Lavoro/mbdyn/flightgear

FGVER = `echo $(FGVERSION) | sed -e 's/\./0/g'`

MODULE_INCLUDE = -D FGVERSION="$(FGVERSION)" -D FGVER=$(FGVER) -I $(FGSRCPATH)/simgear-$(FGVERSION)/ -I $(FGSRCPATH)/flightgear-$(FGVERSION)/src/

