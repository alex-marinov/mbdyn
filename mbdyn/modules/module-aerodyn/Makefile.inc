MODULE_DEPENDENCIES=AeroMods.o AeroSubs.o GenSubs.o AeroDyn.o
# FIXME: needed to link objects compiled with Intel ifc
MODULE_LINK=-L /opt/intel/ia32/lib/ -lF90 -lIEPCF90 -lintrins
