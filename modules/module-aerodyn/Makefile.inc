## $Header$
MODULE_DEPENDENCIES=mbdyn_ad.lo
MODULE_LINK=-L. -lAeroDyn
# FIXME: needed to link objects compiled with Intel ifc
#MODULE_LINK=-L/opt/intel/ia32/lib/ -lF90 -lIEPCF90 -lintrins $(MODULE_LINK)
