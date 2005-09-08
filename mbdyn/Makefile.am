include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = CHANGELOG \
contrib/SimulinkInterface/Makefile \
contrib/SimulinkInterface/MBDynCom.mdl \
contrib/SimulinkInterface/README \
contrib/SimulinkInterface/setup.m \
contrib/SimulinkInterface/sfun_mbdyn_com_read.c \
contrib/SimulinkInterface/sfun_mbdyn_com_write.c \
contrib/SimulinkInterface/sfun_mbdyn_start.c \
	contrib/SimulinkInterface/examples/pendulum \
	contrib/SimulinkInterface/examples/pendulum.mdl \
	contrib/SimulinkInterface/examples/README \
	contrib/SimulinkInterface/tmp/README \
contrib/Patches \
libraries/libboost

dist-hook:
	rm -rf `find $(distdir) -name CVS`



SUBDIRS = include build etc libraries manual mbdyn modules utils var

include $(top_srcdir)/build/bot.mk

