include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

octavedir = "$(datarootdir)/octave"
c81dir = "$(datarootdir)/c81"
awkdir = "$(datarootdir)/awk"

export octavedir c81dir awkdir

codeaster_dir = "$(datarootdir)/Aster"
codeaster_srcdir = $(top_srcdir)/contrib/CodeAster
codeaster_files = \
$(codeaster_srcdir)/cms/cms.py

Pythonpre_srcdit = $(top_srcdir)/contrib/PythonPreprocessor
Pythonpre_files = \
$(Pythonpre_srcdit)/README \
$(Pythonpre_srcdit)/MBDynPreprocess.py \
$(Pythonpre_srcdit)/MBDynLib.py \
$(Pythonpre_srcdit)/examples/CrankPanel_v2.mbd \
$(Pythonpre_srcdit)/examples/include1.dat \
$(Pythonpre_srcdit)/examples/include2.dat \
$(Pythonpre_srcdit)/examples/LateralBuckling2_v2.mbd 

Post_srcdir = $(top_srcdir)/contrib/Post
Post_files = \
$(Post_srcdir)/Makefile \
$(Post_srcdir)/Mbdynpost.flex \
$(Post_srcdir)/Post.C \
$(Post_srcdir)/Post.hh

mls_dir = $(octavedir)
mls_srcdir = $(top_srcdir)/contrib/MLS
mls_files = \
$(mls_srcdir)/COPYING.kd \
$(mls_srcdir)/create_mls_interface.m \
$(mls_srcdir)/example_MLS.m \
$(mls_srcdir)/f_polyn_basis_derivatives.m \
$(mls_srcdir)/f_polyn_basis.m \
$(mls_srcdir)/f_regression_derivatives_kdtree.m \
$(mls_srcdir)/f_weight_function_derivatives.m \
$(mls_srcdir)/f_weight_function_givendist.m \
$(mls_srcdir)/kd_buildtree.m \
$(mls_srcdir)/kd_closestpointfast.m \
$(mls_srcdir)/kd_closestpointgood.m \
$(mls_srcdir)/kd_demo.m \
$(mls_srcdir)/kd_knn.m \
$(mls_srcdir)/kd_nclosestpoints.m \
$(mls_srcdir)/kd_plotbox.m \
$(mls_srcdir)/kd_rangequery.m \
$(mls_srcdir)/README

scicos_dir = "$(datarootdir)/Scicos"
scicos_srcdir = $(top_srcdir)/contrib/ScicosInterface
scicos_files = \
$(scicos_srcdir)/buildblocks.sci \
$(scicos_srcdir)/loadblocks.sci \
$(scicos_srcdir)/sockread_intf.sci \
$(scicos_srcdir)/sockwrite_intf.sci \
$(scicos_srcdir)/sockread_comp.c \
$(scicos_srcdir)/sockwrite_comp.c \
$(scicos_srcdir)/README \
$(scicos_srcdir)/TODO \
$(scicos_srcdir)/examples/pendulum \
$(scicos_srcdir)/examples/Pendulum.cosf

simulink_dir = "$(datarootdir)/Simulink"
simulink_srcdir = $(top_srcdir)/contrib/SimulinkInterface
simulink_files = \
$(simulink_srcdir)/Makefile \
$(simulink_srcdir)/MBDynCom.mdl \
$(simulink_srcdir)/mbdyn.sh \
$(simulink_srcdir)/README \
$(simulink_srcdir)/setup.m \
$(simulink_srcdir)/sfun_mbdyn_com_read.c \
$(simulink_srcdir)/sfun_mbdyn_com_write.c \
$(simulink_srcdir)/sfun_mbdyn_start.c \
	$(simulink_srcdir)/examples/pendulum \
	$(simulink_srcdir)/examples/pendulum.mdl \
	$(simulink_srcdir)/examples/README \
	$(simulink_srcdir)/tmp/README

EXTRA_DIST = CHANGES \
$(top_srcdir)/contrib/Patches \
$(codeaster_files) \
$(Pythonpre_files) \
$(Post_files) \
$(mls_files) \
$(scicos_files) \
$(simulink_files)

install-data-local-codeaster: $(codeaster_files)
	$(mkinstalldirs) $(DESTDIR)$(datarootdir) $(DESTDIR)$(codeaster_dir)
	@for f in $^ ; do \
		case $$f in \
		*.py) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(codeaster_dir)" ; $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(codeaster_dir) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

install-data-local: install-data-local-codeaster

#dist-hook:
#	rm -rf `find $(distdir) -name CVS`

uninstall-local:
	rm -f $(DESTDIR)$(codeaster_dir)/cms.py

SUBDIRS = include build etc libraries manual mbdyn modules utils var

include $(top_srcdir)/build/bot.mk

