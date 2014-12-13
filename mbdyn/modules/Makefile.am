include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

export LTCXXCOMPILE = @LIBTOOL@ --tag=CXX --mode=compile @CXX@ @DEFS@ \
	-I. -I@srcdir@ \
	-I../../include \
	-I../$(top_srcdir)/include \
	-I../$(top_srcdir)/libraries/libmbutil \
	-I../$(top_srcdir)/libraries/libmbmath \
	-I../$(top_srcdir)/libraries/libmbwrap \
	-I../$(top_srcdir)/libraries/libobjs \
	-I../$(top_srcdir)/mbdyn \
	-I../$(top_srcdir)/mbdyn/base \
	-I../$(top_srcdir)/mbdyn/struct \
	-I../$(top_srcdir)/mbdyn/aero \
	-I../$(top_srcdir)/mbdyn/elec \
	-I../$(top_srcdir)/mbdyn/hydr \
	-prefer-pic \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CXXFLAGS) @CXXFLAGS@
CXXLD = @CXX@
export CXXLINK = @LIBTOOL@ --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	@CXXFLAGS@ @LDFLAGS@ -shared @FCLIBS@
export LTCOMPILE = @LIBTOOL@ --mode=compile @CC@ @DEFS@ \
	-I. -I@srcdir@ \
	-I../../include \
	-I../$(top_srcdir)/include \
	-I../$(top_srcdir)/libraries/libmbutil \
	-I../$(top_srcdir)/libraries/libmbmath \
	-I../$(top_srcdir)/libraries/libmbwrap \
	-I../$(top_srcdir)/libraries/libobjs \
	-I../$(top_srcdir)/mbdyn \
	-I../$(top_srcdir)/mbdyn/base \
	-I../$(top_srcdir)/mbdyn/struct \
	-I../$(top_srcdir)/mbdyn/aero \
	-I../$(top_srcdir)/mbdyn/elec \
	-I../$(top_srcdir)/mbdyn/hydr \
	-prefer-pic \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CFLAGS) @CFLAGS@
export LTF77COMPILE = @LIBTOOL@ --tag=F77 --mode=compile @F77@ $(AM_FFLAGS) @FFLAGS@ \
	-prefer-pic
export LTFCCOMPILE = @LIBTOOL@ --tag=FC --mode=compile @FC@ $(AM_FFLAGS) @FFLAGS@ \
	-prefer-pic
export LIBTOOL = @LIBTOOL@
export MODULES_DIRS=@MODULES_DIRS@
export INSTALL_PROGRAM=${INSTALL}
export LIBEXEC_DIR=$(libexecdir)

.PHONY: all all__ clean distclean

# Note: this hack is required because macro "LIBTOOL" may expand
# to sh ../libtool, so, when we change directory to the modules,
# the effective libtool is not found.

all: libtool all__

libtool: 
	-@ln -sf ../libtool .
	
all__:
	@echo "Configured modules: $(MODULES_DIRS)"
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Playing with $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(srcdir)/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f,f90,f95} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f|\*.f90|\*.f95) ;; *) if test ! -f $$i/$$f ; then ( cd $$i && ln -sf ../$$isrc/$$f $$f ) ; fi ; ;; esac ; \
			done ; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			( cd $$i && $(MAKE) MBDYN_MODULE_SRCDIR=../$$isrc ) ; \
			rm -f $$i/Makefile ; \
		done \
	fi

install: all__
	$(mkinstalldirs) $(libexecdir)
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Installing $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(srcdir)/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f,f90,f95} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f|\*.f90|\*.f95) ;; *) if test ! -f $$i/$$f ; then ( cd $$i && ln -sf ../$$isrc/$$f $$f ) ; fi ; ;; esac ; \
			done ; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			( cd $$i && $(MAKE) install ) ; \
			rm -f $$i/Makefile ; \
		done \
	fi

clean clean-so:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Cleaning $$i ..."; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) $@ ; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi

distclean:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Cleaning $$i ..."; \
			rm -rf $$i/*.so $$i/*.lo $$i/*.o ; \
		done \
	fi
	rm -f Makefile libtool

# Add distributed modules below
EXTRA_DIST = Makefile.module Makefile.template

EXTRA_DIST += \
module-aerodyn/NREL_AeroDyn.h \
module-aerodyn/Makefile.inc \
module-aerodyn/mbdyn_ad.f90 \
module-aerodyn/module-aerodyn.cc \
module-aerodyn/README

EXTRA_DIST += \
module-asynchronous_machine/asynchronous_machine.mbdyn \
module-asynchronous_machine/module-asynchronous_machine.cc \
module-asynchronous_machine/module-asynchronous_machine.h \
module-asynchronous_machine/doc/Doxyfile \
module-asynchronous_machine/doc/Makefile

EXTRA_DIST += \
module-autodiff_test/autodiff_test.f \
module-autodiff_test/Makefile.inc \
module-autodiff_test/module-autodiff_test.cc \
module-autodiff_test/module-autodiff_test.h \
module-autodiff_test/examples/deformable_joint_ad.cod \
module-autodiff_test/examples/deformable_joint_ad.elm \
module-autodiff_test/examples/deformable_joint_ad.inv \
module-autodiff_test/examples/deformable_joint_ad.nod \
module-autodiff_test/examples/deformable_joint_ad.ref \
module-autodiff_test/examples/deformable_joint_ad.set \
module-autodiff_test/examples/deformable_joint_ad_f77.mbdyn \
module-autodiff_test/examples/deformable_joint_ad_tmp.mbdyn \
module-autodiff_test/examples/deformable_joint_ad_trad.mbdyn \
module-autodiff_test/examples/inline_joint_ad.mbdyn \
module-autodiff_test/examples/inline_joint_ad_2.mbdyn \
module-autodiff_test/examples/inline_joint_ad_2_elements.elm \
module-autodiff_test/examples/inline_joint_ad_2_parameters.set \
module-autodiff_test/examples/inline_joint_ad_2_tmp.elm \
module-autodiff_test/examples/inline_joint_ad_2_trad.elm

EXTRA_DIST += \
module-constlaw/module-constlaw.cc \
module-constlaw/README

EXTRA_DIST += \
module-constlaw-f90/Makefile.inc \
module-constlaw-f90/module-constlaw-f90.cc \
module-constlaw-f90/README \
module-constlaw-f90/usrsub.f90 \
module-constlaw-f90/usrsub.h

EXTRA_DIST += \
module-constlaw-f95/Makefile.inc \
module-constlaw-f95/module-constlaw-f95.cc \
module-constlaw-f95/README \
module-constlaw-f95/usrsub.f95 \
module-constlaw-f95/usrsub.h

EXTRA_DIST += \
module-cont-contact/module-cont-contact.cc

EXTRA_DIST += \
module-cyclocopter/module-cyclocopter.cc \
module-cyclocopter/module-cyclocopter.h

EXTRA_DIST += \
module-damper-graall/damper_graall.f \
module-damper-graall/Makefile.inc \
module-damper-graall/module-damper-graall.cc

EXTRA_DIST += \
module-damper-hydraulic/module-damper-hydraulic.cc

EXTRA_DIST += \
module-drive/module-drive.cc \
module-drive/README

EXTRA_DIST += \
module-eu2phi/module-eu2phi.cc

EXTRA_DIST += \
module-hid/README \
module-hid/Makefile.inc \
module-hid/module-hid.cc \
module-hid/hid_detect.c \
module-hid/hid.set \
module-hid/hid_ms_sidewinder.set \
module-hid/hid_thrustmaster.set \
module-hid/hid_xxx.set \
module-hid/hid_aircraft \
module-hid/hid_test

EXTRA_DIST += \
module-hydrodynamic_plain_bearing/diffsizes.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force_dv.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force.h \
module-hydrodynamic_plain_bearing/Makefile.inc \
module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.cc \
module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.h \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing.mbdyn

EXTRA_DIST += \
module-hunt-crossley/module-hunt-crossley.cc \
module-hunt-crossley/hunt_crossley

EXTRA_DIST += \
module-imu/module-imu.cc \
module-imu/imu \
module-imu/imu_constraint \
module-imu/imu.tex \
module-imu/README

EXTRA_DIST += \
module-inline_friction/module-inline_friction.cc \
module-inline_friction/module-inline_friction.h \
module-inline_friction/@InLineFriction/AssJac.m \
module-inline_friction/@InLineFriction/AssRes.m \
module-inline_friction/@InLineFriction/ComputeDeltaXP.m \
module-inline_friction/@InLineFriction/ComputeInitialResidual.m \
module-inline_friction/@InLineFriction/ComputeResidual.m \
module-inline_friction/@InLineFriction/ComputeTau.m \
module-inline_friction/@InLineFriction/DescribeDof.m \
module-inline_friction/@InLineFriction/DescribeEq.m \
module-inline_friction/@InLineFriction/dGetPrivData.m \
module-inline_friction/@InLineFriction/display.m \
module-inline_friction/@InLineFriction/GetConnectedNodes.m \
module-inline_friction/@InLineFriction/GetDofType.m \
module-inline_friction/@InLineFriction/GetEqType.m \
module-inline_friction/@InLineFriction/GetInitialStateVector.m \
module-inline_friction/@InLineFriction/GetStateVector.m \
module-inline_friction/@InLineFriction/iGetInitialNumDof.m \
module-inline_friction/@InLineFriction/iGetNumConnectedNodes.m \
module-inline_friction/@InLineFriction/iGetNumDof.m \
module-inline_friction/@InLineFriction/iGetNumPrivData.m \
module-inline_friction/@InLineFriction/iGetPrivDataIdx.m \
module-inline_friction/@InLineFriction/InitialAssJac.m \
module-inline_friction/@InLineFriction/InitialAssRes.m \
module-inline_friction/@InLineFriction/InitialWorkSpaceDim.m \
module-inline_friction/@InLineFriction/InLineFriction.m \
module-inline_friction/@InLineFriction/SetInitialValue.m \
module-inline_friction/@InLineFriction/SetValue.m \
module-inline_friction/@InLineFriction/Update.m \
module-inline_friction/@InLineFriction/WorkSpaceDim.m \
module-inline_friction/examples/inline_friction2_bench.elm \
module-inline_friction/examples/inline_friction2.elm \
module-inline_friction/examples/inline_friction2.mbdyn \
module-inline_friction/examples/inline_friction2.nod \
module-inline_friction/examples/inline_friction.mbdyn

EXTRA_DIST += \
module-inplane_friction/@InPlaneFriction/AssJac.m \
module-inplane_friction/@InPlaneFriction/AssRes.m \
module-inplane_friction/@InPlaneFriction/ComputeInitialResidual.m \
module-inplane_friction/@InPlaneFriction/ComputeResidual.m \
module-inplane_friction/@InPlaneFriction/DescribeDof.m \
module-inplane_friction/@InPlaneFriction/DescribeEq.m \
module-inplane_friction/@InPlaneFriction/dGetPrivData.m \
module-inplane_friction/@InPlaneFriction/display.m \
module-inplane_friction/@InPlaneFriction/GetConnectedNodes.m \
module-inplane_friction/@InPlaneFriction/GetDofType.m \
module-inplane_friction/@InPlaneFriction/GetEqType.m \
module-inplane_friction/@InPlaneFriction/GetInitialStateVector.m \
module-inplane_friction/@InPlaneFriction/GetStateVector.m \
module-inplane_friction/@InPlaneFriction/iGetInitialNumDof.m \
module-inplane_friction/@InPlaneFriction/iGetNumConnectedNodes.m \
module-inplane_friction/@InPlaneFriction/iGetNumDof.m \
module-inplane_friction/@InPlaneFriction/iGetNumPrivData.m \
module-inplane_friction/@InPlaneFriction/iGetPrivDataIdx.m \
module-inplane_friction/@InPlaneFriction/InitialAssJac.m \
module-inplane_friction/@InPlaneFriction/InitialAssRes.m \
module-inplane_friction/@InPlaneFriction/InitialWorkSpaceDim.m \
module-inplane_friction/@InPlaneFriction/InPlaneFriction.m \
module-inplane_friction/@InPlaneFriction/SetInitialValue.m \
module-inplane_friction/@InPlaneFriction/SetValue.m \
module-inplane_friction/@InPlaneFriction/Update.m \
module-inplane_friction/@InPlaneFriction/WorkSpaceDim.m \
module-inplane_friction/examples/inplane_friction.mbdyn

EXTRA_DIST += \
module-mds/mds \
module-mds/module-mds.cc

EXTRA_DIST += \
module-minmaxdrive/minmax.mbdyn \
module-minmaxdrive/module-minmaxdrive.cc \
module-minmaxdrive/module-minmaxdrive.h

EXTRA_DIST += \
module-muscles/muscle_test \
module-muscles/module-muscles.cc

EXTRA_DIST += \
module-namespace/module-namespace.cc \
module-namespace/tns

EXTRA_DIST += \
module-nodedistdrive/module-nodedistdrive.cc \
module-nodedistdrive/module-nodedistdrive.h

EXTRA_DIST += \
module-nonsmooth-node/Makefile.inc \
module-nonsmooth-node/mbdyn_siconos.cc \
module-nonsmooth-node/mbdyn_siconos.h \
module-nonsmooth-node/module-nonsmooth-node.cc

EXTRA_DIST += \
module-octave/Makefile.inc \
module-octave/module-octave.cc \
module-octave/module-octave.h \
module-octave/octave_object.cc \
module-octave/octave_object.h \
module-octave/README \
module-octave/demo/my_func2.m \
module-octave/demo/my_func2.m \
module-octave/demo/my_func3.m \
module-octave/demo/my_func4.m \
module-octave/demo/my_func5.m \
module-octave/demo/my_func.m \
module-octave/demo/my_tpl_drive1D.m \
module-octave/demo/my_tpl_drive3x3D.m \
module-octave/demo/octave_const_law.mbdyn \
module-octave/demo/octave_drive_caller.mbdyn \
module-octave/demo/octave_element.mbdyn \
module-octave/demo/@MyConstLaw1/GetConstLawType.m \
module-octave/demo/@MyConstLaw1/MyConstLaw1.m \
module-octave/demo/@MyConstLaw1/Update.m \
module-octave/demo/@MyElem1/AfterConvergence.m \
module-octave/demo/@MyElem1/AssJac.m \
module-octave/demo/@MyElem1/AssRes.m \
module-octave/demo/@MyElem1/DescribeDof.m \
module-octave/demo/@MyElem1/DescribeEq.m \
module-octave/demo/@MyElem1/dGetPrivData.m \
module-octave/demo/@MyElem1/display.m \
module-octave/demo/@MyElem1/GetConnectedNodes.m \
module-octave/demo/@MyElem1/GetDofType.m \
module-octave/demo/@MyElem1/GetEqType.m \
module-octave/demo/@MyElem1/iGetInitialNumDof.m \
module-octave/demo/@MyElem1/iGetNumConnectedNodes.m \
module-octave/demo/@MyElem1/iGetNumDof.m \
module-octave/demo/@MyElem1/iGetNumPrivData.m \
module-octave/demo/@MyElem1/iGetPrivDataIdx.m \
module-octave/demo/@MyElem1/InitialAssJac.m \
module-octave/demo/@MyElem1/InitialAssRes.m \
module-octave/demo/@MyElem1/InitialWorkSpaceDim.m \
module-octave/demo/@MyElem1/MyElem1.m \
module-octave/demo/@MyElem1/Output.m \
module-octave/demo/@MyElem1/SetInitialValue.m \
module-octave/demo/@MyElem1/SetValue.m \
module-octave/demo/@MyElem1/Update.m \
module-octave/demo/@MyElem1/WorkSpaceDim.m

EXTRA_DIST += \
module-scalarfunc/module-scalarfunc.cc \
module-scalarfunc/README

EXTRA_DIST += \
module-tclpgin/Makefile.inc \
module-tclpgin/module-tclpgin.cc

EXTRA_DIST += \
module-template/module-template.cc

EXTRA_DIST += \
module-template2/module-template2.cc

EXTRA_DIST += \
module-udunits/Makefile.inc \
module-udunits/module-udunits.cc \
module-udunits/README

EXTRA_DIST += \
module-wheel2/module-wheel2.cc \
module-wheel2/module-wheel2.h

EXTRA_DIST += \
module-wheel4/module-wheel4.cc \
module-wheel4/module-wheel4.h \
module-wheel4/axleExampleNoData

include $(top_srcdir)/build/bot.mk
