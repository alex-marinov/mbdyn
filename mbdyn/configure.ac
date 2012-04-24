dnl $Header$
dnl Process this file with autoconf to produce a configure script.
dnl
dnl MBDyn (C) is a multibody analysis code.
dnl http://www.mbdyn.org
dnl
dnl Copyright (C) 1996-2012
dnl
dnl Pierangelo Masarati     <masarati@aero.polimi.it>
dnl Paolo Mantegazza        <mantegazza@aero.polimi.it>
dnl
dnl Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
dnl via La Masa, 34 - 20156 Milano, Italy
dnl http://www.aero.polimi.it
dnl
dnl Changing this copyright notice is forbidden.
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation (version 2 of the License).
dnl 
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
dnl
dnl ----------------------------------------------------------------
dnl Redefine AC_INIT_BINSH to provide copyright notice at top 
dnl of generated configure script.  Prints simple copyright.
define([AC_INIT_BINSH],
[[#! /bin/sh
] ] [

#
# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2012
#
# Pierangelo Masarati     <masarati@aero.polimi.it>
# Paolo Mantegazza        <mantegazza@aero.polimi.it>
#
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
#
# Changing this copyright notice is forbidden.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
#
echo ""
echo "	Copyright (C) 1996-2012"
echo "	"
echo "	Pierangelo Masarati     <masarati@aero.polimi.it>"
echo "	Paolo Mantegazza        <mantegazza@aero.polimi.it>"
echo "	"
echo "	Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano"
echo "	via La Masa, 34 - 20156 Milano, Italy"
echo "	http://www.aero.polimi.it"
echo "	"
echo "	Restrictions apply, see COPYRIGHT and LICENSE files."
echo ""
])dnl
dnl ----------------------------------------------------------------
dnl Disable config.cache!
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl

dnl ================================================================
dnl configure.in for MBDyn
AC_INIT(build/version)

DATE=`date "+%B %e, %Y"`

# set unset (borrowed from autoconf 2.14a)
if (unset FOO) >/dev/null 2>&1; then
	ol_unset=unset
else
	ol_unset=false
fi
# unset CDPATH
$ol_unset CDPATH || test "${CDPATH:+set}" != set || CDPATH=: && export CDPATH

AC_CONFIG_AUX_DIR(build)dnl

MB_VERSION=`cat $ac_aux_dir/version`
if test -z "$MB_VERSION"; then
	AC_MSG_ERROR([could not determine version])
fi
echo "Configuring MBDyn $MB_VERSION ..."

dnl Determine host platform
dnl             we try not to use this for much
AC_CANONICAL_SYSTEM

AM_INIT_AUTOMAKE(mbdyn, [$MB_VERSION])

dnl We use autoconf features new to 2.13.
dnl aclocal.m4 should be built using aclocal from automake 1.4
dnl libtool 1.3.3 should be installed.
AC_PREREQ(2.13)dnl Required Autoconf version

AC_CONFIG_HEADERS([include/mbconfig.h:include/mbconfig.h.in])
AH_TOP([
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */ 
])
AH_BOTTOM([
/* this file must be included by every .c/.cc file */
#include "mbdyn.h"
])

top_builddir=`pwd`
dnl AC_SUBST(top_builddir)dnl

dnl ----------------------------------------------------------------
dnl Checks for args.
dnl
AC_MSG_CHECKING(configure arguments)
AC_PREFIX_DEFAULT(/usr/local/mbdyn)

dnl ----------------------------------------------------------------
dnl --with-subdir
mbdyn_subdir="/mbdyn"

AC_ARG_WITH(subdir,
[  --with-subdir=DIR       change default subdirectory used for installs],
[case "$withval" in
	no) 
		mbdyn_subdir=""
		;;
	yes)
		;;
	/*|\\*)
		mbdyn_subdir="$withval"
		;;
	*)
		mbdyn_subdir="/$withval"
		;;
esac
])dnl
AC_SUBST(mbdyn_subdir)dnl

dnl ----------------------------------------------------------------
dnl General "enable" options
OL_ARG_ENABLE(debug,[  --enable-debug          enable debugging], no)dnl
OL_ARG_WITH(debug_mode,[  --with-debug-mode[={none|mem}] with debug mode {none|mem}],none,[none mem])
OL_ARG_ENABLE(runtime_loading,[  --enable-runtime-loading    enable runtime loading], auto)dnl
OL_ARG_WITH(static_modules,[  --with-static-modules   build (known) modules as static],auto,[no yes auto])
OL_ARG_ENABLE(crypt,[  --enable-crypt          enable crypt (deprecated)], no)dnl
OL_ARG_ENABLE(eig,[  --enable-eig            enable direct eigensolution
                          (needs either LAPACK, ARPACK or JDQZ)],auto,[auto yes no])dnl needs LAPACK, ARPACK or JDQZ
OL_ARG_ENABLE(schur,[  --enable-schur          enable Schur parallel solver
                          (needs MPI and either Metis or Chaco)],auto,[auto yes no])dnl needs MPI and Metis or Chaco
OL_ARG_ENABLE(multithread,[  --enable-multithread    enable multithread assembly],no,[auto yes no force])dnl
OL_ARG_ENABLE(multithread_naive,[  --enable-multithread-naive    enable multithread naive solver],no,[auto yes no force])dnl
OL_ARG_ENABLE(adams,[  --enable-adams          enable MSC.ADAMS output],no)dnl
OL_ARG_ENABLE(motionview,[  --enable-motionview     enable Altair's Motion View output],no)dnl
OL_ARG_ENABLE(netcdf,[  --enable-netcdf         enable NetCDF based binary output],auto)dnl
OL_ARG_ENABLE(python,[  --enable-python         enable Python support],no)dnl
OL_ARG_ENABLE(octave,[  --enable-octave         enable Octave support],no)dnl

dnl ----------------------------------------------------------------
dnl General "with" options
OL_ARG_WITH(g2c,[  --with-g2c              with g2c header],no,[auto yes no])
OL_ARG_WITH(fs,[  --with-fs[={unix|dos}]    filesystem type],unix,[unix dos])
OL_ARG_WITH(mpi,[  --with-mpi              with MPI support (=pmpi for profiling)],auto,[auto yes pmpi no])
OL_ARG_ENABLE(debug_mpi,[  --enable-debug-mpi      enable MPI debugging], no)dnl
OL_ARG_WITH(metis,[  --with-metis            with Metis model partitioning support],auto,[auto yes no])

dnl temporarily disabled until we figure out how to make chaco work...
OL_ARG_WITH(chaco,[  --with-chaco            with Chaco model partitioning support],no,[auto yes no])
dnl ol_with_chaco=no

OL_ARG_WITH(threads,[  --with-threads          with threads],auto,[auto yes no])
OL_ARG_WITH(rt,[  --with-rt               with POSIX realtime support],auto,[auto yes no])
OL_ARG_WITH(rtai,[  --with-rtai             with RTAI support],no,[auto yes no])
OL_ARG_WITH(boost,[  --with-boost            with BOOST support],no,[auto yes no])
OL_ARG_WITH(ann,[  --with-ann              with ANN support],no,[auto yes no])
OL_ARG_WITH(bullet,[  --with-bullet           with bullet support],no,[auto yes no])
OL_ARG_WITH(blas,[
       math libraries:
  --with-blas[={auto|goto|atlas|cblas|blas}]
                          with (C)BLAS math library],auto,[auto goto atlas cblas blas])
AC_ARG_WITH(goto,[  --with-goto=lib(s)      with Goto BLAS implementation],[
	ol_with_goto_lib="$withval"
],[ol_with_goto_lib=])
OL_ARG_WITH(ginac,[  --with-ginac            with GiNaC support
                          (ginac-config must be in $PATH)],auto,[auto yes no])
OL_ARG_WITH(y12,[
       optional linear algebra solvers (naive is enabled by default):
  --with-y12              with Y12 sparse math library],yes,[auto yes no])
OL_ARG_WITH(umfpack,[  --with-umfpack          with Umfpack math library],auto,[auto yes no])
OL_ARG_WITH(klu,[  --with-klu          with KLU math library],auto,[auto yes no])
OL_ARG_WITH(lapack,[  --with-lapack           with LAPACK math library],auto,[auto yes no])
OL_ARG_WITH(harwell,[  --with-hsl              with HSL (Harwell) sparse math library - historical],no,[auto yes no])
OL_ARG_WITH(meschach,[  --with-meschach         with Meschach math library - historical],no,[auto yes no])
OL_ARG_WITH(superlu,[  --with-superlu          with SuperLU math library - eXperimental],no,[auto yes no])
OL_ARG_WITH(wsmp,[  --with-wsmp             with Watson library - eXperimental],no,[auto yes no])
dnl temporarily disabled until we find out how to make it work...
dnl OL_ARG_WITH(taucs,[  --with-taucs            with TAUCS math library - eXperimental],auto,[auto yes no])
ol_with_taucs=no
OL_ARG_WITH(arpack,[  --with-arpack           with ARPACK math library],auto,[auto yes no])
OL_ARG_WITH(jdqz,[  --with-jdqz             with JDQZ math library],auto,[auto yes no])

OL_ARG_WITH(pam,[
       misc security libraries:
  --with-pam              with PAM support],no,[auto yes no])
OL_ARG_WITH(sasl2,[  --with-sasl2            with Cyrus SASL2 support],no,[auto yes no])

dnl module with's
AC_ARG_WITH(module,[
  --with-module=<list>    build listed modules (see modules/)
],[
	ol_arg=
	for ol_val in $withval ; do
		ol_arg="$ol_arg module-$ol_val"
	done
	if test "${ol_arg:-set}" = set ; then
		AC_MSG_ERROR(bad value $withval for --with-module)
	fi
	MODULES_DIRS=$ol_arg
],[MODULES_DIRS=])dnl

OL_ARG_WITH(charm,[  --with-charm            with CHARM WP module],no,[auto yes no])

fool_libtool=yes
case "$host" in
	*-apple-*)
		fool_libtool=no
		;;
	*)
		;;
esac

dnl ----------------------------------------------------------------
dnl
dnl Checks for programs.
dnl
AC_PROG_CC
if test "x$CC" = "x" ; then
	AC_MSG_ERROR([Need a working C compiler])
fi
AC_PROG_CXX
if test "x$CXX" = "x" ; then
	AC_MSG_ERROR([Need a working C++ compiler])
fi
AC_PROG_FC
AC_FC_DUMMY_MAIN
AC_FC_LIBRARY_LDFLAGS
if test "x$FC" = "x" ; then
	AC_PROG_F77
	if test "x$F77" = "x" ; then
		AC_MSG_ERROR([Need a working Fortran compiler])
	fi
	AC_F77_LIBRARY_LDFLAGS
	LIBS="$LIBS $F77LIBS"
else
	F77="$FC"
	F77LIBS="$FCLIBS"
	LIBS="$LIBS $FCLIBS"
fi

if test "$GCC" = "yes" ; then
	dnl check for 3.3
	AC_MSG_CHECKING(for gcc version = 3.3)
	AC_CACHE_VAL(ac_cv_gcc33,
	AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	[AC_TRY_COMPILE(
	[#if !(__GNUC__ == 3 && __GNUC_MINOR__ == 3)
	#error GCC is not 3.3 (good news ;)
	#endif],
	[return 0;],
	[ac_cv_gcc33=yes],
	[ac_cv_gcc33=no])]
	AC_LANG_RESTORE)
	if [[ $ac_cv_gcc33 = no ]]; then
		AC_MSG_RESULT([no])

		dnl then check for 3.2
		AC_MSG_CHECKING(for gcc version = 3.2)
		AC_CACHE_VAL(ac_cv_gcc32,
		AC_LANG_SAVE
		AC_LANG_CPLUSPLUS
		[AC_TRY_COMPILE(
		[#if !(__GNUC__ == 3 && __GNUC_MINOR__ == 2)
		#error GCC is not 3.2 (good news ;)
		#endif],
		[return 0;],
		[ac_cv_gcc32=yes],
		[ac_cv_gcc32=no])]
		AC_LANG_RESTORE)
		if [[ $ac_cv_gcc32 = no ]]; then
			AC_MSG_RESULT([no])
		else
			AC_DEFINE(MBDYN_X_WORKAROUND_GCC_3_2,1,[define to workaround known GCC 3.2 bug])
			AC_MSG_RESULT([ok])
		fi
	else
		AC_DEFINE(MBDYN_X_WORKAROUND_GCC_3_3,1,[define to workaround known GCC 3.3 bug])
		AC_MSG_RESULT([ok])
	fi
fi

dnl AM_DISABLE_SHARED dnl Remove when you want to build shared libs also
AM_PROG_LIBTOOL

dnl ----------------------------------------------------------------
dnl
dnl Initialize vars
dnl
MPI_LIBS=
METIS_LIBS=
CHACO_LIBS=
Y12_LIBS=
UMFPACK_LIBS=
KLU_LIBS=
BLAS_LIBS=
THREAD_LIBS=
ATOMIC_OPS_LIBS=
MESCHACH_LIBS=
SUPERLU_LIBS=
WSMP_LIBS=
HARWELL_LIBS=
LAPACK_LIBS=
TAUCS_LIBS=
ARPACK_LIBS=
JDQZ_LIBS=
SECURITY_LIBS=
NETCDF_LIBS=
UDUNITS_LIBS=
RT_LIBS=
ANN_LIBS=
BULLET_LIBS=
PYTHON_NUMPY_HEADERS=
OCTAVE_LIBS=

dnl GINACLIB_CONFIG=
GINACLIB_LIBS=
GINACLIB_CPPFLAGS=

dnl MODULES_DIRS=
MODULE_LOADPATH=

dnl ----------------------------------------------------------------
dnl
dnl MBDyn requires STDC features
dnl
AM_PROG_CC_STDC
if test "X${am_cv_prog_cc_stdc}" = "Xno" ; then
	AC_MSG_ERROR([I'm afraid MBDyn requires compiler to support STDC constructs.])
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for header files.
dnl
AC_STDC_HEADERS
AC_CHECK_HEADERS(arpa/inet.h)
AC_CHECK_HEADERS(atomic_ops.h)
AC_CHECK_HEADERS(cxxabi.h)
AC_CHECK_HEADERS(fenv.h)
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS(netdb.h)
AC_CHECK_HEADERS(netinet/in.h)
AC_CHECK_HEADERS(pwd.h)
AC_CHECK_HEADERS(pthread.h)
AC_CHECK_HEADERS(semaphore.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_HEADERS(stdint.h)
AC_CHECK_HEADERS(stdio.h)
AC_CHECK_HEADERS(stdlib.h)
AC_CHECK_HEADERS(string.h)
AC_CHECK_HEADERS(sys/ipc.h)
AC_CHECK_HEADERS(sys/mman.h)
AC_CHECK_HEADERS(sys/param.h)
AC_CHECK_HEADERS(sys/poll.h)
AC_CHECK_HEADERS(sys/pstat.h)
AC_CHECK_HEADERS(sys/select.h)
AC_CHECK_HEADERS(sys/shm.h)
AC_CHECK_HEADERS(sys/socket.h)
AC_CHECK_HEADERS(sys/sysinfo.h)
AC_CHECK_HEADERS(sys/times.h)
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(time.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(values.h)


dnl ----------------------------------------------------------------
dnl
dnl Checks for C++ header files (takes care of ANSI/pre-ANSI naming)
dnl
AC_LANG_CPLUSPLUS

AC_CHECK_HEADERS(exception)
AC_CHECK_HEADERS(stdexcept)
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Sort out defines
dnl
dnl !!! set to "no" before release !!!
MBDYN_DEVEL=yes
dnl MBDYN_DEVEL=no
if test x"$MBDYN_DEVEL" = x"yes" ; then 
	AC_DEFINE(MBDYN_DEVEL,1,[define to enable features under development])
fi
AM_CONDITIONAL(MBDYN_DEVEL,[test "$MBDYN_DEVEL" = "yes"])

dnl ----------------------------------------------------------------
dnl
dnl Typedefs
dnl

dnl ----------------------------------------------------------------
dnl
dnl Set symbol table for Math Parser
dnl
AC_DEFINE(USE_TABLE,1,[define this to use external symbol table (default)])

dnl ----------------------------------------------------------------
dnl
dnl Set inter-language naming convention for calling Fortran routines
dnl
AC_DEFINE(USE_UNDERSCORE,1,[append underscore to fortran calls])

dnl ----------------------------------------------------------------
dnl
dnl Check for DEBUG
dnl
if test "$ol_enable_debug_mpi" != no ; then
	AC_DEFINE(DEBUG_MPI,1,[define this to add MPI debugging code])
fi
if test "$ol_enable_debug" != no -o "$ol_enable_debug_mpi" != no ; then
	AC_DEFINE(DEBUG,1,[define this to add debugging code])
	if test "$ol_with_debug_mode" = "mem" ; then
		AC_DEFINE(DEBUG_MEMMANAGER,1,[define this to add memory debugging code])
	fi
fi

dnl try to detect MPI and either Metis or Chaco
use_schur=no
if test "$ol_enable_schur" = yes ; then
	case "$ol_with_mpi" in
	no)
		AC_MSG_ERROR([Schur parallel solver needs MPI])
		;;

	auto)
		ol_with_mpi=yes
		;;
	esac

	dnl add further partitioning package detection code here
	if test "$ol_with_metis" = no \
		-a "$ol_with_chaco" = no \
	; then
		AC_MSG_WARN([No partitioning library available; Schur parallel solver will only allow manual partitioning.])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable ADAMS/View compatible output (no need to detect stuff)
dnl
if test "$ol_enable_adams" != no ; then
	AC_DEFINE(USE_ADAMS,1,[define this to enable ADAMS/View output])
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable MotionView compatible output (detect result library?)
dnl
if test "$ol_enable_motionview" != no ; then
	AC_CHECK_HEADER(mbs_result_api.h,[have_mbs_result_api_h=yes],[have_mbs_result_api_h=no])
	if test $have_mbs_result_api_h = "yes" ; then
		AC_DEFINE(USE_MOTIONVIEW,1,[define this to enable MotionView output])
		AC_LANG_CPLUSPLUS
		AC_TRY_COMPILE([
#include <mbs_result_api.h>
],[
			mrfOpenResult();
],		[mb_motionview_major_version=6],[mb_motionview_major_version=5])
		AC_LANG_C
		AC_DEFINE_UNQUOTED(PROBED_MOTIONVIEW_VER_MAJ,$mb_motionview_major_version,[set to the major version number of MotionView])
		if test $mb_motionview_major_version -lt 6 ; then
			AC_DEFINE(MRFOPENRESULT_NEEDS_ARGS,1,[define if MotionView API func mrfOpenResult() needs args])
		fi
	else
		if test $ol_enable_motionview != auto ; then
			AC_MSG_ERROR([Could not locate mbs_result_api.h header])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable NetCDF binary output
dnl
if test "$ol_enable_netcdf" != no ; then
	have_netcdf=no
	have_libnetcdf=no
	AC_LANG_CPLUSPLUS
	AC_CHECK_HEADER(netcdfcpp.h,[have_netcdfcpp_h=yes],[have_netcdfcpp_h=no])
	AC_LANG_C
	if test $have_netcdfcpp_h = "yes" ; then
		save_LIBS="$LIBS"
		try_NETCDF_LIBS="-lnetcdf_c++ -lnetcdf"
		LIBS="$try_NETCDF_LIBS $LIBS"
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <netcdfcpp.h>
			],[
	/* test for netcdf */
	NcFile nc("foobar", NcFile::Replace); // Create and leave in define mode
			],[
				have_libnetcdf=yes
				NETCDF_LIBS="$try_NETCDF_LIBS"
			])
		AC_LANG_C
		LIBS="$save_LIBS"
	fi

	if test $have_libnetcdf = yes ; then
		AC_DEFINE(USE_NETCDF,1,[define this to enable NetCDF based output])
	elif test $ol_enable_netcdf != auto ; then
		AC_MSG_ERROR([Could not locate netcdfcpp.h header])
	fi

	have_udunits=no
	have_libudunits=no
	AC_CHECK_HEADER(udunits.h,[have_udunits_h=yes],[have_udunits_h=no])
	if test $have_udunits_h = "yes" ; then
		save_LIBS="$LIBS"
		try_UDUNITS_LIBS="-ludunits"
		LIBS="$try_UDUNITS_LIBS $LIBS"
		AC_CHECK_LIB(udunits,utInit,[
			have_libudunits=yes
			UDUNITS_LIBS="$try_UDUNITS_LIBS"
		],[have_libudunits=no],[-lm])
		LIBS="$save_LIBS"
	fi

	if test $have_libudunits = yes ; then
		AC_DEFINE(USE_UDUNITS,1,[define this to enable udunits in NetCDF based output])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable Python support
dnl
use_python=no
if test "$ol_enable_python" != no ; then
	AX_PYTHON_DEVEL
	dnl ----------------------------------------------------------------
	dnl
	dnl detect and configure swig for python bindings
	dnl
	dnl AX_PKG_SWIG(1.3.17, [], [ AC_MSG_ERROR([SWIG required]) ])
	AX_PKG_SWIG

	if test x"$SWIG" = x -a "$ol_enable_python" = yes ; then
		AC_MSG_ERROR([python needs swig])
	fi
	if test "$ol_enable_python" = yes ; then
	        PYTHON_NUMPY_HEADERS=`python -c "import numpy; print numpy.get_include()" 2>/dev/null`
		if test x"$PYTHON_NUMPY_HEADERS" = x ; then
			AC_MSG_ERROR([python needs numpy - including its development packages])
		else
			PYTHON_NUMPY_HEADERS="-I$PYTHON_NUMPY_HEADERS"
		fi
	fi

	AX_SWIG_ENABLE_CXX
	AX_SWIG_MULTI_MODULE_SUPPORT
	AX_SWIG_PYTHON

	use_python=yes
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable Octave support <http://www.octave.org/>
dnl
use_octave=no
if test "$ol_enable_octave" != no ; then
	AC_LANG_CPLUSPLUS
	AC_CHECK_HEADER(octave/oct.h,[have_octave_oct_h=yes],[have_octave_oct_h=no])
	AC_CHECK_HEADER(octave/parse.h,[have_octave_parse_h=yes],[have_octave_parse_h=no],[
#include <octave/oct.h>
		])
	AC_CHECK_HEADER(octave/toplev.h,[have_octave_toplev_h=yes],[have_octave_toplev_h=no],[
#include <octave/oct.h>
		])
	AC_CHECK_HEADER(octave/octave.h,[have_octave_octave_h=yes],[have_octave_octave_h=no],[
#include <octave/oct.h>
		])
	AC_LANG_C

	have_mkoctfile=no
	which mkoctfile > /dev/null 2>&1
	RC=$?
	if test "$RC" = 0 ; then
		have_mkoctfile=yes
	fi

	link_octave="no"
	if test "$have_octave_oct_h" = "yes" \
		-a "$have_octave_parse_h" = "yes" \
		-a "$have_octave_toplev_h" = "yes" \
		-a "$have_octave_octave_h" = "yes" \
		-a "$have_mkoctfile" = "yes" ; \
	then
		save_LIBS="$LIBS"
		try_LIBS="-loctinterp -loctave -lcruft"
		LIBS="$try_LIBS $LIBS"
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>
			],[
		octave_main(0, NULL, 1);
			],[link_octave="yes"])
		AC_LANG_C
		LIBS="$save_LIBS"
		OCTAVE_LIBS="$try_LIBS"
	fi

	if test "$link_octave" = yes ; then
		use_octave=yes
		AC_DEFINE(USE_OCTAVE, 1, [define if octave can be used])

	else
		if test "$ol_enable_octave" = "yes" ; then
			AC_MSG_ERROR([octave support not available])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for working offsetof
dnl
have_offsetof=no
AC_MSG_CHECKING([for offsetof])
AC_TRY_COMPILE([#include <stddef.h>],[
	struct dummy_t { int dummy; };
	(void)offsetof(struct dummy_t, dummy);
],[have_offsetof=yes],[have_offsetof=no])
if test $have_offsetof != no ; then
	AC_DEFINE(HAVE_OFFSETOF,1, [define if offsetof() works])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for socket
dnl
AC_CHECK_FUNC(socket,[have_socket=yes],[have_socket=no])
use_socket=no
if test "$have_socket" = "yes" ; then
	AC_DEFINE(HAVE_SOCKET,1,[define if socket() is available])
	AC_DEFINE(USE_SOCKET,1,[define if socket() can be used])
	use_socket=yes

	if test $have_socket != no ; then
		AC_MSG_CHECKING([for socklen_t])
		AC_TRY_COMPILE([#include <sys/socket.h>],[
			socklen_t l;
		],[have_socklen_t=yes],[have_socklen_t=no])
		if test $have_socklen_t != no ; then
			AC_DEFINE(HAVE_SOCKLEN_T,1,
				[define if typedef socklen_t is defined])
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for sig_atomic_t
dnl
AC_MSG_CHECKING([for sig_atomic_t])
AC_TRY_COMPILE([#include <signal.h>],[
	sig_atomic_t	atom;
	],[
		AC_DEFINE(HAVE_SIG_ATOMIC_T,1,[define if sig_atomic_t is available])
		AC_MSG_RESULT([yes])
	],[
		AC_MSG_RESULT([no])
	])

dnl ----------------------------------------------------------------
dnl
dnl Check for bool
dnl
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether bool is defined])
AC_TRY_COMPILE([],[
	bool b;
],[have_bool=yes],[have_bool=no])
if test $have_bool != no ; then
	AC_DEFINE(HAVE_BOOL,1,[define if bool is defined])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
	AC_MSG_CHECKING([whether bool needs bool.h])
	AC_TRY_COMPILE([],[
	#include <bool.h>
	bool b;
	],[need_bool_h=yes],[need_bool_h=no])
	if test $need_bool_h != no ; then
	        AC_DEFINE(NEED_BOOL_H,1,[define if bool needs bool.h])
		AC_MSG_RESULT([yes])
	else
	        AC_MSG_RESULT([no])
	fi
fi
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Check for exceptions
dnl
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether the C++ compiler handles exceptions])
AC_TRY_COMPILE([],[
	try {
		(void *)0;
	}
	catch (...) {
		(void *)0;
	}
],[have_exceptions=yes],[have_exceptions=no])
if test "$have_exceptions" = "yes" ; then
	AC_MSG_RESULT([yes])
else
	AC_MSG_ERROR([C++ compiler does not handle exceptions])
fi
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Check for std namespace
dnl

AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether the C++ compiler allows to use namespace std])
AC_TRY_COMPILE([
	#include <math.h>
	using namespace std;
],[
	sqrt(2.);
],[allow_std=yes],[allow_std=no])
if test $allow_std != no ; then
	AC_DEFINE(USING_STD,1,[define if C++ compiler allows namespace std])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Check for UNIX fs !?!
dnl
case "$ol_with_fs" in
	unix)
		AC_DEFINE(USE_UNIX_FILE,1,
			[define this to use UNIX-style filenames])
		;;
	dnl dos)	AC_DEFINE(USE_DOS_FILE,1,
	dnl	[define this to use DOS-style filenames])
	dnl	;;
	*)
		AC_MSG_ERROR([only UNIX-stype filenames are allowed])
		;;
esac

dnl ----------------------------------------------------------------
dnl
dnl Check for threads (check anyway...)
dnl
have_threads="no"
if test "$ol_with_threads" != "no" ; then
	AC_MSG_CHECKING([for usable threads])
	if test "$ac_cv_header_pthread_h" = "yes" ; then
		save_LIBS="$LIBS"
		have_threads="no"
		if test "$have_threads" = no ; then
			AC_TRY_LINK([
			/* Check whether threads are linked by default */
#include <pthread.h>
				],[
		return pthread_create(0, 0, 0, 0);
				],[have_threads="yes"])
			LIBS="$save_LIBS"
		fi

		if test "$have_threads" = no ; then
			THREAD_LIBS="-lpthread"
			LIBS="$save_LIBS $THREAD_LIBS"
			AC_TRY_LINK([
			/* Check whether threads need -lpthread */
#include <pthread.h>
				],[
		return pthread_create(0, 0, 0, 0);
				],[have_threads="yes"])
			LIBS="$save_LIBS"
		fi

		if test "$have_threads" = "yes" ; then
			AC_DEFINE(HAVE_THREADS,1,[define if usable threads are available])
			AC_MSG_RESULT([yes])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for loadable
dnl
have_runtime_loading=no
if test "$ol_enable_runtime_loading" != "no" ; then
	AC_CHECK_HEADER(ltdl.h,[have_ltdl_h=yes],[have_ltdl_h=no])
	if test $have_ltdl_h = "yes" ; then
		AC_DEFINE(HAVE_LTDL_H,1,[define if you have <ltdl.h>])
		AC_CHECK_LIB(ltdl, lt_dlinit,[have_lt_dlinit=yes],[have_lt_dlinit=no])
		if test $have_lt_dlinit = "yes" ; then
			LIBS="$LIBS -lltdl"
			have_runtime_loading=yes
			MODULE_LOADPATH=$libexecdir
		fi
	fi

	if test "$have_runtime_loading" = "yes" ; then
		AC_DEFINE(USE_RUNTIME_LOADING,1,[define if you want to load runtime])
	else
		if test "$ol_enable_runtime_loading" = "yes" ; then
			AC_MSG_ERROR([no dynamic linking capability detected])
		else
			AC_MSG_WARN([no dynamic linking capability detected; runtime loading disabled])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl static modules...
dnl
build_static_modules=no
case "$ol_with_static_modules" in
	no)
		if test "$have_runtime_loading" != "yes" -a test "x$MODULES_DIRS" != "x" ; then
			dnl ?!? check
			AC_MSG_ERROR([modules cannot be built as dynamic because no runtime loading capability is available])
		fi
		;;
	auto)
		if test "$have_runtime_loading" != "yes" ; then
			build_static_modules=yes
		fi
		;;
	yes)
		build_static_modules=yes
		;;
	*)
		;;
esac
if test "$build_static_modules" = "yes" ; then
	AC_DEFINE(STATIC_MODULES,1,[define to statically build modules])
fi

build_charm=no
if test "$ol_with_charm" != "no" ; then
	build_charm=yes
	AC_DEFINE(HAVE_CHARM,1,[define to build CHARM WP module])
	dnl add test
	if test 1 = 1 ; then
		AC_DEFINE(NEED_CHARM_SHIPWAKE,1,[define if CHARM WP module needs shipwake])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for multithread assembly (use previous checks)
dnl
ol_link_multithread=no
if test "$ol_enable_multithread" != no ; then
	if test "$ol_enable_multithread" = "yes" -a "$ol_with_thread" = "no" ; then
		AC_MSG_ERROR([Can't enable multithread --without-threads])
	fi
	AC_MSG_CHECKING([whether to enable multithread])

	if test "$ac_cv_header_pthread_h" != "yes" ; then
		AC_MSG_ERROR([need pthreads.h])
	fi

	if test "$ac_cv_header_atomic_ops_h" != "yes" ; then
		AC_MSG_ERROR([need atomic_ops.h (from HP's atomic_ops project)])
	fi

	if test "$have_threads" = no ; then
		case "$ol_enable_multithread" in
		yes)
			AC_MSG_ERROR([No threads available])
			;;
			
		force)
			have_threads=yes
			;;

		*)
			;; dnl do nothing ...
		esac
	fi

	if test "$have_threads" = no ; then
		AC_MSG_RESULT([no])
	else
		ol_link_multithread=yes
		AC_DEFINE(USE_MULTITHREAD,1,[define to enable multithread assembly])
		AC_MSG_RESULT([yes])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for multithread naive solver (use previous checks)
dnl
ol_link_multithread_naive=no
if test "$ol_enable_multithread_naive" != no ; then
	if test "$ol_enable_multithread_naive" = "yes" -a "$ol_with_thread" = "no" ; then
		AC_MSG_ERROR([Can't enable multithread naive solver --without-threads])
	fi
	AC_MSG_CHECKING([whether to enable multithread naive solver])

	if test "$ac_cv_header_pthread_h" != "yes" ; then
		AC_MSG_ERROR([need pthreads.h])
	fi

	if test "$ac_cv_header_atomic_ops_h" != "yes" ; then
		AC_MSG_ERROR([need atomic_ops.h (from HP's atomic_ops project)])
	fi

	if test "$have_threads" = no ; then
		case "$ol_enable_multithread_naive" in
		yes)
			AC_MSG_ERROR([No threads available])
			;;
			
		force)
			have_threads=yes
			;;

		*)
			;; dnl do nothing ...
		esac
	fi

	if test "$have_threads" = no ; then
		AC_MSG_RESULT([no])
	else
		if test "$ac_cv_header_asm_atomic_h" = no \
				-o "$ac_cv_header_asm_bitops_h" = no \
				-o "$ac_cv_header_asm_system_h" = no \
				; then
			AC_MSG_ERROR([missing required system headers])
		fi
		ol_link_multithread_naive=yes
		AC_DEFINE(USE_NAIVE_MULTITHREAD,1,[define to enable multithread naive solver])
		AC_MSG_RESULT([yes])
	fi
fi

ol_link_multithread=no
if test "$ol_enable_multithread" != no ; then
	AC_MSG_CHECKING([for usable threads])

	if test "$ac_cv_header_pthread_h" != "yes" ; then
		AC_MSG_ERROR([need pthreads.h])
	fi

	if test "$ac_cv_header_atomic_ops_h" != "yes" ; then
		AC_MSG_ERROR([need atomic_ops.h (from HP's atomic_ops project)])
	fi

	save_LIBS="$LIBS"
	have_threads="no"
	if test "$have_threads" = no ; then
		AC_TRY_LINK([
		/* Check whether threads are linked by default */
#include <pthread.h>
			],[
	return pthread_create(0, 0, 0, 0);
			],[have_threads="yes"])
		LIBS="$save_LIBS"
	fi

	if test "$have_threads" = no ; then
		THREAD_LIBS="-lpthread"
		LIBS="$save_LIBS $THREAD_LIBS"
		AC_TRY_LINK([
		/* Check whether threads need -lpthread */
#include <pthread.h>
			],[
	return pthread_create(0, 0, 0, 0);
			],[have_threads="yes"])
		LIBS="$save_LIBS"
	fi
	
	if test "$have_threads" = no ; then
		AC_MSG_RESULT([no])
		case "$ol_enable_multithread" in
		yes)
			AC_MSG_ERROR([No threads available])
			;;
			
		force)
			have_threads=yes
			;;

		*)
			;; dnl do nothing ...
		esac
	fi

	if test "$have_threads" != no ; then
		ol_link_multithread=yes
		AC_DEFINE(USE_MULTITHREAD,1,[define to enable multithreaded solution])
		AC_MSG_RESULT([yes])
	fi
fi

if test "$ol_link_multithread" = "yes" -o "$ol_link_multithread_naive" = "yes" ; then
	SAVE_LIBS="$LIBS"
	LIBS="$LIBS -latomic_ops_gpl -latomic_ops"
	AC_TRY_LINK([#include <atomic_ops.h>], [AO_t t; AO_TS_t tt; AO_compare_and_swap_full(&t, t, t); AO_test_and_set_full(&tt)],[have_atomic_ops=yes],[have_atomic_ops=no],[-latomic_ops])	LIBS="$SAVE_LIBS"
	if test $have_atomic_ops = "yes" ; then
		AC_DEFINE(HAVE_ATOMIC_OPS,1,[define to use HP's atomic_ops library])
	else
		AC_MSG_ERROR([Need HP's atomic_ops library for multithreaded features])
		ATOMIC_OPS_LIBS="-latomic_ops_gpl -latomic_ops"
	fi
fi



dnl ----------------------------------------------------------------
dnl
dnl Check for /dev/TASK2CPU
dnl
AC_MSG_CHECKING([for /dev/TASK2CPU])
if test -e "/dev/TASK2CPU" ; then
	AC_DEFINE(HAVE_TASK2CPU,1,[define if /dev/TASK2CPU is available])
	AC_MSG_RESULT(yes)
else
	AC_MSG_RESULT(no)
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for BOOST
dnl
ol_link_boost=no
if test $ol_with_boost != no ; then
	AC_LANG_CPLUSPLUS
	AC_CHECK_HEADERS(boost/config.hpp)
	if test "$ac_cv_header_boost_config_hpp" = "yes" ; then
		AC_DEFINE(USE_BOOST,1,[define to use BOOST])
		AC_CHECK_HEADERS(boost/graph/adjacency_list.hpp)
		AC_CHECK_HEADERS(boost/graph/cuthill_mckee_ordering.hpp)
		AC_CHECK_HEADERS(boost/graph/king_ordering.hpp)
		AC_CHECK_HEADERS(boost/graph/sloan_ordering.hpp)
		AC_CHECK_HEADERS(boost/graph/properties.hpp)
		AC_CHECK_HEADERS(boost/graph/bandwidth.hpp)
		AC_CHECK_HEADERS(boost/graph/wavefront.hpp)
		AC_CHECK_HEADERS(boost/graph/minimum_degree_ordering.hpp)
	else
		if test $ol_with_boost != auto ; then
			AC_MSG_ERROR([BOOST library headers not available])
		fi
	fi
	AC_LANG_C
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for bullet <http://bulletphysics.com>
dnl
ol_link_bullet=no
if test "$ol_with_bullet" != "no" ; then
	AC_LANG_CPLUSPLUS
	have_libbullet=no
	AC_CHECK_HEADER(BulletCollision/CollisionShapes/btCollisionShape.h,[have_bullet_h=yes],[have_bullet_h=no])
	if test "$have_bullet_h" = "yes" ; then
		SAVE_LIBS="$LIBS"
		try_BULLET_LIBS="-lbulletcollision -lbulletmath"
		LIBS="$try_BULLET_LIBS $LIBS"
		AC_CHECK_LIB(bulletcollision,btBulletCollisionProbe,[
			have_libbullet=yes
			BULLET_LIBS="$try_BULLET_LIBS"
		],[have_libbullet=no],[-lbulletmath])
		LIBS="$save_LIBS"
	fi

	if test $have_libbullet = yes ; then
		AC_DEFINE(USE_BULLET,1,[define this to enable bullet <http://bulletphysics.com>])
	else
		if test "$ol_with_bullet" = "yes" ; then
			AC_MSG_ERROR([bullet not available])
		fi
	fi
	AC_LANG_C
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for BLAS
dnl
ol_link_blas=no
if test $ol_with_blas != no ; then
	dnl FIXME: any standard BLAS headers?
	AC_MSG_CHECKING([for BLAS library])

	save_LIBS="$LIBS"
	have_blas="no"

	if test "${ol_with_goto_lib:-set}" = set ; then
		dnl in case one relies on a symbolic link named libgoto.{a,so}
		ol_with_goto_lib="goto"
	fi

	try_blas_LIBS=
	if test "$have_blas" = "no" ; then
		try_blas_LIBS="-l$ol_with_goto_lib -lxerbla $FCLIBS $F77LIBS -lm"
		if test "$ol_with_blas" = "auto" -o "$ol_with_blas" = "goto" ; then
			LIBS="$try_blas_LIBS $LIBS"

			AC_TRY_LINK([
	/* FIXME: some architectures might use different
	 * naming conventions, e.g. uppercase, no underscore
	 */
	extern int dgemm_(void);
			],[
	/* test for blas */
	return dgemm_();
			],[have_blas="goto"],[
				have_blas="no"
				if test "$ol_with_blas" = "goto" ; then
					AC_MSG_ERROR([Goto blas implementation not available])
				fi
			])
			LIBS="$save_LIBS"
		fi
	fi

	if test "$have_blas" = "no" ; then
		try_blas_LIBS="-lcblas -latlas -lm"
		if test "$ol_with_blas" = "auto" -o "$ol_with_blas" = "atlas" ; then
			LIBS="$try_blas_LIBS $LIBS"

			AC_TRY_LINK([
	extern int cblas_dgemm(void);
			],[
	/* test for blas */
	return cblas_dgemm();
			],[have_blas="atlas"
				LIBS="-lcblas -lf77blas -latlas $FCLIBS $F77LIBS -lm $LIBS"
				AC_TRY_LINK([
					/* FIXME: some architectures might use different
					 * naming conventions, e.g. uppercase, no underscore
					 */
					extern int dgemm_(void);],[return dgemm_();
				],[try_blas_LIBS="-lcblas -lf77blas -latlas $FCLIBS $F77LIBS -lm"],[])
			],[
				have_blas="no"
				if test "$ol_with_blas" = "atlas" ; then
					AC_MSG_ERROR([Atlas blas implementation not available])
				fi
			])
			LIBS="$save_LIBS"
		fi
	fi

	if test "$have_blas" = "no" ; then
		try_blas_LIBS="-lcblas -lm"
		if test "$ol_with_blas" = "auto" -o "$ol_with_blas" = "cblas" ; then
			LIBS="$try_blas_LIBS $LIBS"

			AC_TRY_LINK([
	extern int cblas_dgemm(void);
			],[
	/* test for blas */
	return cblas_dgemm();
			],[have_blas="cblas"],[
				have_blas="no"
				if test "$ol_with_blas" = "cblas" ; then
					AC_MSG_ERROR([CBLAS not available])
				fi
			])
			LIBS="$save_LIBS"
		fi
	fi

	if test "$have_blas" = "no" ; then
		try_blas_LIBS="-lblas $FCLIBS $F77LIBS -lm"
		if test "$ol_with_blas" = "auto" -o "$ol_with_blas" = "blas" ; then
			LIBS="$try_blas_LIBS $LIBS"

			AC_TRY_LINK([
	/* FIXME: some architectures might use different
	 * naming conventions, e.g. uppercase, no underscore
	 */
	extern int dgemm_(void);
			],[
	/* test for blas */
	return dgemm_();
			],[have_blas="blas"],[
				have_blas="no"
				if test "$ol_with_blas" = "blas" ; then
					AC_MSG_ERROR([BLAS not available])
				fi
			])
			LIBS="$save_LIBS"
		fi
	fi

	dnl revert to original libraries; the BLAS libraries 
	dnl will be explicitly added by those libraries that require them
	if test "$have_blas" != "no" ; then
		ol_link_blas="yes"
		BLAS_LIBS="$try_blas_LIBS"
		AC_DEFINE(HAVE_BLAS,1,[define if BLAS are available])
	fi

	dnl from now on, add $BLAS_LIBS to libraries that (may) depend on them
	AC_MSG_RESULT([$have_blas])
fi


dnl ----------------------------------------------------------------
dnl
dnl Need either Metis or Chaco with MPI and Taucs
dnl
dnl First we check for Metis
ol_link_metis=no
if test "$ol_with_metis" != "no" ; then
	AC_CHECK_LIB(metis, METIS_PartGraphVKway,
		[ol_link_metis=yes],[ol_link_metis=no],[-lm])
	
	AC_CHECK_HEADERS(metis.h)

	if test "$ac_cv_header_metis_h" = "yes" ; then
		if test "$ol_link_metis" != "no" ; then
			AC_DEFINE(USE_METIS,1,
				[define if you want to use Metis])
		fi
	fi

	if test "$ol_link_metis" = "no" -o "$ac_cv_header_metis_h" != "yes" ; then
		if test "$ol_with_metis" != "auto" ; then
			AC_MSG_ERROR([Could not locate Metis])
		else
			AC_MSG_WARN([Metis not supported])
			dnl set to "no" because we cannot link metis
			dnl without header file
			ol_link_metis="no"
		fi
	else
		METIS_LIBS="-lmetis -lm"
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Now we check for Chaco.
dnl
ol_link_chaco=no
if test "$ol_with_chaco" != no ; then
	AC_CHECK_LIB(chaco,interface,[have_chaco=yes],[have_chaco=no],[-lm])

	if test $have_chaco != no ; then
		CHACO_LIBS="-lchaco -lm"
		AC_DEFINE(USE_CHACO,1,[define if you want to use Chaco])
		ol_link_chaco=yes
	fi

	if test "$ol_link_chaco" = no ; then
                if test $ol_with_chaco != auto ; then
                        AC_MSG_ERROR([Could not locate Chaco])
                else
                        AC_MSG_WARN([Chaco not supported])
                fi
        fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for LAPACK
dnl
ol_link_lapack=no
if test $ol_with_lapack != no ; then
	AC_MSG_CHECKING([for lapack])

	save_LIBS="$LIBS"
	try_lapack_LIBS="-llapack"
	LIBS="$try_lapack_LIBS $BLAS_LIBS $FCLIBS $F77LIBS $save_LIBS"
	AC_TRY_LINK([
	extern double dgesvd(double *);
		],[
	double i = 1., o;
	o = dgesvd(&i);
		],[have_lapack=yes
		   have_lapack_=no],[
		AC_TRY_LINK([
	extern double dgesvd_(double *);
			],[
	double i = 1., o;
	o = dgesvd_(&i);
		],[have_lapack=yes
		   have_lapack_=yes],[have_lapack=no
				      have_lapack_=no])])
	
	if test $have_lapack != no ; then
		AC_DEFINE(USE_LAPACK,1,[define to use LAPACK libraries])
		AC_DEFINE(HAVE_DGESVD,1,[define to use LAPACK's dgesvd in GPC])
		LAPACK_LIBS="$try_lapack_LIBS"
	fi

	if test "$have_lapack_" = "yes" ; then
		AC_TRY_LINK([
	extern double dggev_(double *);
			],[
	double i = 1., o;
	o = dggev_(&i);
			],[AC_DEFINE(HAVE_DGGEV,1,[define to use LAPACK's dggev in Eig])])
	elif test "$have_lapack" = "yes" ; then
		AC_TRY_LINK([
	extern double dggev(double *);
			],[
	double i = 1., o;
	o = dggev(&i);
			],[AC_DEFINE(HAVE_DGGEV_,1,[define to use LAPACK's dggev in Eig])])
	fi

	if test "$have_lapack_" = "yes" ; then
		AC_TRY_LINK([
	extern double dgegv_(double *);
			],[
	double i = 1., o;
	o = dgegv_(&i);
			],[AC_DEFINE(HAVE_DGEGV,1,[define to use LAPACK's dgegv in Eig])])
	elif test "$have_lapack" = "yes" ; then
		AC_TRY_LINK([
	extern double dgegv(double *);
			],[
	double i = 1., o;
	o = dgegv(&i);
			],[AC_DEFINE(HAVE_DGEGV_,1,[define to use LAPACK's dgegv in Eig])])
	fi

	if test "$have_lapack" = "yes" ; then
		ol_link_lapack=yes
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])	
	fi
	LIBS="$save_LIBS"
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Umfpack
dnl
ol_link_umfpack=no
if test "$ol_with_umfpack" != "no" ; then
	AC_CHECK_HEADERS(umfpack.h)

	if test "$ac_cv_header_umfpack_h" = "yes" ; then

		dnl check for AMD
		AC_MSG_CHECKING([for AMD library])
		save_LIBS="$LIBS"
		try_amd_LIBS="-lamd -lm"
		LIBS="$try_amd_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for amd */
	return 0;
		],[have_amd="yes"],[have_amd="no"])
		if test "$have_amd" = "yes" ; then
			try_umfpack_LIBS="$try_amd_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_umfpack_LIBS=
			AC_MSG_RESULT([no])
		fi


		dnl check for COLAMD
		AC_MSG_CHECKING([for COLAMD library])
		try_colamd_LIBS="-lcolamd -lm"
		LIBS="$try_colamd_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for colamd */
	return 0;
		],[have_colamd="yes"],[have_colamd="no"])
		if test "$have_colamd" = "yes" ; then
			try_umfpack_LIBS="$try_colamd_LIBS $try_umfpack_LIBS"
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi


		dnl check for CCOLAMD
		AC_MSG_CHECKING([for CCOLAMD library])
		try_ccolamd_LIBS="-lccolamd -lm"
		LIBS="$try_ccolamd_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for ccolamd */
	return 0;
		],[have_ccolamd="yes"],[have_ccolamd="no"])
		if test "$have_ccolamd" = "yes" ; then
			try_umfpack_LIBS="$try_ccolamd_LIBS $try_umfpack_LIBS"
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi


		dnl check for CAMD
		AC_MSG_CHECKING([for CAMD library])
		try_camd_LIBS="-lcamd -lm"
		LIBS="$try_camd_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for camd */
	return 0;
		],[have_camd="yes"],[have_camd="no"])
		if test "$have_camd" = "yes" ; then
			try_umfpack_LIBS="$try_camd_LIBS $try_umfpack_LIBS"
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi


		dnl check for CHOLMOD
		AC_MSG_CHECKING([for CHOLMOD library])
		try_cholmod_LIBS="-lcholmod $METIS_LIBS -lm"
		LIBS="$try_cholmod_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for cholmod */
	return 0;
		],[have_cholmod="yes"],[have_cholmod="no"])
		if test "$have_cholmod" = "yes" ; then
			try_umfpack_LIBS="$try_cholmod_LIBS $try_umfpack_LIBS"
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi

		have_umfpack_version="no"
		have_umfpack="no"
		AC_MSG_CHECKING([for umfpack library version])
		AC_TRY_COMPILE([
#include <sys/types.h>
#include <umfpack.h>
		],[
#ifndef UMFPACK_MAIN_VERSION
#error "UMFPACK_MAIN_VERSION not defined (less than 5?)"
#endif // UMFPACK_MAIN_VERSION
#if UMFPACK_MAIN_VERSION < 5
#error "UMFPACK_MAIN_VERSION less than 5?"
#endif // UMFPACK_MAIN_VERSION
		],[
			have_umfpack_version="yes"
			AC_MSG_RESULT([5 or greater])
		],[
			AC_MSG_RESULT([less than 5?])
		])

		if test "$have_umfpack_version" = "yes" ; then
			umfpack_use_long="no"
			AC_MSG_CHECKING([if umfpack library needs long])
			AC_LANG_CPLUSPLUS
			AC_TRY_COMPILE([
#include <sys/types.h>
#include <umfpack.h>
			],[
#define SIZE 2
	int32_t sys;
	int32_t Ap[SIZE], Ai[SIZE];
	double Ax[SIZE], X[SIZE], B[SIZE];
	void *Numeric;
	double Control[UMFPACK_CONTROL];
	double Info[UMFPACK_INFO];
	umfpack_di_solve(sys, Ap, Ai, Ax, X, B, &Numeric, Control, Info);
			],[],[
				AC_TRY_COMPILE([
#include <sys/types.h>
#include <umfpack.h>
				],[
#define SIZE 2
	int32_t sys;
	int32_t Ap[SIZE], Ai[SIZE];
	double Ax[SIZE], X[SIZE], B[SIZE];
	void *Numeric;
	double Control[UMFPACK_CONTROL];
	double Info[UMFPACK_INFO];
	umfpack_dl_solve(sys, Ap, Ai, Ax, X, B, &Numeric, Control, Info);
				],[
					AC_DEFINE(USE_UMFPACK_LONG,1,[define if need to use umfpack_dl_* calls])
					umfpack_use_long="yes"
				],[
					AC_MSG_ERROR([umfpack might be broken])
				]) 
			])
			AC_LANG_C

			AC_MSG_RESULT([$umfpack_use_long])

			try_umfpack_LIBS="-lumfpack $try_umfpack_LIBS"
			LIBS="$try_umfpack_LIBS $BLAS_LIBS $save_LIBS"
			AC_MSG_CHECKING([for umfpack library])
			AC_TRY_LINK([
#include <umfpack.h>
			],[
	double t;
	t = umfpack_timer();
			],[
				have_umfpack="yes"
				AC_MSG_RESULT([yes])
			],[
				AC_MSG_RESULT([no])
			])
		fi

		if test "$have_umfpack" = "yes" ; then
			ol_link_umfpack=yes
			UMFPACK_LIBS="$try_umfpack_LIBS"
			if test "${BLAS_LIBS:-set}" = "set" ; then
				AC_DEFINE(NBLAS,1,[define if BLAS/ATLAS are not available to Umfpack])
				AC_MSG_WARN([Umfpack will be used without cblas/atlas])
			fi
		fi
		LIBS="$save_LIBS"
	fi
fi

if test $ol_link_umfpack = no ; then
	case "$ol_with_umfpack" in
	yes) 
		AC_MSG_ERROR([Could not locate Umfpack])
		;;
	auto)
		AC_MSG_WARN([Umfpack not supported])
		;;
	*)
		;;
	esac
else
	AC_DEFINE(USE_UMFPACK,1,[define if you want to use Umfpack])
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for KLU
dnl
ol_link_klu=no
if test "$ol_with_klu" != "no" ; then
	dnl OL_ADD_DIRS(klu)
	
	AC_CHECK_HEADERS(klu.h)

	if test "$ac_cv_header_klu_h" = "yes" ; then

		dnl check for AMD
		AC_MSG_CHECKING([for AMD library])
		save_LIBS="$LIBS"
		try_amd_LIBS="-lamd"
		LIBS="$try_amd_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for amd */
	return 0;
		],[have_amd="yes"],[have_amd="no"])
		if test "$have_amd" = "yes" ; then
			try_klu_LIBS="$try_amd_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_klu_LIBS=
			AC_MSG_RESULT([no])
		fi

		dnl check for COLAMD
		AC_MSG_CHECKING([for COLAMD library])
		try_colamd_LIBS="-lcolamd"
		LIBS="$try_colamd_LIBS $save_LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for colamd */
	return 0;
		],[have_colamd="yes"],[have_colamd="no"])
		if test "$have_colamd" = "yes" ; then
			try_klu_LIBS="$try_klu_LIBS $try_colamd_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_klu_LIBS=
			AC_MSG_RESULT([no])
		fi

		dnl check for BTF
		AC_MSG_CHECKING([for BTF library])
		try_btf_LIBS="-lbtf"
		LIBS="$try_btf_LIBS $save_LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for btf */
	return 0;
		],[have_btf="yes"],[have_btf="no"])
		if test "$have_btf" = "yes" ; then
			try_klu_LIBS="$try_klu_LIBS $try_btf_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_klu_LIBS=
			AC_MSG_RESULT([no])
		fi

		dnl check for KLU
		AC_MSG_CHECKING([for KLU library])
		try_klu_LIBS="-lklu $try_klu_LIBS"
		LIBS="$try_klu_LIBS $save_LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for klu */
	return 0;
		],[have_klu="yes"],[have_klu="no"])
		if test "$have_klu" = "yes" ; then
			AC_MSG_RESULT([yes])
		else
			try_klu_LIBS=
			AC_MSG_RESULT([no])
		fi

		if test "$have_klu" = "yes" ; then
			ol_link_klu=yes
			KLU_LIBS="$try_klu_LIBS"
			if test "${BLAS_LIBS:-set}" = "set" ; then
				AC_DEFINE(NBLAS,1,[define if BLAS/ATLAS are not available to KLU])
				AC_MSG_WARN([KLU will be used without cblas/atlas])
			fi
		fi
		LIBS="$save_LIBS"

	fi
fi
		

if test $ol_link_klu = no ; then
	case "$ol_with_klu" in
	yes) 
		AC_MSG_ERROR([Could not locate KLU])
		;;
	auto)
		AC_MSG_WARN([KLU not supported])
		;;
	*)
		;;
	esac
else
	AC_DEFINE(USE_KLU,1,[define if you want to use KLU])
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Meschach
dnl
ol_link_meschach=no
if test $ol_with_meschach != no ; then
	AC_CHECK_HEADERS(meschach/sparse2.h)

	if test "$ac_cv_header_meschach_sparse2_h" = "yes" ; then

		AC_MSG_CHECKING([for meschach library])
		try_meschach_LIBS="-lmeschach"
		save_LIBS=$LIBS
		LIBS="$LIBS $try_meschach_LIBS"
		AC_TRY_LINK([
#include <meschach/sparse2.h>
],[
	SPMAT* mat;
	mat = sp_get(10, 20, 10);
		],[have_meschach=yes],[have_meschach=no])
		LIBS=$save_LIBS

		if test $have_meschach != no ; then
			MESCHACH_LIBS=$try_meschach_LIBS
			AC_DEFINE(USE_MESCHACH,1,[define if you want to use Meschach])
			ol_link_meschach=yes
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi
	fi

	if test $ol_link_meschach = no ; then
		if test $ol_with_meschach != auto ; then
			AC_MSG_ERROR([Could not locate Meschach])
		else
			AC_MSG_WARN([Meschach not supported])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for SuperLU
dnl
ol_link_superlu=no
if test $ol_with_superlu != no ; then
	AC_CHECK_HEADERS(dsp_defs.h)
	AC_CHECK_HEADERS(pdsp_defs.h)
	AC_CHECK_HEADERS(util.h)

	if test "$ac_cv_header_util_h" = "yes" ; then
		AC_MSG_CHECKING([for SuperLU library])
		try_superlu_LIBS="-lsuperlu"
		save_LIBS="$LIBS"
		have_superlu_mt="no"
		if test "$ac_cv_header_pdsp_defs_h" = "yes" ; then
			LIBS="$try_superlu_LIBS $save_LIBS $BLAS_LIBS $THREAD_LIBS -lm"
			AC_TRY_LINK([
#include <pdsp_defs.h>
#include <util.h>
],[
	SuperMatrix	AC, L, U;
	(void)pdgstrf_thread_init(&AC, &L, &U,
			(void *)0, (void *)0, (void *)0, (void *)0);
			],[have_superlu_mt="yes"])
		fi
		have_superlu="no"
		if test "$ac_cv_header_dsp_defs_h" = "yes" ; then
			LIBS="$try_superlu_LIBS $save_LIBS $BLAS_LIBS -lm"
			AC_TRY_LINK([
#include <dsp_defs.h>
#include <util.h>
],[
	SuperMatrix	B;
	dCreate_Dense_Matrix(&B, 1, 1, (void *)0, 0,
			SLU_DN, SLU_D, SLU_GE);
			],[have_superlu="yes"])
		fi
		LIBS="$save_LIBS"

		if test "$have_superlu_mt" != "no" -o \
				"$have_superlu" != "no" ; then
			SUPERLU_LIBS=$try_superlu_LIBS
			ol_link_superlu=yes
			if test "$have_superlu_mt" != "no" ; then
				AC_DEFINE(USE_SUPERLU_MT,1,[define if you want to use multithreaded SuperLU])
			fi
			AC_DEFINE(USE_SUPERLU,1,[define if you want to use SuperLU])

			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi
	fi

	if test $ol_link_superlu = no ; then
		if test $ol_with_superlu != auto ; then
			AC_MSG_ERROR([Could not locate SuperLU])
		else
			AC_MSG_WARN([SuperLU not supported])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Watson wsmp
dnl
ol_link_wsmp=no
if test $ol_with_wsmp != no ; then
	try_wsmp_LIBS="-lwsmp"
	AC_CHECK_LIB(wsmp,wgsmp_,[have_wsmp=yes],[have_wsmp=no],[$try_wsmp_LIBS])
	if test "$have_wsmp" = "yes" ; then
		AC_DEFINE(HAVE_WSMP,1,[define to use Watson wsmp library])
		WSMP_LIBS="$try_wsmp_LIBS"
	else
		if test "$ol_with_wsmp" = "yes" ; then
			AC_MSG_ERROR([wsmp library not available])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Y12
dnl
ol_link_y12=no
if test $ol_with_y12 != no ; then
	ol_link_y12=yes
	AC_DEFINE(USE_Y12,1,[define to use Y12 sparse LU library])
	Y12_LIBS="${top_builddir}/libraries/liby12/liby12.la -lm"
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Harwell
dnl
ol_link_harwell=no
if test $ol_with_harwell != no ; then
	AC_CHECK_LIB(hsl,ma28ad,[have_libhsl=yes],[have_libhsl=no],["$FCLIBS $F77LIBS -lm"])
	if test $have_libhsl = no ; then
		AC_CHECK_LIB(hsl,ma28ad_,[have_libhsl=yes],[have_libhsl=no],["$FCLIBS $F77LIBS -lm"])
	fi
	if test $have_libhsl = yes ; then
		ol_link_harwell=yes
		AC_DEFINE(USE_HARWELL,1,[define to use Harwell sparse LU library])
		HARWELL_LIBS=-lhsl
	else
		if test $ol_with_harwell != auto ; then
			AC_MSG_ERROR([Could not locate Harwell lib (-lhsl)])
		else
			AC_MSG_WARN([Harwell lib (-lhsl) not supported])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for ARPACK
dnl
ol_link_arpack=no
if test $ol_with_arpack != no ; then
	AC_MSG_CHECKING([for arpack])

	save_LIBS="$LIBS"
	try_arpack_LIBS="-larpack"
	LIBS="$try_arpack_LIBS $LAPACK_LIBS $BLAS_LIBS $FCLIBS $F77LIBS $save_LIBS"
	AC_TRY_LINK([
	extern double dnaupd(double *);
		],[
	double i = 1., o;
	o = dnaupd(&i);
		],[have_arpack=yes
		   have_arpack_=no],[
		AC_TRY_LINK([
	extern double dnaupd_(double *);
			],[
	double i = 1., o;
	o = dnaupd_(&i);
		],[have_arpack=yes
		   have_arpack_=yes],[have_arpack=no
				      have_arpack_=no])])
	
	if test $have_arpack != no ; then
		AC_DEFINE(USE_ARPACK,1,[define to use ARPACK libraries])
		ARPACK_LIBS="$try_arpack_LIBS"
	fi

	if test "$have_arpack" = "yes" ; then
		ol_link_arpack=yes
		AC_MSG_RESULT([yes])
	else
		if test "$ol_with_arpack" = "yes" ; then
			AC_MSG_ERROR([Could not locate arpack])
		fi
		AC_MSG_RESULT([no])
	fi
	LIBS="$save_LIBS"
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for JDQZ
dnl
ol_link_jdqz=no
if test $ol_with_jdqz != no ; then
	AC_MSG_CHECKING([for jdqz])

	save_LIBS="$LIBS"
	try_jdqz_LIBS="-ljdqz"
	LIBS="$try_jdqz_LIBS $LAPACK_LIBS $BLAS_LIBS $FCLIBS $F77LIBS $save_LIBS"
	AC_TRY_LINK([
	double amul(double *x) { return 0.; }
	double bmul(double *x) { return 0.; }
	double precon(double *x) { return 0.; }
	extern double jdqz(double *);
		],[
	double i = 1., o;
	o = jdqz(&i);
		],[have_jdqz=yes
		   have_jdqz_=no],[
		AC_TRY_LINK([
	double amul_(double *x) { return 0.; }
	double bmul_(double *x) { return 0.; }
	double precon_(double *x) { return 0.; }
	extern double jdqz_(double *);
			],[
	double i = 1., o;
	o = jdqz_(&i);
		],[have_jdqz=yes
		   have_jdqz_=yes],[have_jdqz=no
				      have_jdqz_=no])])
	
	if test $have_jdqz != no ; then
		AC_DEFINE(USE_JDQZ,1,[define to use JDQZ library])
		JDQZ_LIBS="$try_jdqz_LIBS"
	fi

	if test "$have_jdqz" = "yes" ; then
		ol_link_jdqz=yes
		AC_MSG_RESULT([yes])
	else
		if test "$ol_with_jdqz" = "yes" ; then
			AC_MSG_ERROR([Could not locate jdqz])
		fi
		AC_MSG_RESULT([no])
	fi
	LIBS="$save_LIBS"
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for eig
dnl
if test $ol_enable_eig != no ; then
	if test $ol_link_arpack = yes -o $ol_link_lapack = yes -o $ol_link_jdqz = yes ; then
		AC_DEFINE(USE_EIG,1,[define if you want to enable direct eigensolution])
	else
		if test $ol_enable_eig = yes ; then
			AC_MSG_ERROR([Need either LAPACK, ARPACK or JDQZ to enable direct eigenanalysis])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for ANN
dnl
ol_link_ann=no
if test "$ol_with_ann" != "no" ; then
	ol_have_ann=no
	AC_LANG_CPLUSPLUS
	ANN_H=
	ANN_ANN_H=
	AC_CHECK_HEADERS(ANN.h, [ANN_H="ANN.h"])
	AC_CHECK_HEADERS(ANN/ANN.h, [ANN_H="ANN/ANN.h"])
	if test "$ac_cv_header_ANN_h" = "yes" -o "$ac_cv_header_ANN_ANN_h" = "yes" ; then
		SAVE_LIBS="$LIBS"
		LIBS="$LIBS -lann"
		AC_TRY_LINK(
	#include <$ANN_H>
			,[
	ANNkd_tree x(0,0);
			],[have_ann=yes],[have_ann=no])
		LIBS="$save_LIBS"
		if test "$have_ann" != "no" ; then
			ANN_LIBS="-lann"
			AC_DEFINE(HAVE_ANN,1,[define to use ANN])
			ol_have_ann="yes"
		fi
	fi

	if test "$ol_have_ann" != "yes" ; then
		if test "$ol_with_ann" = "yes" ; then
			AC_MSG_ERROR([ANN library headers not available])
		fi
	else
		if test "$ol_link_lapack" != "yes" ; then
			AC_MSG_ERROR([ANN library needs LAPACK])
		fi
		LIBS="$SAVE_LIBS"
	fi
	AC_LANG_C

	if test "$ol_have_ann" = "yes" ; then
		ol_link_ann="yes"
		AC_DEFINE(USE_ANN,1,[define to use ANN])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for SASL2
dnl
ol_link_sasl2="no"
if test "$ol_with_sasl2" != "no" ; then
	AC_CHECK_HEADERS(sasl.h sasl/sasl.h)

	have_sasl_h="no"
	if test "$ac_cv_header_sasl_h" = "yes" ; then
		have_sasl_h="yes"
	elif test "$ac_cv_header_sasl_sasl_h" = "yes" ; then
		have_sasl_h="yes"
	fi
	if test "$have_sasl_h" = "yes" ; then
		AC_CHECK_LIB(sasl2,sasl_client_start,[have_sasl2=yes],[have_sasl2=no])
		if test "$have_sasl2" != "no" ; then
			SECURITY_LIBS="$SECURITY_LIBS -lsasl2"
			AC_DEFINE(HAVE_SASL2,1,[define to use Cyrus SASL2])
			ol_link_sasl2="yes"
		fi
	fi

	if test "$ol_link_sasl2" = "no" ; then
		if test "$ol_with_sasl2" != "auto" ; then
			AC_MSG_ERROR([Could not locate SASL2])
		else
			AC_MSG_WARN([SASL2 not supported])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for PAM
dnl
ol_link_pam=no
if test $ol_with_pam != no ; then
	AC_CHECK_HEADERS(security/pam_appl.h pam_misc.h security/pam_misc.h)
	
	if test "$ac_cv_header_security_pam_appl_h" = "yes" ; then
		if test "$ac_cv_header_pam_misc_h" = "yes" \
				-o "$ac_cv_header_security_pam_misc_h" = "yes" ; then
			AC_CHECK_LIB(pam,pam_start,
				[have_pam=yes],[have_pam=no],
				[-lcrypt -lpam_misc -ldl])
			
			if test $have_pam != no ; then
				SECURITY_LIBS="$SECURITY_LIBS -lcrypt -lpam -lpam_misc -ldl"
				AC_DEFINE(USE_PAM,1,
					[define if you want to use PAM])
				ol_link_pam=yes
			fi
		fi
	fi

	if test $ol_link_pam = no ; then
		if test $ol_with_pam != auto ; then
			AC_MSG_ERROR([Could not locate PAM])
		else
			AC_MSG_WARN([PAM not supported])
		fi		
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check if crypt needs libcrypt
dnl
ol_link_crypt=no
if test "$ol_enable_crypt" != no ; then
	AC_MSG_CHECKING([how to link crypt])
	crypt_pool="\"\" -lcrypt -lcrypto"
	for i in $crypt_pool ; do
		save_LIBS=$LIBS
		LIBS="$LIBS $i"
		AC_TRY_LINK([
#include <unistd.h>
#include <pwd.h>
		],[
	crypt(NULL,NULL);
		],[crypt_ok=yes],[crypt_ok=no])
		LIBS=$save_LIBS

		if test $crypt_ok = yes ; then
			if test "x$i" = "x" ; then
				AC_MSG_RESULT([automatically linked])
			else
				AC_MSG_RESULT([with $i])
				SECURITY_LIBS="$SECURITY_LIBS $i"
			fi
			AC_DEFINE(HAVE_CRYPT,1,[define if crypt(3) can be linked])
			break
		fi
	done

	if test $crypt_ok != yes ; then
		AC_MSG_RESULT([no])
		if test "$ol_enable_crypt" = "yes" ; then
			AC_MSG_ERROR([don't know how to link crypt])
		else
			AC_MSG_WARN([don't know how to link crypt...])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for MPI & MPI++
dnl
dnl FIXME: needs to be reworked for MPI2
dnl
ol_link_mpi=no
if test $ol_with_mpi != no ; then
	AC_CHECK_HEADERS(mpi.h)
	AC_LANG_CPLUSPLUS
	AC_CHECK_HEADERS(mpi++.h)
	AC_CHECK_HEADERS(mpicxx.h)
	AC_LANG_C

	have_mpi_cxx_header=no
	if test $ac_cv_header_mpipp_h = yes ; then
		dnl mpich
		have_mpi_cxx_header=yes
		mpi_cxx_header="mpi++.h"
		if test "$ol_with_mpi" = "pmpi" ; then
        		try_mpi_LIBS="-lpmpich++ -lmpich -lpmpich -lnsl"
		else
        		try_mpi_LIBS="-lmpich++ -lmpich -lnsl"
		fi
	elif test $ac_cv_header_mpicxx_h = yes ; then
		dnl openmpi
		have_mpi_cxx_header=yes
		mpi_cxx_header="mpicxx.h"
		try_mpi_LIBS="-lmpi -lmpi_cxx"
	fi
	
	if test $ac_cv_header_mpi_h = yes -a $have_mpi_cxx_header = yes ; then
		AC_MSG_CHECKING([for mpi libraries])
		AC_LANG_CPLUSPLUS
		try_mpi_LIBS=
		save_LIBS=$LIBS
 		AC_TRY_LINK([
#include <$mpi_cxx_header>
		],[
 	int argc = 0;
 	char **argv = 0;
 	MPI::Init(argc, argv);
 		],[have_mpi=yes],[have_mpi=no])
                if test $have_mpi = no ; then
		        LIBS="$LIBS $try_mpi_LIBS"
		        AC_TRY_LINK([
#include <mpi_cxx_header>
			],[
	int argc = 0;
	char **argv = 0;
	MPI::Init(argc, argv);
        		],[have_mpi=yes],[have_mpi=no])
	        	LIBS=$save_LIBS
                fi
		AC_LANG_C

		if test $have_mpi != no ; then
			MPI_LIBS=$try_mpi_LIBS
			AC_DEFINE(USE_MPI,1,[define if you want to use MPI])
			ol_link_mpi=yes
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi
	fi

	if test $ol_link_mpi = no ; then	
		if test $ol_with_mpi != auto ; then
			AC_MSG_ERROR([Could not locate MPI])
		else
			AC_MSG_WARN([MPI not supported])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Need Schur parallel solver with MPI
dnl
if test $ol_enable_schur != no ; then
	if test $ol_link_mpi = yes ; then
		if test $ol_link_metis = yes -o $ol_link_chaco = yes ; then
			use_schur=yes
			AC_DEFINE(USE_SCHUR,1,[define if you want to use Schur parallel solver])
		else
			if test $ol_enable_schur = yes ; then
				AC_MSG_ERROR([Need either Metis or Chaco for Schur parallel solver])
			fi
		fi
	else
		if test $ol_enable_schur = yes ; then
			AC_MSG_ERROR([Need MPI for Schur parallel solver])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for TAUCS
dnl
ol_link_taucs=no
if test $ol_with_taucs != no ; then
	AC_MSG_CHECKING([for taucs])

	if test $ol_link_metis = no ; then
		have_taucs=no
		if test $ol_with_taucs != auto ; then
			AC_MSG_ERROR([Need Metis library for Taucs solver])
		else
                        AC_MSG_WARN([Taucs needs Metis; not supported])
		fi
	elif test $have_lapack = no ; then
		have_taucs=no
		if test $ol_with_taucs != auto ; then
			AC_MSG_ERROR([Need Lapack library for Taucs solver])
		else
                        AC_MSG_WARN([Taucs needs Lapack; not supported])
		fi
	elif test $have_blas = no ; then
		have_taucs=no
		if test $ol_with_taucs != auto ; then
			AC_MSG_ERROR([Need Blas library for Taucs solver])
		else
                        AC_MSG_WARN([Taucs needs BLAS; not supported])
		fi
	else
		save_LIBS="$LIBS"
		try_taucs_LIBS="-ltaucs"
		LIBS="$try_taucs_LIBS $METIS_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS"
		AC_TRY_LINK([#include "taucs.h"
			],[
			taucs_ccs_matrix* A;
			void**            F;
			int		  nrhs;
			void*		  X;
			void*		  B;
			char**		  options;
			void**		  opt_arg;
			int taucs_linsolve(A, F, nrhs, X, B, options, opt_arg);
		],[have_taucs=yes],[have_taucs=no])
		LIBS="$save_LIBS"
	fi
	if test $have_taucs != no ; then
		AC_DEFINE(USE_TAUCS,1,[define to use Taucs libraries])
		TAUCS_LIBS="$try_taucs_LIBS"
	fi
	if test "$have_taucs" = "yes" ; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])	
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for POSIX realtime
dnl
if test $ol_with_rt != no ; then
	save_LIBS="$LIBS"
	try_RT_LIBS="-lrt"
	LIBS="$try_RT_LIBS $LIBS"
	AC_CHECK_LIB(rt,clock_nanosleep,[
		have_librt=yes
		RT_LIBS="$try_RT_LIBS"
	],[have_librt=no])
	LIBS="$save_LIBS"

	if test "$have_librt" = "yes" ; then
		AC_DEFINE(USE_RT,1,[define to use POSIX realtime])
	elif test "$ol_with_rt" = "yes" ; then
		AC_MSG_ERROR([POSIX realtime is not available])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for RTAI
dnl
use_rtai=no
if test $ol_with_rtai != no ; then

	dnl RTAI >= 3.0
	AC_CHECK_HEADERS(rtai_lxrt.h)
	AC_CHECK_HEADERS(rtai_netrpc.h)

	if test "$ac_cv_header_rtai_lxrt_h" = "yes" -a \
			"$ac_cv_header_rtai_netrpc_h" = "yes" ; then
		use_rtai=yes
	else
		dnl RTAI ~ 2.4.13
		AC_CHECK_HEADERS(rtai_lxrt_user.h)
		AC_CHECK_HEADERS(net_rpc.h)
		if test "$ac_cv_header_rtai_lxrt_user_h" = "yes" -a \
				"$ac_cv_header_net_rpc_h" = "yes" ; then
			use_rtai=yes
		else
			rtai_mismatch=no
			if test "$ac_cv_header_rtai_lxrt_h" = "yes" -a \
					"$ac_cv_header_rtai_netrpc_h" = "no" ; then
				rtai_mismatch=yes
			elif test "$ac_cv_header_rtai_lxrt_user_h" = "yes" -a \
					"$ac_cv_header_net_rpc_h" = "no" ; then
				rtai_mismatch=yes
			fi

			if test "$rtai_mismatch" = "yes" ; then
				AC_MSG_WARN([possible RTAI version mismatch; YMMV])
			fi

		fi
	fi

	dnl FIXME: no longer used
	dnl FIXME: is this correct for >= 3.0 as well?
	dnl AC_CHECK_HEADERS(rtai_usp_posix.h)
	dnl AC_CHECK_HEADERS(rtai_posix.h)
	dnl if test "$ac_cv_header_rtai_usp_posix_h" = "no" -o \
	dnl 		"$ac_cv_header_rtai_posix_h" = "no" ; then
	dnl 	use_rtai=no
	dnl fi
	if test "$use_rtai" = "yes" ; then
		AC_DEFINE(USE_RTAI,1,[define to use RTAI])
	elif test "$ol_with_rtai" = "yes" ; then
		AC_MSG_ERROR([RTAI is not available])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks id DeadBeat can be built
dnl
if test $ol_link_lapack != no -o $ol_link_meschach != no ; then
	AC_DEFINE(USE_DBC,1,[define if you want to use GPC/DeadBeat control])
else
	dnl Required for GPC/Deadbeat (should be able to work it 'round...)
	AC_MSG_WARN([GPC/Deadbeat not allowed])
fi

dnl ----------------------------------------------------------------
dnl
dnl Use GiNaC?
dnl
ol_link_ginac=no
if test $ol_with_ginac != no ; then
	AC_LANG_CPLUSPLUS
	AC_LANG_C
	AC_MSG_CHECKING([for GiNaC library])
	have_ginac=no

	ginac-config --version >/dev/null 2>&1
	if test $? = 0 ; then
		GINACLIB_LIBS=`ginac-config --libs`
		GINACLIB_CPPFLAGS=`ginac-config --cppflags`
		have_ginac=yes
	else
		pkg-config ginac
		if test $? = 0; then
			GINACLIB_LIBS=`pkg-config --libs ginac`
			GINACLIB_CPPFLAGS=`pkg-config --cflags ginac`
			have_ginac=yes
		fi
	fi

	if test $have_ginac = yes ; then
		ol_link_ginac=yes
		AC_DEFINE(USE_GINAC,1,[define if GiNaC is available])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		if test $ol_with_ginac != auto; then
			AC_MSG_ERROR([GiNaC not available])
		else
			AC_MSG_WARN([GiNaC not available])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for library functions.
dnl
AC_CHECK_FUNCS([ \
	atoi \
	atol \
	chdir \
	copysign \
	finite \
	get_ncpus \
	get_nprocs \
	get_nprocs_conf \
	getaddrinfo \
	getcwd \
	gethostbyname \
	gethostbyname2 \
	getopt \
	getopt_long \
	getpid \
	inet_aton \
	localtime \
	memset \
	memmove \
	mkstemp \
	mlockall \
	nanosleep \
	pstat_getdynamic \
	pstat_getprocessor \
	putenv \
	setenv \
	shmget \
	signal \
	sleep \
	strftime \
	strtod \
	strtol \
	strtoul \
	time \
	times \
	unsetenv \
	usleep \
])

dnl ----------------------------------------------------------------
dnl
dnl Checks if we can detect how many processors can be used on smp arch
dnl
if test "$ac_cv_func_get_nprocs" = "no" ; then
	if test "$ac_cv_header_sys_sysinfo_h" = "yes" ; then
	        AC_TRY_LINK([
#include <sys/sysinfo.h>
       		],[
        return get_nprocs();
       		],[have_get_nprocs=yes],[have_get_nprocs=no])
	else
	        AC_TRY_LINK([
       		],[
        return get_nprocs();
       		],[have_get_nprocs=yes],[have_get_nprocs=no])
	fi
	
	if test "$have_get_nprocs" = "yes" ; then
		AC_DEFINE(HAVE_GET_NPROCS,1,[define if get_nprocs is available])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Make sure we can use sleep with C++
dnl
if test "$ac_cv_func_sleep" = "yes" ; then
	AC_LANG_CPLUSPLUS
	AC_TRY_LINK([
#include <unistd.h>
       	],[
        sleep(1);
       	],[use_sleep=yes],[use_sleep=no])
	if test "$use_sleep" = "yes" ; then
		AC_DEFINE(USE_SLEEP,1,[define if sleep is usable with C++ (issue with MinGW)])
	fi
	AC_LANG_C
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks whether ios has fmtflags type
dnl

AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether class ios has fmtflags type])
AC_TRY_COMPILE([#include <iostream>],[
	std::ios::fmtflags f;
],[have_fmtflags=yes],[have_fmtflags=no])

if test $have_fmtflags != no ; then
	AC_DEFINE(HAVE_FMTFLAGS_IN_IOS,1,
		[define if class ios has 'fmtflags' type])
	AC_MSG_RESULT([yes])
else
	AC_TRY_COMPILE([#include <iostream.h>],[
		std::ios::fmtflags f;
],[have_fmtflags=yes],[have_fmtflags=no])

	if test $have_fmtflags != no ; then
		AC_DEFINE(HAVE_FMTFLAGS_IN_IOS,1,
			[define if class ios has 'fmtflags' type])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi
fi
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Variable substituion
dnl
AC_SUBST(MPI_LIBS)
AC_SUBST(CHACO_LIBS)
AC_SUBST(METIS_LIBS)
AC_SUBST(UMFPACK_LIBS)
AC_SUBST(KLU_LIBS)
AC_SUBST(Y12_LIBS)
AC_SUBST(SUPERLU_LIBS)
AC_SUBST(WSMP_LIBS)
AC_SUBST(HARWELL_LIBS)
AC_SUBST(MESCHACH_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(TAUCS_LIBS)
AC_SUBST(ARPACK_LIBS)
AC_SUBST(JDQZ_LIBS)
AC_SUBST(THREAD_LIBS)
AC_SUBST(ATOMIC_OPS_LIBS)
AC_SUBST(SECURITY_LIBS)
AC_SUBST(NETCDF_LIBS)
AC_SUBST(UDUNITS_LIBS)
AC_SUBST(RT_LIBS)
AC_SUBST(ANN_LIBS)
AC_SUBST(BULLET_LIBS)
AC_SUBST(PYTHON_NUMPY_HEADERS)
AC_SUBST(OCTAVE_LIBS)

AC_SUBST(MODULES_DIRS)
AC_SUBST(MODULE_LOADPATH)

dnl AC_SUBST(GINACLIB_CONFIG)
AC_SUBST(GINACLIB_CPPFLAGS)
AC_SUBST(GINACLIB_LIBS)

AC_SUBST(DATE)

dnl conditional build declarations
AM_CONDITIONAL(USE_ARPACK,[test "$ol_link_arpack" = yes])
AM_CONDITIONAL(USE_CHACO,[test "$ol_link_chaco" = yes])
AM_CONDITIONAL(USE_GINAC,[test "$ol_link_ginac" = yes])
AM_CONDITIONAL(USE_JDQZ,[test "$ol_link_jdqz" = yes])
AM_CONDITIONAL(USE_METIS,[test "$ol_link_metis" = yes])
AM_CONDITIONAL(USE_OCTAVE,[test "$use_octave" = yes])
AM_CONDITIONAL(USE_PYTHON,[test "$use_python" = yes])
AM_CONDITIONAL(USE_RTAI,[test "$use_rtai" = yes])
AM_CONDITIONAL(USE_SCHUR,[test "$use_schur" = yes])
AM_CONDITIONAL(USE_SOCKET,[test "$use_socket" = yes])
AM_CONDITIONAL(USE_ANN,[test "$ol_link_ann" = yes])

AM_CONDITIONAL(FOOL_LIBTOOL,[test "$fool_libtool" = yes])
AM_CONDITIONAL(BUILD_Y12,[test "$ol_link_y12" = yes])

AM_CONDITIONAL(BUILD_STATIC_MODULES,[test "$build_static_modules" = yes])
AM_CONDITIONAL(BUILD_CHARM,[test "$build_charm" = yes])

dnl ----------------------------------------------------------------
dnl
dnl Prepare for libtool output:
dnl

dnl LTLIBOBJS=`echo "$LIBOBJS" | sed 's/\.o/\.lo/g'`
dnl AC_SUBST(LTLIBOBJS)
dnl LTALLOCA=`echo "$ALLOCA" | sed 's/\.o/\.lo/g'`
dnl AC_SUBST(LTALLOCA)

dnl ----------------------------------------------------------------
dnl
dnl final output
dnl
AC_OUTPUT([ \
	Makefile:Makefile.in
	./libraries/Makefile:libraries/Makefile.in
	./libraries/libmbutil/Makefile:libraries/libmbutil/Makefile.in
	./libraries/libobjs/Makefile:libraries/libobjs/Makefile.in
	./libraries/libcolamd/Makefile:libraries/libcolamd/Makefile.in
	./libraries/libnaive/Makefile:libraries/libnaive/Makefile.in
	./libraries/liby12/Makefile:libraries/liby12/Makefile.in
	./libraries/libmbmath/Makefile:libraries/libmbmath/Makefile.in
	./libraries/libmbwrap/Makefile:libraries/libmbwrap/Makefile.in
	./libraries/libann/Makefile:libraries/libann/Makefile.in
	./libraries/libmbc/Makefile:libraries/libmbc/Makefile.in
	./mbdyn/Makefile:mbdyn/Makefile.in
	./mbdyn/base/Makefile:mbdyn/base/Makefile.in
	./mbdyn/aero/Makefile:mbdyn/aero/Makefile.in
	./mbdyn/elec/Makefile:mbdyn/elec/Makefile.in
	./mbdyn/thermo/Makefile:mbdyn/thermo/Makefile.in
	./mbdyn/hydr/Makefile:mbdyn/hydr/Makefile.in
	./mbdyn/struct/Makefile:mbdyn/struct/Makefile.in
	./modules/Makefile:modules/Makefile.in
	./utils/Makefile:utils/Makefile.in
	./etc/Makefile:etc/Makefile.in
	./etc/modal.d/Makefile:etc/modal.d/Makefile.in
	./etc/pam.d/Makefile:etc/pam.d/Makefile.in
	./var/Makefile:var/Makefile.in
	./manual/Makefile:manual/Makefile.in
	./manual/mbdyn.1:manual/mbdyn.1.in
	./include/Makefile:include/Makefile.in
	./include/ac/Makefile:include/ac/Makefile.in
	./build/Makefile:build/Makefile.in
],[
	date > stamp-h
])

