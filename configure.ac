dnl $Header$
dnl Process this file with autoconf to produce a configure script.
dnl
dnl MBDyn (C) is a multibody analysis code.
dnl http://www.mbdyn.org
dnl
dnl Copyright (C) 1996-2021
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
# Copyright (C) 1996-2017
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
echo "	Copyright (C) 1996-2017"
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
dnl configure.ac for MBDyn

DATE=`date "+%B %e, %Y"`

# set unset (borrowed from autoconf 2.14a)
if (unset FOO) >/dev/null 2>&1; then
	ol_unset=unset
else
	ol_unset=false
fi
# unset CDPATH
$ol_unset CDPATH || test "${CDPATH:+set}" != set || CDPATH=: && export CDPATH

AC_INIT(mbdyn, m4_esyscmd([tr -d '\n' < build/version]))
AC_CONFIG_AUX_DIR(build)dnl
AC_CONFIG_SRCDIR([build/version])

MBDYN_DEVEL=no
if test -f "$ac_aux_dir/DEVEL" ; then
	MBDYN_DEVEL=yes
	AC_DEFINE(MBDYN_DEVEL,1,[define to enable features under development])
	MB_VERSION="Devel"
else
	MB_VERSION=`cat $ac_aux_dir/version`
fi

if test -z "$MB_VERSION"; then
	AC_MSG_ERROR([could not determine version])
fi
echo "Configuring MBDyn $MB_VERSION ..."

dnl Determine host platform
dnl             we try not to use this for much
AC_CANONICAL_SYSTEM

AM_INIT_AUTOMAKE

dnl We use autoconf features new to 2.13.
dnl aclocal.m4 should be built using aclocal from automake 1.4
dnl libtool 1.3.3 should be installed.
AC_PREREQ(2.13)dnl Required Autoconf version

AC_CONFIG_HEADERS([include/mbconfig.h])
AH_TOP([
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

PKG_PROG_PKG_CONFIG

dnl ----------------------------------------------------------------
dnl General "enable" options
OL_ARG_ENABLE(Werror,[  --enable-Werror          enable building with -Werror], yes)dnl
OL_ARG_ENABLE(debug,[  --enable-debug          enable debugging], no)dnl
OL_ARG_WITH(debug_mode,[  --with-debug-mode[={none|mem}] with debug mode {none|mem}],none,[none mem])
OL_ARG_ENABLE(runtime_loading,[  --enable-runtime-loading    enable runtime loading], auto)dnl
OL_ARG_WITH(static_modules,[  --with-static-modules   build (known) modules as static],auto,[no yes auto])
OL_ARG_ENABLE(crypt,[  --enable-crypt          enable crypt (deprecated)], no)dnl
OL_ARG_ENABLE(eig,[  --enable-eig            enable direct eigensolution
                          (needs one of LAPACK, ARPACK or JDQZ)],auto,[auto yes no])dnl needs LAPACK, ARPACK or JDQZ
OL_ARG_ENABLE(schur,[  --enable-schur          enable Schur parallel solver
                          (needs MPI and either Metis or Chaco)],auto,[auto yes no])dnl needs MPI and Metis or Chaco
OL_ARG_ENABLE(multithread,[  --enable-multithread    enable multithread assembly],no,[auto yes no force])dnl
OL_ARG_ENABLE(multithread_naive,[  --enable-multithread-naive    enable multithread naive solver],no,[auto yes no force])dnl
OL_ARG_ENABLE(mbc,[  --enable-mbc            enable MBC - multibody communication library],yes)dnl
OL_ARG_ENABLE(netcdf,[  --enable-netcdf         enable NetCDF based binary output],auto,[auto yes no legacy netcdf4])dnl
OL_ARG_ENABLE(python,[  --enable-python         enable Python support],no)dnl
OL_ARG_ENABLE(octave,[  --enable-octave         enable Octave support],no)dnl
OL_ARG_ENABLE(octave_utils,[  --enable-octave-utils   enable Octave utility functions],no)dnl
AC_ARG_WITH(mkoctfile,[  --with-mkoctfile=exec-path Octave build tool mkoctfile],[
	ol_with_mkoctfile_exec="$withval"
],[ol_with_mkoctfile_exec=mkoctfile])
OL_ARG_ENABLE(autodiff,[  --enable-autodiff       enable automatic differentiation in order to compute the Jacobian matrix],no)dnl
OL_ARG_ENABLE(autodiff_vect, [  --enable-autodiff_vect={0|16|32|64}     enable vectorization for automatic differentiation],0,[0 16 32 64])dnl
OL_ARG_ENABLE(sparse_autodiff, [  --enable-sparse_autodiff                enable automatic differentiation based on sparse algorithms],no)dnl
OL_ARG_ENABLE(install_test_progs, [  --enable-install_test_progs       enable installation of some test programs which are otherwise not installed], no)
dnl ----------------------------------------------------------------
dnl General "with" options
OL_ARG_WITH(g2c,[  --with-g2c              with g2c header],no,[auto yes no])
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
  --with-blas[={auto|goto|atlas|cblas|blas|openblas}]
                          with (C)BLAS math library],auto,[auto goto atlas cblas blas openblas])
AC_ARG_WITH(goto,[  --with-goto=lib(s)      with Goto BLAS implementation],[
	ol_with_goto_lib="$withval"
],[ol_with_goto_lib=])
AC_ARG_WITH(openblas,[  --with-openblas=lib(s)  with OpenBLAS implementation],[
	ol_with_openblas_lib="$withval"
],[ol_with_openblas_lib=])
OL_ARG_WITH(ginac,[  --with-ginac            with GiNaC support
                          (ginac-config must be in $PATH)],auto,[auto yes no])
OL_ARG_WITH(y12,[
       optional linear algebra solvers (naive is enabled by default):
  --with-y12              with Y12 sparse math library],auto,[auto yes no])
OL_ARG_WITH(umfpack,[  --with-umfpack          with Umfpack math library],auto,[auto yes no])
OL_ARG_WITH(klu,[  --with-klu              with KLU math library],auto,[auto yes no])
OL_ARG_WITH(lapack,[  --with-lapack           with LAPACK math library],auto,[auto yes no])
OL_ARG_WITH(harwell,[  --with-hsl              with HSL (Harwell) sparse math library - historical],no,[auto yes no])
OL_ARG_WITH(meschach,[  --with-meschach         with Meschach math library - historical],no,[auto yes no])
OL_ARG_WITH(superlu,[  --with-superlu          with SuperLU math library - eXperimental],no,[auto yes no])
OL_ARG_WITH(wsmp,[  --with-wsmp             with IBM Watson library],no,[auto yes no])
OL_ARG_WITH(pardiso,[  --with-pardiso          with Intel MKL Pardiso library],no,[auto yes no])
OL_ARG_WITH(pastix,[  --with-pastix           with Pastix math library],no,[auto yes no])
with_suitesparseqr_default=auto
OL_ARG_WITH(qrupdate,[  --with-qrupdate         with QRUPDATE library],$with_suitesparseqr_default,[auto yes no])
OL_ARG_WITH(suitesparseqr,[  --with-suitesparseqr    with SuiteSparseQR library],$with_suitesparseqr_default,[auto yes no])
OL_ARG_WITH(strumpack,[  --with-strumpack        with STRUMPACK library],auto,[auto yes no])  
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
build_windows=no
case "$host" in
	*-apple-*)
		fool_libtool=no
		;;
    cygwin*|*mingw*)
        build_windows=yes
        ;;
	*)
		;;
esac

dnl Add a conditional we can use in Makefile.am for MS Windows builds
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])

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
AC_PROG_F77
have_fc_f77=no
if test "x$FC" = "x" ; then
	if test "x$F77" = "x" ; then
		AC_MSG_ERROR([Need a working Fortran compiler])
	fi
	AC_F77_LIBRARY_LDFLAGS
	LIBS="$LIBS $F77LIBS"
	have_fc_f77=yes
else
	F77="$FC"
	F77LIBS="$FCLIBS"
	LIBS="$LIBS $FCLIBS"
	have_fc_f77=yes
fi

if test "$have_fc_f77" = "yes" ; then
	AC_DEFINE(HAVE_FC_F77, 1, [define if a Fortran 77 compiler is available])
	AC_DEFINE(USE_AEROD2_F, 1, [define Fortran 77 aerodynamic routines can be used])
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
PARDISO_LIBS=
WSMP_LIBS=
HARWELL_LIBS=
LAPACK_LIBS=
TAUCS_LIBS=
PASTIX_LIBS=
PASTIX_INC=
QRUPDATE_LIBS=
SUITESPARSEQR_LIBS=
STRUMPACK_LIBS=
ARPACK_LIBS=
JDQZ_LIBS=
SECURITY_LIBS=
NETCDF_LIBS=
UDUNITS_LIBS=
RTAI_LDFLAGS=
RT_LIBS=
ANN_LIBS=
BULLET_LIBS=
PYTHON_NUMPY_HEADERS=
OCTAVE_LIBS=
OCTAVE_INCLUDE=
OCTAVE_MKOCTFILE=
dnl GINACLIB_CONFIG=
GINACLIB_LIBS=
GINACLIB_CPPFLAGS=
wl=
EXPORT_DYNAMIC_FLAG_SPEC=

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
AC_CHECK_HEADERS(sched.h)
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
AC_CHECK_HEADERS(sys/ioctl.h)
AC_CHECK_HEADERS(time.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(values.h)
dnl termios.h used in utils autopilot.c, we don't build if we don't have it (e.g. on windows)
AC_CHECK_HEADER(termios.h, [have_termios=yes], [have_termios=no])

dnl on windows we use winsock
AC_CHECK_HEADERS(winsock2.h)


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

have_environ=no
AC_TRY_COMPILE([
],[
	environ;
], [have_environ=yes], [have_environ=no])
if test $have_environ != yes ; then
	AC_TRY_COMPILE([
#include <unistd.h>
	],[
	environ;
	], [have_environ=yes], [have_environ=no])
	if test $have_environ = yes ; then
		AC_DEFINE(ENVIRON_NEEDS_UNISTD_H,1,[define if variable environ needs unistd.h])
	fi
fi
if test $have_environ != yes ; then
	AC_TRY_COMPILE([
#define _GNU_SOURCE
#include <unistd.h>
	],[
	environ;
	], [have_environ=yes], [have_environ=no])
	if test $have_environ = yes ; then
		AC_DEFINE(ENVIRON_NEEDS_GNU_SOURCE_AND_UNISTD_H,1,[define if variable environ needs _GNU_SOURCE and unistd.h])
	fi
fi
if test $have_environ != yes ; then
	AC_TRY_LINK([
extern char **environ;
	],[
	environ;
	],[have_environ=yes],[have_environ=no])
	if test $have_environ = yes ; then
		AC_DEFINE(ENVIRON_NEEDS_DECLARATION,1,[define if variable environ needs to be declared])
	fi
fi
if test $have_environ = yes; then
	AC_DEFINE(HAVE_ENVIRON,1,[define if variable environ is declared])
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
dnl Enable NetCDF binary output (Classic/legacy or netcdf4)
dnl

if test "$ol_enable_netcdf" != no ; then
	have_netcdf4=no
	have_libnetcdf4=no
	have_netcdfc=no
	have_libnetcdfc=no
	AC_LANG_CPLUSPLUS
	if test "$ol_enable_netcdf" != legacy ; then
dnl checks for netcdf header file (for netcdf4 c++):
		AC_CHECK_HEADER(netcdf,[have_netcdf4=yes],[have_netcdf4=no])
		AC_LANG_C
		if test $have_netcdf4 = "yes" ; then
			AC_DEFINE(USE_NETCDF4,1,[must include <netcdf> to use NetCDF4 C++ support])
			save_LIBS="$LIBS"
			try_NETCDF_LIBS="-lnetcdf_c++4 -lnetcdf"
			LIBS="$try_NETCDF_LIBS $LIBS"
			AC_LANG_CPLUSPLUS
			AC_TRY_LINK([
	#include <netcdf>
				],[
		/* test for netcdf */
		netCDF::NcFile nc("foobar", netCDF::NcFile::replace); // Create and leave in define mode
				],[
					have_libnetcdf4=yes
					NETCDF_LIBS="$try_NETCDF_LIBS"
				],[
				  try_NETCDF_LIBS="-lnetcdf-cxx4 -lnetcdf"
			      LIBS="$try_NETCDF_LIBS $save_LIBS"
				  AC_TRY_LINK([
	#include <netcdf>
				  ],[
		/* test for netcdf */
		netCDF::NcFile nc("foobar", netCDF::NcFile::replace); // Create and leave in define mode
				  ],[
					have_libnetcdf4=yes
					NETCDF_LIBS="$try_NETCDF_LIBS"
				  ])
				])
			AC_TRY_LINK([
	#include <netcdf>
				],[
		/* test for netcdf sync */
		netCDF::NcFile nc("foobar", netCDF::NcFile::replace); // Create and leave in define mode
		nc.sync(); // test for sync option (comes in later version, >=4.3.0)
		nc.close(); // test for close (appears in 4.3.0)
				],[

				], [
					have_libnetcdf4=no
					AC_MSG_WARN([netcdf-cxx4 version too old, either link version 4.3.0 or more recent
					 or otherwise use the legacy version with the --enable-netcdf command
					 falling back to checking for legacy netcdf version ])
				])
			AC_LANG_C
			LIBS="$save_LIBS"
		fi
	fi
	if test "$ol_enable_netcdf" != netcdf4 ; then
		if test $have_libnetcdf4 = no ; then
			AC_LANG_CPLUSPLUS
			AC_CHECK_HEADER(netcdfcpp.h,[have_netcdfcpp_h=yes],[have_netcdfcpp_h=no])
			AC_LANG_C
			if test $have_netcdfcpp_h = "yes" ; then
				AC_DEFINE(HAVE_NETCDFCPP_H,1,[must include <netcdfcpp.h> to use NetCDF C++ support])
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
						have_libnetcdfc=yes
						NETCDF_LIBS="$try_NETCDF_LIBS"
					])
				AC_LANG_C
				LIBS="$save_LIBS"
			else
				dnl apparently, NetCDF 4 only needs header <netcdf.h>
				dnl and thus does not provide <netcdfcpp.h> any longer
				AC_LANG_CPLUSPLUS
				AC_CHECK_HEADER(netcdf.h,[have_netcdf_h=yes],[have_netcdf_h=no])
				AC_LANG_C
				if test $have_netcdf_h = "yes" ; then
					AC_DEFINE(HAVE_NETCDF_H,1,[<netcdf.h> suffices to use NetCDF C++ support])
					save_LIBS="$LIBS"
					try_NETCDF_LIBS="-lnetcdf_c++4 -lnetcdf"
					LIBS="$try_NETCDF_LIBS $LIBS"
					AC_LANG_CPLUSPLUS
					AC_TRY_LINK([
		#include <netcdf.h>
						],[
				/* test for netcdf */
				NcFile nc("foobar", NcFile::Replace); // Create and leave in define mode
						],[
							have_libnetcdfc=yes
							NETCDF_LIBS="$try_NETCDF_LIBS"
						],[
						try_NETCDF_LIBS="-lnetcdf-cxx4 -lnetcdf"
					    LIBS="$try_NETCDF_LIBS $save_LIBS"
                        AC_TRY_LINK([
		#include <netcdf.h>
						  ],[
				/* test for netcdf */
				NcFile nc("foobar", NcFile::Replace); // Create and leave in define mode
						  ],[
							have_libnetcdfc=yes
							NETCDF_LIBS="$try_NETCDF_LIBS"
						  ],[
							have_libnetcdfc=no
						  ])
						])
					AC_LANG_C
					LIBS="$save_LIBS"
				fi
			fi
		fi
	fi
	if test $have_libnetcdfc = yes ; then
		AC_DEFINE(USE_NETCDF,1,[define this to enable NetCDF based output])
		AC_DEFINE(USE_NETCDFC,1,[define this to enable classic/legacy NetCDF based output])
		AC_MSG_NOTICE([Enabling netcdf classic/legacy output capabilities])
	elif test $ol_enable_netcdf != auto ; then
		if test "$ol_enable_netcdf" = legacy ; then
			AC_MSG_ERROR([User requested netcdf classic/legacy but cannot be enabled, see messages above.])
		else
			AC_MSG_WARN([Not using classic/legacy netcdf])
		fi
	fi

	if test $have_libnetcdf4 = yes ; then
		AC_DEFINE(USE_NETCDF4,1,[define this to enable NetCDF4 based output])
		AC_DEFINE(USE_NETCDF,1,[also define this to enable NetCDF4 based output])
		AC_MSG_NOTICE([Enabling netcdf4 output capabilities])
	elif test $ol_enable_netcdf != auto ; then
		if test "$ol_enable_netcdf" = netcdf4 ; then
			AC_MSG_ERROR([User requested netcdf4 but cannot be enabled, see messages above.])
		else
			AC_MSG_WARN([Not using netcdf4])
		fi
	fi

	if test "$have_libnetcdf4" != yes ; then
		if test "$have_libnetcdfc" != yes ; then
			if test $ol_enable_netcdf != auto ; then
				AC_MSG_ERROR([Could not locate neither netcdf4 nor netcdf classic/legacy headers or improper installation or versions installed.])
			fi
		fi
	fi

fi

if test "$ol_enable_netcdf" != no ; then
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
	        PYTHON_NUMPY_HEADERS=`$PYTHON -c "from __future__ import print_function; import numpy; print(\"{}\".format(numpy.get_include()))" 2>/dev/null`
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
use_octave_utils=no
if test "$ol_enable_octave" != no || test "$ol_enable_octave_utils" != no ; then
	have_mkoctfile=no
      	AC_MSG_CHECKING(for mkoctfile)
   	$ol_with_mkoctfile_exec --help > /dev/null 2>&1
	RC=$?
	if test "$RC" = 0 ; then
   	   	AC_MSG_RESULT([yes])
		have_mkoctfile=yes
	   	OCTAVE_MKOCTFILE=$ol_with_mkoctfile_exec
		dnl On Windows systems the output from mkoctfile contains backslashes which cause build failures!
	   	try_OCTAVE_LIBS="`$ol_with_mkoctfile_exec -p OCTAVE_LIBS` `$ol_with_mkoctfile_exec -p LDFLAGS | sed 's/-IC\:/-I\/c/g;s/-LC\:/-L\/c/g;s/\\\\/\//g'` -L`$ol_with_mkoctfile_exec -p OCTLIBDIR | sed 's/-IC\:/-I\/c/g;s/-LC\:/-L\/c/g;s/\\\\/\//g'`"
	   	try_OCTAVE_INCLUDE="`$ol_with_mkoctfile_exec -p INCFLAGS | sed 's/-IC\:/-I\/c/g;s/-LC\:/-L\/c/g;s/\\\\/\//g'`"
       	   	export LD_LIBRARY_PATH="`$ol_with_mkoctfile_exec -p OCTLIBDIR`:${LD_LIBRARY_PATH}"
    else
		AC_MSG_RESULT([no])
    	try_OCTAVE_LIBS="-loctinterp -loctave -lcruft"
    	try_OCTAVE_INCLUDE=
	fi

	AC_LANG_CPLUSPLUS
	save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $try_OCTAVE_INCLUDE"
	AC_CHECK_HEADER(octave/oct.h,[have_octave_oct_h=yes],[have_octave_oct_h=no])
	AC_CHECK_HEADER(octave/parse.h,[have_octave_parse_h=yes],[have_octave_parse_h=no],[
#include <octave/oct.h>
		])
	AC_CHECK_HEADER(octave/toplev.h,[have_octave_toplev_h=yes],[have_octave_toplev_h=no],[
#include <octave/oct.h>
		])
	AC_CHECK_HEADER(octave/interpreter.h,[have_octave_interpreter_h=yes],[have_octave_interpreter_h=no],[
#include <octave/oct.h>
		])
	AC_CHECK_HEADER(octave/octave.h,[have_octave_octave_h=yes],[have_octave_octave_h=no],[
#include <octave/oct.h>
		])
	CPPFLAGS="$save_CPPFLAGS"
	AC_LANG_C

	link_octave="no"
        have_clean_up_and_exit="no"
        have_do_octave_atexit="no"
        
	if test "$have_octave_oct_h" = "yes" \
		-a "$have_octave_parse_h" = "yes" \
		-a "$have_octave_toplev_h" = "yes" \
		-a "$have_octave_octave_h" = "yes" \
		-a "$have_mkoctfile" = "yes" ; \
	then
		OCTAVE_INCLUDE="$try_OCTAVE_INCLUDE"
		save_CPPFLAGS="$CPPFLAGS"
		save_LIBS="$LIBS"
		try_LIBS="$try_OCTAVE_LIBS"
		LIBS="$try_LIBS $LIBS"
		CPPFLAGS="$CPPFLAGS $OCTAVE_INCLUDE"
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>
			],[
		octave_main(0, NULL, 1);
			],[link_octave="yes"])
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>
			],[
		do_octave_atexit();
			],[have_do_octave_atexit="yes"])
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>
			],[
		clean_up_and_exit(0, true);
			],[have_clean_up_and_exit="yes"])
		AC_LANG_C
		LIBS="$save_LIBS"
		CPPFLAGS="$save_CPPFLAGS"
		OCTAVE_LIBS="$try_LIBS"
	fi

	if test "$have_octave_oct_h" = "yes" \
		-a "$have_octave_parse_h" = "yes" \
		-a "$have_octave_interpreter_h" = "yes" \
		-a "$have_octave_octave_h" = "yes" \
		-a "$have_mkoctfile" = "yes" ; \
	then
		OCTAVE_INCLUDE="$try_OCTAVE_INCLUDE"
		save_CPPFLAGS="$CPPFLAGS"
		save_LIBS="$LIBS"
		try_LIBS="$try_OCTAVE_LIBS"
		LIBS="$try_LIBS $LIBS"
		CPPFLAGS="$CPPFLAGS $OCTAVE_INCLUDE"
		AC_LANG_CPLUSPLUS
		AC_TRY_LINK([
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
#include <octave/octave.h>
			],[
                        octave::interpreter interpreter;
                        int status = interpreter.execute();
                        octave_value_list out = octave::feval ("cos", octave_value(0.), 1);
			],[link_octave="yes"])
		AC_LANG_C
		LIBS="$save_LIBS"
		CPPFLAGS="$save_CPPFLAGS"
		OCTAVE_LIBS="$try_LIBS"                        
	fi

	if test "$link_octave" = yes ; then
		if test "$have_do_octave_atexit" = yes; then
			AC_DEFINE(HAVE_DO_OCTAVE_ATEXIT, 1, [define if do_octave_atexit can be used])
		fi
		if test "$have_clean_up_and_exit" = yes; then
			AC_DEFINE(HAVE_CLEAN_UP_AND_EXIT, 1, [define if clean_up_and_exit can be used])
		fi
		use_octave=yes
		AC_DEFINE(USE_OCTAVE, 1, [define if octave can be used])
	else
		if test "$ol_enable_octave" = "yes" ; then
			AC_MSG_ERROR([octave support not available])
		fi
	fi

	if test "$have_octave_oct_h" = "yes" -a "$have_mkoctfile" = "yes"; then
	        use_octave_utils="yes"
	else
	        AC_MSG_ERROR([octave utils support not available])
	fi
fi

if test "$ol_enable_autodiff" != no || test "$ol_enable_sparse_autodiff" != no; then
	AC_MSG_CHECKING([for _aligned_malloc and _aligned_free])
	AC_TRY_LINK([
#include <stdlib.h>
				],[
void* p = _aligned_malloc(1024, 32);
_aligned_free(p);
				],[have_aligned_malloc="yes"])

	if test "$have_aligned_malloc" == "yes"; then
		AC_DEFINE(HAVE_ALIGNED_MALLOC, 1, [define if _aligned_malloc and _aligned_free can be used])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi

	AC_MSG_CHECKING([for posix_memalign])
	AC_TRY_LINK([
#include <stdlib.h>
				],[
void *p;
posix_memalign(&p, 32, 1024);
			],[have_posix_memalign="yes"])

	if test "$have_posix_memalign" == "yes"; then
		AC_DEFINE(HAVE_POSIX_MEMALIGN, 1, [define if posix_memalign can be used])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi

	AC_MSG_CHECKING([for memalign])
	AC_TRY_LINK([
#include <stdlib.h>
				],[
void* p = memalign(32, 1024);
				],[have_memalign="yes"])

	if test "$have_memalign" == "yes"; then
		AC_DEFINE(HAVE_MEMALIGN, 1, [define if memalign can be used])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi

	AC_MSG_CHECKING([for aligned_alloc])
	AC_TRY_LINK([
#include <stdlib.h>
				],[
void* p = aligned_alloc(32, 1024);
free(p);
				],[have_aligned_alloc="yes"])

	if test "$have_aligned_alloc" == "yes"; then
		AC_DEFINE(HAVE_ALIGNED_ALLOC, 1, [define if aligned_alloc can be used])
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi
fi

if test "$ol_enable_autodiff" != no ; then
	AC_DEFINE(USE_AUTODIFF, 1, [define if automatic differentiation can be used in order to compute the Jacobian matrix])

        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_CHECK_HEADERS(dense_hash_map)
        if test "$ac_cv_header_dense_hash_map" = "yes" ; then
                AC_MSG_CHECKING([for dense_hash_map])
                AC_TRY_COMPILE([#include <dense_hash_map>
                                #include <functional>],[
                                google::dense_hash_map<int,
                                int,
                                std::hash<int>,
                                std::equal_to<int>,
                                std::allocator<std::pair<const int, int> > > dhm;
                                dhm.set_empty_key(-1);
                                dhm.set_deleted_key(-2);
                                dhm[1] = 1;
                                for (auto i = dhm.begin(); i != dhm.end(); ++i) {
                                i->first;
                                i->second;
                                }
                                auto i = dhm.find(1);
                                i->first;
                                i->second;
                                dhm.clear();
                ], [have_dense_hash_map=yes], [have_dense_hash_map=no])
                if test $have_dense_hash_map = yes ; then
                        AC_DEFINE(USE_DENSE_HASH_MAP, 1, [define if google::dense_hash_map can be used])
                        AC_MSG_RESULT([yes])
                else
                        AC_MSG_RESULT([no])
                fi
        fi
        AC_LANG_RESTORE

	if test "$ol_enable_autodiff_vect" == 0 ; then
		save_CXXFLAGS="$CXXFLAGS"
		VECTOR_FLAGS="-march=native -mtune=native"
		CXXFLAGS="$CXXFLAGS $VECTOR_FLAGS"

		AC_LANG_CPLUSPLUS
		AC_TRY_RUN([
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <iomanip>
#include <fstream>
using namespace std;

template <typename scalar_type, int VECTOR_SIZE>
class vector
{
public:
    static const int vec_size = VECTOR_SIZE / sizeof(scalar_type);
    typedef scalar_type __attribute__((vector_size(VECTOR_SIZE))) vector_type;

    vector(int n)
    {
        N = n / vec_size;

        if (n % vec_size != 0)
        {
            ++N;
        }

        void* p;

#if HAVE_POSIX_MEMALIGN == 1
    	if (0 != posix_memalign(&p, sizeof(vector_type), vector_size() * sizeof(vector_type)))
    	{
    		p = 0;
    	}
#elif HAVE_ALIGNED_MALLOC == 1
        p = _aligned_malloc(vector_size() * sizeof(vector_type), sizeof(vector_type));
#elif HAVE_MEMALIGN == 1
    	p = memalign(sizeof(vector_type), vector_size() * sizeof(vector_type));
#else
		p = malloc(vector_size() * sizeof(vector_type));
#endif

        v = reinterpret_cast<vector_type*>(p);

        if (v == 0)
            throw std::bad_alloc();
    }

    ~vector()
    {
#if HAVE_ALIGNED_MALLOC == 1
        _aligned_free(v);
#else
    	free(v);
#endif
    }

    int size() const { return N * vec_size; }

    int vector_size() const { return N; }

    vector_type& operator[](int i)
    {
        return v[i];
    }

    vector_type operator[](int i) const
    {
        return v[i];
    }

    scalar_type& operator()(int i)
    {
        return s[i];
    }

    scalar_type operator()(int i) const
    {
        return s[i];
    }
private:
    union {
        scalar_type* s;
        vector_type* v;
    };
    int N;
};

template <typename T, int VECTOR_SIZE>
T func_vec(const T a, const T b, const vector<T,VECTOR_SIZE>& ad, const vector<T,VECTOR_SIZE>& bd,  vector<T,VECTOR_SIZE>& cd)
{
    for (int i = 0; i < ad.vector_size(); ++i)
    {
        cd[i] = ad[i] * b + a * bd[i];
    }

    return a * b;
}

template<typename T, int VECTOR_SIZE>
clock_t testVector(int N, int M)
{
    T a, b, c;

    a = 10;
    b = 5;
    vector<T,VECTOR_SIZE> ad(N), bd(N), cd(N);

    for (int i = 0; i < ad.size(); ++i)
    {
        ad(i) = i + 1;
        bd(i) = 10 * (i + 1);
    }

    const clock_t start = clock();

    for (int i = 0; i < M; ++i)
        c = func_vec(a, b, ad, bd, cd);

    const clock_t dt = clock() - start;

    assert(c == a * b);
    for (int i = 0; i < ad.size(); ++i)
    {
        assert(cd(i) == ad(i) * b + bd(i) * a);
    }

    return dt;
}

int vec_size = 0;
clock_t tmin = std::numeric_limits<clock_t>::max();

template <int i>
void runTest(int N, int M, clock_t t[])
{
    const int iprev = i - 1;
    runTest<iprev>(N, M, t);
    const int K = 1 << (i + 3);

    cout << "checking for vector_size=" << setw(3) << K << " ... ";

    t[i] = testVector<double, K>(N, M);

    cout << "[OK] t=" << setw(5) << t[i] << endl;

    if (t[i] < tmin)
    {
        tmin = t[i];
        vec_size = K;
        ofstream("grad_vector_size.dat") << K;
    }
}

template <>
void runTest<-1>(int N, int M, clock_t t[])
{
}

int main(int argc, char* argv[])
{
    int N = 1024;
    int M = 1000000;

    if (argc > 1)
        N = atoi(argv[1]);

    if (argc > 2)
        M = atoi(argv[2]);

    const int K = 3;
    clock_t t[K];

    ofstream("grad_vector_size.dat") << 0;

    runTest<K>(N, M, t);

    cout << "result: optimum(vector_size)=" << vec_size << endl;
    return 0;
}					],
					[use_vector_flags=enabled],
					[use_vector_flags=disabled; printf "[FAILED]\n"])

		dnl If AVX is working, we must keep -march=native in order to avoid compiler errors.
		
		if test $use_vector_flags != enabled; then
		   CXXFLAGS="$save_CXXFLAGS"
		fi
		
		if test $use_vector_flags == enabled && test -f grad_vector_size.dat; then
			ol_enable_autodiff_vect=`cat grad_vector_size.dat`

			if test -z "$ol_enable_autodiff_vect"; then
				ol_enable_autodiff_vect=0
			fi

			rm grad_vector_size.dat
		fi

		AC_LANG_C
	fi

	case "$ol_enable_autodiff_vect" in
	16)
		AC_DEFINE(GRADIENT_VECTOR_REGISTER_SIZE, 16, [the number of bytes of an avx or sse register])
		AC_MSG_RESULT([vectorization of derivatives will use 16 byte vectors])
		;;
	32)
		AC_DEFINE(GRADIENT_VECTOR_REGISTER_SIZE, 32, [the number of bytes of an avx or sse register])
		AC_MSG_RESULT([vectorization of derivatives will use 32 byte vectors])
		;;
	64)
		AC_DEFINE(GRADIENT_VECTOR_REGISTER_SIZE, 64, [the number of bytes of an avx or sse register])
		AC_MSG_RESULT([vectorization of derivatives will use 64 byte vectors])
		;;
	*)
		AC_DEFINE(GRADIENT_VECTOR_REGISTER_SIZE, 0, [the number of bytes of an avx or sse register])
		AC_MSG_WARN([vectorization of derivatives is not possible])
	esac
fi

if test "$ol_enable_sparse_autodiff" != no ; then
	AC_DEFINE(USE_SPARSE_AUTODIFF, 1, [define if automatic differentiation based on sparse algorithms can be used])
fi

dnl ----------------------------------------------------------------
dnl
dnl Enable blitz++ support <http://oonumerics.org/blitz/>
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_HEADER(blitz/blitz.h,[have_blitz_blitz_h=yes],[have_blitz_blitz_h=no])
AC_CHECK_HEADER(blitz/array.h,[have_blitz_array_h=yes],[have_blitz_array_h=no])
AC_CHECK_HEADER(blitz/tinyvec.h,[have_blitz_tinyvec_h=yes],[have_blitz_tinyvec_h=no])
AC_CHECK_HEADER(blitz/tinyvec-et.h,[have_blitz_tinyvecet_h=yes],[have_blitz_tinyvecet_h=no])
AC_CHECK_HEADER(blitz/tinymat.h,[have_blitz_tinymat_h=yes],[have_blitz_tinymat_h=no])
AC_CHECK_HEADER(blitz/matrix.h,[have_blitz_matrix_h=yes],[have_blitz_matrix_h=no])
AC_LANG_C

if test "$have_blitz_blitz_h" = "yes" \
	-a "$have_blitz_array_h" = "yes" \
	-a "$have_blitz_tinyvec_h" = "yes" \
	-a "$have_blitz_tinyvecet_h" = "yes" \
	-a "$have_blitz_tinymat_h" = "yes" \
	-a "$have_blitz_matrix_h" = "yes"; \
then
	AC_DEFINE(HAVE_BLITZ, 1, [define if blitz++ can be used])
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

	if test "$have_librt" = "yes" ; then
		AC_DEFINE(USE_RT,1,[define to use POSIX realtime])

        dnl ----------------------------------------------------------------
        dnl
        dnl Check for clock_gettime()
        dnl
        AC_TRY_LINK([
        #include <time.h>
                    ],[
            struct timespec tp;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);
                    ],[use_clock_gettime=yes],[use_clock_gettime=no])

        if test $use_clock_gettime = yes; then
            AC_DEFINE(HAVE_CLOCK_GETTIME, 1, [define if clock_gettime can be used])
        fi
	elif test "$ol_with_rt" = "yes" ; then
		AC_MSG_ERROR([POSIX realtime is not available])
	fi

    LIBS="$save_LIBS"
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

	RTAI_CONFIG=`which rtai-config`
	RC=$?
	if test $RC = 0 ; then
		RTAI_LDFLAGS=`$RTAI_CONFIG --lxrt-ldflags`
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for noexcept
dnl
AC_LANG_CPLUSPLUS
have_noexcept=no
AC_MSG_CHECKING([for noexcept])
AC_TRY_COMPILE([int pippo() noexcept {return 1;}],[],[have_noexcept=yes],[have_noexcept=no])
if test $have_noexcept != no ; then
	AC_DEFINE(HAVE_NOEXCEPT,1, [define if the compiler accept noexcept])
	AC_MSG_RESULT([yes])
else
	dnl AC_MSG_RESULT([no])
	AC_MSG_ERROR([The C++ compiler does not support the noexcept keyword. You should upgrade to a more recent compiler supporting c++11.])
fi
AC_LANG_C

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
AC_CHECK_LIB([ws2_32], [main])
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

		if test "$have_socklen_t" = "no" ; then
			AC_TRY_COMPILE([#include <ws2tcpip.h>],[
				socklen_t l;
			],[have_socklen_t=yes],[have_socklen_t=no])
		fi

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
dnl Check for threads (check anyway...)
dnl
have_threads="no"
if test "$ol_with_threads" != "no" ; then
	AX_PTHREAD([have_threads="yes"],[have_threads="no"])
	if test "$ol_with_threads" = "yes" ; then
		if test "$have_thread" = "no" ; then
			AC_MSG_ERROR([no pthread thread capability detected])
		fi
	fi
	if test "$have_threads" = "yes" ; then
		AC_MSG_CHECKING([for threads])
		save_LIBS="$LIBS"
		save_CXXFLAGS="$CXXFLAGS"
		save_CFLAGS="$CFLAGS"
		save_CC="$CC"
		LIBS="$LIBS $PTHREAD_LIBS"
		CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
		CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
		LDFLAGS="$LDFLAGS  $PTHREAD_CFLAGS"
		CC="$PTHREAD_CC"
		AC_TRY_LINK([
				/* Check whether threads really works */
				#include <pthread.h>
			],[
			return pthread_create(0, 0, 0, 0);
			],[
			],[
				LIBS="$save_LIBS"
				CXXFLAGS="$save_CXXFLAGS"
				CFLAGS="$save_CFLAGS"
				CC="$save_CC"
				have_threads = "no"
				if test "$ol_with_threads" = "yes" ; then
					AC_MSG_ERROR([no pthread thread linking capability detected])
				fi
		])
		if test "$have_threads" = "yes" ; then
			AC_DEFINE(HAVE_THREADS,1,[define if usable threads are available])
			AC_MSG_RESULT([yes])
		fi

		AC_MSG_CHECKING([for usable pthread_setaffinity_np])
		AC_TRY_LINK([
				#define _GNU_SOURCE
				#include <pthread.h>
				#include <sched.h>
			],[
				cpu_set_t cpuset;
				pthread_attr_t attr;

				CPU_ZERO(&cpuset);
				CPU_SET(0, &cpuset);

				pthread_getattr_np(pthread_self(), &attr);
				pthread_setaffinity_np(pthread_self(), sizeof(cpuset), &cpuset);
				return pthread_getaffinity_np(pthread_self(), sizeof(cpuset), &cpuset);
			],[
				AC_DEFINE(USE_PTHREAD_SETAFFINITY_NP,1,[define to use pthread_setaffinity_np])
				AC_MSG_RESULT([yes])
			],[
				AC_MSG_RESULT([no])
		])

		AC_MSG_CHECKING([for usable pthread_sigmask])
		AC_TRY_LINK([
				#include <pthread.h>
				#include <signal.h>
			],[
				sigset_t newset;
				sigemptyset(&newset);
				sigaddset(&newset, SIGTERM);
				pthread_sigmask(SIG_BLOCK, &newset, NULL);
			],[
				AC_DEFINE(HAVE_PTHREAD_SIGMASK,1,[define to use pthread_sigmask])
				AC_MSG_RESULT([yes])
			],[
				AC_MSG_RESULT([no])
		])
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
			eval `libtool --tag=CC --config|grep '^wl='`
			LIBS="$LIBS -lltdl"
			EXPORT_DYNAMIC_FLAG_SPEC=`(./libtool --config; echo eval echo \\$export_dynamic_flag_spec) | sh`
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
		if test "$have_runtime_loading" != "yes" -a "x$MODULES_DIRS" != "x" ; then
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
dnl Check for build mbc library
dnl
dnl do not build libmbc if sockets are not available (fixes issue in MinGW)
AM_CONDITIONAL(BUILD_MBC,[test "$ol_enable_mbc" = "yes" -a "$have_socket" = "yes"])

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

	if test "${ol_with_openblas_lib:-set}" = set ; then
		dnl in case one relies on a symbolic link named libopenblas.{a,so}
		ol_with_openblas_lib="openblas"
	fi

	try_blas_LIBS=
	if test "$have_blas" = "no" ; then
		try_blas_LIBS="-l$ol_with_openblas_lib $FCLIBS $F77LIBS -lm"
		if test "$ol_with_blas" = "auto" -o "$ol_with_blas" = "openblas" ; then
			LIBS="$try_blas_LIBS $LIBS"

			AC_TRY_LINK([
	/* FIXME: some architectures might use different
	 * naming conventions, e.g. uppercase, no underscore
	 */
	extern int dgemm_(void);
			],[
	/* test for blas */
	return dgemm_();
			],[have_blas="openblas"],[
				have_blas="no"
				if test "$ol_with_blas" = "openblas" ; then
					AC_MSG_ERROR([OpenBLAS implementation not available])
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
				have_cblas="yes"
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
	#include<cblas.h>
			],[
	/* test for blas */
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 1, 1, 1, 1, NULL, 1, NULL, 1, 0., NULL, 1);
	return 0;
			],[have_blas="cblas"
				have_cblas="yes"
			],[
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
		if test "$have_cblas" == "yes" ; then
			AC_DEFINE(HAVE_CBLAS,1,[define if CBLAS are available])
		fi
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
new_metis_interface=no
if test "$ol_with_metis" != "no" ; then
	AC_CHECK_LIB(metis, METIS_PartGraphVKway,
		[ol_link_metis=yes],[ol_link_metis=no],[-lm])
	if test "$ol_link_metis" != "yes" ; then
		AC_CHECK_LIB(metis, METIS_PartGraphKway,
			[ol_link_metis=yes],[ol_link_metis=no],[-lm])
	fi

	AC_CHECK_HEADERS(metis.h)

	if test "$ac_cv_header_metis_h" = "yes" ; then
		if test "$ol_link_metis" != "no" ; then
			AC_DEFINE(USE_METIS,1,
				[define if you want to use Metis])
			AC_TRY_COMPILE([#include "metis.h"], [idx_t i;], [new_metis_interface=yes],[new_metis_interface=no] )
			if test "$new_metis_interface" != "no" ; then
				AC_DEFINE(NEW_METIS_INTERFACE,1,[define if you want to use new Metis API])
			else
				AC_DEFINE(NEW_METIS_INTERFACE,0,[define if you want to use new Metis API])
			fi
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

        dnl Some optimized BLAS libraries like OpenBLAS already contain LAPACK
        try_lapack_LIBS="$BLAS_LIBS"

        LIBS="$try_lapack_LIBS $save_LIBS"
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

	if test "$have_lapack" == "no"; then
            LIBS=$save_LIBS
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
	fi

    if test "$have_lapack" != "no" ; then
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
                        ],[AC_DEFINE(HAVE_DGGEV,1,[define to use LAPACK's dggev in Eig])])
	fi

	if test "$have_lapack_" = "yes" ; then
		AC_TRY_LINK([
	extern double dgecon_(double *);
			],[
	double i = 1., o;
	o = dgecon_(&i);
                        ],[AC_DEFINE(HAVE_DGECON,1,[define to use LAPACK's dgecon in newton raphson solver])])
	elif test "$have_lapack" = "yes" ; then
		AC_TRY_LINK([
	extern double dgecon(double *);
			],[
	double i = 1., o;
	o = dgecon(&i);
			],[AC_DEFINE(HAVE_DGECON,1,[define to use LAPACK's dgecon in newton raphson solver])])
	fi

	if test "$have_lapack_" = "yes" ; then
		AC_TRY_LINK([
	extern double dlamch_(const char *);
			],[
	double o;
	o = dlamch_("S");
                        ],[AC_DEFINE(HAVE_DLAMCH,1,[define to use LAPACK's dlamch in dgeequ])])
	elif test "$have_lapack" = "yes" ; then
		AC_TRY_LINK([
	extern double dlamch(const char *);
			],[
	double o;
	o = dlamch("S");
			],[AC_DEFINE(HAVE_DLAMCH,1,[define to use LAPACK's dlamch in dgeequ])])
	fi

	if test "$have_lapack_" = "yes" ; then
		AC_TRY_LINK([
	extern double dgetrf_(double *);
			],[
	double i = 1., o;
	o = dgetrf_(&i);
                        ],[AC_DEFINE(HAVE_DGETRF,1,[define to use LAPACK's dgetrf in newton raphson solver])])
	elif test "$have_lapack" = "yes" ; then
		AC_TRY_LINK([
	extern double dgetrf(double *);
			],[
	double i = 1., o;
	o = dgetrf(&i);
			],[AC_DEFINE(HAVE_DGETRF,1,[define to use LAPACK's dgetrf in newton raphson solver])])
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
                        ],[AC_DEFINE(HAVE_DGEGV,1,[define to use LAPACK's dgegv in Eig])])
        fi

        if test "$have_lapack_" = "yes" ; then
                AC_TRY_LINK([
        extern double dtrtrs_(double *);
                        ],[
        double i = 1., o;
        o = dtrtrs_(&i);
                        ],[have_dtrtrs=yes
                           AC_DEFINE(HAVE_DTRTRS,1,[define to use LAPACK's dtrtrs in SPQR])])
        elif test "$have_lapack" = "yes" ; then
                AC_TRY_LINK([
        extern double dtrtrs(double *);
                        ],[
        double i = 1., o;
        o = dtrtrs(&i);
                        ],[have_dtrtrs=yes
                           AC_DEFINE(HAVE_DTRTRS,1,[define to use LAPACK's dtrtrs in SPQR])])
	fi

	if test "$have_lapack" = "yes" ; then
		ol_link_lapack=yes
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
	fi

	LIBS="$save_LIBS"
fi

if test $ol_link_lapack = no ; then
	case "$ol_with_lapack" in
	yes)
		AC_MSG_ERROR([Could not locate Lapack])
		;;
	auto)
		AC_MSG_WARN([Lapack not supported])
		;;
	*)
		;;
	esac
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for Umfpack
dnl
ol_link_umfpack=no
if test "$ol_with_umfpack" != "no" ; then
	AC_CHECK_HEADERS(umfpack.h)

	have_suitesparseconfig="no"
	if test "$ac_cv_header_umfpack_h" = "yes" ; then

		dnl check for SuiteSparseConfig
		AC_MSG_CHECKING([for SuiteSparse_config library])
		save_LIBS="$LIBS"
		try_suitesparseconfig_LIBS="-lsuitesparseconfig -lm $RT_LIBS"
		LIBS="$try_suitesparseconfig_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for suitesparseconfig */
	return 0;
		],[have_suitesparseconfig="yes"],[have_suitesparseconfig="no"])
		if test "$have_suitesparseconfig" = "yes" ; then
			try_umfpack_LIBS="$try_suitesparseconfig_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_umfpack_LIBS=
			LIBS="$save_LIBS"
			AC_MSG_RESULT([no])
		fi

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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
			AC_MSG_RESULT([no])
		fi

		have_umfpack_version="no"
		have_umfpack="no"
		AC_MSG_CHECKING([for umfpack library version])
		AC_TRY_COMPILE([
#if defined(HAVE_STDINT_H)
#include <stdint.h>
#elif defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#else
#error "need int32_t"
#endif
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
#if defined(HAVE_STDINT_H)
#include <stdint.h>
#elif defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#endif
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
#if defined(HAVE_STDINT_H)
#include <stdint.h>
#elif defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#endif
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

	have_suitesparseconfig="no"
	if test "$ac_cv_header_klu_h" = "yes" ; then

		dnl check for SuiteSparseConfig
		AC_MSG_CHECKING([for SuiteSparse_config library])
		save_LIBS="$LIBS"
		try_suitesparseconfig_LIBS="-lsuitesparseconfig $RT_LIBS -lm"
		LIBS="$try_suitesparseconfig_LIBS $LIBS"

		AC_TRY_LINK([
	/* no test: the library must exist, then we trust it */
		],[
	/* test for suitesparseconfig */
	return 0;
		],[have_suitesparseconfig="yes"],[have_suitesparseconfig="no"])
		if test "$have_suitesparseconfig" = "yes" ; then
			try_umfpack_LIBS="$try_suitesparseconfig_LIBS"
			AC_MSG_RESULT([yes])
		else
			try_umfpack_LIBS=
			LIBS="$save_LIBS"
			AC_MSG_RESULT([no])
		fi

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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
			LIBS="$save_LIBS"
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
dnl Check for Intel MKL Pardiso
dnl
ol_link_pardiso=no
if test $ol_with_pardiso != no ; then
        AC_CHECK_HEADER(mkl/mkl_pardiso.h, [have_mkl_pardiso_h=yes], [have_mkl_pardiso_h=no])
        AC_CHECK_HEADER(omp.h, [have_omp_h=yes], [have_omp_h=no])
        if test "$have_mkl_pardiso_h" = "yes" ; then
	        try_PARDISO_LIBS="-lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lm -lgomp"
	        AC_CHECK_LIB(mkl_gf_lp64,pardiso_,[have_mkl_pardiso_lib=yes],[have_mkl_pardiso_lib=no],[$try_PARDISO_LIBS])
                AC_CHECK_LIB(gomp,omp_set_num_threads,[have_omp_set_num_threads=yes],[have_omp_set_num_threads=no],[$try_PARDISO_LIBS])
	        if test "$have_mkl_pardiso_lib" = "yes"; then
		        AC_DEFINE(USE_PARDISO,1,[define to use Intel MKL Pardiso library])
                        ol_link_pardiso=yes
		        PARDISO_LIBS="$try_PARDISO_LIBS"
	        else
                        if test "$ol_with_pardiso" = "yes" ; then
			   AC_MSG_ERROR([Intel MKL pardiso library not available])
		        fi
	        fi

                if test "$have_omp_h" = "yes" && test "$have_omp_set_num_threads" = "yes" ; then
                   AC_DEFINE(USE_OMP_SET_NUM_THREADS,1,[define to use omp_set_num_threads])
                fi
        fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for IBM Watson wsmp
dnl
ol_link_wsmp=no
if test $ol_with_wsmp != no ; then
	try_wsmp_LIBS="-lwsmp64"
	AC_CHECK_LIB(wsmp64,wgsmp_,[have_wsmp=yes],[have_wsmp=no],[$try_wsmp_LIBS])
	if test "$have_wsmp" = "yes" ; then
		AC_DEFINE(USE_WSMP,1,[define to use IBM Watson wsmp library])
                ol_link_wsmp=yes
		WSMP_LIBS="$try_wsmp_LIBS"
	else
		if test "$ol_with_wsmp" = "yes" ; then
			AC_MSG_ERROR([IBM Watson wsmp library not available])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Y12
dnl
ol_link_y12=no
if test $ol_with_y12 != no ; then
	if test "$have_fc_f77" = "no" ; then
		if test $ol_with_y12 = yes ; then
			AC_MSG_ERROR([Y12 sparse LU library needs Fortran compiler])
		fi
	else
		ol_link_y12=yes
		AC_DEFINE(USE_Y12,1,[define to use Y12 sparse LU library])
		Y12_LIBS="${top_builddir}/libraries/liby12/liby12.la -lm"
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Pastix
dnl
ol_link_pastix=no
if test $ol_with_pastix != no ; then
        have_pastix="no"

	if test -z "${PASTIX_INC}"; then
     	        if ! test -z "${PKG_CONFIG}"; then
     	                PASTIX_INC="`${PKG_CONFIG} --cflags pastix`"
     	        fi
	fi

   	if test -z "${PASTIX_LIBS}"; then
     	        if ! test -z "${PKG_CONFIG}"; then
                        PASTIX_LIBS="`${PKG_CONFIG} --libs pastix`"
                else
                        PASTIX_LIBS="-lpastix -lpastix_kernels -lspm -lhwloc"
                fi
        fi

        save_CPPFLAGS="${CPPFLAGS}"
	CPPFLAGS="${PASTIX_INC} ${CPPFLAGS}"

	AC_CHECK_HEADER([pastix.h])

	if test "${ac_cv_header_pastix_h}" = yes; then
	        save_LIBS="$LIBS"
	        LIBS="${PASTIX_LIBS} $save_LIBS -lm"

	        AC_CHECK_FUNCS([pastixInit \
	                        pastixInitParam \
			        pastix_task_analyze \
			        pastix_task_numfact \
			        pastix_task_solve \
				pastix_task_refine \
				pastixFinalize \
				spmNorm \
				spmScalMatrix \
				spmScalVector \
				spmExit \
				spmCheckAxb],[have_pastix=yes],[have_pastix=no])

		LIBS="$save_LIBS"
	fi

	CPPFLAGS="${save_CPPFLAGS}"
	
	if test "${have_pastix}" = "yes"; then
		ol_link_pastix=yes
	        AC_DEFINE([USE_PASTIX],1,[Define if PaStiX can be used])
	else
		PASTIX_LIBS=""
		PASTIX_INC=""
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for QRUPDATE
dnl
ol_link_qrupdate=no
if test $ol_with_qrupdate != no ; then
   if test $ol_with_suitesparseqr = yes ; then
      AC_MSG_NOTICE([suitesparseqr needs qrupdate; promoting qrupdate checking from auto to yes])
      ol_with_qrupdate=yes
   fi
   AC_MSG_CHECKING([for qrupdate library])
   try_qrupdate_LIBS="-lqrupdate"
   save_LIBS="$LIBS"
   have_qrupdate="no"
   LIBS="$try_qrupdate_LIBS $save_LIBS -lm"
   AC_TRY_LINK([
        extern int dqr1up(double*);
   ],[
        double i = 1;
        int o;
        o = dqr1up(&i);
   ],[have_qrupdate="yes"])

   if test "$have_qrupdate" = "no"; then
      AC_TRY_LINK([
        extern int dqr1up_(double*);
      ],[
        double i = 1;
        int o;
        o = dqr1up_(&i);
      ],[have_qrupdate="yes"])
   fi

   LIBS="$save_LIBS"

   if test "$have_qrupdate" = "yes"; then
      AC_MSG_RESULT([yes])
      ol_link_qrupdate=yes
      QRUPDATE_LIBS="$try_qrupdate_LIBS"
      AC_DEFINE(USE_QRUPDATE,1,[define if you want to use QRUPDATE])
   else
      AC_MSG_RESULT([no])
   fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for SuiteSparseQr
dnl
ol_link_suitesparseqr=no
if test $ol_with_suitesparseqr != no ; then
        AC_LANG_CPLUSPLUS
        AC_CHECK_HEADERS(SuiteSparseQR.hpp)

        if test "$ac_cv_header_SuiteSparseQR_hpp" = "yes" ; then
                AC_MSG_CHECKING([for SuiteSparseQR library])
                try_suitesparseqr_LIBS="-lspqr -lcholmod"
                save_LIBS="$LIBS"
                have_suitesparseqr="no"
                LIBS="$try_suitesparseqr_LIBS $save_LIBS -lm"
                AC_TRY_LINK([
                        #include <SuiteSparseQR.hpp>
                        ],[
                        cholmod_sparse *A, *Q, *R;
                        SuiteSparse_long *E;
                        cholmod_common oCommon;
                        cholmod_l_start(&oCommon);
                        A = cholmod_l_speye(100, 100, CHOLMOD_REAL, &oCommon);
                        SuiteSparseQR<double>(SPQR_ORDERING_DEFAULT,
                                              SPQR_NO_TOL,
                                              A->ncol,
                                              A,
                                              &Q,
                                              &R,
                                              &E,
                                              &oCommon);
                        cholmod_l_free_sparse(&A, &oCommon);
                        cholmod_l_free_sparse(&Q, &oCommon);
                        cholmod_l_free_sparse(&R, &oCommon);
                        cholmod_l_free(A->nrow, sizeof(*E), E, &oCommon);
                        cholmod_l_finish(&oCommon);
                ],[have_suitesparseqr="yes"])

                LIBS="$save_LIBS"

                if test "$have_suitesparseqr" = "yes"; then
                   AC_MSG_RESULT([yes])

                   if test "$have_qrupdate" = "yes" && test "$have_dtrtrs" = "yes" ; then
                     ol_link_suitesparseqr=yes
                     SUITESPARSEQR_LIBS="$try_suitesparseqr_LIBS"
                     AC_DEFINE(USE_SUITESPARSE_QR,1,[define if you want to use SuiteSparseQR])
                   fi
                else
                   AC_MSG_RESULT([no])
                fi
        fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for STRUMPACK
dnl
ol_link_strumpack=no
if test $ol_with_strumpack != no ; then
        AC_LANG_CPLUSPLUS
        AC_CHECK_HEADERS(StrumpackSparseSolver.hpp)

        if test "$ac_cv_header_StrumpackSparseSolver_hpp" = "yes" ; then
                AC_MSG_CHECKING([for STRUMPACK library])
                try_strumpack_LIBS="-lstrumpack"
                save_LIBS="$LIBS"
                have_strumpack="no"
                LIBS="$try_strumpack_LIBS $save_LIBS"
                AC_TRY_LINK([
                        #include <StrumpackSparseSolver.hpp>
			 ],[
			   double Ax[] = {1., 2., 3., 4.};
  			   int Ar[] = {0, 1, 2, 3, 4};
			   int Ac[] = {0, 1, 2, 3};
			   double x[4];
			   double b[4];
			   using namespace strumpack;
			   StrumpackSparseSolver<double, int> sp;							 
                           sp.options().set_compression(CompressionType::HSS);
	                   sp.options().set_rel_tol(1e-8);
	                   sp.options().set_abs_tol(1e-8);
	                   sp.options().set_gmres_restart(10);     
			   sp.set_csr_matrix(4, Ar, Ac, Ax);
			   sp.reorder();
			   sp.factor();
			   sp.solve(b, x);
                ],[have_strumpack="yes"])

                LIBS="$save_LIBS"

                if test "$have_strumpack" = "yes"; then
                   AC_MSG_RESULT([yes])
                   ol_link_strumpack=yes
                   STRUMPACK_LIBS="$try_strumpack_LIBS"
                   AC_DEFINE(USE_STRUMPACK,1,[define if you want to use STRUMPACK])
                else
                   AC_MSG_RESULT([no])
                fi
        fi
fi

AC_LANG_C
if test "$have_fc_f77" = "no" ; then
   	if test $ol_with_y12 = yes ; then
           	AC_MSG_ERROR([Y12 sparse LU library needs Fortran compiler])
	fi
else
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
#include <$mpi_cxx_header>
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
		$PKG_CONFIG ginac
		if test $? = 0; then
			GINACLIB_LIBS=`$PKG_CONFIG --libs ginac`
			GINACLIB_CPPFLAGS=`$PKG_CONFIG --cflags ginac`
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
	feenableexcept \
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
	sched_setaffinity \
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
dnl Add EXPORT_DYNAMIC_FLAG_SPEC to LDFLAGS immediately before exit
dnl otherwise all the tests during configure would fail

LDFLAGS="$LDFLAGS $EXPORT_DYNAMIC_FLAG_SPEC"

if test "$ol_enable_Werror" != no ; then
	CFLAGS="$CFLAGS -Werror -Wno-error=cpp -Wno-error=format-truncation"
	CXXFLAGS="$CXXFLAGS -Werror -Wno-error=cpp -Wno-error=format-truncation"
	FFLAGS="$FFLAGS"
fi


dnl ----------------------------------------------------------------
dnl TODO: add a check about whether a fortran (F77, F95) compiler is
dnl needed, and bail out if not available

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
AC_SUBST(PARDISO_LIBS)
AC_SUBST(WSMP_LIBS)
AC_SUBST(HARWELL_LIBS)
AC_SUBST(PASTIX_INC)
AC_SUBST(PASTIX_LIBS)
AC_SUBST(QRUPDATE_LIBS)
AC_SUBST(SUITESPARSEQR_LIBS)
AC_SUBST(STRUMPACK_LIBS)
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
AC_SUBST(RTAI_LDFLAGS)
AC_SUBST(ANN_LIBS)
AC_SUBST(BULLET_LIBS)
AC_SUBST(PYTHON_NUMPY_HEADERS)
AC_SUBST(OCTAVE_LIBS)
AC_SUBST(OCTAVE_INCLUDE)
AC_SUBST(OCTAVE_MKOCTFILE)
AC_SUBST(MODULES_DIRS)
AC_SUBST(MODULE_LOADPATH)
AC_SUBST(wl)

dnl AC_SUBST(GINACLIB_CONFIG)
AC_SUBST(GINACLIB_CPPFLAGS)
AC_SUBST(GINACLIB_LIBS)

AC_SUBST(DATE)

dnl conditional build declarations
AM_CONDITIONAL(MBDYN_DEVEL,[test "$MBDYN_DEVEL" = "yes"])

AM_CONDITIONAL(USE_ARPACK,[test "$ol_link_arpack" = yes])
AM_CONDITIONAL(USE_CHACO,[test "$ol_link_chaco" = yes])
AM_CONDITIONAL(USE_GINAC,[test "$ol_link_ginac" = yes])
AM_CONDITIONAL(USE_JDQZ,[test "$ol_link_jdqz" = yes])
AM_CONDITIONAL(USE_METIS,[test "$ol_link_metis" = yes])
AM_CONDITIONAL(NEW_METIS_INTERFACE,[test "$new_metis_interface" = yes])
AM_CONDITIONAL(USE_NETCDF4,[test "$have_libnetcdf4" = yes])
AM_CONDITIONAL(USE_OCTAVE,[test "$use_octave" = yes])
AM_CONDITIONAL(USE_OCTAVE_UTILS,[test "$use_octave_utils" = yes])
AM_CONDITIONAL(USE_PYTHON,[test "$use_python" = yes])
AM_CONDITIONAL(USE_RTAI,[test "$use_rtai" = yes])
AM_CONDITIONAL(USE_SCHUR,[test "$use_schur" = yes])
AM_CONDITIONAL(USE_SOCKET,[test "$use_socket" = yes])
AM_CONDITIONAL(USE_ANN,[test "$ol_link_ann" = yes])
AM_CONDITIONAL(USE_AUTODIFF,[test "$ol_enable_autodiff" = yes])
AM_CONDITIONAL(USE_SPARSE_AUTODIFF,[test "$ol_enable_sparse_autodiff" = yes])
AM_CONDITIONAL([USE_PASTIX],[test "$ol_link_pastix" = "yes"])

AM_CONDITIONAL(FOOL_LIBTOOL,[test "$fool_libtool" = yes])

AM_CONDITIONAL(USE_FORTRAN,[test "$have_fc_f77" = yes])
AM_CONDITIONAL(BUILD_Y12,[test "$ol_link_y12" = yes])
AM_CONDITIONAL(BUILD_AEROD2_F,[test "$have_fc_f77" = yes])

AM_CONDITIONAL(BUILD_STATIC_MODULES,[test "$build_static_modules" = yes])
AM_CONDITIONAL(BUILD_CHARM,[test "$build_charm" = yes])

AM_CONDITIONAL(HAVE_TERMIOS,[test "$have_termios" = yes])
AM_CONDITIONAL(HAVE_GETPASS,[test "$have_getpass" = yes])
AM_CONDITIONAL(INSTALL_TEST_PROGRAMS,[test "$ol_enable_install_test_progs" = yes])

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
AC_CONFIG_FILES([Makefile \
 	libraries/Makefile \
 	libraries/libmbutil/Makefile \
 	libraries/libobjs/Makefile \
 	libraries/libcolamd/Makefile \
 	libraries/libnaive/Makefile \
 	libraries/liby12/Makefile \
 	libraries/libmbmath/Makefile \
 	libraries/libmbwrap/Makefile \
 	libraries/libann/Makefile \
 	libraries/libmbc/Makefile \
 	mbdyn/Makefile \
 	mbdyn/base/Makefile \
 	mbdyn/aero/Makefile \
 	mbdyn/elec/Makefile \
 	mbdyn/thermo/Makefile \
 	mbdyn/hydr/Makefile \
 	mbdyn/struct/Makefile \
 	modules/Makefile \
 	utils/Makefile \
 	etc/Makefile \
 	etc/modal.d/Makefile \
 	etc/pam.d/Makefile \
 	var/Makefile \
 	manual/Makefile \
 	manual/femgen.1 \
 	manual/mbdyn.1 \
 	include/Makefile \
 	include/ac/Makefile
 	build/Makefile])
AC_OUTPUT
