dnl Process this file with autoconf to produce a configure script.
dnl
dnl MBDyn (C) is a multibody analysis code.
dnl http://www.mbdyn.org
dnl
dnl Copyright (C) 1996-2003
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
# Copyright (C) 1996-2003
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
echo "	Copyright (C) 1996-2003"
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
$ol_unset CDPATH || test "${CDPATH+set}" != set || CDPATH=: && export CDPATH

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

AM_CONFIG_HEADER(include/mbconfig.h)

top_builddir=`pwd`
AC_SUBST(top_builddir)dnl

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
OL_ARG_ENABLE(debug_mpi,[  --enable-debug-mpi      enable MPI debugging], no)dnl
OL_ARG_ENABLE(socket_drives,[  --enable-socket-drives  enable socket drives], auto)dnl
OL_ARG_ENABLE(adams,[  --enable-adams          enable MSC.ADAMS output], no)dnl
OL_ARG_ENABLE(motionview,[  --enable-motionview     enable Altair's Motion View output], no)dnl

dnl ----------------------------------------------------------------
dnl General "with" options
OL_ARG_WITH(debug_mode,[  --with-debug-mode       with debug mode],none,[none memory])
OL_ARG_WITH(tcl,[  --with-tcl              with tcl interpreters],auto,[auto yes no])
OL_ARG_WITH(libf2c,[  --with-libf2c[={f2c|g2c}]  with f2c library],auto,[auto yes no g2c])
OL_ARG_WITH(fs,[  --with-fs[={unix|dos}]  filesystem type],unix,[unix dos])
OL_ARG_WITH(mpi,[  --with-mpi              with MPI support],auto,[auto yes no])
OL_ARG_WITH(metis,[  --with-metis            with Metis support],auto,[auto yes no])
OL_ARG_WITH(rtai,[  --with-rtai              with RTAI support],no,[auto yes no])
OL_ARG_WITH(meschach,[  --with-meschach         with Meschach math library],auto,[auto yes no])
OL_ARG_WITH(umfpack,[  --with-umfpack         with Umfpack math library],auto,[auto yes no])
OL_ARG_WITH(y12,[  --with-y12              with Y12 sparse math library],auto,[auto yes no])
OL_ARG_WITH(harwell,[  --with-hsl          with HSL (Harwell) sparse math library],auto,[auto yes no])
OL_ARG_WITH(lapack,[  --with-lapack           with LAPACK math library],auto,[auto yes no])
OL_ARG_WITH(pam,[  --with-pam              with PAM support],auto,[auto yes no])
OL_ARG_WITH(ginac,[  --with-ginac            with GiNaC support],auto,[auto yes no])

OL_ARG_WITH(aero_output,[  --with-aero-output={std,gauss,node}  aerodynamic output mode],auto,[auto std gauss node])

OL_ARG_WITH(struct,[  --with-struct           with structural elements],yes,[yes no])
OL_ARG_WITH(elec,[  --with-elec             with electric stuff],yes,[yes no])
OL_ARG_WITH(aero,[  --with-aero             with aerodynamic stuff],yes,[yes no])
OL_ARG_WITH(hydr,[  --with-hydr             with hydraulic stuff],yes,[yes no])

dnl module with's
AC_ARG_WITH(module,[  --with-module=modulelist	with listed modules],[
	ol_arg=
	for ol_val in $withval ; do
		ol_arg="$ol_arg module-$ol_val"
	done
	if test "${ol_arg+set}" != set ; then
		AC_MSG_ERROR(bad value $withval for --with-module)
	fi
	MODULES_DIRS=$ol_arg
],[MODULES_DIRS=])dnl

dnl ----------------------------------------------------------------
dnl
dnl Validate general "enable" options
dnl


dnl ----------------------------------------------------------------
dnl
dnl Checks for programs.
dnl
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AM_DISABLE_SHARED dnl Remove when you want to build shared libs also
AM_PROG_LIBTOOL

dnl ----------------------------------------------------------------
dnl
dnl Initialize vars
dnl
F2C_LIBS=
MPI_LIBS=
METIS_LIBS=
Y12_LIBS=
UMFPACK_LIBS=
MESCHACH_LIBS=
HARWELL_LIBS=
LAPACK_LIBS=
SECURITY_LIBS=
TCL_LIBS=

GINACLIB_CONFIG=
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
AC_CHECK_HEADERS(float.h)
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS(math.h)
AC_CHECK_HEADERS(netdb.h)
AC_CHECK_HEADERS(netinet/in.h)
AC_CHECK_HEADERS(pwd.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_HEADERS(string.h)
AC_CHECK_HEADERS(sys/mman.h)
AC_CHECK_HEADERS(sys/socket.h)
AC_CHECK_HEADERS(sys/times.h)
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(values.h)


dnl ----------------------------------------------------------------
dnl
dnl Checks for C++ header files (takes care of ANSI/pre-ANSI naming)
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_HEADERS(iostream)
AC_CHECK_HEADERS(iostream.h)
AC_CHECK_HEADERS(iomanip)
AC_CHECK_HEADERS(iomanip.h)
AC_CHECK_HEADERS(fstream)
AC_CHECK_HEADERS(fstream.h)

AC_CHECK_HEADERS(exception)
AC_CHECK_HEADERS(stdexcept)
AC_LANG_C


dnl ----------------------------------------------------------------
dnl
dnl Sort out defines
dnl


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
dnl Use element iterator
dnl
AC_DEFINE(USE_ELEM_ITER,1,[define this to use dedicated element iterator])


dnl ----------------------------------------------------------------
dnl
dnl Set inter-language naming convention for calling Fortran routines
dnl
AC_DEFINE(USE_UNDERSCORE,1,[append underscore to fortran calls])


dnl ----------------------------------------------------------------
dnl
dnl define packages (will be made configurable)
dnl
dnl AC_DEFINE(USE_STRUCT_NODES,1,[use structural nodes])
dnl AC_DEFINE(USE_AERODYNAMIC_ELEMS,1,[use aerodynamic elements (requires USE_STRUCT_NODES)])
dnl AC_DEFINE(USE_ELECTRIC_NODES,1,[use electric nodes])
dnl AC_DEFINE(USE_HYDRAULIC_NODES,1,[use hydraulic nodes])


dnl ----------------------------------------------------------------
dnl
dnl Check for DEBUG
dnl
if test "$ol_enable_debug_mpi" != no ; then
	AC_DEFINE(DEBUG_MPI,1,[define this to add MPI debugging code])
fi
if test "$ol_enable_debug" != no -o "$ol_enable_debug_mpi" != no ; then
	AC_DEFINE(DEBUG,1,[define this to add debugging code])
	if test "$ol_with_debug_mode" = memory ; then
		AC_DEFINE(DEBUG_MEMMANAGER,1,[define this to add memory debugging code])
	fi
fi

if test "$ol_enable_adams" != no ; then
	AC_DEFINE(USE_ADAMS,1,[define this to enable ADAMS/View output])
fi

if test "$ol_enable_motionview" != no ; then
	AC_CHECK_HEADER(mbs_result_api.h,[have_mbs_result_api_h=yes],[have_mbs_result_api_h=no])
	if test $have_mbs_result_api_h = "yes" ; then
		AC_DEFINE(USE_MOTIONVIEW,1,[define this to enable Motion View output])
	else
		if test $ol_enable_motionview != auto ; then
			AC_MSG_ERROR([Could not locate mbs_result_api.h header])
		fi
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for socket
dnl
if test "$ol_enable_socket_drives" != no ; then
	AC_CHECK_FUNC(socket,[have_socket=yes],[have_socket=no])

	if test $have_socket != no ; then
		AC_DEFINE(USE_SOCKET_DRIVES,1,[define to enable socket drives])

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
		
	else
		if test $ol_enable_socket_drives != auto ; then
			AC_MSG_ERROR([cannot find sockets])
		else
			AC_MSG_WARN([cannot find sockets])
			AC_MSG_WARN([socket drives not supported!])
		fi
	fi
else
	AC_MSG_WARN([socket drives not supported!])
fi

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
if test $have_exceptions != no ; then
	AC_DEFINE(USE_EXCEPTIONS,1,[define if C++ compiler handles exceptions])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
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
dnl Check for aerodynamic output flags
dnl
if test $ol_with_aero_output = gauss ; then
	AC_DEFINE(AEROD_OUTPUT,1,[set to 0,1,2 for std,gauss,node])
elif test $ol_with_aero_output = node ; then
	AC_DEFINE(AEROD_OUTPUT,2,[set to 0,1,2 for std,gauss,node])
else
	AC_DEFINE(AEROD_OUTPUT,0,[set to 0,1,2 for std,gauss,node])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for packages
dnl
if test $ol_with_struct = yes ; then
        AC_DEFINE(USE_STRUCT_NODES,1,[set to 1 to use structural stuff])
fi
if test $ol_with_elec = yes ; then
	AC_DEFINE(USE_ELECTRIC_NODES,1,[set to 1 to use electric stuff])
fi
if test $ol_with_hydr = yes ; then
        AC_DEFINE(USE_HYDRAULIC_NODES,1,[set to 1 to use hydraulic stuff])
fi
if test $ol_with_aero = yes ; then
	if test $ol_with_struct != yes ; then
		AC_MSG_WARN([Need --with-struct to enable aerodynamic stuff])
	else
        	AC_DEFINE(USE_AERODYNAMIC_ELEMS,1,[set to 1 to use electric stuff])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for tcl
dnl
ol_link_tcl=no
if test $ol_with_tcl != no ; then
	OL_ADD_DIRS(tcl)

	AC_CHECK_HEADERS(tcl.h)

	AC_MSG_CHECKING([for tcl library])

	AC_CHECK_LIB(tcl,Tcl_CreateInterp,[have_libtcl=yes],[have_libtcl=no])
	if test $have_libtcl != no ; then
		AC_DEFINE(USE_TCL,1,[define if you want to use tcl])
		TCL_LIBS=-ltcl
	else
		if test $ol_with_tcl != auto ; then
			AC_MSG_ERROR([Could not locate tcl library])
		else
			AC_MSG_WARN([tcl not supported!])
		fi
	fi
else
	AC_MSG_WARN([tcl library not supported!])
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks for sstream
dnl
AC_LANG_CPLUSPLUS
AC_CHECK_HEADERS(sstream,[have_sstream=yes],[have_sstream=no])
AC_CHECK_HEADERS(strstream.h,[have_strstream=yes],[have_strstream=no])
AC_MSG_CHECKING([whether the system uses sstream or strstream.h])
if test $have_sstream != no ; then
	AC_MSG_RESULT([sstream])
else
	if test $have_strstream != no ; then
		AC_MSG_RESULT([strstream.h])
		AC_MSG_WARN([strstream.h is deprecated])
	else
		AC_MSG_ERROR([No sstream or strstream.h found])
	fi
fi
AC_LANG_C

dnl ----------------------------------------------------------------
dnl
dnl Check for f2c libraries
dnl
ol_link_libf2c=no
if test $ol_with_libf2c != no ; then
	AC_CHECK_HEADERS(g2c.h f2c.h)

	AC_MSG_CHECKING([for f2c library])
	AC_MSG_RESULT([])
		
	if test $ol_with_libf2c = g2c ; then
		libf2c_pool="g2c"
	else
		libf2c_pool="f2c g2c"
	fi

	for l in $libf2c_pool ; do
		AC_MSG_CHECKING([for lib$l])

		try_libf2c="$l"
		try_libf2c_LIBS="-l$l"
		save_LIBS=$LIBS
		LIBS="$LIBS $try_libf2c_LIBS -lm"
		AC_TRY_LINK([#include <$top_builddir/include/ac/f2c.h>],[
			doublereal i = 1., o;
			o = d_cos(&i);
		],[have_libf2c=yes],[have_libf2c=no])
		LIBS=$save_LIBS

		if test $have_libf2c != no ; then
			F2C_LIBS=$try_libf2c_LIBS
			AC_MSG_RESULT([yes])
			ol_link_libf2c=yes
			break
		else
			AC_MSG_RESULT([no])
		fi
	done

	if test $ol_link_libf2c = no ; then
	
		if test $ol_with_libf2c != auto ; then
			AC_MSG_ERROR([Could not locate f2c library])
		else
			AC_MSG_WARN([f2c library not supported!])
		fi
	else
		dnl Hack to allow different Fortran compatibility libraries
		dnl that apparently do not implement integer power of double
		AC_CHECK_LIB($try_libf2c,pow_di,[
			AC_DEFINE(HAVE_POW_DI,1,[define if pow_di() function is available])])
	fi
else
	AC_MSG_WARN([f2c library not supported!])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for UNIX fs !?!
dnl
case "$ol_with_fs" in
	unix)
		AC_DEFINE(USE_UNIX_FILE,1,
			[define this to use UNIX-style filenames])
		;;
	dos)	AC_DEFINE(USE_DOS_FILE,1,
		[define this to use DOS-style filenames])
		;;
esac

dnl ----------------------------------------------------------------
dnl
dnl Check for dlopen
dnl
have_loadable=no
AC_CHECK_HEADER(ltdl.h,[have_ltdl_h=yes],[have_ltdl_h=no])
if test $have_ltdl_h = "yes" ; then
	AC_DEFINE(HAVE_LTDL_H,1,[define if you have <ltdl.h>])
	AC_CHECK_LIB(ltdl, lt_dlinit,[have_lt_dlinit=yes],[have_lt_dlinit=no])
	if test $have_lt_dlinit = "yes" ; then
		LIBS="$LIBS -lltdl"
		have_loadable=yes
		MODULE_LOADPATH=$libexecdir
	fi
fi
if test $have_loadable = no ; then
	AC_CHECK_HEADER(dlfcn.h,[have_dlfcn_h=yes],[have_dlfcn_h=no])
	if test $have_dlfcn_h="yes" ; then
		AC_DEFINE(HAVE_DLFCN_H,1,[define if you have <dlfcn.h>])
		AC_CHECK_LIB(dl, dlopen,[have_dlopen_in_ldl=yes],[have_dlopen_in_ldl=no], -ldl)
		if test $have_dlopen_in_ldl = "yes" ; then
			LIBS="$LIBS -ldl"
			have_loadable=yes
		fi
	fi
fi
if test $have_loadable=yes ; then
	AC_DEFINE(HAVE_LOADABLE,1,[define if you want to use loadable modules])
else
	AC_MSG_WARN([no dynamic linking library detected; loadable elements will be disabled])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for Umfpack
dnl
ol_link_umfpack=no
if test $ol_with_umfpack != no ; then
	OL_ADD_DIRS(umfpack)
	
	AC_CHECK_HEADERS(umfpack.h)

	if test $ac_cv_header_umfpack_h = yes ; then

		AC_MSG_CHECKING([for umfpack library])
		try_umfpack_LIBS="-lumfpack -lcblas -latlas"
		save_LIBS=$LIBS
		LIBS="$LIBS $try_umfpack_LIBS"
		AC_TRY_LINK([
#include <umfpack.h>
],[
			double t;
			t = umfpack_timer();
		],[have_umfpack=yes],[have_umfpack=no])
		AC_TRY_LINK([
#include <umfpack.h>
],[
			double Control[UMFPACK_CONTROL];
			umfpack_dl_defaults(Control);
		],[have_umfpack4=yes],[have_umfpack4=no])
		LIBS=$save_LIBS

		if test $have_umfpack != no ; then
			UMFPACK_LIBS="$try_umfpack_LIBS"
			AC_DEFINE(USE_UMFPACK,1,[define if you want to use Umfpack])
			if test $have_umfpack4 != no ; then
				AC_DEFINE(HAVE_UMFPACK4,1,[define if Umfpack4 is available])
			fi
			ol_link_umfpack=yes
			if test $have_umfpack4 != no ; then
				AC_MSG_RESULT([(>= 4.0) yes])
			else
				AC_MSG_RESULT([yes])
			fi
		else
			try_umfpack_LIBS="-lumfpack"
			save_LIBS=$LIBS
			LIBS="$LIBS $try_umfpack_LIBS"
			AC_TRY_LINK([
#include <umfpack.h>
],[
				double t;
				t = umfpack_timer();
			],[have_umfpack=yes],[have_umfpack=no])
			AC_TRY_LINK([
#include <umfpack.h>
],[
				double Control[UMFPACK_CONTROL];
				umfpack_dl_defaults(Control);
			],[have_umfpack4=yes],[have_umfpack4=no])
			LIBS=$save_LIBS
	
			if test $have_umfpack != no ; then
				UMFPACK_LIBS="$try_umfpack_LIBS"
				AC_DEFINE(USE_UMFPACK,1,[define if you want to use Umfpack])
				if test $have_umfpack4 != no ; then
					AC_DEFINE(HAVE_UMFPACK4,1,[define if Umfpack4 is available])
				fi
				ol_link_umfpack=yes
				if test $have_umfpack4 != no ; then
					AC_MSG_RESULT([(>= 4.0) yes])
				else
					AC_MSG_RESULT([yes])
				fi
				AC_DEFINE(NBLAS,1,[define if BLAS/ATLAS are not available to Umfpack])
				AC_MSG_WARN([Umfpack will be used without cblas/atlas])
			else
				AC_MSG_RESULT([no])
			fi
		fi
	fi

	if test $ol_link_umfpack = no ; then
		if test $ol_with_umfpack != auto ; then
			AC_MSG_ERROR([Could not locate Umfpack])
		else
			AC_MSG_WARN([Could not locate Umfpack])
			AC_MSG_WARN([Umfpack not supported!])
		fi
	fi
else
	AC_MSG_WARN([Umfpack not supported!])
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for Meschach
dnl
ol_link_meschach=no
if test $ol_with_meschach != no ; then
	OL_ADD_DIRS(meschach)
	
	AC_CHECK_HEADERS(meschach/sparse2.h)

	if test $ac_cv_header_meschach_sparse2_h = yes ; then

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
			AC_MSG_WARN([Meschach not supported!])
		fi
	fi
else
	AC_MSG_WARN([Meschach not supported!])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for Y12
dnl
ol_link_y12=no
if test $ol_with_y12 != no ; then
	ol_link_y12=yes
	AC_DEFINE(USE_Y12,1,[define to use Y12 sparse LU library])
	Y12_LIBS="../libraries/liby12/liby12.la"
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for Harwell
dnl
ol_link_harwell=no
if test $ol_with_harwell != no ; then
	AC_CHECK_LIB(hsl,ma28ad,[have_libhsl=yes],[have_libhsl=no],[${F2C_LIBS} -lm])
	if test $have_libhsl = no ; then
		AC_CHECK_LIB(hsl,ma28ad_,[have_libhsl=yes],[have_libhsl=no],[${F2C_LIBS} -lm])
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
dnl Check for a solver
dnl
if test $ol_link_meschach = no \
	-a $ol_link_y12 = no \
	-a $ol_link_harwell = no ; then
	AC_MSG_ERROR([Need at least one (sparse) solver])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for LAPACK
dnl
ol_link_lapack=no
if test $ol_with_lapack != no ; then
	OL_ADD_DIRS(lapack)

	AC_LANG_FORTRAN77
	AC_CHECK_LIB(lapack,dgesvd,[have_lapack=yes],[have_lapack=no],[-lblas])
	AC_LANG_C
	
	if test $have_lapack != no ; then
		AC_DEFINE(USE_LAPACK,1,[define to use LAPACK libraries])
		AC_DEFINE(HAVE_DGESVD,1,[define to use LAPACK's dgesvd in GPC])
		LAPACK_LIBS="-llapack -lblas "
	fi

	AC_LANG_FORTRAN77
	AC_CHECK_LIB(lapack,dgegv,[have_lapack=yes],[have_lapack=no],[-lblas])
	AC_LANG_C

	if test $have_lapack != no ; then
		ol_link_lapack=yes
		AC_DEFINE(HAVE_DGEGV,1,[define to use LAPACK's dgegv in Eig])
	fi
fi


dnl ----------------------------------------------------------------
dnl
dnl Check for PAM
dnl
ol_link_pam=no
if test $ol_with_pam != no ; then
	OL_ADD_DIRS(pam)
	AC_CHECK_HEADERS(security/pam_appl.h pam_misc.h security/pam_misc.h)
	
	if test $ac_cv_header_security_pam_appl_h = yes ; then
		if test $ac_cv_header_pam_misc_h = yes \
			-o $ac_cv_header_security_pam_misc_h = yes ; then
			AC_CHECK_LIB(pam,pam_start,
				[have_pam=yes],[have_pam=no],
				[-lcrypt -lpam_misc -ldl])
			
			if test $have_pam != no ; then
				SECURITY_LIBS="-lcrypt -lpam -lpam_misc -ldl"
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
			AC_MSG_WARN([Could not locate PAM])
			AC_MSG_WARN([PAM not supported!])
		fi		
	fi
else
	AC_MSG_WARN([PAM not supported!])
fi


dnl ----------------------------------------------------------------
dnl
dnl Check if crypt needs libcrypt
dnl
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
	AC_MSG_WARN([don't know how to link crypt...])
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for MPI & MPI++
dnl
ol_link_mpi=no
if test $ol_with_mpi != no ; then
	OL_ADD_DIRS(mpi)
	AC_CHECK_HEADER(mpi.h)
	AC_LANG_CPLUSPLUS
	AC_CHECK_HEADER(mpi++.h)
	AC_LANG_C
	
	if test $ac_cv_header_mpi_h = yes -a $ac_cv_header_mpipp_h = yes ; then
		AC_MSG_CHECKING([for mpi libraries])
		AC_LANG_CPLUSPLUS
		try_mpi_LIBS="-lmpich++ -lmpich -lnsl"
		save_LIBS=$LIBS
		LIBS="$LIBS $try_mpi_LIBS"
		AC_TRY_LINK([#include <mpi++.h>],[
			int argc = 0;
			char **argv = 0;
			MPI::Init(argc, argv);
		],[have_mpi=yes],[have_mpi=no])
		LIBS=$save_LIBS
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
			AC_MSG_WARN([Could not locate MPI])
			AC_MSG_WARN([MPI not supported!])
		fi
	fi
else
	AC_MSG_WARN([MPI not supported!])
fi


dnl ----------------------------------------------------------------
dnl
dnl Need Schur parallel solver with MPI
dnl
if test $ol_link_mpi = yes ; then
	AC_DEFINE(USE_SCHUR,1,[define if you want to use Schur parallel solver])
fi


dnl ----------------------------------------------------------------
dnl
dnl Need also Metis with MPI
dnl
ol_link_metis=no
if test $ol_with_metis != no ; then
	OL_ADD_DIRS(metis)
	
	AC_CHECK_HEADERS(metis.h)

	if test $ac_cv_header_metis_h = yes ; then
		AC_CHECK_LIB(metis, METIS_PartGraphVKway,
			[have_metis=yes],[have_metis=no])

		if test $have_metis != no ; then
			METIS_LIBS="-lmetis -lm"
			AC_DEFINE(USE_METIS,1,
				[define if you want to use Metis])
			ol_link_metis=yes
		fi
	fi

	if test $ol_link_metis = no ; then
		if test $ol_with_metis != auto ; then
			AC_MSG_ERROR([Could not locate Metis])
		else
			AC_MSG_WARN([Could not locate Metis])
			AC_MSG_WARN([Metis not supported!])
		fi

		if test $ol_link_mpi = yes ; then
			AC_MSG_ERROR([You NEED Metis \
				to use Schur parallel solver!])
		fi
	fi
else
	AC_MSG_WARN([Metis not supported!])
fi

dnl ----------------------------------------------------------------
dnl
dnl Check for RTAI
dnl
if test $ol_with_rtai != no ; then
	dnl FIXME: LOT TODO ...
	AC_CHECK_HEADERS(rtai_lxrt_user.h)
	AC_CHECK_HEADERS(net_rpc.h)

	AC_DEFINE(USE_RTAI,1,[define to use RTAI])
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks id DeadBeat can be built
dnl
if test $ol_link_lapack != no -o $ol_link_meschach != no ; then
	AC_DEFINE(USE_DBC,1,[define if you want to use GPC/DeadBeat control])
else
	if test $ol_with_elec != no ; then
		AC_MSG_WARN([GPC/Deadbeat not allowed!])
	fi
fi

dnl ----------------------------------------------------------------
dnl
dnl Use GiNaC?
dnl
ol_link_ginac=no
if test $ol_with_ginac != no ; then
	dnl AC_LANG_CPLUSPLUS
	dnl AM_PATH_GINAC(0.7.0,[ol_link_ginac=yes],[ol_link_ginac=no])
	dnl AC_LANG_C
	AC_MSG_CHECKING([for GiNaC library])
	which ginac-config >/dev/null 2>&1
	if test $? = 0 ; then
		GINACLIB_CONFIG="ginac-config"
		GINACLIB_LIBS=`$GINACLIB_CONFIG --libs`
		GINACLIB_CPPFLAGS=`$GINACLIB_CONFIG --cppflags`
		AC_DEFINE(HAVE_GINAC,1,[define if GiNaC is available])
		AC_MSG_RESULT([yes])
	else 
		if test $ol_with_ginac != yes ; then
			AC_MSG_WARN([GiNaC support not available])
		else
			AC_MSG_ERROR([need to have GiNaC installed])
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
	getcwd \
	gethostbyname \
	gethostbyname2 \
	getopt \
	getopt_long \
	inet_aton \
	localtime \
	memset \
	mlockall \
	nice \
	signal \
	strftime \
	strtod \
	strtol \
	strtoul \
	time \
])


dnl ----------------------------------------------------------------
dnl
dnl Checks whether ostream has form method
dnl
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether class ostream has method 'form'])
AC_TRY_COMPILE([#include <iostream.h>],[
	cout.form("%d", 1);
],[have_form=yes],[have_form=no])
AC_LANG_C

if test $have_form != no ; then
	AC_DEFINE(HAVE_FORM_IN_OSTREAM,1,
		[define if class ostream has 'form' method])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi

dnl ----------------------------------------------------------------
dnl
dnl Checks whether cmath defines std::isfinite
dnl
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([whether cmath defines 'std::isfinite'])
AC_TRY_LINK([#include <cmath>],[
	using std::isfinite;
	isfinite(1.);
],[have_std_isfinite=yes],[have_std_isfinite=no])
AC_LANG_C

if test $have_std_isfinite != no ; then
	AC_DEFINE(HAVE_STD_ISFINITE_IN_CMATH,1,
		[define if cmath defines 'std::isfinite'])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
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
AC_SUBST(F2C_LIBS)
AC_SUBST(MPI_LIBS)
AC_SUBST(METIS_LIBS)
AC_SUBST(UMFPACK_LIBS)
AC_SUBST(MESCHACH_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(SECURITY_LIBS)
AC_SUBST(TCL_LIBS)

AC_SUBST(Y12_LIBS)
AC_SUBST(HARWELL_LIBS)
AC_SUBST(MODULES_DIRS)
AC_SUBST(MODULE_LOADPATH)

AC_SUBST(GINACLIB_CONFIG)
AC_SUBST(GINACLIB_CPPFLAGS)
AC_SUBST(GINACLIB_LIBS)

AC_SUBST(DATE)


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
	libraries/Makefile:libraries/Makefile.in
	libraries/libmbutil/Makefile:libraries/libmbutil/Makefile.in
	libraries/libmbmath/Makefile:libraries/libmbmath/Makefile.in
	libraries/libmbwrap/Makefile:libraries/libmbwrap/Makefile.in
	libraries/liby12/Makefile:libraries/liby12/Makefile.in
	libraries/libobjs/Makefile:libraries/libobjs/Makefile.in
	mbdyn/Makefile:mbdyn/Makefile.in
	mbdyn/base/Makefile:mbdyn/base/Makefile.in
	mbdyn/aero/Makefile:mbdyn/aero/Makefile.in
	mbdyn/elec/Makefile:mbdyn/elec/Makefile.in
	mbdyn/hydr/Makefile:mbdyn/hydr/Makefile.in
	mbdyn/struct/Makefile:mbdyn/struct/Makefile.in
	modules/Makefile:modules/Makefile.in
	utils/Makefile:utils/Makefile.in
	etc/Makefile:etc/Makefile.in
	etc/modal.d/Makefile:etc/modal.d/Makefile.in
	etc/pam.d/Makefile:etc/pam.d/Makefile.in
	var/Makefile:var/Makefile.in
	manual/Makefile:manual/Makefile.in
	manual/mbdyn.1:manual/mbdyn.1.in
	include/Makefile:include/Makefile.in
	include/ac/Makefile:include/ac/Makefile.in
	build/Makefile:build/Makefile.in
],[
	date > stamp-h
])

