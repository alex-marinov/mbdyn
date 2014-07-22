dnl $Header$
dnl MBDyn (C) is a multibody analysis code.
dnl http://www.mbdyn.org
dnl
dnl Copyright (C) 1996-2014
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
dnl *** Borrowed from:
dnl
dnl Copyright 1998-2004 The OpenLDAP Foundation,  All Rights Reserved.
dnl COPYING RESTRICTIONS APPLY, See COPYRIGHT file
dnl
dnl OpenLDAP Autoconf Macros
dnl
divert(-1)
builtin(include, build/mbdyn.m4)dnl
builtin(include, build/m4_ax_pkg_swig.m4)dnl
builtin(include, build/m4_ax_python_devel.m4)dnl
builtin(include, build/m4_ax_swig_enable_cxx.m4)dnl
builtin(include, build/m4_ax_swig_multi_module_support.m4)dnl
builtin(include, build/m4_ax_swig_python.m4)dnl
