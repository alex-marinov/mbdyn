# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2015
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
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
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 

%module mbc_py
%include "typemaps.i"
%header
%{
        #include <numpy/arrayobject.h>
        #include "mbc_py_global.h"

        extern int
        mbc_py_nodal_initialize(const char *const path,
                const char *const host, unsigned port,
                int timeout, unsigned verbose, unsigned data_and_next,
                unsigned refnode, unsigned nodes,
                unsigned labels, unsigned rot, unsigned accels);

        extern int
        mbc_py_nodal_negotiate(unsigned id);

        extern int
        mbc_py_nodal_send(unsigned id, int last);

        extern int
        mbc_py_nodal_recv(unsigned id);

        extern int
        mbc_py_nodal_destroy(unsigned id);

        extern int
        mbc_py_modal_initialize(const char *const path,
                const char *const host, unsigned port,
                int timeout, unsigned verbose, unsigned data_and_next,
                unsigned refnode, unsigned modes);

        extern int
        mbc_py_modal_negotiate(unsigned id);

        extern int
        mbc_py_modal_send(unsigned id, int last);

        extern int
        mbc_py_modal_recv(unsigned id);

        extern int
        mbc_py_modal_destroy(unsigned id);
%}

%include "mbc_py_global.h"
int
mbc_py_nodal_initialize(const char *path,
        const char *host, unsigned port,
        int timeout, unsigned verbose, unsigned data_and_next,
        unsigned refnode, unsigned nodes,
        unsigned labels, unsigned rot, unsigned accels);

int
mbc_py_nodal_negotiate(unsigned id);

int
mbc_py_nodal_send(unsigned id, int last);

int
mbc_py_nodal_recv(unsigned id);

int
mbc_py_nodal_destroy(unsigned id);

int
mbc_py_modal_initialize(const char *const path,
        const char *const host, unsigned port,
        int timeout, unsigned verbose, unsigned data_and_next,
        unsigned refnode, unsigned modes);

int
mbc_py_modal_negotiate(unsigned id);

int
mbc_py_modal_send(unsigned id, int last);

int
mbc_py_modal_recv(unsigned id);

int
mbc_py_modal_destroy(unsigned id);

%init
%{
        import_array();
%}

