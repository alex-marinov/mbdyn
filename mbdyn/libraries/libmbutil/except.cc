/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <sstream>

#include "except.h"

MBDynErrBase::MBDynErrBase(MBDYN_EXCEPT_ARGS_DECL_NODEF)
{
	std::stringstream ss;
	ss << "[" << file << ":" << line << ",func=" << func << "]";
	if (!r.empty()) {
		ss << " (" << r << ")";
	}
	s = ss.str();
}

void
MBDynErrBase::Set(const std::string& s)
{
	this->s = s;
}

const char *
MBDynErrBase::what(void) const throw()
{
	return s.c_str();
}

void
ErrIndexOutOfRange::WriteMsg(const char *idx_type, int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL_NODEF)
{
	std::stringstream ss;
	ss << "[" << file << ":" << line << ",func=" << func << "]";
	if (!r.empty()) {
		ss << " (" << r << ")";
	}
	ss << ": " << idx_type << "index=" << idx << " out of range (" << imin << ":" << imax << ")";
	Set(ss.str());
}

ErrIndexOutOfRange::ErrIndexOutOfRange(const char *idx_type, int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL_NODEF)
: ErrOutOfRange(MBDYN_EXCEPT_ARGS_PASSTHRU)
{
	WriteMsg(idx_type, idx, imin, imax, MBDYN_EXCEPT_ARGS);
}

ErrIndexOutOfRange::ErrIndexOutOfRange(int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL_NODEF)
: ErrOutOfRange(MBDYN_EXCEPT_ARGS_PASSTHRU)
{
	WriteMsg("", idx, imin, imax, MBDYN_EXCEPT_ARGS);
}

