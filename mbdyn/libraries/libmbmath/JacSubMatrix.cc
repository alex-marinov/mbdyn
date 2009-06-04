/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003
 *
 * Marco Morandini	<morandini@aero.polimi.it>
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

/* here goes Morandini's copyright */

/* Portions Copyright (C) 2009 Pierangelo Masarati <masarati@aero.polimi.it> */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <iomanip>

#include "myassert.h"
#include "JacSubMatrix.h"

static ExpandableRowVector::ExpandableRow er_Zero;

ExpandableRowVector::ExpandableRowVector(void) {};
ExpandableRowVector::ExpandableRowVector(const integer n) {
	ReDim(n);
}
ExpandableRowVector::~ExpandableRowVector(void) {};
void ExpandableRowVector::ReDim(const integer n) {
	//we have to accept = 0, some elements do ReDim(0,0) (PointForceElement))
	ASSERTMSGBREAK(n>=0, "ExpandableRowVector:ReDim(), n shold be >= 0");
	v.resize(n, er_Zero);
}
void ExpandableRowVector::Zero(void) {
	for (std::vector<ExpandableRow>::iterator i = v.begin(); i != v.end(); i++) {
		i->x = 0.;
	}
}
void ExpandableRowVector::Reset(void) {
	std::fill(v.begin(), v.end(), er_Zero);
}
void ExpandableRowVector::Link(const integer i, const ExpandableRowVector*const xp) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Link() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::Link() overflow");
	ASSERTMSGBREAK(v[i - 1].idx == 0, "ExpandableRowVector::Link() fatal error");
	v[i - 1].xm = xp;
}
void ExpandableRowVector::Set(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Set() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::Set() overflow");
	v[i - 1].x = xx;
}
void ExpandableRowVector::SetIdx(integer i, integer iidx) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::SetIdx() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::SetIdx() overflow");
	v[i - 1].idx = iidx;
}
void ExpandableRowVector::Set(doublereal xx, integer i, integer iidx) {
	Set(xx, i);
	SetIdx(i, iidx);
}
doublereal&
ExpandableRowVector::operator ()(integer i)
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::() overflow");
	return v[i - 1].x;
}
const doublereal&
ExpandableRowVector::operator ()(integer i) const
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::() overflow");
	return v[i - 1].x;
}
void ExpandableRowVector::Add(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Add() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::Add() overflow");
	v[i - 1].x += xx;
}
void ExpandableRowVector::Sub(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Sub() underflow");
	ASSERTMSGBREAK(i <= v.size(), "ExpandableRowVector::Sub() overflow");
	v[i - 1].x -= xx;
}
void ExpandableRowVector::Add(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WorkVec.Add(v[i].idx, c*v[i].x);
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			v[i].xm->Add(WorkVec, c*v[i].x);
		}
	}
}
void ExpandableRowVector::Sub(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WorkVec.Sub(v[i].idx, c*v[i].x);
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
			v[i].xm->Sub(WorkVec, c*v[i].x);
		}
	}
}
void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
	const integer eq,
	const doublereal c) const {
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WM.IncCoef(eq, v[i].idx, c*v[i].x);
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			v[i].xm->Add(WM, eq, c*v[i].x);
		}
	}
}
void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
	const std::vector<integer>& eq,
	const std::vector<doublereal>& cc,
	const doublereal c) const
{
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		integer idx = v[i].idx;
		if (idx != 0) {
			doublereal d = c*v[i].x;
			for (std::vector<integer>::size_type j = 0; j < eq.size(); j++) {
				WM.IncCoef(eq[j], idx, cc[j]*d);
			}
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			v[i].xm->Add(WM, eq, cc, c*v[i].x);
		}
	}
}
void ExpandableRowVector::Sub(FullSubMatrixHandler& WM,
	const integer eq,
	const doublereal c) const {
	for (std::vector<ExpandableRow>::size_type i = 0 ; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WM.DecCoef(eq, v[i].idx, c*v[i].x);
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
			v[i].xm->Sub(WM, eq, c*v[i].x);
		}
	}
}
void ExpandableRowVector::Sub(FullSubMatrixHandler& WM, 
	const std::vector<integer>& eq,
	const std::vector<doublereal>& cc,
	const doublereal c) const
{
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		integer idx = v[i].idx;
		if (idx != 0) {
			doublereal d = c*v[i].x;
			for (std::vector<integer>::size_type j = 0; j < eq.size(); j++) {
				WM.DecCoef(eq[j], idx, cc[j]*d);
			}
		} else {
			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			v[i].xm->Sub(WM, eq, cc, c*v[i].x);
		}
	}
}

std::ostream & ExpandableRowVector::Write(std::ostream &out, const char *sFill) const {
	out << "LocalDof: ";
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		if (v[i].idx != 0) {
			out << sFill << std::setw(12) << v[i].idx;
		} else {
			out << sFill << std::setw(12) << "linked";
		}		
	}
	out << std::endl;
	out << "   Value: ";
	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
		out << sFill << std::setw(12) << v[i].x;
	}
	out << std::endl;
	return out;
}

std::ostream & operator << (std::ostream & s, const ExpandableRowVector & z) {
	return z.Write(s);	
}

