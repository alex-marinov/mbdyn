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


#ifndef JacSubMatrix_hh
#define JacSubMatrix_hh

#include <vector>

#include "ac/f2c.h"
#include "submat.h"

class ExpandableRowVector {
private:
	std::vector<doublereal> x;
	std::vector<const ExpandableRowVector*> xm;
	std::vector<integer> idx;
	ExpandableRowVector & operator = (const ExpandableRowVector &); // not to be implemented
	ExpandableRowVector (const ExpandableRowVector &); // not to be implemented
public:
	ExpandableRowVector();
	ExpandableRowVector(const integer n);
	virtual ~ExpandableRowVector();
	void ReDim(const integer n);
	void Zero();
	void Reset();
	void Link(const integer i, const ExpandableRowVector* const xp);
	void Set(doublereal xx, integer i);
	void SetIdx(integer i, integer iidx);
	void Set(doublereal xx, integer i, integer iidx);
	void Add(doublereal xx, integer i);
	void Sub(doublereal xx, integer i);
	void Add(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
	void Sub(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
	void Add(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
	void Sub(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
	std::ostream & Write(std::ostream &out, 
		const char *sFill="") const;
};

std::ostream & operator << (std::ostream & s, const ExpandableRowVector & z);

#endif /* JacSubMatrix_hh */

