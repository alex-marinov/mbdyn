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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <algorithm>
#include <iomanip>

#include "myassert.h"
#include "JacSubMatrix.h"

static ExpandableRowVector::ExpandableRowElement er_Zero;

ExpandableRowVector::ExpandableRowVector(void) {};
ExpandableRowVector::ExpandableRowVector(const integer n) {
	ReDim(n);
}
ExpandableRowVector::~ExpandableRowVector(void) {};
void ExpandableRowVector::ReDim(const integer n) {
	//we have to accept = 0, some elements do ReDim(0,0) (PointForceElement))
	ASSERTMSGBREAK(n>=0, "ExpandableRowVector:ReDim(), n shold be >= 0");
// 	std::cerr << "ExpandableRowVector::ReDim("  << n << ")" << std::endl;
// 	std::cerr << "----------------------" << std::endl;
	v.resize(n, er_Zero);
}
integer ExpandableRowVector::GetDim(void) const {
	return v.size();
};
void ExpandableRowVector::Zero(void) {
	for (std::vector<ExpandableRowElement>::iterator i = v.begin(); i != v.end(); ++i) {
		i->x = 0.;
	}
}
void ExpandableRowVector::Reset(void) {
	std::fill(v.begin(), v.end(), er_Zero);
}
void ExpandableRowVector::Link(const integer i, const ExpandableRowVector*const xp, const integer rhs_block) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Link() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::Link() overflow");
	ASSERTMSGBREAK(v[i - 1].idx == 0, "ExpandableRowVector::Link() fatal error");
	if (std::vector<ExpandableRowVector *>::size_type(rhs_block) > v[ i - 1].xm.size()) {
		v[ i - 1].xm.resize(rhs_block);
	}
	v[i - 1].xm[rhs_block - 1] = xp;
}
void ExpandableRowVector::Set(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Set() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::Set() overflow");
	v[i - 1].x = xx;
}
void ExpandableRowVector::Set(Vec3 xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Set() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i+3) <= v.size(), "ExpandableRowVector::Set() overflow");
	v[i - 1].x = xx(1);
	v[i - 0].x = xx(2);
	v[i + 1].x = xx(3);
}
void ExpandableRowVector::SetIdx(integer i, integer iidx) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::SetIdx() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::SetIdx() overflow");
// 	std::cerr << "ExpandableRowVector::SetIdx(" << i <<", " <<  iidx << ")" << std::endl;
	v[i - 1].idx = iidx;
}
void ExpandableRowVector::Set(doublereal xx, integer i, integer iidx) {
	Set(xx, i);
	SetIdx(i, iidx);
}
void ExpandableRowVector::Set(Vec3 xx, integer i, integer iidx) {
	Set(xx, i);
	SetIdx(i, iidx);
	SetIdx(i+1, iidx+1);
	SetIdx(i+2, iidx+2);
}
doublereal&
ExpandableRowVector::operator ()(integer i)
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::() overflow");
	return v[i - 1].x;
}
const doublereal&
ExpandableRowVector::operator ()(integer i) const
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::() overflow");
	return v[i - 1].x;
}
void ExpandableRowVector::Add(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Add() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::Add() overflow");
	v[i - 1].x += xx;
}
void ExpandableRowVector::Sub(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Sub() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableRowVector>::size_type(i) <= v.size(), "ExpandableRowVector::Sub() overflow");
	v[i - 1].x -= xx;
}
void ExpandableRowVector::Add(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WorkVec.IncCoef(v[i].idx, c*v[i].x);
		} else {
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Add(WorkVec, c*v[i].x);
			}	
		}
	}
}
void ExpandableRowVector::Sub(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WorkVec.DecCoef(v[i].idx, c*v[i].x);
		} else {
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Sub(WorkVec, c*v[i].x);
			}
		}
	}
}
void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
	const integer eq,
	const doublereal c) const {
// 	std::cerr << "\t\t\tExpandableRowVector::Add" << std::endl;
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
// 		std::cerr << "\t\t\t\ti " << i << std::endl;

		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
// 			std::cerr << "Adding("<< eq << ", " << v[i].idx << ") = " <<
// 			c*v[i].x << std::endl;
// 			std::cerr << "\t\t\t\tdentro if" << std::endl;
			WM.IncCoef(eq, v[i].idx, c*v[i].x);
		} else {
// 			std::cerr << "\t\t\t\tdentro else" << std::endl;
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
// 				std::cerr << "\t\t\t\t\trhs_block " << rhs_block << std::endl;
// 				std::cerr << "\t\t\t\t\t\trhs_block_ptr " << v[i].xm[rhs_block] << std::endl;
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Add(WM, eq, c*v[i].x);
			}
		}
	}
// 	std::cerr << "\t\t\t ---- ExpandableRowVector::Add" << std::endl;
}
void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
	const std::vector<integer>& eq,
	const std::vector<doublereal>& cc,
	const doublereal c) const
{
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
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
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Add(WM, eq, cc, c*v[i].x);
			}
		}
	}
}
void ExpandableRowVector::Sub(FullSubMatrixHandler& WM,
	const integer eq,
	const doublereal c) const {
	for (std::vector<ExpandableRowElement>::size_type i = 0 ; i < v.size(); i++) {
		if (v[i].x == 0.) {
			continue;
		}

		if (v[i].idx != 0) {
			WM.DecCoef(eq, v[i].idx, c*v[i].x);
		} else {
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Sub(WM, eq, c*v[i].x);
			}
		}
	}
}
void ExpandableRowVector::Sub(FullSubMatrixHandler& WM, 
	const std::vector<integer>& eq,
	const std::vector<doublereal>& cc,
	const doublereal c) const
{
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
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
			for (std::vector<const ExpandableRowVector *>::size_type rhs_block = 0; rhs_block < v[i].xm.size(); rhs_block++) {
				ASSERTMSGBREAK(v[i].xm[rhs_block] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
				v[i].xm[rhs_block]->Sub(WM, eq, cc, c*v[i].x);
			}
		}
	}
}

std::ostream & ExpandableRowVector::Write(std::ostream &out, const char *sFill) const {
	out << "LocalDof: ";
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
		if (v[i].idx != 0) {
			out << sFill << std::setw(12) << v[i].idx;
		} else {
			out << sFill << std::setw(12) << "linked";
		}		
	}
	out << std::endl;
	out << "   Value: ";
	for (std::vector<ExpandableRowElement>::size_type i = 0; i < v.size(); i++) {
		out << sFill << std::setw(12) << v[i].x;
	}
	out << std::endl;
	return out;
}

std::ostream & operator << (std::ostream & s, const ExpandableRowVector & z) {
	return z.Write(s);	
}

// -------------------------------------
static ExpandableMatrix::ExpandableColBlock ecb_Zero;


ExpandableMatrix::ExpandableMatrix(void) {};
ExpandableMatrix::ExpandableMatrix(const integer n, const integer m) {
	ReDim(n, m);
}
ExpandableMatrix::~ExpandableMatrix(void) {};
void ExpandableMatrix::ReDim(const integer n, const integer m) {
	//we have to accept = 0, some elements do ReDim(0,0) (PointForceElement))
	ASSERTMSGBREAK(n>=0, "ExpandableMatrix:ReDim(), n shold be >= 0");
	ASSERTMSGBREAK(m>=0, "ExpandableMatrix:ReDim(), m shold be >= 0");
	v.resize(m, ecb_Zero);
	for (std::vector<ExpandableColBlock>::iterator i = v.begin(); i != v.end(); ++i) {
		i->ReDim(n);
	}
}
void ExpandableMatrix::SetBlockDim(const integer block, const integer ncols) {
	ASSERTMSGBREAK(block>0, "ExpandableMatrix:SetBlockDim(), block shold be > 0");
	ASSERTMSGBREAK(block <= GetNBlocks(), "ExpandableMatrix:SetBlockDim(), block shold be "
		"<= GetNBlocks()");
// 	std::cerr << "ExpandableMatrix::SetBlockDim(" << block <<", " << ncols << ")" << std::endl;
	int pippo = 0;
	for (std::vector<ExpandableRowVector>::iterator i = v[block-1].rows.begin();
		i != v[block-1].rows.end(); ++i)
	{
			pippo++;
// 			std::cerr << "\t" << pippo << std::endl;
			i->ReDim(ncols);
	}
}
void ExpandableMatrix::Zero(void) {
	for (std::vector<ExpandableColBlock>::iterator i = v.begin(); i != v.end(); ++i) {
		i->Zero();
	}
}
integer ExpandableMatrix::GetNRows() const {
	if (v.begin() != v.end() ) {
		return v[0].GetBlockNRows(); 
	} else {
		return 0;
	}
};
integer ExpandableMatrix::GetNBlocks() const {
	return v.size();
};
integer ExpandableMatrix::GetBlockNCols(const integer block) const {
	ASSERTMSGBREAK(block >= 0, "ExpandableMatrix:GetBlockNCols(), "
		"block shold be >= 0");
	ASSERTMSGBREAK((unsigned long)block < v.size(), "ExpandableMatrix:GetBlockNCols(), "
		"block shold be < nblocks");
	return v[block].GetBlockNCols();
}
void ExpandableMatrix::Reset(void) {
	for (std::vector<ExpandableColBlock>::iterator i = v.begin();
		i != v.end(); ++i)
	{
		i->Reset();
	}
}

void ExpandableMatrix::Link(const integer i, const ExpandableMatrix*const xp) {
	ASSERTMSGBREAK(i > 0, "ExpandableMatrix::Link() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableColBlock>::size_type(i) <= v.size(), "ExpandableMatrix::Link() overflow");
	//FIXME
	//ASSERTMSGBREAK(v[i - 1].idx == 0, "ExpandableMatrix::Link() fatal error");
	ASSERTMSGBREAK(v[i - 1].GetBlockNCols() == xp->GetNRows(), 
		"ExpandableMatrix::Link() dimension mismatch");
	v[i - 1].Link(xp);
}

void ExpandableMatrix::Link(const integer i, const ExpandableRowVector*const xp) {
	ASSERTMSGBREAK(i > 0, "ExpandableMatrix::Link() underflow");
	ASSERTMSGBREAK(std::vector<ExpandableColBlock>::size_type(i) <= v.size(), "ExpandableMatrix::Link() overflow");
	//FIXME
	//ASSERTMSGBREAK(v[i - 1].idx == 0, "ExpandableMatrix::Link() fatal error");
	ASSERTMSGBREAK(v[i - 1].GetBlockNCols() == 1, 
		"ExpandableMatrix::Link() dimension mismatch");
	v[i - 1].Link(xp);
}

void ExpandableMatrix::Set(doublereal xx, integer eq, integer block, integer block_col) {
	ASSERTMSGBREAK(eq > 0, "ExpandableMatrix::Set() underflow");
	ASSERTMSGBREAK(block > 0, "ExpandableMatrix::Set() underflow");
	ASSERTMSGBREAK(block_col > 0, "ExpandableMatrix::Set() underflow");
	ASSERTMSGBREAK(eq <= GetNRows(), "ExpandableRowVector::Set() overflow");
	ASSERTMSGBREAK(block <= GetNBlocks(), "ExpandableRowVector::Set() overflow");
	ASSERTMSGBREAK(block_col <= GetBlockNCols(block), "ExpandableRowVector::Set() overflow");
	v[block - 1].rows[eq - 1].Set(xx, block_col);
}
void ExpandableMatrix::Set(Vec3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		Set(xx(i + 1), eq + i, block, block_col);
	}
}
void ExpandableMatrix::Set(Mat3x3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		for (integer ii = 0; ii <3; ii++) {
			Set(xx(i + 1, ii + 1), eq + i, block, block_col + ii);
		}
	}
}
void ExpandableMatrix::SetBlockIdx(integer block, integer iidx) {
	ASSERTMSGBREAK(block > 0, "ExpandableMatrix::SetColIdx() underflow");
	ASSERTMSGBREAK(block <= GetNBlocks(), "ExpandableMatrix::SetColIdx() overflow");
// 	std::cerr << "ExpandableMatrix::SetBlockIdx" << std::endl;
	v[block - 1].SetColIdx(iidx);
}
// void ExpandableRowVector::Set(doublereal xx, integer block, integer i, integer iidx) {
// 	Set(xx, i);
// 	SetIdx(i, iidx);
// }
// doublereal&
// ExpandableRowVector::operator ()(integer i)
// {
// 	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
// 	ASSERTMSGBREAK(std::vector<ExpandableRow>::size_type(i) <= v.size(), "ExpandableRowVector::() overflow");
// 	return v[i - 1].x;
// }
// const doublereal&
// ExpandableRowVector::operator ()(integer i) const
// {
// 	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
// 	ASSERTMSGBREAK(std::vector<ExpandableRow>::size_type(i) <= v.size(), "ExpandableRowVector::() overflow");
// 	return v[i - 1].x;
// }
void ExpandableMatrix::Add(doublereal xx, integer eq, integer block, integer block_col) {
	ASSERTMSGBREAK(eq > 0, "ExpandableMatrix::Add() underflow");
	ASSERTMSGBREAK(block > 0, "ExpandableMatrix::Add() underflow");
	ASSERTMSGBREAK(block_col > 0, "ExpandableMatrix::Add() underflow");
	ASSERTMSGBREAK(eq <= GetNRows(), "ExpandableRowVector::Add() overflow");
	ASSERTMSGBREAK(block <= GetNBlocks(), "ExpandableRowVector::Add() overflow");
	ASSERTMSGBREAK(block_col <= GetBlockNCols(block), "ExpandableRowVector::Add() overflow");
	v[block - 1].rows[eq - 1].Add(xx, block_col);
}
void ExpandableMatrix::Sub(doublereal xx, integer eq, integer block, integer block_col) {
	ASSERTMSGBREAK(eq > 0, "ExpandableMatrix::Sub() underflow");
	ASSERTMSGBREAK(block > 0, "ExpandableMatrix::Sub() underflow");
	ASSERTMSGBREAK(block_col > 0, "ExpandableMatrix::Sub() underflow");
	ASSERTMSGBREAK(eq <= GetNRows(), "ExpandableRowVector::Sub() overflow");
	ASSERTMSGBREAK(block <= GetNBlocks(), "ExpandableRowVector::Sub() overflow");
	ASSERTMSGBREAK(block_col <= GetBlockNCols(block), "ExpandableRowVector::Sub() overflow");
	v[block - 1].rows[eq - 1].Sub(xx, block_col);
}
void ExpandableMatrix::Add(Vec3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		Add(xx(i + 1), eq + i, block, block_col);
	}
}
void ExpandableMatrix::Add(Mat3x3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		for (integer ii = 0; ii <3; ii++) {
			Add(xx(i + 1, ii + 1), eq + i, block, block_col + ii);
		}
	}
}
void ExpandableMatrix::Sub(Vec3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		Sub(xx(i + 1), eq + i, block, block_col);
	}
}
void ExpandableMatrix::Sub(Mat3x3 xx, integer eq, integer block, integer block_col) {
	for (integer i = 0; i <3; i++) {
		for (integer ii = 0; ii <3; ii++) {
			Sub(xx(i + 1, ii + 1), eq + i, block, block_col + ii);
		}
	}
}
// void ExpandableRowVector::Add(SubVectorHandler& WorkVec, const doublereal c) const {
// 	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
// 		if (v[i].x == 0.) {
// 			continue;
// 		}
// 
// 		if (v[i].idx != 0) {
// 			WorkVec.Add(v[i].idx, c*v[i].x);
// 		} else {
// 			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
// 			v[i].xm->Add(WorkVec, c*v[i].x);
// 		}
// 	}
// }
// void ExpandableRowVector::Sub(SubVectorHandler& WorkVec, const doublereal c) const {
// 	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
// 		if (v[i].x == 0.) {
// 			continue;
// 		}
// 
// 		if (v[i].idx != 0) {
// 			WorkVec.Sub(v[i].idx, c*v[i].x);
// 		} else {
// 			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
// 			v[i].xm->Sub(WorkVec, c*v[i].x);
// 		}
// 	}
// }
void ExpandableMatrix::Add(FullSubMatrixHandler& WM,
	const integer eq,
	const doublereal c) const {
// 	std::cerr << "ExpandableMatrix::Add" << std::endl;
// 		std::cerr << "\teq " << eq << std::endl;
// 		std::cerr << "\tc " << c << std::endl;
	for (integer block = 0; block < GetNBlocks(); block++) {
// 		std::cerr << "\t\tblock " << block << std::endl;
		v[block].Add(WM, eq, c);
	}
}
// void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
// 	const std::vector<integer>& eq,
// 	const std::vector<doublereal>& cc,
// 	const doublereal c) const
// {
// 	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
// 		if (v[i].x == 0.) {
// 			continue;
// 		}
// 
// 		integer idx = v[i].idx;
// 		if (idx != 0) {
// 			doublereal d = c*v[i].x;
// 			for (std::vector<integer>::size_type j = 0; j < eq.size(); j++) {
// 				WM.IncCoef(eq[j], idx, cc[j]*d);
// 			}
// 		} else {
// 			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
// 			v[i].xm->Add(WM, eq, cc, c*v[i].x);
// 		}
// 	}
// }
void ExpandableMatrix::Sub(FullSubMatrixHandler& WM,
	const integer eq,
	const doublereal c) const {
	for (integer block = 0; block < GetNBlocks(); block++) {
		v[block].Sub(WM, eq, c);
	}
}
// void ExpandableRowVector::Sub(FullSubMatrixHandler& WM, 
// 	const std::vector<integer>& eq,
// 	const std::vector<doublereal>& cc,
// 	const doublereal c) const
// {
// 	for (std::vector<ExpandableRow>::size_type i = 0; i < v.size(); i++) {
// 		if (v[i].x == 0.) {
// 			continue;
// 		}
// 
// 		integer idx = v[i].idx;
// 		if (idx != 0) {
// 			doublereal d = c*v[i].x;
// 			for (std::vector<integer>::size_type j = 0; j < eq.size(); j++) {
// 				WM.DecCoef(eq[j], idx, cc[j]*d);
// 			}
// 		} else {
// 			ASSERTMSGBREAK(v[i].xm != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
// 			v[i].xm->Sub(WM, eq, cc, c*v[i].x);
// 		}
// 	}
// }

std::ostream & ExpandableMatrix::Write(std::ostream &out, const char *sFill) const {
	for (integer block = 0; block < GetNBlocks(); block++) {
		out << "Block: " <<  block << std::endl;
		v[block].Write(out, "\t");
	}
	return out;
}

std::ostream & operator << (std::ostream & s, const ExpandableMatrix & z) {
	return z.Write(s);	
}
