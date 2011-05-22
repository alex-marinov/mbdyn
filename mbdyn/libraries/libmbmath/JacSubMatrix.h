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

/* Portions Copyright (C) 2009 Pierangelo Masarati <masarati@aero.polimi.it> */


#ifndef JacSubMatrix_hh
#define JacSubMatrix_hh

#include <vector>

#include "ac/f2c.h"
#include "submat.h"

class ExpandableRowVector {
public:
	struct ExpandableRowElement {
		doublereal x;
		std::vector<const ExpandableRowVector *> xm;
		integer idx;
		ExpandableRowElement(void) : x(0.), idx(0) {xm.resize(1); x = 0;};
	};
private:
	std::vector<ExpandableRowElement> v;

// 	ExpandableRowVector & operator = (const ExpandableRowVector &); // not to be implemented
// 	ExpandableRowVector (const ExpandableRowVector &); // not to be implemented

public:
	ExpandableRowVector(void);
	ExpandableRowVector(const integer n);
	virtual ~ExpandableRowVector(void);
	void ReDim(const integer n);
	integer GetDim(void) const;
	void Zero(void);
	void Reset(void);
	void Link(const integer i, const ExpandableRowVector* const xp, const integer rhs_block = 1);
	void SetIdx(integer i, integer iidx);
	void Set(doublereal xx, integer i, integer iidx);
	void Set(Vec3 xx, integer i, integer iidx);
	doublereal& operator ()(integer i);
	const doublereal& operator ()(integer i) const;
	void Set(doublereal xx, integer i);
	void Set(Vec3 xx, integer i);
	void Add(doublereal xx, integer i);
	void Sub(doublereal xx, integer i);
	void Add(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
	void Sub(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
	void Add(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
	void Add(FullSubMatrixHandler& WM, const std::vector<integer>& eq,
		const std::vector<doublereal>& cc, const doublereal c = 1.) const;
	void Sub(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
	void Sub(FullSubMatrixHandler& WM, const std::vector<integer>& eq,
		const std::vector<doublereal>& cc, const doublereal c = 1.) const;
	std::ostream & Write(std::ostream &out, const char *sFill = "") const;
};

std::ostream & operator << (std::ostream & s, const ExpandableRowVector & z);

class ExpandableMatrix {
public:
	struct ExpandableColBlock;
private:
	std::vector<ExpandableColBlock> v;

	ExpandableMatrix & operator = (const ExpandableMatrix &); // not to be implemented
	ExpandableMatrix (const ExpandableMatrix &); // not to be implemented

public:
	struct ExpandableColBlock {
		std::vector<ExpandableRowVector> rows;
// 		const ExpandableMatrix *xm;
// 		integer ncols;
// 		integer idx;
		ExpandableColBlock(void) : rows(0) {}; //, idx(0) {};
		void Zero(void) {
			for (std::vector<ExpandableRowVector>::iterator i = rows.begin();
				i != rows.end(); ++i)
			{
				i->Zero();
			}
		};
		void ReDim(const integer nrows) {
			rows.resize(nrows);
		};
		void SetBlockDim(const integer nrows, const integer ncols) {
			for (std::vector<ExpandableRowVector>::iterator i = rows.begin();
				i != rows.end(); ++i)
			{
				i->ReDim(ncols);
			}
		};
		integer GetBlockNRows() const {
			return rows.size();
		};
		integer GetBlockNCols() const {
			if (rows.begin() != rows.end()) {
				return rows[0].GetDim();
			} else {
				return 0;
			}
		};
		void Reset() {
			for (std::vector<ExpandableRowVector>::iterator i = rows.begin();
				i != rows.end(); ++i)
			{
				i->Reset();
			}
		};
// 		void Link(ExpandableRowVector & row, std::vector<ExpandableRowVector>& m1) {
// 			ASSERTMSGBREAK(row.GetDim() == m1.size(), 
// 				"ExpandableMatrix::ExpandableColBlock::Link() dimension mismatch");
// 			for (integer i = 0; i < m1.size(); i++) {
// 				v1.Link(i, m1[i]);
// 			}
// 		}
		void Link(const ExpandableMatrix *const xp) {
// 			xm = xp;
// 			integer rhs_block = 0;
// 			integer rhs_block_nrows = xp->GetNRows();
// 			std::cerr << "ExpandableColBlock::Link" << std::endl;
			for (integer col = 1; col <= GetBlockNCols(); col++) {
// 				std::cerr << "\t col " << col << std::endl;
// 				if (col >= rhs_block_ncols) {
// 					rhs_block++;
// 					rhs_block_ncols += xp->v[rhs_block].GetBlockNCols();
// 				}
				for (std::vector<ExpandableRowVector>::size_type eq = 0; eq < rows.size(); eq++) {
// 					std::cerr << "\t\t eq " << eq << std::endl;
					for (integer rhs_block = 1; rhs_block <=
						xp->GetNBlocks(); rhs_block++) {
// 						std::cerr << "\t\t\t rhs_block " << rhs_block << std::endl;
						rows[eq].Link(col, &(xp->v[rhs_block - 1].rows[col - 1]), rhs_block);
					}
				}
			}
// 			std::cerr << "-------------------" << std::endl;
		};
		void Link(const ExpandableRowVector *const xp) {
// 			xm = xp;
			for (std::vector<ExpandableRowVector>::size_type eq = 0; eq < rows.size(); eq++) {
				rows[eq].Link(1, xp);
			}
		};
		void SetColIdx(integer iidx) {
// 			std::cerr << "\tBloc::SetColIdx" << std::endl;
			for (integer col = 1; col <= GetBlockNCols(); col++, iidx++) {
// 				std::cerr << "\t\tcol " << col << " iidx " << iidx << std::endl;
				for (std::vector<ExpandableRowVector>::size_type eq = 0; eq < rows.size(); eq++) {
// 					std::cerr << "\t\t\teq " << eq <<std::endl;
					rows[eq].SetIdx(col, iidx);
				}
			}
		};
		void Add(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const {
// 			std::cerr << "\t\tBlock::Add " << std::endl;
			for (integer e = 0; e < GetBlockNRows(); e++) {
// 				std::cerr << "\t\t\te " << e << " eq+e " << eq+e << std::endl;
				rows[e].Add(WM, eq + e, c);
			}
		};
		void Sub(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const {
			for (integer e = 0; e < GetBlockNRows(); e++) {
				rows[e].Sub(WM, eq + e, c);
			}
		};
		std::ostream & Write(std::ostream &out, const char *sFill) const {
			for (integer eq = 0; eq < GetBlockNRows(); eq++) {
				out << sFill << "Row: " << eq << std::endl;
				rows[eq].Write(out, "\t\t");
			}
			return out;
		}
	};
public:
	ExpandableMatrix(void);
	ExpandableMatrix(const integer n, const integer m);
	virtual ~ExpandableMatrix(void);
	void ReDim(const integer n, const integer m);
	void SetBlockDim(const integer block, const integer ncols);
	void Zero(void);
	integer GetNRows() const;
	integer GetNBlocks() const;
	integer GetBlockNCols(const integer block) const;
	
	void Reset(void);
	void Link(const integer block, const ExpandableMatrix* const xp);
	void Link(const integer block, const ExpandableRowVector* const xp);
	void SetBlockIdx(integer block, integer iidx);

// 	void Set(doublereal xx, integer block, integer eq, integer iidx);
// 	doublereal& operator ()(integer i);
// 	const doublereal& operator ()(integer i) const;
	void Set(doublereal xx, integer eq, integer block, integer block_col = 1);
	void Set(Vec3 xx, integer eq, integer block, integer block_col = 1);
	void Set(Mat3x3 xx, integer eq, integer block, integer block_col = 1);
	void Add(doublereal xx, integer eq, integer block, integer block_col = 1);
	void Sub(doublereal xx, integer eq, integer block, integer block_col = 1);
	void Add(Vec3 xx, integer eq, integer block, integer block_col = 1);
	void Sub(Vec3 xx, integer eq, integer block, integer block_col = 1);
	void Add(Mat3x3 xx, integer eq, integer block, integer block_col = 1);
	void Sub(Mat3x3 xx, integer eq, integer block, integer block_col = 1);
// 	void Add(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
// 	void Sub(SubVectorHandler& WorkVec, const doublereal c = 1.) const;
	void Add(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
// 	void Add(FullSubMatrixHandler& WM, const std::vector<integer>& eq,
// 		const std::vector<doublereal>& cc, const doublereal c = 1.) const;
	void Sub(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.) const;
// 	void Sub(FullSubMatrixHandler& WM, const std::vector<integer>& eq,
// 		const std::vector<doublereal>& cc, const doublereal c = 1.) const;
	std::ostream & Write(std::ostream &out, const char *sFill = "") const;
};


#endif /* JacSubMatrix_hh */

