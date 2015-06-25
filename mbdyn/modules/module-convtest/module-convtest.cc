/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

/*
 * user-defined
 */
class ConvergenceTest
: virtual public Elem, public UserDefinedElem {
private:
	int numIter; 
	int curIter; 

public:
	ConvergenceTest(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ConvergenceTest(void);

	unsigned int iGetNumDof(void) const;
	DofOrder::Order SetDof(unsigned int i);
	void Output(OutputHandler& OH) const;
	std::ostream& Restart(std::ostream& out) const;
	void
	SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef, 
		const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
	void AfterPredict(VectorHandler& X, VectorHandler& XP);
	void Update( const VectorHandler& X, const VectorHandler& XP);
	unsigned int iGetInitialNumDof(void) const;
	void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);
	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

ConvergenceTest::ConvergenceTest(unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout("Module convtest" << std::endl);
	}

	numIter = HP.GetInt();
	curIter = 0;
}

ConvergenceTest::~ConvergenceTest(void)
{
	NO_OP;
}

unsigned int
ConvergenceTest::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
ConvergenceTest::SetDof(unsigned int i)
{
	return DofOrder::DIFFERENTIAL;
}

void
ConvergenceTest::Output(OutputHandler& OH) const
{
	return;
}

std::ostream&
ConvergenceTest::Restart(std::ostream& out) const
{
	return out << "not implemented yet;" << std::endl;
}

void 
ConvergenceTest::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
   	NO_OP;
}

void
ConvergenceTest::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
ConvergenceTest::AssJac(
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(iNumRows, iNumCols);

	integer iIndex = iGetFirstIndex() + 1;
	WM.PutRowIndex(1, iIndex);
	WM.PutColIndex(1, iIndex);

	/*
	 * set sub-matrix indices and coefs
	 */

	WM.PutCoef(1, 1, dCoef);
	
	return WorkMat;
}

SubVectorHandler& 
ConvergenceTest::AssRes(
		SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);
	WorkVec.PutRowIndex(1, iGetFirstIndex() + 1);

	/*
	 * set sub-vector indices and coefs
	 */

	doublereal d;

	if (curIter < numIter) {
		/* force an error */
		d = 1.;

	} else {
		/* force the exact solution */
		d = - XCurr(iGetFirstIndex() + 1);
	}

	WorkVec(1) = d;
	
	return WorkVec;
}

void
ConvergenceTest::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/* reset counter */
	curIter = 0;
}

void
ConvergenceTest::Update( const VectorHandler& X, const VectorHandler& XP)
{
	/* increment counter */
	curIter++;
}

unsigned int
ConvergenceTest::iGetInitialNumDof(void) const
{
	return 0;
}

void
ConvergenceTest::InitialWorkSpaceDim( integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;   
}

VariableSubMatrixHandler& 
ConvergenceTest::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	
	return WorkMat;
}

SubVectorHandler& 
ConvergenceTest::InitialAssRes(
	SubVectorHandler& WorkVec, 
	const VectorHandler& XCurr)
{
	WorkVec.Resize(0);
	
	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ConvergenceTest>;

	if (!SetUDE("convtest", rf)) {
		delete rf;

		silent_cerr("ConvergenceTest: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

