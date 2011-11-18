/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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

/* joint regularization */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "jointreg.h"
#include "dataman.h"
#include "driven.h"

/* JointRegularization - begin */

JointRegularization::JointRegularization(unsigned int uL,
	const Joint *j,
	flag fOut)
: Elem(uL, fOut),
InitialAssemblyElem(uL, fOut),
pJ(j)
{
	NO_OP;
}

JointRegularization::~JointRegularization(void)
{
	NO_OP;
}

unsigned int
JointRegularization::iGetInitialNumDof(void) const
{
	return 0;
}

/* JointRegularization - end */


/* TikhonovRegularization - begin */

TikhonovRegularization::TikhonovRegularization(unsigned int uL,
	const Joint *j,
	const std::vector<doublereal>& c,
	flag fOut)
: Elem(uL, fOut),
JointRegularization(uL, j, fOut),
dC(c)
{
	NO_OP;
}

TikhonovRegularization::~TikhonovRegularization(void)
{
	NO_OP;
}

JointRegularization::Type
TikhonovRegularization::GetJointRegularizationType(void) const
{
	return JointRegularization::TIKHONOV_REGULARIZATION;
}


void
TikhonovRegularization::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetNumDof();
}

/* assemblaggio residuo */
SubVectorHandler&
TikhonovRegularization::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WorkVec.ResizeReset(iNumDofs);

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[iCnt - 1]*XCurr(iFirstIndex + iCnt)/dCoef);
	}

	return WorkVec;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
TikhonovRegularization::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WM.ResizeReset(iNumDofs, 0);

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[iCnt - 1]/dCoef);
	}

	return WorkMat;
}

/* inverse dynamics capable element */
bool
TikhonovRegularization::bInverseDynamics(void) const
{
	return true;
}

/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
TikhonovRegularization::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	return AssJac(WorkMat, 1., XCurr, XCurr);
}

void
TikhonovRegularization::InitialWorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetInitialNumDof();
}

VariableSubMatrixHandler& 
TikhonovRegularization::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WM.ResizeReset(iNumDofs, 0);
	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[tmp2]);
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}
	}



	return WorkMat;
}

SubVectorHandler& 
TikhonovRegularization::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{	
	
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WorkVec.ResizeReset(iNumDofs);

	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[tmp2]*XCurr(iFirstIndex + iCnt));
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}
	}


	return WorkVec;
}

/* TikhonovRegularization - end */


/* DynamicRegularization - begin */

DynamicRegularization::DynamicRegularization(unsigned int uL,
	const Joint *j,
	const std::vector<doublereal>& c,
	flag fOut)
: Elem(uL, fOut),
JointRegularization(uL, j, fOut),
dC(c),
dLambda(c.size())
{
	NO_OP;
}

DynamicRegularization::~DynamicRegularization(void)
{
	NO_OP;
}

JointRegularization::Type
DynamicRegularization::GetJointRegularizationType(void) const
{
	return JointRegularization::DYNAMIC_REGULARIZATION;
}

void
DynamicRegularization::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetNumDof();
}

void
DynamicRegularization::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	for (std::vector<doublereal>::iterator i = dLambda.begin();
		i != dLambda.end(); ++i )
	{
		*i = 0.;
	}
}

/* assemblaggio residuo */
SubVectorHandler&
DynamicRegularization::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WorkVec.ResizeReset(iNumDofs);

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		unsigned iCntm1 = iCnt - 1;
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[iCntm1]*dLambda[iCntm1]/dCoef);
		dLambda[iCntm1] = XCurr(iFirstIndex + iCnt);
	}

	return WorkVec;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
DynamicRegularization::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WM.ResizeReset(iNumDofs, 0);
	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[iCnt - 1]/dCoef);
	}

	return WorkMat;
}


void
DynamicRegularization::InitialWorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetInitialNumDof();
}

VariableSubMatrixHandler& 
DynamicRegularization::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{		
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WM.ResizeReset(iNumDofs, 0);
	
	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[tmp2]);
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}
	}
	
	return WorkMat;
}

SubVectorHandler& 
DynamicRegularization::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WorkVec.ResizeReset(iNumDofs);

	std::vector<doublereal> dInitialLambda(iNumDofs);

	for (std::vector<doublereal>::iterator i = dInitialLambda.begin();
		i != dInitialLambda.end(); ++i )
	{
		*i = 0.;
	}

	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		unsigned iCntm1 = iCnt - 1;
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[tmp2]*dInitialLambda[iCntm1]);
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}

		dInitialLambda[iCntm1] = XCurr(iFirstIndex + iCnt);
	}


	return WorkVec;
}

/* JacobianRegularization - begin */

JacobianRegularization::JacobianRegularization(unsigned int uL,
	const Joint *j,
	const std::vector<doublereal>& c,
	flag fOut)
: Elem(uL, fOut),
JointRegularization(uL, j, fOut),
dC(c)
{
	NO_OP;
}

JacobianRegularization::~JacobianRegularization(void)
{
	NO_OP;
}

JointRegularization::Type
JacobianRegularization::GetJointRegularizationType(void) const
{
	return JointRegularization::JACOBIAN_REGULARIZATION;
}


void
JacobianRegularization::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetNumDof();
}

/* assemblaggio residuo */
SubVectorHandler&
JacobianRegularization::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
/*
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WorkVec.ResizeReset(iNumDofs);

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[iCnt - 1]*XCurr(iFirstIndex + iCnt)/dCoef);
	}
*/
	return WorkVec;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
JacobianRegularization::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetNumDof();

	WM.ResizeReset(iNumDofs, 0);

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[iCnt - 1]/dCoef);
	}

	return WorkMat;
}

void
JacobianRegularization::InitialWorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = pJ->iGetInitialNumDof();
}

VariableSubMatrixHandler& 
JacobianRegularization::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WM.ResizeReset(iNumDofs, 0);
	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WM.PutItem(iCnt, iFirstIndex + iCnt, iFirstIndex + iCnt,
			- dC[tmp2]);
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}
	}



	return WorkMat;
}

SubVectorHandler& 
JacobianRegularization::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{	
/*	
	integer iFirstIndex = pJ->iGetFirstIndex();
	unsigned iNumDofs = pJ->iGetInitialNumDof();

	WorkVec.ResizeReset(iNumDofs);

	unsigned tmp1 = 0;
	unsigned tmp2 = 0;

	for (unsigned iCnt = 1; iCnt <= iNumDofs; iCnt++) {
		WorkVec.PutItem(iCnt, iFirstIndex + iCnt,
			dC[tmp2]*XCurr(iFirstIndex + iCnt));
		if(tmp1 == 1)	{
			tmp1 = 0;
			tmp2++;
		} else	{
			tmp1++;
		}
	}
*/

	return WorkVec;
}

/* JacobianRegularization - end */

/* JointRegularization - end */


Elem *
ReadJointRegularization(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadRegularizationJoint");

	const char* sKeyWords[] = {
		"tikhonov",
		"dynamic",
		"jacobian",

		0
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		TIKHONOV = 0,
		DYNAMIC,
		JACOBIAN,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di vincolo */
	KeyWords CurrKeyWord = KeyWords(HP.IsKeyWord());

#ifdef DEBUG
	if (CurrKeyWord >= 0) {
		std::cout << "joint regularization type: " << sKeyWords[CurrKeyWord] << std::endl;
	}
#endif // DEBUG

	Elem* pEl = 0;

	switch (CurrKeyWord) {

	case TIKHONOV:
	case DYNAMIC:
	case JACOBIAN:
	{
		Elem *pTmpEl = pDM->pFindElem(Elem::JOINT, uLabel);
		Joint *pJ = dynamic_cast<Joint *>(pTmpEl);
		if (pJ == 0) {
			DrivenElem *pDE = dynamic_cast<DrivenElem *>(pTmpEl);
			if (pDE != 0) {
				pJ = dynamic_cast<Joint *>(pDE->pGetElem());
			}
		}

		if (pJ == 0) {
			silent_cerr(psElemNames[Elem::JOINT_REGULARIZATION]
				<< "(" << uLabel << "): "
				"unable to find "
				<< psElemNames[Elem::JOINT]
				<< "(" << uLabel << ") "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		unsigned iNumDofs = pJ->iGetNumDof();
		if (iNumDofs == 0) {
			silent_cerr(psElemNames[Elem::JOINT_REGULARIZATION]
				<< "(" << uLabel << "): "
				"no dofs for "
				<< psElemNames[Elem::JOINT]
				<< "(" << uLabel << ") "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

		}

		std::vector<doublereal> dC(iNumDofs);
		if (HP.IsKeyWord("list")) {
			for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
				dC[iCnt] = HP.GetReal();
			}

		} else if (HP.IsArg()) {
			doublereal d = HP.GetReal();

			for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
				dC[iCnt] = d;
			}

		} else {
			for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
				dC[iCnt] = 1e-6;
			}
		}

		switch (CurrKeyWord) {
		case TIKHONOV:
			SAFENEWWITHCONSTRUCTOR(pEl,
				TikhonovRegularization,
				TikhonovRegularization(uLabel, pJ, dC, 0));
			break;

		case DYNAMIC:
			SAFENEWWITHCONSTRUCTOR(pEl,
				DynamicRegularization,
				DynamicRegularization(uLabel, pJ, dC, 0));
			break;

		case JACOBIAN:
			SAFENEWWITHCONSTRUCTOR(pEl,
				JacobianRegularization,
				JacobianRegularization(uLabel, pJ, dC, 0));
			break;

		default:
			DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} break;

	/* Aggiungere qui altri vincoli */

	default:
		break;
	}

	/* Se non c'e' il punto e virgola finale */
	if (pEl == 0) {
		DEBUGCERR("");
		silent_cerr("error in allocation of "
			<< psElemNames[Elem::JOINT_REGULARIZATION]
			<< "(" << uLabel << ")" << std::endl);

		throw ErrMemory(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* ReadJointRegularization() */
