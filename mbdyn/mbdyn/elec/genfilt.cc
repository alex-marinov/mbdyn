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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "genfilt.h"

#ifdef USE_LAPACK
extern "C" int
__FC_DECL__(dgebal)(char *JOB, integer *N, doublereal *pdA, integer *LDA,
	integer *ILO, integer *IHI, doublereal *SCALE, integer *INFO);
extern "C" int
__FC_DECL__(dggbal)(char *JOB, integer *N, doublereal *pdA, integer *LDA,
	doublereal *pdB, integer *LDB, integer *ILO, integer *IHI,
	doublereal *LSCALE, doublereal *RSCALE, doublereal *WORK,
	integer *INFO);
#endif // USE_LAPACK

/* GenelStateSpaceSISO - begin */

GenelStateSpaceSISO::GenelStateSpaceSISO(unsigned int uLabel,
	const DofOwner* pDO,
	const ScalarDof& y,
	ScalarValue* u,
	unsigned int Order,
	doublereal* pE,
	doublereal* pA,
	doublereal* pB,
	doublereal* pC,
	doublereal D,
	bool bBalance,
	doublereal *pdX0,
	doublereal *pdXP0,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
SD_y(y), SV_u(u),
iNumDofs(Order),
pdE(pE), pdA(pA), pdB(pB), pdC(pC), dD(D),
pdX(0), pdXP(0)
{
	ASSERT(Order > 0);
	ASSERT(pdA != 0);
	ASSERT(pdB != 0);
	ASSERT(pdC != 0);
	ASSERT(SD_y.iOrder == 0);
	DEBUGCOUT("GenelStateSpaceSISO " << uLabel
		<< ", NumDofs: " << iNumDofs << std::endl);

#ifdef USE_LAPACK
	// try balancing the matrix?
	if (bBalance) {
		char JOB = 'S';
		integer N = Order;
		integer LDA = Order;
		integer LDB = Order;
		integer ILO;
		integer IHI;
		integer INFO;

		std::vector<doublereal> RSCALE(Order);
		std::vector<doublereal> LSCALE(Order);

		if (pdE != 0) {
			std::vector<doublereal> WORK(6*Order);

			(void)__FC_DECL__(dggbal)(&JOB, &N, pdA, &LDA,
				pdE, &LDB, &ILO, &IHI,
				&LSCALE[0], &RSCALE[0], &WORK[0], &INFO);

		} else {
			(void)__FC_DECL__(dgebal)(&JOB, &N, pdA, &LDA,
				&ILO, &IHI, &RSCALE[0], &INFO);

			for (unsigned i = 0; i < Order; i++) {
				LSCALE[i] = 1./RSCALE[i];
			}
		}

		if (INFO != 0) {
			silent_cout("GenelStateSpaceSISO(" << uLabel << "): "
				"balancing failed (ignored)" << std::endl);

		} else {
			for (unsigned i = 0; i < Order; i++) {
				pdB[i] *= RSCALE[i];
				pdC[i] *= LSCALE[i];
			}
		}
	}
#endif // USE_LAPACK

	SAFENEWARR(pdX, doublereal, 2*Order);
	pdXP = pdX + Order;

	for (unsigned i = 0; i < 2*Order; i++) {
		pdX[i] = 0.;
	}

	if (pdX0) {
		for (unsigned i = 0; i < Order; i++) {
			pdX[i] = pdX0[i];
		}

		if (pdXP0) {
			for (unsigned i = 0; i < Order; i++) {
				pdXP[i] = pdXP0[i];
			}
		}
	}
}

GenelStateSpaceSISO::~GenelStateSpaceSISO(void)
{
	if (pdX != 0) {
		SAFEDELETEARR(pdX);
	}

	if (pdC != 0) {
		SAFEDELETEARR(pdC);
	}

	if (pdB != 0) {
		SAFEDELETEARR(pdB);
	}

	if (pdA != 0) {
		SAFEDELETEARR(pdA);
	}

	if (pdE != 0) {
		SAFEDELETEARR(pdE);
	}

	if (SV_u != 0) {
		SAFEDELETE(SV_u);
	}
}

unsigned int
GenelStateSpaceSISO::iGetNumDof(void) const
{
	return iNumDofs;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
GenelStateSpaceSISO::GetDofType(unsigned int i) const
{
	ASSERT(i < iNumDofs);
	return DofOrder::DIFFERENTIAL;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelStateSpaceSISO::Restart(std::ostream& out) const
{
	return out << "GenelStateSpaceSISO: not implemented yet!" << std::endl;
}

/* Dimensioni del workspace */
void
GenelStateSpaceSISO::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = iNumDofs + 1;
	*piNumCols = iNumDofs + 1;

	// inputs may contribute to the Jacobian matrix
	if (dynamic_cast<ScalarDofValue *>(SV_u) != 0) {
		*piNumCols += 1;
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelStateSpaceSISO::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GenelStateSpaceSISO::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex() + 1;
	integer iColIndex_y = SD_y.pNode->iGetFirstColIndex() + 1;
	integer iFirstIndex = iGetFirstIndex();

	WM.PutRowIndex(iNumRows, iRowIndex_y);
	WM.PutColIndex(iNumCols, iColIndex_y);

	// inputs may contribute to the Jacobian matrix
	ScalarDofValue *SDV_u = dynamic_cast<ScalarDofValue *>(SV_u);
	if (SDV_u != 0) {
		ScalarDof& SD_u = dynamic_cast<ScalarDof &>(*SV_u);
		integer iColIndex_u = SD_u.pNode->iGetFirstRowIndex() + 1;
		integer iIdx_u = iNumCols - 1;

		WM.PutColIndex(iIdx_u, iColIndex_u);

		doublereal dd;
		if (SD_u.iOrder == 0) {
			dd = dCoef;
		} else {
			dd = 1.;
		}

		doublereal *pdb = pdB - 1;
		for (unsigned int i = iNumDofs; i > 0; i--) {
			WM.PutCoef(i, iIdx_u, pdb[i]*dd);
		}
	}

	WM.PutCoef(iNumRows, iNumCols, dCoef);

	doublereal* pda = pdA + iNumDofs*iNumDofs - 1;
	doublereal* pdc = pdC - 1;
	if (pdE) {
		doublereal* pde = pdE + iNumDofs*iNumDofs - 1;

		for (unsigned int i = iNumDofs; i > 0; i--) {
			WM.PutRowIndex(i, iFirstIndex + i);
			WM.PutColIndex(i, iFirstIndex + i);
			WM.PutCoef(iNumRows, i, -pdc[i]*dCoef);
			pde -= iNumDofs;
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				/* Attenzione: si assume A orientata per righe:
				 * a_11, a_12, ..., a_1n, a_21, ..., a_2n, ..., a_nn */
				WM.PutCoef(i, j, pde[j] - pda[j]*dCoef);
			}
		}

	} else {
		for (unsigned int i = iNumDofs; i > 0; i--) {
			WM.PutRowIndex(i, iFirstIndex + i);
			WM.PutColIndex(i, iFirstIndex + i);
			WM.PutCoef(iNumRows, i, -pdc[i]*dCoef);
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				/* Attenzione: si assume A orientata per righe:
				 * a_11, a_12, ..., a_1n, a_21, ..., a_2n, ..., a_nn */
				WM.PutCoef(i, j, -pda[j]*dCoef);
			}
			WM.IncCoef(i, i, 1.);
		}
	}

	return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler&
GenelStateSpaceSISO::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelStateSpaceSISO::AssRes()" << std::endl);

	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
	integer iFirstIndex = iGetFirstIndex();

	doublereal y = SD_y.pNode->dGetX();
	doublereal u = SV_u->dGetValue();

	WorkVec.PutRowIndex(iNumRows, iRowIndex_y);

	doublereal* pdx = pdX - 1;
	doublereal* pdxp = pdXP - 1;
	doublereal* pdc = pdC - 1;
	doublereal d = dD*u - y;
	for (unsigned int i = iNumDofs; i > 0; i--) {
		WorkVec.PutRowIndex(i, iFirstIndex + i);
		pdx[i] = XCurr(iFirstIndex + i);
		pdxp[i] = XPrimeCurr(iFirstIndex + i);
		d += pdc[i]*pdx[i];
	}

	WorkVec.PutCoef(iNumRows, d);

	doublereal* pda = pdA + iNumDofs*iNumDofs;
	doublereal* pdb = pdB;
	pdxp = pdXP;
	pdx = pdX;
	if (pdE) {
		doublereal *pde = pdE + iNumDofs*iNumDofs;
		for (unsigned int i = iNumDofs; i-- > 0; ) {
			d = pdb[i]*u;
			pde -= iNumDofs;
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				d += pda[j]*pdx[j] - pde[j]*pdxp[j];
			}
			WorkVec.PutCoef(i + 1, d);
		}

	} else {
		for (unsigned int i = iNumDofs; i-- > 0; ) {
			d = pdb[i]*u - pdxp[i];
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				d += pda[j]*pdx[j];
			}
			WorkVec.PutCoef(i + 1, d);
		}
	}

	return WorkVec;
}

void
GenelStateSpaceSISO::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = iGetFirstIndex() + 1;

	if (pdE == 0) {
		doublereal u = SV_u->dGetValue();
		doublereal* pda = pdA + iNumDofs*iNumDofs;
		doublereal* pdb = pdB;

		for (unsigned int i = iNumDofs; i-- > 0; ) {
			pdXP[i] = pdb[i]*u;
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				pdXP[i] += pda[j]*pdX[j];
			}
		}
	}

	for (unsigned i = 0; i < iNumDofs; i++) {
		X(iFirstIndex + i) = pdX[i];
		XP(iFirstIndex + i) = pdXP[i];
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
GenelStateSpaceSISO::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream &out(OH.Genels());
		out << std::setw(8) << GetLabel();
		for (unsigned int i = 0; i < iNumDofs; i++) {
			out << " " << pdX[i];
		}
		for (unsigned int i = 0; i < iNumDofs; i++) {
			out << " " << pdXP[i];
		}
		out << "\n";
	}
}

void
GenelStateSpaceSISO::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	unsigned iNodes = 1;
	if (dynamic_cast<NodeDof *>(SV_u)) {
		iNodes++;
	}
	connectedNodes.resize(iNodes);
	connectedNodes[0] = SD_y.pNode;
	if (dynamic_cast<NodeDof *>(SV_u)) {
		connectedNodes[1] = dynamic_cast<NodeDof *>(SV_u)->pNode;
	}
}

/* GenelStateSpaceSISO - end */


/* GenelStateSpaceMIMO - begin */

GenelStateSpaceMIMO::GenelStateSpaceMIMO(unsigned int uLabel,
	const DofOwner* pDO,
	unsigned int iNumOut,
	const ScalarDof* y,
	std::vector<ScalarValue *>& u,
	unsigned int Order,
	doublereal* pE,
	doublereal* pA,
	doublereal* pB,
	doublereal* pC,
	doublereal* pD,
	bool bBalance,
	doublereal *pdX0,
	doublereal *pdXP0,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
iNumOutputs(iNumOut), iNumInputs(u.size()),
pvSD_y(const_cast<ScalarDof *>(y)), SV_u(u),
iNumDofs(Order),
pdE(pE), pdA(pA), pdB(pB), pdC(pC), pdD(pD),
pdX(0), pdXP(0)
{
	ASSERT(iNumDofs > 0);
	ASSERT(iNumOutputs > 0);
	ASSERT(pvSD_y != 0);
	for (int i = iNumOutputs; i-- > 0; ) {
		ASSERT(pvSD_y[i].iOrder == 0);
	}
	ASSERT(iNumInputs > 0);
	ASSERT(pdA != 0);
	ASSERT(pdB != 0);
	ASSERT(pdC != 0);
	DEBUGCOUT("GenelStateSpaceMIMO " << uLabel
		<< ", NumDofs: " << iNumDofs << std::endl);

#ifdef USE_LAPACK
	// try balancing the matrix?
	if (bBalance) {
		int rc;
		char JOB = 'S';
		integer N = Order;
		integer LDA = Order;
		integer LDB = Order;
		integer ILO;
		integer IHI;
		integer INFO;

		std::vector<doublereal> RSCALE(Order);
		std::vector<doublereal> LSCALE(Order);

		if (pdE != 0) {
			std::vector<doublereal> WORK(6*Order);

			rc = __FC_DECL__(dggbal)(&JOB, &N, pdA, &LDA,
				pdE, &LDB, &ILO, &IHI,
				&LSCALE[0], &RSCALE[0], &WORK[0], &INFO);

		} else {
			rc = __FC_DECL__(dgebal)(&JOB, &N, pdA, &LDA,
				&ILO, &IHI, &RSCALE[0], &INFO);

			for (unsigned i = 0; i < Order; i++) {
				LSCALE[i] = 1./RSCALE[i];
			}
		}

		if (INFO != 0) {
			silent_cout("GenelStateSpaceMIMO(" << uLabel << "): "
				"balancing failed (ignored)" << std::endl);

		} else {
			doublereal *pd = pdB;

			for (unsigned i = 0; i < Order; i++) {
				doublereal d = RSCALE[i];
				for (unsigned j = 0; j < iNumInputs; j++) {
					*pd++ *= d;
				}
			}
	
			pd = pdC;
			for (unsigned j = 0; j < iNumOutputs; j++) {
				for (unsigned i = 0; i < Order; i++) {
					*pd++ *= LSCALE[i];
				}
			}
		}
	}
#endif // USE_LAPACK

	SAFENEWARR(pdX, doublereal, 2*Order);
	pdXP = pdX + Order;

	for (unsigned i = 0; i < 2*Order; i++) {
		pdX[i] = 0.;
	}

	if (pdX0) {
		for (unsigned i = 0; i < Order; i++) {
			pdX[i] = pdX0[i];
		}

		if (pdXP0) {
			for (unsigned i = 0; i < Order; i++) {
				pdXP[i] = pdXP0[i];
			}
		}
	}
}

GenelStateSpaceMIMO::~GenelStateSpaceMIMO(void)
{
	if (pdX != 0) {
		SAFEDELETEARR(pdX);
	}

	if (pdD != 0) {
		SAFEDELETEARR(pdD);
	}

	if (pdC != 0) {
		SAFEDELETEARR(pdC);
	}

	if (pdB != 0) {
		SAFEDELETEARR(pdB);
	}

	if (pdA != 0) {
		SAFEDELETEARR(pdA);
	}

	if (pdE != 0) {
		SAFEDELETEARR(pdE);
	}

	for (std::vector<ScalarValue *>::iterator i = SV_u.begin();
		i != SV_u.end(); ++i)
	{
		delete *i;
	}

	if (pvSD_y != 0) {
		SAFEDELETEARR(pvSD_y);
	}
}

unsigned int GenelStateSpaceMIMO::iGetNumDof(void) const
{
	return iNumDofs;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order GenelStateSpaceMIMO::GetDofType(unsigned int i) const
{
	ASSERT(i < iNumDofs);
	return DofOrder::DIFFERENTIAL;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream& GenelStateSpaceMIMO::Restart(std::ostream& out) const
{
	return out << "GenelStateSpaceMIMO: not implemented yet!" << std::endl;
}

/* Dimensioni del workspace */
void GenelStateSpaceMIMO::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = iNumDofs + iNumOutputs;
	*piNumCols = iNumDofs + iNumOutputs;

	// inputs may contribute to the Jacobian matrix
	for (unsigned int j = iNumInputs; j-- > 0; ) {
		if (dynamic_cast<ScalarDofValue *>(SV_u[j]) != 0) {
			*piNumCols += 1;
		}
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelStateSpaceMIMO::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GenelStateSpaceMIMO::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	// inputs may contribute to the Jacobian matrix
	integer iIdx_u = iNumDofs + iNumOutputs;
	for (unsigned int j = 0; j < iNumInputs; j++) {
		ScalarDofValue *SDV_u = dynamic_cast<ScalarDofValue *>(SV_u[j]);
		if (SDV_u != 0) {
			ScalarDof& SD_u = dynamic_cast<ScalarDof &>(*SV_u[j]);
			integer iColIndex_u = SD_u.pNode->iGetFirstRowIndex() + 1;
			iIdx_u += 1;
			WM.PutColIndex(iIdx_u, iColIndex_u);

			doublereal dd;
			if (SD_u.iOrder == 0) {
				dd = dCoef;
			} else {
				dd = 1.;
			}

			doublereal *pdb = &pdB[j];
			for (unsigned int i = 1; i <= iNumDofs; i++) {
				WM.PutCoef(i, iIdx_u, pdb[0]*dd);
				pdb += iNumInputs;
			}
		}
	}

	integer iFirstIndex = iGetFirstIndex();
	doublereal* pdc = pdC + iNumOutputs*iNumDofs - 1;
	for (unsigned int i = iNumOutputs; i > 0; i--) {
		integer iRowIndex_y = pvSD_y[i - 1].pNode->iGetFirstRowIndex() + 1;
		integer iColIndex_y = pvSD_y[i - 1].pNode->iGetFirstColIndex() + 1;

		WM.PutRowIndex(iNumDofs + i, iRowIndex_y);
		WM.PutColIndex(iNumDofs + i, iColIndex_y);
		// 1 sulla diagonale
		WM.PutCoef(iNumDofs + i, iNumDofs + i, dCoef);

		pdc -= iNumDofs;
		for (unsigned int j = iNumDofs; j > 0; j--) {
			/* Attenzione: si assume C orientata per righe:
			 * c_11, c_12, ..., c_1n, c_21, ..., c_2n, ..., c_nn */
			WM.PutCoef(iNumDofs + i, j, -pdc[j]*dCoef);
		}
	}

	doublereal* pda = pdA + iNumDofs*iNumDofs - 1;
	if (pdE) {
		doublereal* pde = pdE + iNumDofs*iNumDofs - 1;

		for (unsigned int i = iNumDofs; i > 0; i--) {
			WM.PutRowIndex(i, iFirstIndex + i);
			WM.PutColIndex(i, iFirstIndex + i);
			pde -= iNumDofs;
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				/* Attenzione: si assume A orientata per righe:
				 * a_11, a_12, ..., a_1n, a_21, ..., a_2n,
				 * ..., a_nn */
				WM.PutCoef(i, j, pde[j] - pda[j]*dCoef);
			}
		}

	} else {
		for (unsigned int i = iNumDofs; i > 0; i--) {
			WM.PutRowIndex(i, iFirstIndex + i);
			WM.PutColIndex(i, iFirstIndex + i);
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				/* Attenzione: si assume A orientata per righe:
				 * a_11, a_12, ..., a_1n, a_21, ..., a_2n,
				 * ..., a_nn */
				WM.PutCoef(i, j, -pda[j]*dCoef);
			}
			WM.IncCoef(i, i, 1.);
		}
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelStateSpaceMIMO::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelStateSpaceMIMO::AssRes()" << std::endl);

	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstIndex = iGetFirstIndex();

	doublereal* pdx = pdX-1;
	doublereal* pdxp = pdXP-1;
	for (unsigned int i = iNumDofs; i > 0; i--) {
		WorkVec.PutRowIndex(i, iFirstIndex+i);
		pdx[i] = XCurr(iFirstIndex+i);
		pdxp[i] = XPrimeCurr(iFirstIndex+i);
	}

	doublereal* pdc = pdC + iNumOutputs*iNumDofs - 1;
	if (pdD != 0) {
		doublereal* pdd = pdD + iNumOutputs*iNumInputs - 1;
		for (int i = iNumOutputs; i > 0; i--) {
			integer iRowIndex_y = pvSD_y[i - 1].pNode->iGetFirstRowIndex()+1;
			WorkVec.PutRowIndex(iNumDofs+i, iRowIndex_y);
			doublereal y = pvSD_y[i - 1].pNode->dGetX();
			doublereal d = -y;
			pdc -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				d += pdc[j]*pdx[j];
			}
			pdd -= iNumInputs;
			for (unsigned int j = iNumInputs; j > 0; j--) {
				d += pdd[j]*SV_u[j - 1]->dGetValue();
			}
			WorkVec.PutCoef(iNumDofs + i, d);
		}

	} else {
		for (int i = iNumOutputs; i > 0; i--) {
			integer iRowIndex_y = pvSD_y[i - 1].pNode->iGetFirstRowIndex() + 1;
			WorkVec.PutRowIndex(iNumDofs + i, iRowIndex_y);
			doublereal y = pvSD_y[i - 1].pNode->dGetX();
			doublereal d = -y;
			pdc -= iNumDofs;
			for (unsigned int j = iNumDofs; j > 0; j--) {
				d += pdc[j]*pdx[j];
			}
			WorkVec.PutCoef(iNumDofs + i, d);
		}
	}

	doublereal* pda = pdA + iNumDofs*iNumDofs;
	doublereal* pdb = pdB + iNumDofs*iNumInputs;
	pdxp = pdXP;
	pdx = pdX;
	if (pdE) {
		doublereal* pde = pdE + iNumDofs*iNumDofs;
		for (unsigned int i = iNumDofs; i-- > 0; ) {
			doublereal d = 0.;
			pdb -= iNumInputs;
			for (unsigned int j = iNumInputs; j-- > 0; ) {
				d += pdb[j]*SV_u[j]->dGetValue();
			}
			pde -= iNumDofs;
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				d += pda[j]*pdx[j] - pde[j]*pdxp[j];
			}
			WorkVec.PutCoef(i + 1, d);
		}

	} else {
		for (unsigned int i = iNumDofs; i-- > 0; ) {
			doublereal d = -pdxp[i];
			pdb -= iNumInputs;
			for (unsigned int j = iNumInputs; j-- > 0; ) {
				d += pdb[j]*SV_u[j]->dGetValue();
			}
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				d += pda[j]*pdx[j];
			}
			WorkVec.PutCoef(i + 1, d);
		}
	}

	return WorkVec;
}

void
GenelStateSpaceMIMO::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = iGetFirstIndex() + 1;

	if (pdE == 0) {
		doublereal* pda = pdA + iNumDofs*iNumDofs;
		doublereal* pdb = pdB + iNumDofs*iNumInputs;

		for (unsigned int i = iNumDofs; i-- > 0; ) {
			pdXP[i] = 0.;
			pdb -= iNumInputs;
			for (unsigned int j = iNumInputs; j-- > 0; ) {
				pdXP[i] += pdb[j]*SV_u[j]->dGetValue();
			}
			pda -= iNumDofs;
			for (unsigned int j = iNumDofs; j-- > 0; ) {
				pdXP[i] += pda[j]*pdX[j];
			}
		}
	}

	for (unsigned i = 0; i < iNumDofs; i++) {
		X(iFirstIndex + i) = pdX[i];
		XP(iFirstIndex + i) = pdXP[i];
	}
}

void
GenelStateSpaceMIMO::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out(OH.Genels());
		out << std::setw(8) << GetLabel();
		for (unsigned int i = 0; i < iNumDofs; i++) {
			out << " " << pdX[i];
		}
		for (unsigned int i = 0; i < iNumDofs; i++) {
			out << " " << pdXP[i];
		}
		out << "\n";
	}
}

void
GenelStateSpaceMIMO::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	unsigned i, iNodes = iNumOutputs;
	for (std::vector<ScalarValue *>::const_iterator u = SV_u.begin();
		u != SV_u.end(); ++u)
	{
		if (dynamic_cast<NodeDof *>(*u)) {
			iNodes++;
		}
	}

	connectedNodes.resize(iNodes);
	for (i = 0; i < iNumOutputs; i++) {
		connectedNodes[i] = pvSD_y[i].pNode;
	}

	for (std::vector<ScalarValue *>::const_iterator u = SV_u.begin();
		u != SV_u.end(); ++u)
	{
		NodeDof* ndp = dynamic_cast<NodeDof *>(*u);
		if (ndp) {
			connectedNodes[i++] = ndp->pNode;
		}
	}
}

/* GenelStateSpaceMIMO - end */

