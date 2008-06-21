/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "genfilt.h"

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
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
SD_y(y), SV_u(u),
iNumDofs(Order),
pdE(pE), pdA(pA), pdB(pB), pdC(pC), dD(D),
pdX(NULL), pdXP(NULL)
{
	ASSERT(Order > 0);
	ASSERT(pdA != NULL);
	ASSERT(pdB != NULL);
	ASSERT(pdC != NULL);
	ASSERT(SD_y.iOrder == 0);
	DEBUGCOUT("GenelStateSpaceSISO " << uLabel
		<< ", NumDofs: " << iNumDofs << std::endl);

	SAFENEWARR(pdX, doublereal, 2*Order);
	pdXP = pdX + Order;
}

GenelStateSpaceSISO::~GenelStateSpaceSISO(void)
{
	if (pdX != NULL) {
		SAFEDELETEARR(pdX);
	}

	if (pdC != NULL) {
		SAFEDELETEARR(pdC);
	}

	if (pdB != NULL) {
		SAFEDELETEARR(pdB);
	}

	if (pdA != NULL) {
		SAFEDELETEARR(pdA);
	}

	if (pdE != NULL) {
		SAFEDELETEARR(pdE);
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
		pdx[i] = XCurr.dGetCoef(iFirstIndex + i);
		pdxp[i] = XPrimeCurr.dGetCoef(iFirstIndex + i);
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
	std::vector<const Node *>& connectedNodes)
{
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
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
iNumOutputs(iNumOut), iNumInputs(u.size()),
pvSD_y(const_cast<ScalarDof *>(y)), SV_u(u),
iNumDofs(Order),
pdE(pE), pdA(pA), pdB(pB), pdC(pC), pdD(pD),
pdX(NULL), pdXP(NULL)
{
#ifdef DEBUG
	ASSERT(iNumDofs > 0);
	ASSERT(iNumOutputs > 0);
	ASSERT(pvSD_y != NULL);
	for (int i = iNumOutputs; i-- > 0; ) {
		ASSERT(pvSD_y[i].iOrder == 0);
	}
	ASSERT(iNumInputs > 0);
	ASSERT(pdA != NULL);
	ASSERT(pdB != NULL);
	ASSERT(pdC != NULL);
	DEBUGCOUT("GenelStateSpaceMIMO " << uLabel
		<< ", NumDofs: " << iNumDofs << std::endl);
#endif /* DEBUG */

	SAFENEWARR(pdX, doublereal, 2*Order);
	pdXP = pdX + Order;
}

GenelStateSpaceMIMO::~GenelStateSpaceMIMO(void)
{
	if (pdX != NULL) {
		SAFEDELETEARR(pdX);
	}

	if (pdD != NULL) {
		SAFEDELETEARR(pdD);
	}

	if (pdC != NULL) {
		SAFEDELETEARR(pdC);
	}

	if (pdB != NULL) {
		SAFEDELETEARR(pdB);
	}

	if (pdA != NULL) {
		SAFEDELETEARR(pdA);
	}

	if (pdE != NULL) {
		SAFEDELETEARR(pdE);
	}

	for (std::vector<ScalarValue *>::iterator i = SV_u.begin();
		i != SV_u.end(); i++)
	{
		delete *i;
	}

	if (pvSD_y != NULL) {
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
		pdx[i] = XCurr.dGetCoef(iFirstIndex+i);
		pdxp[i] = XPrimeCurr.dGetCoef(iFirstIndex+i);
	}

	doublereal* pdc = pdC + iNumOutputs*iNumDofs - 1;
	if (pdD != NULL) {
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
	std::vector<const Node *>& connectedNodes)
{
	unsigned i, iNodes = iNumOutputs;
	for (std::vector<ScalarValue *>::const_iterator u = SV_u.begin();
		u != SV_u.end(); u++)
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
		u != SV_u.end(); u++)
	{
		NodeDof* ndp = dynamic_cast<NodeDof *>(*u);
		if (ndp) {
			connectedNodes[i++] = ndp->pNode;
		}
	}
}

/* GenelStateSpaceMIMO - end */

