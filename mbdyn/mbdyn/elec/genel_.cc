/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "dataman.h"
#include "genel_.h"

/* GenelClamp - begin */

GenelClamp::GenelClamp(unsigned int uLabel,
	const DofOwner* pDO,
	const DriveCaller* pDC,
	const ScalarDof& sd,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
DriveOwner(pDC), SD(sd), dRct(0.)
{
	NO_OP;
}

GenelClamp::~GenelClamp(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

unsigned int
GenelClamp::iGetNumDof(void) const
{
	return 1;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
GenelClamp::GetDofType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::ALGEBRAIC;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
GenelClamp::GetEqType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::ALGEBRAIC;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelClamp::Restart(std::ostream& out) const
{
	return out;
}

/* Tipo di Genel */
Genel::Type
GenelClamp::GetGenelType(void) const
{
	return Genel::CLAMP;
}

/* Dimensioni del workspace */
void
GenelClamp::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 2;
}

void
GenelClamp::Output(OutputHandler& OH ) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Genels();
		out << std::setw(8) << GetLabel() << " " << dRct << std::endl;
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelClamp::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelClamp::AssJac()" << std::endl);

	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.Resize(2, 0);

	integer iRowIndex = SD.pNode->iGetFirstRowIndex() + 1;
	integer iColIndex = SD.pNode->iGetFirstColIndex() + 1;
	integer iFirstReactionIndex = iGetFirstIndex() + 1;

	WM.PutItem(1, iRowIndex, iFirstReactionIndex, 1.);
	WM.PutItem(2, iFirstReactionIndex, iColIndex, 1.);

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelClamp::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelClamp::AssRes()" << std::endl);

	WorkVec.Resize(2);

	integer iRowIndex = SD.pNode->iGetFirstRowIndex() + 1;
	integer iFirstReactionIndex = iGetFirstIndex() + 1;

	doublereal dVal = SD.pNode->dGetDofValue(1, SD.iOrder);
	dRct = XCurr(iFirstReactionIndex);

	WorkVec.PutItem(1, iRowIndex, -dRct);

	doublereal dConstr = dGet() - dVal;
	if (SD.iOrder == 0
		&& SD.pNode->GetDofType(0) == DofOrder::DIFFERENTIAL
		&& dCoef != 0.)
	{
		dConstr /= dCoef;
	}
	WorkVec.PutItem(2, iFirstReactionIndex, dConstr);

	return WorkVec;
}

void
GenelClamp::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (SD.iOrder == 0) {
		X.PutCoef(SD.pNode->iGetFirstRowIndex() + 1, dGet());
	} else if (SD.iOrder == 1) {
		XP.PutCoef(SD.pNode->iGetFirstRowIndex() + 1, dGet());
	}
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelClamp::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(1);
	connectedNodes[0] = SD.pNode;
}
/* ************************************************ */

/* GenelClamp - end */


/* GenelDistance - begin */

GenelDistance::GenelDistance(unsigned int uLabel,
	const DofOwner* pDO,
	const DriveCaller* pDC,
	const ScalarDof& sd1,
	const ScalarDof& sd2,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
DriveOwner(pDC), SD1(sd1), SD2(sd2), dRct(0.)
{
	NO_OP;
}

GenelDistance::~GenelDistance(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD1.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}

	pn2s = dynamic_cast<const Node2Scalar *>(SD2.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

unsigned int
GenelDistance::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
GenelDistance::GetDofType(unsigned int i ) const
{
	ASSERT(i == 0);
	return DofOrder::ALGEBRAIC;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelDistance::Restart(std::ostream& out) const
{
	return out;
}

/* Tipo di Genel */
Genel::Type
GenelDistance::GetGenelType(void) const
{
	return Genel::DISTANCE;
}

/* Dimensioni del workspace */
void
GenelDistance::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

void
GenelDistance::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Genels();
		out << std::setw(8) << GetLabel() << " " << dRct << std::endl;
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelDistance::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelDistance::AssJac()" << std::endl);

	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(4, 0);

	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex() + 1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex() + 1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex() + 1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex() + 1;
	integer iFirstReactionIndex = iGetFirstIndex() + 1;

	WM.PutItem(1, iNode1RowIndex, iFirstReactionIndex, -1.);
	WM.PutItem(2, iNode2RowIndex, iFirstReactionIndex, 1.);

	doublereal d = dCoef;
	if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {
		d = 1.;
	}

	if (SD1.iOrder == 1) {
		WM.PutItem(3, iFirstReactionIndex, iNode1ColIndex, -1.);
	} else {
		WM.PutItem(3, iFirstReactionIndex, iNode1ColIndex, -d);
	}

	if (SD2.iOrder == 1) {
		WM.PutItem(4, iFirstReactionIndex, iNode2ColIndex, 1.);
	} else {
		WM.PutItem(4, iFirstReactionIndex, iNode2ColIndex, d);
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelDistance::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelDistance::AssRes()" << std::endl);

	WorkVec.ResizeReset(3);

	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex() + 1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex() + 1;
	integer iFirstReactionIndex = iGetFirstIndex() + 1;

	doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
	doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);
	dRct = XCurr(iFirstReactionIndex);

	WorkVec.PutItem(1, iNode1RowIndex, dRct);
	WorkVec.PutItem(2, iNode2RowIndex, -dRct);

	if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {
		ASSERT(dCoef != 0.);
		WorkVec.PutItem(3, iFirstReactionIndex,
			(dGet() - dVal2 + dVal1)/dCoef);

	} else {
		WorkVec.PutItem(3, iFirstReactionIndex, dGet() - dVal2 + dVal1);
	}

	return WorkVec;
}

void
GenelDistance::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (SD2.iOrder == 0) {
		X.PutCoef(SD2.pNode->iGetFirstRowIndex() + 1,
			dGet() - SD2.pNode->dGetX()
			+ SD1.pNode->dGetDofValue(1, SD1.iOrder));
	} else if (SD2.iOrder == 1) {
		XP.PutCoef(SD2.pNode->iGetFirstRowIndex() + 1,
			dGet() - SD2.pNode->dGetXPrime()
			+ SD1.pNode->dGetDofValue(1, SD1.iOrder));
	}
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelDistance::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(2);
	connectedNodes[0] = SD1.pNode;
	connectedNodes[1] = SD2.pNode;
}
/* ************************************************ */

/* GenelDistance - end */


/* GenelSpring - begin */

GenelSpring::GenelSpring(unsigned int uLabel,
	const DofOwner* pDO,
	const ConstitutiveLaw1D* pCL,
	const ScalarDof& sd1,
	const ScalarDof& sd2,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
ConstitutiveLaw1DOwner(pCL), SD1(sd1), SD2(sd2), dVal(0.)
{
	NO_OP;
}

GenelSpring::~GenelSpring(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD1.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}

	pn2s = dynamic_cast<const Node2Scalar *>(SD2.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelSpring::Restart(std::ostream& out) const
{
	return out;
}

void
GenelSpring::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
}

/* Tipo di Genel */
Genel::Type
GenelSpring::GetGenelType(void) const
{
	return Genel::SPRING;
}

/* Dimensioni del workspace */
void
GenelSpring::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 2;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelSpring::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpring::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(2, 2);

	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex() + 1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex() + 1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex() + 1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNode1RowIndex);
	WM.PutColIndex(1, iNode1ColIndex);
	WM.PutRowIndex(2, iNode2RowIndex);
	WM.PutColIndex(2, iNode2ColIndex);

	doublereal dFDE = GetFDE();

	if (SD1.iOrder == 1) {
		WM.PutCoef(1, 1, dFDE);
		WM.PutCoef(2, 1, -dFDE);
	} else {
		WM.PutCoef(1, 1, dFDE*dCoef);
		WM.PutCoef(2, 1, -dFDE*dCoef);
	}

	if (SD2.iOrder == 1) {
		WM.PutCoef(1, 2, -dFDE);
		WM.PutCoef(2, 2, dFDE);
	} else {
		WM.PutCoef(1, 2, -dFDE*dCoef);
		WM.PutCoef(2, 2, dFDE*dCoef);
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelSpring::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpring::AssRes()" << std::endl);

	WorkVec.ResizeReset(2);

	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;

	doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
	doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);

	dVal = dVal2 - dVal1;
	ConstitutiveLaw1DOwner::Update(dVal, 0.);

	doublereal d = GetF();

	WorkVec.PutItem(1, iNode1RowIndex, d);
	WorkVec.PutItem(2, iNode2RowIndex, -d);

	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelSpring::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(2);
	connectedNodes[0] = SD1.pNode;
	connectedNodes[1] = SD2.pNode;
}
/* ************************************************ */

/* GenelSpring - end */


/* GenelSpringSupport - begin */

GenelSpringSupport::GenelSpringSupport(unsigned int uLabel,
	const DofOwner* pDO,
	const ConstitutiveLaw1D* pCL,
	const ScalarDof& sd, flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
ConstitutiveLaw1DOwner(pCL), SD(sd), dVal(0.)
{
	ASSERT(SD.iOrder == 0);
}

GenelSpringSupport::~GenelSpringSupport(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelSpringSupport::Restart(std::ostream& out) const
{
	return out;
}

void
GenelSpringSupport::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
}

/* Tipo di Genel */
Genel::Type
GenelSpringSupport::GetGenelType(void) const
{
	return Genel::SPRINGSUPPORT;
}

/* Dimensioni del workspace */
void
GenelSpringSupport::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelSpringSupport::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpringSupport::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);

	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex() + 1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);

	WM.PutCoef(1, 1, GetFDE()*dCoef);

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelSpringSupport::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpringSupport::AssRes()" << std::endl);

	WorkVec.ResizeReset(1);

	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex() + 1;

	dVal = SD.pNode->dGetX();
	ConstitutiveLaw1DOwner::Update(dVal, 0.);

	WorkVec.PutItem(1, iNodeRowIndex, -GetF());

	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelSpringSupport::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(1);
	connectedNodes[0] = SD.pNode;
}

unsigned int GenelSpringSupport::iGetNumPrivData(void) const
{
	return 1u;
}

unsigned int GenelSpringSupport::iGetPrivDataIdx(const char *s) const
{
	if (0 == strcmp(s, "F")) {
		return 1u;
	}

	silent_cerr("genel(" << GetLabel() << "): private data \"" << s << "\" not available" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal GenelSpringSupport::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1u:
		return GetF();
	default:
		silent_cerr("genel(" << GetLabel() << "): private data index " << i << " out of range" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}
/* ************************************************ */

/* GenelSpringSupport - end */


/* GenelCrossSpringSupport - begin */

GenelCrossSpringSupport::GenelCrossSpringSupport(unsigned int uLabel,
	const DofOwner* pDO,
	const ConstitutiveLaw1D* pCL,
	const ScalarDof& sdrow,
	const ScalarDof& sdcol,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
ConstitutiveLaw1DOwner(pCL),
SDRow(sdrow), SDCol(sdcol), dVal(0.)
{
	ASSERT(SDCol.iOrder == 0);
}

GenelCrossSpringSupport::~GenelCrossSpringSupport(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SDRow.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}

	pn2s = dynamic_cast<const Node2Scalar *>(SDCol.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelCrossSpringSupport::Restart(std::ostream& out) const
{
	return out;
}

void
GenelCrossSpringSupport::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
}

/* Tipo di Genel */
Genel::Type
GenelCrossSpringSupport::GetGenelType(void) const
{
	return Genel::CROSSSPRINGSUPPORT;
}

/* Dimensioni del workspace */
void
GenelCrossSpringSupport::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelCrossSpringSupport::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelCrossSpringSupport::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);

	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex() + 1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);

	WM.PutCoef(1, 1, GetFDE()*dCoef);

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelCrossSpringSupport::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelCrossSpringSupport::AssRes()" << std::endl);

	WorkVec.ResizeReset(1);

	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;

	dVal = SDCol.pNode->dGetX();
	ConstitutiveLaw1DOwner::Update(dVal, 0.);

	WorkVec.PutItem(1, iNodeRowIndex, -GetF());

	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelCrossSpringSupport::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(2);
	connectedNodes[0] = SDRow.pNode;
	connectedNodes[1] = SDCol.pNode;
}
/* ************************************************ */

/* GenelCrossSpringSupport - end */


/* GenelCrossSpringDamperSupport - begin */

GenelCrossSpringDamperSupport::GenelCrossSpringDamperSupport(
	unsigned int uLabel, const DofOwner* pDO,
	const ConstitutiveLaw1D* pCL,
	const ScalarDof& sdrow,
	const ScalarDof& sdcol,
	flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
ConstitutiveLaw1DOwner(pCL), SDRow(sdrow), SDCol(sdcol),
dVal(0.), dValPrime(0.)
{
	ASSERT(SDCol.iOrder == 0);
}

GenelCrossSpringDamperSupport::~GenelCrossSpringDamperSupport(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SDRow.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}

	pn2s = dynamic_cast<const Node2Scalar *>(SDCol.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelCrossSpringDamperSupport::Restart(std::ostream& out) const
{
	return out;
}

void
GenelCrossSpringDamperSupport::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dVal, dValPrime);
}

/* Tipo di Genel */
Genel::Type
GenelCrossSpringDamperSupport::GetGenelType(void) const
{
	return Genel::CROSSSPRINGDAMPERSUPPORT;
}

/* Dimensioni del workspace */
void
GenelCrossSpringDamperSupport::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelCrossSpringDamperSupport::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);

	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex() + 1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);

	WM.PutCoef(1, 1, GetFDE()*dCoef+GetFDEPrime());

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelCrossSpringDamperSupport::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssRes()" << std::endl);

	WorkVec.ResizeReset(1);

	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex() + 1;

	dVal = SDCol.pNode->dGetX();
	dValPrime = SDCol.pNode->dGetXPrime();
	ConstitutiveLaw1DOwner::Update(dVal, dValPrime);

	WorkVec.PutItem(1, iNodeRowIndex, -GetF());

	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelCrossSpringDamperSupport::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(2);
	connectedNodes[0] = SDRow.pNode;
	connectedNodes[1] = SDCol.pNode;
}
/* ************************************************ */

/* GenelCrossSpringDamperSupport - end */


/* GenelSpringDamperSupport - begin */

GenelSpringDamperSupport::GenelSpringDamperSupport(unsigned int uLabel,
	const DofOwner* pDO,
	const ConstitutiveLaw1D* pCL,
	const ScalarDof& sd, flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
ConstitutiveLaw1DOwner(pCL),
SD(sd), dVal(0.),dValPrime(0.)
{
	ASSERT(sd.pNode->GetDofType(0) == DofOrder::DIFFERENTIAL);
	ASSERT(sd.iOrder == 0);
}

GenelSpringDamperSupport::~GenelSpringDamperSupport(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelSpringDamperSupport::Restart(std::ostream& out) const
{
	return out;
}

void
GenelSpringDamperSupport::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dVal, dValPrime);
}

/* Tipo di Genel */
Genel::Type
GenelSpringDamperSupport::GetGenelType(void) const
{
	return Genel::SPRINGDAMPERSUPPORT;
}

/* Dimensioni del workspace */
void
GenelSpringDamperSupport::WorkSpaceDim(integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelSpringDamperSupport::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpringDamperSupport::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);

	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex() + 1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);

	WM.PutCoef(1, 1, GetFDE()*dCoef + GetFDEPrime());

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelSpringDamperSupport::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelSpringDamperSupport::AssRes()" << std::endl);

	WorkVec.ResizeReset(1);

	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;

	dVal = SD.pNode->dGetX();
	dValPrime = SD.pNode->dGetXPrime();

	ConstitutiveLaw1DOwner::Update(dVal, dValPrime);

	WorkVec.PutItem(1, iNodeRowIndex, -GetF());

	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelSpringDamperSupport::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(1);
	connectedNodes[0] = SD.pNode;
}
/* ************************************************ */

/* GenelSpringDamperSupport - end */


/* GenelMass - begin */

GenelMass::GenelMass(unsigned int uLabel,
	const DofOwner* pDO, const DriveCaller* pDC,
	const ScalarDof& sd, flag fOutput)
: Elem(uLabel, fOutput),
Genel(uLabel, pDO, fOutput),
DriveOwner(pDC), SD(sd)
{
	NO_OP;
}

GenelMass::~GenelMass(void)
{
	const Node2Scalar *pn2s;

	pn2s = dynamic_cast<const Node2Scalar *>(SD.pNode);
	if (pn2s) {
		SAFEDELETE(pn2s);
	}
}

unsigned int
GenelMass::iGetNumDof(void) const
{
	return 1;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
GenelMass::GetDofType(unsigned int i ) const
{
	ASSERT(i == 0);
	return DofOrder::DIFFERENTIAL;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenelMass::Restart(std::ostream& out) const
{
	return out;
}

/* Tipo di Genel */
Genel::Type
GenelMass::GetGenelType(void) const
{
	return Genel::MASS;
}

/* Dimensioni del workspace */
void
GenelMass::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 2;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GenelMass::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelMass::AssJac()" << std::endl);

	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(3, 0);

	integer iRowIndex = SD.pNode->iGetFirstRowIndex() + 1;
	integer iColIndex = SD.pNode->iGetFirstColIndex() + 1;
	integer iDerivativeIndex = iGetFirstIndex() + 1;

	WM.PutItem(1, iRowIndex, iDerivativeIndex, dGet());
	WM.PutItem(2, iDerivativeIndex, iColIndex, -1.);
	WM.PutItem(3, iDerivativeIndex, iDerivativeIndex, dCoef);

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
GenelMass::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering GenelMass::AssRes()" << std::endl);

	WorkVec.ResizeReset(2);

	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iDerivativeIndex = iGetFirstIndex()+1;

	doublereal dVal = SD.pNode->dGetXPrime();
	doublereal dDer = XCurr(iDerivativeIndex);
	doublereal dDerPrime = XPrimeCurr(iDerivativeIndex);

	WorkVec.PutItem(1, iRowIndex, -dGet()*dDerPrime);
	WorkVec.PutItem(2, iDerivativeIndex, dVal - dDer);

	return WorkVec;
}

void
GenelMass::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	X.PutCoef(iGetFirstIndex() + 1, SD.pNode->dGetXPrime());
}

/* *******PER IL SOLUTORE PARALLELO******** */
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
GenelMass::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(1);
	connectedNodes[0] = SD.pNode;
}
/* ************************************************ */

/* GenelMass - end */

