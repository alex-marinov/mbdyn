#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"

#include "module-chrono-interface.h"

ChronoInterfaceBaseElem::ChronoInterfaceBaseElem(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	module-chrono-interface						\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// Read element: obtain information from MBDyn script
	// or create a read function to read data
	// get the node
	pNode = dynamic_cast<const StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL)); // obtain the node
	// get the force
	f = HP.GetVec3();
}

ChronoInterfaceBaseElem::~ChronoInterfaceBaseElem(void)
{
	// destroy private data
	NO_OP;
}

void
ChronoInterfaceBaseElem::Output(OutputHandler& OH) const
{
	// should do something useful
	NO_OP;
}

void
ChronoInterfaceBaseElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 3;
}

VariableSubMatrixHandler& 
ChronoInterfaceBaseElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// similar like the structural force element;
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.ResizeReset(3, 3);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++)
	{
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}
	Vec3 TmpArm(pNode->GetRRef() * Zero3);
	WM.Sub(1, 1, Mat3x3(MatCrossCross, f, TmpArm * dCoef));
	return WorkMat;
}

SubVectorHandler& 
ChronoInterfaceBaseElem::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// refer to structural force element
	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	const Mat3x3& R(pNode->GetRCurr());
	Vec3 F(f);
	Vec3 M((R*Zero3).Cross(F));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);


	return WorkVec;
}

unsigned int
ChronoInterfaceBaseElem::iGetNumPrivData(void) const
{
	return 0;
}

int
ChronoInterfaceBaseElem::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ChronoInterfaceBaseElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ChronoInterfaceBaseElem::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ChronoInterfaceBaseElem::Restart(std::ostream& out) const
{
	return out << "# ModuleChronoInterface: not implemented" << std::endl;
}

unsigned int
ChronoInterfaceBaseElem::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ChronoInterfaceBaseElem::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}

VariableSubMatrixHandler&
ChronoInterfaceBaseElem::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// refer to the structural force element
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.ResizeReset(6, 6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(3+iCnt, iFirstVelocityIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(3+iCnt, iFirstVelocityIndex + iCnt);
	}


	Vec3 TmpArm(pNode->GetRRef()*Zero3);
	Vec3 TmpDir = f;
	const Vec3& Omega(pNode->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	Mat3x3 MTmp(MatCrossCross, TmpDir, TmpArm);

	WM.Sub(1, 1, MTmp);
	WM.Sub(4, 1, Mat3x3(MatCrossCross, TmpDir, Omega)*Mat3x3(MatCross, TmpArm));
	WM.Sub(4, 4, MTmp);

	return WorkMat;
}

SubVectorHandler& 
ChronoInterfaceBaseElem::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// refer to structural force element
	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);


	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex + iCnt);
	}


	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir(f);
	Vec3 TmpArm(R*Zero3);
	const Vec3& Omega(pNode->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm.Cross(TmpDir));
	
	WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir));

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ChronoInterfaceBaseElem>; // or new ChronoInterfaceElemRead;

	if (!SetUDE("ChronoInterface", rf)) {
		delete rf;

		silent_cerr("module-Chrono-Interface: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}