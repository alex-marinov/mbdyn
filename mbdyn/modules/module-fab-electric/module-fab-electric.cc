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

/**
    Library of electric components for "digital fabrication" machines (alpha version) [2013]
    Eduardo Okabe (okabe@unicamp.br)
    Postdoc CNPq at Aero/Polimi
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "module-fab-electric.h"

Resistor::Resistor(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	resistor			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("Resistor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("Resistor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the resistance [Ohms] from .mbd file:
   R1 = HP.GetReal();

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

Resistor::~Resistor(void)
{
	// destroy private data
	NO_OP;
}

void
Resistor::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on resistor
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << std::endl;
   }

}

unsigned int
Resistor::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Resistor::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
Resistor::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();

}

void
Resistor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

VariableSubMatrixHandler&
Resistor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);

   WM.IncCoef(1, 3, 1.);
   WM.DecCoef(2, 3, 1.);

   WM.DecCoef(3, 1, dCoef);
   WM.IncCoef(3, 2, dCoef);

   WM.IncCoef(3, 3, R1);

	return WorkMat;
}

SubVectorHandler&
Resistor::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

	doublereal i = XCurr(iFirstIndex);
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();

   DEBUGCOUT("Resistor::AssRes(), V1, V2, resistor current: " << V1 << ", " << V2 << ", " << i << std::endl);

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, V1 - V2 - R1*i);

	return WorkVec;
}

unsigned int
Resistor::iGetNumPrivData(void) const
{
	return 0;
}

int
Resistor::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
Resistor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pElec1;
	connectedNodes[1] = pElec2;
}

void
Resistor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
Resistor::Restart(std::ostream& out) const
{
	return out << "# Resistor: not implemented" << std::endl;
}

unsigned int
Resistor::iGetInitialNumDof(void) const
{
	return 0;
}

void
Resistor::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Resistor::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
Resistor::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

// Capacitor:

Capacitor::Capacitor(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Capacitor			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("Capacitor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("Capacitor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the Capacitance [F] from .mbd file:
   C1 = HP.GetReal();

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

Capacitor::~Capacitor(void)
{
	// destroy private data
	NO_OP;
}

void
Capacitor::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on Capacitor
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << std::endl;
   }

}

unsigned int
Capacitor::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Capacitor::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
Capacitor::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();

}

void
Capacitor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

VariableSubMatrixHandler&
Capacitor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);

   WM.IncCoef(1, 3, 1.);
   WM.DecCoef(2, 3, 1.);

   WM.IncCoef(3, 1, C1);
   WM.DecCoef(3, 2, C1);

   WM.IncCoef(3, 3, 1.);

	return WorkMat;
}

SubVectorHandler&
Capacitor::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

	doublereal i = XCurr(iFirstIndex);
	doublereal VP1 = pElec1->dGetXPrime();
	doublereal VP2 = pElec2->dGetXPrime();

   DEBUGCOUT("Capacitor::AssRes(), VP1, VP2, Capacitor current: " << VP1 << ", " << VP2 << ", " << i << std::endl);

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, i - C1*(VP1 - VP2));

	return WorkVec;
}

unsigned int
Capacitor::iGetNumPrivData(void) const
{
	return 0;
}

int
Capacitor::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
Capacitor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pElec1;
	connectedNodes[1] = pElec2;
}

void
Capacitor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
Capacitor::Restart(std::ostream& out) const
{
	return out << "# Capacitor: not implemented" << std::endl;
}

unsigned int
Capacitor::iGetInitialNumDof(void) const
{
	return 0;
}

void
Capacitor::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Capacitor::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
Capacitor::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


// Inductor:


Inductor::Inductor(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Inductor			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("Inductor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("Inductor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the Inductance [H] from .mbd file:
   L1 = HP.GetReal();

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

Inductor::~Inductor(void)
{
	// destroy private data
	NO_OP;
}

void
Inductor::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on Inductor
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << std::endl;
   }

}

unsigned int
Inductor::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Inductor::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
Inductor::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();

}

void
Inductor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

VariableSubMatrixHandler&
Inductor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);

   WM.IncCoef(1, 3, dCoef);
   WM.DecCoef(2, 3, dCoef);

   WM.DecCoef(3, 1, dCoef);
   WM.IncCoef(3, 2, dCoef);

   WM.IncCoef(3, 3, L1);

	return WorkMat;
}

SubVectorHandler&
Inductor::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

	doublereal i = XCurr(iFirstIndex);
	doublereal iP = XPrimeCurr(iFirstIndex);
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();

   DEBUGCOUT("Inductor::AssRes(), VP1, VP2, Inductor current: " << V1 << ", " << V2 << ", " << iP << std::endl);

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, V1 - V2 - L1*iP);

	return WorkVec;
}

unsigned int
Inductor::iGetNumPrivData(void) const
{
	return 0;
}

int
Inductor::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
Inductor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pElec1;
	connectedNodes[1] = pElec2;
}

void
Inductor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
Inductor::Restart(std::ostream& out) const
{
	return out << "# Inductor: not implemented" << std::endl;
}

unsigned int
Inductor::iGetInitialNumDof(void) const
{
	return 0;
}

void
Inductor::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Inductor::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
Inductor::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


// Diode model
// ref.: http://en.wikipedia.org/wiki/Diode

Diode::Diode(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Diode			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("Diode (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("Diode (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the forward saturation current [A] from .mbd file:
   IS = HP.GetReal();

	// Read the forward ideality factor [1-2] from .mbd file:
   NF = HP.GetReal();

	// Read the breakdown current [A] from .mbd file:
   IBV = HP.GetReal();

	// Read the breakdown voltage [V] from .mbd file:
   BV = HP.GetReal();

	// Read the reverse ideality factor from .mbd file:
   NR = HP.GetReal();

   if (HP.IsKeyWord("thermal voltage")) {
        DEBUGCOUT("Thermal voltage is supplied" << std::endl);
        VT = HP.GetReal();
   } else {
        VT = 25.85e-3;
   }

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

Diode::~Diode(void)
{
	// destroy private data
	NO_OP;
}

void
Diode::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on Diode
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << std::endl;
   }

}

unsigned int
Diode::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Diode::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
Diode::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();

}

void
Diode::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

VariableSubMatrixHandler&
Diode::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);

	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();
	doublereal dV = V1 - V2;

   doublereal Tmp1 = exp(dV/(NF*VT));
   doublereal Tmp2 = exp(-(dV+BV)/(NR*VT));
   doublereal Tmp3 = IS/(NF*VT)*Tmp1 + IBV/(NR*VT)*Tmp2;

   WM.IncCoef(1, 3, 1.);
   WM.DecCoef(2, 3, 1.);

   WM.IncCoef(3, 1, Tmp3*dCoef);
   WM.DecCoef(3, 2, Tmp3*dCoef);

   WM.DecCoef(3, 3, 1.);

	return WorkMat;
}

SubVectorHandler&
Diode::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

	doublereal i = XCurr(iFirstIndex);
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();
	doublereal dV = V1 - V2;

   DEBUGCOUT("Diode::AssRes(), VP1, VP2, Diode current: " << V1 << ", " << V2 << ", " << i << std::endl);

   doublereal Tmp1 = exp(dV/(NF*VT));
   doublereal Tmp2 = exp(-(dV+BV)/(NR*VT));

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, i - IS*Tmp1 + IBV*Tmp2);

	return WorkVec;
}

unsigned int
Diode::iGetNumPrivData(void) const
{
	return 0;
}

int
Diode::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
Diode::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pElec1;
	connectedNodes[1] = pElec2;
}

void
Diode::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
Diode::Restart(std::ostream& out) const
{
	return out << "# Diode: not implemented" << std::endl;
}

unsigned int
Diode::iGetInitialNumDof(void) const
{
	return 0;
}

void
Diode::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Diode::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
Diode::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



// Switch
// Model not tested!!!!!

Switch::Switch(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Switch			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("Switch (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("Switch (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Get driver caller associated to the switch state (open:S1=0/close:S1=1) from .mbd file:
   S1drv.Set(HP.GetDriveCaller());

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

Switch::~Switch(void)
{
	// destroy private data
	NO_OP;
}

void
Switch::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on switch
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << " " << dS1           // switch state
      << std::endl;
   }

}

unsigned int
Switch::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Switch::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
Switch::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   dS1 = S1drv.dGet();
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();
}

void
Switch::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

VariableSubMatrixHandler&
Switch::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);

	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();
   doublereal S1 = S1drv.dGet();

   WM.IncCoef(1, 3, 1.);
   WM.DecCoef(2, 3, 1.);

   WM.IncCoef(3, 1,-S1*dCoef);
   WM.DecCoef(3, 2,-S1*dCoef);

   WM.DecCoef(3, 3, 1.-S1);

	return WorkMat;
}

SubVectorHandler&
Switch::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

	doublereal i = XCurr(iFirstIndex);
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();
   doublereal S1 = S1drv.dGet();

   DEBUGCOUT("Switch::AssRes(), V1, V2, Switch current: " << V1 << ", " << V2 << ", " << i << std::endl);

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, S1*(V1-V2)+(1.-S1)*i);

	return WorkVec;
}

unsigned int
Switch::iGetNumPrivData(void) const
{
	return 0;
}

int
Switch::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
Switch::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pElec1;
	connectedNodes[1] = pElec2;
}

void
Switch::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
Switch::Restart(std::ostream& out) const
{
	return out << "# Switch: not implemented" << std::endl;
}

unsigned int
Switch::iGetInitialNumDof(void) const
{
	return 0;
}

void
Switch::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Switch::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
Switch::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


// Electrical Source (current or voltage)
// options: no control, current control, voltage control
// Model not tested!!!!!

ElectricalSource::ElectricalSource(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
         "									\n"
         "Module: 	ElectricalSource			\n"
         "									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   if (HP.IsKeyWord("current")) {
      source_type = CURRENTSOURCE;
   } else if (HP.IsKeyWord("voltage")) {
      source_type = VOLTAGESOURCE;
   } else {
      silent_cerr("unknown electrical source type at line "
                  << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   if (HP.IsKeyWord("control")) {
      if (HP.IsKeyWord("current")) {
         if (source_type==CURRENTSOURCE) {
            source_type = CURRENTCONTROLLED_CURRENTSOURCE;
         } else {
            source_type = CURRENTCONTROLLED_VOLTAGESOURCE;
         }
      } else if (HP.IsKeyWord("voltage")) {
         if (source_type==CURRENTSOURCE) {
            source_type = VOLTAGECONTROLLED_CURRENTSOURCE;
         } else {
            source_type = VOLTAGECONTROLLED_VOLTAGESOURCE;
         }
      } else {
      silent_cerr("unknown control type of electrical source at line "
                  << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
   }

   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
      pElecIn1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
      if (!pElecIn1) {
         silent_cerr("ElectricalSource (" << GetLabel() << "): electric node input (-) expected at line "
            << HP.GetLineData() << std::endl);
         throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }

      pElecIn2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
      if (!pElecIn2) {
         silent_cerr("ElectricalSource (" << GetLabel() << "): electric node input (+) expected at line "
            << HP.GetLineData() << std::endl);
         throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
   }

   pElecOut1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecOut1) {
      silent_cerr("ElectricalSource (" << GetLabel() << "): electric node output (-) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecOut2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecOut2) {
      silent_cerr("ElectricalSource (" << GetLabel() << "): electric node output (+) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Get driver caller associated to the Gain, Current or Voltage (depends on the type) from .mbd file:
   Vi1drv.Set(HP.GetDriveCaller());

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

ElectricalSource::~ElectricalSource(void)
{
	// destroy private data
	NO_OP;
}

void
ElectricalSource::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
      std::ostream& out = OH.Loadable();
      if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
         out << std::setw(8) << GetLabel()
            << " " << i_currIn         // input current
            << " " << i_currOut        // output current
            << " " << VoltageIn1       // voltage on node 1 (input)
            << " " << VoltageIn2       // voltage on node 2 (input)
            << " " << VoltageOut1      // voltage on node 1 (output)
            << " " << VoltageOut2      // voltage on node 2 (output)
            << " " << dVi1             // ElectricalSource gain
            << std::endl;
      } else {
         out << std::setw(8) << GetLabel()
            << " " << i_currOut        // current on ElectricalSource
            << " " << VoltageOut1      // voltage on node 1
            << " " << VoltageOut2      // voltage on node 2
            << " " << dVi1             // current or voltage
            << std::endl;
      }
   }

}

unsigned int
ElectricalSource::iGetNumDof(void) const
{
   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) return 2;
      else return 1;
}

DofOrder::Order
ElectricalSource::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
ElectricalSource::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric variables:
   dVi1 = Vi1drv.dGet();

   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
      i_currIn = X(iGetFirstIndex()+1);
      i_currOut = X(iGetFirstIndex()+2);
      VoltageIn1 = pElecIn1->dGetX();
      VoltageIn2 = pElecIn2->dGetX();
      VoltageOut1 = pElecOut1->dGetX();
      VoltageOut2 = pElecOut2->dGetX();
   } else {
      i_currOut = X(iGetFirstIndex()+1);
      VoltageOut1 = pElecOut1->dGetX();
      VoltageOut2 = pElecOut2->dGetX();
   }

}

void
ElectricalSource::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
	   *piNumRows = 6;
	   *piNumCols = 6;
   } else {
	   *piNumRows = 3;
	   *piNumCols = 3;
   }
}

VariableSubMatrixHandler&
ElectricalSource::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecOutNodeFirstIndex1 = pElecOut1->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex2 = pElecOut2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
	   integer iElecInNodeFirstIndex1 = pElecIn1->iGetFirstRowIndex() + 1;
	   integer iElecInNodeFirstIndex2 = pElecIn2->iGetFirstRowIndex() + 1;

	   WM.PutRowIndex(1, iElecInNodeFirstIndex1);
	   WM.PutRowIndex(2, iElecInNodeFirstIndex2);
	   WM.PutRowIndex(3, iElecOutNodeFirstIndex1);
	   WM.PutRowIndex(4, iElecOutNodeFirstIndex2);
	   WM.PutRowIndex(5, iFirstIndex);
	   WM.PutRowIndex(6, iFirstIndex + 1);

	   WM.PutColIndex(1, iElecInNodeFirstIndex1);
	   WM.PutColIndex(2, iElecInNodeFirstIndex2);
	   WM.PutColIndex(3, iElecOutNodeFirstIndex1);
	   WM.PutColIndex(4, iElecOutNodeFirstIndex2);
	   WM.PutColIndex(5, iFirstIndex);
	   WM.PutColIndex(6, iFirstIndex + 1);

   } else {

	   WM.PutRowIndex(1, iElecOutNodeFirstIndex1);
	   WM.PutRowIndex(2, iElecOutNodeFirstIndex2);
	   WM.PutRowIndex(3, iFirstIndex);

	   WM.PutColIndex(1, iElecOutNodeFirstIndex1);
	   WM.PutColIndex(2, iElecOutNodeFirstIndex2);
	   WM.PutColIndex(3, iFirstIndex);

   }

   doublereal G1 = Vi1drv.dGet();

   if (source_type==CURRENTSOURCE || source_type==VOLTAGESOURCE) {
      WM.IncCoef(1, 3, 1.);
      WM.DecCoef(2, 3, 1.);
   } else {
      WM.IncCoef(1, 5, 1.);
      WM.DecCoef(2, 5, 1.);
      WM.IncCoef(3, 6, 1.);
      WM.DecCoef(4, 6, 1.);
   }

   switch (source_type) {
      case VOLTAGESOURCE:
         WM.DecCoef(3, 1, dCoef);
         WM.IncCoef(3, 2, dCoef);
         break;
      case CURRENTSOURCE:
         WM.DecCoef(3, 3, 1.);
         break;
      case VOLTAGECONTROLLED_VOLTAGESOURCE:
         WM.IncCoef(5, 1, dCoef*G1);
         WM.DecCoef(5, 2, dCoef*G1);
         WM.DecCoef(5, 3, dCoef);
         WM.IncCoef(5, 4, dCoef);
         WM.DecCoef(6, 5, 1.);
         break;
      case CURRENTCONTROLLED_VOLTAGESOURCE:
         WM.DecCoef(5, 3, dCoef);
         WM.IncCoef(5, 4, dCoef);
         WM.DecCoef(5, 5, G1);
         WM.IncCoef(6, 1, dCoef);
         WM.DecCoef(6, 2, dCoef);
         break;
      case VOLTAGECONTROLLED_CURRENTSOURCE:
         WM.IncCoef(5, 1, dCoef*G1);
         WM.DecCoef(5, 2, dCoef*G1);
         WM.IncCoef(5, 6, 1.);
         WM.DecCoef(6, 5, 1.);
         break;
      case CURRENTCONTROLLED_CURRENTSOURCE:
         WM.DecCoef(5, 5, G1);
         WM.IncCoef(5, 6, 1.);
         WM.IncCoef(6, 1, dCoef);
         WM.DecCoef(6, 2, dCoef);
         break;
   }

	return WorkMat;
}

SubVectorHandler&
ElectricalSource::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iElecOutNodeFirstIndex1 = pElecOut1->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex2 = pElecOut2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   doublereal V1in = 0.;
   doublereal V2in = 0.;
   doublereal V1out = pElecOut1->dGetX();
	doublereal V2out = pElecOut2->dGetX();
   doublereal G1 = Vi1drv.dGet();

   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {

	   integer iElecInNodeFirstIndex1 = pElecIn1->iGetFirstRowIndex() + 1;
	   integer iElecInNodeFirstIndex2 = pElecIn2->iGetFirstRowIndex() + 1;

	   WorkVec.PutRowIndex(1, iElecInNodeFirstIndex1);
	   WorkVec.PutRowIndex(2, iElecInNodeFirstIndex2);
	   WorkVec.PutRowIndex(3, iElecOutNodeFirstIndex1);
	   WorkVec.PutRowIndex(4, iElecOutNodeFirstIndex2);
	   WorkVec.PutRowIndex(5, iFirstIndex);
	   WorkVec.PutRowIndex(6, iFirstIndex + 1);

      doublereal i_currIn = XCurr(iFirstIndex);
      doublereal i_currOut = XCurr(iFirstIndex + 1);
      doublereal V1in = pElecIn1->dGetX();
	   doublereal V2in = pElecIn2->dGetX();

      DEBUGCOUT("ElectricalSource::AssRes(), V1in, V2in, i_currIn, V1out, V2out, i_currOut, G1: " << V1in
         << ", " << V2in << ", " << i_currIn << ", " << V1out << ", " << V2out << ", "
         << i_currOut  << ", " << G1 << std::endl);

   } else {

	   WorkVec.PutRowIndex(1, iElecOutNodeFirstIndex1);
	   WorkVec.PutRowIndex(2, iElecOutNodeFirstIndex2);
	   WorkVec.PutRowIndex(3, iFirstIndex);

      doublereal i_currOut = XCurr(iFirstIndex);

      DEBUGCOUT("ElectricalSource::AssRes(), V1out, V2out, i_currOut, G1: " <<
         V1out << ", " << V2out << ", " << i_currOut  << ", " << G1 << std::endl);

   }


   if (source_type==CURRENTSOURCE || source_type==VOLTAGESOURCE) {
      WorkVec.DecCoef(1, i_currOut);
      WorkVec.IncCoef(2, i_currOut);
   } else {
      WorkVec.DecCoef(1, i_currIn);
      WorkVec.IncCoef(2, i_currIn);
      WorkVec.DecCoef(3, i_currOut);
      WorkVec.IncCoef(4, i_currOut);
   }

   switch (source_type) {
      case VOLTAGESOURCE:
         WorkVec.IncCoef(3, V2out - V1out - G1);
         break;
      case CURRENTSOURCE:
         WorkVec.IncCoef(3, i_currOut - G1);
         break;
      case VOLTAGECONTROLLED_VOLTAGESOURCE:
         WorkVec.IncCoef(5, G1*(V2in - V1in) - (V2out - V1out));
         WorkVec.IncCoef(6, i_currIn);
         break;
      case CURRENTCONTROLLED_VOLTAGESOURCE:
         WorkVec.IncCoef(5, G1*i_currIn - (V2out - V1out));
         WorkVec.IncCoef(6, V2in - V1in);
         break;
      case VOLTAGECONTROLLED_CURRENTSOURCE:
         WorkVec.IncCoef(5, G1*(V2in - V1in) - i_currOut);
         WorkVec.IncCoef(6, i_currIn);
         break;
      case CURRENTCONTROLLED_CURRENTSOURCE:
         WorkVec.IncCoef(5, G1*i_currIn - i_currOut);
         WorkVec.IncCoef(6, V2in - V1in);
         break;
   }

	return WorkVec;
}

unsigned int
ElectricalSource::iGetNumPrivData(void) const
{
	return 0;
}

int
ElectricalSource::iGetNumConnectedNodes(void) const
{
   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) return 4;
      else return 2;
}

void
ElectricalSource::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
   if (source_type!=CURRENTSOURCE && source_type!=VOLTAGESOURCE) {
	   connectedNodes.resize(4);
	   connectedNodes[0] = pElecIn1;
	   connectedNodes[1] = pElecIn2;
	   connectedNodes[2] = pElecOut1;
	   connectedNodes[3] = pElecOut2;
   } else {
	   connectedNodes.resize(2);
	   connectedNodes[0] = pElecOut1;
	   connectedNodes[1] = pElecOut2;
   }
}

void
ElectricalSource::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ElectricalSource::Restart(std::ostream& out) const
{
	return out << "# ElectricalSource: not implemented" << std::endl;
}

unsigned int
ElectricalSource::iGetInitialNumDof(void) const
{
	return 0;
}

void
ElectricalSource::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ElectricalSource::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
ElectricalSource::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



// Ideal Tranformer
// Model not tested!!!!!

IdealTransformer::IdealTransformer(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
         "									\n"
         "Module: 	Ideal Transformer			\n"
         "									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElecIn1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecIn1) {
      silent_cerr("IdealTransformer (" << GetLabel() << "): electric node input (-) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecIn2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecIn2) {
      silent_cerr("IdealTransformer (" << GetLabel() << "): electric node input (+) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecOut1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecOut1) {
      silent_cerr("IdealTransformer (" << GetLabel() << "): electric node output (-) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecOut2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecOut2) {
      silent_cerr("IdealTransformer (" << GetLabel() << "): electric node output (+) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Get driver caller associated to the transformer ratio (Nout/Nin) from .mbd file:
   G1drv.Set(HP.GetDriveCaller());

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

IdealTransformer::~IdealTransformer(void)
{
	// destroy private data
	NO_OP;
}

void
IdealTransformer::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
      std::ostream& out = OH.Loadable();
      out << std::setw(8) << GetLabel()
         << " " << i_currIn         // input current
         << " " << i_currOut        // output current
         << " " << VoltageIn1       // voltage on node 1 (input)
         << " " << VoltageIn2       // voltage on node 2 (input)
         << " " << VoltageOut1      // voltage on node 1 (output)
         << " " << VoltageOut2      // voltage on node 2 (output)
         << " " << dG1             // IdealTransformer gain
         << std::endl;
   }

}

unsigned int
IdealTransformer::iGetNumDof(void) const
{
   return 2;
}

DofOrder::Order
IdealTransformer::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
IdealTransformer::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletrical variables:
   dG1 = G1drv.dGet();
   i_currIn = X(iGetFirstIndex()+1);
   i_currOut = X(iGetFirstIndex()+2);
   VoltageIn1 = pElecIn1->dGetX();
   VoltageIn2 = pElecIn2->dGetX();
   VoltageOut1 = pElecOut1->dGetX();
   VoltageOut2 = pElecOut2->dGetX();

}

void
IdealTransformer::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = 6;
   *piNumCols = 6;
}

VariableSubMatrixHandler&
IdealTransformer::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

   integer iElecInNodeFirstIndex1 = pElecIn1->iGetFirstRowIndex() + 1;
   integer iElecInNodeFirstIndex2 = pElecIn2->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex1 = pElecOut1->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex2 = pElecOut2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   WM.PutRowIndex(1, iElecInNodeFirstIndex1);
   WM.PutRowIndex(2, iElecInNodeFirstIndex2);
   WM.PutRowIndex(3, iElecOutNodeFirstIndex1);
   WM.PutRowIndex(4, iElecOutNodeFirstIndex2);
   WM.PutRowIndex(5, iFirstIndex);
   WM.PutRowIndex(6, iFirstIndex + 1);

   WM.PutColIndex(1, iElecInNodeFirstIndex1);
   WM.PutColIndex(2, iElecInNodeFirstIndex2);
   WM.PutColIndex(3, iElecOutNodeFirstIndex1);
   WM.PutColIndex(4, iElecOutNodeFirstIndex2);
   WM.PutColIndex(5, iFirstIndex);
   WM.PutColIndex(6, iFirstIndex + 1);

   doublereal G1 = G1drv.dGet();

   WM.IncCoef(1, 5, 1.);
   WM.DecCoef(2, 5, 1.);
   WM.IncCoef(3, 6, 1.);
   WM.DecCoef(4, 6, 1.);
   WM.DecCoef(5, 1, dCoef*G1);
   WM.IncCoef(5, 2, dCoef*G1);
   WM.IncCoef(5, 3, dCoef);
   WM.DecCoef(5, 4, dCoef);
   WM.DecCoef(6, 5, 1.);
   WM.IncCoef(6, 6, G1);

	return WorkMat;
}

SubVectorHandler&
IdealTransformer::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   integer iElecInNodeFirstIndex1 = pElecIn1->iGetFirstRowIndex() + 1;
   integer iElecInNodeFirstIndex2 = pElecIn2->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex1 = pElecOut1->iGetFirstRowIndex() + 1;
	integer iElecOutNodeFirstIndex2 = pElecOut2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   doublereal i_currIn = XCurr(iFirstIndex);
   doublereal i_currOut = XCurr(iFirstIndex + 1);
   doublereal V1in = pElecIn1->dGetX();
   doublereal V2in = pElecIn2->dGetX();
   doublereal V1out = pElecOut1->dGetX();
	doublereal V2out = pElecOut2->dGetX();
   doublereal G1 = G1drv.dGet();

   WorkVec.PutRowIndex(1, iElecInNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecInNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iElecOutNodeFirstIndex1);
   WorkVec.PutRowIndex(4, iElecOutNodeFirstIndex2);
   WorkVec.PutRowIndex(5, iFirstIndex);
   WorkVec.PutRowIndex(6, iFirstIndex + 1);


   DEBUGCOUT("IdealTransformer::AssRes(), V1in, V2in, i_currIn, V1out, V2out, i_currOut, G1: " << V1in
      << ", " << V2in << ", " << i_currIn << ", " << V1out << ", " << V2out << ", "
      << i_currOut  << ", " << G1 << std::endl);

   WorkVec.DecCoef(1, i_currIn);
   WorkVec.IncCoef(2, i_currIn);
   WorkVec.DecCoef(3, i_currOut);
   WorkVec.IncCoef(4, i_currOut);
   WorkVec.IncCoef(5, (V2out - V1out) - G1*(V2in - V1in));
   WorkVec.IncCoef(6, i_currIn - G1*i_currOut);

	return WorkVec;
}

unsigned int
IdealTransformer::iGetNumPrivData(void) const
{
	return 0;
}

int
IdealTransformer::iGetNumConnectedNodes(void) const
{
   return 4;
}

void
IdealTransformer::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
   connectedNodes.resize(4);
   connectedNodes[0] = pElecIn1;
   connectedNodes[1] = pElecIn2;
   connectedNodes[2] = pElecOut1;
   connectedNodes[3] = pElecOut2;
}

void
IdealTransformer::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
IdealTransformer::Restart(std::ostream& out) const
{
	return out << "# IdealTransformer: not implemented" << std::endl;
}

unsigned int
IdealTransformer::iGetInitialNumDof(void) const
{
	return 0;
}

void
IdealTransformer::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
IdealTransformer::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
IdealTransformer::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



// Operational Amplifier
// Model not tested!!!!!

OperationalAmplifier::OperationalAmplifier(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
         "									\n"
         "Module: 	Operational Amplifier			\n"
         "									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   pElecNeg = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecNeg) {
      silent_cerr("OperationalAmplifier (" << GetLabel() << "): electric node input (-) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecPos = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecPos) {
      silent_cerr("OperationalAmplifier (" << GetLabel() << "): electric node input (+) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecRef = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecRef) {
      silent_cerr("OperationalAmplifier (" << GetLabel() << "): electric node reference expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecOut = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecOut) {
      silent_cerr("OperationalAmplifier (" << GetLabel() << "): electric node output (+) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   // Amplifier gain (ideal amplifier has infinite gain):
   Gain = 1.e8;

   if (HP.IsKeyWord("gain")) Gain = HP.GetReal();

   // Amplifier input resistance (ideal amplifier has infinite input impedance):
   Rinput = 1.e8;

   if (HP.IsKeyWord("input resistance")) Rinput = HP.GetReal();

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

OperationalAmplifier::~OperationalAmplifier(void)
{
	// destroy private data
	NO_OP;
}

void
OperationalAmplifier::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
      std::ostream& out = OH.Loadable();
      out << std::setw(8) << GetLabel()
         << " " << i_curr           // output current
         << " " << VoltageNeg       // voltage on node 1 (input)
         << " " << VoltagePos       // voltage on node 2 (input)
         << " " << VoltageOut       // voltage on node 1 (output)
         << " " << VoltageRef       // voltage on node 2 (output)
         << std::endl;
   }

}

unsigned int
OperationalAmplifier::iGetNumDof(void) const
{
   return 2;
}

DofOrder::Order
OperationalAmplifier::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
OperationalAmplifier::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletrical variables:
   i_curr = X(iGetFirstIndex()+2);
   VoltageNeg = pElecNeg->dGetX();
   VoltagePos = pElecPos->dGetX();
   VoltageOut = pElecOut->dGetX();
   VoltageRef = pElecRef->dGetX();

}

void
OperationalAmplifier::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = 6;
   *piNumCols = 6;
}

VariableSubMatrixHandler&
OperationalAmplifier::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

   integer iElecNodeFirstIndexNeg = pElecNeg->iGetFirstRowIndex() + 1;
   integer iElecNodeFirstIndexPos = pElecPos->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexOut = pElecOut->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexRef = pElecRef->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   WM.PutRowIndex(1, iElecNodeFirstIndexNeg);
   WM.PutRowIndex(2, iElecNodeFirstIndexPos);
   WM.PutRowIndex(3, iElecNodeFirstIndexRef);
   WM.PutRowIndex(4, iElecNodeFirstIndexOut);
   WM.PutRowIndex(5, iFirstIndex);
   WM.PutRowIndex(6, iFirstIndex + 1);

   WM.PutColIndex(1, iElecNodeFirstIndexNeg);
   WM.PutColIndex(2, iElecNodeFirstIndexPos);
   WM.PutColIndex(3, iElecNodeFirstIndexRef);
   WM.PutColIndex(4, iElecNodeFirstIndexOut);
   WM.PutColIndex(5, iFirstIndex);
   WM.PutColIndex(6, iFirstIndex + 1);

   WM.IncCoef(1, 5, 1.);
   WM.DecCoef(2, 5, 1.);
   WM.IncCoef(3, 6, 1.);
   WM.DecCoef(4, 6, 1.);
   WM.DecCoef(5, 1, dCoef*Gain);
   WM.IncCoef(5, 2, dCoef*Gain);
   WM.IncCoef(5, 3, dCoef);
   WM.DecCoef(5, 4, dCoef);
   WM.IncCoef(6, 1, dCoef);
   WM.DecCoef(6, 2, dCoef);
   WM.IncCoef(6, 5, Rinput);

	return WorkMat;
}

SubVectorHandler&
OperationalAmplifier::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   integer iElecNodeFirstIndexNeg = pElecNeg->iGetFirstRowIndex() + 1;
   integer iElecNodeFirstIndexPos = pElecPos->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexOut = pElecOut->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexRef = pElecRef->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   doublereal i_currIn = XCurr(iFirstIndex);
   doublereal i_currOut = XCurr(iFirstIndex + 1);
   doublereal Vneg = pElecNeg->dGetX();
   doublereal Vpos = pElecPos->dGetX();
   doublereal Vref = pElecRef->dGetX();
	doublereal Vout = pElecOut->dGetX();

   WorkVec.PutRowIndex(1, iElecNodeFirstIndexNeg);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndexPos);
   WorkVec.PutRowIndex(3, iElecNodeFirstIndexRef);
   WorkVec.PutRowIndex(4, iElecNodeFirstIndexOut);
   WorkVec.PutRowIndex(5, iFirstIndex);
   WorkVec.PutRowIndex(6, iFirstIndex + 1);

   DEBUGCOUT("OperationalAmplifier::AssRes(), Vneg, Vpos, Vref, Vout, i_currIn: " << Vneg
      << ", " << Vpos << ", " << Vref << ", " << Vout << ", " << i_currIn << std::endl);

   WorkVec.DecCoef(1, i_currIn);
   WorkVec.IncCoef(2, i_currIn);
   WorkVec.DecCoef(3, i_currOut);
   WorkVec.IncCoef(4, i_currOut);
   WorkVec.IncCoef(5, Vout - Vref - Gain*(Vpos - Vneg));
   WorkVec.IncCoef(6, Vpos - Vneg - Rinput*i_currIn);

	return WorkVec;
}

unsigned int
OperationalAmplifier::iGetNumPrivData(void) const
{
	return 0;
}

int
OperationalAmplifier::iGetNumConnectedNodes(void) const
{
   return 4;
}

void
OperationalAmplifier::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
   connectedNodes.resize(4);
   connectedNodes[0] = pElecNeg;
   connectedNodes[1] = pElecPos;
   connectedNodes[2] = pElecRef;
   connectedNodes[3] = pElecOut;
}

void
OperationalAmplifier::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
OperationalAmplifier::Restart(std::ostream& out) const
{
	return out << "# OperationalAmplifier: not implemented" << std::endl;
}

unsigned int
OperationalAmplifier::iGetInitialNumDof(void) const
{
	return 0;
}

void
OperationalAmplifier::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
OperationalAmplifier::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
OperationalAmplifier::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



// Bipolar junction transistor
// Model not tested!!!!!

BipolarTransistor::BipolarTransistor(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
         "									\n"
         "Module: 	Operational Amplifier			\n"
         "									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   if (HP.IsKeyWord("npn")) {
      bjt_type = NPN;
   } else if (HP.IsKeyWord("pnp")) {
      bjt_type = PNP;
   } else {
      silent_cerr("unknown bipolar transistor type at line "
                  << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecC = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecC) {
      silent_cerr("BipolarTransistor (" << GetLabel() << "): electric node (collector) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecB = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecB) {
      silent_cerr("BipolarTransistor (" << GetLabel() << "): electric node (base) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   pElecE = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElecE) {
      silent_cerr("BipolarTransistor (" << GetLabel() << "): electric node (emitter) expected at line "
         << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the B-E leakage saturation current [A] from .mbd file:
   Ise = HP.GetReal();

	// Read the B-C leakage saturation current [A] from .mbd file:
   Isc = HP.GetReal();

	// Read ideal maximum forward beta from .mbd file:
   Bf = HP.GetReal();

	// Read ideal maximum reverse beta from .mbd file:
   Br = HP.GetReal();

   if (HP.IsKeyWord("thermal voltage")) {
        DEBUGCOUT("Thermal voltage is supplied" << std::endl);
        Vt = HP.GetReal();
   } else {
        Vt = 25.85e-3;
   }

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

BipolarTransistor::~BipolarTransistor(void)
{
	// destroy private data
	NO_OP;
}

void
BipolarTransistor::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
      std::ostream& out = OH.Loadable();
      out << std::setw(8) << GetLabel()
         << " " << icurrC           // collector current
         << " " << icurrB           // base current
         << " " << icurrE           // emitter current
         << " " << VoltageC         // voltage on collector node
         << " " << VoltageB         // voltage on base node
         << " " << VoltageE         // voltage on emitter node
         << std::endl;
   }

}

unsigned int
BipolarTransistor::iGetNumDof(void) const
{
   return 3;
}

DofOrder::Order
BipolarTransistor::GetDofType(unsigned int i) const
{
	return DofOrder::ALGEBRAIC;
}

void
BipolarTransistor::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletrical variables:
   icurrC = X(iGetFirstIndex()+1);
   icurrB = X(iGetFirstIndex()+2);
   icurrE = X(iGetFirstIndex()+3);
   VoltageC = pElecC->dGetX();
   VoltageB = pElecB->dGetX();
   VoltageE = pElecE->dGetX();
}

void
BipolarTransistor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = 6;
   *piNumCols = 6;
}

VariableSubMatrixHandler&
BipolarTransistor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

   integer iElecNodeFirstIndexC = pElecC->iGetFirstRowIndex() + 1;
   integer iElecNodeFirstIndexB = pElecB->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexE = pElecE->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   WM.PutRowIndex(1, iElecNodeFirstIndexC);
   WM.PutRowIndex(2, iElecNodeFirstIndexB);
   WM.PutRowIndex(3, iElecNodeFirstIndexE);
   WM.PutRowIndex(4, iFirstIndex);
   WM.PutRowIndex(5, iFirstIndex + 1);
   WM.PutRowIndex(6, iFirstIndex + 2);

   WM.PutColIndex(1, iElecNodeFirstIndexC);
   WM.PutColIndex(2, iElecNodeFirstIndexB);
   WM.PutColIndex(3, iElecNodeFirstIndexE);
   WM.PutColIndex(4, iFirstIndex);
   WM.PutColIndex(5, iFirstIndex + 1);
   WM.PutColIndex(6, iFirstIndex + 2);

   doublereal Af = Bf/(Bf+1.);
   doublereal Ar = Br/(Br+1.);
   doublereal Vbe = pElecB->dGetX() - pElecE->dGetX();
   doublereal Vbc = pElecB->dGetX() - pElecC->dGetX();

   WM.DecCoef(1, 4, 1.);
   WM.DecCoef(2, 5, 1.);
   WM.IncCoef(3, 6, 1.);

   if (bjt_type==NPN) {
      doublereal Tmp1 = (Ise/Vt)*exp(Vbe/Vt);
      doublereal Tmp2 = (Isc/Vt)*exp(Vbc/Vt);

      WM.IncCoef(4, 1, Ar*Tmp2*dCoef);
      WM.IncCoef(4, 2, (Tmp1 - Ar*Tmp2)*dCoef);
      WM.DecCoef(4, 3, Tmp1*dCoef);

      WM.IncCoef(5, 1, Tmp2*dCoef);
      WM.IncCoef(5, 2, (Af*Tmp1 - Tmp2)*dCoef);
      WM.DecCoef(5, 3, Af*Tmp1*dCoef);
   } else {
      doublereal Tmp1 = (Ise/Vt)*exp(-Vbe/Vt);
      doublereal Tmp2 = (Isc/Vt)*exp(-Vbc/Vt);

      WM.DecCoef(4, 1, Ar*Tmp2*dCoef);
      WM.DecCoef(4, 2, (Tmp1 - Ar*Tmp2)*dCoef);
      WM.IncCoef(4, 3, Tmp1*dCoef);

      WM.DecCoef(5, 1, Tmp2*dCoef);
      WM.DecCoef(5, 2, (Af*Tmp1 - Tmp2)*dCoef);
      WM.IncCoef(5, 3, Af*Tmp1*dCoef);
   }

   WM.DecCoef(4, 6, 1.);
   WM.DecCoef(5, 4, 1.);

   WM.DecCoef(6, 4, 1.);
   WM.DecCoef(6, 5, 1.);
   WM.IncCoef(6, 6, 1.);

	return WorkMat;
}

SubVectorHandler&
BipolarTransistor::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);


   integer iElecNodeFirstIndexC = pElecC->iGetFirstRowIndex() + 1;
   integer iElecNodeFirstIndexB = pElecB->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndexE = pElecE->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   WorkVec.PutRowIndex(1, iElecNodeFirstIndexC);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndexB);
   WorkVec.PutRowIndex(3, iElecNodeFirstIndexE);
   WorkVec.PutRowIndex(4, iFirstIndex);
   WorkVec.PutRowIndex(5, iFirstIndex + 1);
   WorkVec.PutRowIndex(6, iFirstIndex + 2);

   doublereal Af = Bf/(Bf+1.);
   doublereal Ar = Br/(Br+1.);
   doublereal Vbe = pElecB->dGetX() - pElecE->dGetX();
   doublereal Vbc = pElecB->dGetX() - pElecC->dGetX();

   doublereal iC = XCurr(iFirstIndex);
   doublereal iB = XCurr(iFirstIndex + 1);
   doublereal iE = XCurr(iFirstIndex + 2);

   doublereal Tmp1(0.);
   doublereal Tmp2(0.);

   if (bjt_type==NPN) {
      Tmp1 = Ise*(exp(Vbe/Vt) - 1.);
      Tmp2 = Isc*(exp(Vbc/Vt) - 1.);
   } else {
      Tmp1 = Ise*(exp(-Vbe/Vt) - 1.);
      Tmp2 = Isc*(exp(-Vbc/Vt) - 1.);
   }

   WorkVec.IncCoef(1, iC);
   WorkVec.IncCoef(2, iB);
   WorkVec.DecCoef(3, iE);
   WorkVec.IncCoef(4, iE - Tmp1 + Ar*Tmp2);
   WorkVec.IncCoef(5, iC - Af*Tmp1 + Tmp2);
   WorkVec.IncCoef(6, iB - iE + iC);

	return WorkVec;
}

unsigned int
BipolarTransistor::iGetNumPrivData(void) const
{
	return 0;
}

int
BipolarTransistor::iGetNumConnectedNodes(void) const
{
   return 3;
}

void
BipolarTransistor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
   connectedNodes.resize(4);
   connectedNodes[0] = pElecC;
   connectedNodes[1] = pElecB;
   connectedNodes[2] = pElecE;
}

void
BipolarTransistor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
BipolarTransistor::Restart(std::ostream& out) const
{
	return out << "# BipolarTransistor: not implemented" << std::endl;
}

unsigned int
BipolarTransistor::iGetInitialNumDof(void) const
{
	return 0;
}

void
BipolarTransistor::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
BipolarTransistor::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
BipolarTransistor::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



// Proximity Sensor:
// Model not tested!!!!!

ProximitySensor::ProximitySensor(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	ProximitySensor			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}


	// Read the node 1 from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
          silent_cerr("Proximity Sensor (" << GetLabel() << ") - node 1: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   X1tilde = Vec3(0., 0., 0.);

	// Read the position offset of node 1, if supplied:
   if (HP.IsKeyWord("position")) {
        DEBUGCOUT("Position offset of node 1 is supplied" << std::endl);
        X1tilde = HP.GetVec3();
   }

	// Read the node 2 from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
          silent_cerr("Proximity Sensor (" << GetLabel() << ") - node 2: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   X2tilde = Vec3(0., 0., 0.);

	// Read the position offset of node 2, if supplied:
   if (HP.IsKeyWord("position")) {
        DEBUGCOUT("Position offset of node 2 is supplied" << std::endl);
        X2tilde = HP.GetVec3();
   }


	// Read electrical node 1:
   pElec1 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec1) {
      silent_cerr("ProximitySensor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read electrical node 2:
   pElec2 = pDM->ReadNode<const ElectricNode, Node::ELECTRIC>(HP);
   if (!pElec2) {
      silent_cerr("ProximitySensor (" << GetLabel() << "): electric node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read scalar function that relates distance to tension or to current:
   sFun = dynamic_cast<const DifferentiableScalarFunction*> (ParseScalarFunction(HP, pDM));

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

ProximitySensor::~ProximitySensor(void)
{
	// destroy private data
	NO_OP;
}

void
ProximitySensor::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << i_curr        // current on ProximitySensor
      << " " << Voltage1      // voltage on node 1
      << " " << Voltage2      // voltage on node 2
      << std::endl;
   }

}

unsigned int
ProximitySensor::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
ProximitySensor::GetDofType(unsigned int i) const
{
	return  DofOrder::ALGEBRAIC;
}

void
ProximitySensor::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   // Get eletric current:
   i_curr = X(iGetFirstIndex()+1);
   Voltage1 = pElec1->dGetX();
   Voltage2 = pElec2->dGetX();

}

void
ProximitySensor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3 + 12;
}

VariableSubMatrixHandler&
ProximitySensor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iStrNode1FirstPosIdx = pNode1->iGetFirstPositionIndex();
	integer iStrNode2FirstPosIdx = pNode2->iGetFirstPositionIndex();
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutRowIndex(1, iElecNodeFirstIndex1);
	WM.PutRowIndex(2, iElecNodeFirstIndex2);
	WM.PutRowIndex(3, iFirstIndex);

	WM.PutColIndex(1, iElecNodeFirstIndex1);
	WM.PutColIndex(2, iElecNodeFirstIndex2);
	WM.PutColIndex(3, iFirstIndex);
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
	   WM.PutColIndex(3 + iCnt, iStrNode1FirstPosIdx + iCnt);
	   WM.PutColIndex(9 + iCnt, iStrNode2FirstPosIdx + iCnt);
   }

   // Get info from nodes:
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 X1(pNode1->GetXCurr());
   Vec3 X2(pNode2->GetXCurr());
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();

   // Calculate relative position and velocities:
   Vec3 dX = R2*X2tilde + X2 - R1*X1tilde - X1;
	doublereal i = XCurr(iFirstIndex);
   doublereal dist = dX.Norm();
   doublereal dFdR = sFun->ComputeDiff(dist, 1);

   WM.IncCoef(1, 3, 1.);
   WM.DecCoef(2, 3, 1.);

   WM.IncCoef(3, 1, dCoef);
   WM.DecCoef(3, 2, dCoef);

   Vec3 dRdg1 = Mat3x3(MatCross, R1*X1tilde).MulTV(dX);
   Vec3 dRdg2 = Mat3x3(MatCross, -R2*X2tilde).MulTV(dX);

   doublereal tmp1 = dFdR/dist*dCoef;

   WM.AddT(3, 4, -dX*tmp1);
   WM.AddT(3, 4 + 3, dRdg1*tmp1);
   WM.AddT(3, 4 + 6, dX*tmp1);
   WM.AddT(3, 4 + 9, dRdg2*tmp1);

	return WorkMat;
}

SubVectorHandler&
ProximitySensor::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Resize and reset the working vector */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   // Recover indexes:
	integer iElecNodeFirstIndex1 = pElec1->iGetFirstRowIndex() + 1;
	integer iElecNodeFirstIndex2 = pElec2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

   // Allocate rows in the WorkVec:
   WorkVec.PutRowIndex(1, iElecNodeFirstIndex1);
   WorkVec.PutRowIndex(2, iElecNodeFirstIndex2);
   WorkVec.PutRowIndex(3, iFirstIndex);

   // Get info from nodes:
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 X1(pNode1->GetXCurr());
   Vec3 X2(pNode2->GetXCurr());
	doublereal V1 = pElec1->dGetX();
	doublereal V2 = pElec2->dGetX();

   // Calculate relative position and velocities:
   Vec3 dX = R2*X2tilde + X2 - R1*X1tilde - X1;
	doublereal i = XCurr(iFirstIndex);

   WorkVec.DecCoef(1, i);
   WorkVec.IncCoef(2, i);
   WorkVec.IncCoef(3, V2 - V1 - (*sFun)(dX.Norm()));

	return WorkVec;
}

unsigned int
ProximitySensor::iGetNumPrivData(void) const
{
	return 0;
}

int
ProximitySensor::iGetNumConnectedNodes(void) const
{
	return 4;
}

void
ProximitySensor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(4);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
	connectedNodes[2] = pElec1;
	connectedNodes[3] = pElec2;
}

void
ProximitySensor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ProximitySensor::Restart(std::ostream& out) const
{
	return out << "# ProximitySensor: not implemented" << std::endl;
}

unsigned int
ProximitySensor::iGetInitialNumDof(void) const
{
	return 0;
}

void
ProximitySensor::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ProximitySensor::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
ProximitySensor::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}



extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<Resistor>;

	if (!SetUDE("resistor", rf1)) {
		delete rf1;

		silent_cerr("resistor: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf2 = new UDERead<Capacitor>;

	if (!SetUDE("capacitor", rf2)) {
		delete rf2;

		silent_cerr("capacitor: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf3 = new UDERead<Inductor>;

	if (!SetUDE("inductor", rf3)) {
		delete rf3;

		silent_cerr("inductor: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf4 = new UDERead<Diode>;

	if (!SetUDE("diode", rf4)) {
		delete rf4;

		silent_cerr("diode: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf5 = new UDERead<Switch>;

	if (!SetUDE("switch", rf5)) {
		delete rf5;

		silent_cerr("switch: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf6 = new UDERead<ElectricalSource>;

	if (!SetUDE("electrical" "source", rf6)) {
		delete rf6;

		silent_cerr("electrical source: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf7 = new UDERead<IdealTransformer>;

	if (!SetUDE("ideal" "transformer", rf7)) {
		delete rf7;

		silent_cerr("ideal transformer: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf8 = new UDERead<OperationalAmplifier>;

	if (!SetUDE("operational" "amplifier", rf8)) {
		delete rf8;

		silent_cerr("operational amplifier: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf9 = new UDERead<BipolarTransistor>;

	if (!SetUDE("bjt", rf9)) {
		delete rf9;

		silent_cerr("bjt: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf10 = new UDERead<ProximitySensor>;

	if (!SetUDE("proximity" "sensor", rf10)) {
		delete rf10;

		silent_cerr("proximity sensor: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}



