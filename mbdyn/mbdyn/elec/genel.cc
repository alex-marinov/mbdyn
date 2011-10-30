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

/* Genel */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "genel_.h"
#include "genfilt.h"
#include "swashpl.h"
#include "rottrim.h"
#include "dataman.h"

/* genel - begin */

Genel::Genel(unsigned int uL,
	const DofOwner* pDO,
	flag fOut)
: Elem(uL, fOut),
ElemWithDofs(uL, pDO, fOut)
{
	NO_OP;
}

Genel::~Genel(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Genel::Restart(std::ostream& out) const
{
	return out << "genel: " << GetLabel();
}

/* Genel - end */


/* Legge un genel */

Elem*
ReadGenel(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadGenel()");

	const char* sKeyWords[] = {
		"swashplate",
		"rotor" "trim",
		"clamp",
		"distance",
		"spring",
		"spring" "support",
		"cross" "spring" "support",
		"mass",
		"scalar" "filter",
		"state" "space" "SISO",
		"state" "space" "MIMO",

		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		SWASHPLATE = 0,
		ROTORTRIM,
		CLAMP,
		DISTANCE,
		SPRING,
		SPRINGSUPPORT,
		CROSSSPRINGSUPPORT,
		MASS,
		SCALARFILTER,
		STATESPACESISO,
		STATESPACEMIMO,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di vincolo */
	KeyWords CurrKeyWord = KeyWords(HP.GetWord());

	Elem* pEl = 0;

	switch (CurrKeyWord) {
	/* genel piatto oscillante */
	case SWASHPLATE: {
		/* nodo Collettivo */
		const ScalarDifferentialNode* pCollIn = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pCollIn != 0);

		flag fCollLimits(0);
		doublereal dCollMax(0.);
		doublereal dCollMin(0.);
		if (HP.IsKeyWord("limits")) {
			fCollLimits = flag(1);
			dCollMin = HP.GetReal();
			dCollMax = HP.GetReal();
		}

		/* nodo Longitudinale */
		const ScalarDifferentialNode* pLongIn = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pLongIn != 0);

		flag fForeAftLimits(0);
		doublereal dForeAftMax(0.);
		doublereal dForeAftMin(0.);
		if (HP.IsKeyWord("limits")) {
			fForeAftLimits = flag(1);
			dForeAftMin = HP.GetReal();
			dForeAftMax = HP.GetReal();
		}

		/* nodo Laterale */
		const ScalarDifferentialNode* pLatIn = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pLatIn != 0);

		flag fLatLimits(0);
		doublereal dLatMax(0.);
		doublereal dLatMin(0.);
		if (HP.IsKeyWord("limits")) {
			fLatLimits = flag(1);
			dLatMin = HP.GetReal();
			dLatMax = HP.GetReal();
		}

		/* nodo collegato 1 */
		const ScalarDifferentialNode* pNode1 = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pNode1 != 0);

		/* nodo collegato 2 */
		const ScalarDifferentialNode* pNode2 = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pNode2 != 0);

		/* nodo collegato 3 */
		const ScalarDifferentialNode* pNode3 = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		ASSERT(pNode3 != 0);

		doublereal dDynCoef = 0.;
		if (HP.IsArg()) {
			dDynCoef = HP.GetReal(dDynCoef);
		}

		doublereal dCyclFact = 1.;
		if (HP.IsArg()) {
			dCyclFact = HP.GetReal(dCyclFact);
		}

		doublereal dCollFact = 1.;
		if (HP.IsArg()) {
			dCollFact = HP.GetReal(dCollFact);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			SwashPlate,
			SwashPlate(uLabel, pDO,
				pCollIn,
				pLongIn,
				pLatIn,
				pNode1, pNode2, pNode3,
				dDynCoef,
				dCyclFact,
				dCollFact,
				fCollLimits,
				dCollMin,
				dCollMax,
				fForeAftLimits,
				dForeAftMin,
				dForeAftMax,
				fLatLimits,
				dLatMin,
				dLatMax,
				fOut));
		} break;

	case ROTORTRIM: {
		// "rotor" (traditional)
		const Rotor* pRot = 0;

		// "generic" 
		const StructNode *pStrNode = 0;
		const DriveCaller *pThrust = 0;
		const DriveCaller *pRollMoment = 0;
		const DriveCaller *pPitchMoment = 0;
		const AirProperties *pAP = 0;
		doublereal dRadius = -1.;
		const DriveCaller *pOmega = 0;
		const DriveCaller *pMu = 0;

		if (HP.IsKeyWord("generic")) {
			if (HP.IsKeyWord("reference" "node")) {
				pStrNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
			}

			pThrust = HP.GetDriveCaller();
			pRollMoment = HP.GetDriveCaller();
			pPitchMoment = HP.GetDriveCaller();

			pAP = dynamic_cast<const AirProperties *>(pDM->pFindElem(Elem::AIRPROPERTIES, 1));

			dRadius = HP.GetReal();
			if (dRadius < std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("RotorTrim(" << uLabel << "): "
					"invalid rotor radius at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			pOmega = HP.GetDriveCaller();
			pMu = HP.GetDriveCaller();

		} else {
			if (!HP.IsKeyWord("rotor")) {
				silent_cout("RotorTrim(" << uLabel << "): "
					"keyword \"rotor\" expected at line " << HP.GetLineData() << std::endl);
			}

			pRot = pDM->ReadElem<const InducedVelocity, Elem::INDUCEDVELOCITY, const Rotor>(HP);
			if (pRot == 0) {
				silent_cerr("RotorTrim(" << uLabel << "): "
					"unable to read rotor "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dRadius = pRot->dGetRadius();
			if (dRadius < std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("RotorTrim(" << uLabel << "): "
					"invalid rotor radius for Rotor(" << pRot->GetLabel() << ") "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		const ScalarDifferentialNode* pvNodes[3];
		pvNodes[0] = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		pvNodes[1] = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		pvNodes[2] = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);

		DEBUGCOUT("Rotor trim " << uLabel
			<< " linked to rotor " << pRot->GetLabel() << std::endl
			<< "abstract nodes: "
			<< pvNodes[0]->GetLabel() << ", "
			<< pvNodes[1]->GetLabel() << ", "
			<< pvNodes[2]->GetLabel() << std::endl);

		DriveCaller* pvDrives[3];
		pvDrives[0] = HP.GetDriveCaller();
		pvDrives[1] = HP.GetDriveCaller();
		pvDrives[2] = HP.GetDriveCaller();

		doublereal dGamma = HP.GetReal();
		DEBUGCOUT("Gamma: " << dGamma << std::endl);

		doublereal dP = HP.GetReal();
		DEBUGCOUT("P: " << dP << std::endl);

		doublereal dTau0 = HP.GetReal();
		DEBUGCOUT("Tau0: " << dTau0 << std::endl);

		doublereal dTau1 = HP.GetReal();
		DEBUGCOUT("Tau1: " << dTau1 << std::endl);

		doublereal dKappa0 = HP.GetReal();
		DEBUGCOUT("Kappa0: " << dKappa0 << std::endl);

		doublereal dKappa1 = HP.GetReal();
		DEBUGCOUT("Kappa1: " << dKappa1 << std::endl);

		DriveCaller *pTrigger = 0;
		if (HP.IsKeyWord("trigger")) {
			pTrigger = HP.GetDriveCaller();

		} else {
			SAFENEW(pTrigger, OneDriveCaller);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		if (pRot) {
			SAFENEWWITHCONSTRUCTOR(pEl,
				RotorTrim,
				RotorTrim(uLabel, pDO, pRot,
					pvNodes[0], pvNodes[1], pvNodes[2],
					pvDrives[0], pvDrives[1], pvDrives[2],
					dGamma, dP,
					dTau0, dTau1, dKappa0, dKappa1,
					pTrigger, fOut));

		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
				RotorTrimGeneric,
				RotorTrimGeneric(uLabel, pDO,
					pStrNode,
					pThrust, pRollMoment, pPitchMoment,
					pAP, dRadius, pOmega, pMu,
					pvNodes[0], pvNodes[1], pvNodes[2],
					pvDrives[0], pvDrives[1], pvDrives[2],
					dGamma, dP,
					dTau0, dTau1, dKappa0, dKappa1,
					pTrigger, fOut));
		}
		} break;

	case CLAMP: {
		ScalarDof SD = ReadScalarDof(pDM, HP, true, true);
		if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelClamp(" << uLabel << "): "
				"parameter nodes are not allowed "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD.iOrder > 1) {
			silent_cerr("GenelClamp(" << uLabel << "): "
				"illegal order " << SD.iOrder
				<< " for ScalarDof "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		DriveCaller* pDC = HP.GetDriveCaller();

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelClamp,
			GenelClamp(uLabel, pDO, pDC, SD, fOut));
		} break;

	case DISTANCE: {
		ScalarDof SD1 = ReadScalarDof(pDM, HP, true, true);
		if (SD1.pNode->GetNodeType() == Node::PARAMETER) {
			silent_cerr("GenelDistance(" << uLabel << "): "
				"parameter nodes not allowed "
				"for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD1.iOrder > 1) {
			silent_cerr("GenelDistance(" << uLabel << "): "
				"illegal order " << SD1.iOrder
				<< " for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarDof SD2 = ReadScalarDof(pDM, HP, true, true);
		if (SD2.pNode->GetNodeType() == Node::PARAMETER) {
			silent_cerr("GenelDistance(" << uLabel << "): "
				"parameter nodes not allowed "
				"for ScalarDof 2 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD2.iOrder > 1) {
			silent_cerr("GenelDistance(" << uLabel << "): "
				"illegal order " << SD2.iOrder
				<< " for ScalarDof 2 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		DriveCaller* pDC = HP.GetDriveCaller();

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelDistance,
			GenelDistance(uLabel, pDO, pDC, SD1, SD2, fOut));

		} break;

	case SPRING: {
		ScalarDof SD1 = ReadScalarDof(pDM, HP, true, true);
		if (SD1.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD1.iOrder > 1) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"illegal order " << SD1.iOrder
				<< " for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarDof SD2 = ReadScalarDof(pDM, HP, true, true);
		if (SD2.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof 2 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD2.iOrder > 1) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"illegal order " << SD2.iOrder
				<< " for ScalarDof 2 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ConstLawType::Type CLType = ConstLawType::UNKNOWN;
		ConstitutiveLaw1D* pCL = HP.GetConstLaw1D(CLType);

		if (pCL->iGetNumDof() != 0) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"support dynamic constitutive laws "
				"not supported "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (CLType != ConstLawType::ELASTIC) {
			silent_cerr("GenelSpring(" << uLabel << "): "
				"only elastic constitutive laws allowed "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelSpring,
			GenelSpring(uLabel, pDO, pCL, SD1, SD2, fOut));

		} break;

	case SPRINGSUPPORT: {
		ScalarDof SD = ReadScalarDof(pDM, HP, true, true);
		if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelSpringSupport(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD.iOrder != 0 ||
			SD.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL)
		{
			silent_cerr("GenelSpringSupport(" << uLabel << "): "
				"illegal order " << SD.iOrder
				<< " for ScalarDof; the algebraic value "
				"of a differential node is required "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ConstLawType::Type CLType = ConstLawType::UNKNOWN;
		ConstitutiveLaw1D* pCL = HP.GetConstLaw1D(CLType);
		if (pCL->iGetNumDof() != 0) {
			silent_cerr("GenelSpringSupport(" << uLabel << "): "
				"only elastic constitutive laws allowed "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		switch (CLType) {
		case ConstLawType::ELASTIC:
			SAFENEWWITHCONSTRUCTOR(pEl,
				GenelSpringSupport,
				GenelSpringSupport(uLabel, pDO, pCL,
					SD, fOut));
			break;

		case ConstLawType::VISCOUS:
		case ConstLawType::VISCOELASTIC:
			SAFENEWWITHCONSTRUCTOR(pEl,
				GenelSpringDamperSupport,
				GenelSpringDamperSupport(uLabel, pDO, pCL,
					SD, fOut));
			break;

		default:
			silent_cerr("You shouldn't be here!" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		} break;

	case CROSSSPRINGSUPPORT: {
		ScalarDof SDRow = ReadScalarDof(pDM, HP, true, true);
		if (SDRow.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr(
				"GenelCrossSpringSupport(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SDRow.iOrder > 1) {
			silent_cerr(
				"GenelCrossSpringSupport(" << uLabel << "): "
				"illegal order " << SDRow.iOrder
				<< " for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarDof SDCol = ReadScalarDof(pDM, HP, true, true);
		if (SDCol.iOrder != 0 ||
			SDCol.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL)
		{
			silent_cerr(
				"GenelCrossSpringSupport(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SDCol.iOrder != 0 ||
			SDCol.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL)
		{
			silent_cerr(
				"GenelCrossSpringSupport(" << uLabel << "): "
				"illegal order " << SDCol.iOrder
				<< " for ScalarDof; the algebraic value "
				"of a differential node is required "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ConstLawType::Type CLType = ConstLawType::UNKNOWN;
		ConstitutiveLaw1D* pCL = HP.GetConstLaw1D(CLType);

		if (pCL->iGetNumDof() != 0) {
			silent_cerr(
				"GenelCrossSpringSupport(" << uLabel << "): "
				"dynamic constitutive laws not supported "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		switch (CLType) {
		case ConstLawType::ELASTIC:
			SAFENEWWITHCONSTRUCTOR(pEl,
				GenelCrossSpringSupport,
				GenelCrossSpringSupport(uLabel, pDO, pCL,
					SDRow, SDCol, fOut));
			break;

		case ConstLawType::VISCOUS:
		case ConstLawType::VISCOELASTIC:
			SAFENEWWITHCONSTRUCTOR(pEl,
				GenelCrossSpringDamperSupport,
				GenelCrossSpringDamperSupport(uLabel, pDO, pCL,
					SDRow, SDCol, fOut));
			break;

		default:
			silent_cerr("You shouldn't be here!" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		} break;

	case MASS: {
		ScalarDof SD = ReadScalarDof(pDM, HP, true, true);
		if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelMass(" << uLabel << "): "
				"parameter nodes not allowed for ScalarDof "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD.iOrder > 1) {
			silent_cerr("GenelMass(" << uLabel << "): "
				"illegal order " << SD.iOrder
				<< " for ScalarDof 1 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL) {
			silent_cerr("GenelMass(" << uLabel << "): "
				"only differential dofs allowed for ScalarDof "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		DriveCaller* pDC = HP.GetDriveCaller();

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelMass,
			GenelMass(uLabel, pDO, pDC, SD, fOut));

		} break;

	case STATESPACESISO: {
		ScalarDof SD_y = ReadScalarDof(pDM, HP, true, true);
		if (SD_y.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
				"parameter nodes not allowed "
				"for output ScalarDof (y) "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD_y.iOrder > 1) {
			silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
				"illegal order " << SD_y.iOrder
				<< " for output ScalarDof (y) "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarValue *SV_u = ReadScalarValue(pDM, HP);

		unsigned int Order = HP.GetInt();
		DEBUGCOUT("State Space SISO " << uLabel
			<< " is of order " << Order << std::endl);

		doublereal* pd = 0;

		doublereal* pdE = 0;
		if (HP.IsKeyWord("matrix" "E")) {
			SAFENEWARR(pdE, doublereal, Order*Order);
			pd = pdE;
			for (unsigned int i = 0; i < Order*Order; i++) {
				*pd++ = HP.GetReal();
			}
		}

		if (!HP.IsKeyWord("matrix" "A")) {
			silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
				"matrix A expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdA = 0;
		SAFENEWARR(pdA, doublereal, Order*Order);
		pd = pdA;
		for (unsigned int i = 0; i < Order*Order; i++) {
			*pd++ = HP.GetReal();
		}

		if (!HP.IsKeyWord("matrix" "B")) {
			silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
				"matrix B expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdB = 0;
		SAFENEWARR(pdB, doublereal, Order);
		pd = pdB;
		for (unsigned int i = 0; i < Order; i++) {
			*pd++ = HP.GetReal();
		}

		if (!HP.IsKeyWord("matrix" "C")) {
			silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
				"matrix C expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdC = 0;
		SAFENEWARR(pdC, doublereal, Order);
		pd = pdC;
		for (unsigned int i = 0; i < Order; i++) {
			*pd++ = HP.GetReal();
		}

		doublereal dD = 0.;
		if (HP.IsKeyWord("matrix" "D")) {
			if (pdE) {
				silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
					"warning, matrix D provided "
					"while in descriptor form "
					"at line " << HP.GetLineData()
					<< std::endl);
			}
			dD = HP.GetReal();
		}

		if (HP.IsKeyWord("gain")) {
			doublereal dGain = HP.GetReal();

			pd = pdC;
			for (unsigned int i = 0; i < Order; i++) {
				*pd++ *= dGain;
			}

			dD *= dGain;
		}

		bool bBalance(true);
		if (HP.IsKeyWord("balance")) {
			if (!HP.GetYesNo(bBalance)) {
				silent_cerr("GenelStateSpaceSISO(" << uLabel << "): "
					"unknown balance mode at line "
					<< HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal *pdX0 = 0;
		std::vector<doublereal> dX0;
		doublereal *pdXP0 = 0;
		std::vector<doublereal> dXP0;

		if (HP.IsKeyWord("value")) {
			dX0.resize(Order);
			pdX0 = &dX0[0];

			for (unsigned i = 0; i < Order; i++) {
				dX0[i] = HP.GetReal();
			}

			if (pdE != 0 && HP.IsKeyWord("derivative")) {
				dXP0.resize(Order);
				pdXP0 = &dXP0[0];

				for (unsigned i = 0; i < Order; i++) {
					dXP0[i] = HP.GetReal();
				}
			}
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelStateSpaceSISO,
			GenelStateSpaceSISO(uLabel, pDO, SD_y, SV_u,
				Order,
				pdE, pdA, pdB, pdC, dD, bBalance,
				pdX0, pdXP0, fOut));

		} break;

	case SCALARFILTER: {
		ScalarDof SD_y = ReadScalarDof(pDM, HP, true, true);
		if (SD_y.pNode->GetNodeType() ==  Node::PARAMETER) {
			silent_cerr("ScalarFilter(" << uLabel << "): "
				"parameter nodes not allowed "
				"for output ScalarDof (y) "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (SD_y.iOrder > 1) {
			silent_cerr("ScalarFilter(" << uLabel << "): "
				"illegal order " << SD_y.iOrder
				<< " for output ScalarDof (y) "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarValue *SV_u = ReadScalarValue(pDM, HP);

/*
 * if nn < nd (strictly proper):
 *
 *	n_0 s^nn + n_1 s^(nn-1) + ... + n_nn
 *	------------------------------------
 *	  s^nd + d_1 s^(nd-1) + ... + d_nd
 *
 * d_0 == 1
 *
 * if nn == nd (proper, not strictly proper):
 *
 *	n_0 s^n + n_1 s^(n-1) + ... + n_n
 *	------------------------------------
 *	  s^n + d_1 s^(n-1) + ... + d_n
 *
 *	      (n_1 - d_1*n_0) s^(n-1) + ... (n_n - d_n *n_0)
 *	n_0 + ----------------------------------------------
 *	               s^n + d_1 s^(n-1) + ... + d_n
 *
 * d_0 == 1
 *
 */

		bool bControllable(true);
		if (HP.IsKeyWord("canonical" "form")) {
			if (HP.IsKeyWord("controllable")) {
				bControllable = true;

			} else if (HP.IsKeyWord("observable")) {
				bControllable = false;

			} else {
				silent_cerr("ScalarFilter(" << uLabel << "):"
					"unknown canonical form "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		int nd = HP.GetInt();
		if (nd < 1) {
			silent_cerr("ScalarFilter(" << uLabel << "):"
				"invalid denominator order " << nd << " "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::vector<doublereal> den(nd + 1);
		den[0] = 1.;
		for (int i = 1; i <= nd; i++) {
			den[i] = HP.GetReal();
		}

		int nn = HP.GetInt();
		if (nn < 0 || nn > nd) {
			silent_cerr("ScalarFilter(" << uLabel << "):"
				"invalid numerator order " << nn << " "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::vector<doublereal> num(nd + 1);
		for (int i = 0; i < nd - nn; i++) {
			num[i] = 0.;
		}

		for (int i = nd - nn; i <= nd; i++) {
			num[i] = HP.GetReal();
		}

		if (HP.IsKeyWord("gain")) {
			doublereal dGain = HP.GetReal();

			for (int i = 0; i <= nd; i++) {
				num[i] *= dGain;
			}
		}

		// turn into strictly proper plus direct feedthrough
		if (nn == nd) {
			for (int i = 1; i <= nd; i++) {
				num[i] -= den[i]*num[0];
			}
		}

		int Order = nd;

		doublereal* pd = 0;

		doublereal* pdE = 0;
		doublereal* pdA = 0;
		SAFENEWARR(pdA, doublereal, Order*Order);
		doublereal* pdB = 0;
		SAFENEWARR(pdB, doublereal, Order);
		doublereal* pdC = 0;
		SAFENEWARR(pdC, doublereal, Order);
		doublereal dD = num[0];

		if (bControllable) {
			// matrix A
			pd = pdA;
			for (int i = 1; i <= Order; i++) {
				*pd++ = -den[i];
			}
			for (int i = 0; i < Order - 2; i++) {
				*pd++ = 1.;
				for (int j = 0; j < Order; j++) {
					*pd++ = 0.;
				}
			}

			if (Order > 1) {
				*pd++ = 1.;
				*pd++ = 0.;
			}

			ASSERT(pd == pdA + Order*Order);
	
			// matrix B
			pd = pdB;
			*pd++ = 1.;
			for (int i = 0; i < Order - 1; i++) {
				*pd++ = 0.;
			}
	
			// matrix C
			pd = pdC;
			for (int i = 1; i <= Order; i++) {
				*pd++ = num[i];
			}

		} else {
			// matrix A
			pd = pdA;
			for (int i = 1; i <= Order; i++) {
				*pd++ = -den[i];
				for (int j = 1; j < Order; j++) {
					if (j == i) {
						*pd++ = 1.;

					} else {
						*pd++ = 0.;
					}
				}
			}

			ASSERT(pd == pdA + Order*Order);

			// matrix B
			pd = pdB;
			for (int i = 1; i <= Order; i++) {
				*pd++ = num[i];
			}
	
			// matrix C
			pd = pdC;
			*pd++ = 1.;
			for (int i = 0; i < Order - 1; i++) {
				*pd++ = 0.;
			}
		}
	
		bool bBalance(true);
		if (HP.IsKeyWord("balance")) {
			if (!HP.GetYesNo(bBalance)) {
				silent_cerr("ScalarFilter(" << uLabel << "): "
					"unknown balance mode at line "
					<< HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelStateSpaceSISO,
			GenelStateSpaceSISO(uLabel, pDO, SD_y, SV_u,
				Order,
				pdE, pdA, pdB, pdC, dD, bBalance,
				0, 0, fOut));

		} break;

	case STATESPACEMIMO: {
		int iNumOutputs = HP.GetInt();
		if (iNumOutputs <= 0) {
			silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
				"illegal number of outputs "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ScalarDof* pvSD_y = 0;
		SAFENEWARRNOFILL(pvSD_y, ScalarDof, iNumOutputs);
		for (int i = 0; i < iNumOutputs; i++) {
			pvSD_y[i] = ReadScalarDof(pDM, HP, true, true);
			if (pvSD_y[i].pNode->GetNodeType() ==  Node::PARAMETER)
			{
				silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
					"parameter nodes not allowed "
					"for output ScalarDof "
					"(y[" << i << "]) "
					"at line " << HP.GetLineData()
					<< std::endl);
				SAFEDELETEARR(pvSD_y);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (pvSD_y[i].iOrder > 1) {
				silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
					"illegal order " << pvSD_y[i].iOrder
					<< " for output ScalarDof "
					"(y[" << i << "]) "
					"at line " << HP.GetLineData()
					<< std::endl);
				SAFEDELETEARR(pvSD_y);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		int iNumInputs = HP.GetInt();
		if (iNumInputs <= 0) {
			silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
				"illegal number of inputs "
				"at line " << HP.GetLineData()
				<< std::endl);
			SAFEDELETEARR(pvSD_y);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::vector<ScalarValue *> SV_u(iNumInputs);
		ReadScalarValues(pDM, HP, SV_u);

		unsigned int Order = HP.GetInt();
		DEBUGCOUT("State Space MIMO " << uLabel
			<< " is of order " << Order << std::endl);

		doublereal* pd = 0;
		doublereal* pdE = 0;
		if (HP.IsKeyWord("matrix" "E")) {
			SAFENEWARR(pdE, doublereal, Order*Order);
			pd = pdE;
			for (unsigned int i = 0; i < Order*Order; i++) {
				*pd++ = HP.GetReal();
			}
		}

		if (!HP.IsKeyWord("matrix" "A")) {
			silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
				"matrix A expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			SAFEDELETEARR(pvSD_y);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdA = 0;
		SAFENEWARR(pdA, doublereal, Order*Order);
		pd = pdA;
		for (unsigned int i = 0; i < Order*Order; i++) {
			*pd++ = HP.GetReal();
		}

		if (!HP.IsKeyWord("matrix" "B")) {
			silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
				"matrix B expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			SAFEDELETEARR(pvSD_y);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdB = 0;
		SAFENEWARR(pdB, doublereal, Order*iNumInputs);
		pd = pdB;
		for (unsigned int i = 0; i < Order*iNumInputs; i++) {
			*pd++ = HP.GetReal();
		}

		if (!HP.IsKeyWord("matrix" "C")) {
			silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
				"matrix C expected "
				"at line " << HP.GetLineNumber()
				<< std::endl);
			SAFEDELETEARR(pvSD_y);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal* pdC = 0;
		SAFENEWARR(pdC, doublereal, iNumOutputs*Order);
		pd = pdC;
		for (unsigned int i = 0; i < iNumOutputs*Order; i++) {
			*pd++ = HP.GetReal();
		}

		doublereal* pdD = 0;
		if (HP.IsKeyWord("matrix" "D")) {
			if (pdE) {
				silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
					"warning, matrix D provided "
					"while in descriptor form "
					"at line " << HP.GetLineData()
					<< std::endl);
				SAFEDELETEARR(pvSD_y);
			}
			SAFENEWARR(pdD, doublereal, iNumOutputs*iNumInputs);
			pd = pdD;
			for (int i = 0; i < iNumOutputs*iNumInputs; i++) {
				*pd++ = HP.GetReal();
			}
		}

		if (HP.IsKeyWord("gain")) {
			doublereal dGain = HP.GetReal();

			pd = pdC;
			for (unsigned int i = 0; i < iNumOutputs*Order; i++) {
				*pd++ *= dGain;
			}

			pd = pdD;
			for (int i = 0; i < iNumOutputs*iNumInputs; i++) {
				*pd++ *= dGain;
			}
		}

		bool bBalance(true);
		if (HP.IsKeyWord("balance")) {
			if (!HP.GetYesNo(bBalance)) {
				silent_cerr("GenelStateSpaceMIMO(" << uLabel << "): "
					"unknown balance mode at line "
					<< HP.GetLineData()
					<< std::endl);
				SAFEDELETEARR(pvSD_y);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal *pdX0 = 0;
		std::vector<doublereal> dX0;
		doublereal *pdXP0 = 0;
		std::vector<doublereal> dXP0;

		if (HP.IsKeyWord("value")) {
			dX0.resize(Order);
			pdX0 = &dX0[0];

			for (unsigned i = 0; i < Order; i++) {
				dX0[i] = HP.GetReal();
			}

			if (pdE != 0 && HP.IsKeyWord("derivative")) {
				dXP0.resize(Order);
				pdXP0 = &dXP0[0];

				for (unsigned i = 0; i < Order; i++) {
					dXP0[i] = HP.GetReal();
				}
			}
		}

		flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

		SAFENEWWITHCONSTRUCTOR(pEl,
			GenelStateSpaceMIMO,
			GenelStateSpaceMIMO(uLabel, pDO,
				iNumOutputs, pvSD_y,
				SV_u,
				Order,
				pdE, pdA, pdB, pdC, pdD, bBalance,
				pdX0, pdXP0, fOut));

		} break;

	/* Aggiungere altri genel */
	default:
		silent_cerr("Genel(" << uLabel << "): "
			"unknown type at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* ReadGenel() */

