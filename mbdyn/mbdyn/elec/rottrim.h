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

#ifndef ROTTRIM_H
#define ROTTRIM_H

#include "rotor.h"
#include "genel.h"

/*
	pRotor->iGetPrivDataIdx("Tz");
	pRotor->iGetPrivDataIdx("Mx");
	pRotor->iGetPrivDataIdx("My");
	doublereal dRho = pRotor->dGetAirDensity(pRotor->GetXCurr());
	doublereal dOmega = pRotor->dGetOmega();
	doublereal dRadius = pRotor->dGetRadius();
	doublereal dMu = pRotor->dGetMu();
	doublereal dThrust = pRotor->dGetPrivData(iRotorTz)/d;
	doublereal dRollMoment = pRotor->dGetPrivData(iRotorMx)/d;
	doublereal dPitchMoment = pRotor->dGetPrivData(iRotorMy)/d;
	pRotor->GetConnectedNodes(connectedNodes);

	dThrust, dRollMoment, dPitchMoment // from drives?
	dRho // from AirProperties
	dOmega // from kinematics
	dRadius // constant
	dMu // from drive?
 */

class RotorTrimBase : virtual public Elem, public Genel {
protected:
	doublereal dRadius;
	const ScalarDifferentialNode* pvNodes[3];
	DriveOwner DThrust, DRollMoment, DPitchMoment;
	DriveOwner Trigger;
	
	doublereal dGamma;
	doublereal dP;
	
	doublereal dP2;
	doublereal dC;
	doublereal dC2;
	
	doublereal dTau0;
	doublereal dTau1;
	doublereal dKappa0;
	doublereal dKappa1;

protected:
	virtual void
	GetData(doublereal &dThrust,
		doublereal &dRollMoment,
		doublereal &dPitchMoment,
		doublereal &dRho,
		doublereal &dOmega,
		doublereal &dMu) const = 0;

public:
	RotorTrimBase(unsigned int uL,
		const DofOwner* pDO,
		const ScalarDifferentialNode* pNode1,
		const ScalarDifferentialNode* pNode2,
		const ScalarDifferentialNode* pNode3,
		const DriveCaller* pDThrust,
		const DriveCaller* pDRollMoment,
		const DriveCaller* pDPitchMoment,
		const doublereal& dG,
		const doublereal& dp,
		const doublereal& dT0,
		const doublereal& dT1,
		const doublereal& dK0,
		const doublereal& dK1,
		const DriveCaller *pTrigger,
		flag fOut);

	virtual ~RotorTrimBase(void);
	
	virtual unsigned int iGetNumDof(void) const;
	
	/* Tipo di Genel */
	virtual Genel::Type GetGenelType(void) const { 
		return Genel::ROTORTRIM; 
	};
	
	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal dCoef,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ );
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ );

 	/* *******PER IL SOLUTORE PARALLELO******** */        
   	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

class RotorTrim : virtual public Elem, public RotorTrimBase {
protected:
	const Rotor* pRotor;

	virtual void
	GetData(doublereal &dThrust,
		doublereal &dRollMoment,
		doublereal &dPitchMoment,
		doublereal &dRho,
		doublereal &dOmega,
		doublereal &dMu) const;

public:
	RotorTrim(unsigned int uL,
		const DofOwner* pDO,
		const Rotor* pRot, 
		const ScalarDifferentialNode* pNode1,
		const ScalarDifferentialNode* pNode2,
		const ScalarDifferentialNode* pNode3,
		const DriveCaller* pDThrust,
		const DriveCaller* pDRollMoment,
		const DriveCaller* pDPitchMoment,
		const doublereal& dG,
		const doublereal& dp,
		const doublereal& dT0,
		const doublereal& dT1,
		const doublereal& dK0,
		const doublereal& dK1,
		const DriveCaller *pTrigger,
		flag fOut);

	virtual ~RotorTrim(void);
	
	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

 	/* *******PER IL SOLUTORE PARALLELO******** */        
   	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

class RotorTrimGeneric : virtual public Elem, public RotorTrimBase {
protected:
	const StructNode *pStrNode;
	DriveOwner Thrust, RollMoment, PitchMoment;
	const AirProperties *pAP;
	DriveOwner Omega, Mu;

	virtual void
	GetData(doublereal &dThrust,
		doublereal &dRollMoment,
		doublereal &dPitchMoment,
		doublereal &dRho,
		doublereal &dOmega,
		doublereal &dMu) const;

public:
	RotorTrimGeneric(unsigned int uL,
		const DofOwner* pDO,
		const StructNode *pStrNode,
		const DriveCaller *pThrust,
		const DriveCaller *pRollMoment,
		const DriveCaller *pPitchMoment,
		const AirProperties *pAP,
		doublereal dRadius,
		const DriveCaller *pOmega,
		const DriveCaller *pMu,
		const ScalarDifferentialNode* pNode1,
		const ScalarDifferentialNode* pNode2,
		const ScalarDifferentialNode* pNode3,
		const DriveCaller* pDThrust,
		const DriveCaller* pDRollMoment,
		const DriveCaller* pDPitchMoment,
		const doublereal& dG,
		const doublereal& dp,
		const doublereal& dT0,
		const doublereal& dT1,
		const doublereal& dK0,
		const doublereal& dK1,
		const DriveCaller *pTrigger,
		flag fOut);

	virtual ~RotorTrimGeneric(void);
	
	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
};

#endif // ROTTRIM_H

