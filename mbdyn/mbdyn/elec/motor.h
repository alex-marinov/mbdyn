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

#ifndef MOTOR_H
#define MOTOR_H

/*
 * Electric motor: an internal couple between two structural nodes
 * whose value is represented by an internal state that is a current,
 * according to equations:

	C1 = Gain * i
	
	C2 = - Gain * i

	i1 = i

	i2 = - i
	
	    d i
	L * --- + R * i = - Gain * Omega + V2 - V1
	    d t
 
 */

class Motor : virtual public Elem, public Electric {
private:
	const StructNode *pStrNode1;
	const StructNode *pStrNode2;
	const ElectricNode *pVoltage1;
	const ElectricNode *pVoltage2;

	Vec3 Dir;
	doublereal dGain;
	doublereal dL;
	doublereal dR;

public:
	Motor(unsigned int uL, const DofOwner* pD, 
			const StructNode* pN1, const StructNode* pN2,
			const ElectricNode* pV1, const ElectricNode* pV2,
			const Vec3& TmpDir, doublereal dG,
			doublereal dl, doublereal dr,
			flag fOut);
	virtual ~Motor(void);

	virtual Electric::Type GetElectricType(void) const {
		return Electric::MOTOR;
	};
   
	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
   
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
     	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
   
     	virtual void SetInitialValue(VectorHandler& /* X */ );

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
     	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(4);
		connectedNodes[0] = pStrNode1;
		connectedNodes[1] = pStrNode2;
		connectedNodes[2] = pVoltage1;
		connectedNodes[3] = pVoltage2;
	};
	/* ************************************************ */
};

#endif /* MOTOR_H */

