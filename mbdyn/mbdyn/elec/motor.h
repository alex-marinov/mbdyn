/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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

	virtual inline void* pGet(void) const;

	virtual Electric::Type GetElectric(void) const {
		return Electric::MOTOR;
	};
   
	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
   
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order SetDof(unsigned int i) const;

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
     	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
   
     	virtual void SetInitialValue(VectorHandler& /* X */ ) const;
     	virtual void SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
     	virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps,
			unsigned int* NdLabels) {
		NumNodes = 4;
		NdTyps[0] = pStrNode1->GetNodeType();
		NdLabels[0] = pStrNode1->GetLabel();
		NdTyps[1] = pStrNode2->GetNodeType();
		NdLabels[1] = pStrNode2->GetLabel();
		NdTyps[2] = pVoltage1->GetNodeType();
		NdLabels[2] = pVoltage1->GetLabel();
		NdTyps[3] = pVoltage2->GetNodeType();
		NdLabels[3] = pVoltage2->GetLabel();
	};
	/* ************************************************ */
};

inline void*
Motor::pGet(void) const
{
	return (void *)this;
}

#endif /* MOTOR_H */

