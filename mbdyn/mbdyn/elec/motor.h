/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

	M = M0(Phi_e) - (Gain + M1(Phi_e)) * i
	
	C1 = -M

	C2 = M

	i1 = - i

	i2 = i
	
	    d i
	L * --- + R * i = - Gain * (Omega2 - Omega1) + V2 - V1
	    d t
 
 * In order to take into account the so called cogging torque or ripple torque 
 * of DC motors, the following data can be provided:
 *   M0(Phi_e) ... variation of the motor torque independent of current
 *   M1(Phi_e) ... variation of the motor torque proportional to current  
 *   p ... number of terminal pairs
 *
 *   Phi_m ... mechanical angle of rotation between rotor and stator
 *   Phi_e = p * Phi_m ... electric angle between rotor field and stator field
 */

class Motor : virtual public Elem, public Electric {
private:
	const StructNode *pStrNode1;
	const StructNode *pStrNode2;
	const ElectricNode *pVoltage1;
	const ElectricNode *pVoltage2;

	Mat3x3 Rn;
	doublereal dGain;
	doublereal dL;
	DriveOwner dR;
	integer p;
	DriveOwner M0, M1;

	inline doublereal dGetOmega() const;
	inline doublereal dGetVoltage() const;
	inline doublereal dGetOmega(const Vec3& n) const;
	inline Vec3 GetAxisOfRotation() const;
	inline doublereal dGetPhiMechanical() const;
	inline doublereal dGetPhiElectric(doublereal Phi_m) const;

	doublereal M, i, iP, Phi_m, Phi_e;

	static const int iNumPrivData = 12;

	static const struct PrivData {
	int index;
		char name[8];
	} rgPrivData[iNumPrivData];

public:
	Motor(unsigned int uL, const DofOwner* pD, 
			const StructNode* pN1, const StructNode* pN2,
			const ElectricNode* pV1, const ElectricNode* pV2,
			const Mat3x3& Rn, doublereal dG,
			doublereal dl, DriveCaller* dR, doublereal i0,
			integer p, const DriveCaller* pM0, const DriveCaller* pM1,
			flag fOut);
	virtual ~Motor(void);

	virtual Electric::Type GetElectricType(void) const;

	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

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

	virtual void SetInitialValue(VectorHandler& X);

	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints* h = 0);
        
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
};

#endif /* MOTOR_H */

