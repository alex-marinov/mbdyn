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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <motor.h>

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

Motor::Motor(unsigned int uL, const DofOwner* pD, 
		const StructNode* pN1, const StructNode* pN2,
		const AbstractNode* pV1, const AbstractNode* pV2,
		const Vec3& TmpDir, doublereal dG,
		doublereal dl, doublereal dr,
		flag fOut)
: Elem(uL, Elem::ELECTRIC, fOut), 
Electric(uL, Electric::MOTOR, pDO, fOut),
pStrNode1(pN1), pStrNode2(pN2), pVoltage1(pV1), pVoltage2(pV2),
Dir(TmpDir), dGain(dG), dL(dl), dR(dr)
{
	NO_OP;
}

Motor::~Motor(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
Motor::Restart(std::ostream& out) const
{
	return out;
}
   
unsigned int
iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Motor::SetDof(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 9;
	*piNumCols = 9;
}
      
VariableSubMatrixHandler&
AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	return WorkMat;
}

SubVectorHandler&
AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	return WorkVec;
}

void
SetInitialValue(VectorHandler& /* X */ ) const
{
	NO_OP;
}

void
SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const
{
	NO_OP;
}

