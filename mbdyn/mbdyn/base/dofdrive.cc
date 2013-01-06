/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Drive che usa  i gradi di liberta' */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dofdrive.h"


DofDriveCaller::DofDriveCaller(const DriveHandler* pDH, 
		const DriveCaller* pDC,
		const ScalarDof& sd)
: DriveCaller(pDH), DriveOwner(pDC), SD(sd)
{
	NO_OP;
};     	

DofDriveCaller::~DofDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller*
DofDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;

	SAFENEWWITHCONSTRUCTOR(pDC, 
			DofDriveCaller,
			DofDriveCaller(pDrvHdl,
				pGetDriveCaller()->pCopy(), SD));
   
	return pDC;
}


/* Restart */
std::ostream&
DofDriveCaller::Restart(std::ostream& out) const
{
	out << " dof, " 
		<< SD.pNode->GetLabel() << ", "
		<< psReadNodesNodes[SD.pNode->GetNodeType()];
	if (SD.iOrder == 0) {
		out << ", algebraic, ";
	} else if (SD.iOrder == 1) {	
		out << ", differential, ";
	} else {
		out << ", order, " << SD.iOrder << ", ";
	}
	return DriveOwner::pGetDriveCaller()->Restart(out);
}

