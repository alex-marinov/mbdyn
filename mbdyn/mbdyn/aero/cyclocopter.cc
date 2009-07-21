/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* Elementi di rotore */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <limits>
#include <cmath>

#include "cyclocopter.h"
#include "dataman.h"

CyclocopterKARI::CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, flag fOut)
: Elem(uL, fOut),
InducedVelocity(uL, pDO, pC, ppres, fOut),
pRotor(pR),
RRot(rrot)
{
	NO_OP;
}

CyclocopterKARI::~CyclocopterKARI(void)
{
	NO_OP;
}

InducedVelocity::Type
CyclocopterKARI::GetInducedVelocityType(void) const
{
	return InducedVelocity::CYCLOCOPTER;
}

void
CyclocopterKARI::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterKARI::Output(OutputHandler& OH) const
{
	NO_OP;
}

std::ostream&
CyclocopterKARI::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

void
CyclocopterKARI::SetInitialValue(VectorHandler& X)
{
	NO_OP;
}

SubVectorHandler&
CyclocopterKARI::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	NO_OP;
}

#if 0
void
CyclocopterKARI::AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X)
{
	InducedVelocity::AddForce(uL, F, M, X);
}
#endif

Vec3
CyclocopterKARI::GetInducedVelocity(const Vec3& X) const
{
	return Zero3;
}

void
CyclocopterKARI::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

Elem*
ReadCyclocopter(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO, 
	unsigned int uLabel,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor)
{
	Elem *pEl = 0;

	if (HP.IsKeyWord("kari")) {
     		ResForceSet **ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);

		pEl = new CyclocopterKARI(uLabel, pDO,
			pCraft, rrot, pRotor,
  			ppres, fOut);

	} else {
		silent_cerr("Rotor(" << uLabel << "): "
			"unknown cyclocopter inflow model "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
}

