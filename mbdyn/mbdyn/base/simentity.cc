/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <sstream>

#include "simentity.h"
#include "dataman.h"
#include "tpldrive.h"
#include "hint_impl.h"

/* SimulationEntity - begin */

SimulationEntity::SimulationEntity(void)
{
	NO_OP;
}

SimulationEntity::~SimulationEntity(void)
{
	NO_OP;
}

bool
SimulationEntity::bIsValidIndex(unsigned int i) const
{
	if (i >= 1 && i <= iGetNumDof()) {
		return true;
	}
	return false;
}

void 
SimulationEntity::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
#if 0
	if (*ppX0_Xp0 != NULL) {
		// scrive in X e XP i dati iniziali
	}
#endif
	NO_OP;
}
     
Hint *
SimulationEntity::ParseHint(DataManager *pDM, const char *s) const
{
	return ::ParseHint(pDM, s);
}

void 
SimulationEntity::BeforePredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ ,
		VectorHandler& /* XPrev */ ,
		VectorHandler& /* XPPrev */ ) const
{
	NO_OP;
}
	
void 
SimulationEntity::AfterPredict(VectorHandler& /* X */ , 
		VectorHandler& /* XP */ )
{
	NO_OP;
}

void
SimulationEntity::Update(const VectorHandler& /* XCurr */ , 
		const VectorHandler& /* XPrimeCurr */ )
{
	NO_OP;
}

void
SimulationEntity::DerivativesUpdate(const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	Update(XCurr, XPrimeCurr);
}

/* Inverse Dynamics: */
void
SimulationEntity::Update(const VectorHandler& /* XCurr */ , 
		InverseDynamics::Order /* iOrder */ )
{
	NO_OP;
}

void 
SimulationEntity::AfterConvergence(const VectorHandler& /* X */ , 
		const VectorHandler& /* XP */ )
{
	NO_OP;
}

/* Inverse Dynamics:*/
void 
SimulationEntity::AfterConvergence(const VectorHandler& /* X */,
			const VectorHandler&  /* XP */,
			const VectorHandler&  /* XPP */)
{
	NO_OP;
}

unsigned int
SimulationEntity::iGetNumPrivData(void) const 
{
	return 0;
}

unsigned int
SimulationEntity::iGetPrivDataIdx(const char *s) const 
{
	silent_cerr("no private data available" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
SimulationEntity::dGetPrivData(unsigned int /* i */ ) const
{
	silent_cerr("no private data available" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

std::ostream&
SimulationEntity::OutputAppend(std::ostream& out) const
{
	return out;
}

void SimulationEntity::ReadInitialState(MBDynParser& HP)
{
	NO_OP;
}

/* SimulationEntity - end */

