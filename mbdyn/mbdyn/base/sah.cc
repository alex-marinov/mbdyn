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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "sah.h"

/* SampleAndHold - begin */

SampleAndHold::SampleAndHold(unsigned int uL,
		const DofOwner* pDO,
		const ScalarNode *pN,
		DriveCaller *pDC,
		const doublereal &dSP,
		flag fOut)
: ParameterNode(uL, pDO, 0., fOut),
pNode(pN),
Time(pDC),
dSamplePeriod(dSP)
{
	ASSERT(pN != 0);
	ASSERT(pDC != 0);
	ASSERT(dSP > 0.);

	dSampleTime = Time.dGet() + dSamplePeriod;
}

SampleAndHold::~SampleAndHold(void)
{
	NO_OP;
}

/* Contributo del nodo al file di restart */
std::ostream&
SampleAndHold::Restart(std::ostream& out) const
{
	return out << "# sample and hold parameter to be implemented" 
		<< std::endl;
}

void
SampleAndHold::SetValue(DataManager *pDM,
		VectorHandler& /* X */, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	dSampleTime = Time.dGet() + dSamplePeriod;
	dX = pNode->dGetX();
}

/* Aggiorna i valori interni */   
void
SampleAndHold::Update(const VectorHandler& /* X */,
		const VectorHandler& /* XP */ )
{
	doublereal dT = Time.dGet();
	if (dT >= dSampleTime) {
		dX = pNode->dGetX();
		dSampleTime += dSamplePeriod;
	}
}

/* SampleAndHold - end */

