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

/*
 * Michele Attolico <attolico@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "bufferstream_out_elem.h"

/* BufferStreamElem - begin */

BufferStreamElem::BufferStreamElem(unsigned int uL,
	unsigned int oe,
	StreamContent *pSC, StreamOutEcho *pSOE)
: Elem(uL, flag(0)),
StreamOutElem(uL, 0, oe),
pSC(pSC), pSOE(pSOE)
{
	pSOE->Init("BufferStreamElem", uLabel, pSC->GetNumChannels());
}

BufferStreamElem::~BufferStreamElem(void)
{
	if (pSC != 0) {
		SAFEDELETE(pSC);
	}

	if (pSOE != 0) {
		delete pSOE;
	}
}

const std::vector<doublereal>&
BufferStreamElem::GetBuf(void) const
{
	return buffer;
}

std::ostream&
BufferStreamElem::Restart(std::ostream& out) const
{   	
	return out << "# BufferStreamElem(" << GetLabel() << "): "
		"not implemented yet" << std::endl;
}	

void
BufferStreamElem::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	// do not send "derivatives"
	OutputCounter = -1;
}

void
BufferStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	/* output only every OutputEvery steps */
	OutputCounter++;
	if (OutputCounter != OutputEvery) {
		return;
	}
	OutputCounter = 0;

	// prepare the output buffer
	pSC->Prepare();

	// check whether echo is needed
	const doublereal *pd = (const doublereal *)pSC->GetBuf();
	pSOE->Echo(pd, pSC->GetNumChannels());

	for (unsigned i = 0; i < pSC->GetNumChannels(); i++) {
		buffer[i] = pd[i];
	}
}

void
BufferStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}


Elem *
ReadBufferStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type)
{
	unsigned int OutputEvery = 1;
	if (HP.IsKeyWord("output" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("BufferStreamElem(" << uLabel << "): "
				"invalid output every value " << i << " "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		OutputEvery = (unsigned int)i;
	}

	StreamOutEcho *pSOE = ReadStreamOutEcho(HP);
	StreamContent *pSC = ReadStreamContent(pDM, HP, type);

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("BufferStreamElem(" << uLabel << "): "
			"semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, BufferStreamElem,
		BufferStreamElem(uLabel, OutputEvery, pSC, pSOE));

	return pEl;
}

