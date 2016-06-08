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

/* BufferStreamElem_base - begin */

BufferStreamElem_base::BufferStreamElem_base(unsigned int uL,
	unsigned int oe,
	StreamContent *pSC, StreamOutEcho *pSOE)
: Elem(uL, flag(0)),
StreamOutElem(uL, "buffer", oe),
pSC(pSC), pSOE(pSOE)
{
	if (pSOE != 0) {
		pSOE->Init("BufferStreamElem_base", uLabel, pSC->GetNumChannels());
	}
}

BufferStreamElem_base::~BufferStreamElem_base(void)
{
	if (pSC != 0) {
		SAFEDELETE(pSC);
	}

	if (pSOE != 0) {
		delete pSOE;
	}
}

const integer
BufferStreamElem_base::GetBufSize(void) const
{
	return pSC->GetNumChannels();
}

void
BufferStreamElem_base::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	// do not send "derivatives"
	OutputCounter = -1;
}

void
BufferStreamElem_base::AfterConvergence(const VectorHandler& X, 
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
	if (pSOE != 0) {
		pSOE->Echo(pd, pSC->GetNumChannels());
	}

	doublereal *pBuffer = const_cast<doublereal *>(GetBufRaw());
	for (unsigned i = 0; i < pSC->GetNumChannels(); i++) {
		pBuffer[i] = pd[i];
	}
}

void
BufferStreamElem_base::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

/* BufferStreamElem_base - end */

/* BufferStreamElem - begin */

BufferStreamElem::BufferStreamElem(unsigned int uL,
	unsigned int oe,
	StreamContent *pSC, StreamOutEcho *pSOE)
: Elem(uL, flag(0)),
BufferStreamElem_base(uL, oe, pSC, pSOE),
buffer(pSC->GetNumChannels())
{
	NO_OP;
}

BufferStreamElem::~BufferStreamElem(void)
{
	NO_OP;
}

const doublereal *
BufferStreamElem::GetBufRaw(void) const
{
	// paranoid sanity check: callers of GetBuf() could have altered the size of the buffer...
	ASSERT(buffer.size() == pSC->GetNumChannels());

	return &buffer[0];
}

const std::vector<doublereal>&
BufferStreamElem::GetBuf(void) const
{
	// paranoid sanity check: callers of GetBuf() could have altered the size of the buffer...
	ASSERT(buffer.size() == pSC->GetNumChannels());

	return buffer;
}

std::ostream&
BufferStreamElem::Restart(std::ostream& out) const
{   	
	return out << "# BufferStreamElem(" << GetLabel() << "): "
		"not implemented yet" << std::endl;
}	

/* BufferStreamElem - end */

/* BufferStreamElemRaw - begin */

BufferStreamElemRaw::BufferStreamElemRaw(unsigned int uL,
	unsigned int oe,
	StreamContent *pSC, StreamOutEcho *pSOE,
	bool bOwnsMemory)
: Elem(uL, flag(0)),
BufferStreamElem_base(uL, oe, pSC, pSOE),
m_bOwnsMemory(bOwnsMemory),
m_pBuffer(0)
{
	if (m_bOwnsMemory) {
		m_pBuffer = new doublereal[pSC->GetNumChannels()];
	}
}

BufferStreamElemRaw::~BufferStreamElemRaw(void)
{
	if (m_bOwnsMemory) {
		delete[] m_pBuffer;
	}
}

bool
BufferStreamElemRaw::bOwnsMemory(void) const
{
	return m_bOwnsMemory;
}

void
BufferStreamElemRaw::SetBufRaw(integer n, const doublereal *p)
{
	if (n != integer(pSC->GetNumChannels())) {
		// error
		std::ostringstream os;
		os << "setting buffer pointer in BufferStreamElemRaw of wrong size (original=" << pSC->GetNumChannels() << ", new=" << n << ")";
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, os.str());
	}

	if (m_bOwnsMemory) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "setting buffer pointer in BufferStreamElemRaw that owns its memory");
	}

	if (m_pBuffer != 0) {
		// error; maybe we could simply replace it, couldn't we?
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "setting buffer pointer in BufferStreamElemRaw that has already been set");
	}

	m_pBuffer = p;
}

const doublereal *
BufferStreamElemRaw::GetBufRaw(void) const
{
	return m_pBuffer;
}

std::ostream&
BufferStreamElemRaw::Restart(std::ostream& out) const
{   	
	return out << "# BufferStreamElem(" << GetLabel() << "): "
		"not implemented yet" << std::endl;
}	

/* BufferStreamElemRaw - end */

Elem *
ReadBufferStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type)
{
	enum {
		STL,
		RAW
	} eType = STL;
	bool bOwnsMemory(true);

	if (HP.IsKeyWord("type")) {
		if (HP.IsKeyWord("raw")) {
			eType = RAW;
			if (HP.IsKeyWord("owns" "memory")) {
				bOwnsMemory = HP.GetYesNoOrBool();
			}

		} else if (!HP.IsKeyWord("stl")) {
			silent_cerr("BufferStreamDrive"
				"(" << uLabel << "\"): "
				"invalid type at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

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
	switch (eType) {
	case STL:
		SAFENEWWITHCONSTRUCTOR(pEl, BufferStreamElem,
			BufferStreamElem(uLabel, OutputEvery, pSC, pSOE));
		break;

	case RAW:
		SAFENEWWITHCONSTRUCTOR(pEl, BufferStreamElemRaw,
			BufferStreamElemRaw(uLabel, OutputEvery, pSC, pSOE, bOwnsMemory));
		break;

	default:
		ASSERT(0);
		break;
	}

	return pEl;
}

