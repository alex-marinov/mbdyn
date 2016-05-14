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
 * Author: Michele Attolico <attolico@aero.polimi.it>
 */

#ifndef BUFFERSTREAM_OUT_ELEM_H
#define BUFFERSTREAM_OUT_ELEM_H

#include "streamoutelem.h"

/* BufferStreamElem_base - begin */

class BufferStreamElem_base : public StreamOutElem, virtual public Elem {
protected:
	StreamContent *pSC;
	StreamOutEcho *pSOE;
	
public:
   	BufferStreamElem_base(unsigned int uL, unsigned int oe,
		StreamContent *pSC, StreamOutEcho *pSOE);

   	virtual ~BufferStreamElem_base(void);

	const integer GetBufSize(void) const;
	virtual const doublereal * GetBufRaw(void) const = 0;

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);

	/* Inverse Dynamics */
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP);
};

/* BufferStreamElem_base - end */

/* BufferStreamElem - begin */

class BufferStreamElem : public BufferStreamElem_base, virtual public Elem {
protected:
	std::vector<doublereal> buffer;
	
public:
   	BufferStreamElem(unsigned int uL, unsigned int oe,
		StreamContent *pSC, StreamOutEcho *pSOE);

   	virtual ~BufferStreamElem(void);

	virtual const doublereal * GetBufRaw(void) const;
	const std::vector<doublereal>& GetBuf(void) const;

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* BufferStreamElem - end */

/* BufferStreamElemRaw - begin */

class BufferStreamElemRaw : public BufferStreamElem_base, virtual public Elem {
protected:
	bool bOwnsMemory;
	const doublereal *pBuffer;
	
public:
   	BufferStreamElemRaw(unsigned int uL, unsigned int oe,
		StreamContent *pSC, StreamOutEcho *pSOE,
		bool bOwnsMemory);

   	virtual ~BufferStreamElemRaw(void);

	virtual void SetBufRaw(integer n, const doublereal *p);
	virtual const doublereal * GetBufRaw(void) const;

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* BufferStreamElem - end */

class DataManager;
class MBDynParser;

extern Elem *
ReadBufferStreamElem(DataManager *pDM, MBDynParser& HP,
	unsigned int uLabel, StreamContent::Type type);

/* BufferStreamElem - end */

#endif /* BUFFERSTREAM_OUT_ELEM_H */

