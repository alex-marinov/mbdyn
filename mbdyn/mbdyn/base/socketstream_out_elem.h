/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef SOCKETSTREAM_OUT_ELEM_H
#define SOCKETSTREAM_OUT_ELEM_H

#include "streamoutelem.h"

#ifdef USE_SOCKET

/* SocketStreamElem - begin */

class SocketStreamElem : public StreamOutElem, virtual public Elem {
protected:
	UseSocket *pUS;
	StreamContent *pSC;

	int send_flags;
	bool bSendFirst;
	bool bAbortIfBroken;

	std::string sOutFileName;
	std::ofstream outFile;
	int iPrecision;
	doublereal dShift;
	
public:
   	SocketStreamElem(unsigned int uL, const std::string& name,
		unsigned int oe,
		UseSocket *pUS, StreamContent *pSC,
		int flags, bool bSendFirst, bool bAbortIfBroken,
		const std::string& sOutFileName, int iPrecision,
		doublereal dShift);

   	virtual ~SocketStreamElem(void);

	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);

	/* Inverse Dynamics */
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP);
};

#endif // USE_SOCKET

class DataManager;
class MBDynParser;

extern Elem *
ReadSocketStreamElem(DataManager *pDM, MBDynParser& HP,
	unsigned int uLabel, StreamContent::Type type);

/* SocketStreamElem - end */

#endif /* SOCKETSTREAM_OUT_ELEM_H */

