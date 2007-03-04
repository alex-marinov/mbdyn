/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifndef SOCKETSTREAMMOTIONELEM_H
#define SOCKETSTREAMMOTIONELEM_H

/*include del programma*/

#include <elem.h>
#include <strnode.h>

/* SocketStreamMotionElem - begin */

class SocketStreamMotionElem : virtual public Elem {
protected:
	std::vector<StructNode *> nodes;

	/* Stream buffer */
	int size;
	char *buf;

	/* output with lower ratio */
	unsigned int OutputEvery;
	mutable unsigned int OutputCounter;

	UseSocket *pUS;

	const char *name;
	int send_flags;
	bool bSendFirst;
	
public:
   	SocketStreamMotionElem(unsigned int uL, std::vector<StructNode *>& n,
		unsigned int oe,
		DataManager *pDM,
		const char *h, const char *m, unsigned short int p,
		bool c, int flags, bool bSendFirst);
   	SocketStreamMotionElem(unsigned int uL, std::vector<StructNode *>& n,
		unsigned int oe,
		DataManager *pDM,
		const char *m, const char* const Path,
		bool c, int flags, bool bSendFirst);
			
   	virtual ~SocketStreamMotionElem(void);

	virtual std::ostream& Restart(std::ostream& out) const;
	virtual Elem::Type GetElemType(void) const;
	virtual void WorkSpaceDim(integer* piRows, integer* piCols) const;
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WrokVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP);
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP);

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);
};

class DataManager;
class MBDynParser;

extern Elem *
ReadSocketStreamMotionElem(DataManager *pDM,
	MBDynParser& HP,
	unsigned int uLabel);

/* SocketSteamElem - end */

#endif /* SOCKETSTREAMMOTIONELEM_H */

