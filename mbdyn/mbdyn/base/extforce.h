/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2008
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

/* Forza */

#ifndef EXTFORCE_H
#define EXTFORCE_H

#include <vector>
#include <string>
#include "force.h"

/* ExtFileHandlerBase - begin */

class ExtFileHandlerBase {
public:
	virtual ~ExtFileHandlerBase(void);

	virtual std::ostream& Send_pre(bool bAfterConvergence = false) = 0;
	virtual void Send_post(bool bAfterConvergence = false) = 0;

	virtual std::istream& Recv_pre(void) = 0;
	virtual void Recv_post(void) = 0;
};

/* ExtFileHandlerBase - end */

/* ExtFileHandler - begin */

class ExtFileHandler : public ExtFileHandlerBase {
protected:
	std::string fin, fout, tmpout;
	bool bRemoveIn, bNoClobberOut;
	int iSleepTime, iPrecision;

	std::ifstream inf;
	std::ofstream outf;

public:
	ExtFileHandler(std::string& fin,
		bool bRemoveIn,
	        std::string& fout,
		bool bNoClobberOut,
		int iSleepTime,
		int iPrecision);
	~ExtFileHandler(void);

	virtual std::ostream& Send_pre(bool bAfterConvergence = false);
	virtual void Send_post(bool bAfterConvergence = false);

	virtual std::istream& Recv_pre(void);
	virtual void Recv_post(void);
};

/* ExtFileHandler - end */

/* ExtFileHandlerEDGE moved to extedge.h, extedge.cc */

/* ExtForce - begin */

class ExtForce : virtual public Elem, public Force {
protected:
	ExtFileHandlerBase *pEFH;

	bool bFirstRes;
	int iCoupling, iCouplingCounter;

	void Send(bool bAfterConvergence = false);
	void Recv(void);

	virtual void Send(std::ostream& out, bool bAfterConvergence = false) = 0;
	virtual void Recv(std::istream& in) = 0;
   
public:
	/* Costruttore */
	ExtForce(unsigned int uL,
		ExtFileHandlerBase *pEFH,
		int iCoupling,
		flag fOut);

	virtual ~ExtForce(void);

	virtual void Update(const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* ExtForce - end */

class DataManager;
class MBDynParser;

extern void
ReadExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel,
	ExtFileHandlerBase*& pEFH,
	int& iCoupling);

#endif // EXTFORCE_H

