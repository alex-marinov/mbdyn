/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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

#include "mbsleep.h"
#include "force.h"
#include "converged.h"

#ifdef USE_SOCKET
#include "usesock.h"
#endif // USE_SOCKET
#include "mbc.h"

/* ExtFileHandlerBase - begin */

class ExtFileHandlerBase {
public:
	enum SendWhen {
		SEND_FIRST_TIME,
		SEND_REGULAR,
		SEND_AFTER_CONVERGENCE
	};

	enum Negotiate {
		NEGOTIATE_NO,
		NEGOTIATE_CLIENT,
		NEGOTIATE_SERVER
	};

	enum OutputFlags {
		OUTPUT_KIN = 0x1000,
		OUTPUT_DYN = 0x2000,
		OUTPUT_REF = 0x4000,

		OUTPUT_REF_KIN = (OUTPUT_REF | OUTPUT_KIN),
		OUTPUT_REF_DYN = (OUTPUT_REF | OUTPUT_DYN)
	};

protected:
	std::streamsize Precision;
	mbsleep_t SleepTime;
	bool bOK;

public:
	ExtFileHandlerBase(mbsleep_t SleepTime, std::streamsize Precision);
	virtual ~ExtFileHandlerBase(void);

	// NOTE: returns true if Prepare() must be called
	virtual bool Prepare_pre(void) = 0;
	// NOTE: returns true if peer must request negotiation during Prepare
	virtual Negotiate NegotiateRequest(void) const;
	virtual void Prepare_post(bool ok) = 0;

	virtual void AfterPredict(void) = 0;

	// NOTE: returns true if Send() must be called
	virtual bool Send_pre(SendWhen when) = 0;
	virtual void Send_post(SendWhen when) = 0;

	// NOTE: returns true if Recv() must be called
	virtual bool Recv_pre(void) = 0;
	// NOTE: returns true if converged
	virtual bool Recv_post(void) = 0;

	virtual std::ostream *GetOutStream(void);
	virtual std::istream *GetInStream(void);
	virtual int GetOutFileDes(void);
	virtual int GetSendFlags(void) const;
	virtual int GetInFileDes(void);
	virtual int GetRecvFlags(void) const;
};

/* ExtFileHandlerBase - end */

/* ExtFileHandler - begin */

class ExtFileHandler : public ExtFileHandlerBase {
protected:
	std::string fin, fout, tmpout;
	bool bRemoveIn, bNoClobberOut;

	std::ifstream inf;
	std::ofstream outf;

public:
	ExtFileHandler(std::string& fin,
		bool bRemoveIn,
	        std::string& fout,
		bool bNoClobberOut,
		mbsleep_t SleepTime,
		std::streamsize Precision);
	virtual ~ExtFileHandler(void);

	virtual bool Prepare_pre(void);
	virtual void Prepare_post(bool ok);

	virtual void AfterPredict(void);

	virtual bool Send_pre(SendWhen when);
	virtual void Send_post(SendWhen when);

	virtual bool Recv_pre(void);
	virtual bool Recv_post(void);

	virtual std::ostream *GetOutStream(void);
	virtual std::istream *GetInStream(void);
};

/* ExtFileHandler - end */

/* ExtSocketHandler - begin */

#ifdef USE_SOCKET
class ExtSocketHandler : public ExtFileHandlerBase {
protected:
	UseSocket *pUS;
	unsigned recv_flags;
	unsigned send_flags;
	bool bReadForces;
	bool bLastReadForce;

	ESCmd u2cmd(unsigned u) const;
	const char *cmd2str(ESCmd cmd) const;

public:
	ExtSocketHandler(UseSocket *pUS, mbsleep_t SleepTime,
		int recv_flags, int send_flags);
	virtual ~ExtSocketHandler(void);

	virtual bool Prepare_pre(void);
	virtual Negotiate NegotiateRequest(void) const;
	virtual void Prepare_post(bool ok);

	virtual void AfterPredict(void);

	virtual bool Send_pre(SendWhen when);
	virtual void Send_post(SendWhen when);

	virtual bool Recv_pre(void);
	virtual bool Recv_post(void);

	virtual int GetOutFileDes(void);
	virtual int GetSendFlags(void) const;
	virtual int GetInFileDes(void);
	virtual int GetRecvFlags(void) const;
};
#endif // USE_SOCKET

/* ExtSocketHandler - end */

/* ExtFileHandlerEDGE moved to extedge.h, extedge.cc */

/* ExtForce - begin */

class ExtForce : virtual public Elem, public Force {
protected:
	Converged c;
	ExtFileHandlerBase *pEFH;

	// exchange after predict?
	bool bSendAfterPredict;

public:
	//-2: no coupling, only send kinematics
	//-1: staggered, send kinematics at end of step, use forces at next step
	// 0: loose coupling
	// 1: tight coupling
	// >1: exchange every iCoupling iterations
	enum {
		COUPLING_NONE = -2,
		COUPLING_STAGGERED = -1,
		COUPLING_LOOSE = 0,
		COUPLING_TIGHT = 1
	};

protected:
	int iCoupling;

	// iteration counter
	mutable int iCouplingCounter;

	// whether the current residual is the first or not...
	mutable bool bFirstSend;
	mutable bool bFirstRecv;

	void Send(ExtFileHandlerBase::SendWhen when);
	void Recv(void);

	virtual bool Prepare(ExtFileHandlerBase *pEFH) = 0;
	virtual void Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when) = 0;
	virtual void Recv(ExtFileHandlerBase *pEFH) = 0;
   
public:
	/* Costruttore */
	ExtForce(unsigned int uL,
		DataManager *pDM,
		ExtFileHandlerBase *pEFH,
		bool bSendAfterPredict,
		int iCoupling,
		flag fOut);

	virtual ~ExtForce(void);

	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints* h = 0);
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
	bool& bSendAfterPredict,
	int& iCoupling);

extern void
ReadExtFileParams(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel,
	mbsleep_t& SleepTime,
	std::streamsize& Precision);

#endif // EXTFORCE_H

