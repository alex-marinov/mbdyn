/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/extforce.h,v 1.30 2017/01/12 14:46:09 masarati Exp $ */
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

#ifndef EXTSOCKET_H
#define EXTSOCKET_H

#include <mbconfig.h>

#ifdef USE_SOCKET

#include <vector>
#include <string>

#include "extforce.h"
#include "mbsleep.h"
#include "force.h"
#include "converged.h"
#include "usesock.h"
#include "mbc.h"


class ExtSocketHandler : public ExtRemoteHandler {
protected:
	UseSocket *pUS;
	unsigned recv_flags;
	unsigned send_flags;

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

	virtual ExtFileHandlerBase::Type GetType(void) { return ExtFileHandlerBase::TYPE_SOCKET; };

	virtual int GetOutFileDes(void);
	virtual int GetSendFlags(void) const;
	virtual int GetInFileDes(void);
	virtual int GetRecvFlags(void) const;
};
#endif // USE_SOCKET

class DataManager;
class MBDynParser;

extern ExtFileHandlerBase *
ReadExtSocketHandler(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel);

#endif // EXTSOCKET_H

