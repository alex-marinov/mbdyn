/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2009
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

#ifndef EXTEDGE_H
#define EXTEDGE_H

#include "extforce.h"

/* ExtFileHandlerEDGE - begin */

class ExtFileHandlerEDGE : public ExtFileHandlerBase {
protected:
	std::string fflagname, fdataname;
	bool bReadForces;

	std::ifstream infile;
	std::ofstream outfile;

	enum EDGEcmd {
		// ?
		EDGE_UNKNOWN		= -1,

		// EDGE is initializing; MBDyn waits
		EDGE_INITIALIZING	= 0,

		// EDGE is busy; MBDyn waits
		EDGE_BUSY		= 1,

		// EDGE waits (is ready to read kinematics); MBDyn iterates
		EDGE_READ_READY		= 2,

		// EDGE is computing; MBDyn waits before reading forces
		EDGE_MBDYN_WRITE_DONE	= 3,

		// EDGE converged; MBDyn advances one step
		EDGE_GOTO_NEXT_STEP	= 4,

		// EDGE wants to end simulation
		EDGE_QUIT		= 5,

		// must be the last one
		EDGE_LAST
	};

	const char *EDGEcmd2str(int cmd) const;

	EDGEcmd CheckFlag(int& cnt);

public:
	ExtFileHandlerEDGE(std::string& fflagname,
		std::string& fdataname, int iSleepTime, int iPrecision);
	~ExtFileHandlerEDGE(void);

	virtual void AfterPredict(void);

	virtual std::ostream& Send_pre(SendWhen when);
	virtual void Send_post(SendWhen when);

	virtual std::istream& Recv_pre(void);
	virtual bool Recv_post(void);
};

/* ExtFileHandlerEDGE - end */

class DataManager;
class MBDynParser;

extern ExtFileHandlerBase *
ReadExtFileHandlerEDGE(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel);

#endif // EXTEDGE_H

