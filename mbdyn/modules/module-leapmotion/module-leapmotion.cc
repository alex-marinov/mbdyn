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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "dataman.h"
#include "filedrv.h"

#include "Leap.h"

class LMDrive : public FileDrive {
public:
	LMDrive(unsigned int uL, const DriveHandler* pDH,
		const std::string& sFileName,
		const std::vector<doublereal>& v0);

	virtual ~LMDrive(void);

	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

LMDrive::LMDrive(unsigned int uL, const DriveHandler* pDH,
	const std::string& sFileName,
	const std::vector<doublereal>& v0)
: FileDrive(uL, pDH, sFileName, 1, v0)
{
	NO_OP;
}

LMDrive::~LMDrive(void)
{
	NO_OP;
}

std::ostream&
LMDrive::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

void
LMDrive::ServePending(const doublereal& t)
{
	NO_OP;
}

struct LMDR : public DriveRead {
	virtual Drive *
	Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP);
};

Drive *
LMDR::Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	leapmotion						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"	file: <label> , leapmotion ,					\n"
"           ...                                                         \n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	const std::string sFileName(HP.GetFileName());

	const std::vector<doublereal> v0;

	return new LMDrive(uLabel, pDM->pGetDrvHdl(), sFileName, v0);
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	DriveRead *rf = new LMDR;

	if (!SetDriveData("leapmotion", rf)) {
		delete rf;

		silent_cerr("LMDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

