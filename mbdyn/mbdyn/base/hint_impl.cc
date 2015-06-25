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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "hint_impl.h"

Hint *
ParseHint(DataManager *pDM, const char *s)
{
	if (strncasecmp(s, "drive{", STRLENOF("drive{")) == 0) {
		s += STRLENOF("drive{");
		
		size_t	len = strlen(s);

		if (s[len - 1] != '}') {
			return 0;
		}

		char *sDriveStr = new char[len + 1];
		memcpy(sDriveStr, s, len + 1);
		sDriveStr[len - 1] = ';';

		return new DriveHint(sDriveStr);

	} else if (strncasecmp(s, "drive3{", STRLENOF("drive3{")) == 0) {
		s += STRLENOF("drive3{");
		
		size_t	len = strlen(s);

		if (s[len - 1] != '}') {
			return 0;
		}

		char *sDriveStr = new char[len + 1];
		memcpy(sDriveStr, s, len + 1);
		sDriveStr[len - 1] = ';';

		return new TplDriveHint<Vec3>(sDriveStr);

	} else if (strncasecmp(s, "drive6{", STRLENOF("drive6{")) == 0) {
		s += STRLENOF("drive6{");
		
		size_t	len = strlen(s);

		if (s[len - 1] != '}') {
			return 0;
		}

		char *sDriveStr = new char[len + 1];
		memcpy(sDriveStr, s, len + 1);
		sDriveStr[len - 1] = ';';

		return new TplDriveHint<Vec6>(sDriveStr);
	} 

	return 0;
}

/* ParsableHint - start */

ParsableHint::ParsableHint(const char *s)
: sHint((char *)s)
{
	NO_OP;
}

ParsableHint::~ParsableHint(void)
{
	if (sHint != 0) {
		SAFEDELETEARR(sHint);
	}
}

/* DriveHint - start */

DriveHint::DriveHint(const char *s)
: ParsableHint(s)
{
	NO_OP;
}

DriveHint::~DriveHint(void)
{
	NO_OP;
}

DriveCaller *
DriveHint::pCreateDrive(DataManager *pDM) const
{
	std::istringstream in(sHint);
	InputStream In(in);

	MBDynParser HP(pDM->GetMathParser(), In, "DriveHint::pCreateDrive");
	HP.ExpectArg();
	HP.SetDataManager(pDM);

	return HP.GetDriveCaller(false);
}

/* DriveHint - end */


/* TplDriveHint - start */

TplDriveHint3 tdh3(0);
TplDriveHint6 tdh6(0);

/* TplDriveHint - end */
