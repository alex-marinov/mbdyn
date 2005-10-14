/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

#ifndef HINT_IMPL_H
#define HINT_IMPL_H

/* include del programma */
#include "ac/sstream"
#include "dataman.h"
#include "tpldrive_.h"

extern Hint *
ParseHint(DataManager *pDM, const char *s);

class DriveHint : public Hint {
protected:
	char	*sDriveStr;

public:
	DriveHint(const char *s);
	virtual ~DriveHint(void);

	DriveCaller *pCreateDrive(DataManager *pDM) const;
};

template <class T>
class TplDriveHint : public DriveHint {
protected:
	T	t;

public:
	TplDriveHint(const T& t, const char *s);
	virtual ~TplDriveHint(void);
	TplDriveCaller<T> *pCreateDrive(DataManager *pDM) const;
};

template <class T>
TplDriveHint<T>::TplDriveHint<T>(const T& t, const char *s)
: DriveHint(s), t(t)
{
	NO_OP;
}

template <class T>
TplDriveHint<T>::~TplDriveHint<T>(void)
{
	NO_OP;
}

template <class T>
TplDriveCaller<T> *
TplDriveHint<T>::pCreateDrive(DataManager *pDM) const
{
#if defined(HAVE_SSTREAM)
	std::istringstream in(sDriveStr);
#else /* HAVE_STRSTREAM_H */
	istrstream in(sDriveStr);
#endif /* HAVE_STRSTREAM_H */
	InputStream In(in);

	MBDynParser HP(pDM->GetMathParser(), In, "TplDriveHint::pCreateDrive");
	HP.ExpectArg();

	return ReadTplDrive(pDM, HP, t);
}

typedef TplDriveHint<Vec3> TplDriveHint3;
typedef TplDriveHint<Vec6> TplDriveHint6;

#endif /* HINT_IMPL_H */

