/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef SCALARVALUE_H
#define SCALARVALUE_H

#include "elem.h"

/* ScalarValue - begin */

class ScalarValue {
public:
	virtual ~ScalarValue(void);
	virtual doublereal dGetValue(void) const = 0;
};

class ScalarDofValue : public ScalarValue, public ScalarDof {
public:
	ScalarDofValue(const ScalarDof& sd);
	doublereal dGetValue(void) const;
};

class ScalarDriveValue : public ScalarValue {
protected:
	const DriveCaller *pDC;

public:
	ScalarDriveValue(const DriveCaller *pdc);
	~ScalarDriveValue(void);
	doublereal dGetValue(void) const;
};

class DataManager;
class MBDynParser;

extern ScalarValue *
ReadScalarValue(DataManager *pDM, MBDynParser& HP);

extern void
ReadScalarValues(DataManager *pDM, MBDynParser& HP,
	std::vector<ScalarValue *>& Values);

/* ScalarValue - end */

#endif /* SCALARVALUE_H */

