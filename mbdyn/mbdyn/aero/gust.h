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

#ifndef GUST_H
#define GUST_H

#include "matvec3.h"

// forward declaration
class MBDynParser;
class DataManager;
class AirProperties;

/* Gust - begin */

class Gust {
protected:
	const AirProperties *pAP;

public:
	virtual ~Gust(void);
	void SetAirProperties(const AirProperties *pap);
	virtual Vec3 GetVelocity(const Vec3& X) const;
	virtual bool GetVelocity(const Vec3& X, Vec3& V) const = 0;
	virtual std::ostream& Restart(std::ostream& out) const = 0;
};

/* prototype of the functional object: reads a gust */
struct GustRead {
public:
	virtual ~GustRead(void);
	virtual Gust *
	Read(const DataManager* pDM, MBDynParser& HP) = 0;
};

/* drive caller registration function: call to register one */
extern bool
SetGustData(const char *name, GustRead* rf);

/* create/destroy */
extern void InitGustData(void);
extern void DestroyGustData(void);

extern Gust *
ReadGustData(const DataManager *pDM, MBDynParser& HP);

/* Gust - end */

#endif // GUST_H

