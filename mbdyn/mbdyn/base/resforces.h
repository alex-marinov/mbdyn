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

#ifndef RESFORCES_H
#define RESFORCES_H

#include <cfloat>
#include <iostream>

#include <set>
#include "strnode.h"

/* ResForces - begin */

class ResForces {
protected:
	Vec3 F;
	Vec3 C;

public:
	ResForces(void);
	virtual ~ResForces(void);

	virtual void Reset(void);
	void AddForce(const Vec3& f);
	void AddForce(const Vec3& f, Vec3& x);
	void AddMoment(const Vec3& c);
	void AddForces(const Vec3& f, const Vec3& c, const Vec3& x);
	void PutForce(const Vec3& f);
	void PutMoment(const Vec3& c);
	void PutForces(const Vec3& f, const Vec3& c);
	
	virtual const Vec3& Force(void) const;
	virtual const Vec3& Moment(void) const;
	virtual const Vec3& Pole(void) const = 0;
};

class ExternResForces : public ResForces {
protected:
	Vec3 X;

public:
	ExternResForces(void);
	virtual ~ExternResForces(void);

	void Reset(void);
	void Reset(const Vec3& x);
	void PutPole(const Vec3& x);
	
	const Vec3& Pole(void) const;
};

class NodeResForces : public ResForces {
protected:
	const StructNode *pNode;

public:
	NodeResForces(const StructNode *n = 0);
	virtual ~NodeResForces(void);

	const Vec3& Pole(void) const;
};

class LocalNodeResForces : public NodeResForces {
protected:
	mutable Vec3 Fr;
	mutable Vec3 Cr;

public:
	LocalNodeResForces(const StructNode *n = 0);
	virtual ~LocalNodeResForces(void);

	const Vec3& Force(void) const;
	const Vec3& Moment(void) const;
};

struct ResForceSet : public WithLabel {
	ResForces *pRes;
	std::set<unsigned int> labelSet;
	
	ResForceSet(unsigned int uLabel, ResForces *p);
	virtual ~ResForceSet(void);
	bool is_in(unsigned int uL);
};

/* ResForces - end */

extern ResForceSet *
ReadResSet(DataManager* pDM, MBDynParser& HP, unsigned int uL);
extern ResForceSet **
ReadResSets(DataManager* pDM, MBDynParser& HP);

#endif /* RESFORCES_H */

