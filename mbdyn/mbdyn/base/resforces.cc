/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <resforces.h>
#include <dataman.h>

/* ResForces - begin */

ResForces::ResForces(void) 
: F(0.), C(0.)
{ 
	NO_OP; 
}

ResForces::~ResForces(void) 
{
	NO_OP;
}

void
ResForces::Reset(void)
{ 
	F = 0.; 
	C = 0.;
}

void
ResForces::AddForce(const Vec3& f)
{
	F += f;
}

void
ResForces::AddForce(const Vec3& f, Vec3& x)
{ 
	F += f; 
	C += (x-Pole()).Cross(f); 
}

void
ResForces::AddCouple(const Vec3& c)
{
	C += c;
}

void
ResForces::AddForces(const Vec3& f, const Vec3& c, const Vec3& x)
{
	F += f;
	C += c + (x-Pole()).Cross(f);
}

void
ResForces::PutForce(const Vec3& f)
{
	F = f;
}

void
ResForces::PutCouple(const Vec3& c)
{
	C = c;
}

void
ResForces::PutForces(const Vec3& f, const Vec3& c)
{
	F = f;
	C = c;
}

const Vec3& 
ResForces::Force(void) const
{
	return F;
}

const Vec3&
ResForces::Couple(void) const
{
	return C;
}

ExternResForces::ExternResForces(void)
: X(0.)
{
	NO_OP;
}

ExternResForces::~ExternResForces(void)
{
	NO_OP;
}

void
ExternResForces::Reset(void)
{
	ResForces::Reset();
}

void
ExternResForces::Reset(const Vec3& x)
{
	X = x;
	ResForces::Reset();
}

void
ExternResForces::PutPole(const Vec3& x)
{
	X = x;
}

const Vec3& 
ExternResForces::Pole(void) const
{
	return X;
}

NodeResForces::NodeResForces(StructNode *n)
: pNode(n)
{
	NO_OP;
}

NodeResForces::~NodeResForces(void)
{
	NO_OP;
} 

const Vec3&
NodeResForces::Pole(void) const
{
	ASSERT(pNode); 
	return pNode->GetXCurr(); 
}

LocalNodeResForces::LocalNodeResForces(StructNode *n)
: NodeResForces(n)
{
	NO_OP;
}

LocalNodeResForces::~LocalNodeResForces(void)
{
	NO_OP;
}

const Vec3&
LocalNodeResForces::Force(void) const
{ 
	return (Fr = pNode->GetRCurr().Transpose()*F); 
}

const Vec3&
LocalNodeResForces::Couple(void) const
{ 
	return (Cr = pNode->GetRCurr().Transpose()*C);
}

ResForceSet::ResForceSet(unsigned int uLabel, ResForces *p) 
: WithLabel(uLabel), pRes(p)
{
	ASSERT(pRes);
}

ResForceSet::~ResForceSet(void)
{
	SAFEDELETE(pRes);
}

bool
ResForceSet::is_in(unsigned int uL)
{ 
	return labelSet.find(uL) != labelSet.end(); 
}

ResForceSet *
ReadResSet(DataManager* pDM, MBDynParser& HP, unsigned int uL)
{
	ResForceSet *pset = NULL;
	ResForces *pres = NULL;

	if (HP.IsKeyWord("external")) {
		SAFENEW(pres, ExternResForces);
		
	} else if (HP.IsKeyWord("node")) {
		StructNode *pNode = 
			(StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);
		
		if (pNode == NULL) {
			THROW(DataManager::ErrGeneric());
		}

		if (HP.IsKeyWord("local")) {
			SAFENEWWITHCONSTRUCTOR(pres, LocalNodeResForces, 
					LocalNodeResForces(pNode));

		} else {
			SAFENEWWITHCONSTRUCTOR(pres, NodeResForces, 
					NodeResForces(pNode));
		}
		
	} else {
		std::cerr << "unknown force set type at line " 
			<< HP.GetLineData() << std::endl;
		THROW(DataManager::ErrGeneric());
	}

	unsigned int nItems = HP.GetInt();
	if (nItems < 1) {
		std::cerr << "illegal number of items " << nItems 
			<< " in set at line " << HP.GetLineData() << std::endl;
		SAFEDELETE(pres);
		THROW(DataManager::ErrGeneric());
	}

	SAFENEWWITHCONSTRUCTOR(pset, ResForceSet, ResForceSet(uL, pres));

	for (unsigned int i = 0; i < nItems; i++) {
		unsigned int uLabel = HP.GetInt();
		std::pair<std::set<unsigned int>::iterator, bool> rc = 
			pset->labelSet.insert(uLabel);

		if (!rc.second) {
			std::cerr << "unable to insert item " << uLabel 
				<< " in set at line " << HP.GetLineData() 
				<< std::endl;
			SAFEDELETE(pset);
			THROW(DataManager::ErrGeneric());
		}
	}

	return pset;
}

ResForceSet **
ReadResSets(DataManager* pDM, MBDynParser& HP)
{
	ResForceSet **ppres = NULL;

	if (HP.IsKeyWord("set")) {
		int nSets = HP.GetInt();

		if (nSets < 1) {
			std::cerr << "Illegal number of sets (" << nSets 
				<< ") at line " << HP.GetLineData() 
				<< std::endl;
			THROW(DataManager::ErrGeneric());
		}

		SAFENEWARR(ppres, ResForceSet*, nSets+1);
		for (unsigned int i = 0; i <= (unsigned int)nSets; i++) {
			ppres[i] = NULL;
		}

		try {
			for (unsigned int i = 0; i < (unsigned int)nSets; i++) {
				ppres[i] = ReadResSet(pDM, HP, i+1);
			}

		} catch (...) {
			for (unsigned int i = 0; ppres[i]; i++) {
				SAFEDELETE(ppres[i]);
			}
			SAFEDELETEARR(ppres);
		}
	}

	return ppres;
}

/* ResForces - end */
