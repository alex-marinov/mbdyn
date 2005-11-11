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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "modelns.h"
#include "strnode.h"
#include <matvecexp.h>
#include <Rot.hh>

typedef Real (*ModelFunc_0args_t)(DataManager *);
typedef Real (*ModelFunc_1args_t)(DataManager *, Real);
typedef Real (*ModelFunc_2args_t)(DataManager *, Real, Real);

enum IDX_t {
	IDX1 = 1,
	IDX2 = 2,
	IDX3 = 3,
	NORM,
	SQUARE
};

static const char *
IDX2str(IDX_t IDX)
{
	switch (IDX) {
	case NORM:
		return "";

	case SQUARE:
		return "2";

	case IDX1:
		return "<1>";

	case IDX2:
		return "<2>";

	case IDX3:
		return "<3>";
	}
}

/*
 * Computes the position of a structural node
 */
template <IDX_t IDX>
static Real
position(DataManager *pDM, Real n)
{
	unsigned uLabel = unsigned(n);

	StructNode *pNode = pDM->pFindStructNode(uLabel);
	if (pNode == 0) {
		silent_cerr("position" << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	switch (IDX) {
	case NORM:
		return pNode->GetXCurr().Norm();

	case SQUARE:
		return pNode->GetXCurr().Dot();

	default:
		return pNode->GetXCurr()(IDX);
	}
}

/*
 * Computes the distance between two structural nodes
 */
template <IDX_t IDX>
static Real
distance(DataManager *pDM, Real n1, Real n2)
{
	unsigned uLabel1 = unsigned(n1);
	unsigned uLabel2 = unsigned(n2);

	StructNode *pNode1 = pDM->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("distance" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = pDM->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("distance" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 d = pNode2->GetXCurr() - pNode1->GetXCurr();

	switch (IDX) {
	case NORM:
		return d.Norm();

	case SQUARE:
		return d.Dot();

	default:
		return d(IDX);
	}
}

/*
 * Computes the distance between two structural nodes
 */
template <IDX_t IDX>
static Real
anglerel(DataManager *pDM, Real n1, Real n2)
{
	unsigned uLabel1 = unsigned(n1);
	unsigned uLabel2 = unsigned(n2);

	StructNode *pNode1 = pDM->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("angle" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = pDM->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("angle" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 phi(RotManip::VecRot(pNode2->GetRCurr()*pNode1->GetRCurr().Transpose()));

	switch (IDX) {
	case NORM:
		return phi.Norm();

	default:
		return phi(IDX);
	}
}

/*
 * Computes the velocity of a structural node
 */
template <IDX_t IDX>
static Real
velocity(DataManager *pDM, Real n)
{
	unsigned uLabel = unsigned(n);

	StructNode *pNode = pDM->pFindStructNode(uLabel);
	if (pNode == 0) {
		silent_cerr("velocity" << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	switch (IDX) {
	case NORM:
		return pNode->GetVCurr().Norm();

	case SQUARE:
		return pNode->GetVCurr().Dot();

	default:
		return pNode->GetVCurr()(IDX);
	}
}

/*
 * Computes the relative velocity between two structural nodes
 */
template <IDX_t IDX>
static Real
vrel(DataManager *pDM, Real n1, Real n2)
{
	unsigned uLabel1 = unsigned(n1);
	unsigned uLabel2 = unsigned(n2);

	StructNode *pNode1 = pDM->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("vrel" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = pDM->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("vrel" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 d = pNode2->GetVCurr() - pNode1->GetVCurr();

	switch (IDX) {
	case NORM:
		return d.Norm();

	case SQUARE:
		return d.Dot();

	default:
		return d(IDX);
	}
}

static MathFunc_t	ModelFunc[] = {
	{ "position",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)position<NORM>) },	0,	"" },
	{ "position2",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)position<SQUARE>) },	0,	"" },
	{ "xposition",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)position<IDX1>) },	0,	"" },
	{ "yposition",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)position<IDX2>) },	0,	"" },
	{ "zposition",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)position<IDX3>) },	0,	"" },

	{ "distance",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)distance<NORM>) },	0,	"" },
	{ "distance2",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)distance<SQUARE>) },	0,	"" },
	{ "xdistance",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)distance<IDX1>) },	0,	"" },
	{ "ydistance",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)distance<IDX2>) },	0,	"" },
	{ "zdistance",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)distance<IDX3>) },	0,	"" },

	{ "anglerel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)anglerel<NORM>) },	0,	"" },
	{ "xanglerel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)anglerel<IDX1>) },	0,	"" },
	{ "yanglerel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)anglerel<IDX2>) },	0,	"" },
	{ "zanglerel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)anglerel<IDX3>) },	0,	"" },

	{ "velocity",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)velocity<NORM>) },	0,	"" },
	{ "velocity2",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)velocity<SQUARE>) },	0,	"" },
	{ "xvelocity",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)velocity<IDX1>) },	0,	"" },
	{ "yvelocity",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)velocity<IDX2>) },	0,	"" },
	{ "zvelocity",	1,	{ (MathFunc_0args_t)((ModelFunc_1args_t)velocity<IDX3>) },	0,	"" },

	{ "vrel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)vrel<NORM>) },		0,	"" },
	{ "vrel2",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)vrel<SQUARE>) },	0,	"" },
	{ "xvrel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)vrel<IDX1>) },		0,	"" },
	{ "yvrel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)vrel<IDX2>) },		0,	"" },
	{ "zvrel",	2,	{ (MathFunc_0args_t)((ModelFunc_2args_t)vrel<IDX3>) },		0,	"" },

     /* add more as needed */
	{ 0,		0,	{ (MathFunc_0args_t)0 },		0,	0 }
};

ModelNameSpace::ModelNameSpace(DataManager *pdm)
: MathParser::StaticNameSpace("model", ModelFunc), pDM(pdm)
{
	NO_OP;
}

ModelNameSpace::~ModelNameSpace(void)
{
	NO_OP;
}

TypedValue 
ModelNameSpace::EvalFunc(MathFunc_t *f, Real *d) const
{
	switch (f->nargs) {
	case 0:
		return TypedValue((*((ModelFunc_0args_t)(f->f.f0)))(pDM));

	case 1:
		return TypedValue((*((ModelFunc_1args_t)(f->f.f1)))(pDM, d[0]));

	case 2:
		return TypedValue((*((ModelFunc_2args_t)(f->f.f2)))(pDM, d[0], d[1]));

	default:
		throw ErrGeneric();
	}
}

