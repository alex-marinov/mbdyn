/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

typedef Real (*ModelFunc_0args_t)(DataManager *);
typedef Real (*ModelFunc_1args_t)(DataManager *, Real);
typedef Real (*ModelFunc_2args_t)(DataManager *, Real, Real);

/*
 * Computes the distance between two structural nodes
 */
static Real
distance(DataManager *pDM, Real n1, Real n2)
{
	unsigned uLabel1 = unsigned(n1);
	unsigned uLabel2 = unsigned(n2);

	StructNode *pNode1 = pDM->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("distance(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = pDM->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("distance(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	return (pNode2->GetXCurr() - pNode1->GetXCurr()).Norm();
}

static MathFunc_t	ModelFunc[] = {
	{ "distance",	2,	{ (MathFunc_0args_t)distance },	0,	"" },

     /* add more as needed */
	{ 0,		0,	{ (MathFunc_0args_t)0 },	0,	0 }
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

