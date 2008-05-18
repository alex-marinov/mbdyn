/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/*
 * Michele Attolico <attolico@aero.polimi.it>
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "streamoutelem.h"

/* ScalarValue - begin */

ScalarValue::~ScalarValue(void)
{
	NO_OP;
}

ScalarDofValue::ScalarDofValue(const ScalarDof& sd)
: ScalarDof(sd)
{
	NO_OP;
}

doublereal
ScalarDofValue::dGetValue(void) const
{
	return ScalarDof::dGetValue();
}

ScalarDriveValue::ScalarDriveValue(const DriveCaller *pdc)
: pDC(pdc)
{
	NO_OP;
}

ScalarDriveValue::~ScalarDriveValue(void)
{
	if (pDC) {
		delete pDC;
	}
}

doublereal
ScalarDriveValue::dGetValue(void) const
{
	return pDC->dGet();
}

/* ScalarValue - end */

/* StreamOutElem - begin */

StreamOutElem::StreamOutElem(unsigned int uL,
		std::vector<ScalarValue *>& pn,
		unsigned int oe)
: Elem(uL, flag(0)),
Values(pn), size(-1), buf(0),
OutputEvery(oe), OutputCounter(0)
{
	ASSERT(OutputEvery > 0);

	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*Values.size();
	SAFENEWARR(buf, char, size);
	memset(buf, 0, size);
}

StreamOutElem::~StreamOutElem(void)
{
	if (buf != 0) {
		SAFEDELETEARR(buf);
	}

	for (std::vector<ScalarValue *>::iterator i = Values.begin();
		i != Values.end(); i++)
	{
		delete *i;
	}
}

Elem::Type
StreamOutElem::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
StreamOutElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
StreamOutElem::AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
StreamOutElem::AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

