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

/* socket driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

/* include del programma */

#include <rtai_out_elem.h>

/* RTAIOutElem - begin */

RTAIOutElem::RTAIOutElem(unsigned int uL, unsigned int nmb,
		ScalarDof *& pn)
: Elem(uL, Elem::RTAI_OUTPUT, flag(0)),
NumChannels(nmb), pNodes(pn)
{
	NO_OP;
}

RTAIOutElem::~RTAIOutElem(void)
{
	SAFEDELETEARR(pNodes);
}

std::ostream&
RTAIOutElem::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

Elem::Type
RTAIOutElem::GetElemType(void) const
{
	return Elem::RTAI_OUTPUT;
}

void
RTAIOutElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
RTAIOutElem::AssRes(SubVectorHandler& WorkVec, double dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
RTAIOutElem::AssJac(VariableSubMatrixHandler& WorkMat, double dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

void
RTAIOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	for (unsigned int i; i < NumChannels; i++) {
		/* assign value somewhere into mailbox buffer */
		(void)pNodes[i].pNode->dGetDofValue(1, pNodes[i].iOrder);
	}
}

#endif /* USE_RTAI */

