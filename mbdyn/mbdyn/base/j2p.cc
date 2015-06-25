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

#include "mynewmem.h"
#include "j2p.h"


/* Elem2Param - begin */

Elem2Param::Elem2Param(unsigned int uL, const DofOwner* pDO, flag fOut)
: ParameterNode(uL, pDO, 0., fOut), pElem(0), iNum(0) 
{
	NO_OP;
}

Elem2Param::~Elem2Param(void)
{
	NO_OP;
}


void
Elem2Param::Bind(const Elem* pEl, unsigned int i)
{
	if (pElem != 0) {
		silent_cerr("Elem2Param::Bind(): parameter (" << GetLabel()
			<< ") is already bound to "
			<< psElemNames[pElem->GetElemType()] 
			<< " (" << pElem->GetLabel() << ')' << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (i <= 0) {
		silent_cerr("Elem2Param::Bind(): illegal value " << i 
			<< " for " << psElemNames[pEl->GetElemType()] 
			<< " (" << pEl->GetLabel() << ") private data"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	iNum = i;     

	if (iNum <= 0 || iNum > pEl->iGetNumPrivData()) {
		silent_cerr("Elem2Param::Bind(): "
			"wrong element private data number "
			<< iNum << " for " << psElemNames[pEl->GetElemType()]
			<< " (" << pEl->GetLabel() << ')' << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pElem = pEl;

	dGetDofValue(1, 0);
}

/* Contributo del nodo al file di restart */
std::ostream&
Elem2Param::Restart(std::ostream& out) const
{
	//return out << "# Elem2Param is not implemented yet" << std::endl;
	return out << "  parameter: "
			<< GetLabel() << " ,element;" << std::endl;
}

std::ostream&
Elem2Param::RestartBind(std::ostream& out) const
{
	//return out << "# Elem2Param is not implemented yet" << std::endl;
	return out << "  bind: "
			<< pElem ->GetLabel() <<", "
			<< psReadElemsElems[pElem -> GetElemType()] << ", "
			<< GetLabel() << ", " << iNum << ";"
			<< std::endl;
}

/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void
Elem2Param::SetDofValue(const doublereal& /* dValue */ ,
		unsigned int /* iDof */ ,
		unsigned int /* iOrder */)
{
	NO_OP;
}

void
Elem2Param::SetValue(DataManager *pDM,
		VectorHandler& /* X */, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	if (pElem == 0) {
		silent_cerr("ParameterNode(" << GetLabel() 
				<< "): never bound" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	dX = dGetX();
}

/* Elem2Param - end */


/* StrainGageParam - begin */

StrainGageParam::StrainGageParam(unsigned int uL, const DofOwner* pDO,
		doublereal dy, doublereal dz, flag fOut)
: Elem2Param(uL, pDO, fOut), dY(dy), dZ(dz) 
{
	NO_OP;
}

StrainGageParam::~StrainGageParam(void)
{
	NO_OP;
}

void
StrainGageParam::Bind(const Elem* pEl, unsigned int i)
{
	ASSERT(pEl != 0);
	if (pEl->GetElemType() != Elem::BEAM) {
		silent_cerr("StrainGageParam::Bind(): must bind to a beam"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Nota: ora i == 1 o 2 per punto di valutazione I o II */
	Elem2Param::Bind(pEl, i);
}

/* Contributo del nodo al file di restart */
std::ostream&
StrainGageParam::Restart(std::ostream& out) const
{
	return out << "StrainGageParam is not implemented yet" << std::endl;
}

/* StrainGageParam - end */

