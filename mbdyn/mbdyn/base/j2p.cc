/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#include <mbconfig.h>

#include <mynewmem.h>
#include <j2p.h>


/* Elem2Param - begin */

Elem2Param::Elem2Param(unsigned int uL, const DofOwner* pDO, flag fOut)
: ParameterNode(uL, pDO, 0., fOut), pElem(NULL), iNum(0) 
{
   NO_OP;
}


Elem2Param::~Elem2Param(void)
{
   NO_OP;
}


void Elem2Param::Bind(const Elem* pEl, unsigned int i)
{
   if (pElem != NULL) {
      cerr << "Elem2Param::Bind(): parameter (" << GetLabel()
	<< ") is already bound to " 
	<< psElemNames[pElem->GetElemType()] 
	<< " (" << pElem->GetLabel() << ')' << endl;
      THROW(ErrGeneric());
   }
   
   if (i <= 0) {
      cerr << "Elem2Param::Bind(): illegal value " << i 
	<< " for " << psElemNames[pEl->GetElemType()] 
	<< " (" << pEl->GetLabel() << ") private data" << endl;
      THROW(ErrGeneric());
   }
   
   iNum = i;     
   
   if (iNum <= 0 || iNum > pEl->iGetNumPrivData()) {
      cerr << "Elem2Param::Bind(): wrong element private data number "
	<< iNum << " for " << psElemNames[pEl->GetElemType()]
	<< " (" << pEl->GetLabel() << ')' << endl;
      THROW(ErrGeneric());
   }
   pElem = (Elem*)pEl;
   
   dGetDofValue(1, 0);
}


/* Contributo del nodo al file di restart */
ostream& Elem2Param::Restart(ostream& out) const
{
   return out << "Elem2Param is not implemented yet" << endl;
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void Elem2Param::SetDofValue(const doublereal& /* dValue */ ,
			     unsigned int /* iDof */ ,
			     unsigned int /* iOrder */)
{
   NO_OP;
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


void StrainGageParam::Bind(const Elem* pEl, unsigned int i)
{
   ASSERT(pEl != NULL);
   if (pEl->GetElemType() != ElemType::BEAM) {
      cerr << "StrainGageParam::Bind(): must bind to a beam" << endl;
      THROW(ErrGeneric());
   }
   
   /* Nota: ora i == 1 o 2 per punto di valutazione I o II */
   Elem2Param::Bind(pEl, i);
}


/* Contributo del nodo al file di restart */
ostream& StrainGageParam::Restart(ostream& out) const
{
   return out << "StrainGageParam is not implemented yet" << endl;
}

/* StrainGageParam - end */
