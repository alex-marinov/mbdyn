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

#include <elem.h>

/* Elem - begin */

Elem::Elem(unsigned int uL, Elem::Type T, flag fOut)
: WithLabel(uL), ToBeOutput(fOut), ElemT(T)
{
   ASSERTMSG(uL > 0, "Null label shouldn't be used");
}


Elem::~Elem(void) 
{
   NO_OP;
}


/* assemblaggio matrici per autovalori */
void Elem::AssMats(VariableSubMatrixHandler& /* WorkMatA */ ,
		  VariableSubMatrixHandler& /* WorkMatB */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) 
{
   std::cerr << psElemNames[GetElemType()] << "(" << GetLabel()
	   << "): AssMats() not implemented yet" << std::endl;
}


unsigned int 
Elem::iGetNumDof(void) const
{
	return 0;
}

DofOrder::Order 
Elem::GetDofType(unsigned int) const
{
	std::cerr << psElemNames[GetElemType()] << "(" << GetLabel()
		<< "): GetDofType() is undefined because element "
		"has no degrees of freedom" << std::endl;
	throw ErrGeneric();
}


Elem* Elem::pGetElem(void) const
{
   return (Elem*)this; 
}


ElemWithDofs* Elem::pGetElemWithDofs(void) const
{
   return NULL;
}


ElemGravityOwner* Elem::pGetElemGravityOwner(void) const
{
   return NULL;
}


AerodynamicElem* Elem::pGetAerodynamicElem(void) const
{
   return NULL;
}


InitialAssemblyElem* Elem::pGetInitialAssemblyElem(void) const
{
   return NULL;
}

/* Elem - end */


/* ElemWithDofs - begin */

ElemWithDofs::ElemWithDofs(unsigned int uL, Elem::Type T, 
			   const DofOwner* pDO, flag fOut)
: Elem(uL, T, fOut), DofOwnerOwner((DofOwner*)pDO)
{
   NO_OP;
}


ElemWithDofs::~ElemWithDofs(void)
{
   NO_OP;
}


/* Consente di effettuare un casting sicuro da Elem* a ElemWithDofs* */
ElemWithDofs* ElemWithDofs::pGetElemWithDofs(void) const
{
   return (ElemWithDofs*)this; 
}

/* ElemWithDofs - end */


/* SubjectToInitialAssembly - begin */

SubjectToInitialAssembly::SubjectToInitialAssembly(void) 
{
   NO_OP;
}


SubjectToInitialAssembly::~SubjectToInitialAssembly(void)
{
   NO_OP;
}

/* SubjectToInitialAssembly - end */


/* InitialAssemblyElem - begin */

InitialAssemblyElem::InitialAssemblyElem(unsigned int uL, 
					 Elem::Type T, 
					 flag fOut)
: Elem(uL, T, fOut)
{ 
   NO_OP;
}


InitialAssemblyElem::~InitialAssemblyElem(void)
{
   NO_OP;
}


/* Consente di effettuare un casting sicuro da Elem* a InitialAssemblyElem* */
InitialAssemblyElem* InitialAssemblyElem::pGetInitialAssemblyElem(void) const
{
   return (InitialAssemblyElem*)this; 
}
   
/* InitialAssemblyElem - end */
