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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <elem.h>

/* Elem - begin */

Elem::Elem(unsigned int uL, ElemType::Type T, flag fOut)
: WithLabel(uL), ToBeOutput(fOut), ElemT(T)
{
   ASSERTMSG(uL > 0, "Null label shouldn't be used");
}


Elem::~Elem(void) 
{
   NO_OP;
}


/* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
unsigned int Elem::iGetNumDof(void) const
{
   return 0;
}

      
/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order Elem::SetDof(unsigned int /* i */ ) const
{
   ASSERTMSG(0, "You shouldn't have called this function");
   return DofOrder::UNKNOWN;
}

   
/* assemblaggio matrici per autovalori */
void Elem::AssEig(VariableSubMatrixHandler& /* WorkMatA */ ,
		  VariableSubMatrixHandler& /* WorkMatB */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) 
{
   cerr << "Sorry, not implemented yet" << endl;     
   THROW(ErrNotImplementedYet());
}


/* Setta i valori iniziali delle variabili (e fa altre cose) 
 * prima di iniziare l'integrazione */
void Elem::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{
   NO_OP;
}


/* Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator */
void Elem::BeforePredict(VectorHandler& /* X */ ,
			 VectorHandler& /* XP */ ,
			 VectorHandler& /* XPrev */ ,
			 VectorHandler& /* XPPrev */ ) const 
{
   NO_OP;
}


void Elem::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ ) 
{
   NO_OP;
}

         
/* Aggiorna dati in base alla soluzione */
void Elem::Update(const VectorHandler& /* XCurr */ , 
		  const VectorHandler& /* XPrimeCurr */ ) 
{
   NO_OP;
}


/* Metodi per l'estrazione di dati "privati".
 * Si suppone che l'estrattore li sappia interpretare.
 * Come default non ci sono dati privati estraibili */
unsigned int Elem::iGetNumPrivData(void) const 
{
   return 0;
}


doublereal Elem::dGetPrivData(unsigned int /* i */ ) const
{
   return 0.;
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

ElemWithDofs::ElemWithDofs(unsigned int uL, ElemType::Type T, 
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
					 ElemType::Type T, 
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
