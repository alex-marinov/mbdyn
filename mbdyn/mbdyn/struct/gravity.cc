/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Elemento accelerazione di gravita' */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "gravity.h"

/* Gravity - begin */

Gravity::Gravity(const TplDriveCaller<Vec3>* pDC, 
		 flag fOut)
: Elem(1, fOut), TplDriveOwner<Vec3>(pDC)
{
   Acc = Get();
}


Gravity::~Gravity(void)
{
   NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream& Gravity::Restart(std::ostream& out) const
{
   return out << "  gravity: /*reference, global,*/ ",
     pGetDriveCaller()->Restart(out) << ";" << std::endl;
}


   /* assemblaggio jacobiano */
VariableSubMatrixHandler& Gravity::AssJac(VariableSubMatrixHandler& WorkMat,
					  doublereal /* dCoef */ , 
					  const VectorHandler& /* XCurr */ ,
					  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Gravity::AssJac()" << std::endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}
   

/* assemblaggio residuo */
SubVectorHandler& Gravity::AssRes(SubVectorHandler& WorkVec,
				  doublereal /* dCoef */ ,
				  const VectorHandler& /* XCurr */ , 
				  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Gravity::AssRes()" << std::endl);
   WorkVec.Resize(0);
   
   /* Approfitto del fatto che Gravity viene aggiornato prima 
    * degli altri elementi (vedi l'enum Elem::Type e la sequenza di
    * assemblaggio) per fargli calcolare Acc una volta per tutte.
    * Quindi, quando viene chiamata GetAcceleration(void), 
    * questa restituisce un reference all'accelerazione con il
    * minimo overhead */
   Acc = Get();
   return WorkVec;
}

void
Gravity::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::GRAVITY)) {
			OH.Gravity() << std::setw(8) << GetLabel()
				<< " " << Acc << std::endl;
		}
	}
}

/* Gravity - end */


/* GravityOwner - begin */

GravityOwner::GravityOwner(void)
: pGravity(NULL)
{
   NO_OP;
}


GravityOwner::~GravityOwner(void)
{
   NO_OP;
}


void GravityOwner::PutGravity(const Gravity* pG)
{
   ASSERT(pGravity == NULL);
   (Gravity*&)pGravity = (Gravity*)pG;
}


bool GravityOwner::bGetGravity(const Vec3& X, Vec3& Acc) const
{
   if(pGravity == NULL) {
      return false;
   }
   
   Acc = pGravity->GetAcceleration(X);
   return true;
}

/* GravityOwner - end */


/* ElemGravityOwner - begin */

ElemGravityOwner::ElemGravityOwner(unsigned int uL,
				   flag fOut)
: Elem(uL, fOut), GravityOwner()
{
   NO_OP;
}


ElemGravityOwner::~ElemGravityOwner(void)
{
   NO_OP;
}

/* ElemGravityOwner - end */
