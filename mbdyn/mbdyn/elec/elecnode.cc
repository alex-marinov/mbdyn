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

#include <mynewmem.h>
#include <elecnode.h>
#include <solman.h>

/* AbstractNode - begin */

/* Output del nodo */
void AbstractNode::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {      
#ifdef DEBUG   
      OH.Output() << "Abstract node " << GetLabel() 
	<< ": x = " << dX << ", x' = " << dXP << endl;
#endif
   
      OH.Abstract() << setw(8) << GetLabel() << " "
	<< dX << " " << dXP << endl;
   }
}


/* Aggiorna dati in base alla soluzione */
void AbstractNode::Update(const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
   integer iFirstIndex = iGetFirstIndex()+1;
   dX = XCurr.dGetCoef(iFirstIndex);
   dXP = XPrimeCurr.dGetCoef(iFirstIndex);
}


/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
void AbstractNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   integer iFirstIndex = iGetFirstIndex()+1;
   X.fPutCoef(iFirstIndex, dX);
   XP.fPutCoef(iFirstIndex, dXP);
}

/* AbstractNode - end */



/* ElectricNode - begin */

/* Output del nodo */
void ElectricNode::Output(OutputHandler& OH) const 
{
   if(fToBeOutput()) {      
#ifdef DEBUG   
      OH.Output() << "Electric node " << uLabel << ':' << endl 
	<< "sorry, not implemented yet" << endl;
#endif   
   }   
}

/* ElectricNode - end */


