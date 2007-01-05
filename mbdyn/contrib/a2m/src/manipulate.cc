/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


// LIBRERIA DI MANIPOLAZIONE DELLE CARD DI MBDYN

#include <manipulate.h>

extern MBDyn_deck MBReference;

void Ref_Interp (MBDyn_reference* P, Id& M1, Id& M2)
{
   MBDyn_reference* Q[2];
   Q[0] = (MBDyn_reference*) Find_MBCard (M1,MBReference);
   Q[1] = (MBDyn_reference*) Find_MBCard (M2,MBReference);
   if ((Q[0]==NULL) | (Q[1]==NULL)) {
      cout << "In interpolatory formula beetween markers " << M1 << " "
	<< M2 << ", a marker in question does not exist" << endl;
   }
   else {
      Ref_Interp (P,Q[0],Q[1]);
   }
   return;
}


void Ref_Interp (MBDyn_reference* P, MBDyn_reference* P1,
		 MBDyn_reference* P2)
{
   P->Abs_Pos=Interp(P1->Abs_Pos,P2->Abs_Pos);
   P->Abs_Rot_Matrix=Interp(P1->Abs_Rot_Matrix,P2->Abs_Rot_Matrix);
   P->Abs_Vel=Interp(P1->Abs_Vel,P2->Abs_Vel);
   P->Abs_Omega=Interp(P1->Abs_Omega,P2->Abs_Omega);
   return;
}

void NodeSet (MBDyn_node_structural* NS, Id LNR)
{
   MBDyn_reference* NR;
   NR=(MBDyn_reference*) Find_MBCard (LNR,MBReference);
   if ((NR==NULL)) {
      cout << "DEBUG: [MANIPULATE.CC->NODESET] A reference to marker " <<
	LNR << " was done but it does not exist" << endl;
   }
   else NodeSet (NS,NR);
   return;
}

void NodeSet (MBDyn_node_structural* NS, MBDyn_reference* NR)
{
   /* Routine che preleva i valori da un reference per utilizzarli
    * nella definizione di un nodo. Utilità principale per le routine
    * BEAM */
   // NR->Restart(cout);
   NS->Mode_type=MBDyn_node_structural::DYNAMIC;
   NS->Abs_Pos=NR->Abs_Pos;
   NS->Abs_Rot_Matrix=NR->Abs_Rot_Matrix;
   NS->Abs_Vel=NR->Abs_Vel;
   NS->Abs_Ang_Vel=NR->Abs_Omega;
   NS->position_initial_stiffness=0;
   NS->velocity_initial_stiffness=0;
   NS->Omega_rotates=N;
   return;
}
