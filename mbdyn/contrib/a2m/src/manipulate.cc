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
   //Per ora effettua solo un'interpolazione bruta, priva di riferimento
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
   NS->Abs_pos=NR->Abs_Pos;
   NS->Abs_rot_matrix=NR->Abs_Rot_Matrix;
   NS->Abs_vel=NR->Abs_Vel;
   NS->Abs_ang_vel=NR->Abs_Omega;
   NS->position_initial_stiffness=0;
   NS->velocity_initial_stiffness=0;
   NS->Omega_rotates=N;
   return;
}
