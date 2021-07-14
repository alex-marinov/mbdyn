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

//                                 SFORCE.CC                                  

#include <sforce.h>
#include <output.h>

extern MBDyn_deck MBReference;
extern MBDyn_deck MBForces;
extern MTP_deck Marker_Table;
extern PTN_deck PartToNode_Table;
extern PTR_deck Part_Table;


s_sforce::s_sforce(void) : Mode(_TRANSLATION),Value(0),
                           Actiononly(N),
                           _Actiononly(N),_Node1(N),_Node2(N),
                           _Value(N),_Mode(N)
                           { Node[0]=0,Node[1]=0; }

s_sforce::~s_sforce(void) {}

inline const char* const s_sforce::Gettype(void) const
{
   return "SFORCE";
}

Boolean s_sforce::Test(void)
{
   const int err_before=nerr;
   /* controlla la presenza dei marker */
   if (_Node1==N) out_error (32,"I");
   if (_Node2==N) out_error (32,"J");
   /* controlla che sia stato definita una modalità di funzionamento */
   if (_Mode==N) out_error (33,"");
   if (err_before != nerr) return Y; else return N;
}

ostream& s_sforce::Print (ostream& out) const
{
   out << endl;
   out << "SFORCE:" << label << endl;
   out << "      " << "Node 1 [" << _Node1 << "] = " << Node[0] << endl
       << "      " << "Node 2 [" << _Node2 << "] = " << Node[1] << endl
       << "      " << "Mode [" << _Mode << "] = " << Mode << endl
       << "      " << "Action Only? " << Actiononly << endl
       << "      " << "Function [" << _Value << "] = " << Value << endl;
   out << endl;
   Display_Formula(out);
   out << endl;
   return out;
}

void s_sforce::Translate (ostream& out)
{
   char* comment = new char[600];
   char* title = new char[80];
   
   /* Retrieve formulas from map */
   char* stf;
   p_formula_entry id;
   id=trova(Recipient,"FUNCTION");
   if ((*id).second != NULL) stf=new char[strlen((*id).second)]; 
   stf=(*id).second;
   
   RVec3 Force_Direction,Diff_Vector;
   Vec3  Relative_Arm(0,0,0);
   MBDyn_reference *MarkerI,*MarkerJ,*MI,*MJ;
   // double Magnitude;
   Id I,J,ReftoJ,ReftoI;
   I=Node[0];
   J=Node[1];
   MBDyn_drive_CONST* Magnitude_FA = new MBDyn_drive_CONST (Value);
   MBDyn_drive_CONST* Magnitude_FB = new MBDyn_drive_CONST (-Value);
   Magnitude_FA->Remark("/* drive */ ");
   Magnitude_FB->Remark("/* drive */ ");
   
   /* Resolve references */
   Id WhichPart, WhichNode, WhichReference;
   MTP Ex;
   p_MTP_entry p1;
   p_PTN_entry p2;
   p_PTR_entry p3;
   p1 = Marker_Table.find(I); /* Trova la parte corrispondente al marker I */
   Ex = ((*p1).second);
   WhichPart = Ex.Num; /* Numero ID della parte relativa al marker I */
   p2 = PartToNode_Table.find(Ex.Num); 
   p3 = Part_Table.find (Ex.Num);
   WhichNode = (*p2).second; /* Nodo corrispondente alla parte Adams */
   WhichReference = (*p3).second; /* Reference corrispondente al BCS sopra */
  
   /* Pointer debugging */
   MI=(MBDyn_reference*) Find_MBCard (I,MBReference);
   MJ=(MBDyn_reference*) Find_MBCard (J,MBReference);
   CHECK_AND_DEBUG (MJ,MarkerJ,J,MBReference);
   CHECK_AND_DEBUG (MI,MarkerI,I,MBReference);
   
   if (Actiononly==Y)
     {
	/* Forza nel marker I diretta secondo l'asse Z del Marker J */ 
		
	/* Definisce la direzione della forza */
	RVec3 Zdirection(Vec3(0,0,1),MBDyn_entity::REFERENCED,I);	
	Relative_Arm = MarkerI->Abs_Pos;
	ReftoI = GetFreeLabel (MBForces,label); /* Label per la forza */
	/* If is a torque */
	if (Mode==_ROTATION) {
	MBDyn_force_couple* FI = 
	  new MBDyn_force_couple (ReftoI,MBDyn_force::CONSERVATIVE,
				  WhichNode,Zdirection,Magnitude_FA);
	sprintf (comment,"[Action Only] Torque %d due to Adams SFORCE %d",
		 ReftoI,label);
	if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
	}
	FI->Remark(comment);
	sprintf (title,"Torque due to Adams SFORCE %d",label);
	FI->Title(title);
	   
	MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));

	}
	/* If is a force */
	if (Mode==_TRANSLATION) {
	   MBDyn_force_structural* FI = 
	     new MBDyn_force_structural (ReftoI,
					 MBDyn_force_structural::CONSERVATIVE,
					 WhichNode,Zdirection,Relative_Arm,
					 Magnitude_FA);
	   sprintf (comment,"[Action Only] Force %d due to Adams SFORCE %d",
		    ReftoI,label);
	if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
	}
	sprintf (title,"Force due to Adams SFORCE %d",label);
	FI->Title(title);
	FI->Remark(comment);

	MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));

	}
     }
   else
     {
       /* Forze di azione, reazione marker I e J
       * La direzione è data dalla congiungente i marker I e J
       * I valori della forza da applicare sono ovviamente uguali
       * e opposti ai due marker */
      
      Diff_Vector.Set(0,0,1);
      Relative_Arm = MarkerI->Abs_Pos;
      // Diff_Vector=(MarkerJ->Abs_Pos)-(MarkerI->Abs_Pos);
      // Force_Direction=Direction (Diff_Vector);
      Force_Direction=(Vec3(0,0,1));
      
      /* ACTION FORCE */
      ReftoI = GetFreeLabel (MBForces,label);
      
      if (Mode==_TRANSLATION) {
	 
      MBDyn_force_structural* FI = 
	new MBDyn_force_structural (ReftoI,
				    MBDyn_force_structural::CONSERVATIVE,
				    WhichNode,Force_Direction,Relative_Arm,
				    Magnitude_FA);
      sprintf (comment,"Action Force %d due to Adams SFORCE %i",
	       ReftoI,label);
      sprintf (title,"Force [action] due to Adams SFORCE %i",label);
      strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
      if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
      }
      FI->Remark(comment);
      FI->Title(title);
      /* Inserimento della forza */
      MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));
	 
      }
      if (Mode==_ROTATION) {
	 
	 MBDyn_force_couple* FI =
	   new MBDyn_force_couple (ReftoI,
				   MBDyn_force::CONSERVATIVE,
				   WhichNode,Force_Direction,
				   Magnitude_FA);
	 sprintf (title,"Torque [action] due to Adams SFORCE %i",label);
	 sprintf (comment,"Action Torque %d due to Adams SFORCE %i",
		  ReftoI,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 FI->Remark(comment);
	 FI->Title(title);
	 /* Inserimento della coppia I */
	 MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));
	 
      }

      /* REACTION FORCE */
      p1 = Marker_Table.find(J); /* Trova la parte corrispondente al marker J */
      Ex = ((*p1).second);
      WhichPart = Ex.Num; /* Numero ID della parte relativa al marker J */
      p2 = PartToNode_Table.find(Ex.Num); 
      p3 = Part_Table.find (Ex.Num);
      WhichNode = (*p2).second; /* Nodo corrispondente alla parte Adams */
      WhichReference = (*p3).second; /* Reference corrisp al BCS sopra */
      Relative_Arm = MarkerJ->Abs_Pos;
      
      ReftoJ = GetFreeLabel (MBForces);

      if (Mode==_TRANSLATION) {
	 
	 MBDyn_force_structural* FJ =
	   new MBDyn_force_structural (ReftoJ,
				    MBDyn_force_structural::CONSERVATIVE,
				    WhichNode,Force_Direction,Relative_Arm,
				    Magnitude_FB);
      
	 sprintf (comment,"Reaction Force %d due to Adams SFORCE %i",
		  ReftoJ,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 sprintf (title,"Force [reaction] due to Adams SFORCE %i",label);
	 FJ->Remark(comment);
	 FJ->Title(title);
	 
	 /* Inserimento della forza di reazione */
	 MBForces.insert (MBDyn_entry(ReftoJ,(MBDyn_card*) FJ));

      }

      if (Mode==_ROTATION) {
	 
	 MBDyn_force_couple* FJ =
	   new MBDyn_force_couple (ReftoJ,
				   MBDyn_force_structural::CONSERVATIVE,
				   WhichNode,Force_Direction,
				   Magnitude_FB);
      
	 sprintf (title,"Torque [reaction] due to Adams SFORCE %i",label);
	 sprintf (comment,"Reaction Torque %d due to Adams SFORCE %i",
		  ReftoJ,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 FJ->Title(title);
	 FJ->Remark(comment);
	 /* Inserimento della forza di reazione */
	 MBForces.insert (MBDyn_entry(ReftoJ,(MBDyn_card*) FJ));

      }

   }  
   return;
}
