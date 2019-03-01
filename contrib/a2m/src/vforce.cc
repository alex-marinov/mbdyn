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

//                                  GFORCE.CC                                 

#include <vforce.h>
#include <output.h>

extern MBDyn_deck MBReference;
extern MBDyn_deck MBForces;

extern MTP_deck Marker_Table;
extern PTN_deck PartToNode_Table;
extern PTR_deck Part_Table;

s_vforce::s_vforce(void) :    Node(0),Jfloat(0),Rm(0),Fx(0),Fy(0),Fz(0),
                              Function(0),
                              _Jfloat(N),_Rm(N),_Fx(N),_Fy(N),_Fz(N),
                              _Function(N)
                              {}

s_vforce::~s_vforce(void) {}

inline const char* const s_vforce::Gettype(void) const
{
   return "VFORCE";
}

Boolean s_vforce::Test(void)
{
   const int err_before=nerr;
   /* Controllo dei parametri fondamentali */
   if (_Rm==N) { out_error (37,""); Rm=Node; }
   if (_Node==N) out_error (38,"");
   if (_Jfloat==N) out_error (39,"");
   /* Controlla che sia definito almeno un valore */
   if ( (_Fx==N) & (_Fy==N) & (_Fz==N) )
     if (_Function==N) out_error (40,"");
   /* Controllo dei parametri rindondanti */
   if ( (_Fx==Y) | (_Fy==Y) | (_Fz==Y) )
     if (_Function==Y) out_error (41,"");
   if (err_before != nerr) return Y; else return N;   
}

ostream& s_vforce::Print (ostream& out) const
{
   out << endl;
   out << "VFORCE:" << label << endl;
   out << "      Marker I [" << _Node << "] = " << Node << endl
       << "      Jfloat [" << _Jfloat << "] = " << Jfloat << endl
       << "      Rm [" << _Rm << "] = " << Rm << endl
       << "      Fx [" << _Fx << "] = " << Fx << endl
       << "      Fy [" << _Fy << "] = " << Fy << endl
       << "      Fz [" << _Fz << "] = " << Fz << endl
       << "      Function [" << _Function << "] = " << Function << endl;
   out << endl;
   Display_Formula(out);
   out << endl;
   return out;
}

void s_vforce::Translate (ostream& out)
{
   char *comment = new char [255];
   char *title = new char[80];
   
   /* Retrieve formulas from map */
   char* stf[4];
   p_formula_entry* id = new p_formula_entry[6];
   id[0]=trova(Recipient,"FX"); id[1]=trova(Recipient,"FY");
   id[2]=trova(Recipient,"FZ");
   id[3]=trova(Recipient,"FORMULA");
   for (int k=0;k<4;k++) {
      if ((*id[k]).second != NULL) stf[k]=new char[strlen((*id[k]).second)]; 
      stf[k]=(*id[k]).second;
   }
   
   MBDyn_reference* ACTION_MARKER,*REACTION_MARKER,*REF_MARKER,*P1,*P2,*P3;
   P1=(MBDyn_reference*) Find_MBCard (Node,MBReference);
   P2=(MBDyn_reference*) Find_MBCard (Rm,MBReference);
   P3=(MBDyn_reference*) Find_MBCard (Jfloat,MBReference);

   /* Pointers debugging */
   CHECK_AND_DEBUG (P1,ACTION_MARKER,Node,MBReference);
   CHECK_AND_DEBUG (P2,REF_MARKER,Rm,MBReference);
   CHECK_AND_DEBUG (P3,REACTION_MARKER,Jfloat,MBReference);
   
   RVec3 DF[3];
   double MF[6]={Fx,-Fx,Fy,-Fy,Fz,-Fz};
   MBDyn_force_structural* FORCE[6];
   Id LF[6];
   Id I,J,idx0,idx1;
   I=Node;
   J=Jfloat;
   Vec3 Relative_Arm;
   MBDyn_drive_CONST* F[6];
   for (int i=0;i<6;i++) {
      F[i]=new MBDyn_drive_CONST (MF[i]);
      F[i]->Remark("/* drive */ ");
   }
   
   /* Resolves references for marker I */
   Id WhichPart,WhichNode,WhichReference;
   MTP Ex;
   p_MTP_entry p1;
   p_PTN_entry p2;
   p_PTR_entry p3;
   p1 = Marker_Table.find(I);
   Ex = ((*p1).second);
   WhichPart = Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode = (*p2).second;
   WhichReference = (*p3).second;
   Relative_Arm = P1->Abs_Pos;
   MBDyn_entity RefM (MBDyn_entity::REFERENCED,REF_MARKER->Label);

   /* Resolves offsets */
   for (int i=0;i<3;i++) {
	   (DF[i])[i]=1;
	   DF[i].REF=RefM;
   }

   /* Insert data in action force statements */
   for (int i=0;i<3;i++) {
      
      idx0=(2*i); idx1=idx0+1;
      
      LF[idx0]=GetFreeLabel (MBForces);
      FORCE[idx0]= new 
	MBDyn_force_structural (LF[idx0],
				MBDyn_force::CONSERVATIVE,
				WhichNode,DF[i],Relative_Arm,F[idx0]);
      sprintf (comment,
	"Action Force %d, direction %d of marker %d, relative to Adams VFORCE %d",
	       LF[idx0],i+1,Rm,label);      
      if (stf[i]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[i]);
      }
      if (stf[3]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[3]);
      }
      sprintf (title,"F%c [action] due to Adams VFORCE %i",(88+i),label);
      
      FORCE[idx0]->Title(title);
      FORCE[idx0]->Remark(comment);      
      MBForces.insert (MBDyn_entry(LF[idx0],(MBDyn_card*) FORCE[idx0]));
   }
   
   /* Reaction Force - > Warning ! Floating marker not supported, so
    * forces are applied to a fixed empty marker */ 
   /* Resolve references for Marker J(float) */
   p1 = Marker_Table.find(J);
   Ex = ((*p1).second);
   WhichPart = Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode = (*p2).second;
   WhichReference = (*p3).second;
   Relative_Arm = P1->Abs_Pos;

   /* Insert data in reaction force statements */
   for (int i=0;i<3;i++) {

      idx1 = (2*i)+1;
      
      LF[idx1]=GetFreeLabel (MBForces);
      FORCE[idx1]= new 
	MBDyn_force_structural (LF[idx1],MBDyn_force::CONSERVATIVE,
				WhichNode,DF[i],Relative_Arm,F[idx1]);
      sprintf (title,"F%c [reaction] due to Adams VFORCE %i",(88+i),label);
      sprintf (comment,
	"Reaction Force %d, direction %d of marker %d, relative to Adams VFORCE %d",
	       LF[idx1],i+1,Rm,label);
      strcat (comment,"\n# WARNING! This is translation of reaction force applied on FLOATING marker");
      if (stf[i]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[i]);
      }
      if (stf[3]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[3]);
      }
      FORCE[idx1]->Title(title);
      FORCE[idx1]->Remark(comment);      
      MBForces.insert (MBDyn_entry(LF[idx1],(MBDyn_card*) FORCE[idx1]));
   }

   return;
}
