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
   double F[3]={Fx,Fy,Fz};
   int NumF=3;
   MBDyn_force_structural* FORCE[NumF];
   Id LF[NumF];
   Id I,J;
   I=Node;
   J=Jfloat;
   Vec3 Relative_Arm;
   
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
   for (int i=0;i<3;i++) {
	   (DF[i])[i]=1;
	   DF[i].REF=RefM;
   }

   /* Insert data in action force statements */
   for (int i=0;i<NumF;i++) {
      LF[i]=GetFreeLabel (MBForces);
      FORCE[i]= new MBDyn_force_structural (LF[i],MBDyn_force::CONSERVATIVE,
					    WhichNode,DF[i],Relative_Arm,F[i]);
      
      sprintf (comment,
	"\n# Action Force %d, direction %d of marker %d, relative to Adams VFORCE %d",
	       LF[i],i+1,Rm,label);      
      if (stf[i]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[i]);
      }
      if (stf[3]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[3]);
      }
      FORCE[i]->Remark(comment);      
      MBForces.insert (MBDyn_entry(LF[i],(MBDyn_card*) FORCE[i]));
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
   for (int i=0;i<NumF;i++) {
      LF[i]=GetFreeLabel (MBForces);
      FORCE[i]= new MBDyn_force_structural (LF[i],MBDyn_force::CONSERVATIVE,
					    WhichNode,DF[i],Relative_Arm,-F[i]);
      
      sprintf (comment,
	"\n# Reaction Force %d, direction %d of marker %d, relative to Adams VFORCE %d",
	       LF[i],i+1,Rm,label);
      strcat (comment,"\n# WARNING! This is translation of reaction force applied on FLOATING marker");
      if (stf[i]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[i]);
      }
      if (stf[3]!=NULL) {
	 strcat (comment,"\n# [formula used for magnitude]:");
	 strcat (comment,stf[3]);
      }
      FORCE[i]->Remark(comment);      
      MBForces.insert (MBDyn_entry(LF[i],(MBDyn_card*) FORCE[i]));
   }

   return;
}
