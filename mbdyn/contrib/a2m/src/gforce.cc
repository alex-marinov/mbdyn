//                                  GFORCE.CC                                 

#include <gforce.h>
#include <output.h>

#include <string.h>

extern MBDyn_deck MBForces;
extern MBDyn_deck MBReference;

extern MTP_deck Marker_Table;
extern PTN_deck PartToNode_Table;
extern PTR_deck Part_Table;


s_gforce::s_gforce(void) :    Node(0),Jfloat(0),Rm(0),Fx(0),Fy(0),Fz(0),
                              Tx(0),Ty(0),Tz(0),Function(0),
                              _Jfloat(N),_Rm(N),_Fx(N),_Fy(N),_Fz(N),
                              _Tx(N),_Ty(N),_Tz(N),_Function(N),_Node(N)
                              {}

s_gforce::~s_gforce(void) {}

inline const char* const s_gforce::Gettype(void) const
{
   return "GFORCE";
}

Boolean s_gforce::Test(void)
{
   /* Test del metodo: FUNCTION o FX,FY,FZ,TX,TY,TZ */
   const int err_before=nerr;
   if (_Function==Y) {
      if (_Fx==Y) out_error (30,"FX");
      if (_Fy==Y) out_error (30,"FY");
      if (_Fz==Y) out_error (30,"FZ");
      if (_Tx==Y) out_error (30,"TX");
      if (_Ty==Y) out_error (30,"TY");
      if (_Tz==Y) out_error (30,"TZ");
   }
   /* Controllo dei parametri importanti per la definizione della card */
   if (_Jfloat==N) out_error (31,"JFLOAT");
   if (_Rm==N) out_error (31,"RM");
   if (_Node==N) out_error (31,"I MARKER");
   if (err_before != nerr) return Y; else return N;
}

ostream& s_gforce::Print (ostream& out) const
{
   out << endl;
   out << "GFORCE:" << label << endl;
   out << "      Marker I [" << _Node << "] = " << Node << endl
       << "      Jfloat [" << _Jfloat << "] = " << Jfloat << endl
       << "      Rm [" << _Rm << "] = " << Rm << endl
       << "      Fx [" << _Fx << "] = " << Fx << endl
       << "      Fy [" << _Fy << "] = " << Fy << endl
       << "      Fz [" << _Fz << "] = " << Fz << endl
       << "      Tx [" << _Tx << "] = " << Tx << endl
       << "      Ty [" << _Ty << "] = " << Ty << endl
       << "      Tz [" << _Tz << "] = " << Tz << endl
       << "      Function [" << _Function << "] = " << Function << endl;
   out << endl;
   Display_Formula(out);
   out << endl;
   return out;
}

void s_gforce::Translate(ostream& out)
{
   //
   // IL MARKER J SI MUOVE, (FLOATING) RISPETTO ALLA SUA PARTE, IN MODO
   // CHE LA POSIZIONE SPAZIALE (ASSOLUTA?) COINCIDA CON QUELLA DEL MARKER
   // I - IN PRATICA I PUNTI DI AZIONE E REAZIONE DEVONO COINCIDERE !
   //
   /* Per ora, come nel caso di sforce si confondo i Marker con i nodi
    * strutturali. Poi bisognerà sostituire tutto! Attenzione ! */
   char* comment = new char[600];

   /* Retrieve formulas from map */
   char* stf[7];
   p_formula_entry* id = new p_formula_entry[7];
   id[0]=trova(Recipient,"FX"); id[1]=trova(Recipient,"FY");
   id[2]=trova(Recipient,"FZ"); id[3]=trova(Recipient,"TX");
   id[4]=trova(Recipient,"TY"); id[5]=trova(Recipient,"TZ");
   id[6]=trova(Recipient,"FUNCTION");
   for (int k=0;k<7;k++) {
      if ((*id[k]).second != NULL) stf[k]=new char[strlen((*id[k]).second)]; 
      stf[k]=(*id[k]).second;
   }

   
   Id I,J;
   I = Node;
   J = Jfloat;
   Vec3 Relative_Arm[2];
   /* INSERIMENTO DELLE COMPONENTI DI FORZA */
   double F[3]={Fx,Fy,Fz};
   double T[3]={Tx,Ty,Tz};
   int NumF = 3;
   int NumT = 3;
   Id LF[2*NumF];
   Id LT[2*NumT];
   MBDyn_force_structural* FORCE[2*NumF];
   MBDyn_force_couple* COUPLE[2*NumT];

   MBDyn_reference* ACTION_MARKER,*REACTION_MARKER,*REF_MARKER,*P1,*P2,*P3;
   P1=(MBDyn_reference*) Find_MBCard (Node,MBReference);
   P2=(MBDyn_reference*) Find_MBCard (Rm,MBReference);
   P3=(MBDyn_reference*) Find_MBCard (Jfloat,MBReference);

   /* Pointer debugging */
   CHECK_AND_DEBUG (P1,ACTION_MARKER,Node,MBReference);
   CHECK_AND_DEBUG (P2,REF_MARKER,Rm,MBReference);
   CHECK_AND_DEBUG (P3,REACTION_MARKER,Jfloat,MBReference);
   
   /* The rotation matrix of marker Rm gives the direction of forces,
    * and, therefore, the direction of torque. */
   RVec3 DF[3];

   /* Resolves references */
   MTP Ex;
   p_MTP_entry p1;
   p_PTN_entry p2;
   p_PTR_entry p3;
   Id WhichPart[2],WhichNode[2],WhichReference[2];
   p1 = Marker_Table.find(I);
   Ex = ((*p1).second);
   WhichPart[0] = Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode[0] = (*p2).second;
   WhichReference[0] = (*p3).second;
   Relative_Arm[0] = P1->Abs_Pos;
   
   p1 = Marker_Table.find(J);
   Ex = ((*p1).second);
   WhichPart[1] = Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode[1] = (*p2).second;
   WhichReference[1] = (*p3).second;
   Relative_Arm[1] = P1->Abs_Pos;
   
   MBDyn_entity RefM (MBDyn_entity::REFERENCED,REF_MARKER->Label);
   for (int i=0;i<3;i++) {
	   (DF[i])[i]=1;
	   DF[i].REF=RefM;
   }

   for (int i=0;i<NumF;i++)
     {
	/* FORCE INSERTION */
	
	/* ACTION FORCE INSERT */
	LF[(2*i)]=GetFreeLabel (MBForces);
	FORCE[(2*i)]=new MBDyn_force_structural (LF[(2*i)],MBDyn_force::CONSERVATIVE,
					WhichNode[0],DF[i],Relative_Arm[0],F[i]);
	sprintf (comment,
		 "\n# Action Force %d, direction %d of marker %d, relative to Adams GFORCE %d",
		 LF[i],i+1,Rm,label);
	
	if (stf[i]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i]);
	}
	if (stf[6]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat(comment,stf[6]);
	}
	FORCE[(2*i)]->Remark(comment);
	
	MBForces.insert (MBDyn_entry(LF[(2*i)],(MBDyn_card*) FORCE[(2*i)]));	
	
	/* REACTION FORCE INSERT */
	LF[(2*i)+1]=GetFreeLabel (MBForces);
	FORCE[(2*i)+1]=new MBDyn_force_structural(LF[(2*i)+1],
						  MBDyn_force::CONSERVATIVE,
						  WhichNode[1],DF[i],
						  Relative_Arm[1],-F[i]);
	sprintf (comment,
		 "\n# Reaction Force %d, direction %d of marker %d, relative to Adams GFORCE %d",
		 LF[(2*i)+1],i+1,Rm,label);
	if (stf[i]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i]);
	}
	if (stf[6]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[6]);
	}
	FORCE[(2*i)+1]->Remark(comment);
	
	MBForces.insert (MBDyn_entry(LF[(2*i)+1],(MBDyn_card*) FORCE[(2*i)+1]));
	
     }
   
   for (int i=0;i<NumT;i++)
     {
	/* TORQUE INSERTION */
	
	/* ACTION TORQUE INSERT */
	LT[(2*i)]=GetFreeLabel (MBForces);
	COUPLE[(2*i)]=new MBDyn_force_couple (LT[(2*i)],
					      MBDyn_force::CONSERVATIVE,
					      WhichNode[0],DF[i],T[i]);
	sprintf (comment,
		 "\n# Action Torque %d, direction %d of marker %d, relative to Adams GFORCE %d",
		 LT[(2*i)],i+1,Rm,label);
	if (stf[i+3]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i+3]);
	}
	if (stf[6]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat(comment,stf[6]);
	}
	COUPLE[(2*i)]->Remark(comment);
	
	MBForces.insert (MBDyn_entry(LT[(2*i)],(MBDyn_card*) COUPLE[(2*i)]));
	
	/* REACTION TORQUE INSERT */
	Id idx = (2*i)+1;
	LT[idx]= GetFreeLabel (MBForces);
	COUPLE[idx]=new MBDyn_force_couple (LT[idx],
					    MBDyn_force::CONSERVATIVE,
					    WhichNode[1],DF[i],-T[i]);
	sprintf (comment,
		 "\n# Reaction Torque %d, direction %d of marker %d, relative to Adams GFORCE %d",
		 LT[idx],i+1,Rm,label);
	if (stf[i+3]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i+3]);
	}
	if (stf[6]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat(comment,stf[6]);
	}
	COUPLE[idx]->Remark(comment);
	
	MBForces.insert (MBDyn_entry(LT[idx],(MBDyn_card*) COUPLE[idx]));

     }
   /* END OF FORCES AND COUPLE INSERTION */
   return;
}
