//                                  GFORCE.CC                                 

#include <vtorque.h>
#include <output.h>


extern MBDyn_deck MBReference;
extern MBDyn_deck MBForces;

extern MTP_deck Marker_Table;
extern PTN_deck PartToNode_Table;
extern PTR_deck Part_Table;

s_vtorque::s_vtorque(void) :  Node(0),Jfloat(0),Rm(0),
                              Tx(0),Ty(0),Tz(0),Function(0),
                              _Jfloat(N),_Rm(N),
                              _Tx(N),_Ty(N),_Tz(N),_Function(N)
                              {}

s_vtorque::~s_vtorque(void) {}

inline const char* const s_vtorque::Gettype(void) const
{
   return "VTORQUE";
}

Boolean s_vtorque::Test(void)
{
   const int err_before=nerr;
   /* Controllo dei parametri fondamentali */
   if (_Rm==N) { out_error (42,""); Rm=Node; }
   if (_Node==N) out_error (43,"");
   if (_Jfloat==N) out_error (44,"");
   /* Controlla che sia definito almeno un valore */
   if ( (_Tx==N) & (_Ty==N) & (_Tz==N) )
     if (_Function==N) out_error (45,"");
   /* Controllo dei parametri rindondanti */
   if ( (_Tx==Y) | (_Ty==Y) | (_Tz==Y) )
     if (_Function==Y) out_error (46,"");
   if (err_before != nerr) return Y; else return N;   
}

ostream& s_vtorque::Print (ostream& out) const
{
   out << endl;
   out << "VTORQUE:" << label << endl;
   out << "      Marker I [" << _Node << "] = " << Node << endl
       << "      Jfloat [" << _Jfloat << "] = " << Jfloat << endl
       << "      Rm [" << _Rm << "] = " << Rm << endl
       << "      Tx [" << _Tx << "] = " << Tx << endl
       << "      Ty [" << _Ty << "] = " << Ty << endl
       << "      Tz [" << _Tz << "] = " << Tz << endl
       << "      Function [" << _Function << "] = " << Function << endl;
   out << endl;
   Display_Formula(out);
   return out;
}

void s_vtorque::Translate (ostream& out)
{
   char* comment= new char [255];

   /* Retrieve formulas from map */
   char* stf[4];
   p_formula_entry* id = new p_formula_entry[4];
   id[0]=trova(Recipient,"TX"); id[1]=trova(Recipient,"TY");
   id[2]=trova(Recipient,"TZ");
   id[3]=trova(Recipient,"FUNCTION");
   for (int k=0;k<4;k++) {
      if ((*id[k]).second != NULL) stf[k]=new char[strlen((*id[k]).second)]; 
      stf[k]=(*id[k]).second;
   }

   MBDyn_reference* ACTION_MARKER,*REACTION_MARKER,*REF_MARKER,*P1,*P2,*P3;
   P1=(MBDyn_reference*) Find_MBCard (Node,MBReference);
   P2=(MBDyn_reference*) Find_MBCard (Rm,MBReference);
   P3=(MBDyn_reference*) Find_MBCard (Jfloat,MBReference);

   /* Pointer debugging */
   CHECK_AND_DEBUG (P1,ACTION_MARKER,Node,MBReference);
   CHECK_AND_DEBUG (P2,REF_MARKER,Rm,MBReference);
   CHECK_AND_DEBUG (P3,REACTION_MARKER,Jfloat,MBReference);
         
   RVec3 DF[3];
   for (int i=0;i<3;i++) {
      (DF[i])[i]=1;
      DF[i].REF=MBDyn_entity(MBDyn_entity::REFERENCED,Rm);
   }

   int NumT = 3;
   Id LT[2*NumT];
   Id I,J;
   I=Node;
   J=Jfloat;
   Id idx1,idx2;
   double T[3]={Tx,Ty,Tz};
   MBDyn_force_couple* COUPLE[2*NumT];

   
   /* Resolves references */
   Id WhichPart[2],WhichNode[2],WhichReference[2];
   MTP Ex;
   p_MTP_entry p1;
   p_PTN_entry p2;
   p_PTR_entry p3;
   p1 = Marker_Table.find(I);
   Ex = ((*p1).second);
   WhichPart[0] = Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode[0] = (*p2).second;
   WhichReference[0] = (*p3).second;
   
   p1 = Marker_Table.find(J);
   Ex = ((*p1).second);
   WhichPart[1]= Ex.Num;
   p2 = PartToNode_Table.find(Ex.Num);
   p3 = Part_Table.find (Ex.Num);
   WhichNode[1] = (*p2).second;
   WhichReference[0] = (*p3).second;
   
   for (int i=0;i<NumT;i++)
     {
	idx1 = (2*i);
	idx2 = (2*i)+1;
	
	LT[idx1]=GetFreeLabel (MBForces);
	COUPLE[idx1]=new MBDyn_force_couple (LT[idx1],MBDyn_force::CONSERVATIVE,
					 WhichNode[0],DF[i],T[i]);
	sprintf (comment,"\n# Couple %d, direction %d of marker %d due to Adams VTORQUE %d",
		 LT[idx1],i+1,Rm,label);
	if (stf[i]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i]);
	}
	if (stf[3]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[3]);
	}
	COUPLE[idx1]->Remark(comment);
	MBForces.insert (MBDyn_entry(LT[idx1],(MBDyn_card*) COUPLE[idx1]));
	
	LT[idx2]=GetFreeLabel (MBForces);
	COUPLE[idx2]=new MBDyn_force_couple (LT[idx2],MBDyn_force::CONSERVATIVE,
					     WhichNode[1],DF[i],-T[i]);
	sprintf (comment,"\n# Reaction Couple %d, direction %d of marker %d due to Adams VTORQUE %d",
		 LT[idx2],i+1,Rm,label);
	if (stf[i]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[i]);
	}
	if (stf[3]!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf[3]);
	}
	COUPLE[idx2]->Remark(comment);
	MBForces.insert (MBDyn_entry(LT[idx2],(MBDyn_card*) COUPLE[idx2]));
	
	
     }
   /* END OF COUPLES INSERTION */
   return;
}
