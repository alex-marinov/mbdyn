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
   char* comment = new char[255];

   /* Retrieve formulas from map */
   char* stf;
   p_formula_entry id;
   id=trova(Recipient,"FUNCTION");
   if ((*id).second != NULL) stf=new char[strlen((*id).second)]; 
   stf=(*id).second;
   
   RVec3 Force_Direction,Diff_Vector;
   Vec3  Relative_Arm(0,0,0);
   MBDyn_reference *MarkerI,*MarkerJ,*MI,*MJ;
   double Magnitude;
   Id I,J,ReftoJ,ReftoI;
   I=Node[0];
   J=Node[1];

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
   CHECK_AND_DEBUG (MJ,MarkerJ,I,MBReference);
   CHECK_AND_DEBUG (MI,MarkerI,J,MBReference);
   
   if (Actiononly==Y)
     {
	/* Forza nel marker I diretta secondo l'asse Z del Marker J */ 
		
	/* Definisce la direzione della forza */
	RVec3 Zdirection(Vec3(0,0,1),MBDyn_entity::REFERENCED,I);
	Magnitude=Value;
	Relative_Arm = MarkerI->Abs_Pos;
	ReftoI = GetFreeLabel (MBForces,label); /* Label per la forza */
	/* If is a torque */
	if (Mode==_ROTATION) {
	MBDyn_force_couple* FI = 
	  new MBDyn_force_couple (ReftoI,MBDyn_force::CONSERVATIVE,
				  WhichNode,Zdirection,Magnitude);
	sprintf (comment,"\n# [Action Only] Torque %d due to Adams SFORCE %d",
		 ReftoI,label);
	if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
	}
	FI->Remark(comment);

	MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));

	}
	/* If is a force */
	if (Mode==_TRANSLATION) {
	   MBDyn_force_structural* FI = 
	     new MBDyn_force_structural (ReftoI,
					 MBDyn_force_structural::CONSERVATIVE,
					 WhichNode,Zdirection,Relative_Arm,
					 Magnitude);
	   sprintf (comment,"\n# [Action Only] Force %d due to Adams SFORCE %d",
		    ReftoI,label);
	if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
	}
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
      Magnitude = Value;
      
      /* ACTION FORCE */
      ReftoI = GetFreeLabel (MBForces,label);
      
      if (Mode==_TRANSLATION) {
	 
      MBDyn_force_structural* FI = 
	new MBDyn_force_structural (ReftoI,
				    MBDyn_force_structural::CONSERVATIVE,
				    WhichNode,Force_Direction,Relative_Arm,
				    Magnitude);
      sprintf (comment,"\n# Action Force %d due to Adams SFORCE %i",
	       ReftoI,label);
      strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
      if (stf!=NULL) {
	   strcat (comment,"\n# [formula used for magnitude]:");
	   strcat (comment,stf);
      }
      FI->Remark(comment);
      /* Inserimento della forza */
      MBForces.insert (MBDyn_entry(ReftoI, (MBDyn_card*) FI));
	 
      }
      if (Mode==_ROTATION) {
	 
	 MBDyn_force_couple* FI =
	   new MBDyn_force_couple (ReftoI,
				   MBDyn_force::CONSERVATIVE,
				   WhichNode,Force_Direction,
				   Magnitude);
	 sprintf (comment,"\n# Action Torque %d due to Adams SFORCE %i",
		  ReftoI,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 FI->Remark(comment);
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
				    -Magnitude);
      
	 sprintf (comment,"\n# Reaction Force %d due to Adams SFORCE %i",
		  ReftoJ,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 FJ->Remark(comment);
	 /* Inserimento della forza di reazione */
	 MBForces.insert (MBDyn_entry(ReftoJ,(MBDyn_card*) FJ));

      }

      if (Mode==_ROTATION) {
	 
	 MBDyn_force_couple* FJ =
	   new MBDyn_force_couple (ReftoJ,
				   MBDyn_force_structural::CONSERVATIVE,
				   WhichNode,Force_Direction,
				   -Magnitude);
      
	 sprintf (comment,"\n# Reaction Torque %d due to Adams SFORCE %i",
		  ReftoJ,label);
	 strcat (comment,"\n# SFORCE IN ACTION-REACTION FORM NOT SUPPORTED - SPECIFY DIRECTION MANUALLY!");
	 if (stf!=NULL) {
	    strcat (comment,"\n# [formula used for magnitude]:");
	    strcat (comment,stf);
	 }
	 FJ->Remark(comment);
	 /* Inserimento della forza di reazione */
	 MBForces.insert (MBDyn_entry(ReftoJ,(MBDyn_card*) FJ));

      }

   }  
   return;
}
