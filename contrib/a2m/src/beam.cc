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


//                                 BEAM.CC
//
//Definizioni per la card Beam

#include <beam.h>

extern MBDyn_deck MBBeams;
extern MBDyn_deck MBReference;
extern MBDyn_deck MBNodes;
extern MBDyn_deck MBBodies;
extern MBDyn_deck MBConstLaw;

extern PTR_deck Part_Table;
extern MTP_deck Marker_Table;
extern PTN_deck PartToNode_Table;
extern BTR_deck Beam_reference_Table;


s_beam::s_beam(void)	: Asy(0), Asz(0), Cratio(0), Ixx(0),Iyy(0),Izz(0),
                          Area(0),Length(0),Emodulus(0),Gmodulus(0),
                          Cmatrix (new double[21]),
                          _Node1(N), _Node2(N), _Asy(N), _Asz(N), _Cratio(N) ,
                          _Ixx(N), _Iyy(N), _Izz(N), _Area(N), _Length(N),
                          _Emodulus(N), _Gmodulus(N),_Cmatrix(N)
{
   Node[0]=0; Node[1]=0; 
}

s_beam::~s_beam(void)
{
}

Boolean s_beam::Test()
{
   const int err_before=nerr;
   // CONTROLLA CHE SIANO PRESENTI I PARAMETRI NECESSARI
   if (_Node1==N)
     {
	out_error(6,"OF ID I");
     }
   if (_Node2==N)
     {
	out_error(6,"OF ID J");
     }
   if (_Length==N)
     {
	out_error(7,"");
     }
   if (_Ixx==N)
     {
	out_error(8,"IXX");
     }
   if (_Iyy==N)
     {
	out_error(8,"IYY");
     }
   if (_Izz==N)
     {
	out_error(8,"IZZ");
     }
   if (_Area==N)
     {
	out_error(8,"AREA");
     }
   if (_Emodulus==N)
     {
	out_error (9,"YOUNG MODULUS");
     }
   if (_Gmodulus==N)
     {
	out_error (9,"SHEAR MODULUS");
     }
   if ((_Cmatrix==Y) & (_Cratio==Y))
     {
	out_error (15,"");
     }
   if (err_before != nerr) return Y; else return N;
}

inline const char* const s_beam::Gettype(void) const
{
   return "BEAM";
}

ostream& s_beam::Print (ostream& out) const
{
   out << endl;
   out << "BEAM:" << label << endl ;
   out << "     " << "Node 1 [" << _Node1 << "] = "  << Node[0] << endl
       << "     " << "Node 2 [" << _Node2 << "] = "  << Node[1] << endl;
   out << "     " << "Length [" << _Length << "] = " << Length << endl;
   out << "     " << "Area [" << _Area << "] = " << Area << endl;
   out << "     " << "Asy [" << _Asy << "] = " << Asy << endl
     << "     " << "Asz [" << _Asz << "] = " << Asz << endl;
   out << "     " << "Ixx [" << _Ixx << "] = " << Ixx << endl
     << "     " << "Iyy [" << _Iyy << "] = " << Iyy << endl
     << "     " << "Izz [" << _Izz << "] = " << Izz << endl;
   out << "     " << "Emodulus [" << _Emodulus << "] = " << Emodulus << endl
     << "     " << "Gmodulus [" << _Gmodulus << "] = " << Gmodulus << endl
     << "     " << "Cratio [" << _Cratio << "] = " << Cratio << endl
     << "     " << "Cmatrix [" << _Cmatrix << "]" << endl;
   if (_Cmatrix==Y)
     {
	out << "     ";
	for (int i=0; i< 8; i++) out << Cmatrix[i] << " ";
	out << endl << "     " ;
	for (int i=8; i< 16; i++) out << Cmatrix[i] << " " ;
	out << endl << "     ";
	for (int i=16; i<21; i++) out << Cmatrix[i] << " " ;
	out << endl;
     }
   out << endl;
   return out;
}

void s_beam::Translate(ostream& out)
{
   double norma;
   Vec3 contrib;
   char* comment = new char [160];
   char* title = new char[80];
   // Traduzione : TRAVE ADAMS = STRUCT NODE + BEAM MBDYN
   /* Dati da inserire nella formula della BEAM */
   Vec3 F1,F2,F3;
   Mat3x3 R_parte_I, R_parte_J, R_nodo_C;
   
   /* R0 e R1 sono le matrici di rotazione relative al materiale
    rispetto al sistema globale o cosa? verificare. Per adesso ..*/
   Mat3x3 R[2]; for (int j=0;j<2;j++) R[j]=Eye3;
   
   Id I,J,NodeID;
   MBDyn_reference* MI, *MJ, *MarkerI, *MarkerJ, 
     *BCS_I, *BCS_J, *BCS_parte_I, *BCS_parte_J;
   MBDyn_node_structural* NI, *NJ, *NodeI, *NodeJ;
   Vec3 OffsetL,OffsetC,OffsetR;
   Vec3 X_I,X_J;
   Vec3 OffsetI,OffsetJ,OffsetM;
   Vec3 X_parte_I, X_parte_J, X_nodo_C;
   
   I = Node[0];
   J = Node[1];
   
   /* Resolve references */
   Id WhichPart[2],WhichNode[2],WhichReference[2];
   MTP Ex;
   p_MTP_entry p1;
   p_PTN_entry p2;
   p_PTR_entry p3;
   for (int i=0;i<2;i++)
     {
	p1 = Marker_Table.find (Node[i]);
	Ex = ((*p1).second);
	
	WhichPart[i] = Ex.Num;
	p2 = PartToNode_Table.find(Ex.Num);
	p3 = Part_Table.find (Ex.Num);
	WhichNode[i] = (*p2).second;
	WhichReference[i] = (*p3).second;
     }
   
   /* Load references relative to markers I and J */
   MI = (MBDyn_reference*) Find_MBCard (I,MBReference);
   MJ = (MBDyn_reference*) Find_MBCard (J,MBReference);
   /* Load BCS (Body Coordinate System) for each part */
   BCS_I = (MBDyn_reference*) Find_MBCard (WhichReference[0],MBReference);
   BCS_J = (MBDyn_reference*) Find_MBCard (WhichReference[1],MBReference);
   /* Load nodes */
   NI = (MBDyn_node_structural*) Find_MBCard (WhichNode[0],MBNodes);
   NJ = (MBDyn_node_structural*) Find_MBCard (WhichNode[1],MBNodes);
   /* Debug references and nodes */
   CHECK_AND_DEBUG (MJ,MarkerJ,J,MBReference);
   CHECK_AND_DEBUG (MI,MarkerI,I,MBReference);
   CHECK_AND_DEBUG (NI,NodeI,WhichNode[0],MBNodes);
   CHECK_AND_DEBUG (NJ,NodeJ,WhichNode[1],MBNodes);
   CHECK_AND_DEBUG (BCS_I,BCS_parte_I,WhichReference[0],MBNodes);
   CHECK_AND_DEBUG (BCS_J,BCS_parte_J,WhichReference[1],MBNodes);
   
   /* Posizioni riferite al sistema assoluto delle parti */
   X_parte_I = Unref (BCS_parte_I->Abs_Pos);
   X_parte_J = Unref (BCS_parte_J->Abs_Pos);   
   /* Le matrici di rotazione delle parti nel GCS */
   R_parte_I = Unref (BCS_parte_I->Abs_Rot_Matrix);
   R_parte_J = Unref (BCS_parte_J->Abs_Rot_Matrix);
   
   /* Offset dei marker nei sistemi di riferimento parte */
   OffsetL = Unref (MarkerJ->Abs_Pos);
   OffsetR = Unref (MarkerI->Abs_Pos);

   /* L'offset della parte a sinistra di trave e' gia' definito */
   OffsetJ = OffsetL;

   /* Posizione del marker J nel sistema assoluto */
   X_J = R_parte_J * OffsetL + X_parte_J;
   /* Posizione del marker I nel sistema assoluto */
   X_I = R_parte_I * OffsetR + X_parte_I;
   /* Norma del vettore differenza */
   norma = Module((X_I-X_J));
    
   /* Versore markerI->markerJ - > contributo relativo alla lunghezza */
   contrib = (((X_I-X_J)/norma)*Length);

   /* Offset del Marker I nel sistema di riferimento della parte I */
   OffsetI = R_parte_I.tr() * (X_J + contrib - X_parte_I);
   
   /* La matrice di rotazione del nodo centrale e' l'interpolazione tra
    * la parte I e la parte J , altrettanto dicasi la posizione */
   R_nodo_C = 0.5 * ( R_parte_I + R_parte_J );
   X_nodo_C = 0.5 * ( X_parte_I + X_parte_J );
   /* Offset del nodo centrale nel sistema di riferimento assoluto */
   OffsetC = 0.5 * (R_parte_I * OffsetI + R_parte_J * OffsetJ);
   /* Offset del nodo centrale nel sistema di riferimento nodo centrale */
   OffsetC = R_nodo_C.tr()*OffsetC;
   
   /* Referenziazione al sistema globale di posizione e matrice di rot */
   RVec3 RDefPos (X_nodo_C,MBDyn_entity(MBDyn_entity::NUL));
   RMat3x3 RDefRot (R_nodo_C,MBDyn_entity(MBDyn_entity::NUL));

   
   /* Creazione del reference BCS per il nodo centrale */
   Id CentrID = GetFreeLabel (MBReference);
   MBDyn_reference* RFR = new 
     MBDyn_reference(CentrID,RDefPos,RDefRot,Y,Zero3,Zero3);
   sprintf (title,"BCS for central node of generated BEAM");
   RFR->Title(title);
   sprintf (comment,"Reference %i generated by interpolation between BCS %i and BCS %i",
	    CentrID,WhichReference[0],WhichReference[1]);
   MBReference.insert (MBDyn_entry (CentrID, (MBDyn_card*) RFR));
   RFR->Remark(comment);
   /* Inserisce la label tra i reference generati dalle Beam [ADAMS] */
   Beam_reference_Table.insert (BTR_entry(label,CentrID));
   
   /* Creazione del nodo strutturale intermedio */
   MBDyn_entity RF(MBDyn_entity::REFERENCED,CentrID);
   RVec3 RNulPos(Zero3,RF);
   RMat3x3 RNulRot(Eye3,RF);
   NodeID = GetFreeLabel (MBNodes);
   MBDyn_node_structural* CENTRALNODE = new
     MBDyn_node_structural (NodeID,MBDyn_node_structural::DYNAMIC,
			    RNulPos,RNulRot,RNulPos,RNulPos,0,0,Y);
   sprintf (comment,"Central node automatically generated for Adams BEAM %i",
	    label);
   CENTRALNODE->Remark(comment);
   sprintf (title,"central node for BEAM %i",label);
   CENTRALNODE->Title(title);
   MBNodes.insert (MBDyn_entry (NodeID,(MBDyn_card*) CENTRALNODE));

   /* Referenziazione al sistema nodo dei 3 offset */
   MBDyn_entity STD(MBDyn_entity::NODE);
   RVec3 STDOffsetL (OffsetJ,STD);
   RVec3 STDOffsetR (OffsetI,STD);
   RVec3 STDOffsetC (OffsetC,STD);

   /* Leggi costitutive nelle sezioni */
   
   Mat6x6* pK[2];
   for (int i=0;i<2;i++) {
      pK[i]=new Mat6x6();
      *(pK[i])= KMatrix (Length,Ixx,Iyy,Izz,Area,Emodulus,
			 Gmodulus,Asy,Asz,NULL,0);
   }

   /* Beam attaccata ai nodi */
   Id BeamID = GetFreeLabel (MBBeams,label);
   MBDyn_beam* BEAM = new MBDyn_beam (BeamID,WhichNode[1],NodeID,
				      WhichNode[0],STDOffsetL,STDOffsetC,
				      STDOffsetR,R[1],R[0],pK[1],pK[0]);
   MBBeams.insert(MBDyn_entry(BeamID,(MBDyn_card*) BEAM));
 
   sprintf (title,"Adams BEAM %i",label);
   sprintf (comment,"MBDyn Beam %i relative to Adams BEAM %i",
	    BeamID,label);   
   BEAM->Remark(comment);
   BEAM->Title(title);

   
   /* End of new code */
   /* return; */
   
   return;
}

