/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2000

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

Copyright (C) 1999-2000
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
   char* comment = new char [160];
   char* title = new char[80];
   // Traduzione : TRAVE ADAMS = STRUCT NODE + BEAM MBDYN
   pCL6D pK1=new ConstitutiveLaw6D(Emodulus,Gmodulus);
   pCL6D pK2=new ConstitutiveLaw6D(Emodulus,Gmodulus);
   /* Dati da inserire nella formula della BEAM */
   Vec3 F1,F2,F3;
   
   /* R0 e R1 sono le matrici di rotazione relative al materiale
    rispetto al sistema globale o cosa? verificare. Per adesso ..*/
   Mat3x3 R[2]; for (int j=0;j<2;j++) R[j]=Eye3;
   
   Id I,J;
   MBDyn_reference* MI, *MJ, *MarkerI, *MarkerJ;
   Vec3 OffsetL,OffsetC,OffsetR;
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
   
   /* Load Markers */
   MI = (MBDyn_reference*) Find_MBCard (I,MBReference);
   MJ = (MBDyn_reference*) Find_MBCard (J,MBReference);
   CHECK_AND_DEBUG (MJ,MarkerJ,J,MBReference);
   CHECK_AND_DEBUG (MI,MarkerI,I,MBReference);
   
   OffsetL=Unref(MarkerI->Abs_Pos);
   OffsetR=Unref(MarkerJ->Abs_Pos);
   OffsetC=Interp(OffsetL,OffsetR,0.5);
   /* ATTENZIONE : LA PARTE DI CODICE SOTTOSTANTE E' VALIDA SE :
    * 1. I BCS DELLE PARTI SONO RIFERITI, ATTRAVERSO IL SEMPLICE
    * MARKER "BCS", AL SISTEMA GLOBALE */

   /* Si crea un nodo strutturale interpolato tra i due e si interpola
    * l'offset relativo */
   /* WhichReference sono già i due BCS */
   
   /* Creazione del reference BCS */
   Id CENTR = GetFreeLabel (MBReference);
   MBDyn_reference* RFR = new MBDyn_reference(CENTR);
   Ref_Interp (RFR,WhichReference[0],WhichReference[1]);
   sprintf (title,"BCS for central node of generated BEAM");
   RFR->Title(title);
   sprintf (comment,"Reference %i generated by interpolation between BCS %i and BCS %i",
	    CENTR,WhichReference[0],WhichReference[1]);
   MBReference.insert (MBDyn_entry (CENTR, (MBDyn_card*) RFR));
   RFR->Remark(comment);

   /* Inserisce la label tra i reference generati dalle Beam [ADAMS] */
   Beam_reference_Table.insert (BTR_entry(label,CENTR));
   
   /* Creazione del nodo strutturale intermedio */
   MBDyn_entity M (MBDyn_entity::REFERENCED,CENTR);
   RVec3 RDefPos (Vec3(0,0,0),M);
   RMat3x3 RDefRot (Eye3,M);
   Id NodeID = GetFreeLabel (MBNodes);
   MBDyn_node_structural* CENTRALNODE = new 
     MBDyn_node_structural (NodeID,MBDyn_node_structural::DYNAMIC,
			    RDefPos,RDefRot,RDefPos,RDefPos,0,0,Y);   
   sprintf (comment,"Central node automatically generated for Adams BEAM %i",
	    label);
   CENTRALNODE->Remark(comment);
   sprintf (title,"central node for BEAM %i",label);
   CENTRALNODE->Title(title);
   MBNodes.insert ( MBDyn_entry (NodeID,(MBDyn_card*) CENTRALNODE));

   
   // Creazione della beam attaccata al nodo centrale.
   // Al solito verifica che ci sia una label libera per la Beam.
   Id IDBEAM = GetFreeLabel (MBBeams,label);
   MBDyn_beam* CENTRALBEAM = new
     MBDyn_beam (IDBEAM,WhichNode[0],NodeID,WhichNode[1],
		 OffsetL,OffsetC,OffsetR,R[0],R[1],pK1,pK2);
   MBBeams.insert(MBDyn_entry(IDBEAM,(MBDyn_card*) CENTRALBEAM));

   sprintf (title,"Adams BEAM %i",label);
   sprintf (comment,"MBDyn Beam %i relative to Adams BEAM %i",
	    IDBEAM,label);   
   CENTRALBEAM->Remark(comment);
   CENTRALBEAM->Title(title);
   
   return;
}

