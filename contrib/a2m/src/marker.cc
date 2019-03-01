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

//                               MARKER.CC                                    

// Definizioni per la card Marker

#include <marker.h>

#include <part.h>
#include <pointmass.h>
extern card_deck parts;
extern card_deck pointmasses;

extern MTP_deck Marker_Table;
extern MBDyn_deck MBReference;
extern PTR_deck Part_Table;
extern MTR_deck Reference_Table;
extern PTR_deck Pmass_Table;

s_marker::s_marker(void)  : Part(0), Point_mass(0),_UseXP(N),           
                           _Floating(N),Flex_body(0),Node_id(0),
                           Mode(s_marker::EYEMATRIX),
                           _Part(N), _Point_mass(N),
                           _Flex_body(N), _Node_id(N), 
                           _Qp(N),_Xp(N),_Zp(N),_Reuler(N) {}


s_marker::~s_marker(void)  {}


Boolean s_marker::Test()
{
   const int err_before=nerr;
   if (_Node_id==Y) {
      // CONTROLLO DELLA SINTASSI IN BASE AL METODO "NODE_ID"
      if (_Zp==Y) out_error(12,"ZP");
      if (_Xp==Y) out_error(12,"XP");
      if (_UseXP==Y) out_error (12,"USEXP");
      if (_Floating==Y) out_error (12,"FLOATING");
      if (_Part==Y) out_error (12,"PART");
      if (_Point_mass==Y) out_error (12,"POINT MASS");
   }
   else {
      // CONTROLLO DELLA SINTASSI IN BASE AL METODO PART O POINT_MASS
      if (_Flex_body==Y) out_error(13,"FLEX_BODY");
      // CONTROLLO DEL TIPO DI DEFINIZIONE:
      // FLOATING:
      if (_Floating) {
	 // CON FLOATING GLI ALTRI PARAMETRI DIVENTANO INUTILI
	 if (_Qp==Y) out_error(14,"QP");
	 if (_Zp==Y) out_error(14,"ZP");
	 if (_Xp==Y) out_error(14,"XP");
	 if (_Reuler==Y) out_error(14,"REULER");
	 if (_UseXP==Y) out_error(14,"USEXP");
      } 
      else {
	 // INDIVIDUAZIONE DEI DUE SISTEMI
	 if (_Reuler==Y) {
	    // SISTEMA DI ANGOLI DI EULERO
	    if (_Zp==Y) out_error(16,"ZP");
	    if (_Xp==Y) out_error(16,"XP");
	    if (_UseXP==Y) out_error(16,"USEXP");
	 } 
      }
   }
   /* PARTE DI DEFINIZIONE DEL MODO */
   /* se non è presente Qp setta la coordinata a 0 *
    * se non è presente Reuler setta gli angoli a 0*
    * se almeno un punto è definito setta il modo a*
    * POINT, invece di default è ANGLES            */
   if (_Qp==N) Qp.Set(0,0,0);
   if (_Reuler==N) Reuler.Set(0,RADIANS,0,RADIANS,0,RADIANS);
   if ((_Zp==Y) | (_Xp==Y)) Mode=s_marker::POINT;
   if (_Reuler==Y) Mode=s_marker::ANGLES;
   if (err_before != nerr) return Y; else return N;   
   /* Se la parte non è stata definita, allora riferisci all'ultima scritta */
   if (_Part==N) { /* qui scrivi la routine che recupera l'ultima parte */ }
}

inline const char* const s_marker::Gettype(void) const
{
   return "MARKER";
}


ostream& s_marker::Print (ostream& out) const
{
   out << endl;
   out << "MARKER:" << label << endl;
   out << "     " << "Part [" << _Part << "] = " << Part << endl;
   out << "     " << "Point mass [" << _Point_mass << "] = " << Point_mass << endl;
   out << "     " << "Floating [" << _Floating << "]" << endl;
   out << "     " << "Flex_Body [" << _Flex_body << "] = " << Flex_body << endl;
   out << "     " << "Node_id [" << _Node_id << "] = " << Node_id << endl;
   out << "     " << "Qp [" << _Qp << "] = ", Qp.Write(out,", ") << endl;
   if (Mode==s_marker::POINT) {
      out << "     " << "Zp [" << _Zp << "] = ", Zp.Write(out,", ") << endl;
      out << "     " << "Xp [" << _Xp << "] = ", Xp.Write(out,", ") << endl;
      out << "     " << "UseXp ? = " << _UseXP << endl;
   }
   if (Mode==s_marker::ANGLES)
     out << "     " << "Reuler [" << _Reuler << "] = " << Reuler << endl;
   return out;
}

void s_marker::Translate (ostream& out)
{
   char* title=new char[80];
   Mat3x3 RelRot,DEBUGROT;
   Vec3 V1,V2,VR1,VR2;
   int IA,IB;
   Id ReferencedTo, MarkerTo;
   Boolean EF=N;
   char* add_remark=new char[80];
   add_remark="";
   // INSERIMENTO NELLA TAVOLA DEI REFERENCE
   if (_Part==Y) {
      // Il marker si riferisce a una parte
      MTP Buffer (Part,MTP::PARTTYPE);
      Marker_Table.insert (MTP_entry(label,Buffer));
      p_PTR_entry pTR;
      pTR=Part_Table.find(Part);
      /* DEBUGGING SECTION */
      if ((*pTR).first==0) {
	 cout << "[DEBUG] : MARKER REFERENCED TO UNEXISTENT PART -> DEBUGGING.." << endl;
	 /* La parte chiamata non esiste -> debugging */
	 s_part* temp = new s_part;
	 temp->label = Part;
	 parts.insert (card_entry(label,temp));
	 temp->Translate(out);
	 pTR=Part_Table.find(Part);
      }
      ReferencedTo = (*pTR).second;
      if (ReferencedTo==0) exit (-3);
   }
   else if (_Point_mass==Y) {
      // Il marker si riferisce a un point-mass
      MTP Buffer (Point_mass,MTP::POINTMASSTYPE);
      Marker_Table.insert (MTP_entry(label,Buffer));
      p_PTR_entry pTR;
      pTR=Pmass_Table.find(Point_mass);
      /* DEBUGGING SECTION */
      if ((*pTR).first==0) {
	 cout << "[DEBUG] : MARKER REFERENCED TO UNEXISTENT POINT MASS -> DEBUGGING.." << endl;
	 /* La parte chiamata non esiste -> debugging */
	 s_pointmass* temp = new s_pointmass;
	 temp->label = Point_mass;
	 pointmasses.insert (card_entry(label,temp));
	 temp->Translate(out);
	 pTR=Pmass_Table.find(Point_mass);
      }      
      ReferencedTo = (*pTR).second;
      if (ReferencedTo==0) exit (-3);
   }
   if (_Floating==Y) {
      /* FLOATING TYPE MARKER */
      /* Questo tipo di marker non è supportabile da MBDYN, per ciò
       * verrà creato un marker riferito alla parte di coordinate pari
       * al BCS della stessa parte */
      add_remark = "\n# FLOATING Marker not supported! - empty fixed marker created";
   }
   switch (Mode) {
    case POINT : {
       V1=(Xp-Qp); V2=(Zp-Qp);
       if (_UseXP==Y) { 
	  IA=1; IB=3;
	  VR1=V1; VR2=V2;
       } else { 
	  IA=3; IB=1; 
          VR1=V2; VR2=V1;
       }
       RelRot = MatR2vec (IA,VR1,IB,VR2);
       /* end of point condition */
    }
      break;
    case ANGLES: {
       /* Trova la matrice di rotazione a partire dagli angoli di Eulero (r) */
       Vec3 T3;
       Euler T = Reuler;
       T.ConvToRadians();
       T3[0]=T.Xangle.value;
       T3[1]=T.Yangle.value;
       T3[2]=T.Zangle.value;
       RelRot = RFromEulerAngles(T3);
       /* End of Euler condition */
    }
      break;
    default :  EF=Y; RelRot=Eye3; break;
   }   
   /* Controllo sulla label : se label è già utilizzato, occorre
    * creare un'altra label libera da utilizzare per lo stesso marker */
   /* ATTENZIONE:FINCHE' E' VALIDO IL QUEUEING DEI MARKER I VALORI SONO = */
   MarkerTo=GetFreeLabel (MBReference,label);
   Reference_Table.insert (MTR_entry(label,MarkerTo));
   /* INSERISCE IL MARKER NELL'APPOSITO DECK DI MBDYN */
   // Velocità ed accellerazione angolare nel sistema MBDyn relativo
   // la label è quella del marker ADAMS, quindi V e Vang=0.
   Vec3 ES1(0,0,0); RVec3 RES1(ES1,MBDyn_entity::REFERENCED,ReferencedTo);
   MBDyn_entity M (MBDyn_entity::REFERENCED,ReferencedTo);
   RVec3 RQp(Qp,M); 
   RMat3x3 RRelRot (RelRot,M);
   MBDyn_reference* R= new MBDyn_reference(MarkerTo,RQp,RRelRot,EF,RES1,RES1);
   char* s = new char[255];
   sprintf (s,"\n# Adams MARKER No. %i become Reference No. %i",
		label,MarkerTo);
   if (add_remark!="") strcat (s,add_remark);
   R->_remark_=s;
   sprintf (title,"MARKER Adams no %i",label);
   R->Title(title);
   
   MBReference.insert(MBDyn_entry(MarkerTo,(MBDyn_card*) R));
   return;
}

/* Only for marker section */

/* La funzione queue alloca le label per i marker con il valore della
 * label ADAMS, in modo che quando viene creato un nuovo reference non
 * viene allocata una label utilizzata da ADAMS per definire i suoi
 * MARKER */

void s_marker::queue (void)
{
   RVec3 RQp,REF,RES1;
   RMat3x3 RRelRot;
   Boolean EF;
   MBDyn_reference* VM = new MBDyn_reference (label,RQp,RRelRot,EF,RES1,RES1);
   MBReference.insert (MBDyn_entry (label,(MBDyn_card*) VM));
   return;
}

void s_marker::dequeue (void)
{
   p_MBDyn_entry i = MBReference.find (label);
   MBReference.erase(i);
   return;
}
