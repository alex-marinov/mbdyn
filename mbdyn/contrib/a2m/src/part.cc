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

//                            PART.CC                                     

//Definizioni per la card Part

#include <part.h>
#include <storage.h>

extern MBDyn_deck MBReference;
extern MBDyn_deck MBNodes;
extern MBDyn_deck MBBodies;
extern MBDyn_deck MBJoints;
extern PTR_deck Part_Table;
extern PTR_deck PartToNode_Table;
extern BTP_deck Body_Table;

extern Id GROUND_PART_ID;

s_part::s_part(void)      : Mass(0), Cm(0), Im(0), Material(new char[80]),
                      	   Vm(0),Wm(0),Ground(N), Ip(0,0,0),
                      	   AddIp(0,0,0),Vx(0),Vy(0),Vz(0),Wx(0),Wy(0),Wz(0),
                      	   Reuler(0,RADIANS,0,RADIANS,0,RADIANS),
                      	   Exact(new coord_type[6]),ncoord(0),
                           Mode (s_part::EYEMATRIX),
                           _Mass(N),_Cm(N),_Im(N),_Material(N),_Qg(N),
                           _Zg(N),_Xg(N),_Vm(N),_Wm(N),_Ip(N),
                      	   _AddIp(N),_Vx(N),_Vy(N),_Vz(N),
                      	   _Wx(N),_Wy(N),_Wz(N),_Reuler(N),_Exact(N)
                      	   { Material="STEEL";}

s_part::~s_part (void)     {}

Boolean s_part::Test()
{
   const int err_before=nerr;
   // SE E' DEFINITA LA MASSA DEVE ESSERE DEFINITO CM
   if ((_Mass==Y) & (_Cm==N)) { out_error(10,""); }
   // SE LA MASSA E' NULLA, I PARAMETRI A SEGUITO NON HANNO SENSO
   if (_Mass==N) {
      /* Rem error code is 11 */
      if (_Ip==Y) out_warning(11,"IP");
      if (_Vx==Y) out_warning(11,"VX");
      if (_Vy==Y) out_warning(11,"VY");
      if (_Vz==Y) out_warning(11,"VZ");
      if (_Wx==Y) out_warning(11,"WX");
      if (_Wy==Y) out_warning(11,"WY");
      if (_Wz==Y) out_warning(11,"WZ");
      if (_Vm==Y) out_warning(11,"VM");
      if (_Wm==Y) out_warning(11,"WM");
      if (_Cm==Y) out_warning(11,"CM");
      if (_Im==Y) out_warning(11,"IM");
   }	   
   /* PARTE DI DEFINIZIONE DEL MODO */
   /* se non è presente Qp setta la coordinata a 0 *
    * se non è presente Reuler setta gli angoli a 0*
    * se almeno un punto è definito setta il modo a*
    * POINT, invece di default è ANGLES            */
   if (_Qg==N) Qg.Set(0,0,0);
   if (_Reuler==N) Reuler.Set(0,RADIANS,0,RADIANS,0,RADIANS);
   if ((_Zg==Y) | (_Xg==Y)) Mode=POINT;
   if (_Reuler==Y) Mode=ANGLES;
   /* Verifica l'univocità della parte ground */
   if ((Ground==Y) & (GROUND_PART_ID!=0)) {
     out_error (47,"");
   } else GROUND_PART_ID=label;
   if (err_before != nerr) return Y; else return N;
}

inline const char* const s_part::Gettype(void) const
{
   return "PART";
}

ostream& s_part::Print (ostream& out) const
{
   out << endl;
   out << "PART:" << label;
   if (Ground != Y) { out << endl;
      out << "     " << "Mass [" << _Mass << "] = " << Mass;
      out << "     " << "Cm   [" << _Cm << "] = " << Cm << endl;
      out << "     " << "Im   [" << _Im << "] = " << Im;
      out << "     " << "Material [" << _Material << "] = " << Material << endl;
      out << "     " << "Vm [" << _Vm << "] = " << Vm;
      out << "     " << "Wm [" << _Wm << "] = " << Wm << endl;
      out << "     " << "Qg [" << _Qg << "] = ", Qg.Write(out,", ") << endl;
      out << "     " << "Zg [" << _Zg << "] = ", Zg.Write(out,", ") << endl;
      out << "     " << "Xg [" << _Xg << "] = ", Xg.Write(out,", ") << endl;
      out << "     " << "Ip [" << _Ip << "] = ", Ip.Write(out,", ") << endl;
      out << "     " << "   [" << _AddIp << "] = ", AddIp.Write(out,", ") << endl;
      out << "     " << "Vx Vy Vz [" << _Vx << " " << _Vy << " " << _Vz
       << "] = " << Vx << " " << Vy << " " << Vz << endl;
      out << "     " << "Wx Wy Wz [" << _Wx << " " << _Wy << " " << _Wz
       << "] = " << Wx << " " << Wy << " " << Wz << endl;
      out << "     Reuler [" << _Reuler << "] = " << Reuler;
      out << "     Exact [" << _Exact << "] = ";
      if (_Exact==Y) outvec(out,Exact,ncoord);
      out << endl << endl;
   } else {
      out << " GROUND" << endl;
      out << "     " << "Material [" << _Material << "] = " << Material << endl;
   }
   return out;
}

void s_part::Translate (ostream& out)
{
   char* comment = new char[255];
   char* title = new char[80];
   Boolean _UseXP=Y;
   Id Xmarker,JMarker;
   Id RefID,PartID,BodyID;
   Vec3 Xcm(0,0,0); RVec3 RXcm;
   RMat3x3 RJMatrix;
   int IA,IB;
   Vec3 V1,V2,VR1,VR2;
   // Creazione di un reference MBDYN con la giacitura della parte
   Mat3x3 RelRot;
   Boolean EF=N;
   /* is a GROUND PART ? */
   if (Ground==Y) {
      Id ClampID;
      Mat3x3 DefRot=Eye3;
      RMat3x3 NulRot(DefRot,MBDyn_entity::NUL);
      RVec3 NulPos(Xcm,MBDyn_entity::NUL);
      /* Inserimento del sistema di riferimento "globale" */
      RefID = GetFreeLabel (MBReference,label);
      Part_Table.insert (PTR_entry(label,RefID));
      MBDyn_reference* R = new MBDyn_reference (RefID,NulPos,NulRot,N,
						NulPos,NulPos);
      sprintf(comment,"Reference %i is the GROUND REFERENCED SYSTEM - Ground Part",RefID,label); 
      sprintf(title,"GCS defined by adams PART %i",label);
      R->Remark(comment);
      R->Title(title);
      
      MBReference.insert(MBDyn_entry(RefID,(MBDyn_card*) R));
      MBDyn_entity M (MBDyn_entity::REFERENCED,RefID);
      /* Inserimento del nodo strutturale statico */
      RVec3 RDefPos (Xcm,M);
      RMat3x3 RDefRot (DefRot,M);
      PartID = GetFreeLabel (MBNodes,label);
      PartToNode_Table.insert (PTN_entry(label,PartID));
      MBDyn_node_structural* GROUNDNODE = new
	MBDyn_node_structural (PartID,MBDyn_node_structural::STATIC,
			       RDefPos,RDefRot,RDefPos,RDefPos,0,0,Y);
      MBNodes.insert (MBDyn_entry(PartID,(MBDyn_card*) GROUNDNODE));
      sprintf(comment,"Structural Node %i due to Adams GROUND PART (%i)",
	   PartID,label);
      sprintf(title,"node generated by GROUND part");
      GROUNDNODE->Title(title);
      GROUNDNODE->Remark(comment);
      /* Inserimento del vincolo "clamp" per vincolare a terra il nodo */
      ClampID = GetFreeLabel (MBJoints,label);
      MBDyn_clamp* GROUNDCLAMP = new MBDyn_clamp (ClampID,PartID);
      MBJoints.insert (MBDyn_entry(ClampID,(MBDyn_card*) GROUNDCLAMP));
      sprintf(comment,"Clamp %i is the Clamp for the GROUND PART",
	      ClampID);
      sprintf(title,"clamp generated by GROUND part");
      GROUNDCLAMP->Title(title);
      GROUNDCLAMP->Remark(comment);
      return;
   }
   switch (Mode) {
    case POINT : {
       /* Trova la matrice di rotazione a partire dai vettori sghembi */
       /* RelRot=Plane_ROT (Qp,Zp,Xp,_UseXP); */
       V1=(Xg-Qg); V2=(Zg-Qg);
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
   /* trova un nuovo valore di label per i sistemi di riferimento */
   RefID=GetFreeLabel(MBReference,label);
   /* Inserisce il riferimento : la parte label crea un reference tmplabel */
   Part_Table.insert (PTR_entry(label,RefID));
   Vec3 V(Vx,Vy,Vz); Vec3 W(Wx,Wy,Wz);
   /* Referenziazione delle coordinate */
   /* - Velocità - */
   MBDyn_entity M (MBDyn_entity::NUL);
   if (_Vm!=N) M = MBDyn_entity (MBDyn_entity::REFERENCED,Vm);
   RVec3 RV(V,M);
   if (_Wm!=N) M = MBDyn_entity (MBDyn_entity::REFERENCED,Wm);
   RVec3 RW(W,M);
   /* Baricentro e matrice di rotazione */
   /* Grandezze nel sistema globale, uso NUL altrimenti dovrei usare
    * MBDyn_entity::GLOBAL */
   RVec3 RQg (Qg,MBDyn_entity::NUL);
   RMat3x3 RRelRot (RelRot,MBDyn_entity::NUL);
   /* fine referenziazione -> reference MBDYN nel sistema nodo */
   MBDyn_reference* R = new MBDyn_reference (RefID,RQg,RRelRot,EF,RV,RW);
   MBReference.insert(MBDyn_entry(RefID,(MBDyn_card*) R));

   sprintf(comment,"Reference %i is BCS for Adams PART %i",RefID,label); 
   R->Remark(comment);
   sprintf(title,"BCS for Adams PART %i",label);
   R->Title(title);
   
   MBDyn_entity MII (MBDyn_entity::REFERENCED,RefID);
   RVec3 DRQg (Xcm,MII);
   RMat3x3 DRRelRot (RelRot,MII);
   RV.REF = MII;
   RW.REF = MII;
   // 1 structural node MBDYN attaccato al reference
   /* Recupera una label libera per i nodi */
   PartID=GetFreeLabel(MBNodes,label);
   /* Inserisce il riferimento-> Parte I uguale a nodo strutturale PartId */
   PartToNode_Table.insert (PTN_entry(label,PartID));
   /* Definizione del nodo da inserire */
   MBDyn_node_structural* CENTRALNODE = new
     MBDyn_node_structural (PartID,MBDyn_node_structural::DYNAMIC,
			    DRQg,DRRelRot,
			    RV,RW,0,0,Y);   
   MBNodes.insert(MBDyn_entry(PartID,(MBDyn_card*) CENTRALNODE));

   sprintf(comment,"Structural Node %i due to Adams PART %i",
	   PartID,label);
   CENTRALNODE->Remark(comment);
   sprintf(title,"node generated by adams PART %i",label);
   CENTRALNODE->Title(title);
   
   /* La referenziazione di centro di massa e matrice di inerzia avviene 
    utilizzando l'informazione proveniente da ADAMS ( CM e IM ) */
   
   /* 1 rigid body MBDyn con l'inerzia della parte */
   BodyID=GetFreeLabel (MBBodies,label);
   Mat3x3 JMatrix=MomentToInertia(Ip,AddIp);
   /* Determinazione di centro di massa e riferimento per la Jmatrix */
   if (_Cm==Y) {
      
      /* ATTENZIONE ! Il marker CM ha la stessa label in MBDyn,
       * in virtù dell'operazione di queueing definita prima */
     
      RXcm=RVec3 (Xcm,MBDyn_entity::REFERENCED,Cm);
      RMat3x3 RJM (JMatrix,MBDyn_entity::REFERENCED,Cm);
      RJMatrix=RJM;
   }
     else { 
	RVec3 XCM (Xcm);
	RXcm=XCM;
     }
   if (_Im==Y) {
      
      /* ATTENZIONE ! Il marker IM ha la stessa label in MBDyn,
       * in virtù dell'operazione di queueing definita prima */
      
      RMat3x3 RJM (JMatrix,MBDyn_entity::REFERENCED,Im);
      RJMatrix=RJM;
   }
   MBDyn_body* MB= new MBDyn_body (BodyID,PartID,Mass,RXcm,RJMatrix);
   MBBodies.insert (MBDyn_entry(BodyID,(MBDyn_card*) MB));

   /* A quale nodo strutturale corrisponde il body BodyID? .>inserisci */
   Body_Table.insert (BTP_entry(BodyID,PartID));
   
   sprintf(title,"body generated by Adams PART %i",label);
   MB->Title(title);
   sprintf(comment,"Rigid Body %i due to Adams PART %i ",BodyID,label);
   MB->Remark(comment);

   return;
}
