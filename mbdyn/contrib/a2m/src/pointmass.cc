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

//                                POINT MASS                                  

//Definizioni per la card POINT MASS

#include <pointmass.h>

extern MBDyn_deck MBReference;
extern MBDyn_deck MBNodes;
extern MBDyn_deck MBBodies;
extern PTR_deck Pmass_Table;

s_pointmass::s_pointmass (void) :   Mass(0), Cm (0), Vx(0),
                                    Vy(0),Vz(0),Exact(new coord_type[3]),ncoord(0),
                                    Reuler(0,RADIANS,0,RADIANS,0,RADIANS),
                                    _Mass(N),_Cm(N),_Qg(N),_Zg(N),_Xg(N),
                                    _Vx(N), _Vy(N), _Vz(N), _Exact(N),
                                    _Reuler(N)
                                    { }

s_pointmass::~s_pointmass(void)     {}


Boolean s_pointmass::Test()
{
   const int err_before=nerr;
   // CONTROLLO DEL TIPO DI DEFINIZIONE
   if (_Reuler==Y) {
      if (_Zg==Y) out_error (29,"ZG");
      if (_Xg==Y) out_error (29,"XG");
   }
   if (err_before!=nerr) return Y; else return N;
}

inline const char* const s_pointmass::Gettype(void) const
{
   return "POINT MASS";
}

ostream& s_pointmass::Print (ostream& out) const
{
   out << endl;
   out << "POINT MASS:" << label << endl;
   out << "     Cm [" << _Cm << "] = " << Cm << endl
       << "     Mass [" << _Mass << "] = " << Mass << endl
       << "     Qg [" << _Qg << "] = ", Qg.Write(out,", ") << endl
       << "     Zg [" << _Zg << "] = ", Zg.Write(out,", ") << endl
       << "     Xg [" << _Xg << "] = ", Xg.Write(out,", ") << endl
       << "     Vx [" << _Vx << "] = " << Vx << endl
       << "     Vy [" << _Vy << "] = " << Vy << endl
       << "     Vz [" << _Vz << "] = " << Vz << endl
       << "     Reuler [" << _Reuler << "] = " << Reuler << endl
       << "     Exact [" << _Exact << "] = ";
   if (_Exact==Y) outvec(out,Exact,ncoord);
   out << endl << endl;
   return out;
}

void s_pointmass::Translate (ostream& out)
{
   /* Pointmass è equivalente a massa puntiforme */   
   char *comment = new char[80];
   Id Xmarker,JMarker;
   Vec3 Xcm(0,0,0); RVec3 RXcm;
   Vec3 V1,V2,VR1,VR2;
   int IA,IB;
   // Creazione di un reference MBDYN con la giacitura della parte
   Mat3x3 RelRot;
   Boolean _UseXP=Y;
   Boolean EF=N;
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
   Id RefID,PartID,BodyID;
   /* trova un nuovo valore di label per i sistemi di riferimento */
   RefID=GetFreeLabel(MBReference,label);
   /* Inserisce il riferimento : la parte label crea un reference tmplabel */
   Pmass_Table.insert (PTR_entry(label,RefID));
   /* Referenziazione delle coordinate */
   Vec3 V(Vx,Vy,Vz); Vec3 N (0,0,0);
   /* - Velocità - */
   RVec3 RV(V,MBDyn_entity::GLOBAL);
   RVec3 RW(N,MBDyn_entity::GLOBAL);
   /* Baricentro e matrice di rotazione */
   RVec3 RQg (Qg,MBDyn_entity::GLOBAL);
   RMat3x3 RRelRot (RelRot,MBDyn_entity::GLOBAL);
   /* fine referenziazione */
   MBDyn_reference* R = new MBDyn_reference (RefID,RQg,RRelRot,EF,RV,RW);
   MBReference.insert(MBDyn_entry(RefID,(MBDyn_card*) R));
   
   /* 1 rigid body MBDyn con l'inerzia della parte */
   BodyID=GetFreeLabel (MBBodies,label);
   /* Determinazione di centro di massa e riferimento per la Jmatrix */
   if (_Cm==Y) {
      /* ATTENZIONE ! Il marker CM ha la stessa label in MBDyn,
       * in virtù dell'operazione di queueing */
      RXcm=RVec3 (Xcm,MBDyn_entity::REFERENCED,Cm);
   }
     else { 
	RVec3 XCM (Xcm);
	RXcm=XCM;
     }
   Mat3x3 NM;
   RMat3x3 RJMatrix (NM);
   /* Nodo strutturale */
   PartID = GetFreeLabel (MBNodes,label);
   MBDyn_node_structural* CENTRALNODE = new 
     MBDyn_node_structural (PartID,MBDyn_node_structural::DYNAMIC,
			    RQg,RRelRot,RV,RW,0,0,Y);
   MBNodes.insert (MBDyn_entry(PartID,(MBDyn_card*) CENTRALNODE));
   /* Referenziazione di massa e inerzia */
   BodyID = GetFreeLabel (MBBodies,label);
   MBDyn_body* MB = new MBDyn_body (BodyID,PartID,Mass,RXcm,RJMatrix);
   MBBodies.insert (MBDyn_entry(BodyID,(MBDyn_card*) MB));
   return;
}
