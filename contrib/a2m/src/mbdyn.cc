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

/* MBDYN.CC - SOURCE CODE FOR MBDYN ELEMENTS AND ROUTINES */

#include <mbdyn.h>
#include <map.h>

extern Boolean REMOVE_REMARK;
char ph='"';

// GENERICA CARTA DI MBDYN, possiede una label e un tipo

MBDyn_card::MBDyn_card (Id L, MBDyn_card::Type q) :
                        Label(L),Card_type(q),_remark_(NULL),
                        _title_(NULL)
{}

MBDyn_card::~MBDyn_card ()
{}

MBDyn_card::MBDyn_card () : Label(0),Card_type(_LAST_CARD),
                            _remark_(NULL),_title_(NULL)
{}

ostream& MBDyn_card::Test(ostream& out) const
{
   out << "<- card" << endl;
   return out;
}

void MBDyn_card::Remark (char* R)
{
   _remark_=new char[strlen(R)+4];
   _remark_[0]='\n';
   _remark_[1]='#';
   _remark_[2]=' ';
   for (int i=0;i<strlen(R);i++) _remark_[i+3]=' ';
   for (int i=0;i<strlen(R);i++) {
      _remark_[i+3]=R[i];
   }
   _remark_[strlen(R)+3]='\0';
   return;
}

void MBDyn_card::Title (char* R)
{
   _title_=new char[strlen(R)+3];
   _title_[0]=ph;
   for (int i=0;i<strlen(R);i++) _title_[i+1]=' ';
   for (int i=0;i<strlen(R);i++) {
      _title_[i+1]=R[i];
   }
   _title_[strlen(R)+1]=ph;
   _title_[strlen(R)+2]='\0';
   return;
}

// GENERICO ELEMENTO DI MBDYN, possiede una label e un tipo elemento

MBDyn_elem::MBDyn_elem (Id L,MBDyn_elem::Type q) :  
                        MBDyn_card (L,MBDyn_card::_ELEMENT), 
                        Elem_type(q)
{}

MBDyn_elem::~MBDyn_elem ()
{}

MBDyn_elem::MBDyn_elem () : MBDyn_card (0,_LAST_CARD),
                            Elem_type(MBDyn_elem::UNKNOWN)
{}

ostream& MBDyn_elem::Restart (ostream& out) const
{
   return out;
}

ostream& MBDyn_elem::Test(ostream& out) const
{
   out << "<- Element -<";
   MBDyn_card::Test(out);
   return out;
}

inline const char* const MBDyn_elem::Gettype(void) const
{
   return "ELEMENT";
}

// LEGGE COSTITUTIVA PER TRAVI

ConstitutiveLaw6D::ConstitutiveLaw6D (double Eyoung=0, double Gtang=0):
				     Matrix(6,6), E(Eyoung), G(Gtang)
{
   /* Legge costitutiva del materiale di tipo elastico isotropico */
   double v=(E/(2*G))-1;
   double c=E/((1+v)*(1-2*v));
   for (int i=0; i<3; i++) (*data[i])[i]=(1-v)*c;
   for (int i=3; i<6; i++) (*data[i])[i]=G;
   (*data[0])[1]=v*c;
   (*data[0])[2]=v*c;
   (*data[1])[2]=v*c;
}

ConstitutiveLaw6D::~ConstitutiveLaw6D ()
{
}

ostream& ConstitutiveLaw6D::Restart(ostream& out, const char* indent) const
{
   out << "linear elastic generic," << endl;
   out << "",Write (out,", ",",\n",indent);
   return out;
}

// BEAM

// costruttore
MBDyn_beam::MBDyn_beam ()
{
   Node[0]=0; Node[1]=0; Node[2]=0;
   Label=0;
   pD[0]=new Mat6x6(),pD[1]=new Mat6x6();
}

MBDyn_beam::MBDyn_beam	(Id L, Id N1, Id N2, Id N3, Vec3 F1, Vec3 F2, 
			 Vec3 F3,Mat3x3 R1, Mat3x3 R2,Mat6x6* CS1,
			 Mat6x6* CS2) :
                         MBDyn_elem (L,_BEAM)
{
   Node[0]=N1; Node[1]=N2; Node[2]=N3;
   R1=Eye3; R2=Eye3;
   f[0]=F1; f[1]=F2; f[2]=F3;
   R[0]=R1; R[1]=R2;
   pD[0]=CS1; pD[1]=CS2;
}
		
// Distruttore
MBDyn_beam::~MBDyn_beam () {}

// funzione di Restart
ostream& MBDyn_beam::Restart (ostream& out) const
{
   const char indent[]="         ";
   out << "  beam:  " << Label;
   if (_title_!=NULL) out << ", name, " << _title_;
   for (unsigned int i=0; i <NUMNODES; i++) {
      out << ", " << endl << indent << Node[i] << ", reference, node, ",
      f[i].Write(out,", ");
   }
   for (unsigned int i=0; i< NUMSEZ; i++) {
      out << ", " << endl << indent << "reference, global, "
	<< endl << indent, R[i].RWrite(out,", ","") << 
	", " << endl
	<< indent << "linear elastic generic," << endl
	<< "",pD[i]->Write(out,",",",\n",indent);
   }
   out << ";" << endl;
   return out;
}

// Restituzione del tipo
inline const char* const MBDyn_beam::Gettype(void) const
{
   return "BEAM";
}

ostream& MBDyn_beam::Test(ostream& out) const
{
   out << "<- beam -<";
   MBDyn_elem::Test(out);
   return out;
}

// RIGID BODY

// Costruttore
MBDyn_body::MBDyn_body()
{
}

MBDyn_body::MBDyn_body (Id L,Id NL, double M, RVec3 Xg, RMat3x3 J0) :
                        Node(NL),
                        MBDyn_elem(L,MBDyn_elem::_RIGIDBODY),Mass(M), 
                        Xgc(Xg), Jota(J0), Inertial(N), defNodeInertial(N)
{}

MBDyn_body::MBDyn_body (Id L,Id NL, double M, RVec3 Xg, RMat3x3 J0, Id L2) :
                        Node(NL),
                        Mass(M),Xgc(Xg),Jota(J0),Inertial(Y),
                        NodeInertial(L2),defNodeInertial(Y),
                        MBDyn_elem(L,MBDyn_elem::_RIGIDBODY)
{}

MBDyn_body::MBDyn_body (Id L, Id NL, double M, RVec3 Xg, RMat3x3 J0, 
			Mat3x3 R) :
                        Node(NL),
                        Mass(M),Xgc(Xg),Jota(J0),Inertial(Y),
                        RotInertial(R),defNodeInertial(N),
                        MBDyn_elem(L,MBDyn_elem::_RIGIDBODY)
{}

// Distruttore
MBDyn_body::~MBDyn_body ()
{}

// Funzione di Restart
ostream& MBDyn_body::Restart(ostream& out) const
{
   const char* indent = "         ";
   out << "  body:  " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   out << endl << indent
     << Node << ", " << Mass << ", " << endl << indent,
     Xgc.Write(out,", ") << ", " << endl
     << indent,Jota.Write(out,", ",", ","");
     if (Inertial==Y) {
	out << endl << indent << ", inertial, ";
	if (defNodeInertial==Y)
	    out << NodeInertial;
	else
	    RotInertial.Write(out,", ",", ");
     }
     out << " ;" << endl;
   return out;
}

// Restituzione del tipo
const char* const MBDyn_body::Gettype(void) const
{
   return "MBDYN RIGID BODY";
}

ostream& MBDyn_body::Test(ostream& out) const
{
   out << "<- Rigid body -<";
   MBDyn_elem::Test(out);
   return out;
}

// MBDYN FORCES
MBDyn_force::MBDyn_force (Id L,MBDyn_force::Type q) : 
                          Force_type(q),
                          MBDyn_elem(L,_FORCE) {}

MBDyn_force::MBDyn_force (Id L,MBDyn_force::Type q, MBDyn_force::FType qm) :
                          Force_type(q),Mode_type(qm),
                          MBDyn_elem(L,_FORCE) {}

MBDyn_force::MBDyn_force () {}
MBDyn_force::~MBDyn_force () {}

ostream& MBDyn_force::Restart(ostream& out) const
{
   out << "  force:  " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   return out;
}

inline const char* const MBDyn_force::Gettype(void) const
{
   return "FORCE";
}

inline const char* const MBDyn_force::GetMode(void) const
{
   if (Mode_type==MBDyn_force::CONSERVATIVE)
     return "conservative";
   else if (Mode_type==MBDyn_force::FOLLOWER)
     return "follower";
}

ostream& MBDyn_force::Test(ostream& out) const
{
   out << "<- force -<";
   MBDyn_elem::Test(out);
   return out;
}

// - ABSTRACT FORCE
MBDyn_force_abstract::MBDyn_force_abstract (Id L, Id dof, 
					    MBDyn_drive_caller *M) :
                      Nodedof(dof),Magnitude(M),MBDyn_force(L,ABSTRACT)
{}
MBDyn_force_abstract::~MBDyn_force_abstract ()
{}

ostream& MBDyn_force_abstract::Restart(ostream& out) const
{
   MBDyn_force::Restart(out) << " abstract, " << Nodedof << ", ";
   Magnitude->Restart(out);
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_force_abstract::Gettype(void) const
{
   return "ABSTRACT FORCE";
}

ostream& MBDyn_force_abstract::Test(ostream& out) const
{
   out << "<- abstract force -<";
   MBDyn_force::Test(out);
   return out;
}

// - STRUCTURAL FORCE
MBDyn_force_structural::MBDyn_force_structural (Id L,MBDyn_force::FType q,
			Id N, RVec3 reldir,Vec3 relarm, 
		        MBDyn_drive_caller* magn) :
                        Node(N),Relative_direction(reldir),
                        Relative_arm(relarm),Magnitude(magn),
                        MBDyn_force(L,STRUCTURALFORCE,q)
{}

MBDyn_force_structural::~MBDyn_force_structural ()
{}

ostream& MBDyn_force_structural::Restart(ostream& out) const
{
   const char* indent = "          ";
   MBDyn_force::Restart(out) << MBDyn_force::GetMode() << ", "
     << endl << indent
     << Node << ", ",Relative_direction.Write(out,", ")
     << ", " << endl << indent,Relative_arm.Write(out,", ") << ", "
     << endl << indent;
     Magnitude->Restart(out);
     out << ";" << endl;
   return out;
}

inline const char* const MBDyn_force_structural::Gettype(void) const
{
   return "STRUCTURAL FORCE";
}

ostream& MBDyn_force_structural::Test(ostream& out) const
{
   out << "<- structural force -<";
   MBDyn_force::Test(out);
   return out;
}

// STRUCTURAL COUPLE
MBDyn_force_couple::MBDyn_force_couple (Id L,MBDyn_force::FType q,
					Id N, RVec3 reldir, 
					MBDyn_drive_caller* magn) :
                                        Node(N),
                                        Relative_direction(reldir), 
                                        Magnitude(magn),
                                        MBDyn_force(L,STRUCTURALCOUPLE,q)
{}

MBDyn_force_couple::~MBDyn_force_couple ()
{}

ostream& MBDyn_force_couple::Restart (ostream& out) const
{
   const char* indent= "           ";
   out << "  couple:  " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   out << MBDyn_force::GetMode() << ", " 
     << endl << indent << Node 
     << ", ",Relative_direction.Write(out,", ")
     << ", " << endl << indent; 
   Magnitude->Restart(out);
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_force_couple::Gettype(void) const
{
   return "STRUCTURAL COUPLE";
}

ostream& MBDyn_force_couple::Test(ostream& out) const
{
   out << "<- structural couple -<";
   MBDyn_force::Test(out);
   return out;
}

// REFERENCE

// Costruttore N.1
MBDyn_reference::MBDyn_reference(Id L,RVec3 P, RMat3x3 R1, 
				 Boolean flag, RVec3 V, RVec3 O) :
				 Abs_Pos(P),Eye_Flag(flag),
				 Abs_Vel(V),Abs_Omega(O),
                                 Abs_Rot_Matrix(R1),
                                 MBDyn_card(L,_REFERENCE)
{}

// Costruttore N.2
MBDyn_reference::~MBDyn_reference()  {}

// Costruttore N.3
MBDyn_reference::MBDyn_reference(Id L) : MBDyn_card (L,_REFERENCE) {}


ostream& MBDyn_reference::Restart(ostream& out) const
{
   const char* indent="              ";
   out << "  reference:  " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", " << endl << indent;
   out << "",Abs_Pos.Write(out,", ") << ", " << endl,
     Abs_Rot_Matrix.RWrite(out,", ","\n",indent);
   out << ", " << endl << indent,
     Abs_Vel.Write(out,", ") << ", " << endl << indent,
     Abs_Omega.Write(out,", ") 
       << ";" << endl;
   return out;
}

inline const char* const MBDyn_reference::Gettype(void) const
{
   return "Mbdyn reference system";
}

ostream& MBDyn_reference::Test(ostream& out) const
{
   out << "<- reference -<";
   return out;
}

// JOINT  

MBDyn_joint::MBDyn_joint (Id N, MBDyn_joint::Type q) : 
                          Joint_type(q), MBDyn_elem (N,_JOINT) {}
MBDyn_joint::MBDyn_joint () : Joint_type(UNKNOWN), MBDyn_elem (0,_JOINT) {}

MBDyn_joint::~MBDyn_joint () {}

ostream& MBDyn_joint::Restart (ostream& out) const
{
   out << "  joint:  " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   return out;
}

inline const char* const MBDyn_joint::Gettype(void) const
{
   return "";
}

ostream& MBDyn_joint::Test(ostream& out) const
{
   out << "<- joint -<";
   MBDyn_elem::Test(out);
   return out;
}


// SPHERICAL JOINT

MBDyn_spherhinge::MBDyn_spherhinge (Id N, Id N1,Id N2,Vec3 D1,
					Mat3x3 R1,Vec3 D2,Mat3x3 R2) :
                                        MBDyn_joint (N,SPHERICALHINGE)
{
   Node[0]=N1; Node[1]=N2;
   d[0]=D1; d[1]=D2;
   Rh[0]=R1; Rh[1]=R2;
}

MBDyn_spherhinge::~MBDyn_spherhinge ()
{}

ostream& MBDyn_spherhinge::Restart(ostream& out) const
{
   const char *indent = "          ";
   Vec3 Null(0,0,0);
   Mat3x3 E=Eye3;
   MBDyn_joint::Restart(out);
   out << "spherical hinge, ";
   for (int i = 0; i<2; i++) {
      out << endl << indent << Node[i];
      if (d[i]!=Null) out << ", reference, node, ", 
	d[i].Write(out,", ");
      if (Rh[i]!=E) out << ", " << endl 
	<< indent << "hinge, reference, node, ", Rh[i].RWrite(out,", ","");
      if (i==0) out << ", ";
   }
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_spherhinge::Gettype(void) const
{
   return "SPHERICAL HINGE";
}

ostream& MBDyn_spherhinge::Test(ostream& out) const
{
   out << "<- spherical hinge -<";
   MBDyn_joint::Test(out);
   return out;
}

// PIN

MBDyn_pin::MBDyn_pin (Id N,Id ND,Vec3 R, Vec3 A) :
                                  Node(ND),Relative_offset(R),
                                  Absolute_pos(A),
                                  MBDyn_joint (N,PIN)
{}

MBDyn_pin::~MBDyn_pin ()
{}

ostream& MBDyn_pin::Restart(ostream& out) const
{
   MBDyn_joint::Restart(out);
   out << " pin, " << Node << ", reference, node, ",
       Relative_offset.Write(out,", ")
       << ", reference, global, ",Absolute_pos.Write(out,", ")
       << ";" << endl;
   return out;
}
   
inline const char* const MBDyn_pin::Gettype(void) const
{
   return "PIN";
}

ostream& MBDyn_pin::Test(ostream& out) const
{
   out << "<- pin -<";
   MBDyn_joint::Test(out);
   return out;
}

// CLAMP

MBDyn_clamp::MBDyn_clamp (Id L,Id ND) :
                   Node(ND),NodeRef1(Y),NodeRef2(Y),
                   MBDyn_joint (L,CLAMP)
{}

MBDyn_clamp::MBDyn_clamp (Id L,Id ND,Mat3x3 ARM) :
                   Node(ND),NodeRef1(Y),Absolute_rot_matrix(ARM),
                   MBDyn_joint (L,CLAMP),NodeRef2(N)
{}

MBDyn_clamp::MBDyn_clamp (Id L,Id ND,Vec3 AP) :
                   Node(ND),Absolute_pos(AP),NodeRef2(Y),
                   MBDyn_joint (L,CLAMP),NodeRef1(N)
{}

MBDyn_clamp::MBDyn_clamp (Id L,Id ND,Vec3 AP,Mat3x3 ARM) :
                   Node(ND),Absolute_pos(AP),Absolute_rot_matrix(ARM),
                   MBDyn_joint (L,CLAMP),NodeRef1(N),NodeRef2(N)
{}

MBDyn_clamp::~MBDyn_clamp ()
{}


ostream& MBDyn_clamp::Restart(ostream& out) const
{
   MBDyn_joint::Restart(out);
   out << " clamp, " << Node;
   if (NodeRef1==Y) out << ", node";
   else {
      out << ", reference, node, ",Absolute_pos.Write(out,", ");
   }
   if (NodeRef2==Y) out << ", node";
   else {
      out << ", reference, node, 1, ",
      (Absolute_rot_matrix.GetVec(1)).Write(out,", ") << ", 2, ",
      (Absolute_rot_matrix.GetVec(2)).Write(out,", ");
   }
   out << ";" << endl;
   return out;
}


inline const char* const MBDyn_clamp::Gettype(void) const
{
   return "CLAMP";
}

ostream& MBDyn_clamp::Test(ostream& out) const
{
   out << "<- clamp -<";
   MBDyn_joint::Test(out);
   return out;
}

// UNIVERSAL HINGE

MBDyn_univhinge::MBDyn_univhinge (Id N, Id N1,Id N2,Vec3 D1,
					Mat3x3 R1,Vec3 D2,Mat3x3 R2) :
                                        MBDyn_joint (N,UNIVERSALHINGE)
{
   Node[0]=N1; Node[1]=N2;
   d[0]=D1; d[1]=D2;
   Rh[0]=R1; Rh[1]=R2;
}

MBDyn_univhinge::~MBDyn_univhinge ()
{}

ostream& MBDyn_univhinge::Restart(ostream& out) const
{
   MBDyn_joint::Restart(out);
   out << " universal hinge ";
   for (int i = 0; i<2; i++) {
      out << ", " << Node[i] << ", reference, node, ",
       d[i].Write(out,", ") << ", hinge, reference, node, 1, ",
       (Rh[i].GetVec(1)).Write(out,", ")
       << ", 2, ",
       (Rh[i].GetVec(2)).Write(out,", ");
   }
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_univhinge::Gettype(void) const
{
   return "UNIVERSAL HINGE";
}

ostream& MBDyn_univhinge::Test(ostream& out) const
{
   out << "<- universal hinge -<";
   MBDyn_joint::Test(out);
   return out;
}

// UNIVERSAL PIN

MBDyn_univpin::MBDyn_univpin (Id L,Id NR, Vec3 RO, Mat3x3 RRM,
			      Vec3 AP, Mat3x3 APR) :
                              MBDyn_joint (L,UNIVERSALPIN),
                              NodeRef(NR), Relative_offset(RO),
                              Relative_rot_matrix(RRM),
                              Absolute_pos (AP),
                              Hinge1(Y),Hinge2(Y),
                              Absolute_pin_rot (APR)
{}

MBDyn_univpin::MBDyn_univpin (Id L, Id NR, Vec3 RO, Vec3 AP,
			      Mat3x3 APR) :
                              MBDyn_joint (L,UNIVERSALPIN),
                              NodeRef(NR), Relative_offset(RO),
                              Absolute_pos(AP),
                              Hinge1(N),Hinge2(Y),
                              Absolute_pin_rot (APR)
{}

MBDyn_univpin::MBDyn_univpin (Id L, Id NR, Vec3 RO, Mat3x3 RRM,
			      Vec3 AP) :
                              MBDyn_joint (L,UNIVERSALPIN),
                              NodeRef(NR), Relative_offset (RO),
                              Relative_rot_matrix (RRM),
                              Hinge1(Y),Hinge2(N),
                              Absolute_pos (AP)
{}

MBDyn_univpin::MBDyn_univpin (Id L, Id NR, Vec3 RO, Vec3 AP) :
                              MBDyn_joint (L,UNIVERSALPIN),
                              NodeRef(NR),Relative_offset(RO),
                              Hinge1(N),Hinge2(N),
                              Absolute_pos (AP)
{}

MBDyn_univpin::~MBDyn_univpin ()
{}

ostream& MBDyn_univpin::Restart (ostream& out) const
{
   MBDyn_joint::Restart(out);
   out << " universal pin, " << NodeRef << ", reference, node, ",
   Relative_offset.Write(out,", ");
   if (Hinge1==Y) {
      out << ", hinge, reference, node, 1, ", 
      (Relative_rot_matrix.GetVec(1)).Write(out,", ") << ", 2, ",
      (Relative_rot_matrix.GetVec(2)).Write(out,", ");
   }
   out << ", reference, global, ", Absolute_pos.Write(out,", ");
   if (Hinge2==Y) {
      out << ", hinge, reference, node, 1, ",
      (Absolute_pin_rot.GetVec(1)).Write(out,", ") << ", 2, ",
      (Absolute_pin_rot.GetVec(2)).Write(out,", ");
   }
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_univpin::Gettype(void) const
{
   return "UNIVERSAL PIN";
}

ostream& MBDyn_univpin::Test(ostream& out) const
{
   out << "<- universal pin -<";
   MBDyn_joint::Test(out);
   return out;
}

// PLANE PIN

MBDyn_planepin::MBDyn_planepin (Id L,Id NR, Vec3 RO, Mat3x3 RRM,
			      Vec3 AP, Mat3x3 APR) :
                              MBDyn_joint (L,PLANEPIN),
                              NodeRef(NR), Relative_offset(RO),
                              Relative_rot_matrix(RRM),
                              Absolute_pos (AP),
                              Hinge1(Y),Hinge2(Y),
                              Absolute_pin_rot (APR)
{}

MBDyn_planepin::MBDyn_planepin (Id L, Id NR, Vec3 RO, Vec3 AP,
			      Mat3x3 APR) :
                              MBDyn_joint (L,PLANEPIN),
                              NodeRef(NR), Relative_offset(RO),
                              Absolute_pos(AP),
                              Hinge1(N),Hinge2(Y),
                              Absolute_pin_rot (APR)
{}

MBDyn_planepin::MBDyn_planepin (Id L, Id NR, Vec3 RO, Mat3x3 RRM,
			      Vec3 AP) :
                              MBDyn_joint (L,PLANEPIN),
                              NodeRef(NR), Relative_offset (RO),
                              Relative_rot_matrix (RRM),
                              Hinge1(Y),Hinge2(N),
                              Absolute_pos (AP)
{}

MBDyn_planepin::MBDyn_planepin (Id L, Id NR, Vec3 RO, Vec3 AP) :
                              MBDyn_joint (L,PLANEPIN),
                              NodeRef(NR),Relative_offset(RO),
                              Hinge1(N),Hinge2(N),
                              Absolute_pos (AP)
{}

MBDyn_planepin::~MBDyn_planepin ()
{}

ostream& MBDyn_planepin::Restart (ostream& out) const
{
   MBDyn_joint::Restart(out);
   out << " plane pin, " << NodeRef << ", reference, node, ",
   Relative_offset.Write(out,", ");
   if (Hinge1==Y) {
      out << ", hinge, reference, node, 1, ", 
      (Relative_rot_matrix.GetVec(1)).Write(out,", ") << ", 2, ",
      (Relative_rot_matrix.GetVec(2)).Write(out,", ");
   }
   out << ", reference, global, ", Absolute_pos.Write(out,", ");
   if (Hinge2==Y) {
      out << ", hinge, reference, node, 1, ",
      (Absolute_pin_rot.GetVec(1)).Write(out,", ") << ", 2, ",
      (Absolute_pin_rot.GetVec(2)).Write(out,", ");
   }
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_planepin::Gettype(void) const
{
   return "PLANE PIN";
}

ostream& MBDyn_planepin::Test(ostream& out) const
{
   out << "<- plane pin -<";
   MBDyn_joint::Test(out);
   return out;
}

// PLANE HINGE

MBDyn_planehinge::MBDyn_planehinge (Id N, Id N1,Id N2,Vec3 D1,
					Mat3x3 R1,Vec3 D2,Mat3x3 R2) :
                                        MBDyn_joint (N,PLANEHINGE)
{
   Node[0]=N1; Node[1]=N2;
   d[0]=D1; d[1]=D2;
   Rh[0]=R1; Rh[1]=R2;
}

MBDyn_planehinge::~MBDyn_planehinge ()
{}

ostream& MBDyn_planehinge::Restart(ostream& out) const
{
   const char *indent = "          ";
   Vec3 Null(0,0,0);
   Mat3x3 E = Eye3;
   MBDyn_joint::Restart(out);
   out << "plane hinge, ";
   for (int i = 0; i<2; i++) {
      out << endl << indent << Node[i];
      out << ", reference, node, ", d[i].Write(out,", ");
      if (Rh[i]!=E) out << ", " << endl 
	<< indent << "hinge, reference, node, ", Rh[i].RWrite(out,", ","");
      if (i==0) out << ", ";
   }
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_planehinge::Gettype(void) const
{
   return "PLANE HINGE";
}

ostream& MBDyn_planehinge::Test(ostream& out) const
{
   out << "<- plane hinge -<";
   MBDyn_joint::Test(out);
   return out;
}

// IN PLANE JOINT

MBDyn_inplane::MBDyn_inplane (Id N, Id N1, Id N2, Vec3 RP,
			      Vec3 RND, Vec3 RO) :
                              MBDyn_joint (N,INPLANESIMPLE)
{
   Node[0]=N1; Node[1]=N2;
   Relative_position=RP;
   Relative_norm_direction=RND;
   Relative_offset=RO;
}

MBDyn_inplane::MBDyn_inplane (Id N, Id N1, Id N2, Vec3 RP,
			      Vec3 RND) : 
                              MBDyn_joint (N, INPLANESIMPLE)
{
   Vec3 Null(0,0,0);
   Node[0]=N1; Node[1]=N2;
   Relative_position = RP;
   Relative_norm_direction=RND;
   Relative_offset = Null;
}

MBDyn_inplane::~MBDyn_inplane ()
{}

ostream& MBDyn_inplane::Restart (ostream& out) const
{
   Vec3 Null (0,0,0);
   const char *indent = "          ";
   MBDyn_joint::Restart (out);
   out << "in plane, "
     << endl << indent << Node[1] << ", "
     << endl << indent,Relative_position.Write(out,",") << ", "
     << endl << indent << Node[2];
     if (Relative_offset!=Null) out << endl << indent
     << "offset, ", Relative_offset.Write(out,", ");
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_inplane::Gettype (void) const
{
   return "IN PLANE";
}

ostream& MBDyn_inplane::Test (ostream& out) const
{
   out << "<- in plane -<";
   MBDyn_joint::Test(out);
   return out;
}


// INLINE JOINT 

MBDyn_inline::MBDyn_inline (Id N, Id N1, Id N2, Vec3 RLP,
			    Mat3x3 ROR, Vec3 ROS) :
                            MBDyn_joint (N, IN_LINE)
{
   Node[0]=N1; Node[1]=N2;
   Relative_line_position=RLP;
   Relative_orientation=ROR;
   Relative_offset = ROS;
}

MBDyn_inline::MBDyn_inline (Id N, Id N1, Id N2, Vec3 RLP,
			    Mat3x3 RO) : 
                            MBDyn_joint (N, IN_LINE)
{
   Vec3 Null(0,0,0);
   Node[0]=N1; Node[1]=N2;
   Relative_line_position=RLP;
   Relative_orientation=RO;
   Relative_offset = Null;
}

MBDyn_inline::~MBDyn_inline ()
{}

ostream& MBDyn_inline::Restart (ostream& out) const
{
   Vec3 Null(0,0,0);
   const char* indent = "          ";
   MBDyn_joint::Restart(out);
   out << "in line, "
     << endl << indent << Node[1] << ", "
     << endl << indent, Relative_line_position.Write(out,", ")
       << endl << indent, Relative_orientation.RWrite(out,", "," ,")
	 << endl << indent << Node[2];
   if (Relative_offset!=Null) out << ", " << endl << indent
     << "offset, ",Relative_offset.Write(out,", ");
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_inline::Gettype (void) const
{
   return "IN LINE";
}

ostream& MBDyn_inline::Test (ostream& out) const
{
   out << "<- in line -<";
   MBDyn_joint::Test (out);
   return out;
}
  
// MBDYN BULK ELEMENT


MBDyn_bulk::MBDyn_bulk (Id N,MBDyn_bulk::Type q) : 
                        Bulk_type(q),MBDyn_elem(N,_BULK)
{}

MBDyn_bulk::~MBDyn_bulk ()
{}

MBDyn_StiffnessSpring::MBDyn_StiffnessSpring (Id N, Id ND, double S) :
                       MBDyn_bulk(N,STIFFNESS_SPRING),
                       Stiffness(S),NodeDof(ND)
{}

MBDyn_StiffnessSpring::~MBDyn_StiffnessSpring ()
{}

inline const char* const MBDyn_bulk::Gettype(void) const
{
   return "BULK";
}

inline const char* const MBDyn_StiffnessSpring::Gettype(void) const
{
   return "STIFFNESS SPRING";
}

ostream& MBDyn_bulk::Restart(ostream& out) const
{
   out << " bulk : " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   return out;
}

ostream& MBDyn_bulk::Test(ostream& out) const
{
   out << "<- bulk -<";
   MBDyn_elem::Test(out);
   return out;
}

ostream& MBDyn_StiffnessSpring::Restart (ostream& out) const
{
   MBDyn_bulk::Restart(out);
   out << " stiffness spring, " << NodeDof << ", " << Stiffness 
       << ";" << endl;
   return out;
}

ostream& MBDyn_StiffnessSpring::Test (ostream& out) const
{
   out << "<- stiffness spring -<";
   MBDyn_bulk::Test(out);
   return out;
}

// NODE 

MBDyn_node::MBDyn_node (Id N, MBDyn_node::Type q) : 
                        Node_type(q),MBDyn_card (N,MBDyn_card::_NODE) 
{}

MBDyn_node::MBDyn_node () : Node_type(UNKNOWN),
                            MBDyn_card (0,MBDyn_card::_NODE)
{}

MBDyn_node::~MBDyn_node () {}


ostream& MBDyn_node::Restart (ostream& out) const
{
   return out;
}

inline const char* const MBDyn_node::Gettype (void) const
{
   return "";
}

ostream& MBDyn_node::Test(ostream& out) const
{
   out << "<- node -<";
   MBDyn_card::Test(out);
   return out;
}


// STRUCTURAL NODE

// Structural node - STATIC,DYNAMIC definition
MBDyn_node_structural::MBDyn_node_structural (Id N,
		       MBDyn_node_structural::Type q,
		       RVec3 AP, RMat3x3 ARM, 
		       RVec3 AV,RVec3 AAV, double PIS, 
		       double VIS, Boolean O) :
        Mode_type(q),Abs_Pos (AP), Abs_Rot_Matrix (ARM), Abs_Vel(AV), 
        Abs_Ang_Vel(AAV),position_initial_stiffness(PIS), 
        velocity_initial_stiffness(VIS), Omega_rotates(O),
        MBDyn_node(N,STRUCTURAL)
{}

// Stuctural node - DUMMY definition
MBDyn_node_structural::MBDyn_node_structural (Id N,
		       MBDyn_node_structural::Type DUMMY,
		       Id sn, RVec3 RO, RMat3x3 RRM ) :
        MBDyn_node(N,STRUCTURAL),StrNodeLabel(sn),
        Mode_type(DUMMY),
        Relative_offset(RO), Relative_rot_matrix(RRM)
{}

MBDyn_node_structural::MBDyn_node_structural (Id N) :
                       MBDyn_node(N,STRUCTURAL),StrNodeLabel(0),
                       Mode_type(DYNAMIC)
{}

MBDyn_node_structural::~MBDyn_node_structural() {}

ostream& MBDyn_node_structural::Restart(ostream& out) const
{
   const char indent[]="                 ";
   out << "  structural: " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   if (Mode_type!=DUMMY) {
     out << GetMode() << ", " << endl << indent,
	Abs_Pos.Write(out,", ") << ", " 
	<< endl, Abs_Rot_Matrix.RWrite(out,", ","\n",indent) << ", "
	<< endl << indent,
	Abs_Vel.Write(out,", ") << ", " 
	<< endl << indent,
	Abs_Ang_Vel.Write(out,", ") << ", " << endl << indent
	<< "assembly, "
	<< position_initial_stiffness << ", "
	<< velocity_initial_stiffness << ", "
	<< Omega_rotates << ";" << endl;      
   }
   else {
      out << "dummy, " << StrNodeLabel << ", offset, " << endl
	<< "reference, relative, ",
      Relative_offset.Write(out,", ") << ", "
      << " reference, relative, " << endl,
      Relative_rot_matrix.RWrite(out,", ","\n",indent)
      << ";" << endl;
   }
   return out;
}

inline const char* const MBDyn_node_structural::Gettype(void) const
{
   return "STRUCTURAL NODE";
}

inline const char* const MBDyn_node_structural::GetMode(void) const
{
   if (Mode_type==MBDyn_node_structural::STATIC)
     return "static";
   else if (Mode_type==MBDyn_node_structural::DYNAMIC)
     return "dynamic";
}

ostream& MBDyn_node_structural::Test(ostream& out) const
{
   out << "<- structural node -<";
   MBDyn_node::Test(out);
   return out;
}


// ELECTRIC NODE


MBDyn_node_electric::MBDyn_node_electric (Id N,double x, double xp) :
                     Dx(x),Dxp(xp),MBDyn_node(N,ELECTRIC)
{}

MBDyn_node_electric::~MBDyn_node_electric () {}


ostream& MBDyn_node_electric::Restart(ostream& out) const
{
   const char* indent = "            ";
   out << " electric : " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   out << endl << indent << Dx << ", "
       << Dxp << ";" << endl;
   return out;
}

inline const char* const MBDyn_node_electric::Gettype(void) const
{
   return "ELECTRIC NODE";
}

ostream& MBDyn_node_electric::Test(ostream& out) const
{
   out << "<- electric node -<";
   MBDyn_node::Test(out);
   return out;
}

// ABSTRACT NODE

MBDyn_node_abstract::MBDyn_node_abstract (Id N,double x,double xp) :
                     Dx(x),Dxp(xp),MBDyn_node(N,ABSTRACT)
{}

ostream& MBDyn_node_abstract::Restart (ostream& out) const
{
   const char* indent = "            ";
   out << " abstract : " << Label << ", ";
   if (_title_!=NULL) out << "name, " << _title_ << ", ";
   out << endl << indent << Dx << ", "
       << Dxp << ";" << endl;
   return out;
}

inline const char* const MBDyn_node_abstract::Gettype(void) const
{
   return "ABSTRACT NODE";
};

MBDyn_node_abstract::~MBDyn_node_abstract ()
{}

ostream& MBDyn_node_abstract::Test(ostream& out) const
{
   out << "<- abstract node -<";
   MBDyn_node::Test(out);
   return out;
}

// GRAVITY

MBDyn_gravity::MBDyn_gravity (Vec3 A,MBDyn_drive_CONST* Iacc) : Acc(A), 
MBDyn_elem(0,MBDyn_elem::_GRAVITY), Intensity(Iacc)
{}

MBDyn_gravity::~MBDyn_gravity ()
{}

ostream& MBDyn_gravity::Restart(ostream& out) const
{
   out << "  gravity:  ",Acc.Write(out,", ") << ", ";
   Intensity->Restart(out);
   out << ";" << endl;
   return out;
}

inline const char* const MBDyn_gravity::Gettype(void) const
{
   return "GRAVITY ELEMENT";
}

ostream& MBDyn_gravity::Test(ostream& out) const
{
   out << "<- gravity element -<";
   MBDyn_elem::Test(out);
   return out;
}

// GENERIC RELATED FUNCTIONS RELATIVE TO MBDYN OUTPUT AND TRANSLATE

// DRIVES AND DRIVE CALLER

MBDyn_drive_caller::MBDyn_drive_caller (MBDyn_drive_caller::Type T) :
                                        Drive_Type (T), _remark_(NULL)
{}

MBDyn_drive_caller::~MBDyn_drive_caller () {}

inline const char* const MBDyn_drive_caller::Gettype (void) const
{
   return "<-drive caller-<";
}

void MBDyn_drive_caller::Remark (char* R)
{
   _remark_=new char[strlen(R)+1];
   for (int i=0;i<strlen(R);i++) _remark_[i]=' ';
   for (int i=0;i<strlen(R);i++) {
      _remark_[i]=R[i];
   }
   _remark_[strlen(R)]='\0';
   return;
}

ostream& MBDyn_drive_caller::Restart (ostream& out) const
{
   if (REMOVE_REMARK==N) out << _remark_;
   return out;
}

// Constant Drive

MBDyn_drive_CONST::MBDyn_drive_CONST (double cf) : MBDyn_drive_caller(CONST),
							       const_coef(cf)
{}

MBDyn_drive_CONST::~MBDyn_drive_CONST ()
{}

ostream& MBDyn_drive_CONST::Restart (ostream& out) const
{
   MBDyn_drive_caller::Restart(out);
   out << "const, " << const_coef;
   return out;
}

ostream& MBDyn_drive_CONST::Test (ostream& out) const
{
   /*NullOP*/
   return out;
}

inline const char* const MBDyn_drive_CONST::Gettype(void) const
{
   return "<-const-<";
}

// Linear drive

MBDyn_drive_LINEAR::MBDyn_drive_LINEAR (double cc, double sc) :
                                        MBDyn_drive_caller (LINEAR),
                                        const_coef (cc), slope_coef (sc)
{}

MBDyn_drive_LINEAR::~MBDyn_drive_LINEAR ()
{}

ostream& MBDyn_drive_LINEAR::Restart (ostream& out) const
{
   MBDyn_drive_caller::Restart(out);
   out << "linear, " << const_coef << ", " << slope_coef;
   return out;
}

ostream& MBDyn_drive_LINEAR::Test (ostream& out) const
{
   /*NullOP*/
   return out;
}

inline const char* const MBDyn_drive_LINEAR::Gettype (void) const
{
   return "<-linear-<";
}

// Parabolic drive

MBDyn_drive_PARABOLIC::MBDyn_drive_PARABOLIC (double cc, double lc,
					      double pc) :
                                              MBDyn_drive_caller (PARABOLIC),
                                              const_coef (cc),
                                              linear_coef (lc),
                                              parabolic_coef (pc)
{}

MBDyn_drive_PARABOLIC::~MBDyn_drive_PARABOLIC ()
{}

ostream& MBDyn_drive_PARABOLIC::Restart (ostream& out) const
{
   MBDyn_drive_caller::Restart(out);
   out << "parabolic, " << const_coef << ", " <<
     linear_coef << ", " << parabolic_coef;
   return out;
}

ostream& MBDyn_drive_PARABOLIC::Test (ostream& out) const
{
   /*NullOP*/
   return out;
}

inline const char* const MBDyn_drive_PARABOLIC::Gettype (void) const
{
   return "<-parabolic-<";
}

// Cubic drive

MBDyn_drive_CUBIC::MBDyn_drive_CUBIC (double cc, double lc,
					      double pc, double cbc) :
                                              MBDyn_drive_caller (CUBIC),
                                              const_coef (cc),
                                              linear_coef (lc),
                                              parabolic_coef (pc),
                                              cubic_coef (cbc)
{}

MBDyn_drive_CUBIC::~MBDyn_drive_CUBIC ()
{}

ostream& MBDyn_drive_CUBIC::Restart (ostream& out) const
{
   MBDyn_drive_caller::Restart(out);
   out << "cubic, " << const_coef << ", " <<
     linear_coef << ", " << parabolic_coef <<
     ", " << cubic_coef;
   return out;
}

ostream& MBDyn_drive_CUBIC::Test (ostream& out) const
{
   /*NullOP*/
   return out;
}

inline const char* const MBDyn_drive_CUBIC::Gettype (void) const
{
   return "<-cubic-<";
}

// END OF DRIVERS 


					
