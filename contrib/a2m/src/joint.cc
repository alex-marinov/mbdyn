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


//                                JOINT.CC                                    

//Definizioni per la card joint

#include <joint.h>
#include <output.h>

extern MBDyn_deck MBJoints;
extern MBDyn_deck MBReference;
extern MBDyn_deck MBNodes;
extern PTN_deck PartToNode_Table;
extern MTP_deck Marker_Table;
extern MTR_deck Reference_Table;

  
s_joint::s_joint(void)    : Joint_Card(NULL),  
                           _Node1(N),_Node2(N),_Joint_Card(N),_Joint_Type(N)
                          { Node[0]=0; Node[1]=0; }

s_joint::~s_joint(void)   {}

Joint_ref::Joint_ref(void): Maximum_deformation(0), Mu_dyn_rot(0),
                             Inner_radius(0),Outer_radius(0),
                             Mu_stat_rot(0),Preload_rad(0),Preload_axial(0),
                             Width(0),Pitch(0),Mu_dyn_trans(0),Mu_stat_trans(0),
                             Delta_v(0),Preload_x(0),Preload_y(0),
                             Height(0),Friction_type(_OFF),
                             Max_fric_rot(0),

                             _Maximum_deformation(N),_Mu_dyn_rot(N),
                             _Ic_rot(N),_Ic_tran(N),_Ic(N),
                             _Inner_radius(N),_Outer_radius(N),
                             _Mu_stat_rot(N), _Preload_rad(N), _Preload_axial(N),
                             _Width(N),_Pitch(N),_Mu_dyn_trans(N),_Mu_stat_trans(N),
                             _Delta_v(N),_Preload_x(N),_Preload_y(N),_Height(N),
                             _Friction_type(N),_Max_fric_rot(N)
                             {}

Joint_ref::~Joint_ref(void)  {}

Boolean s_joint::Test()
{
   const int err_before=nerr;
   Joint_ref* a = Joint_Card;
   // Controllo preliminare sui nodi
   if (_Node1==0) out_error (19,"NODE I");
   if (_Node2==0) out_error (19,"NODE J");
   // Controllo del tipo di vincolo
   if (Joint_Card != 0) {
      // Controlla che i tipi di vincolo semplice non abbiano dati
      // rindondanti derivanti dalla card di estensione
      if (Joint_Type==_CONVEL) out_error (20,"CONVEL");
      if (Joint_Type==_FIXED) out_error (20,"FIXED");
      if (Joint_Type==_HOOKE) out_error (20,"HOOKE");
      if (Joint_Type==_PLANAR) out_error (20,"PLANAR");
      if (Joint_Type==_SPHERICAL) out_error (20,"SPHERICAL");
      if (Joint_Type==_UNIVERSAL) out_error (20,"UNIVERSAL");
      // CYLINDRICAL
      if (Joint_Type==_CYLINDRICAL) {
	 // Controllo parametri rindondanti
	 if (a->_Pd==Y) out_error (21,"PD");
	 if (a->_Delta_v==Y) out_error (21,"DELTA_V");
	 if (a->_Inner_radius==Y) out_error (21,"INNER RADIUS");
	 if (a->_Friction_type==Y) out_error (21,"FRICTION");
	 if (a->_Maximum_deformation) out_error (21,"MAXIMUM DEFORMATION");
	 if (a->_Mu_dyn_rot==Y) out_error (21,"MU_DYN_ROT");
	 if (a->_Mu_stat_rot==Y) out_error (21,"MU_STAT_ROT");
	 if (a->_Outer_radius==Y) out_error (21,"OUTER_RADIUS");
	 if (a->_Preload_rad==Y) out_error (21,"PRELOAD_RAD");
	 if (a->_Preload_axial==Y) out_error (21,"PRELOAD_AXIAL");
	 if (a->_Height==Y) out_error (21,"HEIGHT");
	 if (a->_Mu_dyn_trans) out_error (21,"MU_DYN_TRANS");
	 if (a->_Mu_stat_trans) out_error (21,"MU_STAT_TRANS");
	 if (a->_Preload_x) out_error (21,"PRELOAD_X");
	 if (a->_Preload_y) out_error (21,"PRELOAD_Y");
	 if (a->_Width) out_error(21,"WIDTH");
	 if (a->_Pitch) out_error(21,"PITCH");
	 if (a->_Ic) out_error (21,"IC");
      }
      if (Joint_Type==_RACKPIN) {
	 // Controllo parametri rindondanti
	 if (a->_Ic_tran==Y) out_error (22,"ICTRAN");
	 if (a->_Ic_rot==Y) out_error (22,"ICROT");
	 if (a->_Delta_v==Y) out_error (22,"DELTA_V");
	 if (a->_Inner_radius==Y) out_error (22,"INNER RADIUS");
	 if (a->_Friction_type==Y) out_error (22,"FRICTION");
	 if (a->_Maximum_deformation) out_error (22,"MAXIMUM DEFORMATION");
	 if (a->_Mu_dyn_rot==Y) out_error (22,"MU_DYN_ROT");
	 if (a->_Mu_stat_rot==Y) out_error (22,"MU_STAT_ROT");
	 if (a->_Outer_radius==Y) out_error (22,"OUTER_RADIUS");
	 if (a->_Preload_rad==Y) out_error (22,"PRELOAD_RAD");
	 if (a->_Preload_axial==Y) out_error (22,"PRELOAD_AXIAL");
	 if (a->_Height==Y) out_error (22,"HEIGHT");
	 if (a->_Mu_dyn_trans) out_error (22,"MU_DYN_TRANS");
	 if (a->_Mu_stat_trans) out_error (22,"MU_STAT_TRANS");
	 if (a->_Preload_x) out_error (22,"PRELOAD_X");
	 if (a->_Preload_y) out_error (22,"PRELOAD_Y");
	 if (a->_Width) out_error(22,"WIDTH");
	 if (a->_Pitch) out_error(22,"PITCH");
	 if (a->_Ic) out_error (22,"IC");
      }
      if (Joint_Type==_REVOLUTE) {
	 // Controllo parametri rindondanti
	 if (a->_Ic_tran) out_error (23,"ICTRAN");
	 if (a->_Ic_rot) out_error (23,"ICROT");
	 if (a->_Pd) out_error (23,"PD");
	 if (a->_Pitch) out_error (23,"PITCH");
	 if (a->_Height) out_error (23,"HEIGHT");
	 if (a->_Width) out_error (23,"WIDHT");
	 if (a->_Preload_x) out_error (23,"PRELOAD X");
	 if (a->_Preload_y) out_error (23,"PRELOAD Y");
	 if (a->_Mu_dyn_trans) out_error (23,"MU_DYN_TRANS");
	 if (a->_Mu_stat_trans) out_error (23,"MU_STAT_TRANS");
      }
      if (Joint_Type==_SCREW) {
	 // Controllo parametri rindondanti
	 if (a->_Ic_tran==Y) out_error (24,"ICTRAN");
	 if (a->_Ic_rot==Y) out_error (24,"ICROT");
	 if (a->_Delta_v==Y) out_error (24,"DELTA_V");
	 if (a->_Inner_radius==Y) out_error (24,"INNER RADIUS");
	 if (a->_Friction_type==Y) out_error (24,"FRICTION");
	 if (a->_Maximum_deformation) out_error (24,"MAXIMUM DEFORMATION");
	 if (a->_Mu_dyn_rot==Y) out_error (24,"MU_DYN_ROT");
	 if (a->_Mu_stat_rot==Y) out_error (24,"MU_STAT_ROT");
	 if (a->_Outer_radius==Y) out_error (24,"OUTER_RADIUS");
	 if (a->_Preload_rad==Y) out_error (24,"PRELOAD_RAD");
	 if (a->_Preload_axial==Y) out_error (24,"PRELOAD_AXIAL");
	 if (a->_Height==Y) out_error (24,"HEIGHT");
	 if (a->_Mu_dyn_trans) out_error (24,"MU_DYN_TRANS");
	 if (a->_Mu_stat_trans) out_error (24,"MU_STAT_TRANS");
	 if (a->_Preload_x) out_error (24,"PRELOAD_X");
	 if (a->_Preload_y) out_error (24,"PRELOAD_Y");
	 if (a->_Width) out_error(24,"WIDTH");
	 if (a->_Ic) out_error (24,"IC");
	 if (a->_Pd) out_error (24,"PD");
      }
      if (Joint_Type==_TRANSLATIONAL) {
	 // Controllo parametri rindondanti
	 if (a->_Ic_tran==Y) out_error (25,"ICTRAN");
	 if (a->_Ic_rot==Y) out_error (25,"ICROT");
	 if (a->_Pd==Y) out_error (25,"PD");
	 if (a->_Inner_radius==Y) out_error(25,"INNER RADIUS");
	 if (a->_Mu_dyn_rot==Y) out_error(25,"MU_DYN_ROT");
	 if (a->_Mu_stat_rot==Y) out_error(25,"MU_STAT_ROT");
	 if (a->_Outer_radius==Y) out_error (25,"OUTER_RADIUS");
	 if (a->_Preload_rad==Y) out_error(25,"PRELOAD RAD");
	 if (a->_Preload_axial==Y) out_error(25,"PRELOAD AXIAL");
	 if (a->_Pitch==Y) out_error (25,"PITCH");
      }
   }
   else {
      // Controlla che sia presente l'estensione per i vincoli che
      // necessitano di più parametri
      if (Joint_Type==_RACKPIN) out_error (26,"RACKPIN");
      if (Joint_Type==_SCREW) out_error (26,"SCREW");
   }
   if (err_before != nerr) return Y; else return N;
}


inline const char* const s_joint::Gettype(void) const
{
   return "JOINT";
}


ostream& s_joint::Print (ostream& out) const
{
   out << endl;
   out << "JOINT:" << label << "        Type:" << Joint_Type;
   out << endl;
   out << "     Node1 [" << _Node1 << "] = " << Node[0] << endl
       << "     Node2 [" << _Node2 << "] = " << Node[1] << endl;
   if (Joint_Card != NULL) Joint_Card->Print(out);
   return out;
}

ostream& Joint_ref::Print (ostream& out) const
{
   out << "     Ic [" << _Ic << "] = ",Ic.Write(out,", ") << endl;
   out << "     Ictran [" << _Ic_tran << "] = ",Ic_tran.Write(out,", ") << endl;
   out << "     Icrot  [" << _Ic_rot << "] = ",Ic_rot.Write(out,", ") << endl;
   out << "     Inner radius [" << _Inner_radius << "] = " << Inner_radius << endl;
   out << "     Outer radius [" << _Outer_radius << "] = " << Outer_radius << endl;
   out << "     Maximum deformation [" << _Maximum_deformation << "] = " << Maximum_deformation << endl;
   out << "     Delta V [" << _Delta_v << "] = " << Delta_v << endl;
   out << "     Friction [" << _Friction_type << "] = " << Friction_type;
   out << "     Mu dyn rot [" << _Mu_dyn_rot << "] = "  << Mu_dyn_rot << endl;
   out << "     Mu stat rot [" << _Mu_stat_rot << "] = " << Mu_stat_rot << endl;
   out << "     Max fric rot [" << _Max_fric_rot << "] = " << Max_fric_rot << endl;
   out << "     Preload rad ["   << _Preload_rad << "] = " << Preload_rad << endl;
   out << "     Preload axial [" << _Preload_axial << "] = " << Preload_axial << endl;
   out << "     Preload x [" << _Preload_x << "] = " << Preload_x << endl;
   out << "     Preload y [" << _Preload_y << "] = " << Preload_y << endl;
   out << "     Mu dyn trans [" << _Mu_dyn_trans << "] = " << Mu_dyn_trans << endl;
   out << "     Mu stat trans [" << _Mu_stat_trans << "] = " << Mu_stat_trans << endl;
   out << "     Pd [" << _Pd << "] = " << Pd << endl;
   out << "     Width [" << _Width << "] = " << Width << endl;
   out << "     Pitch [" << _Pitch << "] = " << Pitch << endl;
   out << "     Height [" << _Height << "] = " << Height << endl;
   out << endl;
   return out;
}

void s_joint::Translate (ostream& out)
{
   /* Common query of data */
   char* comment= new char[255];
   char* title=new char[80];
   Id JointID;
   Id REFNODE[2];
   Id NODE[2];
   Id MARKER[2];
   Vec3 RefPos[2];
   Mat3x3 RefMat[2];
   MBDyn_reference *MBR[2];
   MBDyn_node_structural *MNS[2];
   /* Build the relations */
   for (int i = 0; i<2; i++) {
      
      /* debuggin' code (14 September 2000) */
      
      /* DETERMINA LA PARTE ADAMS CORRISPONDENTE  AL MARKER ADAMS NODE [I] */
      /* E NE PRELEVA SOLO IL NUMERO (IL TIPO E' PART O POINTMASS)*/
      
      REFNODE[i]= ((*(Marker_Table.find(Node[i]))).second).Num;
      /* TROVA A QUALE NODO STRUTTURALE MBDYN CORRISPONDE LA PARTE TROVATA */
      NODE[i]=(*(PartToNode_Table.find(REFNODE[i]))).second;
      
      /* DETERMINA IL REFERENCE MBDYN CORRISPONDENTE AL MARKER ADAMS I     */
      /* FINTANTO CHE L'OPERAZIONE DI QUEUEING VALE LE LABEL SONO =        */
      MARKER[i]=(*(Reference_Table.find(Node[i]))).second;
      
      /* PUNTATORE AL MARKER ADAMS - REF MBDYN CON LABEL I                 */
      MBR[i]=(MBDyn_reference*) Find_MBCard (MARKER[i],MBReference);
      
      /* PUNTATORE AL NODO STRUTTURALE MBDYN                               */
      MNS[i]=(MBDyn_node_structural*) Find_MBCard (NODE[i],MBNodes);
      
      /* DEREFERENZIAZIONE DELLA ORIGINE E DELLA MATRICE DI ROTAZIONE DEL
       * SISTEMA MARKER RISPETTO ALLA SUA PARTE */      
      RefPos[i]=Unref(MBR[i]->Abs_Pos);
      RefMat[i]=Unref(MBR[i]->Abs_Rot_Matrix);

   }
   
   /* CREA LA LABEL PER IL NUOVO JOINT */
   JointID = GetFreeLabel (MBJoints,label);
   
   switch ( Joint_Type ) {
    case _SPHERICAL : { 
       /* Spherical hinge */
       MBDyn_spherhinge* HJ = new
	 MBDyn_spherhinge (JointID, NODE[0], NODE[1], RefPos[0], RefMat[0],
			   RefPos[1], RefMat[1]);
       MBJoints.insert (MBDyn_entry(JointID, (MBDyn_card*) HJ));
       
       sprintf (comment,
		"Spherical joint %d is related to Adams SPHERICAL JOINT %d",
		JointID,label);
       sprintf (title,"Adams SPHERICAL joint %d",label);
       HJ->Title(title);
       HJ->Remark(comment);
       
    }
      break;
    case _REVOLUTE : {
       /* Plane and universal Hinge and Pin */
       /* Pin = Hinge when a node is grounded. */
       MBDyn_planehinge* HJ = new
	 MBDyn_planehinge(JointID,NODE[0],NODE[1], RefPos[0], RefMat[0],
			RefPos[1],RefMat[1]);
       MBJoints.insert (MBDyn_entry(JointID, (MBDyn_card*) HJ));
       
       sprintf (comment,
		"Plane hinge joint %d is related to Adams REVOLUTE %d",
		JointID,label);
       sprintf (title,"Adams REVOLUTE joint %d",label);
       HJ->Title(title);
       HJ->Remark(comment);
       
    }
      break;
    case _PLANAR : {
       /* Planar is equivalent to in plane joint */
       Vec3 RPP,RND,ROFF;
       RVec3 RNormal_direction = (MBR[1]->Abs_Rot_Matrix.GetVec(3));
       RND = Unref (RNormal_direction);
       RPP.Set (0,0,0);
       RND.Set (0,0,0);
       ROFF.Set (0,0,0);
       MBDyn_inplane* HJ = new
	 MBDyn_inplane (JointID,NODE[0],NODE[1],RPP,RND,ROFF);
       MBJoints.insert (MBDyn_entry(JointID, (MBDyn_card*) HJ));
       
       sprintf (comment,
		"In plane joint %d is related to Adams PLANAR JOINT %d",
		JointID,label);
       sprintf (title,"Adams PLANAR joint %d",label);
       HJ->Title(title);
       HJ->Remark(comment);
       
    }
      break;
    case _TRANSLATIONAL : {
       /* Translational joint is equivalent to inline joint */
       Vec3 RLP,ROFF;
       Mat3x3 ROR;
       MBDyn_inline* HJ = new
	 MBDyn_inline (JointID,NODE[0],NODE[1],RLP,ROR,ROFF);
       MBJoints.insert (MBDyn_entry(JointID, (MBDyn_card*) HJ));
       
       sprintf (comment,
		"Inline joint %d is related to Adams TRANSLATIONAL JOINT %d",
		JointID,label);
       sprintf (title,"Adams TRANSLATIONAL joint %d", label);
       
       HJ->Title(title);
       HJ->Remark(comment);
       
    }
      break;
    case _FIXED : {
       /* Fixed joint is equivalent to "INCASTRO" [ */
    }
      break;
   }  
   return;
}

void Joint_ref::Translate(ostream& out)
{
   return;
}

