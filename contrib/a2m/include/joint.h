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

#ifndef JOINT_H
#define JOINT_H

#include <common.h>

struct Joint_ref {
   // QUESTA E' UNA CARD DI ESTENSIONE DELLA CARD JOINT 
   Vec2 Ic;
   double Delta_v;
   double Inner_radius;
   Friction Friction_type;
   double Maximum_deformation;
   double Mu_dyn_rot;
   double Mu_stat_rot;
   double Outer_radius;
   double Preload_rad;
   double Preload_axial;
   Vec2 Ic_rot;
   Vec2 Ic_tran;
   double Width;
   double Pitch;
   double Preload_x;
   double Preload_y;
   double Mu_dyn_trans;
   double Mu_stat_trans;
   double Height;
   double Pd;
   double Max_fric_rot;

   Boolean _Maximum_deformation,_Mu_dyn_rot,_Mu_stat_rot;
   Boolean _Outer_radius,_Preload_rad,_Preload_axial,_Ic_rot;
   Boolean _Ic_tran,_Width,_Pitch,_Preload_x,_Preload_y;
   Boolean _Mu_dyn_trans,_Mu_stat_trans,_Pd;
   Boolean _Delta_v,_Inner_radius;
   Boolean _Ic,_Height,_Friction_type;
   Boolean _Max_fric_rot;

   Joint_ref(void);
   ~Joint_ref(void);
   ostream& Print (ostream& out) const;
   void Translate(ostream& out);
   void Comment (const char*);
   
   template <class T>
     void Set_Part (int param, T value) {
	switch (param) {
	 case HEIGHT: Set_Param(Height,value,_Height); break;
	 case IC: Set_Param (Ic,value,_Ic); break;
	 case PRELOAD_RADIAL: Set_Param (Preload_rad,value,_Preload_rad); break;
	 case PRELOAD_AXIAL: Set_Param (Preload_axial,value,_Preload_axial); break;
	 case ICROT: Set_Param (Ic_rot, value, _Ic_rot); break;
	 case ICTRAN: Set_Param (Ic_tran, value, _Ic_tran); break;
	 case WIDTH: Set_Param (Width, value, _Width); break;
	 case PITCH: Set_Param (Pitch, value, _Pitch); break;
	 case MU_STAT_ROT: Set_Param (Mu_stat_rot, value,_Mu_stat_rot); break;
	 case MU_DYN_ROT: Set_Param (Mu_dyn_rot, value, _Mu_dyn_rot); break;
	 case PRELOAD_X: Set_Param (Preload_x, value, _Preload_x); break;
	 case PRELOAD_Y: Set_Param (Preload_y, value, _Preload_y); break;
	 case OUTER_RADIUS : Set_Param (Outer_radius,value,_Outer_radius); break;
	 case MAXIMUM_DEFORMATION: Set_Param (Maximum_deformation,value,_Maximum_deformation); break;
	 case PD: Set_Param (Pd, value, _Pd); break;
	 case DELTA_V: Set_Param (Delta_v, value, _Delta_v); break;
	 case INNER_RADIUS: Set_Param (Inner_radius, value, _Inner_radius); break;
	 case MU_DYN_TRANS: Set_Param (Mu_dyn_trans, value, _Mu_dyn_trans); break;
	 case MU_STAT_TRANS: Set_Param (Mu_stat_trans, value, _Mu_stat_trans); break;
	 case FRICTION: Set_Param (Friction_type,value,_Friction_type); break;
	 case MAX_FRIC_ROT:Set_Param (Max_fric_rot, value, _Max_fric_rot); break;
	 default : out_error(3,Find_Token(param)); break;
	}
     }    
};


struct s_joint : public s_card {
   Id Node[2];
   Joint Joint_Type;
   Joint_ref* Joint_Card;
   // For joint type specification see defs.h header file
   
   Boolean _Node1,_Node2,_Joint_Type, _Joint_Card;
   
   s_joint(void);
   ~s_joint(void);
   ostream& Print (ostream& out) const; 
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case I     : Set_Param (Node[0], value, _Node1); break;
	 case J     : Set_Param (Node[1], value, _Node2); break;
	 case _SET_ : Set_Param(Joint_Type, value,_Joint_Type); break;
	 default : {
	    if (Joint_Card==NULL) Joint_Card=new Joint_ref; 
	    Joint_Card->Set_Part (param, value);
	    } break;
	}
     }
};

#endif
