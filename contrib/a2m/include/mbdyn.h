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


#ifndef MBDYN_H
#define MBDYN_H

#include <defs.h>
#include <map.h>
#include <mathem.h>

// DRIVERS AND DRIVE CALLER

struct MBDyn_drive_caller {
   
   enum Type {
      CONST,
	LINEAR,
	PARABOLIC,
	CUBIC,
	STEP,
	DOUBLE_STEP,
	RAMP,
	DOUBLE_RAMP,
	SINE,
	COSINE,
	FREQUENCY_SWEEP,
	EXPONENTIAL,
	RANDOM,
	STRING,
	DOF,
	ARRAY,
	TEMPLATE_DRIVE
   };
      
   Type Drive_Type;
   char* _remark_;
   char* _title_;
   
   MBDyn_drive_caller (Type);
   ~MBDyn_drive_caller ();
   virtual ostream& Restart (ostream& out) const;
   virtual ostream& Test(ostream& out) const=0;
   inline const char* const Gettype (void) const;
   void Remark(char*);
};

struct MBDyn_drive_CONST : public MBDyn_drive_caller {
  
   double const_coef;
   
   MBDyn_drive_CONST (double coef);
   ~MBDyn_drive_CONST ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
};

struct MBDyn_drive_LINEAR : public MBDyn_drive_caller {
   
   double const_coef,slope_coef;
   
   MBDyn_drive_LINEAR (double, double);
   ~MBDyn_drive_LINEAR ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
};

struct MBDyn_drive_PARABOLIC : public MBDyn_drive_caller {
   
   double const_coef, linear_coef, parabolic_coef;
   
   MBDyn_drive_PARABOLIC (double, double, double);
   ~MBDyn_drive_PARABOLIC ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
};

struct MBDyn_drive_CUBIC : public MBDyn_drive_caller {
   
   double const_coef, linear_coef, parabolic_coef, cubic_coef;
   
   MBDyn_drive_CUBIC (double, double, double, double);
   ~MBDyn_drive_CUBIC ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
};


// SCHEDA GENERICA DI MBDYN

struct MBDyn_card {
   unsigned int Label;
   char* _remark_;
   char* _title_;
   
   enum Type {
      _ELEMENT,
      _REFERENCE,
      _NODE,
      _CONSTLAW,
      _DRIVE,
      
      _LAST_CARD
   };

   Type Card_type;
   
   MBDyn_card ();
   MBDyn_card (Id, MBDyn_card::Type);
   ~MBDyn_card ();
   virtual ostream& Restart(ostream& out) const=0;
   virtual inline const char* const Gettype() const=0;
   virtual ostream& Test(ostream& out) const;
   void Remark (char*);
   void Title (char*);
};

struct MBDyn_elem : public MBDyn_card {
   enum Type {
      UNKNOWN = 0,
      _BEAM,
      _JOINT,
      _RIGIDBODY,
      _BULK,
      _FORCE,
      _GRAVITY
   };
   
   Type Elem_type;
   
   MBDyn_elem (Id,MBDyn_elem::Type);
   MBDyn_elem ();
   ~MBDyn_elem ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype() const;
};


// LEGGE COSTITUTIVA DEL MATERIALE
struct ConstitutiveLaw6D : public Matrix {
   ConstitutiveLaw6D (double,double);
   ~ConstitutiveLaw6D ();
   double E;
   double G;
   double v;
   
   ostream& Restart(ostream& out, const char*) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype() const;
};

typedef ConstitutiveLaw6D* pCL6D;

//

// TRAVE DI MBDYN

struct MBDyn_beam : public MBDyn_elem {
   MBDyn_beam();
   MBDyn_beam (Id,Id,Id,Id,Vec3,Vec3,Vec3,Mat3x3,Mat3x3,Mat6x6*,Mat6x6*);
   ~MBDyn_beam();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype() const;

   enum NumNodes { NODE1=0, NODE2=1, NODE3=2, NUMNODES=3 };
   enum NumSections { S_I=0, S_II=1, NUMSEZ=2 };
   
   Id Node[3];      			// Posizione dei nodi
   Vec3 f[3];			        // Offset dei nodi
   Mat3x3 R[2];			        // Matrici di rotazione
   Mat6x6* pD[2];		        // Legge costitutiva nella sezione
};

// CORPO RIGIDO DI MBDYN

struct MBDyn_body : public MBDyn_elem {
   MBDyn_body ();
   MBDyn_body (Id,Id,double,RVec3,RMat3x3);
   MBDyn_body (Id,Id,double,RVec3,RMat3x3,Id);
   MBDyn_body (Id,Id,double,RVec3,RMat3x3,Mat3x3);
   ~MBDyn_body ();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype() const;
   
   Id Node;
   double Mass;			       // Valore della massa
   RVec3 Xgc;			       // Posizione del centro di massa
   RMat3x3 Jota;                       // Matrice delle proprietà inerziali
   Boolean Inertial;                   // Nel caso sia definito
   Id NodeInertial;
   Boolean defNodeInertial;            // Si definisce in uno dei 2
   Mat3x3 RotInertial;
};

// FORZA E COPPIA IN MBDYN

struct MBDyn_force : public MBDyn_elem {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   inline const char* const GetMode(void) const;
   
   enum Type {
      UNKNOWN=-1,
      ABSTRACT=0,
      STRUCTURALFORCE,
      STRUCTURALCOUPLE,
     
      LAST_FORCE_TYPE
   };
   
   enum FType {
      CONSERVATIVE,
      FOLLOWER
   };
   
   Type Force_type;
   FType Mode_type;
   
   MBDyn_force();
   MBDyn_force(Id,MBDyn_force::Type);
   MBDyn_force(Id,MBDyn_force::Type,MBDyn_force::FType);
   ~MBDyn_force();
};

struct MBDyn_force_abstract : public MBDyn_force {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Nodedof;
   MBDyn_drive_caller* Magnitude;
   
   MBDyn_force_abstract (Id,Id,MBDyn_drive_caller*);
   ~MBDyn_force_abstract ();
};

struct MBDyn_force_structural : public MBDyn_force {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node;
   RVec3 Relative_direction;
   Vec3 Relative_arm;
   MBDyn_drive_caller* Magnitude;
   
   MBDyn_force_structural (Id,FType,Id,RVec3,Vec3,MBDyn_drive_caller*);
   ~MBDyn_force_structural ();
};

struct MBDyn_force_couple : public MBDyn_force {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node;
   RVec3 Relative_direction;
   MBDyn_drive_caller* Magnitude;

   MBDyn_force_couple (Id,FType,Id,RVec3,MBDyn_drive_caller*);
   ~MBDyn_force_couple ();
   
};

// JOINT DI MBDYN
struct MBDyn_joint : public MBDyn_elem {
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype() const;
   
   enum Type {      
        UNKNOWN = -1,
	DISTANCE = 0,
	DISTANCEWITHOFFSET,
	CLAMP,
	SPHERICALHINGE,
	PIN,
	UNIVERSALHINGE,
	UNIVERSALPIN,
	PLANEHINGE,
	PLANEPIN,
	AXIALROTATION,
	INPLANESIMPLE,
	INPLANECONTACT,
	IN_LINE,
	ROD,
	PRISMATIC,
	DRIVEHINGE,
	MODAL,
	
	LASTJOINTTYPE
   };
   
   Type Joint_type;

   MBDyn_joint ();
   MBDyn_joint (Id,MBDyn_joint::Type q);
   ~MBDyn_joint ();
};

// SOTTOTIPI JOINT DI MBDYN
// SPHERICAL JOINT

struct MBDyn_spherhinge : public MBDyn_joint {
   MBDyn_spherhinge(Id,Id,Id,Vec3,Mat3x3,Vec3,Mat3x3);
   ~MBDyn_spherhinge();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node[2];
   Vec3 d[2];
   Mat3x3 Rh[2];
};

// PIN JOINT

struct MBDyn_pin : public MBDyn_joint {
   MBDyn_pin (Id,Id,Vec3,Vec3);
   ~MBDyn_pin ();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Id Node;
   Vec3 Relative_offset;
   Vec3 Absolute_pos;
};

// CLAMP

struct MBDyn_clamp : public MBDyn_joint {
   MBDyn_clamp (Id,Id);
   MBDyn_clamp (Id,Id,Vec3);
   MBDyn_clamp (Id,Id,Mat3x3);
   MBDyn_clamp (Id,Id,Vec3,Mat3x3);
   ~MBDyn_clamp ();
   
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node;
   Boolean NodeRef1;
   Boolean NodeRef2;
   Vec3 Absolute_pos;
   Mat3x3 Absolute_rot_matrix;
};

// UNIVERSAL HINGE

struct MBDyn_univhinge : public MBDyn_joint {
   MBDyn_univhinge(Id,Id,Id,Vec3,Mat3x3,Vec3,Mat3x3);
   ~MBDyn_univhinge();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node[2];
   Vec3 d[2];
   Mat3x3 Rh[2];
};

// UNIVERSAL PIN

struct MBDyn_univpin : public MBDyn_joint {
   MBDyn_univpin (Id,Id,Vec3,Mat3x3,Vec3,Mat3x3);
   MBDyn_univpin (Id,Id,Vec3,Vec3,Mat3x3);
   MBDyn_univpin (Id,Id,Vec3,Mat3x3,Vec3);
   MBDyn_univpin (Id,Id,Vec3,Vec3);
   ~MBDyn_univpin ();
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Id NodeRef;
   Vec3 Relative_offset;
   Mat3x3 Relative_rot_matrix;
   Vec3 Absolute_pos;
   Mat3x3 Absolute_pin_rot;
   Boolean Hinge1,Hinge2;
};

// PLANE HINGE

struct MBDyn_planehinge : public MBDyn_joint {
   MBDyn_planehinge(Id,Id,Id,Vec3,Mat3x3,Vec3,Mat3x3);
   ~MBDyn_planehinge();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   Id Node[2];
   Vec3 d[2];
   Mat3x3 Rh[2];
};

// PLANE PIN

struct MBDyn_planepin : public MBDyn_joint {
   MBDyn_planepin (Id,Id,Vec3,Mat3x3,Vec3,Mat3x3);
   MBDyn_planepin (Id,Id,Vec3,Vec3,Mat3x3);
   MBDyn_planepin (Id,Id,Vec3,Mat3x3,Vec3);
   MBDyn_planepin (Id,Id,Vec3,Vec3);
   ~MBDyn_planepin ();
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Id NodeRef;
   Vec3 Relative_offset;
   Mat3x3 Relative_rot_matrix;
   Vec3 Absolute_pos;
   Mat3x3 Absolute_pin_rot;
   Boolean Hinge1,Hinge2;
};

// IN PLANE

struct MBDyn_inplane : public MBDyn_joint {
   MBDyn_inplane (Id,Id,Id,Vec3,Vec3,Vec3);
   MBDyn_inplane (Id,Id,Id,Vec3,Vec3);
   ~MBDyn_inplane ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Id Node[2];
   Vec3 Relative_position;
   Vec3 Relative_norm_direction;
   Vec3 Relative_offset;
};

// IN LINE

struct MBDyn_inline : public MBDyn_joint {
   MBDyn_inline (Id,Id,Id,Vec3,Mat3x3,Vec3);
   MBDyn_inline (Id,Id,Id,Vec3,Mat3x3);
   ~MBDyn_inline ();
   ostream& Restart (ostream& out) const;
   ostream& Test (ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Id Node[2];
   Vec3 Relative_line_position;
   Mat3x3 Relative_orientation;
   Vec3 Relative_offset;
};

// BULK ELEMENTS

struct MBDyn_bulk : public MBDyn_elem {
   enum Type {
      STIFFNESS_SPRING,
      
      LAST_BULK_TYPE
   };
   
   Type Bulk_type;
   
   MBDyn_bulk (Id,MBDyn_bulk::Type);
   ~MBDyn_bulk ();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
};

struct MBDyn_StiffnessSpring : public MBDyn_bulk {
   MBDyn_StiffnessSpring (Id, Id, double);
   ~MBDyn_StiffnessSpring ();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;

   double Stiffness;
   Id NodeDof;
};


// REFERENCE DI MBDYN

struct MBDyn_reference : public MBDyn_card {
   MBDyn_reference(Id,RVec3,RMat3x3,Boolean,RVec3,RVec3);
   MBDyn_reference();
   MBDyn_reference(Id);
   ~MBDyn_reference();
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype(void) const;
   
   // Posizione, Matrice rotazione, Velocità, Velocità angolare 
   // Il tutto nel sistema di riferimento assoluto
   RVec3 Abs_Pos;				
   RMat3x3 Abs_Rot_Matrix;
   Boolean Eye_Flag;
   RVec3 Abs_Vel;
   RVec3 Abs_Omega;
};

// NODO DI MBDYN

struct MBDyn_node : public MBDyn_card {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   enum Type {
      UNKNOWN = -1,
      ABSTRACT = 0,
      STRUCTURAL,
      ELECTRIC,
      PARAMETER,
      HYDRAULIC,
      
      LAST_NODE_TYPE 
   };

   Type Node_type;
   
   MBDyn_node ();
   MBDyn_node (Id,MBDyn_node::Type);
   ~MBDyn_node ();
};

// NODI CORRELATI

struct MBDyn_node_structural : public MBDyn_node {
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   inline const char* const GetMode (void) const;
   
   enum Type {
      STATIC,
      DYNAMIC,
      DUMMY
   };
   
   Type Mode_type;
   RVec3 Abs_Pos;
   RMat3x3 Abs_Rot_Matrix;
   RVec3 Abs_Vel;
   RVec3 Abs_Ang_Vel;
   double position_initial_stiffness;
   double velocity_initial_stiffness;
   Boolean Omega_rotates;
   // DUMMY CASE
   Id StrNodeLabel;
   RVec3 Relative_offset;
   RMat3x3 Relative_rot_matrix;
   
   MBDyn_node_structural (Id,Type,RVec3,RMat3x3,RVec3,RVec3,
			  double,double,Boolean);
   MBDyn_node_structural (Id,MBDyn_node_structural::Type,Id,RVec3,RMat3x3);
   MBDyn_node_structural (Id);
   ~MBDyn_node_structural ();

};

struct MBDyn_node_electric : public MBDyn_node {
   MBDyn_node_electric (Id,double,double);
   ~MBDyn_node_electric ();
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;

   double Dx;
   double Dxp;
};

struct MBDyn_node_abstract : public MBDyn_node {
   MBDyn_node_abstract (Id,double,double);
   ~MBDyn_node_abstract ();
   ostream& Restart (ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   double Dx;
   double Dxp;
};

// GRAVITY

struct MBDyn_gravity : public MBDyn_elem {
   MBDyn_gravity (Vec3,MBDyn_drive_CONST*);
   ~MBDyn_gravity ();
   ostream& Restart(ostream& out) const;
   ostream& Test(ostream& out) const;
   inline const char* const Gettype (void) const;
   
   Vec3 Acc;
   MBDyn_drive_CONST* Intensity;
};

#endif
