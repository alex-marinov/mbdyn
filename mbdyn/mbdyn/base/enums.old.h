/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* enums ed altre dichiarazioni comuni piu' o meno a tutti */


#ifndef ENUMS_H
#define ENUMS_H


#include "myf2c.h"


/* Tipi di Elem. Lasciare sempre UNKNOWN = -1, cosi' il primo elemento
 * ha tipo zero, e l'ultima entry dell'enum, LAST...TYPE, e' uguale
 * al numero di tipi definiti, quindi puo' essere usata come costante nel 
 * dimensionamento degli arrays e come flag di fine tipi. */
class ElemType {
 public:
   enum Type {
      UNKNOWN = -1,
	FORCE = 0,

	AUTOMATICSTRUCTURAL,
	GRAVITY,
	BODY,
	JOINT,
	BEAM,
	PLATE,

	AIRPROPERTIES,
	ROTOR,
	AERODYNAMIC,

	ELECTRICBULK,
	ELECTRIC,
	GENEL,

	HYDRAULIC,
	
	BULK,			
	DRIVEN,
	
	LASTELEMTYPE
   };      
};

extern const char* psElemNames[];
extern const char* psReadControlElems[];


/* Tipi di Joint */
class JointType {
 public:
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
	INPLANE,
	INPLANECONTACT,
	ROD,
	DEFORMABLEHINGE,
	LINEARVELOCITY,
	ANGULARVELOCITY,
	PRISMATIC,
	DRIVEHINGE,
	
	LASTJOINTTYPE
   };
};

extern const char* psJointNames[];


/* Tipi di rods */
class RodType {
 public:
   enum Type {
      UNKNOWN = -1,
	ELASTIC = 0,
	VISCOELASTIC,
	VISCOELASTICWITHOFFSET,
	
	LASTRODTYPE
   };
};

extern const char* psRodNames[];


/* Tipi di cerniere deformabili */
class DefHingeType {
 public:
   enum Type {
      UNKNOWN = -1,
	ELASTIC = 0,
	VISCOUS,
	VISCOELASTIC,
	
	LASTDEFHINGETYPE
   };
};

extern const char* psDefHingeNames[];


/* Tipi di Genel */
class GenelType {
 public:
   enum Type {
      UNKNOWN = -1,
	SWASHPLATE = 0,
	ROTORTRIM,
	CLAMP,
	DISTANCE,
	SPRING,
	SPRINGSUPPORT,
	CROSSSPRINGSUPPORT,
	SPRINGDAMPER,
	SPRINGDAMPERSUPPORT,
	CROSSSPRINGDAMPERSUPPORT,
	MASS,
	SCALARFILTER,
	STATESPACESISO,
	STATESPACEMIMO,
	
	LASTGENELTYPE
   };
};

extern const char* psGenelNames[];


/* Tipi di Force */
class ForceType {
 public:
   enum Type {
      UNKNOWN = -1,
	ABSTRACTFORCE = 0,
	
	CONSERVATIVEFORCE,
	FOLLOWERFORCE,
	CONSERVATIVECOUPLE,
	FOLLOWERCOUPLE,
	
	LASTFORCETYPE
   };
};

extern const char* psForceNames[];


/* Tipi di elementi elettrici */
class ElectricType {
 public:
   enum Type {
      UNKNOWN = -1,
	ACCELEROMETER = 0,
	DISCRETECONTROL,
	
	LASTELECTRICTYPE
   };
};
      
extern const char* psElectricNames[];


/* Tipi di elementi idraulici */
class HydraulicType {
 public:
   enum Type {
      UNKNOWN = -1,
	PIPE = 0,	
	LASTHYDRAULICTYPE
   };
};
      
extern const char* psHydraulicNames[];


/* Tipi di travi */
class BeamType {
 public:
   enum Type {
      UNKNOWN = -1,
	ELASTIC = 0,
	VISCOELASTIC,
	PIEZOELECTRIC,
	
	LASTBEAMTYPE
   };
};

extern const char* psBeamNames[];


/* Tipi di elementi aerodinamici */
class AeroType {
 public:
   enum Type {
      UNKNOWN = -1,
	ROTOR = 0,
	AERODYNAMICBODY,
	AERODYNAMICBEAM,
	
	LASTAEROTYPE
   };
};

extern const char* psAeroNames[];


/* Tipi di rotori */
class RotorType {
 public:
   enum Type {
      UNKNOWN = -1,
	NO = 0,
	UNIFORM,
	GLAUERT,
	MANGLER,
	DYNAMICINFLOW,
	
	LASTROTORTYPE
   };
};

extern const char* psRotorNames[];


/* Tipi di bulk */
class BulkType {
 public:
   enum Type {
      UNKNOWN = -1,
	SPRINGSUPPORT = 0,
	SPRING,
	
	LASTBULKTYPE
   };
};

extern const char* psBulkNames[];


/* Tipi di Drive */
class DriveType {
 public:
   enum Type {
      UNKNOWN = -1,
	FILE = 0,
	
	LASTDRIVETYPE
   };
   
   enum Func {
      UNKNOWNFUNC = -1,
	CONST = 0,
	STEP,
	DOUBLESTEP,
	RAMP,	
	DOUBLERAMP,
	SINE,
	COSINE,
	
	LASTFUNCTYPE
   };
};      
   
extern const char* psDriveNames[];
extern const char* psReadControlDrivers[];

extern const char* psFuncNames[];


/* Tipi di Node */
class NodeType {
 public:
   enum Type {
      UNKNOWN = -1,
	ABSTRACT = 0,

	STRUCTURAL,

	ELECTRIC,

	PARAMETER,

	HYDRAULIC,

	LASTNODETYPE
   };
};    
   
extern const char* psNodeNames[];
extern const char* psReadControlNodes[];
extern const char* psReadNodesNodes[];


/* Tipi di StructNode */
class StructNodeType {
 public:
   enum Type {
      UNKNOWN = -1,
	DYNAMIC = 0,
	STATIC,
	
	LASTSTRUCTNODETYPE
   };
};    
   
extern const char* psStructNodeNames[];


/* Tipi di DofOwner */
class DofType {
 public:
   enum Type {
      UNKNOWN = -1,
	STRUCTURALNODE = 0,
	ELECTRICNODE,
	ABSTRACTNODE,
	HYDRAULICNODE,
	
	JOINT,
	GENEL,
	ROTOR,
	UNSTEADYAERO,
	ELECTRICBULK,
	ELECTRIC,
	HYDRAULIC,
	
	LASTDOFTYPE
   };
};
   
extern const char* psDofOwnerNames[];   


/* ordine dei dof */
class DofOrder {
 public:
   enum Order {
      UNKNOWN = -1,
	ALGEBRAIC = 0,
	DIFFERENTIAL,
	
	LASTORDER
   };
};


#endif // ENUMS_H
