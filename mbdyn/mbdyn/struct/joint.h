/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

/* vincoli, tipo: Elem::Type JOINT */


#ifndef JOINT_H
#define JOINT_H

/* include per derivazione della classe */

#include "elem.h"
#include "strnode.h"

extern "C" {
#include <float.h>
}

extern const char* psJointNames[];


/* Joint - begin */

class Joint 
: virtual public Elem, public ElemWithDofs, public InitialAssemblyElem {
 public:
   /* Tipi di Joint */
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
      PLANEDISP,
      PLANEDISPPIN,
      INPLANE,
      INPLANECONTACT,
      INLINE,
      ROD,
      DEFORMABLEHINGE,
      LINEARVELOCITY,
      ANGULARVELOCITY,
      LINEARACCELERATION,
      ANGULARACCELERATION,
      PRISMATIC,
      DRIVEHINGE,
      IMPOSEDKINEMATICS,

      BEAMSLIDER,
      
      MODAL,
      
      LASTJOINTTYPE
   };

 public: 
   class ErrGeneric {};
   
 private:
   Joint::Type JointT;
   
 public:
   Joint(unsigned int uL, Joint::Type T, const DofOwner* pD, flag fOut);
   virtual ~Joint(void);

   /* Derivate da Elem */
   
   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual Elem::Type GetElemType(void) const { 
      return Elem::JOINT; 
   };   
   
   /* Tipo di joint */
   virtual Joint::Type GetJointType(void) const {
      return JointT;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out << "  joint: " << GetLabel();
   };
   
   /* Output dei vincoli */
   virtual void Output(OutputHandler& OH) const;
   
   /* Output specifico dei vincoli */
   std::ostream& Output(std::ostream& out, const char* sJointName,
	           unsigned int uLabel, 
		   const Vec3& FLocal, const Vec3& MLocal,
		   const Vec3& FGlobal, const Vec3& MGlobal) const;

   /* Derivate da ElemWith Dofs */
   
   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& /* X */ ) const { 
      NO_OP;
   };
   
   virtual void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const { 
			    NO_OP;
   };

   /* per la lettura dei dati dell'elemento modale */

   friend Joint* ReadModal(DataManager* pDM,MBDynParser& HP, const DofOwner* pD0, unsigned int uLabel, const StructNode* pModalNode);

};

/* Joint - end */


/* Lettura Joints */
class DataManager;
class MBDynParser;

extern Elem* ReadJoint(DataManager* pDM,
		       MBDynParser& HP,
		       const DofOwner* pDO,
		       unsigned int uLabel);

#endif // JOINT_H
