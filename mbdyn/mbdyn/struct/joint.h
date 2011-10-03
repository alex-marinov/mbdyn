/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include <cfloat>

#include "strnode.h"
#include "elem.h"
#include "gravity.h"
#include "hint_impl.h"

extern const char* psJointNames[];

/* Joint - begin */

class Joint
: virtual public Elem, public ElemGravityOwner,
	public ElemWithDofs, public InitialAssemblyElem
{
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
		UNIVERSALROTATION,
		UNIVERSALPIN,
		PLANEHINGE,
		PLANEROTATION,
		PLANEPIN,
		AXIALROTATION,
		PLANEDISP,
		PLANEDISPPIN,
		INPLANE,
		INPLANECONTACT,
		J_INLINE,
		ROD,
		DEFORMABLEHINGE,
		DEFORMABLEDISPJOINT,
		DEFORMABLEJOINT,
		DEFORMABLEAXIALJOINT,
		VISCOUSBODY,
		LINEARVELOCITY,
		ANGULARVELOCITY,
		LINEARACCELERATION,
		ANGULARACCELERATION,
		PRISMATIC,
		DRIVEHINGE,
		DRIVEDISP,
		DRIVEDISPPIN,
		IMPOSEDORIENTATION,
		IMPOSEDDISP,
		IMPOSEDDISPPIN,
		IMPOSEDKINEMATICS,
		BEAMSLIDER,
		BRAKE,
		GIMBAL,
		POINT_SURFACE_CONTACT,
		TOTALJOINT,
		TOTALPINJOINT,
		TOTALEQUATION,
		TOTALREACTION,

		MODAL,
		SCREWJOINT,

		LASTJOINTTYPE
	};

public:
	struct JointHint : public Hint {
		virtual ~JointHint(void) {};
	};

	template <int i>
	struct OffsetHint : public Joint::JointHint {};

	template <int i>
	struct HingeHint : public Joint::JointHint {};
	template <int i>
	struct PositionHingeHint : public Joint::HingeHint<i> {};
	template <int i>
	struct OrientationHingeHint : public Joint::HingeHint<i> {};

	template <class T>
	struct JointDriveHint : public Joint::JointHint {
		TplDriveHint<T> *pTDH;
		JointDriveHint(TplDriveHint<T> *pTDH) : pTDH(pTDH) {};
		~JointDriveHint(void) {
			if (pTDH) {
				delete pTDH;
			}
		};
	};
	template <class T>
	struct PositionDriveHint : public Joint::JointDriveHint<T> {
		PositionDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};
	template <class T>
	struct VelocityDriveHint : public Joint::JointDriveHint<T> {
		VelocityDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};
	template <class T>
	struct AccelerationDriveHint : public Joint::JointDriveHint<T> {
		AccelerationDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};
	template <class T>
	struct OrientationDriveHint : public Joint::JointDriveHint<T> {
		OrientationDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};
	template <class T>
	struct AngularVelocityDriveHint : public Joint::JointDriveHint<T> {
		AngularVelocityDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};
	template <class T>
	struct AngularAccelerationDriveHint : public Joint::JointDriveHint<T> {
		AngularAccelerationDriveHint(TplDriveHint<T> *pTDH) : Joint::JointDriveHint<T>(pTDH) {};
	};

	struct ReactionsHint : public Joint::JointHint {};
	struct ForcesHint : public Joint::ReactionsHint {};
	struct MomentsHint : public Joint::ReactionsHint {};

public:
 	class ErrGeneric : public MBDynErrBase {
  	public:
 		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

public:
	Joint(unsigned int uL, const DofOwner* pD, flag fOut);
	virtual ~Joint(void);

	/* Derivate da Elem */

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::JOINT;
	};

	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const = 0;

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "  joint: " << GetLabel();
	};

	/* Output specifico dei vincoli */
	std::ostream&
	Output(std::ostream& out, const char* sJointName,
		unsigned int uLabel,
		const Vec3& FLocal, const Vec3& MLocal,
		const Vec3& FGlobal, const Vec3& MGlobal) const;

	/* Derivate da ElemWith Dofs */

	/* Setta il valore iniziale delle proprie variabili */
	virtual void SetInitialValue(VectorHandler& /* X */ ) {
		NO_OP;
	};

	virtual void SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* Xp */ ,
		SimulationEntity::Hints *ph = 0)
	{
		NO_OP;
	};

	/* inverse dynamics update */
	virtual void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);
	
	/* per la lettura dei dati dell'elemento modale */

	friend Joint *
	ReadModal(DataManager* pDM,MBDynParser& HP, const DofOwner* pD0,
		unsigned int uLabel, const StructNode* pModalNode);
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
