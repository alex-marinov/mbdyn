/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifndef STRINTERP_H
#define STRINTERP_H

#include <nestedelem.h>
#include <node.h>

/* Geometry */
struct Geometry {
	unsigned uLabel;

	// kimenatics
	Vec3 X;
	Mat3x3 R;
	Vec3 V;
	Vec3 W;

	// optional kinematics
	Vec3 XPP;
	Vec3 WP;

	// optional forces
	Vec3 F;
	Vec3 M;
};

struct GeometryData {
	enum Flags {
		X			= 0x01U,
		R			= 0x02U,
		V			= 0x04U,
		W			= 0x08U,

		XPP			= 0x10U,
		WP			= 0x20U,

		ACCELERATIONS_MASK	= (XPP | WP ),

		F			= 0x40U,
		M			= 0x80U,

		FORCES_MASK		= (F | M)
	};
	unsigned uFlags;

	std::vector<Geometry> data;
};

/* StructOutputManip - begin */

class StructOutputManip {
public:
   	StructOutputManip(void);

	virtual ~StructOutputManip(void);

	virtual void Manipulate(const GeometryData& data) = 0;
};

/* StructOutputManip - end */

/* StructOutputEnd - begin */

class StructOutputEnd : virtual public Elem, public StructOutputManip {
public:
   	StructOutputEnd(unsigned uLabel, flag fOut);

	virtual ~StructOutputEnd(void);

	virtual Elem::Type GetElemType(void) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
};

/* StructOutputEnd - end */

/* StructOutputStart - begin */

class StructOutputStart : virtual public Elem, public NestedElem {
public:
   	StructOutputStart(const Elem *pE);

	virtual ~StructOutputStart(void);

	virtual void Manipulate(const GeometryData& data);
};

/* StructOutputStart - end */

/* StructOutput - begin */

class StructOutput : virtual public Elem, public NestedElem, public StructOutputManip {
public:
   	StructOutput(const Elem *pE);

	virtual ~StructOutput(void);

	virtual void Manipulate(const GeometryData& data);
};

/* StructOutput - end */

/* StructOutputCollectBase - begin */

class StructOutputCollectBase : virtual public Elem, public StructOutputStart {
protected:
	std::vector<const StructNode *> Nodes;
	GeometryData data;

	// Copiare da strext.cc:Send()
	virtual void Manipulate_int(void) = 0;

public:
   	StructOutputCollectBase(const Elem *pE,
		unsigned uFlags,
		const std::vector<const StructNode *>& nodes);

	virtual ~StructOutputCollectBase(void);

	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
};

/* StructOutputCollectBase - end */

/* StructOutputCollect - begin */

class StructOutputCollect : virtual public Elem, public StructOutputCollectBase {
protected:
	// Copiare da strext.cc:Send()
	virtual void Manipulate_int(void);

public:
   	StructOutputCollect(const Elem *pE,
		unsigned uFlags,
		const std::vector<const StructNode *>& nodes);

	virtual ~StructOutputCollect(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputCollect - end */

/* StructOutputCollectRelative - begin */

class StructOutputCollectRelative : virtual public Elem, public StructOutputCollectBase {
protected:
	const StructNode *pRefNode;

	// Copiare da strext.cc:Send()
	virtual void Manipulate_int(void);

public:
   	StructOutputCollectRelative(const Elem *pE,
		unsigned uFlags,
		const StructNode *pRefNode,
		const std::vector<const StructNode *>& nodes);

	virtual ~StructOutputCollectRelative(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputCollectRelative - end */

/* StructOutputInterpBase - begin */

class StructOutputInterpBase : virtual public Elem, public StructOutput {
protected:
	GeometryData data;

	virtual void InterpInit(void);
	virtual void InterpInit_int(void) = 0;
	virtual void Manipulate(const GeometryData& data);

public:
   	StructOutputInterpBase(const Elem *pE);

	virtual ~StructOutputInterpBase(void);
};

/* StructOutputInterpBase - end */

/* StructOutputInterp - begin */

class StructOutputInterp : virtual public Elem, public StructOutputInterpBase {
protected:
	virtual void InterpInit_int(void);

public:
   	StructOutputInterp(const Elem *pE);

	virtual ~StructOutputInterp(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputInterp - end */

/* StructOutputInterpOP2 - begin */

class StructOutputInterpOP2 : virtual public Elem, public StructOutputInterpBase {
protected:
	virtual void InterpInit_int(void);

public:
   	StructOutputInterpOP2(const Elem *pE);

	virtual ~StructOutputInterpOP2(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputInterp - end */

/* StructOutputWriteBase - begin */

class StructOutputWriteBase : virtual public Elem, public StructOutputEnd {
protected:
	const std::string outfilename;
	bool bNoClobberOut;

public:
   	StructOutputWriteBase(unsigned uLabel,
		const std::string& outfilename,
		bool bNoClobberOut,
		flag fOut);

	virtual ~StructOutputWriteBase(void);
};

/* StructOutputWriteBase - end */

/* StructOutputWrite - begin */

class StructOutputWrite : virtual public Elem, public StructOutputWriteBase {
protected:
	int iPrecision;

	// Scrive su file "nativo"
	virtual void Manipulate(const GeometryData& data);

public:
   	StructOutputWrite(unsigned uLabel,
		const std::string& outfilename,
		bool bNoClobberOut,
		int iPrecision,
		flag fOut);

	virtual ~StructOutputWrite(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputWrite - end */

/* StructOutputWriteNASTRAN - begin */

class StructOutputWriteNASTRAN : virtual public Elem, public StructOutputWriteBase {
protected:

	// Scrive su file bulk NASTRAN
	virtual void Manipulate(const GeometryData& data);

public:
   	StructOutputWriteNASTRAN(unsigned uLabel,
		const std::string& outfilename,
		bool bNoClobberOut,
		flag fOut);

	virtual ~StructOutputWriteNASTRAN(void);

	virtual std::ostream& Restart(std::ostream& out) const;
};

/* StructOutputWriteNASTRAN - end */

class DataManager;
class MBDynParser;

extern Elem *
ReadStructOutput(DataManager *pDM, MBDynParser& HP, unsigned int uLabel);

#endif // STRINTERP_H

