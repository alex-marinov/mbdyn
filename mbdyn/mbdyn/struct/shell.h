/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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

/*
 * Authors: Pierangelo Masarati, Marco Morandini, Riccardo Vescovini
 */

/*
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#ifndef SHELL_H
#define SHELL_H

#include "myassert.h"
#include "except.h"

#include "strnode.h"
#include "elem.h"

#include "constltp.h"

extern const char* psShellNames[];

// Forward declaration
class DataManager;
class MBDynParser;

// Shell - begin

class Shell
: virtual public Elem,
// public ElemGravityOwner,
public InitialAssemblyElem
{
public:
	// Shell types
	enum Type {
		UNKNOWN = -1,
		ELASTIC = 0,
		VISCOELASTIC,

		LASTSHELLTYPE
	};

	Shell(unsigned uLabel, flag fOut);
	virtual ~Shell(void);
};

// Shell - end

// Shell4 - begin

class Shell4
: virtual public Elem,
public Shell
{
#if 0
protected:
 	static const unsigned int iNumPrivData =
		+3		//  0 ( 1 ->  3) - strain
		+3		//  3 ( 4 ->  6) - curvature
		+3		//  6 ( 7 ->  9) - force
		+3		//  9 (10 -> 12) - moment
		+3		// 12 (13 -> 15) - position
		+3		// 15 (16 -> 18) - orientation vector
		+3		// 18 (19 -> 21) - angular velocity
		+3		// 21 (22 -> 24) - strain rate
		+3		// 24 (25 -> 27) - curvature rate
	;

	static unsigned int iGetPrivDataIdx_int(const char *s,
		ConstLawType::Type type);
#endif

public:
	// numbered according to
	//
	//         ^
	// 4 o-----+-----o 3
	//   | 1_2 | 2_2 |
	// --+-----+-----+->
	//   | 1_1 | 2_1 |
	// 1 o-----+-----o 2
	//
	enum IntegrationPoint {
		IP_1_1 = 0,
		IP_1_2 = 1,
		IP_2_1 = 2,
		IP_2_2 = 3,

		NUMIP = 4
	};

	// numbered according to the side they are defined on
	enum ShearStrainEvaluationPoint {
		SSEP_1 = 0,
		SSEP_2 = 1,
		SSEP_3 = 2,
		SSEP_4 = 3,

		NUMSSEP = 4
	};

	enum NodeName {
		NODE1 = 0,
		NODE2 = 1,
		NODE3 = 2,
		NODE4 = 3,

		NUMNODES = 4
	};

	enum Deformations {
		STRAIN = 0,
		CURVAT = 1,

		NUMDEFORM = 2
	};

protected:
	// Pointers to nodes
	const StructNode* pNode[NUMNODES];

#if 0
	// Node offsets - TODO: offsets
	const Vec3 f[NUMNODES];
	Vec3 fRef[NUMNODES];
#endif

	const Mat3x3 RNode[NUMNODES];

	// Orientation matrix of the integration points
	Mat3x3 R[NUMIP];
	Mat3x3 RRef[NUMIP];
	Mat3x3 RPrev[NUMIP];

	// Orientation matrix of the shear strain evaluation points
	Mat3x3 R[NUMSSEP];
	Mat3x3 RRef[NUMSSEP];
	Mat3x3 RPrev[NUMSSEP];

	// Constitutive law handlers at integration points
	ConstitutiveLaw6DOwner* pD[NUMIP];
	// Reference constitutive law tangent matrices
	Mat6x6 DRef[NUMIP];

#if 0
	// Angular velocity of the sections - TODO: viscoelastic
	Vec3 Omega[NUMSEZ];
	Vec3 OmegaRef[NUMSEZ];
#endif

	// Temporary data - TODO
#if 0
	Vec6 Az[NUMSEZ];
	Vec6 AzRef[NUMSEZ];
	Vec6 AzLoc[NUMSEZ];
	Vec6 DefLoc[NUMSEZ];
	Vec6 DefLocRef[NUMSEZ];
	Vec6 DefLocPrev[NUMSEZ];
#endif

	Vec3 p[NUMIP];
	Vec3 g[NUMIP];
	Vec3 L0[NUMIP];
	Vec3 L[NUMIP];
	Vec3 LRef[NUMIP];

	doublereal dsdxi[NUMIP];

	// Is first residual
	bool bFirstRes;

	// Helpers
	static Vec3
	InterpState(const Vec3& v1,
		const Vec3& v2,
		const Vec3& v3,
		const Vec3& v4,
		enum Section Sec);
	Vec3
	InterpDeriv(const Vec3& v1,
		const Vec3& v2,
		const Vec3& v3,
		const Vec3& v4,
		enum Section Sec);

	// Compute matrices
	virtual void
	AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void
	AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// for derived elements that add external contributions
	// to internal forces
	virtual void
	AddInternalForces(Vec6& /* AzLoc */ , unsigned int /* iSez */ ) {
		NO_OP;
	};

	// Initialize shape functions derivatives; compute geometric curvature
	virtual void DsDxi(void);

	// Initial angular velocity of the sections
	// TODO: can be improved
	virtual void Omega0(void);

	// Internal restart
	virtual std::ostream& Restart_(std::ostream& out) const;

	// Data Initialization
	void Init(void);

public:
	Shell4(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const StructNode* pN3, const StructNode* pN4,
		const Vec3& f1, const Vec3& f2,
		const Vec3& f3, const Vec3& f4,
		const Mat3x3& R1, const Mat3x3& R2,
		const Mat3x3& R3, const Mat3x3& R4,
		const Mat3x3& r_I, const Mat3x3& rII,
		// TODO: add
		const ConstitutiveLaw6D* pD_I, const ConstitutiveLaw6D* pDII,
		// TODO: add
		flag fOut);

	virtual ~Shell4(void);

	// Shell type
	virtual Shell::Type GetShellType(void) const {
		return Shell::ELASTIC;
	};

	// Element type
	virtual Elem::Type GetElemType(void) const {
		return Elem::PLATE;
	};

	// Contribution to restart file
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// Inverse dynamics
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP);

	// Workspace size
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 24;
		*piNumCols = 24;
	};

	// Initial settings
	void
	SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph = 0);

	// Prepares reference parameters after prediction
	virtual void
	AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ );

	// Residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Inverse dynamics
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler&  XCurr ,
		const VectorHandler&  XPrimeCurr ,
		const VectorHandler&  XPrimePrimeCurr ,
		int iOrder = -1);

	// Jacobian matrix assembly
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Matrix assembly for eigenvalues
	void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;

	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 24;
		*piNumCols = 24;
	};

	virtual void SetInitialValue(VectorHandler& /* X */ ) {
		NO_OP;
	};

	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	// Access to private data
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	// Access to nodes
	virtual const StructNode* pGetNode(unsigned int i) const;

	/******** PER IL SOLUTORE PARALLELO *********/
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(NUMNODES);
		for (int i = 0; i < NUMNODES; i++) {
			connectedNodes[i] = pNode[i];
		}
	};
	/**************************************************/
};

// Shell4 - end

extern Elem*
ReadShell4(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif // SHELL_H

