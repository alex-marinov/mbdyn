/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2011-2017
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
 * Tommaso Solcia	<solcia@aero.polimi.it>
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
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#ifndef MEMBRANEEAS_H
#define MEMBRANEEAS_H

#include "myassert.h"
#include "except.h"

#include "strnode.h"
#include "elem.h"

#include "constltp.h"

/* da spostare in .cc */
#include "Rot.hh"
#include "joint.h"

#include "membrane.h"


// Membrane4EAS - begin

class Membrane4EAS
: virtual public Elem,
public Membrane
{
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
		IP_1_3 = 2,
		IP_2_1 = 3,
		IP_2_2 = 4,
		IP_2_3 = 5,
		IP_3_1 = 6,
		IP_3_2 = 7,
		IP_3_3 = 8,

		NUMIP = 4
	};
	
	static doublereal xi_i[NUMIP][2];
	static doublereal w_i[NUMIP];

	enum NodeName {
		NODE1 = 0,
		NODE2 = 1,
		NODE3 = 2,
		NODE4 = 3,

		NUMNODES = 4
	};

	static doublereal xi_n[NUMNODES][2];
	
	static doublereal xi_0[2];

	enum Deformations {
		STRAIN = 0,

		NUMDEFORM = 1
	};

protected:
	// Pointers to nodes
	const StructDispNode* pNode[NUMNODES];

	// nodal positions (0: initial; otherwise current)
	Vec3 xa_0[NUMNODES];
	Vec3 xa[NUMNODES];

	// Orientation matrix 
	//    .. in reference configuration
	Mat3x3 T_0_0;
	Mat3x3 T_0_i[NUMIP];
	//    .. in current configuration
	Mat3x3 T_0;
	Mat3x3 T_i[NUMIP];

	// linear deformation vectors
	//    .. in reference configuration
	Vec3 eps_tilde_0_i[NUMIP];
	//    .. in current configuration
	Vec3 eps_tilde_i[NUMIP];


	// Temporary data - TODO
#if 0
	Vec6 Az[NUMSEZ];
	Vec6 AzRef[NUMSEZ];
	Vec6 AzLoc[NUMSEZ];
	Vec6 DefLoc[NUMSEZ];
	Vec6 DefLocRef[NUMSEZ];
	Vec6 DefLocPrev[NUMSEZ];
#endif

protected:
	fmh S_alpha_beta_0;
	doublereal alpha_0;
	doublereal alpha_i[NUMIP];
	vfmh L_alpha_beta_i;

	vfmh B_overline_i;

	vfmh P_i;
	
	Vec3 y_i_1[NUMIP];
	Vec3 y_i_2[NUMIP];
	
	vh beta;
	vh epsilon_hat;
	vh epsilon;

#ifdef USE_CL_IN_MEMBRANE
	// Constitutive law handlers at integration points
	ConstitutiveLawOwner<vh, fmh>* pD[NUMIP];
#else // ! USE_CL_IN_MEMBRANE
	vvh FRef;
	bool bPreStress;
	vh PreStress;
#endif // ! USE_CL_IN_MEMBRANE

	// Reference constitutive law tangent matrices
	vfmh DRef;

	//stress
	vvh stress_i;

	// Is first residual
	bool bFirstRes;

	// for derived elements that add external contributions
	// to internal forces
	virtual void
	AddInternalForces(Vec6& /* AzLoc */ , unsigned int /* iSez */ ) {
		NO_OP;
	};

private:
	void UpdateNodalAndAveragePos();
	void ComputeInitialIptOrientation();

public:
	Membrane4EAS(unsigned int uL,
		const DofOwner* pDO,
		const StructDispNode* pN1, const StructDispNode* pN2,
		const StructDispNode* pN3, const StructDispNode* pN4,
#ifdef USE_CL_IN_MEMBRANE
		const ConstitutiveLaw<vh, fmh>** pDTmp, 
#else // ! USE_CL_IN_MEMBRANE
		const fmh& pDTmp,
		const vh& PreStress,
#endif // ! USE_CL_IN_MEMBRANE
		flag fOut);

	virtual ~Membrane4EAS(void);

	// Membrane type
	virtual Membrane::Type GetMembraneType(void) const {
		return Membrane::ELASTIC;
	};

	virtual unsigned int iGetNumDof(void) const { 
		//return 14;
		//return 13;
		return 7;
	};

	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i <= iGetNumDof());
		return DofOrder::ALGEBRAIC;
	};

	virtual DofOrder::Order GetEqType(unsigned int i) const {
		ASSERT(i >= 0 && i <= iGetNumDof());
		return DofOrder::ALGEBRAIC;
	};

	// Contribution to restart file
	virtual std::ostream& Restart(std::ostream& out) const;

#if 0
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// Inverse dynamics
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP);
#endif

	// Workspace size
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3*4 + iGetNumDof();
		*piNumCols = 3*4 + iGetNumDof();
	};

	// Initial settings
	void
	SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph = 0);

#if 0
	// Prepares reference parameters after prediction
	virtual void
	AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ );
#endif

	// Residual assembly
	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Inverse dynamics
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler&  XCurr ,
		const VectorHandler&  XPrimeCurr ,
		const VectorHandler&  XPrimePrimeCurr ,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);
#endif

	// Jacobian matrix assembly
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Matrix assembly for eigenvalues
	void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
#endif

	virtual void Output(OutputHandler& OH) const;

	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3*4;
		*piNumCols = 3*4;
	};

	virtual void SetInitialValue(VectorHandler& /* X */ ) {
		NO_OP;
	};

	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#if 0
	// Access to private data
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
#endif

	// Access to nodes
	virtual const StructDispNode* pGetNode(unsigned int i) const;

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

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

// Memebrane4EAS - end

#endif // MEMBRANEEAS_H

