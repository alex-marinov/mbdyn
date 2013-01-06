/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
 * Author: Matteo Fancello <matteo.fancello@gmail.com>
 * Nonsmooth dynamics element;
 * uses SICONOS <http://siconos.gforge.inria.fr/>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <iostream>
#include <cfloat>
#include "dataman.h"
#include "userelem.h"
#include "stlvh.h"

#include "mbdyn_siconos.h"


enum Constraining {
	// the interaction with MBDyn results
	// in imposing the position of the node associated
	// to the nonsmooth integration and in the force Lambda exchanged
	POSITION_ONLY,

	// the interaction with MBDyn results
	// in imposing the velocity of the node associated
	// to the nonsmooth integration and in the force Lambda exchanged
	VELOCITY_ONLY,

	// the interaction with MBDyn results
	// in imposing the position and velocity of the node associated
	// to the nonsmooth integration and in the force Lambda exchanged
	POSITION_AND_VELOCITY
};

struct plane {
	// Definizione del piano di vincolo:
	// nodo alla cui rotazione è ancorato, punto per cui passa, rotazione relativa al sys di rife del nodo
	const StructNode *AttNode;
	Vec3 Point;
	Mat3x3 RotRel;

	// H: relazione tra posizione e reazione, al passo k e k+1
	Mat3x3 H_k;
	Mat3x3 H_kp1;

	doublereal Gap;
	doublereal e_rest;
	doublereal mu;
};

struct NS_subsys {
	// Contains all the parameters of the nonsmooth system and its actual status

	// ---- NONSMOOTH SUBSYSTEM STATUS--------
	Vec3 Xns_k, Xns_kp1;		// NonSmooth node global position at time k, k+1
	Vec3 Vns_k, Vns_kp1;		// NonSmooth node global velocity at time k
	Vec3 Lambda_k, Lambda_kp1;	// Lagrange multiplier for node 1 at step k and k+1

	Vec3 GravityAcceleration;
	bool bGravity;
	doublereal hstep;

	// ---- NONSMOOTH SUBSYSTEM INPUT PARAMETERS -----
	doublereal mass_ns, radius_ns;
	int Np;					// number of planes which constitute a unilateral constraint each
	std::vector<plane> constr;

	doublereal theta_ts, gamma_pred;	// integration step
	solver_parameters solparam;			// solver parameters: tolerance, max iterations,

	Constraining constraint_type;
	int iterations;
	bool bVerbose;

	// ---- NONSMOOTH SUBSYSTEM DEBUG AND OUTPUT -----

	double pkp1; 	// norm of the reaion impulse
//	double bLCP, WNN, pkp1, Gap, Pkplus1, Un_k, Un_kp1, vfree;
	Vec3 mu;
	int ia;		// size of the vector of the active constraints

	Vec3 Piter;
	int niter;
	int niterLCPmax;
	Vec3 Poutput;
};

class ModuleNonsmoothNode
	: virtual public Elem, public UserDefinedElem
{
private:
	const DataManager *m_pDM;

	// Node associated to the nonsmooth subsystem
	const StructDispNode *m_pNode;

	// all the parameters of the nonsmooth subsystem and the updated variables
	NS_subsys  NS_data;

	// It allows to know if step changes (for the bVerbose output mainly)
	bool bStepToggle;

	// Constraint relaxation variable
	Vec3 Mu;

	// Input parameter to set a frictional problem
	bool bFrictional;

	// iteration number
	int iIter;

	// function timestepping a frictionless subsys
	void mbs_get_force(NS_subsys& );

	// function timestepping a frictional subsys
	void mbs_get_force_frictional(NS_subsys& );

public:
	ModuleNonsmoothNode(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleNonsmoothNode(void);

	/* inherited from SimulationEntity */
	virtual unsigned int iGetNumDof(void) const;

#if 0 // TODO
	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "", bool bInitial = false) const;
	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false, int i = -1) const;
	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "", bool bInitial = false) const;
	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false, int i = -1) const;
#endif

	virtual DofOrder::Order GetDofType(unsigned int) const;

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};


ModuleNonsmoothNode::ModuleNonsmoothNode( unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pDM(pDM),
m_pNode(0),
bStepToggle(false),
Mu(::Zero3),
bFrictional(false),
iIter(0)
{
	NS_data.bVerbose = false;
	NS_data.niter = 0;
	NS_data.Piter = ::Zero3;

	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module: 	nonsmooth-node\n"
"Author: 	Matteo Fancello <matteo.fancello@gmail.com>\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale\n"
"		Politecnico di Milano\n"
"		http://www.aero.polimi.it/\n"
"\n"
"	All rights reserved\n"
"\n"
"	Defines a unilateral constraint in form of a contact\n"
"	between a node and one or more planes, optionally with friction.\n"
"\n"
"	# in the \"control data\" block: increase by 1 the count of \"loadable elements\"\n"
"\n"
"	# in the \"elements\" block: add a user-defined element according to the syntax\n"
"	user defined :\n"
"		<element_label> , nonsmooth node ,\n"
"		[ frictional , { yes | no | (bool) <frictional> } , ]\n"
"		(StructDispNode) <NonsmoothNODELABEL> ,\n"
"		mass , (real) <mass> , radius , (real) <radius> ,\n"
"		planes , (int) <number_of_planes> ,\n"
"			<PlaneDefiningNODELABEL> ,\n"
"			position ,\n"
"			(Vec3) <relative_plane_position> ,\n"
"			rotation orientation ,\n"
"			(OrientationMatrix) <rel_rot_orientation_1> ,\n"
"			restitution , (real) <rest_coeff>\n"
"			[ , frictional coefficient , <mu> ]\n"
"			[ , <other planes> ]\n"
"		[ , constraint type , { position | velocity | both } ] # default: both\n"
"		[ , theta , <theta> ] [ , gamma , <gamma> ]\n"
"		[ , LCP solver , <solver> ]\n"
"		[ , tolerance , <tolerance> ][ , max iterations , <num_iter> ]\n"
"			# these options depend on LCP solver support, see\n"
"			# http://siconos.gforge.inria.fr/Numerics/LCPSolvers.html\n"
"		[ , limit iterations , <niterations> ]\n"
"		[ , limit LCP iterations , <niterations> ]\n"
"		[ , verbose , { yes | no | (bool) <verbose> } ] ;\n"
"\n"
"        <solver> ::= lexico lemke # the default\n"
"                   | rpgs\n"
"                   | qp\n"
"                   | cpg\n"
"                   | pgs\n"
"                   | psor\n"
"                   | nsqp\n"
"                   | latin\n"
"                   | latin_w\n"
"                   | newton_min\n"
"                   | newton_FB\n"
"\n"
"	OUTPUT:\n"
"		1:     element label\n"
"		2-4:   impulse on  NonSmooth node in global ref. frame\n"
"		5-7:   position of NonSmooth node in global ref. frame\n"
"		8-10:  velocity of Nonsmooth node in global ref. frame\n"
"		11-13: Lambda: constr. reaction between mbs node and NS node\n"
"		14:    p: norm of the impulse reaction in global ref. frame\n"
"		15:    Ia: number of active constraints in this step\n"
"\n"
"	if verbose, also:\n"
"		16-18: Mu: constraint relaxation factor\n"
"		19:    LCPsolver status: 0 -> convergence, 1 -> iter=maxiter, >1 -> fail\n"
"		       (only for some solvers)\n"
"		20:    LCPsolver resulting_error\n"
"		       (only for some solvers)\n"
			<< std::endl);


		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_pNode = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);

	if (HP.IsKeyWord("frictional")) {
		bFrictional = HP.GetYesNoOrBool(true);
	}


	if (!HP.IsKeyWord("mass")) {
		silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"mass\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	NS_data.mass_ns   = HP.GetReal();

	if (!HP.IsKeyWord("radius")) {
		silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"radius\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	NS_data.radius_ns = HP.GetReal();

	if (!HP.IsKeyWord("planes")) {
		silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"planes\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	NS_data.Np  = HP.GetInt();
	NS_data.constr.resize(NS_data.Np);

	for (int i = 0; i < NS_data.Np; i++) {
		NS_data.constr[i].AttNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		ReferenceFrame RFnode(NS_data.constr[i].AttNode);

		if (!HP.IsKeyWord("position")) {
			silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"position\" expected at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		NS_data.constr[i].Point = HP.GetPosRel(RFnode);

		if (!HP.IsKeyWord("rotation" "orientation")) {
			silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"rotation orientation \" expected at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		NS_data.constr[i].RotRel = HP.GetRotRel(RFnode);

		if (!HP.IsKeyWord("restitution")) {
			silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"restitution \" expected at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		NS_data.constr[i].e_rest  = HP.GetReal();

		if (bFrictional) {
			if (!HP.IsKeyWord("friction" "coefficient")) {
				silent_cerr("ModuleNonsmoothNode(" << uLabel << "): \"friction coefficient\" expected at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			NS_data.constr[i].mu = HP.GetReal();
		}
	}

	NS_data.constraint_type = POSITION_AND_VELOCITY;

	if (HP.IsKeyWord("constraint" "type")) {
		if (HP.IsKeyWord("position")) {
			NS_data.constraint_type = POSITION_ONLY;
		} else if (HP.IsKeyWord("velocity")) {
			NS_data.constraint_type = VELOCITY_ONLY;
		} else if (HP.IsKeyWord("both")) {
			NS_data.constraint_type = POSITION_AND_VELOCITY;
		} else {
			silent_cerr("ModuleNonsmoothNode: invalid constraint_type"
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("theta")) {
		NS_data.theta_ts = HP.GetReal();
	} else {
		NS_data.theta_ts = 0.5;
	}

	if (HP.IsKeyWord("gamma")) {
		NS_data.gamma_pred = HP.GetReal();
	} else {
		NS_data.gamma_pred = 1.;
	}

	// default: LEXICO_LEMKE
	NS_data.solparam.solver = LEXICO_LEMKE;
	if (HP.IsKeyWord("LCP" "solver")) {
		if (HP.IsKeyWord("lexico" "lemke")) {
			NS_data.solparam.solver = LEXICO_LEMKE;
		} else if (HP.IsKeyWord("rpgs")) {
			NS_data.solparam.solver = RPGS;
		} else if (HP.IsKeyWord("qp")) {
			NS_data.solparam.solver = QP;
		} else if (HP.IsKeyWord("cpg")) {
			NS_data.solparam.solver = CPG;
		} else if (HP.IsKeyWord("pgs")) {
			NS_data.solparam.solver = PGS;
		} else if (HP.IsKeyWord("psor")) {
			NS_data.solparam.solver = PSOR;
		} else if (HP.IsKeyWord("nsqp")) {
			NS_data.solparam.solver = NSQP;
		} else if (HP.IsKeyWord("latin")) {
			NS_data.solparam.solver = LATIN;
		} else if (HP.IsKeyWord("latin_w")) {
			NS_data.solparam.solver = LATIN_W;
		} else if (HP.IsKeyWord("newton_min")) {
			NS_data.solparam.solver = NEWTON_MIN;
		} else if (HP.IsKeyWord("newton_FB")) {
			NS_data.solparam.solver = NEWTON_FB;
		} else {
			silent_cerr("ModuleNonsmoothNode(" << GetLabel() << "): invalid lcp solver type"
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	NS_data.solparam.solvertol = 1.e-14;
	if (HP.IsKeyWord("tolerance")) {
		NS_data.solparam.solvertol = HP.GetReal();
	}

	NS_data.solparam.solveritermax = 100;
	if (HP.IsKeyWord("max" "iterations")) {
		NS_data.solparam.solveritermax = HP.GetInt();
	}

	NS_data.iterations = 0;
	if (HP.IsKeyWord("limit" "iterations")) {
		NS_data.iterations = HP.GetInt();
	}

	NS_data.niterLCPmax = 0;
	if (HP.IsKeyWord("limit" "LCP" "iterations")) {
		NS_data.niterLCPmax = HP.GetInt();
	}

	if (HP.IsKeyWord("verbose")) {
		NS_data.bVerbose = HP.GetYesNoOrBool(true);
	}


	//----------------------------------------------
	//		CHECK GAP
	//----------------------------------------------

	for (int i = 0; i < NS_data.Np; i++) {
		NS_data.constr[i].H_kp1 = NS_data.constr[i].AttNode->GetRCurr() * NS_data.constr[i].RotRel;
		Vec3 Plane_point = NS_data.constr[i].AttNode->GetXCurr() + NS_data.constr[i].AttNode->GetRCurr() * NS_data.constr[i].Point;
		double Gap = (m_pNode->GetXCurr() - Plane_point) * NS_data.constr[i].H_kp1.GetVec(3) - NS_data.radius_ns;

		if (NS_data.bVerbose) {
			silent_cout("ModuleNonsmoothNode(" << GetLabel() << ")" << std::endl
				<< "  initial nonsmooth node position " << m_pNode->GetXCurr() << std::endl
				<< "  Plane_point " << Plane_point << std::endl
				<< "  normal " << NS_data.constr[i].H_kp1.GetVec(3) << std::endl
				<< "  Gap " << Gap << std::endl
				<< "  RFnode " << NS_data.constr[i].AttNode->GetRCurr() << std::endl);
		}

		if (Gap <= 0.) {
			silent_cerr("ModuleNonsmoothNode(" << GetLabel() << "): "
				"unilateral constraint " << i << " violated before simulation; Gap=" << Gap << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	//Inizializzazione del sottosistema NonSmooth (struct) nel costruttore:
	// integration step
	NS_data.hstep = m_pDM->pGetDrvHdl()->dGetTimeStep();

	// Inizializzazione dello stato
	NS_data.Lambda_kp1 = Zero3;
	NS_data.Lambda_k = Zero3;
	NS_data.Xns_kp1 = m_pNode->GetXCurr();	//NS state at step k +1
	NS_data.Vns_kp1 = m_pNode->GetVCurr();	//NS velocity at step k	+ 1
	NS_data.Xns_k = NS_data.Xns_kp1;		//NS state at step k
	NS_data.Vns_k = NS_data.Vns_kp1;		//NS velocity at step k

	NS_data.solparam.info = 0;
	NS_data.solparam.resulting_error = 0;

	if (NS_data.bVerbose) {
		silent_cout("ModuleNonsmoothNode(" << GetLabel() << "): initial nonsmooth node position " << NS_data.Xns_kp1 << std::endl);
	}
}

ModuleNonsmoothNode::~ModuleNonsmoothNode(void)
{
	// destroy private data
	NO_OP;
}

unsigned int
ModuleNonsmoothNode::iGetNumDof(void) const
{
	switch (NS_data.constraint_type) {
	case POSITION_ONLY:
		return 3;

	case VELOCITY_ONLY:
		return 3;

	case POSITION_AND_VELOCITY:
		return 6;
	}

	// impossible
	ASSERT(0);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}


#if 0
std::ostream&
ModuleNonsmoothNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
}

void
ModuleNonsmoothNode::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
}

std::ostream&
ModuleNonsmoothNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
}

void
ModuleNonsmoothNode::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
}
#endif

DofOrder::Order
ModuleNonsmoothNode::GetDofType(unsigned int) const
{
	return DofOrder::ALGEBRAIC;
}

void
ModuleNonsmoothNode::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	switch (NS_data.constraint_type) {
	case POSITION_ONLY:
		*piNumRows = 6;
		*piNumCols = 6;
		break;

	case VELOCITY_ONLY:
		*piNumRows = 6;
		*piNumCols = 6;
		break;

	case POSITION_AND_VELOCITY:
		*piNumRows = 9;
		*piNumCols = 9;
		break;

	default:
		ASSERT(0);
	}
}

void
ModuleNonsmoothNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	// called right after prediction, before any iteration within the time step
	// estrai Lambda1 e Lambda2 da XCurr e "mettili via" per il non-smooth
	integer iFirstIndex = iGetFirstIndex();
	NS_data.Lambda_kp1 = Vec3(X, iFirstIndex + 1);
	// Inizializzazione dello stato
	NS_data.Xns_kp1 = m_pNode->GetXCurr();	//NS state at step k +1
	NS_data.Vns_kp1 = m_pNode->GetVCurr();	//NS velocity at step k	+ 1
}

void
ModuleNonsmoothNode::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	// called after a step converged
	// store lambda at step k
	NS_data.Lambda_k = NS_data.Lambda_kp1;
	bStepToggle = false;
	iIter = 0;
	NS_data.niter = 0;
	NS_data.Piter = Zero3;
}

VariableSubMatrixHandler&
ModuleNonsmoothNode::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	integer iNode1FirstPositionIndex = m_pNode->iGetFirstPositionIndex();
	integer iNode1FirstMomentumIndex = m_pNode->iGetFirstMomentumIndex();
	integer iFirstIndex = iGetFirstIndex();

	switch (NS_data.constraint_type) {
	case POSITION_ONLY:
		WM.ResizeReset(6, 1);
		// add	  + I Dl = -lambda      to existing momentum equation
		WM.PutItem(1, iNode1FirstMomentumIndex + 1, iFirstIndex + 1, 1.);
		WM.PutItem(2, iNode1FirstMomentumIndex + 2, iFirstIndex + 2, 1.);
		WM.PutItem(3, iNode1FirstMomentumIndex + 3, iFirstIndex + 3, 1.);
		// constraint equation on position: I dq  = q-q
		WM.PutItem(4, iFirstIndex + 1, iNode1FirstPositionIndex + 1, 1.);
		WM.PutItem(5, iFirstIndex + 2, iNode1FirstPositionIndex + 2, 1.);
		WM.PutItem(6, iFirstIndex + 3, iNode1FirstPositionIndex + 3, 1.);
		break;

	case VELOCITY_ONLY:
		WM.ResizeReset(6, 1);
		// add	  + I Dl = -lambda      to existing momentum equation
		WM.PutItem(1, iNode1FirstMomentumIndex + 1, iFirstIndex + 1, 1.);
		WM.PutItem(2, iNode1FirstMomentumIndex + 2, iFirstIndex + 2, 1.);
		WM.PutItem(3, iNode1FirstMomentumIndex + 3, iFirstIndex + 3, 1.);
		// constraint equation on velocity:		I Dv = v-v
		WM.PutItem(4, iFirstIndex + 1, iNode1FirstPositionIndex + 1, 1.);
		WM.PutItem(5, iFirstIndex + 2, iNode1FirstPositionIndex + 2, 1.);
		WM.PutItem(6, iFirstIndex + 3, iNode1FirstPositionIndex + 3, 1.);
		break;

	case POSITION_AND_VELOCITY:
		WM.ResizeReset(12, 1);
		// add	  + I Dl = -lambda      to existing momentum equation
		WM.PutItem(1, iNode1FirstMomentumIndex + 1, iFirstIndex + 1, 1.);
		WM.PutItem(2, iNode1FirstMomentumIndex + 2, iFirstIndex + 2, 1.);
		WM.PutItem(3, iNode1FirstMomentumIndex + 3, iFirstIndex + 3, 1.);
		// constraint equation on position, relaxed: I dq + I dmu = q-q -mu
		WM.PutItem(4, iFirstIndex + 1, iNode1FirstPositionIndex + 1, dCoef);
		WM.PutItem(5, iFirstIndex + 2, iNode1FirstPositionIndex + 2, dCoef);
		WM.PutItem(6, iFirstIndex + 3, iNode1FirstPositionIndex + 3, dCoef);
		WM.PutItem(7, iFirstIndex + 1, iFirstIndex + 4, 1.);
		WM.PutItem(8, iFirstIndex + 2, iFirstIndex + 5,	1.);
		WM.PutItem(9, iFirstIndex + 3, iFirstIndex + 6, 1.);
		// constraint equation on velocity:		I Dv = v-v
		WM.PutItem(10, iFirstIndex + 4, iNode1FirstPositionIndex + 1, 1.);
		WM.PutItem(11, iFirstIndex + 5, iNode1FirstPositionIndex + 2, 1.);
		WM.PutItem(12, iFirstIndex + 6, iNode1FirstPositionIndex + 3, 1.);
		break;

	default:
		// impossible
		ASSERT(0);
	}

	return WorkMat;
}

SubVectorHandler&
ModuleNonsmoothNode::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	integer iNode1FirstIndex = m_pNode->iGetFirstMomentumIndex();
	integer iFirstIndex = iGetFirstIndex();
	// 	Variables needed for the iteration:
	//	timestep value, gravity acceleration (if it is set),
	//	position and velocity of the node coincident with the nonsmooth node,
	//	forces applied on that node by the rest of the model.
	NS_data.hstep = m_pDM->pGetDrvHdl()->dGetTimeStep();
	const Vec3& Xn = m_pNode->GetXCurr();
	NS_data.bGravity = GravityOwner::bGetGravity(Xn, NS_data.GravityAcceleration);

	NS_data.Xns_k = m_pNode->GetXPrev();
	NS_data.Vns_k = m_pNode->GetVPrev();

	NS_data.Xns_kp1 = m_pNode->GetXCurr();
//	std::cout << "NS_data.Xns_kp1="<<NS_data.Xns_kp1;
	NS_data.Vns_kp1 = m_pNode->GetVCurr();
	NS_data.Lambda_kp1 = Vec3(XCurr, iFirstIndex + 1);

	// integra un'iterazione del systema NonSmooth e restituisci il nuovo Lambda
	// e Xns1 e Vns1 all'interno di NS_data

	if (((NS_data.iterations > 0) && (iIter < NS_data.iterations)) || (NS_data.iterations == 0)) {
		if (bFrictional) {
			mbs_get_force_frictional(NS_data);
		} else {
			mbs_get_force(NS_data);
		}
		iIter++;
	}

	Vec3 DX1(NS_data.Xns_kp1 - m_pNode->GetXCurr());
	Vec3 DV1(NS_data.Vns_kp1 - m_pNode->GetVCurr());

	if (NS_data.bVerbose) {
		silent_cout("ModuleNonsmoothNode(" << GetLabel() << ") step=" << m_pDM->pGetDrvHdl()->iGetStep() << std::endl
			<< "  NS_data.Xns_kp1={" << NS_data.Xns_kp1 << "} DX={" << DX1 << "}" << std::endl
			<< "  NS_data.Vns_kp1={" << NS_data.Vns_kp1 << "} DV={" << DV1 << "}" << std::endl
			<< "  NS_data.Lambda_kp1={" << NS_data.Lambda_kp1 << "}" << std::endl);
	}

	switch (NS_data.constraint_type) {
	case POSITION_ONLY:
		WorkVec.ResizeReset(6);
		WorkVec.PutItem(1, iNode1FirstIndex + 1, - NS_data.Lambda_kp1(1));
		WorkVec.PutItem(2, iNode1FirstIndex + 2, - NS_data.Lambda_kp1(2));
		WorkVec.PutItem(3, iNode1FirstIndex + 3, - NS_data.Lambda_kp1(3));
		WorkVec.PutItem(4, iFirstIndex + 1, DX1(1)/dCoef);
		WorkVec.PutItem(5, iFirstIndex + 2, DX1(2)/dCoef);
		WorkVec.PutItem(6, iFirstIndex + 3, DX1(3)/dCoef);
		break;

	case VELOCITY_ONLY:
		WorkVec.ResizeReset(6);
		WorkVec.PutItem(1, iNode1FirstIndex + 1, - NS_data.Lambda_kp1(1));
		WorkVec.PutItem(2, iNode1FirstIndex + 2, - NS_data.Lambda_kp1(2));
		WorkVec.PutItem(3, iNode1FirstIndex + 3, - NS_data.Lambda_kp1(3));
		WorkVec.PutItem(4, iFirstIndex + 1, DV1(1));
		WorkVec.PutItem(5, iFirstIndex + 2, DV1(2));
		WorkVec.PutItem(6, iFirstIndex + 3, DV1(3));
		break;

	case POSITION_AND_VELOCITY:
		WorkVec.ResizeReset(9);
		Mu = Vec3(XCurr, iFirstIndex + 4);
		// residuo per le equazioni di congruenza del nodo 1
		DX1 -= Mu;
		NS_data.mu = Mu;
		WorkVec.PutItem(1, iNode1FirstIndex + 1, - NS_data.Lambda_kp1(1));
		WorkVec.PutItem(2, iNode1FirstIndex + 2, - NS_data.Lambda_kp1(2));
		WorkVec.PutItem(3, iNode1FirstIndex + 3, - NS_data.Lambda_kp1(3));
		WorkVec.PutItem(4, iFirstIndex + 1, DX1(1));
		WorkVec.PutItem(5, iFirstIndex + 2, DX1(2));
		WorkVec.PutItem(6, iFirstIndex + 3, DX1(3));
		WorkVec.PutItem(7, iFirstIndex + 4, DV1(1));
		WorkVec.PutItem(8, iFirstIndex + 5, DV1(2));
		WorkVec.PutItem(9, iFirstIndex + 6, DV1(3));
		break;
	}

	if (NS_data.bVerbose) {
		silent_cout("ModuleNonsmoothNode(" << GetLabel() << ") workvec:" << std::endl << WorkVec << std::endl);
	}

	return WorkVec;
}

unsigned int
ModuleNonsmoothNode::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleNonsmoothNode::iGetNumConnectedNodes(void) const
{
	// silent_cerr("start iGetNumConnectedNodes " << std::endl);
	return 1 + NS_data.constr.size();
}

void
ModuleNonsmoothNode::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(iGetNumConnectedNodes());
	connectedNodes[0] = m_pNode;
	unsigned i = 0;
	for (; i < NS_data.constr.size(); i++) {
		connectedNodes[i + 1] = NS_data.constr[i].AttNode;
	}
}

void
ModuleNonsmoothNode::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleNonsmoothNode::Restart(std::ostream& out) const
{
	return out << "# ModuleNonsmoothNode: not implemented" << std::endl;
}

unsigned int
ModuleNonsmoothNode::iGetInitialNumDof(void) const
{
	return 0;
}

void
ModuleNonsmoothNode::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleNonsmoothNode::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
ModuleNonsmoothNode::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);
	WorkVec.ResizeReset(0);
	return WorkVec;
}

void
ModuleNonsmoothNode::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();
		out << std::setw(8) << GetLabel()
			<< " " << NS_data.Poutput
			<< " " << NS_data.Xns_kp1
			<< " " << NS_data.Vns_kp1
			<< " " << NS_data.Lambda_kp1
			<< " " << NS_data.pkp1
			<< " " << NS_data.ia;

		if (NS_data.bVerbose) {
			out
				<< " " << NS_data.mu
				<< " " << NS_data.solparam.info
				<< " " << NS_data.solparam.resulting_error;
		}

		out << std::endl;
	}

	// TODO: use NetCDF?
}

void
ModuleNonsmoothNode::mbs_get_force(NS_subsys& NS)
{
	// Some Output initialization
	Vec3 p_kp1 = Zero3;
	NS.pkp1 = 0.;
	NS.ia = 0;

//	Mat3x3 M_ns = Eye3 * NS.mass_ns;
//	Mat3x3 K_ns = Zero3x3;
//	Mat3x3 C_ns = Zero3x3;
//	Mat3x3 Mhat = M_ns + C_ns * (NS.hstep * NS.theta_ts) + K_ns * (NS.hstep * NS.hstep * NS.theta_ts * NS.theta_ts);
//	Mat3x3 invMhat = Mhat.Inv();
	Mat3x3 M_ns = Eye3 * NS.mass_ns;
	Mat3x3 invMhat = Eye3 * (1./NS.mass_ns);

	// Determine the predicted NS state at the step k+1, to verify if the gap is closed
	Vec3 x_NS_predicted = NS.Xns_k + NS.Vns_k * (NS.hstep *  NS.gamma_pred);

	// Create index of active constraints: Ia, through evaluation of the gap distance.
	std::vector<int> Ia;
	for (int i = 0; i < NS.Np; i++) {
		// Store the orientation matrix H of the reference frame with Z aligned with contact direction
		NS.constr[i].H_kp1 = NS.constr[i].AttNode->GetRCurr() * NS.constr[i].RotRel;
		NS.constr[i].H_k = NS.constr[i].AttNode->GetRPrev() * NS.constr[i].RotRel;

		Vec3 Plane_point = NS.constr[i].AttNode->GetXCurr() + NS.constr[i].AttNode->GetRCurr() * NS.constr[i].Point;
		NS.constr[i].Gap = (x_NS_predicted - Plane_point) * NS.constr[i].H_kp1.GetVec(3) - NS.radius_ns;

		// Valuto gap(xpredicted) e creo indice vincoli attivi
		if (NS.constr[i].Gap <= 0) {
			int indice = i;
			Ia.push_back(indice);
		}
	}

	// External forces: gravity and the reaction passed from the rest of the system integrated by MBDYN
	Vec3 Fg = NS.bGravity ? NS.GravityAcceleration * NS.mass_ns : Zero3;
	Vec3 F_k = NS.Lambda_k + Fg;
	Vec3 F_kp1 = NS.Lambda_kp1 + Fg;

	// Velocity of the non-smooth subsystem calculated without the contribute of contact impulse
	// Vec3 vfree_kp1= NS.Vns_k + invMhat * ( - C_ns*NS.Vns_k*NS.hstep - K_ns*NS.Xns_k*NS.hstep - K_ns*NS.Vns_k*(NS.hstep*NS.hstep*NS.theta_ts) + (F_kp1*NS.theta_ts + F_k* (1.-NS.theta_ts) ) *NS.hstep);
	Vec3 vfree_kp1 = NS.Vns_k + invMhat * (F_kp1*NS.theta_ts + F_k* (1. - NS.theta_ts)) * NS.hstep;

	Vec3 VfreeRel = vfree_kp1 - NS.constr[0].AttNode->GetVCurr();		// VfreeMinusVplane
	Vec3 VnskRel = NS.Vns_k - NS.constr[0].AttNode->GetVCurr();

	// Composing and solving the Linear Complementarity System to determine the contact impulse p_kp1
	if (Ia.size() > 0) {
		Mat3xN Hfull_kp1(Ia.size());
		Mat3xN Hfull_k(Ia.size());

		for (unsigned i = 0; i < Ia.size(); i++) {
			Hfull_kp1.PutVec(i + 1, NS.constr[Ia[i]].H_kp1.GetVec(3));
			Hfull_k.PutVec(i + 1, NS.constr[Ia[i]].H_k.GetVec(3));
		}

		MatNx3 Hfull_kp1T(Ia.size());
		Hfull_kp1T.Transpose(Hfull_kp1);
		MatNx3 Hfull_kT(Ia.size());
		Hfull_kT.Transpose(Hfull_k);
		MatNx3 HTW(Ia.size());
		HTW.RightMult(Hfull_kp1T, invMhat);
		MatNxN W_delassus(Ia.size(), 0.);
		W_delassus.Mult(HTW, Hfull_kp1);

		// Attenzione, matrici riempite alla FORTRAN! Passo a Siconos le matrici come le vuole lui.
		double W_passed[Ia.size()*Ia.size()];
		int indice = 0;
		for (unsigned i = 1; i < Ia.size() + 1; i++) {
			for (unsigned j = 1; j < Ia.size() + 1; j++) {
				W_passed[indice] = W_delassus(j, i);
				indice++;
			}
		}

		double bLCP_new[Ia.size()];
		for (unsigned i = 0; i < Ia.size(); i++) {
			bLCP_new[i] = Hfull_kp1T.GetVec(i + 1) * VfreeRel + Hfull_kT.GetVec(i + 1) * VnskRel * NS.constr[Ia[i]].e_rest;
		}

		VecN P(Ia.size(), 0.);
		VecN U(Ia.size(), 0.);

		if ((NS.niter < NS.niterLCPmax) || (NS.niterLCPmax == 0)) {
			// Solving the LCP with one of Siconos solvers
			mbdyn_siconos_LCP_call(Ia.size(), W_passed, bLCP_new, P.pGetVec(), U.pGetVec(), NS.solparam);

			// Reazioni nel sistema di riferimento globale
			p_kp1 = Hfull_kp1 * P;

			NS.Piter = p_kp1;

			NS.niter++;
		}
	}

	// Updating: state actualization with the contribution of the contact impulse
	NS.Vns_kp1 = vfree_kp1 + invMhat * NS.Piter;
	NS.Xns_kp1  = NS.Xns_k + (NS.Vns_kp1 * NS.theta_ts + NS.Vns_k* (1. - NS.theta_ts)) * NS.hstep;
	// Reaction force to be exchanged with the rest of the system in MBDYN
	NS.Lambda_kp1 = - (Fg + p_kp1 / NS.hstep - M_ns * (NS.Vns_kp1 - NS.Vns_k) / NS.hstep + NS.Lambda_k*(1. - NS.theta_ts)) / NS.theta_ts;

	// Store output
	NS.ia = Ia.size();
	NS.pkp1 = p_kp1.Norm();
	NS.Poutput = p_kp1;
	// p_kp1 = NS.Piter;

	// Indicatore di nuovo step, per evitare troppo output con bVerbose.
	bStepToggle = true;
}

// Creo matrici fullmh IPa e IRa contenenti i vettori delle direzioni delle facets del friction cone outer discretized
static const doublereal dIPa[2][2] = {
	{ 1., 0. },
	{ 0., 1. }
};

static const doublereal dIRa[2][14] = {
	{  0.92388,  0.70711,  0.38268, -0.38268, -0.70711, -0.92388, -1.00000, -0.92388, -0.70711, -0.38268, -0.00000,  0.38268,  0.70711,  0.92388 },
	{  0.38268,  0.70711,  0.92388,  0.92388,  0.70711,  0.38268,  0.00000, -0.38268, -0.70711, -0.92388, -1.00000, -0.92388, -0.70711, -0.38268 }
};

static FullMatrixHandler IRa;
static FullMatrixHandler IRat;
static FullMatrixHandler IPa;
static bool bIRa(false);

void
ModuleNonsmoothNode::mbs_get_force_frictional(NS_subsys& NS)
{
	// FIXME
	// vanno messe nel costruttore, omega changeable ?-----------------------------

	// number of facets of the discretized friction cone
	int omegan = 16;

	if (!bIRa) {
		bIRa = true;

		IRa.Resize(2, 14);
		for (int r = 0; r < 2; r++) {
			for (int c = 0; c < 14; c++) {
				IRa(r + 1, c + 1) = dIRa[r][c];
			}
		}

		IRat.Resize(14, 2);
		for (int r = 0; r < 14; r++) {
			for (int c = 0; c < 2; c++) {
				IRat(r + 1, c + 1) = dIRa[c][r];
			}
		}

		IPa.Resize(2, 2);
		for (int r = 0; r < 2; r++) {
			for (int c = 0; c < 2; c++) {
				IPa(r + 1, c + 1) = dIPa[r][c];
			}
		}
	}

//------------------------------------------------------------------------------

//	Mat3x3 M_ns = Eye3 * NS.mass_ns;
//	Mat3x3 K_ns = Zero3x3;
//	Mat3x3 C_ns = Zero3x3;
//	Mat3x3 Mhat = M_ns + C_ns * (NS.hstep * NS.theta_ts) + K_ns * (NS.hstep * NS.hstep * NS.theta_ts * NS.theta_ts);
//	Mat3x3 invMhat = Mhat.Inv();
	Mat3x3 M_ns = Eye3 * NS.mass_ns;
	Mat3x3 invMhat = Eye3 * (1./NS.mass_ns);

	Vec3 p_kp1 = Zero3;
	NS.Poutput = Zero3;
	NS.ia = 0;

	// Determine the predicted NS state at the step k+1, to verify if the gap is closed
	Vec3 x_NS_predicted = NS.Xns_k + NS.Vns_k * NS.hstep * NS.gamma_pred;

	// Create index of active constraints: Ia, through evaluation of the gap distance.
	std::vector<int> Ia;
	for (int i = 0; i < NS.Np; i++) {
		// Store the orientation matrix H of the reference frame with Z aligned with contact direction
		NS.constr[i].H_kp1 = NS.constr[i].AttNode->GetRCurr() * NS.constr[i].RotRel;
		NS.constr[i].H_k = NS.constr[i].AttNode->GetRPrev() * NS.constr[i].RotRel;

		Vec3 Plane_point = NS.constr[i].AttNode->GetXCurr() + NS.constr[i].AttNode->GetRCurr() * NS.constr[i].Point;
		NS.constr[i].Gap = (x_NS_predicted -  Plane_point) * NS.constr[i].H_kp1.GetVec(3) - NS.radius_ns;

		// Valuto gap(xpredicted) e creo indice vincoli attivi
		if (NS.constr[i].Gap <= 0) {
			int indice = i;
			Ia.push_back(indice);
		}
	}

	// External forces: gravity and the reaction passed from the rest of the system integrated by MBDYN
	Vec3 Fg = NS.bGravity ? NS.GravityAcceleration * NS.mass_ns : Zero3;
	Vec3 F_k = NS.Lambda_k + Fg;
	Vec3 F_kp1 = NS.Lambda_kp1 + Fg;

	// Velocity of the non-smooth subsystem calculated without the contribute of contact impulse
	// Vec3 vfree_kp1= NS.Vns_k + invMhat * ( - C_ns*NS.Vns_k*NS.hstep - K_ns*NS.Xns_k*NS.hstep - K_ns*NS.Vns_k*(NS.hstep*NS.hstep*NS.theta_ts) + (F_kp1*NS.theta_ts + F_k* (1.-NS.theta_ts) ) *NS.hstep);
	Vec3 vfree_kp1 = NS.Vns_k + invMhat * (F_kp1*NS.theta_ts + F_k * (1. - NS.theta_ts)) * NS.hstep;
	Vec3 VfreeRel = vfree_kp1 - NS.constr[0].AttNode->GetVCurr();	// VfreeMinusVplane

	if (Ia.size() > 0) {
		int ac = Ia.size();			// active constraints !!!
		FullMatrixHandler H_NN(3, ac);
		FullMatrixHandler H_TT(3, 2*ac);
		FullMatrixHandler H_NNt(ac, 3);
		FullMatrixHandler H_TTt(2*ac, 3);

		FullMatrixHandler Htot(3, 3*ac);

		FullMatrixHandler I_P(2*ac, 2*ac);
		FullMatrixHandler I_R(2*ac, (omegan - 2)*ac);
		FullMatrixHandler I_Rt((omegan - 2)*ac, 2*ac);
		FullMatrixHandler mu_P(2*ac, ac);
		FullMatrixHandler mu_R((omegan - 2)*ac, ac);
		FullMatrixHandler en(1, ac);

		for (unsigned i = 0; i < Ia.size(); i++) {
			en(1, i + 1) = NS.constr[Ia[i]].e_rest;

			H_NN.Put(1, i + 1, NS.constr[Ia[i]].H_kp1.GetVec(3));
			H_TT.Put(1, 2*i + 1, NS.constr[Ia[i]].H_kp1.GetVec(1));
			H_TT.Put(1, 2*i + 2, NS.constr[Ia[i]].H_kp1.GetVec(2));
			H_NNt.PutT(i + 1, 1, NS.constr[Ia[i]].H_kp1.GetVec(3));
			H_TTt.PutT(2*i + 1, 1, NS.constr[Ia[i]].H_kp1.GetVec(1));
			H_TTt.PutT(2*i + 2, 1, NS.constr[Ia[i]].H_kp1.GetVec(2));

			I_P.Put(2*(i + 1 - 1) + 1, 2*(i + 1 - 1) + 1, IPa);
			I_R.Put(2*(i + 1 - 1) + 1, (omegan - 2)*(i + 1 - 1) + 1, IRa);
			I_Rt.Put((omegan - 2)*(i + 1 - 1) + 1, 2*(i + 1 - 1) + 1, IRat);

			FullMatrixHandler mu_Piw(2, 1);
			mu_Piw(1, 1) = NS.constr[Ia[i]].mu;
			mu_Piw(2, 1) = NS.constr[Ia[i]].mu;

			FullMatrixHandler mu_Riw(omegan - 2, 1);

			for (int in = 1; in < omegan - 2 + 1; in++) {
				mu_Riw(in, 1) = NS.constr[Ia[i]].mu;
			}

			mu_P.Put(2*(i + 1 - 1) + 1, i + 1, mu_Piw);
			mu_R.Put((omegan - 2)*(i + 1 - 1) + 1, i + 1, mu_Riw);
		}

		Htot.Put(1, 1, H_TT);
		Htot.Put(1, 2*ac + 1, H_NN);

		FullMatrixHandler v_free(3, 1);
		v_free.Put(1, 1, VfreeRel);
		FullMatrixHandler U_Nfree(1*ac, 1);
		FullMatrixHandler U_Tfree(2*ac, 1);
		U_Nfree.MatMul(H_NNt, v_free);
		U_Tfree.MatMul(H_TTt, v_free);

		FullMatrixHandler iMhat(3, 3);
		iMhat.Put(1, 1, invMhat);

		FullMatrixHandler WNNl(ac, 3);
		WNNl.MatMul(H_NNt, iMhat);
		FullMatrixHandler WNN(ac, ac);
		WNN.MatMul(WNNl, H_NN);

		FullMatrixHandler WTTl(2*ac, 3);
		WTTl.MatMul(H_TTt, iMhat);
		FullMatrixHandler WTT(2*ac, 2*ac);
		WTT.MatMul(WTTl, H_TT);

		FullMatrixHandler WNTl(ac, 3);
		WNTl.MatMul(H_NNt, iMhat);
		FullMatrixHandler WNT(ac, 2*ac);
		WNT.MatMul(WNTl, H_TT);

		FullMatrixHandler WTNl(2*ac, 3);
		WTNl.MatMul(H_TTt, iMhat);
		FullMatrixHandler WTN(2*ac, ac);
		WTN.MatMul(WTNl, H_NN);

		// FIXME: adesso è facile perchè eye(2*ac)
		FullMatrixHandler IPmT(2*ac, 2*ac);
		for (int ipmti = 1; ipmti < 2*ac + 1; ipmti++) {
			IPmT(ipmti, ipmti) = 1;
		}

		// MLCP construction
		int ProbDim = (3 + omegan - 2)*ac;
		FullMatrixHandler Mlcp(ProbDim, ProbDim);
		FullMatrixHandler mlcp11(ac, ac);
		FullMatrixHandler mlcp12(ac, 2*ac);
		FullMatrixHandler mlcp21(2*ac, ac);
		FullMatrixHandler mlcp22(2*ac, 2*ac);
		FullMatrixHandler mlcp23(2*ac, (omegan - 2)*ac);
		FullMatrixHandler mlcp31((omegan - 2)*ac, ac);
		FullMatrixHandler mlcp32((omegan - 2)*ac, 2*ac);

		FullMatrixHandler WNT_IPmT(ac, 2*ac);
		WNT_IPmT.MatMul(WNT, IPmT);
		FullMatrixHandler WNT_IPmT_muP(ac, ac);
		WNT_IPmT_muP.MatMul(WNT_IPmT, mu_P);

		mlcp11.Put(1, 1, WNN);
		mlcp11.Add(1, 1, WNT_IPmT_muP);
		mlcp12.Sub(1, 1, WNT_IPmT);

		FullMatrixHandler WTT_IPmT(2*ac, 2*ac);
		WTT_IPmT.MatMul(WTT, IPmT);
		FullMatrixHandler WTT_IPmT_muP(2*ac, ac);
		WTT_IPmT_muP.MatMul(WTT_IPmT, mu_P);
		WTT_IPmT_muP.Add(1, 1, WTN);
		mlcp21.Sub(1, 1, WTT_IPmT_muP);

		mlcp22.Put(1, 1, WTT);

		mlcp23.Sub(1, 1, I_R);

		FullMatrixHandler I_Rt_mu_P((omegan - 2)*ac, ac);
		I_Rt_mu_P.MatMul(I_Rt, mu_P);

		mlcp31.Add(1, 1, mu_R);
		mlcp31.Sub(1, 1, I_Rt_mu_P);
		mlcp32.Add(1, 1, I_Rt);

		Mlcp.Put(1, 1, mlcp11);
		Mlcp.Put(1, ac + 1, mlcp12);
		Mlcp.Put(ac + 1, 1, mlcp21);
		Mlcp.Put(ac + 1, ac + 1, mlcp22);
		Mlcp.Put(ac + 1, 3*ac + 1, mlcp23);
		Mlcp.Put(3*ac + 1, 1, mlcp31);
		Mlcp.Put(3*ac + 1, ac + 1, mlcp32);

		// FIXME: use STLVectorHandler?
		FullMatrixHandler qlcp(ProbDim, 1);
		FullMatrixHandler qlcp1(ac, 1);

		STLVectorHandler VnskMinusVplane(3);
		VnskMinusVplane.Put(1, NS.Vns_k - NS.constr[0].AttNode->GetVCurr());

		STLVectorHandler U_N(ac);
		H_NNt.MatVecMul(U_N, VnskMinusVplane);

		for (int ien = 1; ien < ac + 1; ien++) {
			qlcp1(ien, 1) = U_N(ien)*en(1, ien);
		}

		qlcp1.Add(1, 1, U_Nfree);
		qlcp.Put(1, 1, qlcp1);
		qlcp.Sub(ac + 1, 1, U_Tfree);

#if 0
		// SOME OUTPUT TO DEBUG

		std::cout<<" H_NN" << std::endl << H_NN << std::endl;
		std::cout<<" H_TT" << std::endl << H_TT << std::endl;
		std::cout<<" I_P" << std::endl << I_P << std::endl;
		std::cout<<" I_R" << std::endl << I_R << std::endl;
		std::cout<<" I_Rt" << std::endl << I_Rt << std::endl;
		std::cout<<" mu_P" << std::endl << mu_P << std::endl;
		std::cout<<" mu_R" << std::endl << mu_R << std::endl;
		std::cout<<" WNN" << std::endl << WNN << std::endl;
		std::cout<<" WNT" << std::endl << WNT << std::endl;
		std::cout<<" WTN" << std::endl << WTN << std::endl;
		std::cout<<" WTT" << std::endl << WTT << std::endl;
		std::cout<<" IPmT" << std::endl << IPmT << std::endl;
//		std::cout<< "MLCP" << std::endl << MLCP << std::endl;
//		std::cout<< "qLCP" << std::endl << qLCP << std::endl;
#endif

		double Pkp1[ProbDim];
		double Unkp1[ProbDim];

		// NOTE: resetting Pkp1, Unkp1 is mandatory because some times
		// they are not (r)eset in mbdyn_siconos_LCP_call()
		memset(Pkp1, 0, sizeof(Pkp1));
		memset(Unkp1, 0, sizeof(Unkp1));

		if ((NS.niter < NS.niterLCPmax) || (NS.niterLCPmax == 0)) {
			// Solving the LCP with one of Siconos solvers
			mbdyn_siconos_LCP_call(ProbDim, Mlcp.pdGetMat(), qlcp.pdGetMat(), &Pkp1[0], &Unkp1[0], NS.solparam);

			STLVectorHandler PN(ac);
			STLVectorHandler sigmaP(2*ac);
			STLVectorHandler PT(2*ac);
			STLVectorHandler Ptot(3*ac);
			STLVectorHandler ptot(3);

			for (int i = 1; i <= ac; i++) {
				PN(i) = Pkp1[i - 1];
			}
			for (int i = ac + 1; i <= 3*ac; i++) {
				sigmaP(i - ac) = Pkp1[i - 1];
			}
			mu_P.MatVecMul(PT, PN);
			PT -= sigmaP;
			for (int i = 1; i <= 2*ac; i++) {
				Ptot(i) = PT(i);
			}
			for (int i = 2*ac + 1; i <= 3*ac; i++) {
				Ptot(i) = PN(i - (2*ac));
			}
			Htot.MatVecMul(ptot, Ptot);

			p_kp1(1) = ptot(1);
			p_kp1(2) = ptot(2);
			p_kp1(3) = ptot(3);

			NS.Piter = p_kp1;

			NS.niter++;
		}

//		// Solving the LCP with one of Siconos solvers
//		mbdyn_siconos_LCP_call(ProbDim, Mlcp_passed, qlcp_passed,  Pkp1,  Unkp1, NS.solparam);

//		FullMatrixHandler PN(ac, 1);
//		FullMatrixHandler sigmaP(2*ac, 1);
//		FullMatrixHandler PT(2*ac, 1);
//		FullMatrixHandler Ptot(3*ac, 1);
//		FullMatrixHandler ptot(3, 1);
//
//		for (int i = 1; i < ac + 1; i++) {
//			PN(i, 1) = Pkp1[i - 1];
//		}
//		for (int i = ac; i < 3*ac; i++) {
//			sigmaP(i - ac + 1, 1) = Pkp1[i];
//		}
//		PT.MatMul(mu_P, PN);
//		PT.Sub(1, 1, sigmaP);
//		for (int i = 1; i < 2*ac + 1; i++) {
//			Ptot(i, 1) = PT(i, 1);
//		}
//		for (int i = 2*ac + 1; i < 3*ac + 1; i++) {
//			Ptot(i, 1) = PN(i - (2*ac), 1);
//		}
//
//		ptot.MatMul(Htot, Ptot);
//
//		p_kp1(1) = ptot(1, 1);
//		p_kp1(2) = ptot(2, 1);
//		p_kp1(3) = ptot(3, 1);
	}

	// Updating: state actualization with the contribution of the contact impulse
//	NS.Vns_kp1 = vfree_kp1 + invMhat * p_kp1;
	NS.Vns_kp1 = vfree_kp1 + invMhat * NS.Piter;
	NS.Xns_kp1 = NS.Xns_k + (NS.Vns_kp1 * NS.theta_ts + NS.Vns_k* (1. - NS.theta_ts)) * NS.hstep;

	// Reaction force to be exchanged with the rest of the system in MBDYN
	NS.Lambda_kp1 = - (Fg + p_kp1 / NS.hstep - M_ns * (NS.Vns_kp1 - NS.Vns_k) / NS.hstep + NS.Lambda_k*(1. - NS.theta_ts)) / NS.theta_ts;

	// Store Output
	NS.ia = Ia.size();
	NS.pkp1 = p_kp1.Norm();
	NS.Poutput = p_kp1;

	// Indicatore di nuovo step, per evitare troppo output con bVerbose.
	bStepToggle = true;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ModuleNonsmoothNode>;

	if (!SetUDE("nonsmooth" "node", rf)) {
		delete rf;
		silent_cerr("module-nonsmooth-node: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

