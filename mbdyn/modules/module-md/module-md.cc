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
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Tingnan Zhang <tingnan1986@gatech.edu>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "dataman.h"
#include "userelem.h"

#include "md_shape.h"
#include "md_sand.h"

// data structure and container declarations
struct NodeData {
	const StructNode *pNode;
	Vec3 F;
	Vec3 M;
};

typedef std::vector<NodeData> NodeContainer;

struct Data {
	NodeContainer m_nodeData;

	// MD-specific data
	MDSand m_sandbox;
	std::vector<MDSphere> m_sphere;
	std::vector<MDCylinder> m_cylinder;
	MDShapeParser m_shape_parser;

	std::string m_sandFileName;
	bool m_bSettle;
	std::string m_fileBallName;
	std::ofstream m_fileBall;

	RealVec m_pos_vec;
	RealVec m_vel_vec;
	RealVec m_ang_vec;
	rmat m_ori_mat;
	RealVec m_foc_vec;
	RealVec m_mom_vec;
	double m_m_ufac;
	double m_l_ufac;
	int m_steps;

	Data(void)
	: m_sandbox(24, 12, 15, 124270), // FIXME: do these data need to be configurable?
	m_bSettle(false),
	m_foc_vec(0, 0, 0), m_mom_vec(0, 0, 0),
	m_m_ufac(1e3), m_l_ufac(1e2), m_steps(0)
	{
		m_ori_mat.resize(3);
		for (int i = 0; i < 3; i++) {
			m_ori_mat[i].resize(3);
		}
	};
};

// examples of MD-specific handlers (they could also be member functions
// of the MBDynMD module class
static int
MD_init(Data& data)
{
	unsigned num_nodes = data.m_nodeData.size();

	data.m_sphere.resize(num_nodes);
	data.m_cylinder.resize(num_nodes);

	// using a structure parser with node number inside it.
	data.m_shape_parser.pshape = new MDShape* [num_nodes];
	data.m_shape_parser.n = num_nodes;
	for (unsigned i = 0; i < num_nodes; i++) {
		// initialize the shape parser
		// shape_parser.pshape[i] = sphere + i;
		// arch[i].ArchInit(1.0, math_pi/2, 0.1, 2.1, 0.1);
		// sphere[i].SphereInit(1.0, 0.01);
		data.m_shape_parser.pshape[i] = &data.m_cylinder[i];
		data.m_cylinder[i].CylinderInit(1.0, 2, 0.01); // FIXME: configurable?
	}

	if (true) {
		// read existing sand bed and began calculation
		// MD_GEN_SAND is used when you initialize some objects
		// inside sandbox. You need to dig those overlapping particles out.
		data.m_sandbox.LoadSand(data.m_sandFileName.c_str());
		silent_cout(std::endl << "MD: loading sand from file \"" << data.m_sandFileName << "\"" << std::endl);
		data.m_sandbox.Init();

	} else {
		data.m_sandbox.GenerateSand();
		data.m_sandbox.Init();
	}
    
	// let the sand bed naturally settle down by runing number of steps
	// not needed if the sand file is already settled.
	if (data.m_bSettle) {
		data.m_sandbox.SwitchOutput(false);
		data.m_sandbox.Run(8e4, data.m_shape_parser, false);
		data.m_sandbox.SwitchOutput(true);
		data.m_sandbox.OutputFile();
		throw NoErr(MBDYN_EXCEPT_ARGS);
	}

	data.m_fileBall.open(data.m_fileBallName.c_str());

	return 0;
}

static int
MD_get_loads(Data& data)
{
	unsigned num_nodes = data.m_nodeData.size();
	for (unsigned n = 0; n < num_nodes; n++) {
		// get position for all nodes
		const Vec3& X(data.m_nodeData[n].pNode->GetXCurr());
		data.m_pos_vec.x = X(1) * data.m_l_ufac;
		data.m_pos_vec.y = X(2) * data.m_l_ufac;
		data.m_pos_vec.z = X(3) * data.m_l_ufac;

		// get velocity for all nodes
		const Vec3& XP(data.m_nodeData[n].pNode->GetVCurr());
		data.m_vel_vec.x = XP(1) * data.m_l_ufac;
		data.m_vel_vec.y = XP(2) * data.m_l_ufac;
		data.m_vel_vec.z = XP(3) * data.m_l_ufac;

		// get orientation for all nodes
		const Mat3x3& R(data.m_nodeData[n].pNode->GetRCurr());
		data.m_ori_mat[0][0] = R(1, 1);
		data.m_ori_mat[0][1] = R(1, 2);
		data.m_ori_mat[0][2] = R(1, 3);
		data.m_ori_mat[1][0] = R(2, 1);
		data.m_ori_mat[1][1] = R(2, 2);
		data.m_ori_mat[1][2] = R(2, 3);
		data.m_ori_mat[2][0] = R(3, 1);
		data.m_ori_mat[2][1] = R(3, 2);
		data.m_ori_mat[2][2] = R(3, 3);

		// get the angular velocity
		const Vec3& Omega(data.m_nodeData[n].pNode->GetWCurr());
		data.m_ang_vec.x = Omega(1);
		data.m_ang_vec.y = Omega(2);
		data.m_ang_vec.z = Omega(3);

		// store those info to shape nodes for interaction with sands
		data.m_shape_parser.pshape[n]->SetNodePosition(data.m_pos_vec);
		data.m_shape_parser.pshape[n]->SetNodeVelocity(data.m_vel_vec);
		data.m_shape_parser.pshape[n]->SetNodeAngularVel(data.m_ang_vec);
		data.m_shape_parser.pshape[n]->SetNodeOrientation(data.m_ori_mat);
	}

	// the sandbox.Run() will also update the current system time
	// as well as kinematics of all particles; We only execute
	// it once per step!
	data.m_sandbox.Run(1, data.m_shape_parser, true);

	for (unsigned n = 0; n < num_nodes; n++) {
		// get calculated force from MD, stored in forceVec
		data.m_shape_parser.pshape[n]->GetNodeForce(data.m_foc_vec);
		data.m_shape_parser.pshape[n]->GetNodeMoment(data.m_mom_vec);

		double tmp_mass;
		data.m_shape_parser.pshape[n]->GetNodeMass(tmp_mass);
		double d = 1./(data.m_l_ufac*data.m_m_ufac);

		if (data.m_steps < 10) {
			d *= double(data.m_steps)/10.0;
		}

		data.m_nodeData[n].F = Vec3(data.m_foc_vec.x*d, data.m_foc_vec.y*d, data.m_foc_vec.z*d);

		d /= data.m_l_ufac;

		data.m_nodeData[n].M = Vec3(data.m_mom_vec.x*d, data.m_mom_vec.y*d, data.m_mom_vec.z*d);

		data.m_fileBall << n
			<< "\t" << data.m_foc_vec.x << "\t" << data.m_foc_vec.y << "\t" << data.m_foc_vec.z
			<< "\t" << data.m_mom_vec.x << "\t" << data.m_mom_vec.y << "\t" << data.m_mom_vec.z
			<< std::endl;
	}

	data.m_steps++;

	silent_cout("MD: end of steps: " << data.m_steps << std::endl);

	return 0;
}

static void
MD_destroy(Data& data)
{
	delete[] data.m_shape_parser.pshape;
}

class MBDynMD
: virtual public Elem, public UserDefinedElem {
private:
	int m_iCoupling;
	int m_iCouplingCounter;

	Data m_data;

public:
	MBDynMD(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~MBDynMD(void);

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
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
	virtual void Output(OutputHandler& OH) const;
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
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
};

MBDynMD::MBDynMD(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_iCoupling(0),
m_iCouplingCounter(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        md\n"
"Author:        Tingnan Zhang <tingnan1986@gatech.edu>\n"
"               Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:	\"Crab Lab\"\n"
"		Georgia Institute of Technology\n"
"		<http://www.physics.gatech.edu/research/goldman/>\n"
"             	Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               <http://www.aero.polimi.it/>\n"
"\n"
"               All rights reserved\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("tight")) {
			m_iCoupling = 1;

		} else if (HP.IsKeyWord("loose")) {
			m_iCoupling = 0;

		} else {
			m_iCoupling = HP.GetInt();
			if (m_iCoupling < 0) {
				silent_cerr("MBDynMD(" << uLabel << "): invalid coupling "
					" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	if (HP.IsKeyWord("sand" "file")) {
		m_data.m_sandFileName = HP.GetFileName();

	} else {
#if MD_GEN_SAND == 1
		m_data.m_sandFileName = "san_low.dat";
#else
		m_data.m_sandFileName = "../saninit.dat";
#endif
		silent_cout("MD: using default sand file \"" << m_data.m_sandFileName << "\"" << std::endl);
	}

	if (HP.IsKeyWord("settle")) {
		m_data.m_bSettle = HP.GetYesNoOrBool();
	}

	if (HP.IsKeyWord("file" "ball")) {
		m_data.m_fileBallName = HP.GetFileName();

	} else {
		m_data.m_fileBallName = "force.txt";
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("MBDynMD(" << uLabel << "): invalid node number " << n
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_data.m_nodeData.resize(n);
	for (unsigned i = 0; i < unsigned(n); i++) {
		m_data.m_nodeData[i].pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		for (unsigned c = 0; c < i; c++) {
			if (m_data.m_nodeData[c].pNode == m_data.m_nodeData[i].pNode) {
				silent_cerr("MBDynMD(" << uLabel << "): "
					"StructNode(" << m_data.m_nodeData[i].pNode->GetLabel() << "), "
					"#" << i << "/" << n << ", "
					"already provided as #" << c << "/" << n
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	if (MD_init(m_data) != 0) {
		silent_cerr("MBDynMD(" << uLabel << "): "
			"MD_init() failed"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

MBDynMD::~MBDynMD(void)
{
	// destroy private data
	MD_destroy(m_data);
}

void
MBDynMD::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (NodeContainer::const_iterator n = m_data.m_nodeData.begin();
			n != m_data.m_nodeData.end(); ++n)
		{
			// format:
			// - for each node
			//   - element label "@" node label
			//   - three components of force
			//   - three components of moment
			out << GetLabel() << "@" << n->pNode->GetLabel()
				<< " " << n->F
				<< " " << n->M
				<< std::endl;
		}
	}
}

void
MBDynMD::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	m_iCouplingCounter = 0;
}

void
MBDynMD::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	// 6 rows for each node; 1 column since only residual is assembled
	*piNumRows = 6*m_data.m_nodeData.size();
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
MBDynMD::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// no contribution to Jacobian matrix
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBDynMD::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(6*m_data.m_nodeData.size());

	// two possible approaches:
	//   1) loose coupling: pass predicted kinematics,
	//	receive forces once, reuse them for all iterations
	//   2) tight coupling: pass actual kinematics,
	//      receive forces at each iteration

	if ((m_iCoupling == 0 && m_iCouplingCounter == 0)
		|| (m_iCoupling > 0 && (m_iCouplingCounter%m_iCoupling) == 0))
	{
		// get loads from MD
		// this helper fills F and M of each node
		if (MD_get_loads(m_data) != 0) {
			silent_cerr("MBDynMD(" << uLabel << "): "
				"MD_get_loads() failed"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	m_iCouplingCounter++;

	integer r = 0;
	for (NodeContainer::const_iterator n = m_data.m_nodeData.begin();
		n != m_data.m_nodeData.end(); ++n)
	{
		integer iFirstIndex = n->pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(r + iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Add(r + 1, n->F);
		WorkVec.Add(r + 4, n->M);

		r += 6;
	}

	return WorkVec;
}

unsigned int
MBDynMD::iGetNumPrivData(void) const
{
	return 0;
}

int
MBDynMD::iGetNumConnectedNodes(void) const
{
	return m_data.m_nodeData.size();
}

void
MBDynMD::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(m_data.m_nodeData.size());

	for (unsigned n = 0; n < m_data.m_nodeData.size(); n++) {
		connectedNodes[n] = m_data.m_nodeData[n].pNode;
	}
}

void
MBDynMD::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	// initialize X and XP according to the initial state of MD as needed
}

std::ostream&
MBDynMD::Restart(std::ostream& out) const
{
	// don't worry about "soft" restart by now
	return out << "# MBDynMD: not implemented" << std::endl;
}

unsigned int
MBDynMD::iGetInitialNumDof(void) const
{
	return 0;
}

void 
MBDynMD::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
MBDynMD::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBDynMD::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<MBDynMD>;

	if (!SetUDE("md", rf)) {
		delete rf;

		silent_cerr("module-md: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

