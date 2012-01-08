/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <cstring>
#include <sstream>

#include "mynewmem.h"
#include "strnode.h"
#include "body.h"
#include "autostr.h"
#include "dataman.h"

#include "matvecexp.h"
#include "Rot.hh"

/* StructNode - begin */

/* Costruttore definitivo */
StructNode::StructNode(unsigned int uL,
	const DofOwner* pDO,
	const Vec3& X0,
	const Mat3x3& R0,
	const Vec3& V0,
	const Vec3& W0,
	const StructNode *pRN,
	const RigidBodyKinematics *pRBK,
	doublereal dPosStiff,
	doublereal dVelStiff,
	bool bOmRot,
	OrientationDescription ood,
	flag fOut)
: Node(uL, pDO, fOut),
RPrev(R0),
RRef(R0),
RCurr(R0),
gRef(Zero3),
gCurr(Zero3),
gPRef(Zero3),
gPCurr(Zero3),
XPrev(X0),
XCurr(X0),
VPrev(V0),
VCurr(V0),
WPrev(W0),
WRef(W0),
WCurr(W0),
XPPCurr(Zero3),
WPCurr(Zero3),
XPPPrev(Zero3),
WPPrev(Zero3),
pRefNode(pRN),
#ifdef USE_NETCDF
Var_X(0),
Var_Phi(0),
Var_XP(0),
Var_Omega(0),
Var_XPP(0),
Var_OmegaP(0),
#endif /* USE_NETCDF */
od(ood),
dPositionStiffness(dPosStiff),
dVelocityStiffness(dVelStiff),
bOmegaRot(bOmRot),
pRefRBK(pRBK),
bOutputAccels(false)
{
	NO_OP;
}

/* Distruttore (per ora e' banale) */
StructNode::~StructNode(void)
{
	NO_OP;
}

/* Tipo di nodo */
Node::Type
StructNode::GetNodeType(void) const
{
	return Node::STRUCTURAL;
}

/* rigid-body kinematics */
const RigidBodyKinematics *
StructNode::pGetRBK(void) const
{
	return pRefRBK;
}

const Vec3&
StructNode::GetX(void) const
{
	return GetXCurr();
}

const Mat3x3&
StructNode::GetR(void) const
{
	return GetRCurr();
}

const Vec3&
StructNode::GetV(void) const
{
	return GetVCurr();
}

const Vec3&
StructNode::GetW(void) const
{
	return GetWCurr();
}

const Vec3&
StructNode::GetXPP(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

const Vec3&
StructNode::GetWP(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

std::ostream&
StructNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position [px,py,pz]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
			"orientation parameters [gx,gy,gz]" << std::endl;

	if (bInitial) {
		iIndex += 6;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"linear velocity [vx,vy,vz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"angular velocity [wx,wy,wz]" << std::endl;
	}

	return out;
}

static const char xyz[] = "xyz";
static const char *dof[] = {
	"position P",
	"incremental rotation parameter g",
	"momentum B",
	"momenta moment G"
};
static const char *eq[] = {
	"momentum definition B",
	"momenta moment definition G",
	"force equilibrium F",
	"moment equilibrium M"
};
static const char *modal_eq[] = {
	"linear velocity definition v",
	"angular velocity definition w",
	"force equilibrium F",
	"moment equilibrium M"
};
static const char *initial_dof[] = {
	"position P",
	"incremental rotation parameter g",
	"velocity v",
	"angular velocity w"
};
static const char *initial_eq[] = {
	"position constraint P",
	"orientation constraint g",
	"position constraint derivative v",
	"orientation constraint derivative w"
};

void
StructNode::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (i == -1) {
		if (bInitial) {
			desc.resize(12);

		} else {
			desc.resize(6);
		}

	} else {
		desc.resize(1);
	}

	std::ostringstream os;
	os << "StructNode(" << GetLabel() << ")";

	// always uses initial_dof[] becuase dof[]
	// and initial_dof[] are the same up to 6
	int iend = bInitial ? 12 : 6;
	if (i == -1) {
		std::string name = os.str();

		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << initial_dof[i/3] << xyz[i%3];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= iend) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		os << ": " << initial_dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

std::ostream&
StructNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	if (bInitial) {
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"position [Px,Py,Pz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"orientation [gx,gy,gz]" << std::endl
			<< prefix << iIndex + 7 << "->" << iIndex + 9 << ": "
				"linear velocity [vx,vy,vz]" << std::endl
			<< prefix << iIndex + 10 << "->" << iIndex + 12 << ": "
				"angular velocity [wx,wy,wz]" << std::endl;
	} else {
		if (dynamic_cast<const DynamicStructNode*>(this) != 0
				|| dynamic_cast<const ModalNode*>(this) != 0)
		{
			iIndex += 6;
		}

		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"force equilibrium [Fx,Fy,Fz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"moment equilibrium [Mx,My,Mz]" << std::endl;
	}

	return out;
}

void
StructNode::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (i == -1) {
		if (bInitial) {
			desc.resize(12);

		} else {
			desc.resize(6);
		}

	} else {
		desc.resize(1);
	}

	std::ostringstream os;
	os << "StructNode(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		if (bInitial) {
			for (i = 0; i < 12; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": " << initial_eq[i/3] << xyz[i%3];
				desc[i] = os.str();
			}

		} else {
			for (i = 0; i < 6; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": " << eq[2 + i/3] << xyz[i%3];
				desc[i] = os.str();
			}
		}

	} else {
		if (bInitial) {
			if (i < 0 || i >= 12) {
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			os << ": " << initial_eq[i/3] << xyz[i%3];

		} else {
			if (i < 0 || i >= 6) {
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			os << ": " << eq[2 + i/3] << xyz[i%3];
		}
		desc[0] = os.str();
	}
}

/* Contributo del nodo strutturale al file di restart */
std::ostream&
StructNode::Restart(std::ostream& out) const
{
	out << "  structural: " << GetLabel() << ", ";
	if (GetStructNodeType() == StructNode::DYNAMIC) {
		out << "dynamic";
	} else if (GetStructNodeType() == StructNode::STATIC) {
		out << "static";
	}
	out << ", reference, global, ";
	XCurr.Write(out, ", ")
		<< ", reference, global, 1, ", (RCurr.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (RCurr.GetVec(2)).Write(out, ", ")
		<< ", reference, global, ",
		VCurr.Write(out, ", ")
		<< ", reference, global, ",
		WCurr.Write(out, ", ") << ", assembly, "
		<< dPositionStiffness << ", "
		<< dVelocityStiffness << ", "
		<< bOmegaRot
		<< ", scale, " << pGetDofOwner()->dGetScale() << ';' << std::endl;

	return out;
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
StructNode::dGetDofValue(int iDof, int iOrder) const
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder == 0 || iOrder == 1);
	if (iDof >= 1 && iDof <= 3) {
		if (iOrder == 0) {
			return XCurr.dGet(iDof);
		} else if (iOrder == 1) {
			return VCurr.dGet(iDof);
		}
	} else if (iDof >= 4 && iDof <= 6) {
		if (iOrder == 1) {
			return WCurr.dGet(iDof - 3);
		} else if (iOrder == 0) {
			silent_cerr("StructNode(" << GetLabel() << "): "
				"unable to return angles" << std::endl);
			throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	} else {
		silent_cerr("StructNode(" << GetLabel() << "): "
			"required dof " << iDof << " (order " << iOrder << ") "
			"is not available." << std::endl);
		throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* dummy return value to workaround compiler complains */
	static doublereal dmy = 0.;
	return dmy;
}

/* Restituisce il valore del dof iDof al passo precedente;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
StructNode::dGetDofValuePrev(int iDof, int iOrder) const
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder == 0 || iOrder == 1);
	if (iDof >= 1 && iDof <= 3) {
		if (iOrder == 0) {
			return XPrev.dGet(iDof);
		} else if (iOrder == 1) {
			return VPrev.dGet(iDof);
		}
	} else if (iDof >= 4 && iDof <= 6) {
		if (iOrder == 1) {
			return WPrev.dGet(iDof - 3);
		} else if (iOrder == 0) {
			silent_cerr("StructNode(" << GetLabel() << "): "
				"unable to return angles" << std::endl);
			throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	} else {
		silent_cerr("StructNode(" << GetLabel() << "): "
			"required dof " << iDof << " (order " << iOrder << ") "
			"is not available." << std::endl);
		throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* dummy return value to workaround compiler complains */
	static doublereal dmy = 0.;
	return dmy;
}

/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void
StructNode::SetDofValue(const doublereal& dValue,
	unsigned int iDof,
	unsigned int iOrder /* = 0 */ )
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder == 0 || iOrder == 1);
	if (iDof >= 1 && iDof <= 3) {
		if (iOrder == 0) {
			XCurr.Put(iDof, dValue);

		} else if (iOrder == 1) {
			VCurr.Put(iDof, dValue);
		}

	} else if (iDof >= 4 && iDof <= 6) {
		if (iOrder == 1) {
			WCurr.Put(iDof - 3, dValue);

		} else if (iOrder == 0) {
			silent_cerr("StructNode(" << GetLabel() << "): "
				"unable to set angles" << std::endl);
			throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		silent_cerr("StructNode(" << GetLabel() << "): "
			"required dof " << iDof << " (order " << iOrder << ") "
			"is not available." << std::endl);
		throw StructNode::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

DofOrder::Order
StructNode::GetDofType(unsigned int i) const
{
	ASSERT(i >= 0 && i < iGetNumDof());
	return DofOrder::DIFFERENTIAL;
}

void
StructNode::OutputPrepare(OutputHandler &OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::STRNODES)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			// node
			const char *type;
			switch (GetStructNodeType()) {
			case STATIC:
				type = "static";
				break;

			case DYNAMIC:
				type = "dynamic";
				break;

			case MODAL:
				type = "modal";
				break;

			case DUMMY:
				type = "dummy";
				break;

			default:
				pedantic_cerr("StructNode::OutputPrepare(" << GetLabel() << "): "
					"warning, unknown node type?" << std::endl);
				type = "unknown";
				break;
			}

			std::ostringstream os;
			os << "node.struct." << GetLabel();
			(void)OH.CreateVar(os.str(), type);

			// node sub-data
			os << '.';
			std::string name = os.str();

			Var_X = OH.CreateVar<Vec3>(name + "X", "m",
				"global position vector (X, Y, Z)");

			Var_Phi = OH.CreateRotationVar(name, "", od, "global");

			Var_XP = OH.CreateVar<Vec3>(name + "XP", "m/s",
				"global velocity vector (v_X, v_Y, v_Z)");

			Var_Omega = OH.CreateVar<Vec3>(name + "Omega", "radian/s",
				"global angular velocity vector (omega_X, omega_Y, omega_Z)");

			// accelerations
			if (bOutputAccels) {
				Var_XPP = OH.CreateVar<Vec3>(name + "XPP", "m/s^2",
					"global acceleration vector (a_X, a_Y, a_Z)");

				Var_OmegaP = OH.CreateVar<Vec3>(name + "OmegaP", "radian/s^2",
					"global angular acceleration vector (omegaP_X, omegaP_Y, omegaP_Z)");
			}
		}
#endif // USE_NETCDF
	}
}

/* Output del nodo strutturale (da mettere a punto) */
void
StructNode::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RCurr)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RCurr)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RCurr)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RCurr);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::STRNODES)) {
			Var_X->put_rec(XCurr.pGetVec(), OH.GetCurrentStep());
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
				break;

			case ORIENTATION_MATRIX:
				Var_Phi->put_rec(RCurr.pGetMat(), OH.GetCurrentStep());
				break;

			default:
				/* impossible */
				break;
			}
			Var_XP->put_rec(VCurr.pGetVec(), OH.GetCurrentStep());
			Var_Omega->put_rec(WCurr.pGetVec(), OH.GetCurrentStep());

			if (bOutputAccels) {
				Var_XPP->put_rec(XPPCurr.pGetVec(), OH.GetCurrentStep());
				Var_OmegaP->put_rec(WPCurr.pGetVec(), OH.GetCurrentStep());
			}
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::STRNODES)) {
			std::ostream& out = OH.StrNodes();
			out
				<< std::setw(8) << GetLabel()
				<< " " << XCurr << " ";
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.StrNodes() << E;
				break;

			case ORIENTATION_MATRIX:
				OH.StrNodes() << RCurr;
				break;

			default:
				/* impossible */
				break;
			}
			OH.StrNodes() << " " << VCurr << " " << WCurr;

			if (bOutputAccels) {
				out
					<< " " << XPPCurr
					<< " " << WPCurr;
			}
			out << std::endl;
		}
	}
}


/* Output della soluzione perturbata (modi ...) */
void
StructNode::Output(OutputHandler& OH,
	const VectorHandler& X,
	const VectorHandler& XP) const
{
	if (fToBeOutput()) {
		integer iFirstIndex = iGetFirstIndex();
		Vec3 DX(X, iFirstIndex + 1);
		Vec3 Dg(X, iFirstIndex + 4);
		Mat3x3 DR(CGR_Rot::MatR, Dg);

		OH.StrNodes() << std::setw(8) << GetLabel()
			<< " " << (XCurr + DX)
			<< " " << MatR2EulerAngles(DR*RCurr)*dRaDegr
			<< " " << "#" << std::endl;
	}
}

/* Aggiorna dati in base alla soluzione */
void
StructNode::Update(const VectorHandler& X, const VectorHandler& XP)
{
	integer iFirstIndex = iGetFirstIndex();

	XCurr = Vec3(X, iFirstIndex + 1);
	VCurr = Vec3(XP, iFirstIndex + 1);

	/* Nota: i g, gP non vengono incrementati */
	gCurr = Vec3(X, iFirstIndex + 4);
	gPCurr = Vec3(XP, iFirstIndex + 4);

#if 0
	// test amplitude of orientation increment
	if (gCurr.Norm() > 1.) {
		silent_cerr("StructNode(" << GetLabel() << "): "
			"incremental rotation too large, YMMV" << std::endl);
	}
#endif

	/* Matrice RDelta, incremento di rotazione da predetto a corrente;
	 * Questo e' piu' efficiente */
	Mat3x3 RDelta(CGR_Rot::MatR, gCurr);

#if 0
	/* Questo e' meno efficiente anche se sembra piu' elegante.
	 * Il problema e' che per scrivere il manipolatore in forma
	 * elegante bisogna aggiungere alla matrice le informazioni
	 * di memorizzazione della funzione di manipolazione.
	 * Oppure occorre un operatore ternario */
	RDelta = MatR << gCurr;
#endif

	/* La matrice di rotazione corrente e' data dalla matrice predetta
	 * (costante) moltiplicata per l'incremento totale occorso;
	 * la velocita' angolare e' data dalla parte incrementale totale
	 * piu' il contributo della velocita' di riferimento (costante) */
	RCurr = RDelta*RRef;
	WCurr = Mat3x3(CGR_Rot::MatG, gCurr)*gPCurr + RDelta*WRef;

#if 0
	/* Nuovo manipolatore (forse e' meno efficiente) */
	WCurr = (CGR_Rot::MatG << gCurr)*gPCurr+RDelta*WRef;
#endif
}


/* Aggiorna dati in base alla soluzione */
void
StructNode::DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP)
{
	integer iFirstIndex = iGetFirstIndex();

	/* Forza configurazione e velocita' al valore iniziale */
	const_cast<VectorHandler &>(X).Put(iFirstIndex + 1, XCurr);
	const_cast<VectorHandler &>(X).Put(iFirstIndex + 4, Zero3);
	const_cast<VectorHandler &>(XP).Put(iFirstIndex + 1, VCurr);
	const_cast<VectorHandler &>(XP).Put(iFirstIndex + 4, Zero3);
}


/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
void
StructNode::InitialUpdate(const VectorHandler& X)
{
	integer iFirstIndex = iGetFirstIndex();

	XCurr = Vec3(X, iFirstIndex + 1);
	VCurr = Vec3(X, iFirstIndex + 7);

	/* Nota: g viene incrementato */
	gCurr = Vec3(X, iFirstIndex + 4);

	Mat3x3 RDelta(CGR_Rot::MatR, gCurr);

	RCurr = RDelta*RRef;
	WCurr = Vec3(X, iFirstIndex + 10);
}

/* Inverse Dynamics: */
void 
StructNode::Update(const VectorHandler& X, InverseDynamics::Order iOrder)
{
	integer iFirstIndex = iGetFirstIndex();
	switch (iOrder)	{
	case InverseDynamics::POSITION: {
		XCurr = Vec3(X, iFirstIndex + 1);
		gCurr = Vec3(X, iFirstIndex + 4);
		Mat3x3 RDelta(CGR_Rot::MatR, gCurr);
		RCurr = RDelta*RRef;
		} break;
		
	case InverseDynamics::VELOCITY: {
		VCurr = Vec3(X, iFirstIndex + 1);
#if 0
		gPCurr = Vec3(X, iFirstIndex + 4);
		Mat3x3 RDelta(CGR_Rot::MatR, gCurr);
		WCurr = Mat3x3(CGR_Rot::MatG, gCurr)*gPCurr + RDelta*WRef;
#endif
		WCurr = Vec3(X, iFirstIndex + 4);
		} break;

	case InverseDynamics::ACCELERATION: {
		XPPCurr = Vec3(X, iFirstIndex + 1);
		WPCurr = Vec3(X, iFirstIndex + 4);
		} break;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
void
StructNode::SetInitialValue(VectorHandler& X)
{
	/* FIXME: why is this called? */
	integer iIndex = iGetFirstIndex();

	X.Put(iIndex + 1, XCurr);
	X.Put(iIndex + 4, Zero3);
	X.Put(iIndex + 7, VCurr);
	X.Put(iIndex + 10, WCurr);
}


void
StructNode::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
#ifdef MBDYN_X_RELATIVE_PREDICTION
	if (pRefNode) {
		Vec3 Xtmp = XPrev - pRefNode->GetXCurr();
		const Mat3x3& R0 = pRefNode->GetRCurr();
		const Vec3& V0 = pRefNode->GetVCurr();
		const Vec3& W0 = pRefNode->GetWCurr();

		XPrev = R0.MulTV(Xtmp);
		RPrev = R0.MulTM(RCurr);
		VPrev = R0.MulTV(VCurr - V0 - W0.Cross(Xtmp));
		WPrev = R0.MulTV(WCurr - W0);

#if 0
		std::cout << "StructNode(" << GetLabel() << "): "
			"SetValue: X=" << XPrev
			<< ", R=" << RPrev
			<< ", V=" << VPrev
			<< ", W=" << WPrev
			<< std::endl;
#endif

	} else
#endif /* MBDYN_X_RELATIVE_PREDICTION */
	{
		/* FIXME: in any case, we start with Crank-Nicolson ... */
		XPrev = XCurr;
		RPrev = RCurr;
		VPrev = VCurr;
		WPrev = WCurr;
	}

	integer iFirstIndex = iGetFirstIndex();
	X.Put(iFirstIndex + 1, XPrev);
	X.Put(iFirstIndex + 4, Zero3);
	gRef = gCurr = gPRef = gPCurr = Zero3;
	XP.Put(iFirstIndex + 1, VPrev);
	XP.Put(iFirstIndex + 4, WPrev);
}


void
StructNode::BeforePredict(VectorHandler& X,
	VectorHandler& XP,
	VectorHandler& XPr,
	VectorHandler& XPPr) const
{
	integer iFirstPos = iGetFirstIndex();

#ifdef MBDYN_X_RELATIVE_PREDICTION
	/* If pRefNode is defined, the prediction is made
	 * on the data in the reference frame it provides */
	if (pRefNode) {

		/*
		   x_r = R_0^T * ( x - x_0 )
		   R_r = R_0^T * R
		   v_r = R_0^T * ( v - v_0 - omega_0 \times ( x - x_0 ) )
		   omega_r = R_0^T * ( omega - omega_0 )
		 */
		Vec3 Xtmp = XCurr - pRefNode->GetXCurr();
		const Mat3x3& R0 = pRefNode->GetRCurr();
		const Vec3& V0 = pRefNode->GetVCurr();
		const Vec3& W0 = pRefNode->GetWCurr();

		XCurr = R0.MulTV(Xtmp);
		RCurr = R0.MulTM(RCurr);
		VCurr = R0.MulTV(VCurr - V0 - W0.Cross(Xtmp));
		WCurr = R0.MulTV(WCurr - W0);

		/* update state vectors with relative position and velocity */
		X.Put(iFirstPos + 1, XCurr);
		XP.Put(iFirstPos + 1, VCurr);
		XPr.Put(iFirstPos + 1, XPrev);
		XPPr.Put(iFirstPos + 1, VPrev);

#if 0
		std::cout << "StructNode(" << GetLabel() << "): "
			"BeforePredict: X=" << XCurr
			<< ", R=" << RCurr
			<< ", V=" << VCurr
			<< ", W=" << WCurr
			<< std::endl;
#endif
	}
#endif /* MBDYN_X_RELATIVE_PREDICTION */

	/* Questa e' la predizione "consistente", ovvero usa come gdl
	 * di rotazione i parametri di rotazione "totali" per predire
	 * la configurazione al nuovo passo, quindi ritorna in forma
	 * incrementale */

	/* Calcolo la matrice RDelta riferita a tutto il passo trascorso
	 * all'indietro */
	Mat3x3 RDelta(RPrev.MulMT(RCurr));

	/* Mi assicuro che g al passo corrente sia nullo */
	X.Put(iFirstPos + 4, Zero3);

	/* Calcolo g al passo precedente attraverso la matrice RDelta riferita
	 * a tutto il passo. Siccome RDelta e' calcolata all'indietro,
	 * i parametri sono gia' con il segno corretto */
	Vec3 gPrev(CGR_Rot::Param, RDelta);
	XPr.Put(iFirstPos + 4, gPrev);

	/* Calcolo gP al passo precedente attraverso la definizione
	 * mediante le Omega. Siccome i parametri sono con il segno meno
	 * e la matrice RDelta e' gia' calcolata all'indietro, l'insieme
	 * e' consistente */
	XPPr.Put(iFirstPos + 4, Mat3x3(CGR_Rot::MatGm1, gPrev)*WPrev);

	/* Metto Omega al passo corrente come gP (perche' G(0) = I) */
	XP.Put(iFirstPos + 4, WCurr);

#if 0
	std::cout
		<< "  " << std::setw(16) << "prev" << std::setw(16) << "curr" << std::setw(16) << GetLabel() << std::endl
		<< "x:" << std::setw(16) << XPrev(1) << std::setw(16) << XCurr(1) << std::endl
		<< "  " << std::setw(16) << XPrev(2) << std::setw(16) << XCurr(2) << std::endl
		<< "  " << std::setw(16) << XPrev(3) << std::setw(16) << XCurr(3) << std::endl
		<< "v:" << std::setw(16) << VPrev(1) << std::setw(16) << VCurr(1) << std::endl
		<< "  " << std::setw(16) << VPrev(2) << std::setw(16) << VCurr(2) << std::endl
		<< "  " << std::setw(16) << VPrev(3) << std::setw(16) << VCurr(3) << std::endl
		<< "g:" << std::setw(16) << gPrev(1) << std::setw(16) << 0 << std::endl
		<< "  " << std::setw(16) << gPrev(2) << std::setw(16) << 0 << std::endl
		<< "  " << std::setw(16) << gPrev(3) << std::setw(16) << 0 << std::endl
		<< "w:" << std::setw(16) << XP(iFirstPos+4) << std::setw(16) << WCurr(1) << std::endl
		<< "  " << std::setw(16) << XP(iFirstPos+5) << std::setw(16) << WCurr(2) << std::endl
		<< "  " << std::setw(16) << XP(iFirstPos+6) << std::setw(16) << WCurr(3) << std::endl;
#endif

	XPrev = XCurr;
	VPrev = VCurr;

	/* Pongo la R al passo precedente uguale a quella corrente
	 * mi servira' se devo ripetere il passo con un diverso Delta t
	 * e per la rettifica dopo la predizione */
	RPrev = RCurr;

	/* Pongo le Omega al passo precedente uguali alle Omega al passo corrente
	 * mi servira' per la correzione dopo la predizione */
	WPrev = WCurr;
}

void
StructNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	integer iFirstIndex = iGetFirstIndex();

	/* Spostamento e velocita' aggiornati */
	XCurr = Vec3(X, iFirstIndex + 1);
	VCurr = Vec3(XP, iFirstIndex + 1);

	/* Ottengo il g predetto */
	gRef = Vec3(X, iFirstIndex + 4);

	/* Calcolo la matrice RDelta derivante dalla predizione */
	Mat3x3 RDelta(CGR_Rot::MatR, gRef);

	/* Calcolo la R corrente in base alla predizione */
	RCurr = RDelta*RPrev;

	/* Calcolo la Omega corrente in base alla predizione (gP "totale") */
	gPRef = Vec3(XP, iFirstIndex + 4);

	/* Calcolo il nuovo Omega */
	WCurr = Mat3x3(CGR_Rot::MatG, gRef)*gPRef;

	/* Resetto i parametri di rotazione e le derivate, g e gP */
	X.Put(iFirstIndex + 4, Zero3);
	XP.Put(iFirstIndex + 4, Zero3);

	gCurr = gPCurr = Zero3;

#ifdef MBDYN_X_RELATIVE_PREDICTION
	if (pRefNode) {

		/*
		   x = x_0 + R_0 * x_r
		   R = R_0 * R_r
		   v = v_0 + omega_0 \times ( R_0 * x_r ) + R_0 * v_r
		   omega = omega_0 + R_0 * omega_r
		 */
		Vec3 X0 = pRefNode->GetXCurr();
		Mat3x3 R0 = pRefNode->GetRCurr();
		Vec3 V0 = pRefNode->GetVCurr();
		Vec3 W0 = pRefNode->GetWCurr();

		XCurr = R0*XCurr;	/* temporary */
		RCurr = R0*RCurr;
		VCurr = V0 + W0.Cross(XCurr) + R0*VCurr;
		WCurr = W0 + R0*WCurr;
		XCurr += X0;		/* plus reference */

		/* alcuni usano anche le predizioni dei parametri
		 * di rotazione e delle loro derivate come riferimento
		 * (approccio updated-updated); quindi calcolo
		 * i parametri di riferimento come i parametri
		 * che danno una predizione pari alla variazione
		 * di R0 piu' l'incremento relativo, e le derivate
		 * dei parametri corrispondenti */
		gRef = Vec3(CGR_Rot::Param, R0*RDelta.MulMT(pRefNode->GetRPrev()));
		gPRef = Mat3x3(CGR_Rot::MatGm1, gRef)*WCurr;

		/* to be safe, the correct values are put back
		 * in the state vectors */
		X.Put(iFirstIndex + 1, XCurr);
		XP.Put(iFirstIndex + 1, VCurr);

#if 0
		std::cout << "StructNode(" << GetLabel() << "): "
			"AfterPredict: X=" << XCurr
			<< ", R=" << RCurr
			<< ", V=" << VCurr
			<< ", W=" << WCurr
			<< std::endl;
#endif
	}
#endif /* MBDYN_X_RELATIVE_PREDICTION */

	RRef = RCurr;
	WRef = WCurr;

#if 0
	/* Ortho check */
	Mat3x3 RRT = RCurr.MulTM(RCurr);
	RRT(1, 1) -= 1.;
	RRT(2, 2) -= 1.;
	RRT(3, 3) -= 1.;
	doublereal dmax = 0.;
	for (int r = 1; r <= 3; r++) {
		for (int c = 1; c <= 3; c++) {
			dmax = std::max(dmax, fabs(RRT(r, c)));
		}
	}
	silent_cout("### StructNode(" << GetLabel() << ") " << dmax << std::endl);
#endif
}

/* Inverse Dynamics: */
void
StructNode::AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP, 
			const VectorHandler& XPP)
{

/* Right now, AfterConvergence is performed only on position 
 * to reset orientation parameters. XPrime and XPrimePrime are 
 * left for compatibility with the virtual method in 
 * class SimulationEntity */

	integer iFirstIndex = iGetFirstIndex();
	
	
	/* Orientation Parameters:
	 * Get g and impose it as gRef: successive iterations 
	 * use gRef as reference and the solution is a perturbation
	 * from it */
	gRef = Vec3(X, iFirstIndex + 4);
	gCurr = Zero3;
	RRef = RCurr;

	XPrev = XCurr;
	RPrev = RCurr;
	VPrev = VCurr;
	WPrev = WCurr;
	XPPPrev = XPPCurr;
	WPPrev = WPCurr;
}

bool
StructNode::ComputeAccelerations(bool b)
{
	return false;
}

/*
 * Metodi per l'estrazione di dati "privati".
 * Si suppone che l'estrattore li sappia interpretare.
 * Come default non ci sono dati privati estraibili
 */
unsigned int
StructNode::iGetNumPrivData(void) const
{
	unsigned i =
		3	// X
		+ 3	// x (R^T * X)
		+ 3	// Phi
		+ 3	// XP
		+ 3	// xP (R^T * XP)
		+ 3	// Omega
		+ 3	// omega (R^T * Omega)
		+ 3	// Euler angles (123)
		+ 3	// Euler angles (313)
		+ 3	// Euler angles (321)
		+ 4;	// Euler parameters

	if (bComputeAccelerations()) {
		i +=
			3	// XPP
			+ 3	// xPP (R^T * XPP)
			+ 3	// OmegaP
			+ 3;	// omegaP (R^T * OmegaP)
	}

	return i;
}

/*
 * Maps a string (possibly with substrings) to a private data;
 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
 * in case of unrecognized data; error must be handled by caller
 */
unsigned int
StructNode::iGetPrivDataIdx(const char *s) const
{
	long	idx;
	char	*next;
	std::string sDataName(s);

	const char	*brk = std::strchr(s, '[' /*]*/ );
	if (brk == 0) {
		return 0;
	}

	size_t	len = brk - s;;
	brk++;

	errno = 0;
	idx = strtol(brk, &next, 10);
	int save_errno = errno;
	if (next == brk || strcmp(next, /*[*/ "]") != 0) {
		return 0;
	}

	if (save_errno == ERANGE) {
		silent_cerr("StructNode(" << GetLabel() << "): "
			"warning, private data index "
			<< std::string(brk, next - brk)
			<< " overflows" << std::endl);
		return 0;
	}

	/*
		X		 0 + idx	idx = {1,3}
		x		 3 + idx	idx = {1,3}
		Phi		 6 + idx	idx = {1,3}
		XP		 9 + idx	idx = {1,3}
		x		12 + idx	idx = {1,3}
		Omega		15 + idx	idx = {1,3}
		omega		18 + idx	idx = {1,3}
		E | E123	21 + idx	idx = {1,3}
		E313		24 + idx	idx = {1,3}
		E321		27 + idx	idx = {1,3}
		PE		31 + idx	idx = {0,3}
		-------------------------------------------
		XPP		34 + idx	idx = {1,3}
		xPP		37 + idx	idx = {1,3}
		OmegaP		40 + idx	idx = {1,3}
		omegaP		43 + idx	idx = {1,3}
	 */

	if (strncmp(s, "PE", len) == 0) {
		if (idx < 0 || idx > 3) {
			return 0;
		}

		return 31 + idx;
	}

	if (idx < 1 || idx > 3) {
		return 0;
	}

	if (strncmp(s, "X", len) == 0) {
		return 0 + idx;
	}

	if (strncmp(s, "x", len) == 0) {
		return 3 + idx;
	}

	if (strncmp(s, "Phi", len) == 0) {
		return 6 + idx;
	}

	if (strncmp(s, "XP", len) == 0) {
		return 9 + idx;
	}

	if (strncmp(s, "xP", len) == 0) {
		return 12 + idx;
	}

	if (strncmp(s, "Omega", len) == 0) {
		return 15 + idx;
	}

	if (strncmp(s, "omega", len) == 0) {
		return 18 + idx;
	}

	if (strncmp(s, "E", len) == 0
		|| strncmp(s, "E123", len) == 0)
	{
		return 21 + idx;
	}

	if (strncmp(s, "E313", len) == 0) {
		return 24 + idx;
	}

	if (strncmp(s, "E321", len) == 0) {
		return 27 + idx;
	}

	bool bca = false;
	unsigned i;
	if (strncmp(s, "XPP", len) == 0) {
		bca = true;
		i = 34 + idx;

	} else if (strncmp(s, "xPP", len) == 0) {
		bca = true;
		i = 37 + idx;

	} else if (strncmp(s, "OmegaP", len) == 0) {
		bca = true;
		i = 40 + idx;

	} else if (strncmp(s, "omegaP", len) == 0) {
		bca = true;
		i = 43 + idx;

	} else {
		// error
		return 0;
	}

	// NOTE: bComputeAccels is set only if iGetPrivDataIdx() is called
	// first; it is not when the (deprecated) idx is directly used.
	if (bca) {
		if (!const_cast<StructNode *>(this)->ComputeAccelerations(true)) {
			silent_cerr("StructNode(" << GetLabel() << "): "
				"request to compute accelerations failed, requested by private data \"" << sDataName << "\""
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	return i;
}

/*
 * Returns the current value of a private data
 * with 0 < i <= iGetNumPrivData()
 */
doublereal
StructNode::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3:
		return XCurr(i);

	case 4:
	case 5:
	case 6:
		return RCurr.GetVec(i - 3)*XCurr;

	case 7:
	case 8:
	case 9: {
		/* TODO */
		Vec3 Phi(RotManip::VecRot(RCurr));
		return Phi(i - 6);
	}

	case 10:
	case 11:
	case 12:
		return VCurr(i - 9);

	case 13:
	case 14:
	case 15:
		return RCurr.GetVec(i - 12)*VCurr;

	case 16:
	case 17:
	case 18:
		return WCurr(i - 15);

	case 19:
	case 20:
	case 21:
		return RCurr.GetVec(i - 18)*WCurr;

	case 22:
	case 23:
	case 24: {
		Vec3 Phi(MatR2EulerAngles123(RCurr));
		return Phi(i - 21);
	}

	case 25:
	case 26:
	case 27: {
		Vec3 Phi(MatR2EulerAngles313(RCurr));
		return Phi(i - 24);
	}

	case 28:
	case 29:
	case 30: {
		Vec3 Phi(MatR2EulerAngles321(RCurr));
		return Phi(i - 27);
	}

	case 31:
	case 32:
	case 33:
	case 34: {
		/* TODO */
		Vec3 e;
		doublereal e0;
		MatR2EulerParams(RCurr, e0, e);
		if (i == 31) {
			return e0;
		}
		return e(i - 31);
	}

	case 35:
	case 36:
	case 37:
		ASSERT(bComputeAccelerations() == true);
		return XPPCurr(i - 34);

	case 38:
	case 39:
	case 40:
		ASSERT(bComputeAccelerations() == true);
		return RCurr.GetVec(i - 37)*XPPCurr;

	case 41:
	case 42:
	case 43:
		ASSERT(bComputeAccelerations() == true);
		return WPCurr(i - 40);

	case 44:
	case 45:
	case 46:
		ASSERT(bComputeAccelerations() == true);
		return RCurr.GetVec(i - 43)*WPCurr;
	}

	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* StructNode - end */


/* DynamicStructNode - begin */

DynamicStructNode::DynamicStructNode(unsigned int uL,
	const DofOwner* pDO,
	const Vec3& X0,
	const Mat3x3& R0,
	const Vec3& V0,
	const Vec3& W0,
	const StructNode *pRN,
	const RigidBodyKinematics *pRBK,
	doublereal dPosStiff,
	doublereal dVelStiff,
	bool bOmRot,
	OrientationDescription ood,
	flag fOut)
: StructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot,
	ood, fOut),
bComputeAccels((fOut & OUTPUT_ACCELERATIONS) == OUTPUT_ACCELERATIONS),
pAutoStr(0)
{
	bOutputAccels = bComputeAccels;
}


/* Distruttore (per ora e' banale) */
DynamicStructNode::~DynamicStructNode(void)
{
	NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type
DynamicStructNode::GetStructNodeType(void) const
{
	return StructNode::DYNAMIC;
}

/* rigid-body kinematics */
const Vec3&
DynamicStructNode::GetXPP(void) const
{
	return GetXPPCurr();
}

const Vec3&
DynamicStructNode::GetWP(void) const
{
	return GetWPCurr();
}

std::ostream&
DynamicStructNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	StructNode::DescribeDof(out, prefix, bInitial);

	if (bInitial == false) {
		out
			<< prefix << iIndex + 7 << "->" << iIndex + 9 << ": "
				"momentum [Bx,By,Bz]" << std::endl
			<< prefix << iIndex + 10 << "->" << iIndex + 12 << ": "
				"momenta moment [Gx,Gy,Gz]" << std::endl;
	}

	return out;
}

void
DynamicStructNode::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (bInitial || i == -1 || (i >= 0 && i < 6)) {
		StructNode::DescribeDof(desc, bInitial, i);

		if (bInitial || (i >= 0 && i < 6)) {
			return;
		}
	}

	if (i == -1) {
		desc.resize(12);

	} else {
		desc.resize(1);
	}
	
	std::ostringstream os;
	os << "StructNode(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();

		for (i = 6; i < 12; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i/3] << xyz[i%3];
			desc[i] = os.str();
		}

	} else {
		if (i < 6 || i >= 12) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		os << ": " << dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

std::ostream&
DynamicStructNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	if (bInitial == false) {
		integer iIndex = iGetFirstIndex();

		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"momentum definition [Bx,By,Bz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"momenta moment definition [Gx,Gy,Gz]" << std::endl;
	}

	StructNode::DescribeEq(out, prefix, bInitial);

	return out;
}

void
DynamicStructNode::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (bInitial || i == -1 || (i >= 6 && i < 12)) {
		int new_i = i;
		if (!bInitial && i != -1) {
			new_i = i - 6;
		}
		StructNode::DescribeEq(desc, bInitial, new_i);

		if (bInitial || (i >= 6 && i < 12)) {
			return;
		}
	}

	if (i == -1) {
		desc.resize(12);
		for (int j = 0; j < 6; j++) {
			desc[6 + j] = desc[j];
		}

	} else {
		desc.resize(1);
	}
	
	std::ostringstream os;
	os << "StructNode(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (i = 0; i < 6; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i/3] << xyz[i%3];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= 6) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		os << ": " << eq[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
 * al posto giusto */
integer
DynamicStructNode::iGetFirstRowIndex(void) const
{
	return iGetFirstMomentumIndex();
}

/* delegate to autostr node */
void
DynamicStructNode::AddInertia(const doublereal& dm, const Vec3& dS,
	const Mat3x3& dJ) const
{
	/* FIXME: do it only if to be output... */
	if (bComputeAccelerations()) {
		pAutoStr->AddInertia(dm, dS, dJ);
	}
}

/* Accesso ai suoi dati */
const Vec3&
DynamicStructNode::GetBCurr(void) const
{
	return pAutoStr->GetBCurr();
}

const Vec3&
DynamicStructNode::GetGCurr(void) const
{
	return pAutoStr->GetGCurr();
}

const Vec3&
DynamicStructNode::GetBPCurr(void) const
{
	return pAutoStr->GetBPCurr();
}

const Vec3&
DynamicStructNode::GetGPCurr(void) const
{
	return pAutoStr->GetGPCurr();
}

bool
DynamicStructNode::ComputeAccelerations(bool b)
{
	bComputeAccels = b;
	return true;
}

void
DynamicStructNode::SetOutputFlag(flag f)
{
	if (f & StructNode::OUTPUT_ACCELERATIONS) {
		// ignore result
		ComputeAccelerations(true);
	}
	ToBeOutput::SetOutputFlag(f);
}

void
DynamicStructNode::Update(const VectorHandler& X, const VectorHandler& XP)
{
	StructNode::Update(X, XP);
	if (bComputeAccelerations()) {
		/* FIXME: pAutoStr is 0 in ModalNode */
		ASSERT(pAutoStr != 0);

		// FIXME: based on values set during previous
		// of AutomaticStructural::AssRes()
		pAutoStr->ComputeAccelerations(XPPCurr, WPCurr);
	}
}

void
DynamicStructNode::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	if (bComputeAccelerations()) {
		/* FIXME: pAutoStr is 0 in ModalNode */
		ASSERT(pAutoStr != 0);

		// FIXME: based on values set during previous
		// of AutomaticStructural::AssRes()
		pAutoStr->ComputeAccelerations(XPPCurr, WPCurr);
	}
}

void
DynamicStructNode::BeforePredict(VectorHandler& X,
	VectorHandler& XP,
	VectorHandler& XPr,
	VectorHandler& XPPr) const
{
	if (bComputeAccelerations()) {
		XPPPrev = XPPCurr;
		WPPrev = WPCurr;
	}

	StructNode::BeforePredict(X, XP, XPr, XPPr);
}

/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
DynamicStructNode::dGetDofValue(int iDof, int iOrder) const
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder >= 0 && iOrder <= 2);

	if (iOrder == 2) {
		/* FIXME: should not happen */
		ASSERT(bComputeAccelerations());
		if (!bComputeAccelerations()) {
			silent_cerr("DynamicStructNode::dGetDofValue("
				<< iDof << "," << iOrder << "): "
				"accelerations are not computed while they should"
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 1
		/* FIXME: might need to compute them in order to be
		 * as up to date as possible; however, elements that contribute
		 * to inertia should assemble first...
		 */
		pAutoStr->ComputeAccelerations(XPPCurr, WPCurr);
#endif

		if (iDof >= 1 && iDof <= 3) {
			return XPPCurr.dGet(iDof);
		} else {
			return WPCurr.dGet(iDof - 3);
		}

	} else {
		return StructNode::dGetDofValue(iDof, iOrder);
	}
}

/* Restituisce il valore del dof iDof al passo precedente;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
DynamicStructNode::dGetDofValuePrev(int iDof, int iOrder) const
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder == 0 || iOrder == 1);

	if (iOrder == 2) {
		/* FIXME: should not happen */
		ASSERT(bComputeAccelerations());
		if (!bComputeAccelerations()) {
			silent_cerr("DynamicStructNode::dGetDofValuePrev("
				<< iDof << "," << iOrder << "): "
				"accelerations are not computed while they should"
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (iDof >= 1 && iDof <= 3) {
			return XPPPrev.dGet(iDof);
		} else {
			return WPPrev.dGet(iDof - 3);
		}
	} else {
		return StructNode::dGetDofValuePrev(iDof, iOrder);
	}
}

/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void
DynamicStructNode::SetDofValue(const doublereal& dValue,
	unsigned int iDof,
	unsigned int iOrder /* = 0 */ )
{
	ASSERT(iDof >= 1 && iDof <= 6);
	ASSERT(iOrder == 0 || iOrder == 1);

	if (iOrder == 2) {
		/* FIXME: should not happen */
		ASSERT(bComputeAccelerations());
		if (!bComputeAccelerations()) {
			silent_cerr("DynamicStructNode::SetDofValue("
				<< dValue << "," << iDof << "," << iOrder << "): "
				"accelerations are not computed while they should"
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (iDof >= 1 && iDof <= 3) {
			XPPCurr.Put(iDof, dValue);

		} else {
			WPCurr.Put(iDof - 3, dValue);
		}

	} else {
		StructNode::SetDofValue(iDof, iOrder);
	}
}

/* DynamicStructNode - end */


/* StaticStructNode - begin */

/* Costruttore definitivo */
StaticStructNode::StaticStructNode(unsigned int uL,
	const DofOwner* pDO,
	const Vec3& X0,
	const Mat3x3& R0,
	const Vec3& V0,
	const Vec3& W0,
	const StructNode *pRN,
	const RigidBodyKinematics *pRBK,
	doublereal dPosStiff,
	doublereal dVelStiff,
	bool bOmRot,
	OrientationDescription ood,
	flag fOut)
: StructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot,
	ood, fOut)
{
	NO_OP;
}


/* Distruttore (per ora e' banale) */
StaticStructNode::~StaticStructNode(void)
{
	NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type
StaticStructNode::GetStructNodeType(void) const
{
	return StructNode::STATIC;
}

/* StaticStructNode - end */


/* ModalNode - begin */

ModalNode::ModalNode(unsigned int uL,
	const DofOwner* pDO,
	const Vec3& X0,
	const Mat3x3& R0,
	const Vec3& V0,
	const Vec3& W0,
	const RigidBodyKinematics *pRBK,
	doublereal dPosStiff,
	doublereal dVelStiff,
	bool bOmRot,
	OrientationDescription ood,
	flag fOut)
: DynamicStructNode(uL, pDO, X0, R0, V0, W0, 0, pRBK,
	dPosStiff, dVelStiff, bOmRot, ood, fOut)
{
	/* XPP and WP are not known in ModalNode */
	ComputeAccelerations(false);
}


/* Distruttore (per ora e' banale) */
ModalNode::~ModalNode(void)
{
	NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type
ModalNode::GetStructNodeType(void) const
{
	return StructNode::MODAL;
}


/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
 * al posto giusto */
integer
ModalNode::iGetFirstRowIndex(void) const
{
	return iGetFirstMomentumIndex();
}

std::ostream&
ModalNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	StructNode::DescribeDof(out, prefix, bInitial);

	if (bInitial == false) {
		integer iIndex = iGetFirstIndex();

		out
			<< prefix << iIndex + 7 << "->" << iIndex + 9 << ": "
				"velocity [vx,vy,vz]" << std::endl
			<< prefix << iIndex + 10 << "->" << iIndex + 12 << ": "
				"angular velocity [wx,wy,wz]" << std::endl;
	}

	return out;
}

void
ModalNode::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (bInitial || i == -1 || (i >= 0 && i < 6)) {
		StructNode::DescribeDof(desc, bInitial, i);

		if (bInitial || (i >= 0 && i < 6)) {
			for (size_t ii = 0; ii < desc.size(); ii++) {
				desc[ii] = "Modal" + desc[ii];
			}
			return;
		}
	}

	if (i == -1) {
		desc.resize(12);

	} else {
		desc.resize(1);
	}
	
	std::ostringstream os;
	os << "ModalStructNode(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();

		for (i = 0; i < 6; i++) {
			desc[i] = "Modal" + desc[i];
		}

		for (i = 6; i < 12; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << initial_dof[i/3] << xyz[i%3];
			desc[i] = os.str();
		}

	} else {
		if (i < 6 || i >= 12) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		os << ": " << initial_dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
	
}

std::ostream&
ModalNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	if (bInitial == false) {
		integer iIndex = iGetFirstIndex();

		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"linear velocity definition [vx,vy,vz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"angular velocity definition [wx,wy,wz]" << std::endl;
	}

	StructNode::DescribeEq(out, prefix, bInitial);

	return out;
}

void
ModalNode::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	if (bInitial || i == -1 || (i >= 6 && i < 12)) {
		int new_i = i;
		if (!bInitial && i != -1) {
			new_i = i - 6;
		}
		StructNode::DescribeEq(desc, bInitial, new_i);

		if (bInitial || (i >= 6 && i < 12)) {
			for (size_t ii = 0; ii < desc.size(); ii++) {
				desc[ii] = "Modal" + desc[ii];
			}
			return;
		}
	}

	if (i == -1) {
		desc.resize(12);
		for (int j = 0; j < 6; j++) {
			desc[6 + j] = "Modal" + desc[j];
		}

	} else {
		desc.resize(1);
	}
	
	std::ostringstream os;
	os << "ModalStructNode(" << GetLabel() << ")";

	const char **xeq = modal_eq;
	if (bInitial) {
		xeq = initial_eq;
	}

	if (i == -1) {
		std::string name = os.str();

		for (i = 0; i < 6; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << xeq[i/3] << xyz[i%3];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= 6) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		os << ": " << xeq[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

/* Aggiorna dati in base alla soluzione */
void
ModalNode::Update(const VectorHandler& X, const VectorHandler& XP)
{
	StructNode::Update(X, XP);

	integer iFirstIndex = iGetFirstIndex();

	/* aggiorno XPP e WP (servono solo a modal.cc) */
	XPPCurr = Vec3(XP, iFirstIndex + 7);
	WPCurr = Vec3(XP, iFirstIndex + 10);
}

/* ModalNode - end */


/* DummyStructNode - begin */

/* Costruttore definitivo */
DummyStructNode::DummyStructNode(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN,
	OrientationDescription ood,
	flag fOut)
: StructNode(uL, pDO, Zero3, Zero3x3, Zero3, Zero3, 0, 0, 0., 0., 0, ood, fOut),
pNode(pN)
{
	ASSERT(pNode != NULL);
}


/* Distruttore (per ora e' banale) */
DummyStructNode::~DummyStructNode(void)
{
	NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type
DummyStructNode::GetStructNodeType(void) const
{
	return StructNode::DUMMY;
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
DummyStructNode::dGetDofValue(int iDof, int iOrder) const
{
	silent_cerr("DummyStructNode(" << GetLabel() << ") has no dofs"
		<< std::endl);
	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}


/* Restituisce il valore del dof iDof al passo precedente;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
DummyStructNode::dGetDofValuePrev(int iDof, int iOrder) const
{
	silent_cerr("DummyStructNode(" << GetLabel() << ") has no dofs"
		<< std::endl);
	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void
DummyStructNode::SetDofValue(const doublereal& dValue,
	unsigned int iDof, unsigned int iOrder)
{
	silent_cerr("DummyStructNode(" << GetLabel() << ") has no dofs"
		<< std::endl);
	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}


/* Aggiorna dati durante l'iterazione fittizia iniziale */
void
DummyStructNode::DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP)
{
	/* posso farlo perche' in genere i dummy nodes si limitano
	 * a copiare i valori di altri nodi, quindi non alterano
	 * le variabili cinematiche */
	Update(X, XP);
}


/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
void
DummyStructNode::InitialUpdate(const VectorHandler& /* X */ )
{
	NO_OP;
}


/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
void
DummyStructNode::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}


void
DummyStructNode::SetValue(DataManager *pDM,
	VectorHandler& X,
	VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	Update(X, XP);
}


/* Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator */
void
DummyStructNode::BeforePredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ ,
	VectorHandler& /* XPrev */ ,
	VectorHandler& /* XPPrev */ ) const
{
	NO_OP;
}


void
DummyStructNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Update(X, XP);
}

bool
DummyStructNode::ComputeAccelerations(bool b)
{
	return const_cast<StructNode *>(pNode)->ComputeAccelerations(b);
}

/* DummyStructNode - end */


/* OffsetDummyStructNode - begin */

/* Costruttore definitivo */
OffsetDummyStructNode::OffsetDummyStructNode(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN,
	const Vec3& f,
	const Mat3x3& R,
	OrientationDescription ood,
	flag fOut)
: DummyStructNode(uL, pDO, pN, ood, fOut), f(f), R(R)
{
	if (pNode->bOutputAccelerations()) {
		bOutputAccels = true;

	} else {
		bOutputAccels = false;
	}

	/* forzo la ricostruzione del nodo strutturale sottostante */
	Update_int();
}


/* Distruttore (per ora e' banale) */
OffsetDummyStructNode::~OffsetDummyStructNode(void)
{
	NO_OP;
}


/* update - interno */
void
OffsetDummyStructNode::Update_int(void)
{
	Vec3 fCurr(pNode->GetRCurr()*f);
	XCurr = pNode->GetXCurr() + fCurr;
	RCurr = pNode->GetRCurr()*R;
	WCurr = pNode->GetWCurr();
	VCurr = pNode->GetVCurr() + WCurr.Cross(fCurr);

	if (bComputeAccelerations()) {
		WPCurr = pNode->GetWPCurr();
		XPPCurr = pNode->GetXPPCurr()
			+ WCurr.Cross(WCurr.Cross(fCurr))
			+ WPCurr.Cross(fCurr);
	}
}


/* Tipo di nodo dummy */
DummyStructNode::Type
OffsetDummyStructNode::GetDummyType(void) const
{
	return DummyStructNode::OFFSET;
}


/* Aggiorna dati in base alla soluzione */
void
OffsetDummyStructNode::Update(const VectorHandler& /* X */ ,
			      const VectorHandler& /* XP */ )
{
	Update_int();
}

/* OffsetDummyStructNode - end */


/* RelFrameDummyStructNode - begin */

/* Costruttore definitivo */
RelFrameDummyStructNode::RelFrameDummyStructNode(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN,
	const StructNode* pNR,
	const Vec3& fh,
	const Mat3x3& Rh,
	OrientationDescription ood,
	flag fOut)
: DummyStructNode(uL, pDO, pN, ood, fOut),
pNodeRef(pNR),
RhT(Rh.Transpose()),
fhT(RhT*fh)
{
	ASSERT(pNodeRef != NULL);

	/*
	 * Note: Rh is transposed from the beginning because it is
	 *       never used directly;
	 *       fh is premultiplied by Rh.Transpose() for the same reason
	 *
	 * Formulas:
	 *
	 * R = RhT * RrT * Rn
	 * X = RhT * RrT * (Xn - Xr)
	 * W = RhT * RrT * (Wn - Wr)
	 * V = RhT * RrT * (Vn - Vr - Wr x (Xn - Xr))
	 *
	 * by defining
	 *
	 * Rn = Rr * Rh * R
	 * Xn = Xr + Rr * (fh + Rh * X)
	 *
	 * and differentiating with respect to time
	 */

	if (pNode->bOutputAccelerations()
		&& pNodeRef->bOutputAccelerations())
	{
		bOutputAccels = true;

	} else {
		bOutputAccels = false;
	}

	/* forzo la ricostruzione del nodo strutturale sottostante */
	Update_int();
}


/* Distruttore (per ora e' banale) */
RelFrameDummyStructNode::~RelFrameDummyStructNode(void)
{
	NO_OP;
}


/* update - interno */
void
RelFrameDummyStructNode::Update_int(void)
{
	Mat3x3 RT(RhT.MulMT(pNodeRef->GetRCurr()));
	Vec3 XRel(pNode->GetXCurr() - pNodeRef->GetXCurr());

	RCurr = RT*pNode->GetRCurr();
	XCurr = RT*XRel - fhT;
	WCurr = RT*(pNode->GetWCurr() - pNodeRef->GetWCurr());

	VCurr = RT*(pNode->GetVCurr()
		- pNodeRef->GetVCurr()
		- pNodeRef->GetWCurr().Cross(XRel));

	if (bComputeAccelerations()) {
		WPCurr = RT*(pNode->GetWPCurr() - pNodeRef->GetWPCurr()
			- pNodeRef->GetWCurr().Cross(pNode->GetWCurr()));
		XPPCurr = RT*(pNode->GetXPPCurr() - pNodeRef->GetXPPCurr()
			- pNodeRef->GetWPCurr().Cross(XRel)
			- (pNodeRef->GetWCurr()*2.).Cross(pNode->GetVCurr() - pNodeRef->GetVCurr())
			+ pNodeRef->GetWCurr().Cross(pNodeRef->GetWCurr().Cross(XRel)));
	}
}


/* Tipo di nodo dummy */
DummyStructNode::Type
RelFrameDummyStructNode::GetDummyType(void) const
{
	return DummyStructNode::RELATIVEFRAME;
}


/* Aggiorna dati in base alla soluzione */
void
RelFrameDummyStructNode::Update(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	Update_int();
}

bool
RelFrameDummyStructNode::ComputeAccelerations(bool b)
{
	bool ok = true;

	if (!const_cast<StructNode *>(pNode)->ComputeAccelerations(b)) {
		ok = false;
	}

	if (!const_cast<StructNode *>(pNodeRef)->ComputeAccelerations(b)) {
		ok = false;
	}

	return ok;
}

/* RelFrameDummyStructNode - end */


/* PivotRelFrameDummyStructNode - begin */

/* Costruttore definitivo */
PivotRelFrameDummyStructNode::PivotRelFrameDummyStructNode(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN,
	const StructNode* pNR,
	const Vec3& fh,
	const Mat3x3& Rh,
	const StructNode* pNR2,
	const Vec3& fh2,
	const Mat3x3& Rh2,
	OrientationDescription ood,
	flag fOut)
: RelFrameDummyStructNode(uL, pDO, pN, pNR, fh, Rh, ood, fOut),
pNodeRef2(pNR2), Rh2(Rh2), fh2(fh2)
{
	ASSERT(pNodeRef2 != NULL);

	/*
	 * Note: Rh is transposed from the beginning because it is
	 *       never used directly;
	 *       fh is premultiplied by Rh.Transpose() for the same reason
	 *
	 * Formulas:
	 *
	 * R = RhT * RrT * Rn
	 * X = RhT * RrT * (Xn - Xr)
	 * W = RhT * RrT * (Wn - Wr)
	 * V = RhT * RrT * (Vn - Vr - Wr x (Xn - Xr))
	 *
	 * by defining
	 *
	 * Rn = Rr * Rh * R
	 * Xn = Xr + Rr * (fh + Rh * X)
	 *
	 * and differentiating with respect to time
	 */

	if (pNode->bOutputAccelerations()
		&& pNodeRef->bOutputAccelerations()
		&& pNodeRef2->bOutputAccelerations())
	{
		bOutputAccels = true;

	} else {
		bOutputAccels = false;
	}

	/* forzo la ricostruzione del nodo strutturale sottostante */
	Update_int();
}


/* Distruttore (per ora e' banale) */
PivotRelFrameDummyStructNode::~PivotRelFrameDummyStructNode(void)
{
	NO_OP;
}


/* update - interno */
void
PivotRelFrameDummyStructNode::Update_int(void)
{
	RelFrameDummyStructNode::Update_int();

	Mat3x3 R2(pNodeRef2->GetRCurr()*Rh2);

	XCurr = pNodeRef2->GetRCurr()*(Rh2*XCurr + fh2);

	if (bComputeAccelerations()) {
		WPCurr = pNodeRef2->GetWPCurr()
			+ pNodeRef2->GetWCurr().Cross(R2*WCurr)
			+ R2*WPCurr;
		XPPCurr = pNodeRef2->GetXPPCurr()
			+ pNodeRef2->GetWPCurr().Cross(XCurr)
			+ pNodeRef2->GetWCurr().Cross(pNodeRef2->GetWCurr().Cross(XCurr))
			+ (pNodeRef2->GetWCurr()*2.).Cross(R2*VCurr)
			+ R2*XPPCurr;
	}

	WCurr = pNodeRef2->GetWCurr() + R2*WCurr;
	VCurr = pNodeRef2->GetVCurr()
		+ pNodeRef2->GetWCurr().Cross(XCurr)
		+ R2*VCurr;
	RCurr = R2*RCurr;
	XCurr += pNodeRef2->GetXCurr();
}


/* Tipo di nodo dummy */
DummyStructNode::Type
PivotRelFrameDummyStructNode::GetDummyType(void) const
{
	return DummyStructNode::PIVOTRELATIVEFRAME;
}


/* Aggiorna dati in base alla soluzione */
void
PivotRelFrameDummyStructNode::Update(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	Update_int();
}

bool
PivotRelFrameDummyStructNode::ComputeAccelerations(bool b)
{
	bool ok = true;

	if (!const_cast<StructNode *>(pNode)->ComputeAccelerations(b)) {
		ok = false;
	}

	if (!const_cast<StructNode *>(pNodeRef)->ComputeAccelerations(b)) {
		ok = false;
	}

	if (!const_cast<StructNode *>(pNodeRef2)->ComputeAccelerations(b)) {
		ok = false;
	}

	return ok;
}

/* RelFrameDummyStructNode - end */


/* Legge un nodo strutturale */

Node*
ReadStructNode(DataManager* pDM,
	MBDynParser& HP,
	DofOwner* pDO,
	unsigned int uLabel)
{
	DEBUGCOUT("Entering ReadStructNode(" << uLabel << ")" << std::endl);

	const char* sKeyWords[] = {
		"static",
		"dynamic",
		"modal",
		"dummy",
			"offset",
			"relative" "frame",
		0
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		STATIC = 0,
		DYNAMIC,
		MODAL,
		DUMMY,

		OFFSET,
		RELATIVEFRAME,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura dati specifici */
	KeyWords CurrType((KeyWords)HP.IsKeyWord());

	/*
	 * explicit node type required; default is no longer "DYNAMIC"
	 */
	if (CurrType == UNKNOWN) {
		silent_cerr("StructNode(" << uLabel << "): "
			"missing node type at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef DEBUG
	if (CurrType == STATIC) {
		std::cout << "Static structural node" << std::endl;
	} else if (CurrType == DYNAMIC) {
		std::cout << "Dynamic structural node" << std::endl;
	} else if (CurrType == DUMMY) {
		std::cout << "Dummy structural node" << std::endl;
	} else if (CurrType == MODAL) {
		std::cout << "Modal node" << std::endl;
	} else {
		std::cout << "Unknown structural node" << std::endl;
	}
#endif /* DEBUG */

	StructNode* pNd = NULL;
	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	KeyWords DummyType = UNKNOWN;
	if (CurrType == DUMMY) {
		StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

		DummyType = KeyWords(HP.GetWord());
		switch (DummyType) {
		case OFFSET: {
			ReferenceFrame RF(pNode);
			Vec3 f(HP.GetPosRel(RF));
			Mat3x3 R(HP.GetRotRel(RF));

			od = ReadOptionalOrientationDescription(pDM, HP);

			flag fOut = pDM->fReadOutput(HP, Node::STRUCTURAL);
			SAFENEWWITHCONSTRUCTOR(pNd,
				OffsetDummyStructNode,
				OffsetDummyStructNode(uLabel, pDO, pNode,
					f, R, od, fOut));
		} break;

		case RELATIVEFRAME: {
			StructNode* pNodeRef = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

			ReferenceFrame RF(pNodeRef);

			Vec3 fh(Zero3);
			if (HP.IsKeyWord("position")) {
				fh = HP.GetPosRel(RF);
			}

			Mat3x3 Rh(Eye3);
			if (HP.IsKeyWord("orientation")) {
				Rh = HP.GetRotRel(RF);
			}

			od = ReadOptionalOrientationDescription(pDM, HP);

			StructNode *pNodeRef2 = 0;
			Vec3 fh2(Zero3);
			Mat3x3 Rh2(Eye3);
			if (HP.IsKeyWord("pivot" "node")) {
				pNodeRef2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

				if (HP.IsKeyWord("position")) {
					fh2 = HP.GetPosRel(RF);
				}

				if (HP.IsKeyWord("orientation")) {
					Rh2 = HP.GetRotRel(RF);
				}
			}

			flag fOut = pDM->fReadOutput(HP, Node::STRUCTURAL);

			if (pNodeRef2) {
				SAFENEWWITHCONSTRUCTOR(pNd,
					PivotRelFrameDummyStructNode,
					PivotRelFrameDummyStructNode(uLabel, pDO,
						pNode, pNodeRef, fh, Rh,
						pNodeRef2, fh2, Rh2, od, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pNd,
					RelFrameDummyStructNode,
					RelFrameDummyStructNode(uLabel, pDO,
						pNode, pNodeRef, fh, Rh, od, fOut));
			}
		} break;

		default:
			silent_cerr("StructNode(" << uLabel << "): "
				"unknown dummy node type "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		/* posizione (vettore di 3 elementi) */
		if (!HP.IsKeyWord("position")) {
			pedantic_cerr("StructNode(" << uLabel << "): "
				"missing keyword \"position\" at line "
				<< HP.GetLineData() << std::endl);
		}
		Vec3 X0(HP.GetPosAbs(AbsRefFrame));
		DEBUGCOUT("X0 =" << std::endl << X0 << std::endl);

		/* sistema di riferimento (trucco dei due vettori) */
		if (!HP.IsKeyWord("orientation")) {
			pedantic_cerr("StructNode(" << uLabel << "): "
				"missing keyword \"orientation\" at line "
				<< HP.GetLineData() << std::endl);
		}
		Mat3x3 R0(HP.GetRotAbs(AbsRefFrame));
		DEBUGCOUT("R0 =" << std::endl << R0 << std::endl);

		/* Velocita' iniziali (due vettori di 3 elementi, con la possibilita'
		 * di usare "null" per porli uguali a zero) */
		if (!HP.IsKeyWord("velocity")) {
			pedantic_cerr("StructNode(" << uLabel << "): "
				"missing keyword \"velocity\" at line "
				<< HP.GetLineData() << std::endl);
		}
		Vec3 XPrime0(HP.GetVelAbs(AbsRefFrame, X0));

		if (!HP.IsKeyWord("angular" "velocity")) {
			pedantic_cerr("StructNode(" << uLabel << "): "
				"missing keyword \"angular velocity\" at line "
				<< HP.GetLineData() << std::endl);
		}
		Vec3 Omega0(HP.GetOmeAbs(AbsRefFrame));
		DEBUGCOUT("Xprime0 =" << std::endl << XPrime0 << std::endl
			<< "Omega0 =" << std::endl << Omega0 << std::endl);

		const StructNode *pRefNode = 0;
		if (HP.IsKeyWord("prediction" "node")) {
			switch (CurrType) {
			case STATIC:
			case DYNAMIC:
				break;

			default:
				silent_cerr("StructNode(" << uLabel << "): "
					"prediction node allowed "
					"for static and dynamic nodes only, "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pRefNode = dynamic_cast<const StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

#ifndef MBDYN_X_RELATIVE_PREDICTION
			silent_cerr("warning, relative prediction disabled; "
				"absolute prediction will be used" << std::endl);
#endif /* ! MBDYN_X_RELATIVE_PREDICTION */
		}

		const RigidBodyKinematics *pRBK = pDM->pGetRBK();

		/* Rigidezza in assemblaggio diversa da quella di default
		 * e flag di output */
		doublereal dPosStiff = pDM->dGetInitialPositionStiffness();
		doublereal dVelStiff = pDM->dGetInitialVelocityStiffness();
		bool bOmRot = pDM->bDoesOmegaRotate();

		if (HP.IsArg()) {
			if (HP.IsKeyWord("assembly")) {
				dPosStiff = HP.GetReal(dPosStiff);
				dVelStiff = HP.GetReal(dVelStiff);

				bOmRot = HP.GetYesNoOrBool(bOmRot);

				DEBUGCOUT("Initial position stiffness: " << dPosStiff << std::endl);
				DEBUGCOUT("Initial velocity stiffness: " << dVelStiff << std::endl);
				DEBUGCOUT("Omega rotates? : " << (bOmRot ? "yes" : "no") << std::endl);
			}
		}

		pDO->SetScale(pDM->dReadScale(HP, DofOwner::STRUCTURALNODE));

		od = ReadOptionalOrientationDescription(pDM, HP);

		flag fOut = pDM->fReadOutput(HP, Node::STRUCTURAL);
		if (CurrType == DYNAMIC || CurrType == MODAL) {
			bool bGotAccels(false);
			bool bAccels(pDM->bOutputAccelerations());

			bool bGotInertia(false);
			bool bInertia(false);

			while (HP.IsArg()) {
				if (HP.IsKeyWord("accelerations")) {
					if (bGotAccels) {
						silent_cerr("StructNode(" << uLabel << "): "
							"\"accelerations\" already set, "
							"repeated at line " << HP.GetLineData()
							<< std::endl);
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					bGotAccels = true;

					if (HP.IsArg()) {
						if (HP.GetYesNoOrBool(false)) {
							bAccels = true;

						} else {
							bAccels = false;
						}

					} else {
						// deprecated
						silent_cout("StructNode(" << uLabel << "): "
							"warning, \"accelerations\" needs \"yes\" or \"no\" "
							"at line " << HP.GetLineData() << std::endl);
						bAccels = true;
					}

				} else if (HP.IsKeyWord("output" "inertia")) {
					if (bGotInertia) {
						silent_cerr("StructNode(" << uLabel << "): "
							"\"inertia\" already set, "
							"repeated at line " << HP.GetLineData()
							<< std::endl);
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					bGotInertia = true;

					if (!HP.IsArg()) {
						silent_cerr("StructNode(" << uLabel << "): "
							"\"inertia\" needs \"yes\" or \"no\" "
							"at line " << HP.GetLineData() << std::endl);
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (HP.GetYesNoOrBool(true)) {
						bInertia = true;

					} else {
						bInertia = false;
					}

				} else {
					break;
				}
			}

			if (fOut) {
				// restore legacy behavior
				if (!bGotInertia) {
					bInertia = true;
				}

				if (bInertia) {
					fOut |= StructNode::OUTPUT_INERTIA;
				}

				if (bAccels) {
					fOut |= StructNode::OUTPUT_ACCELERATIONS;
				}
			}
		}

		if (CurrType == DYNAMIC && pDM->bIsStaticModel()) {
			pedantic_cout("DynamicStructNode(" << uLabel << ") turned into static" << std::endl);
			CurrType = STATIC;
		}

		/* Se non c'e' il punto e virgola finale */
		if (HP.IsArg()) {
			silent_cerr("ReadStructNode(" << uLabel << "): semicolon expected "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* costruzione del nodo */
		if (CurrType == STATIC) {
			SAFENEWWITHCONSTRUCTOR(pNd, StaticStructNode,
				StaticStructNode(uLabel, pDO,
					X0, R0,
					XPrime0, Omega0,
					pRefNode,
					pRBK,
					dPosStiff, dVelStiff,
					bOmRot, od, fOut));

		} else if (CurrType == DYNAMIC) {
			SAFENEWWITHCONSTRUCTOR(pNd, DynamicStructNode,
				DynamicStructNode(uLabel, pDO,
					X0, R0,
					XPrime0, Omega0,
					pRefNode,
					pRBK,
					dPosStiff, dVelStiff,
					bOmRot, od, fOut));

			/* Incrementa il numero di elementi automatici dei nodi dinamici */
			pDM->IncElemCount(Elem::AUTOMATICSTRUCTURAL);

		} else if (CurrType == MODAL) {
			SAFENEWWITHCONSTRUCTOR(pNd, ModalNode,
				ModalNode(uLabel, pDO,
					X0, R0,
					XPrime0, Omega0,
					pRBK,
					dPosStiff, dVelStiff,
					bOmRot, od, fOut));
		}
	}

	std::ostream& out = pDM->GetLogFile();

	const char *description = "structural node: ";

	switch (CurrType){
	case DUMMY:
		switch (DummyType){
		case RELATIVEFRAME:
			description = "relative frame structural node: ";
			break;

		default:
			break;
		}
		break;

	default:
		break;
	}

	out << description << uLabel
		<< " ", pNd->GetXCurr().Write(out, " ")
		<< " ";
	switch (od) {
	case EULER_123:
		out << "euler123 ",
			(MatR2EulerAngles123(pNd->GetRCurr())*dRaDegr).Write(out, " ");
		break;

	case EULER_313:
		out << "euler313 ",
			(MatR2EulerAngles313(pNd->GetRCurr())*dRaDegr).Write(out, " ");
		break;

	case EULER_321:
		out << "euler321 ",
			(MatR2EulerAngles321(pNd->GetRCurr())*dRaDegr).Write(out, " ");
		break;

	case ORIENTATION_VECTOR:
		out << "phi ",
			RotManip::VecRot(pNd->GetRCurr()).Write(out, " ");
		break;

	case ORIENTATION_MATRIX:
		out << "mat ",
			pNd->GetRCurr().Write(out, " ");
		break;

	default:
		/* impossible */
		break;
	}

	out << std::endl;

	ASSERT(pNd != NULL);

	return pNd;
} /* End of ReadStructNode() */

