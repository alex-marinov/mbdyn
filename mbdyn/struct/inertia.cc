/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* inertia element */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <set>
#include <limits>

#include "inertia.h"
#include "dataman.h"
#include "Rot.hh"
#include <sstream> // not sure why this is now needed...

/* CenterOfMass - begin */

std::ostream&
CenterOfMass::Output_int(std::ostream& out) const
{
	return out
		<< "    mass:        " << dMass << std::endl
		<< "    J:           " << J << std::endl
		<< "    Xcg:         " << X_cm << std::endl
		<< "    Jcg:         " << J_cm << std::endl
		<< "    Vcg:         " << V_cm << std::endl
		<< "    Wcg:         " << Omega_cm << std::endl;
}

void
CenterOfMass::Collect_int(void)
{
	dMass = 0.;
	S = Zero3;
	J = Zero3x3;

	Vec3 B(Zero3), G(Zero3);

	for (std::set<const ElemGravityOwner *>::const_iterator i = elements.begin();
		i != elements.end(); ++i)
	{
		dMass += (*i)->dGetM();
		S += (*i)->GetS();
		J += (*i)->GetJ();

		B += (*i)->GetB();
		G += (*i)->GetG();
	}

	J_cm = J;
	if (dMass < std::numeric_limits<doublereal>::epsilon()) {
		X_cm = Zero3;
		V_cm = Zero3;
		Omega_cm = Zero3;
		J_cm = J_cm.Symm();

	} else {
		X_cm = S/dMass;
		V_cm = B/dMass;

		/*
		 * FIXME: should also rotate it in the principal
		 * reference frame, and log the angles
		 */
		J_cm += Mat3x3(MatCrossCross, S, X_cm);
		J_cm = J_cm.Symm();

		ASSERT(J_cm.IsSymmetric()); // NOTE: should be a run time test
		Omega_cm = J_cm.LDLSolve(G - X_cm.Cross(B));
	}
}

/* Costruttore definitivo (da mettere a punto) */
CenterOfMass::CenterOfMass(std::set<const ElemGravityOwner *>& elements) :
elements(elements)
{
	NO_OP;
}

CenterOfMass::~CenterOfMass(void)
{
	NO_OP;
}

/* CenterOfMass - end */

/* Inertia - begin */

	/* momento statico */
Vec3
Inertia::GetS_int(void) const
{
	return S;
}

/* momento d'inerzia */
Mat3x3
Inertia::GetJ_int(void) const
{
	return J;
}

std::ostream&
Inertia::Output_int(std::ostream& out) const
{
	out
		<< "inertia: " << GetLabel()
		<< ( GetName().empty() ? "" : ( std::string(" \"") + GetName() + "\"" ) )
		<< std::endl;
	Vec3 DX(X_cm - X0);
	Mat3x3 JX(J_cm - Mat3x3(MatCrossCross, DX, DX*dMass));
	CenterOfMass::Output_int(out)
		<< "    Xcg-X:       " << DX << std::endl
		<< "    R^T*(Xcg-X): " << R0.MulTV(DX) << std::endl
		<< "    J(X):        " << JX << std::endl
		<< "    R^T*J(X)*R:  " << R0.MulTM(JX)*R0 << std::endl
		<< "    Rp:          " << R_princ << std::endl
		<< "    Thetap:      " << RotManip::VecRot(R_princ) << std::endl
		<< "    Jp:          " << J_princ << std::endl;

	//printf("\n\nmassa %1.16e\n",dMass);
	//Vec3 tmp = R0.MulTV(X_cm-X0);
	//printf("xcg %1.16e %1.16e %1.16e\n", tmp(1), tmp(2), tmp(3) );
	//printf("RR %1.16e %1.16e %1.16e\n", R_princ(1,1), R_princ(1,2), R_princ(1,3) );
	//printf("RR %1.16e %1.16e %1.16e\n", R_princ(2,1), R_princ(2,2), R_princ(2,3) );
	//printf("RR %1.16e %1.16e %1.16e\n", R_princ(3,1), R_princ(3,2), R_princ(3,3) );
	//printf("JJ %1.16e %1.16e %1.16e\n\n\n", J_princ(1), J_princ(2), J_princ(3) );
	return out;
}

/* Costruttore definitivo (da mettere a punto) */
Inertia::Inertia(unsigned int uL, const std::string& sN, std::set<const ElemGravityOwner *>& elements,
		const Vec3& x0, const Mat3x3& r0, std::ostream& log, flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
CenterOfMass(elements),
flags(0), X0(x0), R0(r0)
#ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
,
Var_dMass(0),
Var_X_cm(0),
Var_V_cm(0),
Var_Omega_cm(0),
Var_DX(0),
Var_dx(0),
Var_Jp(0),
Var_Phip(0)
#endif // USE_NETCDFC
{
	this->PutName(sN);

	if (!X0.IsNull()) {
		flags |= 1;
	}

	if (!(R0 - Eye3).IsNull()) {
		flags |= 2;
	}

	bool done(false);
	if ((fOut & OUTPUT_OUT) && (!silent_output)) {
		if (!done) {
			Collect_int();
			done = true;
		}
		Output_int(std::cout);
	}

	if (fOut & OUTPUT_LOG) {
		if (!done) {
			Collect_int();
			done = true;
		}
		Output_int(log);
	}
}

Inertia::~Inertia(void)
{
	NO_OP;
}

/* massa totale */
doublereal
Inertia::dGetM(void) const
{
	return dMass;
}

/* Tipo dell'elemento (usato solo per debug ecc.) */
Elem::Type
Inertia::GetElemType(void) const
{
	return Elem::INERTIA;
}

/* Numero gdl durante l'assemblaggio iniziale */
unsigned int
Inertia::iGetInitialNumDof(void) const
{
	return 0;
}

void
Inertia::Collect_int(void)
{
	CenterOfMass::Collect_int();

	if (dMass < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Inertia(" << GetLabel() << "): "
			"mass is null" << std::endl);

		R_princ = Eye3;
		J_princ = Zero3;

	} else {
		J_cm.EigSym(J_princ, R_princ);
	}

	if (flags) {
		if (flags & 1) {
			Vec3 DX = X_cm - X0;
			J0 = J + Mat3x3(MatCrossCross, X0, X0*dMass)
				+ Mat3x3(MatCrossCross, DX, X0*dMass)
				+ Mat3x3(MatCrossCross, X0, DX*dMass);
		}

		if (flags & 2) {
			if (flags & 1) {
				J0 = R0*J0.MulMT(R0);
			} else {
				J0 = R0*J.MulMT(R0);
			}
		}

	} else {
		J0 = J;
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Inertia::Restart(std::ostream& out) const
{
	return out;
}

void
Inertia::Output(OutputHandler& OH) const
{
	if ((fToBeOutput() & Inertia::OUTPUT_ALWAYS)  && (!silent_output)) {
		if (OH.UseText(OutputHandler::INERTIA)) {
			Output_int(std::cout);
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA)) {

			Vec3 DX(X_cm - X0);
			Vec3 dx = R0.MulTV(DX);
			Vec3 Phip = RotManip::VecRot(R_princ);

			OH.WriteNcVar(Var_dMass, dMass);
			OH.WriteNcVar(Var_X_cm, X_cm);
			OH.WriteNcVar(Var_V_cm, V_cm);
			OH.WriteNcVar(Var_Omega_cm, Omega_cm);

			OH.WriteNcVar(Var_DX, DX);
			OH.WriteNcVar(Var_dx, dx);
			OH.WriteNcVar(Var_Jp, J_princ);
			OH.WriteNcVar(Var_Phip, Phip);
		}
#endif // USE_NETCDF
	}
}

void
Inertia::OutputPrepare_int(OutputHandler &OH, std::string& name)
{
#ifdef USE_NETCDF
	ASSERT(OH.IsOpen(OutputHandler::NETCDF));

	std::ostringstream os;
	os << "elem.inertia." << GetLabel();
	(void)OH.CreateVar(os.str(), "inertia");

	os << ".";
	name = os.str();
#endif // USE_NETCDF
}

void
Inertia::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA))
		{
			std::string name;
			OutputPrepare_int(OH, name);

			Var_dMass = OH.CreateVar<doublereal>(name + "M", "kg",
				"total mass");
			Var_X_cm = OH.CreateVar<Vec3>(name + "X_cm", "m",
				"center of mass position (x, y, z)");
			Var_V_cm = OH.CreateVar<Vec3>(name + "V_cm", "m/s",
				"center of mass velocity (x, y, z)");
			Var_Omega_cm = OH.CreateVar<Vec3>(name + "Omega_cm", "rad/s",
				"center of mass angular velocity (x, y, z)");

			Var_DX = OH.CreateVar<Vec3>(name + "DX", "m",
				"relative center of mass position, global frame (x, y, z)");
			Var_dx = OH.CreateVar<Vec3>(name + "dx", "m",
			 	"relative center of mass position, local frame (x, y, z)");
			Var_Jp = OH.CreateVar<Vec3>(name + "Jp", "kg*m^2",
				"global inertia matrix, w.r.t. principal axes");
			Var_Phip = OH.CreateVar<Vec3>(name + "Phip", "-",
				"orientation vector of principal axes, global frame");
		}
#endif // USE_NETCDF
	}
}

void
Inertia::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Inertia::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
Inertia::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkVec.Resize(0);

	Collect_int();

	return WorkVec;
}

bool 
Inertia::bInverseDynamics(void) const
{
	return true;
}

/* Inverse Dynamics residual assembly */
SubVectorHandler&
Inertia::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  /* XPrimePrimeCurr */,
		InverseDynamics::Order iOrder)
{
	WorkVec.Resize(0);

	Collect_int();

	return WorkVec;
}

/* inverse dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
Inertia::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* Dimensione del workspace durante l'assemblaggio iniziale.
 * Occorre tener conto del numero di dof che l'elemento definisce
 * in questa fase e dei dof dei nodi che vengono utilizzati.
 * Sono considerati dof indipendenti la posizione e la velocita'
 * dei nodi */
void
Inertia::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
Inertia::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
Inertia::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	WorkVec.Resize(0);

	Collect_int();

	return WorkVec;
}

/* Usata per inizializzare la quantita' di moto */
void
Inertia::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

unsigned int
Inertia::iGetNumPrivData(void) const
{
	return
			+3		// X[1-3]
			+3		// Phi[1-3]
			+3		// V[1-3]
			+3		// Omega[1-3]
			+3		// principal inertia moments
			+3*3		// inertia tensor
			+1		// mass
		;
}

unsigned int
Inertia::iGetPrivDataIdx(const char *s) const
{
	unsigned int idx = 0;

	switch (s[0]) {
	case 'X':
		break;

	case 'P':
		if (strncasecmp(s, "Phi", STRLENOF("Phi")) != 0) {
			return 0;
		}
		s += STRLENOF("Phi") - 1;
		idx = 3;
		break;

	case 'V':
		idx = 6;
		break;

	case 'O':
		if (strncasecmp(s, "Omega", STRLENOF("Omega")) != 0) {
			return 0;
		}
		s += STRLENOF("Omega") - 1;
		idx = 9;
		break;

	case 'J':
		if (s[1] == 'P') {
			s++;
			idx = 12;
			break;
		}
		idx = 15;
		break;

	case 'm':
		if (s[1] != '\0') {
			return 0;
		}
		return 25;

	default:
		return 0;
	}

	s++;
	if (s[0] != '[') {
		return 0;
	}

	s++;
	switch (s[0]) {
	case '1':
	case '2':
	case '3':
		idx += s[0] - '0';
		s++;
		if (idx > 15) {
			if (s[0] != ',') {
				return 0;
			}
			s++;
			switch (s[0]) {
			case '1':
			case '2':
			case '3':
				break;

			default:
				return 0;
			}
			idx += 3*(s[0] - '1');
			s++;
		}
		break;

	default:
		return 0;
	}

	//s++;
	if (strcmp(s, "]") != 0) {
		return 0;
	}
	return idx;
}

doublereal
Inertia::dGetPrivData(unsigned int i) const
{

	unsigned int what = (i - 1)/3;
	unsigned int which = (i - 1)%3 + 1;

	switch (what) {
	case 0:
		// center of mass position
		return X_cm(which);

	case 1:
		// principal inertia axes Euler vector components
		return 0.;

	case 2:
		// center of mass velocity
		return V_cm(which);

	case 3:
		// angular velocity
		return Omega_cm(which);

	case 4:
		// principal inertia moments
		return J_princ(which);
	}

	// mass
	if (i == 25) {
		return dMass;
	}

	// inertia tensor
	int ir = (i - 15 - 1)%3 + 1;
	int ic = (i - 15 - 1)/3 + 1;
	
	return J0(ir, ic);

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}
