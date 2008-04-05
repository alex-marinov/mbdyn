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

/* inertia element */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <float.h>
#include <set>

#include "inertia.h"
#include "dataman.h"
#include "Rot.hh"

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
	S = 0.;
	J = 0.;

	Vec3 B(0.), G(0.);

	for (std::set<const ElemGravityOwner *>::const_iterator i = elements.begin();
		i != elements.end();
		i++)
	{
		dMass += (*i)->dGetM();
		S += (*i)->GetS();
		J += (*i)->GetJ();

		B += (*i)->GetB();
		G += (*i)->GetG();
	}

	J_cm = J;
	if (dMass < DBL_EPSILON) {
		X_cm = 0.;
		V_cm = 0.;
		Omega_cm = 0.;

	} else {
		X_cm = S/dMass;
		V_cm = B/dMass;

		/*
		 * FIXME: should also rotate it in the principal
		 * reference frame, and log the angles
		 */
		J_cm += Mat3x3(S, X_cm);

		Omega_cm = J_cm.Solve(G);
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
		<< " (" << ( GetName() ? GetName() : "unnamed" ) << ")"
		<< std::endl;
	CenterOfMass::Output_int(out)
		<< "    Xcg-X:       " << (X_cm - X0) << std::endl
		<< "    R^T*(Xcg-X): " << R0.MulTV(X_cm - X0) << std::endl
		<< "    Rp:          " << R_princ << std::endl
		<< "    Thetap:      " << RotManip::VecRot(R_princ) << std::endl
		<< "    Jp:          " << J_princ << std::endl;
	return out;
}

/* Costruttore definitivo (da mettere a punto) */
Inertia::Inertia(unsigned int uL, std::set<const ElemGravityOwner *>& elements,
		const Vec3& x0, const Mat3x3& r0, std::ostream& log, flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
CenterOfMass(elements),
flags(0), X0(x0), R0(r0)
{
	if (!X0.IsNull()) {
		flags |= 1;
	}

	if (!(R0 - Eye3).IsNull()) {
		flags |= 2;
	}

	bool done(false);
	if ((fOut & 0x1) & !silent_output) {
		if (!done) {
			Collect_int();
			done = true;
		}
		Output_int(std::cout);
	}

	if (fOut & 0x2) {
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

	if (dMass < DBL_EPSILON) {
		silent_cerr("Inertia(" << GetLabel() << "): "
			"mass is null" << std::endl);

		R_princ = 0.;
		J_princ = 0.;

	} else {
		J_cm.EigSym(J_princ, R_princ);
	}

	if (flags) {
		if (flags & 1) {
			Vec3 DX = X_cm - X0;
			J0 = J + Mat3x3(X0, X0*dMass) + Mat3x3(DX, X0*dMass) + Mat3x3(X0, DX*dMass);
		}

		if (flags & 2) {
			if (flags & 1) {
				J0 = R0*J0*R0.Transpose();
			} else {
				J0 = R0*J*R0.Transpose();
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
	if (fToBeOutput() & (0x1|0x8)) {
		Output_int(std::cout);
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
		break;

	default:
		return 0;
	}

	s++;
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
		return X_cm(which);

	case 1:
		return 0.;

	case 2:
		return V_cm(which);

	case 3:
		return Omega_cm(which);
	}

	throw ErrGeneric();
}
