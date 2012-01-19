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

#include <cmath>
#include <cfloat>

#include "autostr.h"

/* AutomaticStructDispElem - begin */

/* Costruttore */
AutomaticStructDispElem::AutomaticStructDispElem(const DynamicStructDispNode* pN)
: Elem(pN->GetLabel(), pN->fToBeOutput()),
pNode(const_cast<DynamicStructDispNode *>(pN)), B(Zero3), BP(Zero3),
m(0.)
{
	pNode->SetAutoStr(this);
}

void
AutomaticStructDispElem::ComputeAccelerations(Vec3& XPP) const
{
	if (m == 0.) {
		XPP = Zero3;
		return;
	}

	XPP = BP/m;
}

void
AutomaticStructDispElem::AddInertia(const doublereal& dm)
{
	m += dm;
}

/* inizializza i dati */
void
AutomaticStructDispElem::Init(const Vec3& b, const Vec3& bp)
{
	B = b;
	BP = bp;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AutomaticStructDispElem::Restart(std::ostream& out) const
{
	out << "automatic structural displacement: " << GetLabel() << ", "
		"reference, global, ", B.Write(out, ", ") << ", "
		"reference, global, ", BP.Write(out, ", ") << ";" << std::endl;

	return out;
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
AutomaticStructDispElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssJac");

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iCoefs = 6;
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		iCoefs += 6;
	}
	WM.ResizeReset(iCoefs, 0);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i
	 */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutItem(iCnt, iFirstPositionIndex + iCnt,
			iFirstMomentumIndex + iCnt, -dCoef);
		WM.PutItem(3 + iCnt, iFirstMomentumIndex + iCnt,
			iFirstMomentumIndex + iCnt, 1.);
	}

	// relative frame dynamics contribution
	// (see tecman, "Dynamics in a Relative Reference Frame")
	if (pRBK) {
		const Vec3& W0 = pRBK->GetW();

		WM.PutCross(7, iFirstMomentumIndex,
			iFirstMomentumIndex, W0*(2.*dCoef));
	}

	return WorkMat;
}


/* assemblaggio autoval */
void
AutomaticStructDispElem::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssMats");

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WMA = WorkMatA.SetSparse();
	SparseSubMatrixHandler& WMB = WorkMatB.SetSparse();

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(3, 0);
	WMB.ResizeReset(3, 0);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i
	 */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutItem(iCnt, iFirstPositionIndex + iCnt,
			iFirstMomentumIndex + iCnt, -1.);
		WMB.PutItem(iCnt, iFirstMomentumIndex + iCnt,
			iFirstMomentumIndex + iCnt, 1.);
	}
}


/* assemblaggio residuo */
SubVectorHandler&
AutomaticStructDispElem::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssRes");

	WorkVec.ResizeReset(6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Collects data */
	B = Vec3(XCurr, iFirstMomentumIndex + 1);
	BP = Vec3(XPrimeCurr, iFirstMomentumIndex + 1);

	/*
	 * Momentum and momenta moment (about node):
	 *
	 * B = m V + W /\ S
	 *
	 * G = S /\ V + J W
	 *
	 * Bp = F
	 *
	 * Gp + V /\ B = M
	 */
	WorkVec.Add(1, B);
	WorkVec.Sub(4, BP);

	// relative frame dynamics contribution
	// (see tecman, "Dynamics in a Relative Reference Frame")
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		const Vec3& W0 = pRBK->GetW();

		WorkVec.Sub(4, W0.Cross(B*2.));
	}

	// reset instantaneous inertia properties
	m = 0.;

	return WorkVec;
}


void
AutomaticStructDispElem::OutputPrepare(OutputHandler &OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			/* get a pointer to binary NetCDF file
			 * -->  pDM->OutHdl.BinFile */
			NcFile *pBinFile = OH.pGetBinFile();
			char buf[BUFSIZ];

			// NOTE: we hijack the structural nodes namespace
			int l = snprintf(buf, sizeof(buf), "node.struct.%lu.",
				(unsigned long)GetLabel());
			// NOTE: "BP" is the longest var name
			if (l < 0 || l >= int(sizeof(buf) - STRLENOF("BP"))) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			/* Add NetCDF (output) variables to the BinFile object
			 * and save the NcVar* pointer returned from add_var
			 * as handle for later write accesses.
			 * Define also variable attributes */

			strcpy(&buf[l], "B");
			Var_B = pBinFile->add_var(buf, ncDouble,
				OH.DimTime(), OH.DimV3());
			if (Var_B == 0) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_B->add_att("units", "kg m/s")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_B->add_att("type", "Vec3")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_B->add_att("description",
				"momentum (X, Y, Z)"))
			{
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			strcpy(&buf[l], "G");
			Var_G = pBinFile->add_var(buf, ncDouble,
				OH.DimTime(), OH.DimV3());
			if (Var_G == 0) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_G->add_att("units", "kg m^2/s")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_G->add_att("type", "Vec3")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_G->add_att("description",
				"momenta moment (X, Y, Z)"))
			{
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			strcpy(&buf[l], "BP");
			Var_BP = pBinFile->add_var(buf, ncDouble,
				OH.DimTime(), OH.DimV3());
			if (Var_BP == 0) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_BP->add_att("units", "kg m/s^2")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_BP->add_att("type", "Vec3")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_BP->add_att("description",
				"momentum derivative (X, Y, Z)"))
			{
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			strcpy(&buf[l], "GP");
			Var_GP = pBinFile->add_var(buf, ncDouble,
				OH.DimTime(), OH.DimV3());
			if (Var_GP == 0) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_GP->add_att("units", "kg m^2/s^2")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_GP->add_att("type", "Vec3")) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!Var_GP->add_att("description",
				"momenta moment derivative (X, Y, Z)"))
			{
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
#endif // USE_NETCDF
	}
}

void
AutomaticStructDispElem::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA)) {
			Var_B->put_rec(B.pGetVec(), OH.GetCurrentStep());
			Var_G->put_rec(::Zero3.pGetVec(), OH.GetCurrentStep());
			Var_BP->put_rec(BP.pGetVec(), OH.GetCurrentStep());
			Var_GP->put_rec(::Zero3.pGetVec(), OH.GetCurrentStep());
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::INERTIA)) {
			OH.Inertia() << std::setw(8) << GetLabel()
				<< " " << B
				<< " " << ::Zero3
				<< " " << BP
				<< " " << ::Zero3
				<< std::endl;
		}
	}
}

/* Setta i valori iniziali delle variabili (e fa altre cose)
 * prima di iniziare l'integrazione */
void
AutomaticStructDispElem::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iIndex = pNode->iGetFirstMomentumIndex();

	XP.Put(iIndex + 1, BP);
}

/* Dati privati */
unsigned int
AutomaticStructDispElem::iGetNumPrivData(void) const
{
	return 12;
}

unsigned int
AutomaticStructDispElem::iGetPrivDataIdx(const char *s) const
{
	/*
	 * beta[1]
	 * beta[2]
	 * beta[3]
	 * gamma[1]
	 * gamma[2]
	 * gamma[3]
	 * betaP[1]
	 * betaP[2]
	 * betaP[3]
	 * gammaP[1]
	 * gammaP[2]
	 * gammaP[3]
	 */
	unsigned int idx = 0;
	if (strncmp(s, "beta", STRLENOF("beta")) == 0) {
		s += STRLENOF("beta");
	} else if (strncmp(s, "gamma", STRLENOF("gamma")) == 0) {
		s += STRLENOF("gamma");
		idx += 3;
	} else {
		return 0;
	}

	if (s[0] == 'P') {
		s++;
		idx += 6;
	}

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
		break;

	default:
		return 0;
	}

	if (s[0] != ']' && s[1] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
AutomaticStructDispElem::dGetPrivData(unsigned int i) const
{
	unsigned int der = (i - 1)/6;
	i -= 6*der;
	unsigned int type = (i - 1)/3;
	i -= 3*type;

	if (der) {
		if (type) {
			return 0.;
		}
		return BP(i);

	} else {
		if (type) {
			return 0.;
		}
		return B(i);
	}
}

/* AutomaticStructDispElem - end */


/* AutomaticStructElem - begin */

/* Costruttore */
AutomaticStructElem::AutomaticStructElem(const DynamicStructNode* pN)
: Elem(pN->GetLabel(), (pN->fToBeOutput() & StructDispNode::OUTPUT_INERTIA) == StructDispNode::OUTPUT_INERTIA),
AutomaticStructDispElem(pN),
G(Zero3), GP(Zero3),
S(Zero3), J(Zero3x3)
{
	ASSERT(dynamic_cast<const DynamicStructNode *>(pNode) != 0);
	pNode->SetAutoStr(this);
}

void
AutomaticStructElem::ComputeAccelerations(Vec3& XPP, Vec3& WP) const
{
	if (m == 0.) {
		XPP = Zero3;
		WP = Zero3;
		return;
	}

	Vec3 Xcg = S/m;
	Mat3x3 Jcg = J + Mat3x3(MatCrossCross, Xcg, S);
	const Vec3& V = pNode->GetVCurr();
	const Vec3& W = dynamic_cast<const DynamicStructNode *>(pNode)->GetWCurr();
	WP = Jcg.LDLSolve(GP - Xcg.Cross(BP) - W.Cross(Jcg*W) + V.Cross(B));
	XPP = (BP - WP.Cross(S) - W.Cross(W.Cross(S)))/m;
}

void
AutomaticStructElem::AddInertia(const doublereal& dm, const Vec3& dS,
		const Mat3x3& dJ)
{
	m += dm;
	S += dS;
	J += dJ;
}

/* inizializza i dati */
void
AutomaticStructElem::Init(const Vec3& b, const Vec3& g,
			  const Vec3& bp, const Vec3& gp)
{
	B = b;
	G = g;
	BP = bp;
	GP = gp;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AutomaticStructElem::Restart(std::ostream& out) const
{
	out << "automatic structural: " << GetLabel() << ", "
		"reference, global, ", B.Write(out, ", ") << ", "
		"reference, global, ", G.Write(out, ", ") << ", "
		"reference, global, ", BP.Write(out, ", ") << ", "
		"reference, global, ", GP.Write(out, ", ") << ";" << std::endl;

	return out;
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
AutomaticStructElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssJac");

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iCoefs = 24;
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		iCoefs += 12;
	}
	WM.ResizeReset(iCoefs, 0);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i
	 */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutItem(iCnt, iFirstPositionIndex + iCnt,
			iFirstMomentumIndex + iCnt, -dCoef);
		WM.PutItem(6+iCnt, iFirstMomentumIndex + iCnt,
			iFirstMomentumIndex + iCnt, 1.);
	}

	WM.PutCross(13, iFirstMomentumIndex + 3,
		iFirstMomentumIndex, pNode->GetVCurr()*dCoef);
	WM.PutCross(19, iFirstMomentumIndex + 3,
		iFirstPositionIndex, -B);

	// relative frame dynamics contribution
	// (see tecman, "Dynamics in a Relative Reference Frame")
	if (pRBK) {
		const Vec3& W0 = pRBK->GetW();

		WM.PutCross(25, iFirstMomentumIndex,
			iFirstMomentumIndex, W0*(2.*dCoef));
		WM.PutCross(31, iFirstMomentumIndex + 3,
			iFirstMomentumIndex + 3, W0*dCoef);
	}

	return WorkMat;
}


/* assemblaggio autoval */
void
AutomaticStructElem::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssMats");

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WMA = WorkMatA.SetSparse();
	SparseSubMatrixHandler& WMB = WorkMatB.SetSparse();

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(12, 0);
	WMB.ResizeReset(12, 0);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i
	 */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutItem(iCnt, iFirstPositionIndex + iCnt,
			iFirstMomentumIndex + iCnt, -1.);
		WMB.PutItem(iCnt, iFirstMomentumIndex + iCnt,
			iFirstMomentumIndex + iCnt, 1.);
	}

	WMA.PutCross(7, iFirstMomentumIndex + 3, iFirstMomentumIndex,
		pNode->GetVCurr());
	WMB.PutCross(7, iFirstMomentumIndex + 3, iFirstPositionIndex, -B);
}


/* assemblaggio residuo */
SubVectorHandler&
AutomaticStructElem::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AutomaticStructElem::AssRes");

	WorkVec.ResizeReset(12);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 12; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	/* Collects data */
	B = Vec3(XCurr, iFirstMomentumIndex + 1);
	G = Vec3(XCurr, iFirstMomentumIndex + 4);
	BP = Vec3(XPrimeCurr, iFirstMomentumIndex + 1);
	GP = Vec3(XPrimeCurr, iFirstMomentumIndex + 4);

	/*
	 * Momentum and momenta moment (about node):
	 *
	 * B = m V + W /\ S
	 *
	 * G = S /\ V + J W
	 *
	 * Bp = F
	 *
	 * Gp + V /\ B = M
	 */
	WorkVec.Add(1, B);
	WorkVec.Add(4, G);
	WorkVec.Sub(7, BP);
	WorkVec.Sub(10, GP + pNode->GetVCurr().Cross(B));

	// relative frame dynamics contribution
	// (see tecman, "Dynamics in a Relative Reference Frame")
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		const Vec3& W0 = pRBK->GetW();

		WorkVec.Sub(7, W0.Cross(B*2.));
		WorkVec.Sub(10, W0.Cross(G));
	}

	// reset instantaneous inertia properties
	m = 0.;
	S = Zero3;
	J = Zero3x3;

	return WorkVec;
}


void
AutomaticStructElem::OutputPrepare(OutputHandler &OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			std::ostringstream os;
			os << "node.struct." << GetLabel() << ".";

			std::string name(os.str());

			Var_B = OH.CreateVar<Vec3>(name + "B", "kg m/s", "momentum (X, Y, Z)");
			Var_G = OH.CreateVar<Vec3>(name + "G", "kg m^2/s", "momenta moment (X, Y, Z)");
			Var_BP = OH.CreateVar<Vec3>(name + "BP", "kg m/s^2", "momentum derivative (X, Y, Z)");
			Var_GP = OH.CreateVar<Vec3>(name + "GP", "kg m^2/s^2", "momenta moment derivative (X, Y, Z)");
		}
#endif // USE_NETCDF
	}
}

void
AutomaticStructElem::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::INERTIA)) {
			Var_B->put_rec(B.pGetVec(), OH.GetCurrentStep());
			Var_G->put_rec(G.pGetVec(), OH.GetCurrentStep());
			Var_BP->put_rec(BP.pGetVec(), OH.GetCurrentStep());
			Var_GP->put_rec(GP.pGetVec(), OH.GetCurrentStep());
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::INERTIA)) {
			OH.Inertia() << std::setw(8) << GetLabel()
				<< " " << B
				<< " " << G
				<< " " << BP
				<< " " << GP
				<< std::endl;
		}
	}
}

/* Setta i valori iniziali delle variabili (e fa altre cose)
 * prima di iniziare l'integrazione */
void
AutomaticStructElem::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iIndex = pNode->iGetFirstMomentumIndex();

	XP.Put(iIndex + 1, BP);
	XP.Put(iIndex + 4, GP);
}

/* Dati privati */
unsigned int
AutomaticStructElem::iGetNumPrivData(void) const
{
	return 12;
}

unsigned int
AutomaticStructElem::iGetPrivDataIdx(const char *s) const
{
	/*
	 * beta[1]
	 * beta[2]
	 * beta[3]
	 * gamma[1]
	 * gamma[2]
	 * gamma[3]
	 * betaP[1]
	 * betaP[2]
	 * betaP[3]
	 * gammaP[1]
	 * gammaP[2]
	 * gammaP[3]
	 */
	unsigned int idx = 0;
	if (strncmp(s, "beta", STRLENOF("beta")) == 0) {
		s += STRLENOF("beta");
	} else if (strncmp(s, "gamma", STRLENOF("gamma")) == 0) {
		s += STRLENOF("gamma");
		idx += 3;
	} else {
		return 0;
	}

	if (s[0] == 'P') {
		s++;
		idx += 6;
	}

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
		break;

	default:
		return 0;
	}

	if (s[0] != ']' && s[1] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
AutomaticStructElem::dGetPrivData(unsigned int i) const
{
	unsigned int der = (i - 1)/6;
	i -= 6*der;
	unsigned int type = (i - 1)/3;
	i -= 3*type;

	if (der) {
		if (type) {
			return GP(i);
		}
		return BP(i);

	} else {
		if (type) {
			return G(i);
		}
		return B(i);
	}
}

/* AutomaticStructElem - end */

