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

/* Elementi aerodinamici bidimensionali */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "aerodata.h"
#include "aeroelem.h"

#include "shapefnc.h"

#include "aerodc81.h"
#include "c81data.h"
#include "Rot.hh"

#include <sstream>

/* AerodynamicOutput - begin */

AerodynamicOutput::AerodynamicOutput(flag f, int iNP,
		unsigned uFlags, OrientationDescription ood)
: uOutputFlags(uFlags),
od(ood)
{
	SetOutputFlag(f, iNP);
}

AerodynamicOutput::~AerodynamicOutput(void)
{
	NO_OP;
}

void
AerodynamicOutput::SetOutputFlag(flag f, int iNP)
{
	m_eOutput = f;

	if (IsOutput() && IsPGAUSS()) {
		pOutput.resize(iNP);

	} else {
		pOutput.resize(0);
	}
}

void
AerodynamicOutput::ResetIterator(void)
{
	if (IsOutput() && IsPGAUSS()) {
		ASSERT(!pOutput.empty());
		pTmpOutput = pOutput.begin();
	}

#ifdef USE_NETCDF
	if (!pNetCDFOutput.empty()) {
		pTmpNetCDFOutput = pNetCDFOutput.begin();
	}
#endif // USE_NETCDF
}

void
AerodynamicOutput::SetData(const Vec3& v, const doublereal* pd,
	const Vec3& X, const Mat3x3& R, const Vec3& V, const Vec3& W, const Vec3& F, const Vec3& M)
{
	if (IsPGAUSS()) {
		ASSERT(!pOutput.empty());
		ASSERT(pTmpOutput >= pOutput.begin());
		ASSERT(pTmpOutput < pOutput.end());

		pTmpOutput->alpha = 180./M_PI*atan2(-v(2), v(1));
		pTmpOutput->f = Vec3(pd[1], pd[0], pd[5]);

		// move iterator forward
		pTmpOutput++;
	}

#ifdef USE_NETCDF
	if (!pNetCDFOutput.empty()) {
		ASSERT(pTmpNetCDFOutput >= pNetCDFOutput.begin());
		ASSERT(pTmpNetCDFOutput < pNetCDFOutput.end());

		if (pTmpNetCDFOutput->Var_X) pTmpNetCDFOutput->X = X;
		if (pTmpNetCDFOutput->Var_Phi) pTmpNetCDFOutput->R = R;
		if (pTmpNetCDFOutput->Var_V) pTmpNetCDFOutput->V = V;
		if (pTmpNetCDFOutput->Var_W) pTmpNetCDFOutput->W = W;
		if (pTmpNetCDFOutput->Var_F) pTmpNetCDFOutput->F = F;
		if (pTmpNetCDFOutput->Var_M) pTmpNetCDFOutput->M = M;

		pTmpNetCDFOutput++;
	}
#endif // USE_NETCDF
}

AerodynamicOutput::eOutput
AerodynamicOutput::GetOutput(void) const
{
	return eOutput(m_eOutput & AEROD_OUT_MASK);
}

bool
AerodynamicOutput::IsOutput(void) const
{
	return (m_eOutput & 0x1);
}

bool
AerodynamicOutput::IsSTD(void) const
{
	return GetOutput() == AEROD_OUT_STD;
}

bool
AerodynamicOutput::IsPGAUSS(void) const
{
	return GetOutput() == AEROD_OUT_PGAUSS;
}

bool
AerodynamicOutput::IsNODE(void) const
{
	return GetOutput() == AEROD_OUT_NODE;
}

/* AerodynamicOutput - end */


/* Aerodynamic2DElem - begin */

static const bool bDefaultUseJacobian = false;

template <unsigned iNN>
Aerodynamic2DElem<iNN>::Aerodynamic2DElem(unsigned int uLabel,
	const DofOwner *pDO,
	InducedVelocity* pR, bool bPassive,
	const Shape* pC, const Shape* pF,
	const Shape* pV, const Shape* pT,
	const Shape* pTL,
	integer iN, AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	unsigned uFlags, OrientationDescription ood,
	flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut),
DriveOwner(pDC),
AerodynamicOutput(fOut, iNN*iN, uFlags, ood),
aerodata(a),
pIndVel(pR),
bPassiveInducedVelocity(bPassive),
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
TipLoss(pTL),
GDI(iN),
OUTA(iNN*iN, outa_Zero),
bJacobian(bUseJacobian)
{
	DEBUGCOUTFNAME("Aerodynamic2DElem::Aerodynamic2DElem");

	ASSERT(iNN >= 1 && iNN <= 3);
	ASSERT(aerodata != 0);

#ifdef DEBUG
	if (pIndVel != 0) {
		ASSERT(pIndVel->GetElemType() == Elem::INDUCEDVELOCITY);
	}
#endif /* DEBUG */
}

template <unsigned iNN>
Aerodynamic2DElem<iNN>::~Aerodynamic2DElem(void)
{
	DEBUGCOUTFNAME("Aerodynamic2DElem::~Aerodynamic2DElem");

	SAFEDELETE(aerodata);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::SetOutputFlag(flag f)
{
	DEBUGCOUTFNAME("Aerodynamic2DElem::SetOutputFlag");
	ToBeOutput::SetOutputFlag(f);
	AerodynamicOutput::SetOutputFlag(f, iNN*GDI.iGetNum());
}

/* Dimensioni del workspace */
template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = *piNumCols = iNN*6 + iGetNumDof();
}

/* inherited from SimulationEntity */
template <unsigned iNN>
unsigned int
Aerodynamic2DElem<iNN>::iGetNumDof(void) const
{
	return aerodata->iGetNumDof()*iNN*GDI.iGetNum();
}

template <unsigned iNN>
std::ostream&
Aerodynamic2DElem<iNN>::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(bInitial == false);

	integer iFirstIndex = iGetFirstIndex();
	integer iNumDof = aerodata->iGetNumDof();

	for (unsigned b = 0; b < iNN; b++) {
		for (integer i = 0; i < GDI.iGetNum(); i++) {
			out 
				<< prefix << iFirstIndex + 1 << "->" << iFirstIndex + iNumDof << ": "
				"internal states at point #" << i << " of " << GDI.iGetNum() << ", "
				"block #" << b << " of " << iNN << std::endl;

			iFirstIndex += iNumDof;
		}
	}

	return out;
}

static const char *elemnames[] = {
	"AerodynamicBody",
	"AerodynamicBeam2",
	"AerodynamicBeam3",
	0
};

template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	if (i < -1) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(bInitial == false);

	std::ostringstream os;
	os << elemnames[iNN - 1] << "(" << GetLabel() << ")";

	integer iNumDof = aerodata->iGetNumDof();
	if (i == -1) {
		integer iTotDof = iNumDof*GDI.iGetNum();
		desc.resize(iTotDof);

		std::string name(os.str());

		for (integer s = 0; s < iTotDof; s++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": internal state #" << s % iNumDof << " of " << iNumDof
				<< " at point #" << s/iNumDof << " of " << GDI.iGetNum() << ", "
				"block #" << s/(iNumDof*GDI.iGetNum()) << " of " << iNN;
			desc[s] = os.str();
		}

		return;
	}

	desc.resize(1);
	os << ": internal state #" << i % iNumDof << " of " << iNumDof
		<< " at point #" << i/iNumDof << " of " << GDI.iGetNum() << ", "
		"block #" << i/(iNumDof*GDI.iGetNum()) << " of " << iNN;
	desc[0] = os.str();
}

template <unsigned iNN>
std::ostream&
Aerodynamic2DElem<iNN>::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(bInitial == false);

	integer iFirstIndex = iGetFirstIndex();
	integer iNumDof = aerodata->iGetNumDof();

	for (unsigned b = 0; b < iNN; b++) {
		for (integer i = 0; i < GDI.iGetNum(); i++) {
			out 
				<< prefix << iFirstIndex + 1 << "->" << iFirstIndex + iNumDof
				<< ": dynamics equations at point #" << i << " of " << GDI.iGetNum() << ", "
				"block #" << b << " of " << iNN << std::endl;

			iFirstIndex += iNumDof;
		}
	}

	return out;
}

template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	if (i < -1) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(bInitial == false);

	std::ostringstream os;
	os << elemnames[iNN - 1] << "(" << GetLabel() << ")";

	integer iNumDof = aerodata->iGetNumDof();
	if (i == -1) {
		integer iTotDof = iNumDof*GDI.iGetNum();
		desc.resize(iTotDof);

		std::string name(os.str());

		for (integer s = 0; s < iTotDof; s++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": internal equation #" << s % iNumDof << " of " << iNumDof
				<< " at point #" << s/iNumDof << " of " << GDI.iGetNum() << ", "
				"block " << s/(iNumDof*GDI.iGetNum()) << " of " << iNN;
			desc[s] = os.str();
		}

		return;
	}

	desc.resize(1);
	os << ": internal equation #" << i % iNumDof << " of " << iNumDof
		<< " at point #" << i/iNumDof << " of " << GDI.iGetNum() << ", "
		"block " << i/(iNumDof*GDI.iGetNum()) << " of " << iNN;
	desc[0] = os.str();
}

template <unsigned iNN>
DofOrder::Order
Aerodynamic2DElem<iNN>::GetDofType(unsigned int i) const
{
	ASSERT(aerodata->iGetNumDof() > 0);

	return aerodata->GetDofType(i % aerodata->iGetNumDof());
}

template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints* h)
{
	integer iNumDof = aerodata->iGetNumDof();
	if (iNumDof) {
		vx.Resize(6*iNN);
		wx.Resize(6*iNN);
		fq.Resize(iNumDof);
		cq.Resize(iNumDof);

		vx.Reset();
		wx.Reset();
		fq.Reset();
		cq.Reset();
	}
}

template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::AfterConvergence(
	const VectorHandler& X,
	const VectorHandler& XP)
{
	/* Memoria in caso di forze instazionarie */
	switch (aerodata->Unsteady()) {
	case AeroData::STEADY:
		break;

	case AeroData::HARRIS:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case AeroData::BIELAWA:
		for (unsigned i = 0; i < iNN*GDI.iGetNum(); i++) {
	 		aerodata->Update(i);
		}
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (unsigned i = 0; i < iNN*GDI.iGetNum(); i++) {
		aerodata->AfterConvergence(i, X, XP);
	}
}

/* Dimensioni del workspace */
template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*iNN;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
template <unsigned iNN>
VariableSubMatrixHandler&
Aerodynamic2DElem<iNN>::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */)
{
	DEBUGCOUTFNAME("Aerodynamic2DElem::InitialAssJac");
	WorkMat.SetNullMatrix();
	return WorkMat;
}

template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::OutputPrepare(OutputHandler &OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::AERODYNAMIC)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			/* get a pointer to binary NetCDF file
			 * -->  pDM->OutHdl.BinFile */
			NcFile *pBinFile = OH.pGetBinFile();
			char buf[BUFSIZ];

			int l = snprintf(buf, sizeof(buf), "elem.aerodynamic.%lu",
				(unsigned long)GetLabel());

			// X_XX
			// R_XX
			// Phi_XX
			// V_XX
			// Omega_XX
			// F_XX
			// M_XX
			
			// NOTE: "Omega_XX" is the longest var name
			if (l < 0 || l >= int(sizeof(buf) - STRLENOF(".Omega_XX"))) {
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			NcVar *Var_Type = pBinFile->add_var(buf, ncChar, OH.DimV1());

			if (!Var_Type->add_att("type", elemnames[iNN - 1])) {
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			// add var name separator
			buf[l++] = '.';

			int totgp = iNN*GDI.iGetNum();
			pNetCDFOutput.resize(totgp);

			int j = 0;
			for (std::vector<AeroNetCDFOutput>::iterator i = pNetCDFOutput.begin();
				i != pNetCDFOutput.end(); i++, j++)
			{
				/* Add NetCDF (output) variables to the BinFile object
				 * and save the NcVar* pointer returned from add_var
				 * as handle for later write accesses.
				 * Define also variable attributes */
				i->Var_X = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_X) {
					strcpy(&buf[l], "X_");
					snprintf(&buf[l + STRLENOF("X_")], sizeof(buf) - l - STRLENOF("X_"), "%d", j);
					i->Var_X = pBinFile->add_var(buf, ncDouble,
						OH.DimTime(), OH.DimV3());
					if (i->Var_X == 0) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_X->add_att("units", "m")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_X->add_att("type", "Vec3")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					char descbuf[BUFSIZ];
					snprintf(descbuf, sizeof(descbuf),
						"global position vector (X, Y, Z) of Gauss point #%d/%d", j, totgp);
					if (!i->Var_X->add_att("description", descbuf)) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				i->Var_Phi = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_R) {
					char descbuf[BUFSIZ];

					switch (od) {
					case ORIENTATION_MATRIX:
						strcpy(&buf[l], "R_");
						snprintf(&buf[l + STRLENOF("R_")],
							sizeof(buf) - l - STRLENOF("R_"), "%d", j);
						i->Var_Phi = pBinFile->add_var(buf, ncDouble,
							OH.DimTime(), OH.DimV3(), OH.DimV3());
						if (i->Var_Phi == 0) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("units", "-")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("type", "Mat3x3")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						snprintf(descbuf, sizeof(descbuf),
							"global orientation matrix (R11, R21, R31, R12, R22, R32, R13, R23, R33) of Gauss point #%d/%d", j, totgp);
						if (!i->Var_Phi->add_att("description", descbuf)) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						break;

					case ORIENTATION_VECTOR:
						strcpy(&buf[l], "Phi_");
						snprintf(&buf[l + STRLENOF("Phi_")],
							sizeof(buf) - l - STRLENOF("Phi_"), "%d", j);
						i->Var_Phi = pBinFile->add_var(buf, ncDouble,
							OH.DimTime(), OH.DimV3());
						if (i->Var_Phi == 0) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("units", "radian")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("type", "Vec3")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						snprintf(descbuf, sizeof(descbuf),
							"global orientation vector (Phi_X, Phi_Y, Phi_Z) of Gauss point #%d/%d", j, totgp);
						if (!i->Var_Phi->add_att("description", descbuf)) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						break;

					case EULER_123:
					case EULER_313:
					case EULER_321:
						{
						strcpy(&buf[l], "E_");
						snprintf(&buf[l + STRLENOF("E_")],
							sizeof(buf) - l - STRLENOF("E_"), "%d", j);
						i->Var_Phi = pBinFile->add_var(buf, ncDouble,
							OH.DimTime(), OH.DimV3());
						if (i->Var_Phi == 0) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("units", "radian")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (!i->Var_Phi->add_att("type", "Vec3")) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						std::string desc;
						switch (od) {
						case EULER_123:
							desc = "global orientation Euler angles (123) (E_X, E_Y, E_Z) ";
							break;

						case EULER_313:
							desc = "global orientation Euler angles (313) (E_Z, E_X, E_Z') ";
							break;

						case EULER_321:
							desc = "global orientation Euler angles (321) (E_Z, E_Y, E_X) ";
							break;

						default:
							ASSERT(0);
							break;
						}

						char pgbuf[sizeof("#XX/XX")];
						snprintf(pgbuf, sizeof(pgbuf), "#%d/%d", j, totgp);
						desc += pgbuf;

						if (!i->Var_Phi->add_att("description", desc.c_str())) {
							throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						} break;

					default:
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				i->Var_V = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_V) {
					strcpy(&buf[l], "V_");
					snprintf(&buf[l + STRLENOF("V_")], sizeof(buf) - l - STRLENOF("V_"), "%d", j);
					i->Var_V = pBinFile->add_var(buf, ncDouble,
						OH.DimTime(), OH.DimV3());
					if (i->Var_V == 0) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_V->add_att("units", "m/s")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_V->add_att("type", "Vec3")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					char descbuf[BUFSIZ];
					snprintf(descbuf, sizeof(descbuf),
						"velocity in global frame (F_X, F_Y, F_Z) of Gauss point #%d/%d", j, totgp);
					if (!i->Var_V->add_att("description", descbuf)) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				i->Var_W = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_W) {
					strcpy(&buf[l], "Omega_");
					snprintf(&buf[l + STRLENOF("Omega_")], sizeof(buf) - l - STRLENOF("Omega_"), "%d", j);
					i->Var_W = pBinFile->add_var(buf, ncDouble,
						OH.DimTime(), OH.DimV3());
					if (i->Var_W == 0) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_W->add_att("units", "radian/s")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_W->add_att("type", "Vec3")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					char descbuf[BUFSIZ];
					snprintf(descbuf, sizeof(descbuf),
						"angular velocity in global frame (F_X, F_Y, F_Z) of Gauss point #%d/%d", j, totgp);
					if (!i->Var_W->add_att("description", descbuf)) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				i->Var_F = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_F) {
					strcpy(&buf[l], "F_");
					snprintf(&buf[l + STRLENOF("F_")], sizeof(buf) - l - STRLENOF("F_"), "%d", j);
					i->Var_F = pBinFile->add_var(buf, ncDouble,
						OH.DimTime(), OH.DimV3());
					if (i->Var_F == 0) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_F->add_att("units", "N")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_F->add_att("type", "Vec3")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					char descbuf[BUFSIZ];
					snprintf(descbuf, sizeof(descbuf),
						"force in local frame (F_X, F_Y, F_Z) of Gauss point #%d/%d", j, totgp);
					if (!i->Var_F->add_att("description", descbuf)) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				i->Var_M = 0;
				if (uOutputFlags & AerodynamicOutput::OUTPUT_GP_M) {
					strcpy(&buf[l], "M_");
					snprintf(&buf[l + STRLENOF("M_")], sizeof(buf) - l - STRLENOF("M_"), "%d", j);
					i->Var_M = pBinFile->add_var(buf, ncDouble,
						OH.DimTime(), OH.DimV3());
					if (i->Var_M == 0) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_M->add_att("units", "Nm")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (!i->Var_M->add_att("type", "Vec3")) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					char descbuf[BUFSIZ];
					snprintf(descbuf, sizeof(descbuf),
						"moment in local frame (M_X, M_Y, M_Z) of Gauss point #%d/%d", j, totgp);
					if (!i->Var_M->add_att("description", descbuf)) {
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}
		}
#endif // USE_NETCDF
	}
}

/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::Output_int(OutputHandler &OH) const
{
#ifdef USE_NETCDF
	if (OH.UseNetCDF(OutputHandler::AERODYNAMIC)) {
		for (std::vector<AeroNetCDFOutput>::const_iterator i = pNetCDFOutput.begin();
			i != pNetCDFOutput.end(); i++)
		{
			if (i->Var_X) {
				i->Var_X->put_rec(i->X.pGetVec(), OH.GetCurrentStep());
			}

			if (i->Var_Phi) {
				Vec3 E;
				switch (od) {
				case EULER_123:
					E = MatR2EulerAngles123(i->R)*dRaDegr;
					break;

				case EULER_313:
					E = MatR2EulerAngles313(i->R)*dRaDegr;
					break;

				case EULER_321:
					E = MatR2EulerAngles321(i->R)*dRaDegr;
					break;

				case ORIENTATION_VECTOR:
					E = RotManip::VecRot(i->R);
					break;

				case ORIENTATION_MATRIX:
					break;

				default:
					/* impossible */
					break;
				}

				switch (od) {
				case EULER_123:
				case EULER_313:
				case EULER_321:
				case ORIENTATION_VECTOR:
					i->Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
					break;

				case ORIENTATION_MATRIX:
					i->Var_Phi->put_rec(i->R.pGetMat(), OH.GetCurrentStep());
					break;

				default:
					/* impossible */
					break;
				}
			}

			if (i->Var_V) {
				i->Var_V->put_rec(i->V.pGetVec(), OH.GetCurrentStep());
			}

			if (i->Var_W) {
				i->Var_W->put_rec(i->W.pGetVec(), OH.GetCurrentStep());
			}

			if (i->Var_F) {
				i->Var_F->put_rec(i->F.pGetVec(), OH.GetCurrentStep());
			}

			if (i->Var_M) {
				i->Var_M->put_rec(i->M.pGetVec(), OH.GetCurrentStep());
			}
		}
	}
#endif /* USE_NETCDF */
}

// only send forces if:
// 1) an induced velocity model is defined
// 2) this element is not "passive" (i.e. contributes to induced velocity)
// 3) the induced velocity model does not require sectional forces
template <unsigned iNN>
void 
Aerodynamic2DElem<iNN>::AddForce_int(const Vec3& F,
	const Vec3& M, const Vec3& X) const
{
	if (pIndVel != 0 && !bPassiveInducedVelocity
		&& !pIndVel->bSectionalForces())
	{
		pIndVel->AddForce(GetLabel(), F, M, X);
	}
}

// only send forces if:
// 1) an induced velocity model is defined
// 2) this element is not "passive" (i.e. contributes to induced velocity)
// 3) the induced velocity model requires sectional forces
template <unsigned iNN>
void
Aerodynamic2DElem<iNN>::AddSectionalForce_int(unsigned uPnt,
	const Vec3& F, const Vec3& M, doublereal dW,
	const Vec3& X, const Mat3x3& R,
	const Vec3& V, const Vec3& W) const
{
	if (pIndVel != 0 && !bPassiveInducedVelocity
		&& pIndVel->bSectionalForces())
	{
		pIndVel->AddSectionalForce(GetElemType(), GetLabel(), uPnt,
			F, M, dW, X, R, V, W);
	}
}

/* Aerodynamic2DElem - end */


/* AerodynamicBody - begin */

AerodynamicBody::AerodynamicBody(unsigned int uLabel,
	const DofOwner *pDO,
	const StructNode* pN, InducedVelocity* pR, bool bPassive,
	const Vec3& fTmp, doublereal dS,
	const Mat3x3& RaTmp,
	const Shape* pC, const Shape* pF,
	const Shape* pV, const Shape* pT,
	const Shape* pTL,
	integer iN, AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	unsigned uFlags, OrientationDescription ood,
	flag fOut)
: Elem(uLabel, fOut),
Aerodynamic2DElem<1>(uLabel, pDO, pR, bPassive, pC, pF, pV, pT, pTL, iN,
	a, pDC, bUseJacobian, uFlags, ood, fOut),
pNode(pN),
f(fTmp),
dHalfSpan(dS/2.),
Ra(RaTmp),
Ra3(RaTmp.GetVec(3)),
F(0.),
M(0.)
{
	DEBUGCOUTFNAME("AerodynamicBody::AerodynamicBody");

	ASSERT(pNode != 0);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);

}

AerodynamicBody::~AerodynamicBody(void)
{
	DEBUGCOUTFNAME("AerodynamicBody::~AerodynamicBody");
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBody::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBody::Restart");

	out << "  aerodynamic body: " << GetLabel() << ", "
		<< pNode->GetLabel();
	if (pIndVel != 0) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f.Write(out, ", ")
		<< ", " << dHalfSpan*2.
		<< ", reference, node, 1, ", (Ra.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< ", " << GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

VariableSubMatrixHandler&
AerodynamicBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering AerodynamicBody::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	/* Equations start here... */
	doublereal dW[6];

	/* Dati del nodo */
	const Vec3& Xn(pNode->GetXCurr());
	const Mat3x3& Rn(pNode->GetRCurr());
	const Vec3& Vn(pNode->GetVCurr());
	const Vec3& Wn(pNode->GetWCurr());

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR(Rn*Ra);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso la posizione del nodo perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 6;

		integer iOffset = iFirstEq - 6;

		for (int iCnt = 6 + 1; iCnt <= iNumRows; iCnt++) {
			WM.PutRowIndex(iCnt, iOffset + iCnt);
			WM.PutColIndex(iCnt, iOffset + iCnt);
		}
	}

	/* Ciclo sui punti di Gauss */
	PntWght PW = GDI.GetFirst();
	int iPnt = 0;
	do {
		doublereal dCsi = PW.dGetPnt();
		Vec3 Xr(Rn*(f + Ra3*(dHalfSpan*dCsi)));
		Vec3 Xnr = Xn + Xr;
		Vec3 Vr(Vn + Wn.Cross(Xr));

		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xnr)) {
			Vr -= VTmp;
		}

		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pIndVel != 0) {
			Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xnr);
		}

		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(dCsi) + dGet();
		aerodata->SetSectionData(dCsi,
			Chord.dGet(dCsi),
			ForcePoint.dGet(dCsi),
			VelocityPoint.dGet(dCsi),
			dTw,
			dOmega);

		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		Mat3x3 RRloc;
		if (dTw != 0.) {
			doublereal dCosT = cos(dTw);
			doublereal dSinT = sin(dTw);
			/* Assumo lo svergolamento positivo a cabrare */
			Mat3x3 RTw( dCosT, dSinT, 0.,
				-dSinT, dCosT, 0.,
				0.,    0.,    1.);
			RRloc = RR*RTw;

		} else {
			RRloc = RR;
		}

		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		VTmp = RRloc.MulTV(Vr);
		VTmp.PutTo(dW);

		Vec3 WTmp = RRloc.MulTV(Wn);
		WTmp.PutTo(&dW[3]);

		/* Funzione di calcolo delle forze aerodinamiche */
		doublereal Fa0[6];
		Mat6x6 JFa;

		/* Jacobian Assembly... */
		doublereal cc = dHalfSpan*PW.dGetWght();

		if (iNumDof) {
			// prepare (v/dot{x} + dCoef*v/x) and so
			Mat3x3 RRlocT(RRloc.Transpose());

			vx.PutMat3x3(1, RRlocT);
			vx.PutMat3x3(4, RRloc.MulTM(Mat3x3((Vr + Xr.Cross(Wn))*dCoef - Xr)));
			wx.PutMat3x3(4, RRlocT);

			// equations from iFirstEq on are dealt with by aerodata
			aerodata->AssJac(WM, dCoef, XCurr, XPrimeCurr,
                		iFirstEq, iFirstSubEq,
				vx, wx, fq, cq, iPnt, dW, Fa0, JFa, OUTA[iPnt]);

			// deal with (f/dot{q} + dCoef*f/q) and so
			integer iOffset = 6 + iPnt*iNumDof;
			for (integer iCol = 1; iCol <= iNumDof; iCol++) {
				Vec3 fqTmp((RRloc*fq.GetVec(iCol))*cc);
				Vec3 cqTmp(Xr.Cross(fqTmp) + (RRloc*cq.GetVec(iCol))*cc);

				WM.Sub(1, iOffset + iCol, fqTmp);
				WM.Sub(4, iOffset + iCol, cqTmp);
			}

			// first equation
			iFirstEq += iNumDof;
			iFirstSubEq += iNumDof;

		} else {
			aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, OUTA[iPnt]);
		}

		// rotate force, couple and Jacobian matrix in absolute frame
		Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

		Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
		Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + Xr.Cross(fTmp));

		// compute submatrices (see tecman.pdf)
		Mat3x3 Mat21(Xr.Cross(JFaR.GetMat11()) + JFaR.GetMat21());

		Mat3x3 Mat12(JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(Xr));

		Mat3x3 Mat22(Xr.Cross(JFaR.GetMat12()) + JFaR.GetMat22() - Mat21*Mat3x3(Xr));

		Mat3x3 MatV((Vr + Xr.Cross(Wn))*dCoef);

		Mat12 += JFaR.GetMat11()*MatV - Mat3x3(fTmp);

		Mat22 += Mat21*MatV - Mat3x3(cTmp);

		// add (actually, sub) contributions, scaled by weight
		WM.Sub(1, 1, JFaR.GetMat11());
		WM.Sub(1, 4, Mat12);
		WM.Sub(4, 1, Mat21);
		WM.Sub(4, 4, Mat22);

		iPnt++;

	} while (GDI.fGetNext(PW));

	return WorkMat;
}


SubVectorHandler&
AerodynamicBody::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::AssRes");

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstIndex = pNode->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	AssVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;
}


SubVectorHandler&
AerodynamicBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::InitialAssRes");

	if (aerodata->iGetNumDof() > 0) {
		WorkVec.ResizeReset(0);
		return WorkVec;
	}

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstIndex = pNode->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	AssVec(WorkVec, 1., XCurr, XCurr);

	return WorkVec;
}


/* assemblaggio residuo */
void
AerodynamicBody::AssVec(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::AssVec");

	doublereal dTng[6];
	doublereal dW[6];

	/* Dati del nodo */
	const Vec3& Xn(pNode->GetXCurr());
	const Mat3x3& Rn(pNode->GetRCurr());
	const Vec3& Vn(pNode->GetVCurr());
	const Vec3& Wn(pNode->GetWCurr());

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR(Rn*Ra);

	/*
	 * Se l'elemento e' collegato ad un rotore,
 	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/* Resetta i dati */
	F.Reset();
	M.Reset();

	/*
	 * Dati "permanenti" (uso la posizione del nodo perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 6;

		integer iOffset = iFirstEq - 6;
		integer iNumRows = WorkVec.iGetSize();
		for (int iCnt = 6 + 1; iCnt <= iNumRows; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iOffset + iCnt);
		}
	}

	/* Ciclo sui punti di Gauss */
	PntWght PW = GDI.GetFirst();
	int iPnt = 0;
	do {
		doublereal dCsi = PW.dGetPnt();
		Vec3 Xr(Rn*(f + Ra3*(dHalfSpan*dCsi)));
		Vec3 Xnr = Xn + Xr;
		Vec3 Vr(Vn + Wn.Cross(Xr));

		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xnr)) {
	 		Vr -= VTmp;
		}

		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pIndVel != 0) {
	 		Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xnr);
		}

		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(dCsi) + dGet();
		aerodata->SetSectionData(dCsi,
			Chord.dGet(dCsi),
			ForcePoint.dGet(dCsi),
			VelocityPoint.dGet(dCsi),
			dTw,
			dOmega);

		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		Mat3x3 RRloc;
		if (dTw != 0.) {
			doublereal dCosT = cos(dTw);
			doublereal dSinT = sin(dTw);
			/* Assumo lo svergolamento positivo a cabrare */
			Mat3x3 RTw( dCosT, dSinT, 0.,
				   -dSinT, dCosT, 0.,
				    0.,    0.,    1.);
			RRloc = RR*RTw;

		} else {
			RRloc = RR;
		}

		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		VTmp = RRloc.MulTV(Vr);
		VTmp.PutTo(dW);

		Vec3 WTmp = RRloc.MulTV(Wn);
		WTmp.PutTo(&dW[3]);

		/* Funzione di calcolo delle forze aerodinamiche */
		if (iNumDof) {
			aerodata->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr,
                		iFirstEq, iFirstSubEq, iPnt, dW, dTng, OUTA[iPnt]);

			// first equation
			iFirstEq += iNumDof;
			iFirstSubEq += iNumDof;

		} else {
			aerodata->GetForces(iPnt, dW, dTng, OUTA[iPnt]);
		}

		/* Dimensionalizza le forze */
		doublereal dWght = dHalfSpan*PW.dGetWght();
		dTng[1] *= TipLoss.dGet(dCsi);
		Vec3 FTmp(RRloc*(Vec3(&dTng[0])));
		Vec3 MTmp(RRloc*(Vec3(&dTng[3])));

		// Se e' definito il rotore, aggiungere il contributo alla trazione
		AddSectionalForce_int(iPnt, FTmp, MTmp, dWght, Xnr, RRloc, Vr, Wn);

		FTmp *= dWght;
		MTmp *= dWght;

		F += FTmp;
		M += MTmp;
		M += Xr.Cross(FTmp);

		// specific for Gauss points force output
		if (fToBeOutput()) {
			SetData(VTmp, dTng, Xr, RRloc, Vr, Wn, FTmp, MTmp);
		}

		iPnt++;

	} while (GDI.fGetNext(PW));

	// Se e' definito il rotore, aggiungere il contributo alla trazione
	AddForce_int(F, M, Xn);

	/* Sommare il termine al residuo */
	WorkVec.Add(1, F);
	WorkVec.Add(4, M);
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBody::Output(OutputHandler& OH) const
{
	/* Output delle forze aerodinamiche F, M su apposito file */
	if (fToBeOutput()) {
		Aerodynamic2DElem<1>::Output_int(OH);

		if (OH.UseText(OutputHandler::AERODYNAMIC)) {
			std::ostream& out = OH.Aerodynamic()
				<< std::setw(8) << GetLabel();

			switch (GetOutput()) {
			case AEROD_OUT_NODE:
				out << " " << std::setw(8) << pNode->GetLabel()
					<< " ", F.Write(out, " ") << " ", M.Write(out, " ");
				break;

			case AEROD_OUT_PGAUSS:
				ASSERT(!pOutput.empty());
				for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
					i != pOutput.end(); i++)
				{
	 				out << " " << i->alpha
						<< " " << i->f;
				}
				break;

			case AEROD_OUT_STD:
				for (int i = 0; i < GDI.iGetNum(); i++) {
	 				out
						<< " " << OUTA[i].alpha
						<< " " << OUTA[i].gamma
						<< " " << OUTA[i].mach
						<< " " << OUTA[i].cl
						<< " " << OUTA[i].cd
						<< " " << OUTA[i].cm
						<< " " << OUTA[i].alf1
						<< " " << OUTA[i].alf2;
				}
				break;

			default:
				ASSERT(0);
				break;
			}

			out << std::endl;
	}	}
}

/* AerodynamicBody - end */

static bool
ReadInducedVelocity(DataManager *pDM, MBDynParser& HP,
	unsigned uLabel, const char *sElemType,
	InducedVelocity*& pIndVel, bool &bPassive)
{
	bool bReadIV(false);
	bool bReadUDIV(false);
	if (HP.IsKeyWord("rotor")) {
		silent_cerr(sElemType << "(" << uLabel << "): "
			"\"rotor\" keyword is deprecated; "
			"use \"induced velocity\" instead "
			"at line " << HP.GetLineData()
			<< std::endl);

		bReadIV = true;

	} else if (HP.IsKeyWord("induced" "velocity")) {
		bReadIV = true;

	} else if (HP.IsKeyWord("user" "defined" "induced" "velocity")) {
		bReadIV = true;
		bReadUDIV = true;
	}

	if (bReadIV) {
		unsigned int uIV = (unsigned int)HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Linked to InducedVelocity(" << uIV << ")" << std::endl);

		bPassive = false;
		if (HP.IsKeyWord("passive")) {
			bPassive = true;
		}

		/*
		 * verifica di esistenza del rotore
		 * NOTA: ovviamente il rotore deve essere definito
		 * prima dell'elemento aerodinamico
		 */
		Elem* p;

		if (bReadUDIV) {
			p = pDM->pFindElem(Elem::LOADABLE, uIV);
			if (p == 0) {
				silent_cerr(sElemType << "(" << uLabel << "): "
					"user-defined InducedVelocity(" << uIV << ") not defined "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else {
	 		p = pDM->pFindElem(Elem::INDUCEDVELOCITY, uIV);
			if (p == 0) {
				// try a user-defined one?
				p = pDM->pFindElem(Elem::LOADABLE, uIV);
				if (p == 0 || !dynamic_cast<InducedVelocity *>(p)) {
					silent_cerr(sElemType << "(" << uLabel << "): "
						"InducedVelocity(" << uIV << ") not defined "
						"at line " << HP.GetLineData()
						<< std::endl);

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				silent_cerr(sElemType << "(" << uLabel << "): "
					"InducedVelocity(" << uIV << ") not defined; using user-defined InducedVelocity(" << uIV << ") "
					"at line " << HP.GetLineData()
					<< std::endl);
			}
		}

		pIndVel = dynamic_cast<InducedVelocity *>(p);
		ASSERT(pIndVel != 0);
	}

	return (pIndVel != 0);
}

void
ReadAerodynamicCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	unsigned& uFlags, OrientationDescription& od)
{
	uFlags = AerodynamicOutput::OUTPUT_NONE;

	while (HP.IsArg()) {
		unsigned uFlag;

		if (HP.IsKeyWord("position")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_X;

		} else if (HP.IsKeyWord("orientation")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_R;

		} else if (HP.IsKeyWord("velocity")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_V;

		} else if (HP.IsKeyWord("angular" "velocity")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_W;

		} else if (HP.IsKeyWord("configuration")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_CONFIGURATION;

		} else if (HP.IsKeyWord("force")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_F;

		} else if (HP.IsKeyWord("moment")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_M;

		} else if (HP.IsKeyWord("forces")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_FORCES;

		} else if (HP.IsKeyWord("all")) {
			uFlag = AerodynamicOutput::OUTPUT_GP_ALL;

		} else {
			break;
		}

		if (uFlags & uFlag) {
			silent_cerr("AerodynamicElement(" << uLabel << "): "
				"duplicate custom output "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (uFlag & AerodynamicOutput::OUTPUT_GP_R) {
			od = ReadOptionalOrientationDescription(pDM, HP);
		}

		uFlags |= uFlag;
	}
}

void
ReadOptionalAerodynamicCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	unsigned& uFlags, OrientationDescription& od)
{
	pDM->GetOutput(Elem::AERODYNAMIC, uFlags, od);
	if (HP.IsKeyWord("custom" "output")) {
		ReadAerodynamicCustomOutput(pDM, HP, uLabel, uFlags, od);
	}
}

Elem *
ReadAerodynamicBody(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBody");

	/* Nodo */
	StructNode* pNode = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));

	InducedVelocity* pIndVel = 0;
	bool bPassive(false);
	(void)ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBody",
		pIndVel, bPassive);

	ReferenceFrame RF(pNode);
	Vec3 f(HP.GetPosRel(RF));

	DEBUGLCOUT(MYDEBUG_INPUT, "Offset: " << f << std::endl);

	Mat3x3 Ra(HP.GetRotRel(RF));

	doublereal dSpan = HP.GetReal();
	DEBUGLCOUT(MYDEBUG_INPUT, "Span: " << dSpan << std::endl);

	Shape* pChord = 0;
	Shape* pForce = 0;
	Shape* pVelocity = 0;
	Shape* pTwist = 0;
	Shape* pTipLoss = 0;

	integer iNumber = 0;
	DriveCaller* pDC = 0;
	AeroData* aerodata = 0;

	ReadAeroData(pDM, HP, 1,
		&pChord, &pForce, &pVelocity, &pTwist, &pTipLoss,
		&iNumber, &pDC, &aerodata);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = HP.GetYesNoOrBool(bDefaultUseJacobian);
	}

	if (aerodata->iGetNumDof() > 0 && !bUseJacobian) {
		silent_cerr("AerodynamicBody(" << uLabel << "): "
			"aerodynamic model needs \"jacobian, yes\" at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	unsigned uFlags = AerodynamicOutput::OUTPUT_NONE;
	ReadOptionalAerodynamicCustomOutput(pDM, HP, uLabel, uFlags, od);

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBody(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBody,
		AerodynamicBody(uLabel, pDO, pNode, pIndVel, bPassive,
			f, dSpan, Ra,
			pChord, pForce, pVelocity, pTwist, pTipLoss,
			iNumber, aerodata, pDC, bUseJacobian, uFlags, od, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 Ra3 = Ra.GetVec(3);
	doublereal dCm1 = pChord->dGet(-1.);
	doublereal dPm1 = pForce->dGet(-1.);
	doublereal dTm1 = pTwist->dGet(-1.);
	Vec3 Ram1 = Ra*Vec3(std::cos(dTm1), std::sin(dTm1), 0.);

	doublereal dCp1 = pChord->dGet(1.);
	doublereal dPp1 = pForce->dGet(1.);
	doublereal dTp1 = pTwist->dGet(1.);
	Vec3 Rap1 = Ra*Vec3(std::cos(dTp1), std::sin(dTp1), 0.);

	std::ostream& out = pDM->GetLogFile();
	out << "aero0: " << uLabel
		<< " " << pNode->GetLabel()
		<< " ", (f - Ra3*(dSpan/2.) + Ram1*(dPm1 - dCm1*3./4.)).Write(out, " ")
		<< " ", (f - Ra3*(dSpan/2.) + Ram1*(dPm1 + dCm1/4.)).Write(out, " ")
		<< " ", (f + Ra3*(dSpan/2.) + Rap1*(dPp1 - dCp1*3./4.)).Write(out, " ")
		<< " ", (f + Ra3*(dSpan/2.) + Rap1*(dPp1 + dCp1/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of ReadAerodynamicBody() */


/* AerodynamicBeam - begin */

AerodynamicBeam::AerodynamicBeam(unsigned int uLabel,
	const DofOwner *pDO,
	const Beam* pB, InducedVelocity* pR, bool bPassive,
	const Vec3& fTmp1,
	const Vec3& fTmp2,
	const Vec3& fTmp3,
	const Mat3x3& Ra1Tmp,
	const Mat3x3& Ra2Tmp,
	const Mat3x3& Ra3Tmp,
	const Shape* pC, const Shape* pF,
	const Shape* pV, const Shape* pT,
	const Shape* pTL,
	integer iN, AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	unsigned uFlags, OrientationDescription ood,
	flag fOut)
: Elem(uLabel, fOut),
Aerodynamic2DElem<3>(uLabel, pDO, pR, bPassive,
	pC, pF, pV, pT, pTL, iN, a, pDC, bUseJacobian, uFlags, ood, fOut),
pBeam(pB),
f1(fTmp1),
f2(fTmp2),
f3(fTmp3),
Ra1(Ra1Tmp),
Ra2(Ra2Tmp),
Ra3(Ra3Tmp),
Ra1_3(Ra1Tmp.GetVec(3)),
Ra2_3(Ra2Tmp.GetVec(3)),
Ra3_3(Ra3Tmp.GetVec(3))
{
	DEBUGCOUTFNAME("AerodynamicBeam::AerodynamicBeam");

	ASSERT(pBeam != 0);
	ASSERT(pBeam->GetElemType() == Elem::BEAM);

	pNode1 = pBeam->pGetNode(1);
	pNode2 = pBeam->pGetNode(2);
	pNode3 = pBeam->pGetNode(3);

	ASSERT(pNode1 != 0);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != 0);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode3 != 0);
	ASSERT(pNode3->GetNodeType() == Node::STRUCTURAL);
}

AerodynamicBeam::~AerodynamicBeam(void)
{
	DEBUGCOUTFNAME("AerodynamicBeam::~AerodynamicBeam");
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBeam::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Restart");
	out << "  aerodynamic beam: " << GetLabel()
		<< ", " << pBeam->GetLabel();
	if (pIndVel != 0) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f1.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra1.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra1.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f2.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra2.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra2.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f3.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra3.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra3.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

static const doublereal d13 = 1./sqrt(3.);
static const doublereal pdsi3[] = { -1., -d13, d13 };
static const doublereal pdsf3[] = { -d13, d13, 1. };

/* Jacobian assembly */
VariableSubMatrixHandler&
AerodynamicBeam::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering AerodynamicBeam::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstMomentumIndex();
	integer iNode3FirstPosIndex = pNode3->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
	}

	doublereal dW[6];

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[3];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Xn[NODE3] = pNode3->GetXCurr();
	const Mat3x3& Rn3(pNode3->GetRCurr());
	const Vec3& Vn3(pNode3->GetVCurr());
	const Vec3& Wn3(pNode3->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);
	Vec3 f3Tmp(Rn3*f3);

	Vec3 fTmp[3];
	fTmp[NODE1] = f1Tmp;
	fTmp[NODE2] = f2Tmp;
	fTmp[NODE3] = f3Tmp;

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);
	Vec3 X3Tmp(Xn[NODE3] + f3Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));
	Vec3 V3Tmp(Vn3 + Wn3.Cross(f3Tmp));

	Vec3 Omega1Crossf1(Wn1.Cross(f1Tmp));
	Vec3 Omega2Crossf2(Wn2.Cross(f2Tmp));
	Vec3 Omega3Crossf3(Wn3.Cross(f3Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);
	Mat3x3 RR3(Rn3*Ra3);

	/*
	 * Parametri di rotazione dai nodi 1 e 3 al nodo 2 (nell'ipotesi
	 * che tale trasformazione non dia luogo ad una singolarita')
	 */

	Vec3 g1(ER_Rot::Param, RR2.MulTM(RR1));
	Vec3 g3(ER_Rot::Param, RR2.MulTM(RR3));

	Mat3x3 GammaInv1(ER_Rot::MatGm1, g1);
	Mat3x3 GammaInv3(ER_Rot::MatGm1, g3);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn[NODE2], rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	int iPnt = 0;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 18;

		integer iOffset = iFirstEq - 18;

		for (int iCnt = 18 + 1; iCnt <= iNumRows; iCnt++) {
			WM.PutRowIndex(iCnt, iOffset + iCnt);
			WM.PutColIndex(iCnt, iOffset + iCnt);
		}
	}

	for (int iNode = 0; iNode < LASTNODE; iNode++) {
		doublereal dsi = pdsi3[iNode];
		doublereal dsf = pdsf3[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		Mat3x3 WM_F[6], WM_M[6];
		WM_F[DELTAx1].Reset();
		WM_F[DELTAg1].Reset();
		WM_F[DELTAx2].Reset();
		WM_F[DELTAg2].Reset();
		WM_F[DELTAx3].Reset();
		WM_F[DELTAg3].Reset();
		WM_M[DELTAx1].Reset();
		WM_M[DELTAg1].Reset();
		WM_M[DELTAx2].Reset();
		WM_M[DELTAg2].Reset();
		WM_M[DELTAx3].Reset();
		WM_M[DELTAg3].Reset();

		//unsigned int iDelta_x1, iDelta_g1, iDelta_x2, iDelta_g2, iDelta_x3, iDelta_g3;
		//iDelta_x1 = 0;, iDelta_g1, iDelta_x2, iDelta_g2, iDelta_x3, iDelta_g3;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi3N(ds,
				Xn[NODE1], Xn[NODE2], Xn[NODE3]);

			doublereal dN1 = ShapeFunc3N(ds, 1);
			doublereal dN2 = ShapeFunc3N(ds, 2);
			doublereal dN3 = ShapeFunc3N(ds, 3);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2 + X3Tmp*dN3);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2 + V3Tmp*dN3);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2 + Wn3*dN3);
			Vec3 gr(g1*dN1 + g3*dN3);
			Mat3x3 Gamma(ER_Rot::MatG, gr);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != 0) {
				Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR2*Mat3x3(ER_Rot::MatR, gr));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);
			/* Funzione di calcolo delle forze aerodinamiche */
			doublereal Fa0[6];
			Mat6x6 JFa;

			doublereal cc = dXds*dsdCsi*PW.dGetWght();

			Vec3 d(Xr - Xn[iNode]);

			Mat3x3 Theta1(RR2*Gamma*GammaInv1.MulMT(RR2*dN1));
			Mat3x3 Theta3(RR2*Gamma*GammaInv3.MulMT(RR2*dN3));
			Mat3x3 Theta2(Eye3 - Theta1 - Theta3);

			Vec3 Vrc(Vr*dCoef);
			Mat3x3 Bv1(Vrc.Cross(Theta1) - Mat3x3(Omega1Crossf1*(dN1*dCoef)));
			Mat3x3 Bv2(Vrc.Cross(Theta2) - Mat3x3(Omega2Crossf2*(dN2*dCoef)));
			Mat3x3 Bv3(Vrc.Cross(Theta3) - Mat3x3(Omega3Crossf3*(dN3*dCoef)));

			Vec3 Wrc(Wr*dCoef);
			Mat3x3 Bw1(Wrc.Cross(Theta1) - Mat3x3(Wn1*(dN1*dCoef)));
			Mat3x3 Bw2(Wrc.Cross(Theta2) - Mat3x3(Wn2*(dN2*dCoef)));
			Mat3x3 Bw3(Wrc.Cross(Theta3) - Mat3x3(Wn3*(dN3*dCoef)));

			if (iNumDof) {
				// prepare (v/dot{x} + dCoef*v/x) and so
				Mat3x3 RRlocT(RRloc.Transpose());
	
				vx.PutMat3x3(1, RRlocT*dN1);
				vx.PutMat3x3(4, RRloc.MulTM(Bv1 - Mat3x3(f1Tmp*dN1)));

				vx.PutMat3x3(6 + 1, RRlocT*dN2);
				vx.PutMat3x3(6 + 4, RRloc.MulTM(Bv2 - Mat3x3(f2Tmp*dN2)));

				vx.PutMat3x3(12 + 1, RRlocT*dN3);
				vx.PutMat3x3(12 + 4, RRloc.MulTM(Bv3 - Mat3x3(f3Tmp*dN3)));

				wx.PutMat3x3(4, RRlocT + Bw1);
				wx.PutMat3x3(6 + 4, RRlocT + Bw2);
				wx.PutMat3x3(12 + 4, RRlocT + Bw3);
	
				// equations from iFirstEq on are dealt with by aerodata
				aerodata->AssJac(WM, dCoef, XCurr, XPrimeCurr,
		         		iFirstEq, iFirstSubEq,
					vx, wx, fq, cq, iPnt, dW, Fa0, JFa, OUTA[iPnt]);
	
				// deal with (f/dot{q} + dCoef*f/q) and so
				integer iOffset = 18 + iPnt*iNumDof;
				for (integer iCol = 1; iCol <= iNumDof; iCol++) {
					Vec3 fqTmp((RRloc*fq.GetVec(iCol))*cc);
					Vec3 cqTmp(d.Cross(fqTmp) + (RRloc*cq.GetVec(iCol))*cc);

					WM.Sub(6*iNode + 1, iOffset + iCol, fqTmp);
					WM.Sub(6*iNode + 4, iOffset + iCol, cqTmp);
				}

				// first equation
				iFirstEq += iNumDof;
				iFirstSubEq += iNumDof;

			} else {
				aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, OUTA[iPnt]);
			}

			// rotate force, couple and Jacobian matrix in absolute frame
			Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

			// force and moment about the node
			Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
			Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + d.Cross(fTmp));

			Mat3x3 WM_F2[6];

			// f <-> x
			WM_F2[DELTAx1] = JFaR.GetMat11()*dN1;

			WM_F2[DELTAx2] = JFaR.GetMat11()*dN2;

			WM_F2[DELTAx3] = JFaR.GetMat11()*dN3;

			doublereal delta;

			// c <-> x
			delta = (iNode == NODE1) ? 1. : 0.;
			WM_M[DELTAx1] += JFaR.GetMat21()*dN1 - Mat3x3(fTmp*(dN1 - delta));

			delta = (iNode == NODE2) ? 1. : 0.;
			WM_M[DELTAx2] += JFaR.GetMat21()*dN2 - Mat3x3(fTmp*(dN2 - delta));

			delta = (iNode == NODE3) ? 1. : 0.;
			WM_M[DELTAx3] += JFaR.GetMat21()*dN3 - Mat3x3(fTmp*(dN3 - delta));

			// f <-> g
			WM_F2[DELTAg1] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f1Tmp))*dN1;
			WM_F2[DELTAg1] += JFaR.GetMat11()*Bv1 + JFaR.GetMat12()*Bw1;
			WM_F2[DELTAg1] -= fTmp.Cross(Theta1);
		
			WM_F2[DELTAg2] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f2Tmp))*dN2;
			WM_F2[DELTAg2] += JFaR.GetMat11()*Bv2 + JFaR.GetMat12()*Bw2;
			WM_F2[DELTAg2] -= fTmp.Cross(Theta2);
		
			WM_F2[DELTAg3] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f3Tmp))*dN3;
			WM_F2[DELTAg3] += JFaR.GetMat11()*Bv3 + JFaR.GetMat12()*Bw3;
			WM_F2[DELTAg3] -= fTmp.Cross(Theta3);

			// c <-> g
			WM_M[DELTAg1] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f1Tmp))*dN1;
			WM_M[DELTAg1] += JFaR.GetMat21()*Bv1 + JFaR.GetMat22()*Bw1;
			WM_M[DELTAg1] -= cTmp.Cross(Theta1);
			WM_M[DELTAg1] += Mat3x3(fTmp, f1Tmp*dN1);

			WM_M[DELTAg2] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f2Tmp))*dN2;
			WM_M[DELTAg2] += JFaR.GetMat21()*Bv2 + JFaR.GetMat22()*Bw2;
			WM_M[DELTAg2] -= cTmp.Cross(Theta2);
			WM_M[DELTAg2] += Mat3x3(fTmp, f2Tmp*dN2);

			WM_M[DELTAg3] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f3Tmp))*dN3;
			WM_M[DELTAg3] += JFaR.GetMat21()*Bv3 + JFaR.GetMat22()*Bw3;
			WM_M[DELTAg3] -= cTmp.Cross(Theta3);
			WM_M[DELTAg3] += Mat3x3(fTmp, f3Tmp*dN3);

			for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
				WM_F[iCnt] += WM_F2[iCnt];
				WM_M[iCnt] += d.Cross(WM_F2[iCnt]);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
			WM.Sub(6*iNode + 1, 3*iCnt + 1, WM_F[iCnt]);
			WM.Sub(6*iNode + 4, 3*iCnt + 1, WM_M[iCnt]);
		}
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
AerodynamicBeam::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssRes");

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
	}

	AssVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;
}

/* assemblaggio iniziale residuo */
SubVectorHandler&
AerodynamicBeam::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam::InitialAssRes");
	WorkVec.ResizeReset(18);

	integer iNode1FirstIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstPositionIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
	}

	AssVec(WorkVec, 1., XCurr, XCurr);

	return WorkVec;
}


/* assemblaggio residuo */
void
AerodynamicBeam::AssVec(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssVec");

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[3];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Xn[NODE3] = pNode3->GetXCurr();
	const Mat3x3& Rn3(pNode3->GetRCurr());
	const Vec3& Vn3(pNode3->GetVCurr());
	const Vec3& Wn3(pNode3->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);
	Vec3 f3Tmp(Rn3*f3);

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);
	Vec3 X3Tmp(Xn[NODE3] + f3Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));
	Vec3 V3Tmp(Vn3 + Wn3.Cross(f3Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);
	Mat3x3 RR3(Rn3*Ra3);

	/*
	 * Parametri di rotazione dai nodi 1 e 3 al nodo 2 (nell'ipotesi
	 * che tale trasformazione non dia luogo ad una singolarita')
	 */
	Vec3 g1(ER_Rot::Param, RR2.MulTM(RR1));
	Vec3 g3(ER_Rot::Param, RR2.MulTM(RR3));

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn[NODE2], rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	int iPnt = 0;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 18;

		integer iOffset = iFirstEq - 18;
		integer iNumRows = WorkVec.iGetSize();
		for (int iCnt = 18 + 1; iCnt <= iNumRows; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iOffset + iCnt);
		}
	}

	for (int iNode = 0; iNode < LASTNODE; iNode++) {

		/* Resetta le forze */
		F[iNode].Reset();
		M[iNode].Reset();

		doublereal dsi = pdsi3[iNode];
		doublereal dsf = pdsf3[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi3N(ds,
				Xn[NODE1], Xn[NODE2], Xn[NODE3]);

			doublereal dN1 = ShapeFunc3N(ds, 1);
			doublereal dN2 = ShapeFunc3N(ds, 2);
			doublereal dN3 = ShapeFunc3N(ds, 3);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2 + X3Tmp*dN3);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2 + V3Tmp*dN3);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2 + Wn3*dN3);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != 0) {
				Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR2*Mat3x3(ER_Rot::MatR, g1*dN1 + g3*dN3));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			doublereal dW[6];
			doublereal dTng[6];

			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);

			/* Funzione di calcolo delle forze aerodinamiche */
			if (iNumDof) {
				aerodata->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr,
                			iFirstEq, iFirstSubEq, iPnt, dW, dTng, OUTA[iPnt]);

				// first equation
				iFirstEq += iNumDof;
				iFirstSubEq += iNumDof;

			} else {
				aerodata->GetForces(iPnt, dW, dTng, OUTA[iPnt]);
			}

			/* Dimensionalizza le forze */
			doublereal dWght = dXds*dsdCsi*PW.dGetWght();
			dTng[1] *= TipLoss.dGet(dCsi);
			Vec3 FTmp(RRloc*(Vec3(&dTng[0])));
			Vec3 MTmp(RRloc*(Vec3(&dTng[3])));

			// Se e' definito il rotore, aggiungere il contributo alla trazione
			AddSectionalForce_int(iPnt, FTmp, MTmp, dWght, Xr, RRloc, Vr, Wr);

			FTmp *= dWght;
			MTmp *= dWght;
			F[iNode] += FTmp;
			M[iNode] += MTmp;
			M[iNode] += (Xr - Xn[iNode]).Cross(FTmp);

			// specific for Gauss points force output
			if (fToBeOutput()) {
				SetData(VTmp, dTng, Xr, RRloc, Vr, Wr, FTmp, MTmp);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		// Se e' definito il rotore, aggiungere il contributo alla trazione
		AddForce_int(F[iNode], M[iNode], Xn[iNode]);

		/* Somma il termine al residuo */
		WorkVec.Add(6*iNode + 1, F[iNode]);
		WorkVec.Add(6*iNode + 4, M[iNode]);
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBeam::Output(OutputHandler& OH) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Output");

	if (fToBeOutput()) {
		Aerodynamic2DElem<3>::Output_int(OH);

		if (OH.UseText(OutputHandler::AERODYNAMIC)) {
			std::ostream& out = OH.Aerodynamic() << std::setw(8) << GetLabel();

			switch (GetOutput()) {
			case AEROD_OUT_NODE:
				out << " " << std::setw(8) << pBeam->GetLabel()
					<< " ", F[NODE1].Write(out, " ") << " ", M[NODE1].Write(out, " ")
					<< " ", F[NODE2].Write(out, " ") << " ", M[NODE2].Write(out, " ")
					<< " ", F[NODE3].Write(out, " ") << " ", M[NODE3].Write(out, " ");
				break;
	
			case AEROD_OUT_PGAUSS:
				ASSERT(!pOutput.empty());
	
				for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
					i != pOutput.end(); i++)
				{
					out << " " << i->alpha
						<< " " << i->f;
				}
				break;
	
			case AEROD_OUT_STD:
				for (int i = 0; i < 3*GDI.iGetNum(); i++) {
					out
						<< " " << OUTA[i].alpha
						<< " " << OUTA[i].gamma
						<< " " << OUTA[i].mach
						<< " " << OUTA[i].cl
						<< " " << OUTA[i].cd
						<< " " << OUTA[i].cm
						<< " " << OUTA[i].alf1
						<< " " << OUTA[i].alf2;
				}
				break;
	
			default:
				ASSERT(0);
				break;
			}
	
			out << std::endl;
		}
	}
}

/* AerodynamicBeam - end */


/* Legge un elemento aerodinamico di trave */

Elem *
ReadAerodynamicBeam(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBeam");

	/* Trave */
	unsigned int uBeam = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << std::endl);

	/* verifica di esistenza della trave */
	Elem* p = pDM->pFindElem(Elem::BEAM, uBeam);
	if (p == 0) {
		silent_cerr("Beam3(" << uBeam << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	Beam *pBeam = dynamic_cast<Beam *>(p);
	if (pBeam == 0) {
		silent_cerr("Beam(" << uBeam << ") is not a Beam3 "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	ASSERT(pBeam != 0);

	/* Eventuale rotore */
	InducedVelocity* pIndVel = 0;
	bool bPassive(false);
	(void)ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBeam3",
		pIndVel, bPassive);

	/* Nodo 1: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode1 = pBeam->pGetNode(1);

	ReferenceFrame RF(pNode1);
	Vec3 f1(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

	Mat3x3 Ra1(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 1 rotation matrix: " << std::endl << Ra1 << std::endl);

	/* Nodo 2: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode2 = pBeam->pGetNode(2);

	RF = ReferenceFrame(pNode2);
	Vec3 f2(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);

	Mat3x3 Ra2(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 2 rotation matrix: " << std::endl << Ra2 << std::endl);

	/* Nodo 3: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode3 = pBeam->pGetNode(3);

	RF = ReferenceFrame(pNode3);
	Vec3 f3(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 3 offset: " << f3 << std::endl);

	Mat3x3 Ra3(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		"Node 3 rotation matrix: " << std::endl << Ra3 << std::endl);

	Shape* pChord = 0;
	Shape* pForce = 0;
	Shape* pVelocity = 0;
	Shape* pTwist = 0;
	Shape* pTipLoss = 0;

	integer iNumber = 0;
	DriveCaller* pDC = 0;
	AeroData* aerodata = 0;

	ReadAeroData(pDM, HP, 3,
		&pChord, &pForce, &pVelocity, &pTwist, &pTipLoss,
		&iNumber, &pDC, &aerodata);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = HP.GetYesNoOrBool(bDefaultUseJacobian);
	}

	if (aerodata->iGetNumDof() > 0 && !bUseJacobian) {
		silent_cerr("AerodynamicBeam3(" << uLabel << "): "
			"aerodynamic model needs \"jacobian, yes\" at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	unsigned uFlags = AerodynamicOutput::OUTPUT_NONE;
	ReadOptionalAerodynamicCustomOutput(pDM, HP, uLabel, uFlags, od);

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBeam3(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = 0;

	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBeam,
		AerodynamicBeam(uLabel, pDO, pBeam, pIndVel, bPassive,
			f1, f2, f3, Ra1, Ra2, Ra3,
			pChord, pForce, pVelocity, pTwist, pTipLoss,
			iNumber, aerodata, pDC, bUseJacobian, uFlags, od, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ostream& out = pDM->GetLogFile();
	out << "aero3: " << uLabel;

	Vec3 ra1 = Ra1.GetVec(1);
	Vec3 ra3 = Ra1.GetVec(3);
	doublereal dC = pChord->dGet(-1.);
	doublereal dP = pForce->dGet(-1.);
	out
		<< " " << pNode1->GetLabel()
		<< " ", (f1 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f1 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra2.GetVec(1);
	ra3 = Ra2.GetVec(3);
	dC = pChord->dGet(0.);
	dP = pForce->dGet(0.);
	out
		<< " " << pNode2->GetLabel()
		<< " ", (f2 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f2 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra3.GetVec(1);
	ra3 = Ra3.GetVec(3);
	dC = pChord->dGet(1.);
	dP = pForce->dGet(1.);
	out
		<< " " << pNode3->GetLabel()
		<< " ", (f3 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f3 + ra1*(dP + dC/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of ReadAerodynamicBeam() */


/* AerodynamicBeam2 - begin */

AerodynamicBeam2::AerodynamicBeam2(
	unsigned int uLabel,
	const DofOwner* pDO,
	const Beam2* pB,
	InducedVelocity* pR, bool bPassive,
	const Vec3& fTmp1,
	const Vec3& fTmp2,
	const Mat3x3& Ra1Tmp,
	const Mat3x3& Ra2Tmp,
	const Shape* pC,
	const Shape* pF,
	const Shape* pV,
	const Shape* pT,
	const Shape* pTL,
	integer iN,
	AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	unsigned uFlags, OrientationDescription ood,
	flag fOut
)
: Elem(uLabel, fOut),
Aerodynamic2DElem<2>(uLabel, pDO, pR, bPassive,
	pC, pF, pV, pT, pTL, iN, a, pDC, bUseJacobian, uFlags, ood, fOut),
pBeam(pB),
f1(fTmp1),
f2(fTmp2),
Ra1(Ra1Tmp),
Ra2(Ra2Tmp),
Ra1_3(Ra1Tmp.GetVec(3)),
Ra2_3(Ra2Tmp.GetVec(3))
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AerodynamicBeam2");

	ASSERT(pBeam != 0);
	ASSERT(pBeam->GetElemType() == Elem::BEAM);

	pNode1 = pBeam->pGetNode(1);
	pNode2 = pBeam->pGetNode(2);

	ASSERT(pNode1 != 0);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != 0);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

AerodynamicBeam2::~AerodynamicBeam2(void)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::~AerodynamicBeam2");
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBeam2::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBeam2::Restart");
	out << "  aerodynamic beam2: " << GetLabel()
		<< ", " << pBeam->GetLabel();
	if (pIndVel != 0) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f1.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra1.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra1.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f2.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra2.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra2.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

static const doublereal pdsi2[] = { -1., 0. };
static const doublereal pdsf2[] = { 0., 1. };

/* Jacobian assembly */
VariableSubMatrixHandler&
AerodynamicBeam2::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering AerodynamicBeam2::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	doublereal dW[6];

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[2];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);

	Vec3 fTmp[2];
	fTmp[NODE1] = f1Tmp;
	fTmp[NODE2] = f2Tmp;

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));

	Vec3 Omega1Crossf1(Wn1.Cross(f1Tmp));
	Vec3 Omega2Crossf2(Wn2.Cross(f2Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);

	/*
	 * half of relative rotation vector from node 1 to node 2
	 */
	Vec3 overline_theta(Vec3(ER_Rot::Param, RR1.MulTM(RR2))/2.);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	Vec3 Xmid = (Xn[NODE2] + Xn[NODE1])/2.;
	doublereal rho, c, p, T;
	GetAirProps(Xmid, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	int iPnt = 0;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 12;

		integer iOffset = iFirstEq - 12;

		for (int iCnt = 12 + 1; iCnt <= iNumRows; iCnt++) {
			WM.PutRowIndex(iCnt, iOffset + iCnt);
			WM.PutColIndex(iCnt, iOffset + iCnt);
		}
	}

	for (int iNode = 0; iNode < LASTNODE; iNode++) {
		doublereal dsi = pdsi2[iNode];
		doublereal dsf = pdsf2[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		Mat3x3 WM_F[4], WM_M[4];
		WM_F[DELTAx1].Reset();
		WM_F[DELTAg1].Reset();
		WM_F[DELTAx2].Reset();
		WM_F[DELTAg2].Reset();
		WM_M[DELTAx1].Reset();
		WM_M[DELTAg1].Reset();
		WM_M[DELTAx2].Reset();
		WM_M[DELTAg2].Reset();

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi2N(ds,
				Xn[NODE1], Xn[NODE2]);

			doublereal dN1 = ShapeFunc2N(ds, 1);
			doublereal dN2 = ShapeFunc2N(ds, 2);

			// note: identical to dN1 and dN2; see tecman.pdf
#if 0
			doublereal dNN1 = (1. + dN1 - dN2)/2.;
			doublereal dNN2 = (1. + dN2 - dN1)/2.;
#endif

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2);
			Vec3 thetar(overline_theta*dN2);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != 0) {
				Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR1*Mat3x3(ER_Rot::MatR, thetar));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);
			/* Funzione di calcolo delle forze aerodinamiche */
			doublereal Fa0[6];
			Mat6x6 JFa;

			doublereal cc = dXds*dsdCsi*PW.dGetWght();

			Vec3 d(Xr - Xn[iNode]);

			Mat3x3 Bv1((Vr - Omega1Crossf1)*(dN1*dCoef));
			Mat3x3 Bv2((Vr - Omega2Crossf2)*(dN2*dCoef));

			Mat3x3 Bw1((Wr - Wn1)*(dN1*dCoef));
			Mat3x3 Bw2((Wr - Wn2)*(dN2*dCoef));

			if (iNumDof) {
				// prepare (v/dot{x} + dCoef*v/x) and so
				Mat3x3 RRlocT(RRloc.Transpose());
	
				vx.PutMat3x3(1, RRlocT*dN1);
				vx.PutMat3x3(4, RRloc.MulTM(Bv1 - Mat3x3(f1Tmp*dN1)));

				vx.PutMat3x3(6 + 1, RRlocT*dN2);
				vx.PutMat3x3(6 + 4, RRloc.MulTM(Bv2 - Mat3x3(f2Tmp*dN2)));

				wx.PutMat3x3(4, RRlocT + Bw1);
				wx.PutMat3x3(6 + 4, RRlocT + Bw2);
	
				// equations from iFirstEq on are dealt with by aerodata
				aerodata->AssJac(WM, dCoef, XCurr, XPrimeCurr,
		         		iFirstEq, iFirstSubEq,
					vx, wx, fq, cq, iPnt, dW, Fa0, JFa, OUTA[iPnt]);
	
				// deal with (f/dot{q} + dCoef*f/q) and so
				integer iOffset = 12 + iPnt*iNumDof;
				for (integer iCol = 1; iCol <= iNumDof; iCol++) {
					Vec3 fqTmp((RRloc*fq.GetVec(iCol))*cc);
					Vec3 cqTmp(d.Cross(fqTmp) + (RRloc*cq.GetVec(iCol))*cc);

					WM.Sub(6*iNode + 1, iOffset + iCol, fqTmp);
					WM.Sub(6*iNode + 4, iOffset + iCol, cqTmp);
				}

				// first equation
				iFirstEq += iNumDof;
				iFirstSubEq += iNumDof;

			} else {
				aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, OUTA[iPnt]);
			}

			// rotate force, couple and Jacobian matrix in absolute frame
			Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

			// force and moment about the node
			Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
			Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + d.Cross(fTmp));

			Mat3x3 WM_F2[4];

			// f <-> x
			WM_F2[DELTAx1] = JFaR.GetMat11()*dN1;

			WM_F2[DELTAx2] = JFaR.GetMat11()*dN2;

			doublereal delta;

			// c <-> x
			delta = (iNode == NODE1) ? 1. : 0.;
			WM_M[DELTAx1] += JFaR.GetMat21()*dN1 - Mat3x3(fTmp*(dN1 - delta));

			delta = (iNode == NODE2) ? 1. : 0.;
			WM_M[DELTAx2] += JFaR.GetMat21()*dN2 - Mat3x3(fTmp*(dN2 - delta));

			// f <-> g
			WM_F2[DELTAg1] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f1Tmp))*dN1;
			WM_F2[DELTAg1] += JFaR.GetMat11()*Bv1 + JFaR.GetMat12()*Bw1;
			WM_F2[DELTAg1] -= Mat3x3(fTmp*dN1);
		
			WM_F2[DELTAg2] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f2Tmp))*dN2;
			WM_F2[DELTAg2] += JFaR.GetMat11()*Bv2 + JFaR.GetMat12()*Bw2;
			WM_F2[DELTAg2] -= Mat3x3(fTmp*dN2);
		
			// c <-> g
			WM_M[DELTAg1] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f1Tmp))*dN1;
			WM_M[DELTAg1] += JFaR.GetMat21()*Bv1 + JFaR.GetMat22()*Bw1;
			WM_M[DELTAg1] -= Mat3x3(cTmp*dN1);
			WM_M[DELTAg1] += Mat3x3(fTmp, f1Tmp*dN1);

			WM_M[DELTAg2] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f2Tmp))*dN2;
			WM_M[DELTAg2] += JFaR.GetMat21()*Bv2 + JFaR.GetMat22()*Bw2;
			WM_M[DELTAg2] -= Mat3x3(cTmp*dN2);
			WM_M[DELTAg2] += Mat3x3(fTmp, f2Tmp*dN2);

			for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
				WM_F[iCnt] += WM_F2[iCnt];
				WM_M[iCnt] += d.Cross(WM_F2[iCnt]);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
			WM.Sub(6*iNode + 1, 3*iCnt + 1, WM_F[iCnt]);
			WM.Sub(6*iNode + 4, 3*iCnt + 1, WM_M[iCnt]);
		}
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
AerodynamicBeam2::AssRes(
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AssRes");

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
	}

	AssVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;
}

/* assemblaggio iniziale residuo */
SubVectorHandler&
AerodynamicBeam2::InitialAssRes( SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::InitialAssRes");
	WorkVec.ResizeReset(12);

	integer iNode1FirstIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
	}

	AssVec(WorkVec, 1., XCurr, XCurr);

	return WorkVec;
}

/* assemblaggio residuo */
void
AerodynamicBeam2::AssVec(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AssVec");

	doublereal dTng[6];
	doublereal dW[6];

	Vec3 Xn[LASTNODE];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);

	/*
	 * half of relative rotation vector from node 1 to node 2
	 */
	Vec3 overline_theta(Vec3(ER_Rot::Param, RR1.MulTM(RR2))/2.);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != 0) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione di mezzo perche'
	 * non dovrebbero cambiare "molto")
	 */
	Vec3 Xmid = (Xn[NODE2] + Xn[NODE1])/2.;
	doublereal rho, c, p, T;
	GetAirProps(Xmid, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	int iPnt = 0;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 12;

		integer iOffset = iFirstEq - 12;
		integer iNumRows = WorkVec.iGetSize();
		for (int iCnt = 12 + 1; iCnt <= iNumRows; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iOffset + iCnt);
		}
	}

	for (int iNode = 0; iNode < LASTNODE; iNode++) {

		/* Resetta i dati */
		F[iNode].Reset();
		M[iNode].Reset();

		doublereal dsi = pdsi2[iNode];
		doublereal dsf = pdsf2[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi2N(ds, Xn[NODE1], Xn[NODE2]);

			doublereal dN1 = ShapeFunc2N(ds, 1);
			doublereal dN2 = ShapeFunc2N(ds, 2);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2);
			Vec3 thetar(overline_theta*((1. + dN2 - dN1)/2.));

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != 0) {
				Vr += pIndVel->GetInducedVelocity(GetElemType(),
				GetLabel(), iPnt, Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			dTw += dGet(); /* Contributo dell'eventuale sup. mobile */
			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
			 * lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR1*Mat3x3(ER_Rot::MatR, thetar));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw( dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(dW);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);

			/* Funzione di calcolo delle forze aerodinamiche */
			if (iNumDof) {
				aerodata->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr,
                			iFirstEq, iFirstSubEq, iPnt, dW, dTng, OUTA[iPnt]);

				// first equation
				iFirstEq += iNumDof;
				iFirstSubEq += iNumDof;

			} else {
				aerodata->GetForces(iPnt, dW, dTng, OUTA[iPnt]);
			}

			/* Dimensionalizza le forze */
			doublereal dWght = dXds*dsdCsi*PW.dGetWght();
			dTng[1] *= TipLoss.dGet(dCsi);
			Vec3 FTmp(RRloc*(Vec3(&dTng[0])));
			Vec3 MTmp(RRloc*(Vec3(&dTng[3])));

			// Se e' definito il rotore, aggiungere il contributo alla trazione
			AddSectionalForce_int(iPnt, FTmp, MTmp, dWght, Xr, RRloc, Vr, Wr);

			FTmp *= dWght;
			MTmp *= dWght;
			F[iNode] += FTmp;
			M[iNode] += MTmp;
			M[iNode] += (Xr - Xn[iNode]).Cross(FTmp);

			// specific for Gauss points force output
			if (fToBeOutput()) {
				SetData(VTmp, dTng, Xr, RRloc, Vr, Wr, FTmp, MTmp);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		// Se e' definito il rotore, aggiungere il contributo alla trazione
		AddForce_int(F[iNode], M[iNode], Xn[iNode]);

		/* Somma il termine al residuo */
		WorkVec.Add(6*iNode + 1, F[iNode]);
		WorkVec.Add(6*iNode + 4, M[iNode]);
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBeam2::Output(OutputHandler& OH ) const
{
	DEBUGCOUTFNAME("AerodynamicBeam2::Output");

	if (fToBeOutput()) {
		Aerodynamic2DElem<2>::Output_int(OH);

		if (OH.UseText(OutputHandler::AERODYNAMIC)) {
			std::ostream& out = OH.Aerodynamic() << std::setw(8) << GetLabel();

			switch (GetOutput()) {
	
			case AEROD_OUT_NODE:
				out << " " << std::setw(8) << pBeam->GetLabel()
					<< " ", F[NODE1].Write(out, " ") << " ", M[NODE1].Write(out, " ")
					<< " ", F[NODE2].Write(out, " ") << " ", M[NODE2].Write(out, " ");
				break;
	
			case AEROD_OUT_PGAUSS:
				ASSERT(!pOutput.empty());
	
				for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
					i != pOutput.end(); i++)
				{
					out << " " << i->alpha
						<< " " << i->f;
				}
				break;
	
			case AEROD_OUT_STD:
				for (int i = 0; i < 2*GDI.iGetNum(); i++) {
		 			out
						<< " " << OUTA[i].alpha
						<< " " << OUTA[i].gamma
						<< " " << OUTA[i].mach
						<< " " << OUTA[i].cl
						<< " " << OUTA[i].cd
						<< " " << OUTA[i].cm
						<< " " << OUTA[i].alf1
						<< " " << OUTA[i].alf2;
				}
				break;
	
			default:
				ASSERT(0);
				break;
			}
	
			out << std::endl;
		}
	}
}
	
/* AerodynamicBeam2 - end */


/* Legge un elemento aerodinamico di trave a due nodi */

Elem *
ReadAerodynamicBeam2(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBeam2");

	/* Trave */
	unsigned int uBeam = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << std::endl);

	/* verifica di esistenza della trave */
	Elem *p = pDM->pFindElem(Elem::BEAM, uBeam);
	if (p == 0) {
		silent_cerr("Beam2(" << uBeam << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	Beam2* pBeam = dynamic_cast<Beam2 *>(p);
	if (pBeam == 0) {
		silent_cerr("Beam(" << uBeam << ") is not a Beam2 "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Eventuale rotore */
	InducedVelocity* pIndVel = 0;
	bool bPassive(false);
	(void)ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBeam2",
		pIndVel, bPassive);

	/* Nodo 1: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode1 = pBeam->pGetNode(1);

	ReferenceFrame RF(pNode1);
	Vec3 f1(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

	Mat3x3 Ra1(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 1 rotation matrix: " << std::endl << Ra1 << std::endl);

	/* Nodo 2: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode2 = pBeam->pGetNode(2);

	RF = ReferenceFrame(pNode2);
	Vec3 f2(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);

	Mat3x3 Ra2(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		"Node 2 rotation matrix: " << std::endl << Ra2 << std::endl);

	Shape* pChord = 0;
	Shape* pForce = 0;
	Shape* pVelocity = 0;
	Shape* pTwist = 0;
	Shape* pTipLoss = 0;

	integer iNumber = 0;
	DriveCaller* pDC = 0;
	AeroData* aerodata = 0;

	ReadAeroData(pDM, HP, 2,
		&pChord, &pForce, &pVelocity, &pTwist, &pTipLoss,
		&iNumber, &pDC, &aerodata);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = HP.GetYesNoOrBool(bDefaultUseJacobian);
	}

	if (aerodata->iGetNumDof() > 0 && !bUseJacobian) {
		silent_cerr("AerodynamicBeam2(" << uLabel << "): "
			"aerodynamic model needs \"jacobian, yes\" at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	unsigned uFlags = AerodynamicOutput::OUTPUT_NONE;
	ReadOptionalAerodynamicCustomOutput(pDM, HP, uLabel, uFlags, od);

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBeam2(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = 0;

	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBeam2,
		AerodynamicBeam2(uLabel, pDO, pBeam, pIndVel, bPassive,
			f1, f2, Ra1, Ra2,
			pChord, pForce, pVelocity, pTwist, pTipLoss,
			iNumber, aerodata, pDC, bUseJacobian, uFlags, od, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ostream& out = pDM->GetLogFile();
	out << "aero2: " << uLabel;

	Vec3 ra1 = Ra1.GetVec(1);
	Vec3 ra3 = Ra1.GetVec(3);
	doublereal dC = pChord->dGet(-1.);
	doublereal dP = pForce->dGet(-1.);
	out
		<< " " << pNode1->GetLabel()
		<< " ", (f1 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f1 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra2.GetVec(1);
	ra3 = Ra2.GetVec(3);
	dC = pChord->dGet(1.);
	dP = pForce->dGet(1.);
	out
		<< " " << pNode2->GetLabel()
		<< " ", (f2 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f2 + ra1*(dP + dC/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of ReadAerodynamicBeam2() */

