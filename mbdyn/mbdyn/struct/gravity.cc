/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* Elemento accelerazione di gravita' */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "gravity.h"

/* Gravity - begin */

Gravity::Gravity(flag fOut)
: Elem(1, fOut)
{
	NO_OP;
}

Gravity::~Gravity(void)
{
	NO_OP;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
Gravity::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering Gravity::AssJac()" << std::endl);
	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
Gravity::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering Gravity::AssRes()" << std::endl);
	WorkVec.Resize(0);
	return WorkVec;
}

/* Gravity - end */


/* UniformGravity - begin */

UniformGravity::UniformGravity(const TplDriveCaller<Vec3>* pDC, flag fOut)
: Elem(1, fOut), Gravity(fOut), TplDriveOwner<Vec3>(pDC)
{
	Acc = Get();
}

UniformGravity::~UniformGravity(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
UniformGravity::Restart(std::ostream& out) const
{
	return out << "  gravity: uniform, /* reference, global, */ ",
		pGetDriveCaller()->Restart(out) << ";" << std::endl;
}

/* assemblaggio residuo */
SubVectorHandler&
UniformGravity::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering UniformGravity::AssRes()" << std::endl);

	/* Approfitto del fatto che Gravity viene aggiornato prima
	 * degli altri elementi (vedi l'enum Elem::Type e la sequenza di
	 * assemblaggio) per fargli calcolare Acc una volta per tutte.
	 * Quindi, quando viene chiamata GetAcceleration(void),
	 * questa restituisce un reference all'accelerazione con il
	 * minimo overhead
	 */
	Acc = Get();
	return Gravity::AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
}

void
UniformGravity::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::GRAVITY)) {
			OH.Gravity() << std::setw(8) << GetLabel()
				<< " " << Acc << std::endl;
		}
	}
}

/* UniformGravity - end */


/* CentralGravity - begin */

CentralGravity::CentralGravity(const Vec3& X0,
	doublereal dM, doublereal dG, flag fOut)
: Elem(1, fOut), Gravity(fOut), m_X0(X0), m_dM(dM), m_dG(dG)
{
	NO_OP;
}

CentralGravity::~CentralGravity(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
CentralGravity::Restart(std::ostream& out) const
{
	return out << "  gravity: central, "
		"origin, ", m_X0.Write(out, ", ") << ", "
		"mass, " << m_dM << ", "
		"G, " << m_dG << ";"
		<< std::endl;
}

void
CentralGravity::Output(OutputHandler& OH) const
{
#if 0
	// nothing to output...
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::GRAVITY)) {
			OH.Gravity() << std::setw(8) << GetLabel()
				<< std::endl;
		}
	}
#endif
	NO_OP;
}

Vec3
CentralGravity::GetAcceleration(const Vec3& X) const
{
	Vec3 D = m_X0 - X;
	doublereal dD = D.Norm();
	return D*(m_dM*m_dG/(dD*dD*dD));
}

/* CentralGravity - end */


/* GravityOwner - begin */

GravityOwner::GravityOwner(void)
: pGravity(0)
{
	NO_OP;
}


GravityOwner::~GravityOwner(void)
{
	NO_OP;
}


void
GravityOwner::PutGravity(const Gravity* pG)
{
	ASSERT(pGravity == 0);
	pGravity = const_cast<Gravity *>(pG);
}

bool
GravityOwner::bGetGravity(const Vec3& X, Vec3& Acc) const
{
	if (pGravity == 0) {
		return false;
	}

	Acc = pGravity->GetAcceleration(X);
	return true;
}

/* GravityOwner - end */


/* ElemGravityOwner - begin */

ElemGravityOwner::ElemGravityOwner(unsigned int uL, flag fOut)
: Elem(uL, fOut), GravityOwner()
{
	NO_OP;
}

ElemGravityOwner::~ElemGravityOwner(void)
{
	NO_OP;
}

/* ElemGravityOwner - end */

Elem *
ReadGravity(DataManager* pDM, MBDynParser& HP)
{
	Elem *pE = 0;

	if (HP.IsKeyWord("central")) {
		Vec3 X0 = ::Zero3;
		if (HP.IsKeyWord("origin")) {
			X0 = HP.GetPosAbs(::AbsRefFrame);
		}

		if (!HP.IsKeyWord("mass")) {
			silent_cerr("Gravity: \"mass\" keyword expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dM = HP.GetReal();
		if (dM <= std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Gravity: mass " << dM << " is negative or too small at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("G")) {
			silent_cerr("Gravity: \"G\" keyword expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		// 6.673 84 x 10-11 m3 kg-1 s-2, http://physics.nist.gov/cuu/Constants/index.html
		doublereal dG = 6.67384e-11;
		if (!HP.IsKeyWord("si")) {
			dG = HP.GetReal();
			if (dG <= std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("Gravity: gravity constant " << dG << " is negative or too small at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		SAFENEWWITHCONSTRUCTOR(pE, CentralGravity, CentralGravity(X0, dM, dG, 0));

	} else {
		if (!HP.IsKeyWord("uniform")) {
			silent_cerr("Gravity: "
				"<type> expected "
				"at line " << HP.GetLineData() << "; "
				"assuming \"uniform\""
				<< std::endl);
		}

		TplDriveCaller<Vec3>* pDC = ReadDC3D(pDM, HP);

		flag fOut = pDM->fReadOutput(HP, Elem::GRAVITY);

		SAFENEWWITHCONSTRUCTOR(pE, UniformGravity, UniformGravity(pDC, fOut));
	}

	return pE;
}

