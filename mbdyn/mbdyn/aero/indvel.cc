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

/* Rotore */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>

#ifdef USE_MPI
#include "mysleep.h"
const int mysleeptime = 300;

#ifdef MPI_PROFILING
extern "C" {
#include <mpe.h>
}
#include <cstdio>
#endif /* MPI_PROFILING */
#include "mbcomm.h"
#endif /* USE_MPI */

#include "indvel.h"
#include "dataman.h"

/* InducedVelocity - begin */

InducedVelocity::InducedVelocity(unsigned int uL,
	const StructNode *pCraft,
	ResForceSet **ppres, flag fOut)
: Elem(uL, fOut),
#ifdef USE_MPI
is_parallel(false),
pBlockLenght(0),
pDispl(0),
ReqV(MPI::REQUEST_NULL),
pIndVelDataType(0),
#endif /* USE_MPI */
pCraft(pCraft),
ppRes(ppres)
{
#ifdef USE_MPI
	iForcesVecDim = 6;
	for (int i = 0; ppRes && ppRes[i]; i++) {
		iForcesVecDim += 6;
	}
	SAFENEWARR(pTmpVecR, doublereal, iForcesVecDim);
	SAFENEWARR(pTmpVecS, doublereal, iForcesVecDim);
	if (MPI::Is_initialized()) {
		is_parallel = true;
   		IndVelComm = MBDynComm.Dup();
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_init(&induced_velocity_mutex, NULL);
	pthread_cond_init(&induced_velocity_cond, NULL);
	pthread_mutex_init(&forces_mutex, NULL);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}

InducedVelocity::~InducedVelocity(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pTmpVecR);
	SAFEDELETEARR(pTmpVecS);
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_destroy(&induced_velocity_mutex);
	pthread_cond_destroy(&induced_velocity_cond);
	pthread_mutex_destroy(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}

bool
InducedVelocity::bSectionalForces(void) const
{
	return false;
}

unsigned int
InducedVelocity::iGetNumPrivData(void) const {
	return 6;
}

unsigned int
InducedVelocity::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != 0);

	unsigned int idx = 0;

	// sanity check
	if (s[0] == '\0' || s[1] == '\0' || s[2] != '\0' ) {
		return 0;
	}

	switch (s[0]) {
	case 'M':
		idx += 3;
		// fallthru

	case 'T':
		switch (s[1]) {
		case 'x':
			return idx + 1;

		case 'y':
			return idx + 2;

		case 'z':
			return idx + 3;
		}
	}

	return 0;
}

doublereal
InducedVelocity::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= 6);

	switch (i) {
	case 1:
	case 2:
	case 3:
		return Res.Force()(i);

	case 4:
	case 5:
	case 6:
		return Res.Moment()(i - 3);
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
InducedVelocity::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	ASSERT(bDone);
	bDone = false;
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
}

/* assemblaggio jacobiano (nullo per tutti tranne che per il DynamicInflow) */
VariableSubMatrixHandler&
InducedVelocity::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering InducedVelocity::AssJac()" << std::endl);
	WorkMat.SetNullMatrix();

	return WorkMat;
}

#ifdef USE_MPI
void
InducedVelocity::ExchangeLoads(flag fWhat)
{
#ifdef MPI_PROFILING
	MPE_Log_event(33, 0, "start Induced Velocity Loads Exchange ");
#endif /* MPI_PROFILING */
	/* Se il rotore è connesso ad una sola macchina non è necessario scambiare messaggi */
	if (is_parallel && IndVelComm.Get_size() > 1) {
		if (fWhat) {
			/* Scambia F e M */
			Res.Force().PutTo(&pTmpVecS[0]);
			Res.Moment().PutTo(&pTmpVecS[3]);
			for (int i = 0; ppRes && ppRes[i]; i++) {
				ppRes[i]->pRes->Force().PutTo(&pTmpVecS[6 + 6*i]);
				ppRes[i]->pRes->Moment().PutTo(&pTmpVecS[9 + 6*i]);
			}
			IndVelComm.Allreduce(pTmpVecS, pTmpVecR, iForcesVecDim, MPI::DOUBLE, MPI::SUM);
			Res.PutForces(Vec3(&pTmpVecR[0]), Vec3(&pTmpVecR[3]));
			for (int i = 0; ppRes && ppRes[i]; i++) {
				ppRes[i]->pRes->PutForces(Vec3(&pTmpVecR[6 + 6*i]),  Vec3(&pTmpVecR[9 + 6*i]));
			}

		} else {
			IndVelComm.Allreduce(Res.Force().pGetVec(), pTmpVecR, 3, MPI::DOUBLE, MPI::SUM);
			Res.PutForce(Vec3(pTmpVecR));
		}
	}
#ifdef MPI_PROFILING
	MPE_Log_event(34, 0, "end Induced Velocity Loads Exchange ");
#endif /* MPI_PROFILING */
}

void
InducedVelocity::InitializeIndVelComm(MPI::Intracomm* iv)
{
	ASSERT(is_parallel);
  	IndVelComm = *iv;
}

void
InducedVelocity::ExchangeVelocity(void)
{
#define ROTDATATYPELABEL	100
	if (is_parallel && IndVelComm.Get_size() > 1) {
		if (IndVelComm.Get_rank() == 0) {
			for (int i = 1; i < IndVelComm.Get_size(); i++) {
				IndVelComm.Send(MPI::BOTTOM, 1, *pIndVelDataType,
						i, ROTDATATYPELABEL);
			}
		} else {
			ReqV = IndVelComm.Irecv((void *)MPI::BOTTOM, 1,
					*pIndVelDataType, 0, ROTDATATYPELABEL);
		}
	}
}
#endif // USE_MPI

// Somma alla trazione il contributo di forza di un elemento generico
void
InducedVelocity::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
	for (int i = 0; ppRes && ppRes[i]; i++) {
		if (ppRes[i]->is_in(pEl->GetLabel())) {
			ppRes[i]->pRes->AddForces(F, M, X);
		}
	}
}

// Somma alla trazione il contributo di forza di un elemento generico
// usando la forza e il momento per unita' di apertura e il peso
void
InducedVelocity::AddSectionalForce(Elem::Type type,
	const Elem *pEl, unsigned uPnt,
	const Vec3& F, const Vec3& M, doublereal dW,
	const Vec3& X, const Mat3x3& R,
	const Vec3& V, const Vec3& W)
{
	ASSERT(bSectionalForces() == true);
	InducedVelocity::AddForce(pEl, 0, F*dW, M*dW, X);
}

void
InducedVelocity::ResetForce(void)
{
	Res.Reset(pCraft->GetXCurr());
	for (int i = 0; ppRes && ppRes[i]; i++) {
		ppRes[i]->pRes->Reset();
	}
}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
void
InducedVelocity::Wait(void) const
{
	pthread_mutex_lock(&induced_velocity_mutex);
	if (!bDone) {
		pthread_cond_wait(&induced_velocity_cond,
				&induced_velocity_mutex);
	}
	pthread_mutex_unlock(&induced_velocity_mutex);
}

void
InducedVelocity::Done(void) const
{
	pthread_mutex_lock(&induced_velocity_mutex);
	ASSERT(!bDone);
	bDone = true;
	pthread_cond_broadcast(&induced_velocity_cond);
	pthread_mutex_unlock(&induced_velocity_mutex);
}
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

/* InducedVelocity - end */

/* InducedVelocityElem - begin */

InducedVelocityElem::InducedVelocityElem(unsigned int uL, const DofOwner* pDO,
	const StructNode *pCraft,
	ResForceSet **ppres, flag fOut)
: Elem(uL, fOut),
AerodynamicElem(uL, pDO, fOut),
InducedVelocity(uL, pCraft, ppres, fOut)
{
	NO_OP;
}

InducedVelocityElem::~InducedVelocityElem(void)
{
	NO_OP;
}

/* Tipo dell'elemento (usato per debug ecc.) */
Elem::Type
InducedVelocityElem::GetElemType(void) const
{
	return Elem::INDUCEDVELOCITY;
}

/* Tipo dell'elemento (usato per debug ecc.) */
AerodynamicElem::Type
InducedVelocityElem::GetAerodynamicElemType(void) const
{
	return AerodynamicElem::INDUCEDVELOCITY;
}

/* InducedVelocity - end */

