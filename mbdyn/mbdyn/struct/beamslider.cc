/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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
 * Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <constltp.h>
#include <shapefnc.h>
#include <beamslider.h>
#include <dataman.h>

/* Costruttore non banale */
BeamSliderJoint::BeamSliderJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN, enum Type iT,
		unsigned int nB, const Beam** ppB,
		const Vec3& fTmp, const Mat3x3& RTmp, flag fOut)
: Elem(uL, Elem::JOINT, fOut),
Joint(uL, Joint::BEAMSLIDER, pDO, fOut),
nRotConstr(0), nBeams(nB), iCurrBeam(0), iType(iT),
pNode(pN), ppBeam(ppB),
f(fTmp), R(RTmp),
F(0.), M(0.)
{
	ASSERT(pNode != NULL);
	ASSERT(nBeams > 0);
	ASSERT(ppBeam != NULL);

	switch (iType) {
	case CLASSIC:
		nRotConstr = 2;
		break;
	case SPLINE:
		nRotConstr = 3;
		break;
	default:
		break;
	}

	for (unsigned int iCnt = 0; iCnt < nBeams; iCnt++) {
		ASSERT(ppBeam[iCnt] != NULL);
	}
}

/* Assemblaggio residuo */
SubVectorHandler& 
BeamSliderJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering BeamSliderJoint::AssRes()" << endl);
	
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);

	/* Indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	const StructNode *pBeamNode[Beam::NUMNODES];
	
	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.fPutRowIndex(iCnt, iNodeFirstMomIndex+iCnt);
	}
	for (int nCnt = 1; nCnt <= Beam::NUMNODES; nCnt++) {
		pBeamNode[nCnt-1] = ppBeam[iCurrBeam]->pGetNode(nCnt);
		integer iBeamFirstMomIndex = 
			pBeamNode[nCnt-1]->iGetFirstMomentumIndex();

		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.fPutRowIndex(6*nCnt+iCnt, 
					iBeamFirstMomIndex+iCnt);
		}
	}
	
	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WorkVec.fPutRowIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
	}
	
	/* Aggiorna i dati propri */
	F = Vec3(XCurr, iFirstReactionIndex+1);
	M = Vec3(0.);
	for (unsigned int iCnt = 1; iCnt <= nRotConstr; iCnt++) {
		M.Put(iCnt, XCurr.dGetCoef(iFirstReactionIndex+3+iCnt));
	}
	
	return WorkVec;
}

