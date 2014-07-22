/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/*
 * Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "constltp.h"
#include "shapefnc.h"
#include "beamslider.h"


/* BeamConn - begin */
BeamConn::BeamConn(const Beam *pB, 
		const Vec3& f1, const Vec3& f2, const Vec3& f3,
		const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3)
: m_pBeam(pB)
{
	m_f[0] = f1;
	m_f[1] = f2;
	m_f[2] = f3;

	m_R[0] = R1;
	m_R[1] = R2;
	m_R[2] = R3;
}

BeamConn::~BeamConn(void)
{
	NO_OP;
}
	
/* BeamConn - end */


/* BeamSliderJoint - begin */

/* Punto di valutazione */
const doublereal dS = 1./sqrt(3.);

/* Costruttore non banale */
BeamSliderJoint::BeamSliderJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN, enum Type iT,
		unsigned int nB, const BeamConn *const * ppB,
		unsigned int uIB, unsigned int uIN,
		doublereal dl,
		const Vec3& fTmp, const Mat3x3& RTmp, flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
nRotConstr(0), nBeams(nB), iCurrBeam(0), iType(iT),
pNode(pN), ppBeam(ppB),
f(fTmp), R(RTmp),
F(Zero3), M(Zero3),
sRef(0.), s(0.),
dL(dl),
x(Zero3), l(Zero3)
{
	ASSERT(pNode != NULL);
	ASSERT(nBeams > 0);
	ASSERT(ppBeam != NULL);
	ASSERT(uIB > 0 && uIB <= nBeams);
	ASSERT(uIN > 0 && uIN <= 3);

	iCurrBeam = uIB-1;

	switch (uIN) {
	case 1:
		activeNode = Beam::NODE1;
		break;

	case 2:
		activeNode = Beam::NODE2;
		break;

	case 3:
		activeNode = Beam::NODE3;
		break;
	}

	sRef = 2.*iCurrBeam + activeNode - 2.;
	s = sRef;

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

#ifdef DEBUG
	for (unsigned int i = 0; i < nBeams; i++) {
		ASSERT(ppBeam[i] != NULL);
	}
#endif /* DEBUG */
}

BeamSliderJoint::~BeamSliderJoint(void)
{
	NO_OP;
}

std::ostream& 
BeamSliderJoint::Restart(std::ostream& out) const
{
	return out << "# beam slider not implemented yet" << std::endl
		<< "beam slider;" << std::endl;
}

void 
BeamSliderJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 RTmp(pNode->GetRCurr()*R);
		Mat3x3 RTmpT(RTmp.Transpose());
		
		Joint::Output(OH.Joints(), "BeamSlider", GetLabel(),
				RTmpT*F, M, F, RTmp*M)
			<< " " << ppBeam[iCurrBeam]->pGetBeam()->GetLabel()
			<< " " << sRef << " " << l << std::endl;
	}
}

#if 0
unsigned int
BeamSliderJoint::iGetNumPrivData(void) const
{
	return 0;
}

unsigned int
BeamSliderJoint::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

doublereal 
BeamSliderJoint::dGetPrivData(unsigned int i) const
{
	return 0.;
}
#endif

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
BeamSliderJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering BeamSliderJoint::AssJac()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(iNumRows, iNumCols);

	/* Indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	const StructNode *pBeamNode[Beam::NUMNODES];

	/*
	 *  1 =>  6:	nodo corpo
	 *  7 => 12:	nodo 1 trave
	 * 13 => 18:	nodo 2 trave
	 * 19 => 24:	nodo 3 trave
	 *       25:	l^T F = 0 (s)
	 * 26 => 28:	vincolo posizione (F)
	 * 29 => 31:	vincoli rotazione, se presenti
	 */
	
	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex+iCnt);
	}
	
	for (int nCnt = 1; nCnt <= Beam::NUMNODES; nCnt++) {
		pBeamNode[nCnt-1] = ppBeam[iCurrBeam]->pGetNode(nCnt);
		integer iBeamFirstMomIndex = 
			pBeamNode[nCnt-1]->iGetFirstMomentumIndex();
		integer iBeamFirstPosIndex = 
			pBeamNode[nCnt-1]->iGetFirstPositionIndex();

		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(6*nCnt+iCnt, 
					iBeamFirstMomIndex+iCnt);
			WM.PutColIndex(6*nCnt+iCnt, 
					iBeamFirstPosIndex+iCnt);
		}
	}
	
	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WM.PutRowIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
		WM.PutColIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
	}

	/* vincolo in posizione */
	for (unsigned int i = 1; i <= 3; i++) {
		/* l^T F = 0 : Delta F */
		doublereal d = l.dGet(i)/dCoef;
		WM.DecCoef(6*(1+Beam::NUMNODES)+1, 
				6*(1+Beam::NUMNODES)+1+i, d);

		/* xc - x = 0: l Delta s */
		WM.IncCoef(6*(1+Beam::NUMNODES)+1+i, 
				6*(1+Beam::NUMNODES)+1, d);

		/* xc - x = 0: Delta x_b */
		WM.DecCoef(6*(1+Beam::NUMNODES)+1+i, i, 1.);
	}

	Vec3 lp(Zero3);
	for (unsigned int iN = 0; iN < Beam::NUMNODES; iN++) {
		Vec3 Tmp(fTmp[iN].Cross(F));

		/* l^T F = 0 : Delta s */
		lp += xTmp[iN]*dNpp[iN];

		for (unsigned int i = 1; i <= 3; i++) {
			/* l^T F = 0 : Delta x */
			doublereal d = F.dGet(i)*dNp[iN];
			WM.DecCoef(6*(1+Beam::NUMNODES)+1, 
					6*(1+iN)+i, d);

			/* l^T F = 0 : Delta g */
			d = Tmp.dGet(i)*dNp[iN];
			WM.DecCoef(6*(1+Beam::NUMNODES)+1, 
					6*(1+iN)+3+i, d);

			/* xc - x = 0: Delta x */
			WM.IncCoef(6*(1+Beam::NUMNODES)+1+i, 
					6*(1+iN)+i, dN[iN]);
		}

		/* xc - x = 0: Delta g */
		WM.Sub(6*(1+Beam::NUMNODES)+1+1, 
				6*(1+iN)+3+1, Mat3x3(MatCross, fTmp[iN]*dN[iN]));
	}

	/* l^T F = 0 : Delta s */
	WM.DecCoef(6*(1+Beam::NUMNODES)+1, 
			6*(1+Beam::NUMNODES)+1, (F*lp)/dCoef);

	/* reazioni vincolari */
	for (unsigned int i = 1; i <= 3; i++) {
		/* corpo: Delta F */
		WM.DecCoef(i, 6*(1+Beam::NUMNODES)+1+i, 1.);

		/* trave: Delta F */
		WM.IncCoef(6*activeNode+i, 6*(1+Beam::NUMNODES)+1+i, dW[0]);
	}

	/* corpo: Delta F (momento) */
	Mat3x3 MTmp(MatCross, fb);
	WM.Sub(3+1, 6*(1+Beam::NUMNODES)+1+1, MTmp);

	/* vincolo posizione: Delta gb */
	WM.Add(6*(1+Beam::NUMNODES)+1+1, 3+1, MTmp);

	/* corpo: Delta gb (momento) */
	Mat3x3 Ffb(MatCrossCross, F, fb*dCoef);
	WM.Sub(3+1, 3+1, Ffb);

	/* trave: Delta gb (momento) */
	WM.Add(6*activeNode+3+1, 3+1, Ffb*dW[0]);

	/* trave: Delta F (momento) */
	MTmp = Mat3x3(MatCross, F*(dCoef*dW[0]));
	WM.Add(6*activeNode+3+1, 6*(1+Beam::NUMNODES)+1+1, 
			Mat3x3(MatCross, (xc - xNod[activeNode-1])*dW[0]));
	WM.Sub(6*activeNode+3+1, 1, MTmp);
	WM.Add(6*activeNode+3+1, 6*activeNode+1, MTmp);

	if (dW[1] != 0.) {
		/* trave: Delta gb (momento) */
		WM.Add(6*(activeNode+1)+3+1, 3+1, Ffb*dW[1]);

/* NOTE: these terms are questionable: the perturbation
 * related to the amplitude of the smearing should go 
 * in the jacobian, but there is no evidence of improvements
 * in convergence */
#define DELTADW
#ifdef DELTADW
		Vec3 m1(M+(xc-xNod[activeNode-1]).Cross(F));
		Vec3 m2(M+(xc-xNod[activeNode]).Cross(F));
#endif /* DELTADW */

		/* reazioni vincolari */
		for (unsigned int i = 1; i <= 3; i++) {
			/* trave: Delta F */
			WM.IncCoef(6*(activeNode+1)+i, 
					6*(1+Beam::NUMNODES)+1+i, dW[1]);

#ifdef DELTADW
			/* trave: Delta s (Delta dW forza) */
			doublereal d = F(i)/(2.*dL);
			WM.DecCoef(6*activeNode+i, 
					6*(1+Beam::NUMNODES)+1, d);
			WM.IncCoef(6*(activeNode+1)+i, 
					6*(1+Beam::NUMNODES)+1, d);

			/* trave: Delta s (Delta dW momento) */
			d = m1(i)/(2.*dL);
			WM.DecCoef(6*activeNode+3+i, 
					6*(1+Beam::NUMNODES)+1, d);
			d = m2(i)/(2.*dL);
			WM.IncCoef(6*(activeNode+1)+3+i, 
					6*(1+Beam::NUMNODES)+1, d);
#endif /* DELTADW */
		}

		/* trave: Delta F (momento) */
		Mat3x3 MTmp(MatCross, F*(dCoef*dW[1]));
		WM.Add(6*(activeNode+1)+3+1, 6*(1+Beam::NUMNODES)+1+1, 
				Mat3x3(MatCross, (xc-xNod[activeNode]))*dW[1]);
		WM.Sub(6*(activeNode+1)+3+1, 1, MTmp);
		WM.Add(6*(activeNode+1)+3+1, 6*(activeNode+1)+1, MTmp);
	}

	/* Vincolo in rotazione */
	if (iType != SPHERICAL) {
		Vec3 eb2 = Rb.GetVec(2);
		Vec3 eb3 = Rb.GetVec(3);

		Vec3 mm(eb2*m(2) + eb3*m(3));

		doublereal d;

		for (unsigned int iN = 0; iN < Beam::NUMNODES; iN++) {
			Vec3 Tmpf2(fTmp[iN].Cross(eb2));
			Vec3 Tmpf3(fTmp[iN].Cross(eb3));

			for (unsigned int i = 1; i <= 3; i++) {
				/* Vincolo in rotazione: Delta x */
				d = eb2.dGet(i)*dNp[iN];
				WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+1, 
						6*(1+iN)+i, d);

				d = eb3.dGet(i)*dNp[iN];
				WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+2, 
						6*(1+iN)+i, d);

				/* Vincolo in rotazione: Delta g */
				d = Tmpf2.dGet(i)*dNp[iN];
				WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+1,
						6*(1+iN)+3+i, d);

				d = Tmpf3.dGet(i)*dNp[iN];
				WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+2,
						6*(1+iN)+3+i, d);

			}

			Vec3 mmTmp(mm*(dNp[iN]*dCoef));
			Mat3x3 mmTmp2(MatCross, mmTmp);
			Mat3x3 mmTmp3(MatCrossCross, mmTmp, fTmp[iN]);

#if 0
			Vec3 MTmp(M*(dNp[iN]*dCoef));
			Mat3x3 MTmp2(MTmp);
			Mat3x3 MTmp3(MTmp, fTmp[iN]);
#endif
			
			/* Reazione vincolare: Delta x */
			WM.Sub(3+1, 6*(1+iN)+1, mmTmp2);

			/* Reazione vincolare: Delta g */
			WM.Add(3+1, 6*(1+iN)+3+1, mmTmp3);
		
			if (dW[1] == 0.) {
				/* Reazione vincolare: Delta x */
				WM.Add(6*activeNode+3+1, 6*(1+iN)+1, mmTmp2);

				/* Reazione vincolare: Delta g */
				WM.Sub(6*activeNode+3+1, 6*(1+iN)+3+1, mmTmp3);
			} else {
				/* Reazione vincolare: Delta x */
				WM.Add(6*activeNode+3+1, 6*(1+iN)+1, 
						mmTmp2*dW[0]);
				WM.Add(6*(activeNode+1)+3+1, 6*(1+iN)+1, 
						mmTmp2*dW[1]);

				/* Reazione vincolare: Delta g */
				WM.Sub(6*activeNode+3+1, 6*(1+iN)+3+1, 
						mmTmp3*dW[0]);
				WM.Sub(6*(activeNode+1)+3+1, 6*(1+iN)+3+1, 
						mmTmp3*dW[1]);
			}
		}

		Vec3 Tmpl2(eb2.Cross(l));
		Vec3 Tmpl3(eb3.Cross(l));
		Vec3 Tmpmmlp(mm.Cross(lp));
		
		for (unsigned int i = 1; i <= 3; i++) {

			/* Vincolo in rotazione: Delta gb */
			d = Tmpl2.dGet(i);
			WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+1, 3+i, d);

			/* Reazione vincolare: Delta M */
			WM.DecCoef(3+i, 6*(1+Beam::NUMNODES)+1+3+1, d);
			WM.IncCoef(6*activeNode+3+i,
					6*(1+Beam::NUMNODES)+1+3+1, d*dW[0]);
			if (dW[1] != 0) {
				WM.IncCoef(6*(activeNode+1)+3+i,
						6*(1+Beam::NUMNODES)+1+3+1, 
						d*dW[1]);
			}

			/* Vincolo in rotazione: Delta gb */
			d = Tmpl3.dGet(i);
			WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+2, 3+i, d);

			/* Reazione vincolare: Delta M */
			WM.DecCoef(3+i, 6*(1+Beam::NUMNODES)+1+3+2, d);
			WM.IncCoef(6*activeNode+3+i,
					6*(1+Beam::NUMNODES)+1+3+2, d*dW[0]);
			if (dW[1] != 0) {
				WM.IncCoef(6*(activeNode+1)+3+i,
						6*(1+Beam::NUMNODES)+1+3+2, 
						d*dW[1]);
			}

			/* Reazione vincolare: Delta s */
			d = Tmpmmlp(i);
			WM.DecCoef(3+i, 6*(1+Beam::NUMNODES)+1, d);
			WM.IncCoef(6*activeNode+3+i, 
					6*(1+Beam::NUMNODES)+1, d*dW[0]);
			if (dW[1] != 0) {
				WM.IncCoef(6*(activeNode+1)+3+i, 
						6*(1+Beam::NUMNODES)+1, 
						d*dW[1]);
			}
		}

		/* Vincolo in rotazione: Delta s */
		d = (eb2*lp)/dCoef;
		WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+1,
				6*(1+Beam::NUMNODES)+1, d);

		d = (eb3*lp)/dCoef;
		WM.DecCoef(6*(1+Beam::NUMNODES)+1+3+2,
				6*(1+Beam::NUMNODES)+1, d);

		/* Reazione vincolare: Delta gb */
		Mat3x3 mmTmp(MatCrossCross, l, mm*dCoef);
		WM.Sub(3+1, 3+1, mmTmp);
		if (dW[1] == 0) {
			WM.Add(6*activeNode+3+1, 3+1, mmTmp);
		} else {
			WM.Add(6*activeNode+3+1, 3+1, mmTmp*dW[0]);
			WM.Add(6*(activeNode+1)+3+1, 3+1, mmTmp*dW[1]);
		}
		
	}
	
	return WorkMat;
}	

/* Assemblaggio residuo */
SubVectorHandler& 
BeamSliderJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering BeamSliderJoint::AssRes()" << std::endl);

	/*
	 * Nota: posso risparmiare tutte le righe dei nodi della trave
	 * che non sono attivi ...
	 */
	
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	const StructNode *pBeamNode[Beam::NUMNODES];
	
	/* Aggiorna i dati propri */
	sRef = XCurr(iFirstReactionIndex+1);
	F = Vec3(XCurr, iFirstReactionIndex+1+1);
	switch (iType) {
		/*
		 * m(2), m(3) are the moments about the axes 
		 * orthogonal to the tangent to the reference line
		 */
	case CLASSIC:
		m.Put(2, XCurr(iFirstReactionIndex+1+3+1));
		m.Put(3, XCurr(iFirstReactionIndex+1+3+2));
		break;
		
		/*
		 * m(1) is the moment about the tangent to the 
		 * reference line;
		 * m(2), m(3) are the moments about the axes 
		 * orthogonal to the tangent to the reference line
		 */
	case SPLINE:
		m = Vec3(XCurr, iFirstReactionIndex+1+3+1);
		break;

		/*
		 * No moment
		 */
	default:
		/* M is set to zero by someone else ... */
		break;
	}

	/*
	 * in base al valore di s decide su quale trave sta operando
	 * (da studiare e implementare ...)
	 *
	 * Nota: passando da una trave all'altra non e' detto che la
	 * metrica sia la stessa (se hanno lunghezze diverse o i nodi 
	 * non sono equispaziati, cambia).
	 * In prima approssimazione faccio finta che sia la stessa;
	 * un raffinamento si potra' avere considerando il rapporto
	 * tra le metriche.
	 */
	s = sRef - 2*iCurrBeam;
	if (s < -1.) {
		/* passo alla trave precedente */
		if (iCurrBeam > 0) {
			s += 2.;
			iCurrBeam--;
		}

	} else if (s > 1.) {
		/* passo alla trave successiva */
		if (iCurrBeam < nBeams-1) {
			s -= 2.;
			iCurrBeam++;
		}
	}

	/* Cerco il tratto di trave a cui le forze si applicano ... */
	/* Primo tratto */
	if (s < -dS - dL) {
		activeNode = 1;

		dW[0] = 1.;
		dW[1] = 0.;

	} else if ( s < -dS + dL) {
		activeNode = 1;

		doublereal d = .5*(dS + s)/dL;
		dW[0] = .5 - d;
		dW[1] = .5 + d;

	/* Ultimo tratto */
	} else if (s > dS + dL) {
		activeNode = 3;

		dW[0] = 1.;
		dW[1] = 0.;

	} else if (s > dS - dL) {
		activeNode = 2;

		doublereal d = .5*(dS - s)/dL;
		dW[0] = .5 + d;
		dW[1] = .5 - d;

	/* Tratto centrale */
	} else {
		activeNode = 2;

		dW[0] = 1.;
		dW[1] = 0.;
	}

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstMomIndex+iCnt);
	}
	for (int nCnt = 1; nCnt <= Beam::NUMNODES; nCnt++) {
		pBeamNode[nCnt-1] = ppBeam[iCurrBeam]->pGetNode(nCnt);
		integer iBeamFirstMomIndex = 
			pBeamNode[nCnt-1]->iGetFirstMomentumIndex();

		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(6*nCnt+iCnt, 
					iBeamFirstMomIndex+iCnt);
		}
	}
	
	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WorkVec.PutRowIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
	}
	
	/*
	 * Recupero dati
	 */
	x = Vec3(Zero3);
	l = Vec3(Zero3);
	for (unsigned int i = 0; i < Beam::NUMNODES; i++) {
		xNod[i] = pBeamNode[i]->GetXCurr();
		fTmp[i] = pBeamNode[i]->GetRCurr()*ppBeam[iCurrBeam]->Getf(i+1);
		xTmp[i] = xNod[i]+fTmp[i];

		dN[i] = ShapeFunc3N(s, i+1);
		dNp[i] = ShapeFunc3N(s, i+1, ORD_D1);
		dNpp[i] = ShapeFunc3N(s, i+1, ORD_D2);
		x += xTmp[i]*dN[i];
		l += xTmp[i]*dNp[i];
	}
	
	Rb = pNode->GetRCurr()*R;
	fb = pNode->GetRCurr()*f;
	xc = pNode->GetXCurr()+fb;

	Vec3 eb2 = Rb.GetVec(2);
	Vec3 eb3 = Rb.GetVec(3);

	/*
	 * vincoli di posizione
	 *
	 * FIXME: togliere la scalatura da F*l ?????
	 */
	WorkVec.PutCoef(6*(1+Beam::NUMNODES)+1, (F*l)/dCoef);
	WorkVec.Add(6*(1+Beam::NUMNODES)+1+1, (xc-x)/dCoef);

	/*
	 * Vincoli di rotazione
	 */
	if (iType != SPHERICAL) {
		/* 2 vincoli di rotazione */
		WorkVec.PutCoef(6*(1+Beam::NUMNODES)+1+3+1, (eb2*l)/dCoef);
		WorkVec.PutCoef(6*(1+Beam::NUMNODES)+1+3+2, (eb3*l)/dCoef);

		/* calcolo momento */
		M = eb2.Cross(l*m(2))+eb3.Cross(l*m(3));

		if (iType == SPLINE) {
			/* FIXME: vincolo spline */
			NO_OP;
		}
	}

	/*
	 * reazioni vincolari
	 */
	WorkVec.Add(1, F);
	WorkVec.Add(3+1, M+fb.Cross(F));

	WorkVec.Sub(6*activeNode+1, F*dW[0]);
	WorkVec.Sub(6*activeNode+3+1, 
			(M+(xc-xNod[activeNode-1]).Cross(F))*dW[0]);

	if (dW[1] != 0.) {
		WorkVec.Sub(6*(activeNode+1)+1, F*dW[1]);
		WorkVec.Sub(6*(activeNode+1)+3+1, 
				(M+(xc-xNod[activeNode]).Cross(F))*dW[1]);
	}

	return WorkVec;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler &
BeamSliderJoint::InitialAssJac(
		VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr
)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler &
BeamSliderJoint::InitialAssRes(
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr
)
{
	WorkVec.Resize(0);

	return WorkVec;
}

/* BeamSliderJoint - end */

