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
		const Vec3& fTmp, const Mat3x3& RTmp, flag fOut)
: Elem(uL, Elem::JOINT, fOut),
Joint(uL, Joint::BEAMSLIDER, pDO, fOut),
nRotConstr(0), nBeams(nB), iCurrBeam(0), iType(iT),
pNode(pN), ppBeam(ppB),
f(fTmp), R(RTmp),
F(0.), M(0.),
s(0.), activeNode(Beam::NODE2),
x(0.), l(0.)
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

ostream& 
BeamSliderJoint::Restart(ostream& out) const
{
	return out << "# beam slider not implemented yet" << endl
		<< "beam slider;" << endl;
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
			<< " " << s << endl;
	}

	cerr << "################" << endl;
}

doublereal 
BeamSliderJoint::dGetPrivData(unsigned int i) const
{
	return 0.;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
BeamSliderJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering BeamSliderJoint::AssJac()" << endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(iNumRows, iNumCols, 0.);

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
		WM.fPutRowIndex(iCnt, iNodeFirstMomIndex+iCnt);
		WM.fPutColIndex(iCnt, iNodeFirstPosIndex+iCnt);
	}
	
	for (int nCnt = 1; nCnt <= Beam::NUMNODES; nCnt++) {
		pBeamNode[nCnt-1] = ppBeam[iCurrBeam]->pGetNode(nCnt);
		integer iBeamFirstMomIndex = 
			pBeamNode[nCnt-1]->iGetFirstMomentumIndex();
		integer iBeamFirstPosIndex = 
			pBeamNode[nCnt-1]->iGetFirstPositionIndex();

		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.fPutRowIndex(6*nCnt+iCnt, 
					iBeamFirstMomIndex+iCnt);
			WM.fPutColIndex(6*nCnt+iCnt, 
					iBeamFirstPosIndex+iCnt);
		}
	}
	
	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WM.fPutRowIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
		WM.fPutColIndex(6*(1+Beam::NUMNODES)+iCnt, 
				iFirstReactionIndex+iCnt);
	}

	/* vincolo in posizione */
	for (unsigned int i = 1; i <= 3; i++) {
		/* l^T v = 0 : Delta v */
		doublereal d = l.dGet(i)/dCoef;
		WM.fDecCoef(6*(1+Beam::NUMNODES)+1, 
				6*(1+Beam::NUMNODES)+1+i, d);

		/* xb - x = 0: l Delta s */
		WM.fIncCoef(6*(1+Beam::NUMNODES)+1+i, 
				6*(1+Beam::NUMNODES)+1, d);

		/* xb - x = 0: Delta x_b */
		WM.fDecCoef(6*(1+Beam::NUMNODES)+1+i, i, 1.);
	}

	Vec3 lp(0.);
	for (unsigned int iN = 0; iN < Beam::NUMNODES; iN++) {
		Vec3 Tmp(fTmp[iN].Cross(F));

		/* l^T v = 0 : Delta s */
		lp += xTmp[iN]*dNpp[iN];

		for (unsigned int i = 1; i <= 3; i++) {
			/* l^T v = 0 : Delta x */
			doublereal d = F.dGet(i)*dNp[iN];
			WM.fDecCoef(6*(1+Beam::NUMNODES)+1, 
					6*(1+iN)+i, d);

			/* l^T v = 0 : Delta g */
			d = Tmp.dGet(i)*dNp[iN];
			WM.fDecCoef(6*(1+Beam::NUMNODES)+1, 
					6*(1+iN)+3+i, d);

			/* xb - x = 0: Delta x */
			WM.fIncCoef(6*(1+Beam::NUMNODES)+1+i, 
					6*(1+iN)+i, dN[iN]);
		}

		/* xb - x = 0: Delta g */
		WM.Sub(6*(1+Beam::NUMNODES)+1+1, 
				6*(1+iN)+3+1, Mat3x3(fTmp[iN]*dN[iN]));
	}

	/* l^T v = 0 : Delta s */
	WM.fDecCoef(6*(1+Beam::NUMNODES)+1, 
			6*(1+Beam::NUMNODES)+1, (F*lp)/dCoef);

	/* reazioni vincolari */
	for (unsigned int i = 1; i <= 3; i++) {
		/* corpo: Delta v */
		WM.fDecCoef(i, 6*(1+Beam::NUMNODES)+1+i, 1.);

		/* trave: Delta v */
		WM.fIncCoef(6*activeNode+i, 6*(1+Beam::NUMNODES)+1+i, 1.);
	}

	/* trave: Delta v (momento) */
	Mat3x3 MTmp(F*dCoef);
	WM.Add(6*activeNode+3+1, 6*(1+Beam::NUMNODES)+1+1, 
			Mat3x3(pNode->GetXCurr()-xTmp[activeNode-1]));
	WM.Sub(6*activeNode+3+1, 1, MTmp);
	WM.Add(6*activeNode+3+1, 6*activeNode+1, MTmp);

	cerr << WM << endl;
	
	return WorkMat;
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
	
	/* Aggiorna i dati propri */
	s = XCurr.dGetCoef(iFirstReactionIndex+1);
	F = Vec3(XCurr, iFirstReactionIndex+2);
	switch (iType) {
		/*
		 * M(2), M(3) are the moments about the axes 
		 * orthogonal to the tangent to the reference line
		 */
	case CLASSIC:
		M.Put(1, 0.);
		M.Put(2, XCurr.dGetCoef(iFirstReactionIndex+5));
		M.Put(3, XCurr.dGetCoef(iFirstReactionIndex+6));
		break;
		
		/*
		 * M(1) is the moment about the tangent to the 
		 * reference line;
		 * M(2), M(3) are the moments about the axes 
		 * orthogonal to the tangent to the reference line
		 */
	case SPLINE:
		M = Vec3(XCurr, iFirstReactionIndex+5);
		break;

		/*
		 * No moment
		 */
	default:
		M = Vec3(0.);
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
	
	/*
	 * Recupero dati
	 */
	x = Vec3(0.);
	l = Vec3(0.);
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
	
	/*
	 * vincoli di posizione 
	 */
	WorkVec.fPutCoef(6*(1+Beam::NUMNODES)+1, (F*l)/dCoef);
	WorkVec.Add(6*(1+Beam::NUMNODES)+1+1, (pNode->GetXCurr()-x)/dCoef);

	/*
	 * reazioni vincolari
	 */
	WorkVec.Add(1, F);
	WorkVec.Add(4, M);

	/* Cerco il tratto di trave a cui le forze si applicano ... */
	/* Primo tratto */
	if (s < -dS) {
		activeNode = 1;

	/* Ultimo tratto */
	} else if (s > dS) {
		activeNode = 3;

	/* Tratto centrale */
	} else {
		activeNode = 2;
	}

	WorkVec.Sub(6*activeNode+1, F);
	WorkVec.Sub(6*activeNode+3+1, 
			M+(pNode->GetXCurr()-xNod[activeNode-1]).Cross(F));

	cerr << "	iCurrBeam: " << iCurrBeam 
		<< "; activeNode: " << activeNode 
		<< "; s: " << s << endl
		<< "	F: " << F << "; M: " << M << endl
		<< "	x: " << x << "; l: " << l
		<< "; F*l: " << F*l << endl;

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

