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

#include <limits>
#include <cfloat>
#include <limits>

#include "dataman.h"
#include "constltp.h"
#include "shapefnc.h"
#include "beam.h"
#include "pzbeam.h"
#include "Rot.hh"

/*
 * Nota: non e' ancora stato implementato il contributo
 * della ViscoElasticBeam all'assemblaggio iniziale
 */

/*
 * Nota: la parte viscoelastica va rivista in accordo con la piu'
 * recente formulazione delle derivate delle deformazioni nel sistema
 * materiale
 */

/* Beam - begin */

/* Costruttore normale */
Beam::Beam(unsigned int uL,
	const StructNode* pN1,
	const StructNode* pN2,
	const StructNode* pN3,
	const Vec3& F1,
	const Vec3& F2,
	const Vec3& F3,
	const Mat3x3& R1,
	const Mat3x3& R2,
	const Mat3x3& R3,
	const Mat3x3& r_I,
	const Mat3x3& rII,
	const ConstitutiveLaw6D* pD_I,
	const ConstitutiveLaw6D* pDII,
	OrientationDescription ood,
	flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
od(ood),
f(),
RNode(),
bConsistentInertia(false),
dMass_I(0.),
S0_I(Zero3),
J0_I(Zero3x3),
dMassII(0.),
S0II(Zero3),
J0II(Zero3x3),
bFirstRes(false)
{
	ASSERT(pN1 != NULL);
	ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN2 != NULL);
	ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN3 != NULL);
	ASSERT(pN3->GetNodeType() == Node::STRUCTURAL);

	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	pNode[NODE3] = pN3;

	const_cast<Vec3&>(f[NODE1]) = F1;
	const_cast<Vec3&>(f[NODE2]) = F2;
	const_cast<Vec3&>(f[NODE3]) = F3;
	const_cast<Mat3x3&>(RNode[NODE1]) = R1;
	const_cast<Mat3x3&>(RNode[NODE2]) = R2;
	const_cast<Mat3x3&>(RNode[NODE3]) = R3;
	RRef[S_I] = R[S_I] = (Mat3x3&)r_I;
	RRef[SII] = R[SII] = (Mat3x3&)rII;

	pD[S_I] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[S_I],
		ConstitutiveLaw6DOwner,
		ConstitutiveLaw6DOwner(pD_I));
	pD[SII] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[SII],
		ConstitutiveLaw6DOwner,
		ConstitutiveLaw6DOwner(pDII));

	Init();
}


/* Costruttore per la trave con forze d'inerzia consistenti */
Beam::Beam(unsigned int uL,
	const StructNode* pN1,
	const StructNode* pN2,
	const StructNode* pN3,
	const Vec3& F1,
	const Vec3& F2,
	const Vec3& F3,
	const Mat3x3& R1,
	const Mat3x3& R2,
	const Mat3x3& R3,
	const Mat3x3& r_I,
	const Mat3x3& rII,
	const ConstitutiveLaw6D* pD_I,
	const ConstitutiveLaw6D* pDII,
	doublereal dM_I,
	const Vec3& s0_I,
	const Mat3x3& j0_I,
	doublereal dMII,
	const Vec3& s0II,
	const Mat3x3& j0II,
	OrientationDescription ood,
	flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
od(ood),
f(),
RNode(),
bConsistentInertia(true),
dMass_I(dM_I),
S0_I(s0_I),
J0_I(j0_I),
dMassII(dMII),
S0II(s0II),
J0II(j0II),
bFirstRes(false)
{
	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	pNode[NODE3] = pN3;
	const_cast<Vec3&>(f[NODE1]) = F1;
	const_cast<Vec3&>(f[NODE2]) = F2;
	const_cast<Vec3&>(f[NODE3]) = F3;
	const_cast<Mat3x3&>(RNode[NODE1]) = R1;
	const_cast<Mat3x3&>(RNode[NODE2]) = R2;
	const_cast<Mat3x3&>(RNode[NODE3]) = R3;
	RPrev[S_I] = RRef[S_I] = R[S_I] = const_cast<Mat3x3&>(r_I);
	RPrev[SII] = RRef[SII] = R[SII] = const_cast<Mat3x3&>(rII);

	pD[S_I] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[S_I],
		ConstitutiveLaw6DOwner,
		ConstitutiveLaw6DOwner(pD_I));
	pD[SII] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[SII],
		ConstitutiveLaw6DOwner,
		ConstitutiveLaw6DOwner(pDII));

	Init();
}

void
Beam::Init(void)
{
	for (unsigned i = 0; i < NUMSEZ; i++) {
#ifdef USE_NETCDF
		Var_X[i] = 0;
		Var_Phi[i] = 0;
		Var_F[i] = 0;
		Var_M[i] = 0;
		Var_Nu[i] = 0;
		Var_K[i] = 0;
		Var_NuP[i] = 0;
		Var_KP[i] = 0;
#endif /* USE_NETCDF */

		Omega[i] = Zero3;
		Az[i] = Zero6;
		AzRef[i] = Zero6;
		AzLoc[i] = Zero6;
		DefLoc[i] = Zero6;
		DefLocRef[i] = Zero6;
		DefLocPrev[i] = Zero6;
		p[i] = Zero3;
		g[i] = Zero3;
		L0[i] = Zero3;
		L[i] = Zero3;
		LRef[i] = Zero3;
	}

	DsDxi();

	Vec3 xTmp[NUMNODES];

	for (unsigned int i = 0; i < NUMNODES; i++) {
		xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
	}

	/* Aggiorna le grandezze della trave nei punti di valutazione */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		p[iSez] = InterpState(xTmp[NODE1],
			xTmp[NODE2],
			xTmp[NODE3],
			Beam::Section(iSez));
	}
}

Beam::~Beam(void)
{
	ASSERT(pD[S_I] != 0);
	if (pD[S_I] != 0) {
		SAFEDELETE(pD[S_I]);
	}

	ASSERT(pD[SII] != 0);
	if (pD[SII] != 0) {
		SAFEDELETE(pD[SII]);
	}
}

static const unsigned int iNumPrivData =
		+3		//  0 ( 1 ->  3) - "e" strain
		+3		//  3 ( 4 ->  6) - "k" curvature
		+3		//  6 ( 7 ->  9) - "F" force
		+3		//  9 (10 -> 12) - "M" moment
		+3		// 12 (13 -> 15) - "X" position
		+3		// 15 (16 -> 18) - "Phi" orientation vector
		+3		// 18 (19 -> 21) - "Omega" angular velocity
		+3		// 21 (22 -> 24) - "eP" strain rate
		+3		// 24 (25 -> 27) - "kP" curvature rate
	;

/* Accesso ai dati privati */
unsigned int
Beam::iGetNumPrivData(void) const
{
	return 2*iNumPrivData;
}

unsigned int
Beam::iGetPrivDataIdx_int(const char *s, ConstLawType::Type type)
{
	ASSERT(s != NULL);

	unsigned int idx = 0;

	switch (s[0]) {
	case 'k':
		idx += 3;
	case 'e':
		if (s[1] == 'P') {
			if (type != ConstLawType::VISCOUS) {
				return 0;
			}
			s++;
			idx += 21;
			break;
		}
		break;

	case 'F':
		idx += 6;
		break;

	case 'M':
		idx += 9;
		break;

	case 'X':
		idx += 12;
		break;

	case 'P':
		if (strncasecmp(s, "Phi", STRLENOF("Phi")) != 0) {
			return 0;
		}
		s += STRLENOF("Phi") - 1;
		idx += 15;
		break;

	case 'O':
		if (strncasecmp(s, "Omega", STRLENOF("Omega")) != 0) {
			return 0;
		}
		s += STRLENOF("Omega") - 1;
		idx += 18;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		idx += 1;
		break;

	case 'y':
		idx += 2;
		break;

	case 'z':
		idx += 3;
		break;

	default:
		return 0;
	}

	if (s[2] != '\0') {
		return 0;
	}

	return idx;
}

unsigned int
Beam::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int p_idx = 0;

	if (strncmp(s, "pI", STRLENOF("pI")) != 0) {
		return 0;
	}
	s += STRLENOF("pI");

	if (s[0] == 'I') {
		p_idx += iNumPrivData;
		s++;
	}

	if (s[0] != '.') {
		return 0;
	}
	s++;

	ConstLawType::Type type = ConstLawType::ELASTIC;
	if (dynamic_cast<const ViscoElasticBeam *>(this)) {
		type = ConstLawType::VISCOUS;
	}

	unsigned int idx = iGetPrivDataIdx_int(s, type);
	if (idx != 0) {
		idx += p_idx;
	}

	return idx;
}

doublereal
Beam::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= iGetNumPrivData());

	int sez = (i - 1)/iNumPrivData;
	int idx = (i - 1)%iNumPrivData + 1;

	switch (idx) {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:

		return DefLoc[sez].dGet(idx);

	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
	case 12:
		return AzLoc[sez].dGet(idx - 6);

	case 13:
	case 14:
	case 15:
		return p[sez].dGet(idx - 12);

	case 16:
	case 17:
	case 18:
		return RotManip::VecRot(R[sez]).dGet(idx - 15);

	case 19:
	case 20:
	case 21:
		return Omega[sez].dGet(idx - 18);

	default:
		silent_cerr("Beam(" << GetLabel() << "): "
			"illegal private data " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3
Beam::InterpState(const Vec3& v1,
	const Vec3& v2,
	const Vec3& v3,
	enum Section Sec)
{
	doublereal* pv1 = const_cast<doublereal *>(v1.pGetVec());
	doublereal* pv2 = const_cast<doublereal *>(v2.pGetVec());
	doublereal* pv3 = const_cast<doublereal *>(v3.pGetVec());

	return Vec3(pv1[0]*dN3[Sec][0] + pv2[0]*dN3[Sec][1] + pv3[0]*dN3[Sec][2],
		pv1[1]*dN3[Sec][0] + pv2[1]*dN3[Sec][1] + pv3[1]*dN3[Sec][2],
		pv1[2]*dN3[Sec][0] + pv2[2]*dN3[Sec][1] + pv3[2]*dN3[Sec][2]);
}


Vec3
Beam::InterpDeriv(const Vec3& v1,
	const Vec3& v2,
	const Vec3& v3,
	enum Section Sec)
{
	doublereal* pv1 = const_cast<doublereal *>(v1.pGetVec());
	doublereal* pv2 = const_cast<doublereal *>(v2.pGetVec());
	doublereal* pv3 = const_cast<doublereal *>(v3.pGetVec());
	return Vec3((pv1[0]*dN3P[Sec][0] + pv2[0]*dN3P[Sec][1]
			+ pv3[0]*dN3P[Sec][2])*dsdxi[Sec],
		(pv1[1]*dN3P[Sec][0] + pv2[1]*dN3P[Sec][1]
			+ pv3[1]*dN3P[Sec][2])*dsdxi[Sec],
		(pv1[2]*dN3P[Sec][0] + pv2[2]*dN3P[Sec][1]
			+ pv3[2]*dN3P[Sec][2])*dsdxi[Sec]);
}


void
Beam::DsDxi(void)
{
	/* Validazione dati */
	ASSERT(pNode[NODE1] != NULL);
	ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode[NODE2] != NULL);
	ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode[NODE3] != NULL);
	ASSERT(pNode[NODE3]->GetNodeType() == Node::STRUCTURAL);

	/* Calcola il ds/dxi e le deformazioni iniziali */
	Mat3x3 RNod[NUMNODES];
	Vec3 xTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		RNod[i] = pNode[i]->GetRCurr();
		xTmp[i] = pNode[i]->GetXCurr() + RNod[i]*f[i];
	}

	Vec3 xGrad[NUMSEZ];
	for (unsigned int i = 0; i < NUMSEZ; i++) {
		/* temporary */
		dsdxi[i] = 1.;
		xGrad[i] = InterpDeriv(xTmp[NODE1],
			xTmp[NODE2],
			xTmp[NODE3],
			Beam::Section(i));
		doublereal d = xGrad[i].Dot();
		if (d > std::numeric_limits<doublereal>::epsilon()) {
			/* final */
			dsdxi[i] = 1./std::sqrt(d);

		} else {
			silent_cerr("warning, beam " << GetLabel()
				<< " has singular metric; aborting..." << std::endl);

			throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
		}
	}

	/* Calcola le deformazioni iniziali */
	for (unsigned int i = 0; i < NUMSEZ; i++) {
		L0[i] = R[i].MulTV(InterpDeriv(xTmp[NODE1],
			xTmp[NODE2],
			xTmp[NODE3],
			Beam::Section(i)));
		pD[i]->Update(Zero6);
		DRef[i] = MultRMRt(pD[i]->GetFDE(), R[i]);
	}
}


/* Calcola la velocita' angolare delle sezioni a partire da quelle dei nodi */
void
Beam::Omega0(void)
{
	/* Validazione dati */
	ASSERT(pNode[NODE1] != NULL);
	ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode[NODE2] != NULL);
	ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode[NODE3] != NULL);
	ASSERT(pNode[NODE3]->GetNodeType() == Node::STRUCTURAL);

	/* Modo consistente: */
	Mat3x3 RNod[NUMNODES];
	Vec3 w[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		RNod[i] = pNode[i]->GetRCurr()*RNode[i];
		w[i] = pNode[i]->GetWCurr();
	}

	/* Calcolo i parametri di rotazione dei nodi di estremo rispetto a quello
	 * centrale, nell'ipotesi (realistica) che siano limitati */
	Vec3 g1(CGR_Rot::Param, RNod[NODE2].MulTM(RNod[NODE1]));
	Vec3 g3(CGR_Rot::Param, RNod[NODE2].MulTM(RNod[NODE3]));

	/* Calcolo le derivate dei parametri di rotazione ai nodi di estremo; in
	 * quello centrale si assume che sia uguale alla velocita' di rotazione */
	Vec3 g1P(Mat3x3(CGR_Rot::MatGm1, g1)*w[NODE1]);
	Vec3 g3P(Mat3x3(CGR_Rot::MatGm1, g3)*w[NODE3]);

	for (unsigned int i = 0; i < NUMSEZ; i++) {
		Vec3 gTmp(g1*dN3[i][NODE1]+g3*dN3[i][NODE3]);
		Vec3 gPTmp(g1P*dN3[i][NODE1]+w[NODE2]*dN3[i][NODE2]+g3P*dN3[i][NODE3]);
		Omega[i] = Mat3x3(CGR_Rot::MatG, gTmp)*gPTmp;
	}

#if 0
	/* Modo brutale: interpolo le velocita' dei nodi */
	Vec3 w[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		w[i] = pNode[i]->GetWCurr();
	}
	for (unsigned int i = 0; i < NUMSEZ; i++) {
	Omega[i] = (w{NODE1]*dN3[i][NODE1]
		+ w[NODE2]*dN3[i][NODE2]
		+ w[NODE3]*dN3[i][NODE3]);
	}
#endif /* 0 */
}


/* Contributo al file di restart */
std::ostream&
Beam::Restart(std::ostream& out) const
{
	return Restart_(out)<< ';' << std::endl;
}

std::ostream&
Beam::Restart_(std::ostream& out) const
{
	out << "  beam: " << GetLabel();
	for (unsigned int i = 0; i < NUMNODES; i++) {
		out << ", " << pNode[i]->GetLabel() << ", reference, node, ",
			f[i].Write(out, ", ");
	}

	for (unsigned int i = 0; i < NUMSEZ; i++) {
		out << ", reference, global, 1, ",
			(R[i].GetVec(1)).Write(out, ", ")
			<< ", 2, ", (R[i].GetVec(2)).Write(out, ", ") << ", ",
			pD[i]->pGetConstLaw()->Restart(out);
	}

	return out;
}

void
Beam::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	for (unsigned i = 0; i < NUMSEZ; i++) {
		RPrev[i] = R[i];
		DefLocPrev[i] = DefLoc[i];
		pD[i]->AfterConvergence(DefLoc[i]);
	}
}

/* Inverse Dynamics: */
void
Beam::AfterConvergence(const VectorHandler& X, const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

/* Assembla la matrice */
void
Beam::AssStiffnessMat(FullSubMatrixHandler& WMA,
			   FullSubMatrixHandler& /* WMB */ ,
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("Beam::AssStiffnessMat");

	/* La matrice arriva gia' dimensionata e con gli indici di righe e colonne
	 * a posto */

	/* offset nel riferimento globale */
	Vec3 fTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
	}

	Mat6x6 AzTmp[NUMSEZ][NUMNODES];

	/* per ogni punto di valutazione: */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		for (unsigned int i = 0; i < NUMNODES; i++) {
			/* Delta - deformazioni */
			AzTmp[iSez][i] = Mat6x6(mb_deye<Mat3x3>(dN3P[iSez][i]*dsdxi[iSez]*dCoef),
				Zero3x3,
				Mat3x3(MatCross, L[iSez]*(dN3[iSez][i]*dCoef) - fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez]*dCoef)),
				mb_deye<Mat3x3>(dN3P[iSez][i]*dsdxi[iSez]*dCoef));
			/* Delta - azioni interne */
			AzTmp[iSez][i] = DRef[iSez]*AzTmp[iSez][i];
			/* Correggo per la rotazione da locale a globale */
			AzTmp[iSez][i].SubMat12(Mat3x3(MatCross, Az[iSez].GetVec1()*(dN3[iSez][i]*dCoef)));
			AzTmp[iSez][i].SubMat22(Mat3x3(MatCross, Az[iSez].GetVec2()*(dN3[iSez][i]*dCoef)));
		}
	} /* end ciclo sui punti di valutazione */

	Vec3 bTmp[2];

	/* ciclo sulle equazioni */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		bTmp[0] = p[iSez] - pNode[iSez]->GetXCurr();
		bTmp[1] = p[iSez] - pNode[iSez + 1]->GetXCurr();

		unsigned int iRow1 = iSez*6 + 1;
		unsigned int iRow2 = iRow1 + 6;

		for (unsigned int i = 0; i < NUMNODES; i++) {
			/* Equazione all'indietro: */
			WMA.Sub(iRow1, 6*i + 1, AzTmp[iSez][i].GetMat11());
			WMA.Sub(iRow1, 6*i + 4, AzTmp[iSez][i].GetMat12());

			WMA.Sub(iRow1 + 3, 6*i + 1,
				AzTmp[iSez][i].GetMat21()
				- Mat3x3(MatCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
				+ bTmp[0].Cross(AzTmp[iSez][i].GetMat11()));
			WMA.Sub(iRow1 + 3, 6*i + 4,
				AzTmp[iSez][i].GetMat22()
				- Mat3x3(MatCrossCross, Az[iSez].GetVec1()*(-dCoef*dN3[iSez][i]), fTmp[i])
				+ bTmp[0].Cross(AzTmp[iSez][i].GetMat12()));

			/* Equazione in avanti: */
			WMA.Add(iRow2, 6*i+1, AzTmp[iSez][i].GetMat11());
			WMA.Add(iRow2, 6*i+4, AzTmp[iSez][i].GetMat12());

			WMA.Add(iRow2 + 3, 6*i + 1,
				AzTmp[iSez][i].GetMat21()
				- Mat3x3(MatCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
				+ bTmp[1].Cross(AzTmp[iSez][i].GetMat11()));
			WMA.Add(iRow2 + 3, 6*i + 4,
				AzTmp[iSez][i].GetMat22()
				+ Mat3x3(MatCrossCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]), fTmp[i])
				+ bTmp[1].Cross(AzTmp[iSez][i].GetMat12()));
		}

		/* correzione alle equazioni */
		Mat3x3 FTmp(MatCross, Az[iSez].GetVec1()*dCoef);
		WMA.Sub(iRow1 + 3, 6*iSez + 1, FTmp);
		WMA.Add(iRow2 + 3, 6*iSez + 7, FTmp);

	} /* end ciclo sui punti di valutazione */
}


/* Assembla il residuo */
void
Beam::AssStiffnessVec(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("Beam::AssStiffnessVec");

	/* Riceve il vettore gia' dimensionato e con gli indici a posto
	 * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */

	/* Per la trattazione teorica, il riferimento e' il file ul-travi.tex
	 * (ora e' superato) */

	if (bFirstRes) {
		bFirstRes = false; /* AfterPredict ha gia' calcolato tutto */
	} else {
		Vec3 gNod[NUMNODES];
		Vec3 xTmp[NUMNODES];

		for (unsigned int i = 0; i < NUMNODES; i++) {
			gNod[i] = pNode[i]->GetgCurr();
			xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRCurr()*f[i];
		}

		Mat3x3 RDelta[NUMSEZ];
		Vec3 gGrad[NUMSEZ];

		/* Aggiorna le grandezze della trave nei punti di valutazione */
		for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

			/* Posizione */
			p[iSez] = InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

			/* Matrici di rotazione */
			g[iSez] = InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
			RDelta[iSez] = Mat3x3(CGR_Rot::MatR, g[iSez]);
			R[iSez] = RDelta[iSez]*RRef[iSez];

			/* Derivate della posizione */
			L[iSez] = InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

			/* Derivate dei parametri di rotazione */
			gGrad[iSez] = InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));

			/* Calcola le deformazioni nel sistema locale nei punti di valutazione */
			DefLoc[iSez] = Vec6(R[iSez].MulTV(L[iSez]) - L0[iSez],
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]) + DefLocRef[iSez].GetVec2());

			/* Calcola le azioni interne */
			pD[iSez]->Update(DefLoc[iSez]);
			AzLoc[iSez] = pD[iSez]->GetF();

			/* corregge le azioni interne locali (piezo, ecc) */
			AddInternalForces(AzLoc[iSez], iSez);

			/* Porta le azioni interne nel sistema globale */
			Az[iSez] = MultRV(AzLoc[iSez], R[iSez]);
		}
	}

	WorkVec.Add(1, Az[S_I].GetVec1());
	WorkVec.Add(4, (p[S_I] - pNode[NODE1]->GetXCurr()).Cross(Az[S_I].GetVec1())
		+ Az[S_I].GetVec2());
	WorkVec.Add(7, Az[SII].GetVec1() - Az[S_I].GetVec1());
	WorkVec.Add(10, Az[SII].GetVec2() - Az[S_I].GetVec2()
		+ (p[SII] - pNode[NODE2]->GetXCurr()).Cross(Az[SII].GetVec1())
		- (p[S_I] - pNode[NODE2]->GetXCurr()).Cross(Az[S_I].GetVec1()));
	WorkVec.Sub(13, Az[SII].GetVec1());
	WorkVec.Add(16, Az[SII].GetVec1().Cross(p[SII] - pNode[NODE3]->GetXCurr())
		- Az[SII].GetVec2());
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
Beam::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("Beam::AssJac => AssStiffnessMat");

	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
	integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	if (bConsistentInertia) {
		WM.ResizeReset(36, 18);
	} else {
		WM.ResizeReset(18, 18);
	}

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode3FirstMomIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
	}

	AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);

	if (bConsistentInertia) {
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(18 + iCnt, iNode1FirstPosIndex + iCnt);
			WM.PutRowIndex(24 + iCnt, iNode2FirstPosIndex + iCnt);
			WM.PutRowIndex(30 + iCnt, iNode3FirstPosIndex + iCnt);
		}

		AssInertiaMat(WM, WM, dCoef, XCurr, XPrimeCurr);
	}

	return WorkMat;
}


/* assemblaggio matrici per autovalori */
void
Beam::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("Beam::AssMats => AssStiffnessMat");

	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
	integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	if (bConsistentInertia) {
		WMA.ResizeReset(36, 18);
		WMB.ResizeReset(36, 18);
	} else  {
		WMA.ResizeReset(18, 18);
		WorkMatB.SetNullMatrix();
	}

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WMA.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
		WMA.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
		WMA.PutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
		WMA.PutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
	}

	AssStiffnessMat(WMA, WMA, 1., XCurr, XPrimeCurr);

	if (bConsistentInertia) {
		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
			WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
			WMB.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
			WMB.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
			WMB.PutRowIndex(12 + iCnt, iNode3FirstMomIndex + iCnt);
			WMB.PutColIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
		}

		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WMA.PutRowIndex(18 + iCnt, iNode1FirstPosIndex + iCnt);
			WMA.PutRowIndex(24 + iCnt, iNode2FirstPosIndex + iCnt);
			WMA.PutRowIndex(30 + iCnt, iNode3FirstPosIndex + iCnt);

			WMB.PutRowIndex(18 + iCnt, iNode1FirstPosIndex + iCnt);
			WMB.PutRowIndex(24 + iCnt, iNode2FirstPosIndex + iCnt);
			WMB.PutRowIndex(30 + iCnt, iNode3FirstPosIndex + iCnt);
		}

		AssInertiaMat(WMA, WMB, 1., XCurr, XPrimeCurr);
	}
}


/* assemblaggio residuo */
SubVectorHandler&
Beam::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("Beam::AssRes => AssStiffnessVec");

	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	if (bConsistentInertia) {
		WorkVec.ResizeReset(36);
	} else {
		WorkVec.ResizeReset(18);
	}

	for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstMomIndex + iCnt);
	}

	AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	if (bConsistentInertia) {
		integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
		integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
		integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(18 + iCnt, iNode1FirstPosIndex + iCnt);
			WorkVec.PutRowIndex(24 + iCnt, iNode2FirstPosIndex + iCnt);
			WorkVec.PutRowIndex(30 + iCnt, iNode3FirstPosIndex + iCnt);
		}

		AssInertiaVec(WorkVec, dCoef, XCurr, XPrimeCurr);
	}

	return WorkVec;
}

/* Inverse Dynamics residual assembly */
SubVectorHandler&
Beam::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler&  XCurr ,
	const VectorHandler&  XPrimeCurr ,
	const VectorHandler&  XPrimePrimeCurr ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("Beam::AssRes => AssStiffnessVec");

	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstPositionIndex();
	integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstPositionIndex();

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	WorkVec.ResizeReset(18);

	for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstMomIndex + iCnt);
	}

	AssStiffnessVec(WorkVec, 1., XCurr, XPrimeCurr);

	return WorkVec;
}

/* inverse dynamics capable element */
bool
Beam::bInverseDynamics(void) const
{
	return true;
}

/* Settings iniziali, prima della prima soluzione */
void
Beam::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		RRef[iSez] = R[iSez];
		LRef[iSez] = L[iSez];
		DefLocRef[iSez] = DefLoc[iSez];
		AzRef[iSez] = Az[iSez];

		/* Aggiorna il legame costitutivo di riferimento
		 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo) */
		DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), RRef[iSez]);
	}

	bFirstRes = true;
}


/* Prepara i parametri di riferimento dopo la predizione */
void
Beam::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */

	/* Recupera i dati dei nodi */
	Vec3   gNod[NUMNODES];
	Vec3   xTmp[NUMNODES];

	for (unsigned int i = 0; i < NUMNODES; i++) {
		gNod[i] = pNode[i]->GetgRef();
		xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRRef()*f[i];
	}

	Mat3x3 RDelta[NUMSEZ];
	Vec3 gGrad[NUMSEZ];

	/* Aggiorna le grandezze della trave nei punti di valutazione */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

		/* Posizione */
		p[iSez] = InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

		/* Matrici di rotazione */
		g[iSez] = InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
		RDelta[iSez] = Mat3x3(CGR_Rot::MatR, g[iSez]);
		R[iSez] = RRef[iSez] = RDelta[iSez]*RPrev[iSez];

		/* Derivate della posizione */
		L[iSez] = LRef[iSez]
			= InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

		/* Derivate dei parametri di rotazione */
		gGrad[iSez]
			= InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));

		/* Calcola le deformazioni nel sistema locale nei punti di valutazione */
		DefLoc[iSez] = DefLocRef[iSez]
			= Vec6(R[iSez].MulTV(L[iSez]) - L0[iSez],
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]) + DefLocPrev[iSez].GetVec2());

		/* Calcola le azioni interne */
		pD[iSez]->Update(DefLoc[iSez]);
		AzLoc[iSez] = pD[iSez]->GetF();

		/* corregge le azioni interne locali (piezo, ecc) */
		AddInternalForces(AzLoc[iSez], iSez);

		/* Porta le azioni interne nel sistema globale */
		Az[iSez] = AzRef[iSez] = MultRV(AzLoc[iSez], R[iSez]);

		/* Aggiorna il legame costitutivo di riferimento */
		DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), RRef[iSez]);
	}

	bFirstRes = true;
}

void
Beam::OutputPrepare(OutputHandler &OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::BEAMS)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			const char *type = 0;
			switch (GetBeamType()) {
			case Beam::ELASTIC:
				type = "elastic";
				break;

			case Beam::VISCOELASTIC:
				type = "viscoelastic";
				break;

			case Beam::PIEZOELECTRICELASTIC:
				type = "piezoelectric elastic";
				break;

			case Beam::PIEZOELECTRICVISCOELASTIC:
				type = "piezoelectric viscoelastic";
				break;

			default:
				type = "unknown";
				break;
			}

			std::ostringstream os;
			os << "elem.beam." << GetLabel();

			(void)OH.CreateVar(os.str(), type);

			os << '.';
			std::string name(os.str());

			unsigned uOutputFlags = (fToBeOutput() & ToBeOutput::OUTPUT_PRIVATE_MASK);

			static const char *sez[] = { "I", "II" };
			for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
				os.str("");
				os << "evaluation point " << sez[iSez] << " ";
				std::string ep(os.str());

				if (uOutputFlags & Beam::OUTPUT_EP_X) {
					Var_X[iSez] = OH.CreateVar<Vec3>(name + "X_" + sez[iSez], "m",
						ep + "global position vector (X, Y, Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_R) {
					Var_Phi[iSez] = OH.CreateRotationVar(name,
						std::string("_") + sez[iSez], od, ep + "global");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_F) {
					Var_F[iSez] = OH.CreateVar<Vec3>(name + "F_" + sez[iSez], "N",
						ep + "internal force in local frame (F_X, F_Y, F_Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_M) {
					Var_M[iSez] = OH.CreateVar<Vec3>(name + "M_" + sez[iSez], "Nm",
						ep + "internal moment in local frame (M_X, M_Y, M_Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_NU) {
					Var_Nu[iSez] = OH.CreateVar<Vec3>(name + "nu_" + sez[iSez], "-",
						ep + "linear strain in local frame (nu_X, nu_Y, nu_Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_K) {
					Var_K[iSez] = OH.CreateVar<Vec3>(name + "k_" + sez[iSez], "1/m",
						ep + "angular strain in local frame (K_X, K_Y, K_Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_NUP) {
					Var_NuP[iSez] = OH.CreateVar<Vec3>(name + "nuP_" + sez[iSez], "1/s",
						ep + "linear strain rate in local frame (nuP_X, nuP_Y, nuP_Z)");
				}

				if (uOutputFlags & Beam::OUTPUT_EP_KP) {
					Var_KP[iSez] = OH.CreateVar<Vec3>(name + "kP_" + sez[iSez], "1/ms",
						ep + "angular strain rate in local frame (KP_X, KP_Y, KP_Z)");
				}
			}
		}
#endif // USE_NETCDF
	}
}

/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void
Beam::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::BEAMS)) {
			for (unsigned iSez = 0; iSez < NUMSEZ; iSez++) {
				if (Var_X[iSez]) {
					Var_X[iSez]->put_rec(p[iSez].pGetVec(), OH.GetCurrentStep());
				}

				if (Var_Phi[iSez]) {
					Vec3 E;
					switch (od) {
					case EULER_123:
						E = MatR2EulerAngles123(R[iSez])*dRaDegr;
						break;

					case EULER_313:
						E = MatR2EulerAngles313(R[iSez])*dRaDegr;
						break;

					case EULER_321:
						E = MatR2EulerAngles321(R[iSez])*dRaDegr;
						break;

					case ORIENTATION_VECTOR:
						E = RotManip::VecRot(R[iSez]);
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
						Var_Phi[iSez]->put_rec(E.pGetVec(), OH.GetCurrentStep());
						break;

					case ORIENTATION_MATRIX:
						Var_Phi[iSez]->put_rec(R[iSez].pGetMat(), OH.GetCurrentStep());
						break;

					default:
						/* impossible */
						break;
					}
				}

				if (Var_F[iSez]) {
					Var_F[iSez]->put_rec(AzLoc[iSez].GetVec1().pGetVec(), OH.GetCurrentStep());
				}

				if (Var_M[iSez]) {
					Var_M[iSez]->put_rec(AzLoc[iSez].GetVec2().pGetVec(), OH.GetCurrentStep());
				}

				if (Var_Nu[iSez]) {
					Var_Nu[iSez]->put_rec(DefLoc[iSez].GetVec1().pGetVec(), OH.GetCurrentStep());
				}

				if (Var_K[iSez]) {
					Var_K[iSez]->put_rec(DefLoc[iSez].GetVec2().pGetVec(), OH.GetCurrentStep());
				}

				if (Var_NuP[iSez]) {
					Var_NuP[iSez]->put_rec(DefPrimeLoc[iSez].GetVec1().pGetVec(), OH.GetCurrentStep());
				}

				if (Var_KP[iSez]) {
					Var_KP[iSez]->put_rec(DefPrimeLoc[iSez].GetVec2().pGetVec(), OH.GetCurrentStep());
				}
			}
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::BEAMS)) {
			OH.Beams() << std::setw(8) << GetLabel()
				<< " " << AzLoc[S_I].GetVec1()
				<< " " << AzLoc[S_I].GetVec2()
				<< " " << AzLoc[SII].GetVec1()
				<< " " << AzLoc[SII].GetVec2()
				<< std::endl;
		}
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
Beam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("Beam::InitialAssJac => AssStiffnessMat");

	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(18, 18);

	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
	}

	AssStiffnessMat(WM, WM, 1., XCurr, XCurr);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler& Beam::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("Beam::InitialAssRes => AssStiffnessVec");

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	WorkVec.ResizeReset(18);

	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
	}

	AssStiffnessVec(WorkVec, 1., XCurr, XCurr);

	return WorkVec;
}


const StructNode*
Beam::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 3);
	switch (i) {
	case 1:
	case 2:
	case 3:
		return pNode[i-1];
	default:
		throw Beam::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}


void
Beam::GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& r) const
{
	ASSERT(part == 1 || part == 2);
	part--;

	x = p[part];
	r = R[part];
}

void
Beam::GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const
{
	ASSERT(part == 1 || part == 2);
	part--;

	/* FIXME */
	v = Zero3;
	w = Zero3;
}

#ifdef USE_ADAMS
std::ostream&
Beam::WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part, unsigned int firstId) const
{
	Vec3 xTmp[NUMNODES];

	part--;

	for (unsigned int i = part; i <= part+1; i++) {
		xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
	}

	out << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << 1+part << std::endl
		<< firstId << " "
		<< p[part] << " "
		<< MatR2EulerAngles(R[part])*dRaDegr << " "
		<< R[part].MulTV(xTmp[part]-p[part]) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[part]->GetRCurr())*dRaDegr */ << " "
		<< R[part].MulTV(xTmp[1+part]-p[part]) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[1+part]->GetRCurr())*dRaDegr */
		<< std::endl;

	return out;
}
#endif /* USE_ADAMS */

/* Beam - end */


/* ViscoElasticBeam - begin */

/* Costruttore normale */
ViscoElasticBeam::ViscoElasticBeam(
	unsigned int uL,
	const StructNode* pN1,
	const StructNode* pN2,
	const StructNode* pN3,
	const Vec3& F1,
	const Vec3& F2,
	const Vec3& F3,
	const Mat3x3& R1,
	const Mat3x3& R2,
	const Mat3x3& R3,
	const Mat3x3& r_I,
	const Mat3x3& rII,
	const ConstitutiveLaw6D* pD_I,
	const ConstitutiveLaw6D* pDII,
	OrientationDescription ood,
	flag fOut)
: Elem(uL, fOut),
Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut)
{
	Init();
}


/* Costruttore per la trave con forze d'inerzia consistenti */
ViscoElasticBeam::ViscoElasticBeam(
	unsigned int uL,
	const StructNode* pN1,
	const StructNode* pN2,
	const StructNode* pN3,
	const Vec3& F1,
	const Vec3& F2,
	const Vec3& F3,
	const Mat3x3& R1,
	const Mat3x3& R2,
	const Mat3x3& R3,
	const Mat3x3& r_I, const Mat3x3& rII,
	const ConstitutiveLaw6D* pD_I,
	const ConstitutiveLaw6D* pDII,
	doublereal dM_I,
	const Vec3& s0_I, const Mat3x3& j0_I,
	doublereal dMII,
	const Vec3& s0II, const Mat3x3& j0II,
	OrientationDescription ood,
	flag fOut)
: Elem(uL, fOut),
Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII,
	  dM_I, s0_I, j0_I, dMII, s0II, j0II, ood, fOut)
{
	Init();
}

void
ViscoElasticBeam::Init(void)
{
	gPrime[S_I] = gPrime[SII] = Zero3;

	LPrime[S_I] = LPrime[SII] = Zero3;
	LPrimeRef[S_I] = LPrimeRef[SII] = Zero3;

	DefPrimeLoc[S_I] = DefPrimeLoc[SII] = Zero6;
	DefPrimeLocRef[S_I] = DefPrimeLocRef[SII] = Zero6;

	/* Nota: DsDxi() viene chiamata dal costruttore di Beam */
	Beam::Omega0();

	OmegaRef[S_I] = Omega[S_I];
	OmegaRef[SII] = Omega[SII];

	const_cast<Mat6x6&>(ERef[S_I]) = MultRMRt(pD[S_I]->GetFDEPrime(), RRef[S_I]);
	const_cast<Mat6x6&>(ERef[SII]) = MultRMRt(pD[SII]->GetFDEPrime(), RRef[SII]);
}

void
ViscoElasticBeam::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	for (unsigned i = 0; i < NUMSEZ; i++) {
		RPrev[i] = R[i];
		DefLocPrev[i] = DefLoc[i];
		pD[i]->AfterConvergence(DefLoc[i], DefPrimeLoc[i]);
	}
}

/* Inverse Dynamics */
void
ViscoElasticBeam::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

/* Assembla la matrice */
void
ViscoElasticBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
				       FullSubMatrixHandler& WMB,
				       doublereal dCoef,
				       const VectorHandler& /* XCurr */ ,
				       const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticBeam::AssStiffnessMat");

	/* La matrice arriva gia' dimensionata e con gli indici di righe e colonne
	 * a posto */

	/* offset nel riferimento globale */
	Vec3 fTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
	}

	Mat6x6 AzTmp[NUMSEZ][NUMNODES];
	Mat6x6 AzPrimeTmp[NUMSEZ][NUMNODES];

	/* per ogni punto di valutazione: */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		for (unsigned int i = 0; i < NUMNODES; i++) {
			/* Delta - deformazioni */
			AzTmp[iSez][i] = AzPrimeTmp[iSez][i]
				= Mat6x6(mb_deye<Mat3x3>(dN3P[iSez][i]*dsdxi[iSez]),
					Zero3x3,
					Mat3x3(MatCross, L[iSez]*dN3[iSez][i] - fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez])),
					mb_deye<Mat3x3>(dN3P[iSez][i]*dsdxi[iSez]));
			AzTmp[iSez][i] = DRef[iSez]*AzTmp[iSez][i]*dCoef;
			AzTmp[iSez][i] +=
				ERef[iSez]*Mat6x6(Mat3x3(MatCross, Omega[iSez]*(-dN3P[iSez][i]*dsdxi[iSez]*dCoef)),
					Zero3x3,
					(Mat3x3(MatCross, LPrime[iSez]) - Mat3x3(MatCrossCross, Omega[iSez], L[iSez]))*(dN3[iSez][i]*dCoef)
						+ Mat3x3(MatCrossCross, Omega[iSez], fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez]*dCoef))
						+ Mat3x3(MatCross, fTmp[i].Cross(pNode[i]->GetWCurr()*(dN3P[iSez][i]*dsdxi[iSez]*dCoef))),
					Mat3x3(MatCross, Omega[iSez]*(-dN3P[iSez][i]*dsdxi[iSez]*dCoef)));
			AzPrimeTmp[iSez][i] = ERef[iSez]*AzPrimeTmp[iSez][i];
	
			/* Correggo per la rotazione da locale a globale */
			AzTmp[iSez][i].SubMat12(Mat3x3(MatCross, Az[iSez].GetVec1()*(dN3[iSez][i]*dCoef)));
			AzTmp[iSez][i].SubMat22(Mat3x3(MatCross, Az[iSez].GetVec2()*(dN3[iSez][i]*dCoef)));
		}
	} /* end ciclo sui punti di valutazione */


	Vec3 bTmp[2];

	/* ciclo sulle equazioni */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		bTmp[0] = p[iSez] - pNode[iSez]->GetXCurr();
		bTmp[1] = p[iSez] - pNode[iSez + 1]->GetXCurr();

		unsigned int iRow1 = iSez*6 + 1;
		unsigned int iRow2 = iRow1 + 6;

		for (unsigned int i = 0; i < NUMNODES; i++) {
			/* Equazione all'indietro: */
			WMA.Sub(iRow1, 6*i + 1, AzTmp[iSez][i].GetMat11());
			WMA.Sub(iRow1, 6*i + 4, AzTmp[iSez][i].GetMat12());

			WMA.Sub(iRow1 + 3, 6*i + 1,
				AzTmp[iSez][i].GetMat21()
				- Mat3x3(MatCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
				+ bTmp[0].Cross(AzTmp[iSez][i].GetMat11()));
			WMA.Sub(iRow1 + 3, 6*i + 4,
				 AzTmp[iSez][i].GetMat22()
				- Mat3x3(MatCrossCross, Az[iSez].GetVec1()*(-dCoef*dN3[iSez][i]), fTmp[i])
				+ bTmp[0].Cross(AzTmp[iSez][i].GetMat12()));

			/* Equazione in avanti: */
			WMA.Add(iRow2, 6*i + 1, AzTmp[iSez][i].GetMat11());
			WMA.Add(iRow2, 6*i + 4, AzTmp[iSez][i].GetMat12());

			WMA.Add(iRow2 + 3, 6*i + 1,
				AzTmp[iSez][i].GetMat21()
				- Mat3x3(MatCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
				+ bTmp[1].Cross(AzTmp[iSez][i].GetMat11()));
			WMA.Add(iRow2 + 3, 6*i + 4,
				AzTmp[iSez][i].GetMat22()
				+ Mat3x3(MatCrossCross, Az[iSez].GetVec1()*(dCoef*dN3[iSez][i]), fTmp[i])
				+ bTmp[1].Cross(AzTmp[iSez][i].GetMat12()));


			/* Equazione viscosa all'indietro: */
			WMB.Sub(iRow1, 6*i + 1, AzPrimeTmp[iSez][i].GetMat11());
			WMB.Sub(iRow1, 6*i + 4, AzPrimeTmp[iSez][i].GetMat12());

			WMB.Sub(iRow1 + 3, 6*i + 1,
				AzPrimeTmp[iSez][i].GetMat21()
				+ bTmp[0].Cross(AzPrimeTmp[iSez][i].GetMat11()));
			WMB.Sub(iRow1 + 3, 6*i + 4,
				AzPrimeTmp[iSez][i].GetMat22()
				+ bTmp[0].Cross(AzPrimeTmp[iSez][i].GetMat12()));

			/* Equazione viscosa in avanti: */
			WMB.Add(iRow2, 6*i + 1, AzPrimeTmp[iSez][i].GetMat11());
			WMB.Add(iRow2, 6*i + 4, AzPrimeTmp[iSez][i].GetMat12());

			WMB.Add(iRow2 + 3, 6*i + 1,
				AzPrimeTmp[iSez][i].GetMat21()
				+ bTmp[1].Cross(AzPrimeTmp[iSez][i].GetMat11()));
			WMB.Add(iRow2 + 3, 6*i + 4,
				AzPrimeTmp[iSez][i].GetMat22()
				+ bTmp[1].Cross(AzPrimeTmp[iSez][i].GetMat12()));
		}

		/* correzione alle equazioni */
		Mat3x3 FTmp(MatCross, Az[iSez].GetVec1()*dCoef);
		WMA.Sub(iRow1 + 3, 6*iSez + 1, FTmp);
		WMA.Add(iRow2 + 3, 6*iSez + 7, FTmp);

	} /* end ciclo sui punti di valutazione */
};


/* Assembla il residuo */
void
ViscoElasticBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticBeam::AssStiffnessVec");

	/* Riceve il vettore gia' dimensionato e con gli indici a posto
	 * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */

	/* Per la trattazione teorica, il riferimento e' il file ul-travi.tex
	 * (ora e' superato) */

	if (bFirstRes) {
		bFirstRes = false; /* AfterPredict ha gia' calcolato tutto */
	} else {
		Vec3 gNod[NUMNODES];
		Vec3 xTmp[NUMNODES];

		Vec3 gPrimeNod[NUMNODES];
		Vec3 xPrimeTmp[NUMNODES];

		for (unsigned int i = 0; i < NUMNODES; i++) {
			gNod[i] = pNode[i]->GetgCurr();
			Vec3 fTmp = pNode[i]->GetRCurr()*f[i];
			xTmp[i] = pNode[i]->GetXCurr() + fTmp;
			gPrimeNod[i] = pNode[i]->GetgPCurr();
			xPrimeTmp[i] = pNode[i]->GetVCurr()+pNode[i]->GetWCurr().Cross(fTmp);
		}

		Mat3x3 RDelta[NUMSEZ];
		Vec3 gGrad[NUMSEZ];
		Vec3 gPrimeGrad[NUMSEZ];

		/* Aggiorna le grandezze della trave nei punti di valutazione */
		for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

			/* Posizione */
			p[iSez] = InterpState(xTmp[NODE1],
				xTmp[NODE2],
				xTmp[NODE3], Beam::Section(iSez));

			/* Matrici di rotazione */
			g[iSez] = InterpState(gNod[NODE1],
				gNod[NODE2],
				gNod[NODE3], Beam::Section(iSez));
			RDelta[iSez] = Mat3x3(CGR_Rot::MatR, g[iSez]);
			R[iSez] = RDelta[iSez]*RRef[iSez];

			/* Velocita' angolare della sezione */
			gPrime[iSez] = InterpState(gPrimeNod[NODE1],
				gPrimeNod[NODE2],
				gPrimeNod[NODE3], Beam::Section(iSez));
			Omega[iSez] = Mat3x3(CGR_Rot::MatG, g[iSez])*gPrime[iSez]
				+ RDelta[iSez]*OmegaRef[iSez];

			/* Derivate della posizione */
			L[iSez] = InterpDeriv(xTmp[NODE1],
				xTmp[NODE2],
				xTmp[NODE3], Beam::Section(iSez));

			/* Derivate della velocita' */
			LPrime[iSez] = InterpDeriv(xPrimeTmp[NODE1],
				xPrimeTmp[NODE2],
				xPrimeTmp[NODE3], Beam::Section(iSez));

			/* Derivate dei parametri di rotazione */
			gGrad[iSez] = InterpDeriv(gNod[NODE1],
				gNod[NODE2],
				gNod[NODE3], Beam::Section(iSez));

			/* Derivate delle derivate spaziali dei parametri di rotazione */
			gPrimeGrad[iSez] = InterpDeriv(gPrimeNod[NODE1],
				gPrimeNod[NODE2],
				gPrimeNod[NODE3], Beam::Section(iSez));

			/* Calcola le deformazioni nel sistema locale nei punti di valutazione */
			DefLoc[iSez] = Vec6(R[iSez].MulTV(L[iSez]) - L0[iSez],
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]) + DefLocRef[iSez].GetVec2());

			/* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
			DefPrimeLoc[iSez] = Vec6(R[iSez].MulTV(LPrime[iSez] + L[iSez].Cross(Omega[iSez])),
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gPrimeGrad[iSez]
				+ (Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]).Cross(Omega[iSez]))
				+ DefPrimeLocRef[iSez].GetVec2());

			/* Calcola le azioni interne */
			pD[iSez]->Update(DefLoc[iSez], DefPrimeLoc[iSez]);
			AzLoc[iSez] = pD[iSez]->GetF();

			/* corregge le azioni interne locali (piezo, ecc) */
			AddInternalForces(AzLoc[iSez], iSez);

			/* Porta le azioni interne nel sistema globale */
			Az[iSez] = MultRV(AzLoc[iSez], R[iSez]);
		}
	}

	WorkVec.Add(1, Az[S_I].GetVec1());
	WorkVec.Add(4, (p[S_I] - pNode[NODE1]->GetXCurr()).Cross(Az[S_I].GetVec1())
		+ Az[S_I].GetVec2());
	WorkVec.Add(7, Az[SII].GetVec1() - Az[S_I].GetVec1());
	WorkVec.Add(10, Az[SII].GetVec2() - Az[S_I].GetVec2()
		+ (p[SII] - pNode[NODE2]->GetXCurr()).Cross(Az[SII].GetVec1())
		- (p[S_I] - pNode[NODE2]->GetXCurr()).Cross(Az[S_I].GetVec1()));
	WorkVec.Sub(13, Az[SII].GetVec1());
	WorkVec.Add(16, Az[SII].GetVec1().Cross(p[SII] - pNode[NODE3]->GetXCurr())
		- Az[SII].GetVec2());
}


/* Settings iniziali, prima della prima soluzione */
void
ViscoElasticBeam::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	Beam::SetValue(pDM, X, XP, ph);

	/* Aggiorna le grandezze della trave nei punti di valutazione */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		OmegaRef[iSez] = Omega[iSez];
		LPrimeRef[iSez] = LPrime[iSez];
		DefPrimeLocRef[iSez] = DefPrimeLoc[iSez];

		/* Aggiorna il legame costitutivo di riferimento
		 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo) */
		ERef[iSez] = MultRMRt(pD[iSez]->GetFDEPrime(), RRef[iSez]);
	}
	ASSERT(bFirstRes == true);
}


/* Prepara i parametri di riferimento dopo la predizione */
void
ViscoElasticBeam::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */

	/* Recupera i dati dei nodi */
	Vec3   gNod[NUMNODES];
	Vec3   xTmp[NUMNODES];

	Vec3   gPrimeNod[NUMNODES];
	Vec3   xPrimeTmp[NUMNODES];

	for (unsigned int i = 0; i < NUMNODES; i++) {
		gNod[i] = pNode[i]->GetgRef();
		Vec3 fTmp = pNode[i]->GetRRef()*f[i];
		xTmp[i] = pNode[i]->GetXCurr() + fTmp;
		gPrimeNod[i] = pNode[i]->GetgPRef();
		xPrimeTmp[i] = pNode[i]->GetVCurr() + pNode[i]->GetWRef().Cross(fTmp);
	}

	Mat3x3 RDelta[NUMSEZ];
	Vec3 gGrad[NUMSEZ];
	Vec3 gPrimeGrad[NUMSEZ];

	/* Aggiorna le grandezze della trave nei punti di valutazione */
	for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
		/* Posizione */
		p[iSez] = InterpState(xTmp[NODE1],
		xTmp[NODE2],
		xTmp[NODE3], Beam::Section(iSez));

		/* Matrici di rotazione */
		g[iSez] = InterpState(gNod[NODE1],
			gNod[NODE2],
			gNod[NODE3], Beam::Section(iSez));
		RDelta[iSez] = Mat3x3(CGR_Rot::MatR, g[iSez]);
		R[iSez] = RRef[iSez] = RDelta[iSez]*RPrev[iSez];

		/* Velocita' angolare della sezione */
		gPrime[iSez] = InterpState(gPrimeNod[NODE1],
			gPrimeNod[NODE2],
			gPrimeNod[NODE3], Beam::Section(iSez));
		Omega[iSez] = OmegaRef[iSez] = Mat3x3(CGR_Rot::MatG, g[iSez])*gPrime[iSez];

		/* Derivate della posizione */
		L[iSez] = LRef[iSez] = InterpDeriv(xTmp[NODE1],
			xTmp[NODE2],
			xTmp[NODE3], Beam::Section(iSez));

		/* Derivate della velocita' */
		LPrime[iSez] = LPrimeRef[iSez]
			= InterpDeriv(xPrimeTmp[NODE1],
				xPrimeTmp[NODE2],
				xPrimeTmp[NODE3], Beam::Section(iSez));

		/* Derivate dei parametri di rotazione */
		gGrad[iSez] = InterpDeriv(gNod[NODE1],
			gNod[NODE2],
			gNod[NODE3], Beam::Section(iSez));

		/* Derivate delle derivate spaziali dei parametri di rotazione */
		gPrimeGrad[iSez] = InterpDeriv(gPrimeNod[NODE1],
			gPrimeNod[NODE2],
			gPrimeNod[NODE3], Beam::Section(iSez));

		/* Calcola le deformazioni nel sistema locale nei punti di valutazione */
		DefLoc[iSez] = DefLocRef[iSez]
			= Vec6(R[iSez].MulTV(L[iSez]) - L0[iSez],
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]) + DefLocPrev[iSez].GetVec2());

		/* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
		DefPrimeLoc[iSez] = DefPrimeLocRef[iSez]
			= Vec6(R[iSez].MulTV(LPrime[iSez] + L[iSez].Cross(Omega[iSez])),
				R[iSez].MulTV(Mat3x3(CGR_Rot::MatG, g[iSez])*gPrimeGrad[iSez]
					+ (Mat3x3(CGR_Rot::MatG, g[iSez])*gGrad[iSez]).Cross(Omega[iSez])));

		/* Calcola le azioni interne */
		pD[iSez]->Update(DefLoc[iSez], DefPrimeLoc[iSez]);
		AzLoc[iSez] = pD[iSez]->GetF();

		/* corregge le azioni interne locali (piezo, ecc) */
		AddInternalForces(AzLoc[iSez], iSez);

		/* Porta le azioni interne nel sistema globale */
		Az[iSez] = AzRef[iSez] = MultRV(AzLoc[iSez], R[iSez]);

		/* Aggiorna il legame costitutivo */
		DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), R[iSez]);
		ERef[iSez] = MultRMRt(pD[iSez]->GetFDEPrime(), R[iSez]);
	}

	bFirstRes = true;
}

doublereal
ViscoElasticBeam::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= iGetNumPrivData());

	int sez = (i - 1)/iNumPrivData;
	int idx = (i - 1)%iNumPrivData + 1;

	switch (idx) {
	case 22:
	case 23:
	case 24:
	case 25:
	case 26:
	case 27:
		return DefPrimeLoc[sez].dGet(idx - 21);

	default:
		return Beam::dGetPrivData(i);
	}
}

/* ViscoElasticBeam - end */

void
ReadBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od)
{
	uFlags = Beam::OUTPUT_NONE;

	while (HP.IsArg()) {
		unsigned uFlag;

		if (HP.IsKeyWord("position")) {
			uFlag = Beam::OUTPUT_EP_X;

		} else if (HP.IsKeyWord("orientation")) {
			uFlag = Beam::OUTPUT_EP_R;

		} else if (HP.IsKeyWord("configuration")) {
			uFlag = Beam::OUTPUT_EP_CONFIGURATION;

		} else if (HP.IsKeyWord("force")) {
			uFlag = Beam::OUTPUT_EP_F;

		} else if (HP.IsKeyWord("moment")) {
			uFlag = Beam::OUTPUT_EP_M;

		} else if (HP.IsKeyWord("forces")) {
			uFlag = Beam::OUTPUT_EP_FORCES;

		} else if (HP.IsKeyWord("linear" "strain")) {
			uFlag = Beam::OUTPUT_EP_NU;

		} else if (HP.IsKeyWord("angular" "strain")) {
			uFlag = Beam::OUTPUT_EP_K;

		} else if (HP.IsKeyWord("strains")) {
			uFlag = Beam::OUTPUT_EP_STRAINS;

		} else if (HP.IsKeyWord("linear" "strain" "rate")) {
			uFlag = Beam::OUTPUT_EP_NUP;

		} else if (HP.IsKeyWord("angular" "strain" "rate")) {
			uFlag = Beam::OUTPUT_EP_KP;

		} else if (HP.IsKeyWord("strain rates")) {
			uFlag = Beam::OUTPUT_EP_STRAIN_RATES;

		} else if (HP.IsKeyWord("all")) {
			uFlag = Beam::OUTPUT_EP_ALL;

		} else if (HP.IsKeyWord("none")) {
			uFlag = Beam::OUTPUT_NONE;

		} else {
			break;
		}

		if (uFlags & uFlag) {
			silent_cerr("Beam(" << uLabel << "): "
				"duplicate custom output "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (uFlag & Beam::OUTPUT_EP_STRAIN_RATES) {
			if (BT == Beam::ELASTIC) {
				silent_cerr("Beam(" << uLabel << "): "
					"strain rates only allowed for viscoelastic beams (ignored) "
					"at line " << HP.GetLineData()
					<< std::endl);
				uFlag &= ~Beam::OUTPUT_EP_STRAIN_RATES;
			}
		}

		if (uFlag & Beam::OUTPUT_EP_R) {
			od = ReadOptionalOrientationDescription(pDM, HP);
		}

		if (uFlag == Beam::OUTPUT_NONE) {
			// reset all
			uFlags = Beam::OUTPUT_NONE;
		}

		uFlags |= uFlag;
	}
}


void
ReadOptionalBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od)
{
	pDM->GetOutput(Elem::BEAM, uFlags, od);
	if (HP.IsKeyWord("custom" "output")) {
		ReadBeamCustomOutput(pDM, HP, uLabel, BT, uFlags, od);
	}
}


/* Legge una trave */

Elem *
ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadBeam");

	/* Per ogni nodo: */

	/* Nodo 1 */
	const StructNode* pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	const Mat3x3& R1(pNode1->GetRCurr());
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
	Mat3x3 Rn1(Eye3);
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		Rn1 = HP.GetRotRel(ReferenceFrame(pNode1));
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): "
		<< f1 << std::endl
		<< "(global frame): "
		<< pNode1->GetXCurr()+pNode1->GetRCurr()*f1 << std::endl);

	/* Nodo 2 */
	const StructNode* pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	Mat3x3 R2(pNode2->GetRCurr());
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
	Mat3x3 Rn2(Eye3);
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		Rn2 = HP.GetRotRel(ReferenceFrame(pNode2));
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "node 2 offset (node reference frame): "
		<< f2 << std::endl
		<< "(global frame): "
		<< pNode2->GetXCurr()+pNode2->GetRCurr()*f2 << std::endl);

	/* Nodo 3 */
	const StructNode* pNode3 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	Mat3x3 R3(pNode3->GetRCurr());
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f3(HP.GetPosRel(ReferenceFrame(pNode3)));
	Mat3x3 Rn3(Eye3);
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		Rn3 = HP.GetRotRel(ReferenceFrame(pNode3));
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "node 3 offset (node reference frame): "
		<< f3 << std::endl
		<< "(global frame): "
		<< pNode3->GetXCurr()+pNode3->GetRCurr()*f3 << std::endl);

	/*   Per ogni punto: */

	/* Punto I */

	/*     Matrice R */
	Mat3x3 R_I;
	bool b_I(false);
	if (HP.IsKeyWord("from" "nodes") || HP.IsKeyWord("node")) {
		b_I = true;
	} else {
		R_I = HP.GetRotAbs(::AbsRefFrame);
	}

	/*     Legame costitutivo */
	ConstLawType::Type CLType_I = ConstLawType::UNKNOWN;
	ConstitutiveLaw6D* pD_I = HP.GetConstLaw6D(CLType_I);

	if (pD_I->iGetNumDof() != 0) {
		silent_cerr("line " << HP.GetLineData()
			<< ": Beam(" << uLabel << ") "
			"does not support dynamic constitutive laws yet"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef DEBUG
	Mat6x6 MTmp(pD_I->GetFDE());
	Mat3x3 D11(MTmp.GetMat11());
	Mat3x3 D12(MTmp.GetMat12());
	Mat3x3 D21(MTmp.GetMat21());
	Mat3x3 D22(MTmp.GetMat22());

	DEBUGLCOUT(MYDEBUG_INPUT,
		"First point matrix D11: " << std::endl << D11 << std::endl
		<< "First point matrix D12: " << std::endl << D12 << std::endl
		<< "First point matrix D21: " << std::endl << D21 << std::endl
		<< "First point matrix D22: " << std::endl << D22 << std::endl);
#endif


	/* Punto II */

	/* Matrice R */
	Mat3x3 RII;
	bool bII(false);
	if (HP.IsKeyWord("same")) {
		if (b_I) {
			bII = true;
		} else {
			RII = R_I;
		}
	} else {
		if (HP.IsKeyWord("from" "nodes") || HP.IsKeyWord("node")) {
			bII = true;
		} else {
			RII = HP.GetRotAbs(::AbsRefFrame);
		}
	}

	/* Legame costitutivo */
	ConstLawType::Type CLTypeII = ConstLawType::UNKNOWN;
	ConstitutiveLaw6D* pDII = NULL;

	if (HP.IsKeyWord("same")) {
		pDII = pD_I->pCopy();
		CLTypeII = CLType_I;

	} else {
		pDII = HP.GetConstLaw6D(CLTypeII);

		if (pDII->iGetNumDof() != 0) {
			silent_cerr("line " << HP.GetLineData()
				<< ": Beam(" << uLabel << ") "
				"does not support dynamic constitutive laws yet"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	Beam::Type Type;
	if (CLType_I == ConstLawType::ELASTIC && CLTypeII == ConstLawType::ELASTIC) {
		Type = Beam::ELASTIC;
	} else {
		Type = Beam::VISCOELASTIC;
	}

#ifdef DEBUG
	MTmp = pDII->GetFDE();
	D11 = MTmp.GetMat11();
	D12 = MTmp.GetMat12();
	D21 = MTmp.GetMat21();
	D22 = MTmp.GetMat22();

	DEBUGLCOUT(MYDEBUG_INPUT,
		"Second point matrix D11: " << std::endl << D11 << std::endl
		<< "Second point matrix D12: " << std::endl << D12 << std::endl
		<< "Second point matrix D21: " << std::endl << D21 << std::endl
		<< "Second point matrix D22: " << std::endl << D22 << std::endl);
#endif

	flag fPiezo(0);
	Mat3xN PiezoMat[2][2];
	integer iNumElec = 0;
	const ScalarDifferentialNode **pvElecDofs = 0;
	if (HP.IsKeyWord("piezoelectric" "actuator")) {
		fPiezo = flag(1);
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Piezoelectric actuator beam is expected" << std::endl);

		iNumElec = HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT,
			"piezo actuator " << uLabel << " has " << iNumElec
			<< " electrodes" << std::endl);
		if (iNumElec <= 0) {
			silent_cerr("PiezoElectricBeam(" << uLabel << "): "
				"illegal number of electric nodes " << iNumElec
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		SAFENEWARR(pvElecDofs, const ScalarDifferentialNode *, iNumElec);

		for (integer i = 0; i < iNumElec; i++) {
			pvElecDofs[i] = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
		}

		PiezoMat[0][0].Resize(iNumElec);
		PiezoMat[1][0].Resize(iNumElec);
		PiezoMat[0][1].Resize(iNumElec);
		PiezoMat[1][1].Resize(iNumElec);

		/* leggere le matrici (6xN sez. 1, 6xN sez. 2) */
		HP.GetMat6xN(PiezoMat[0][0], PiezoMat[1][0], iNumElec);
		if (HP.IsKeyWord("same")) {
			PiezoMat[0][1].Copy(PiezoMat[0][0]);
			PiezoMat[1][1].Copy(PiezoMat[1][0]);
		} else {
			HP.GetMat6xN(PiezoMat[0][1], PiezoMat[1][1], iNumElec);
		}

#if 0
		DEBUGLCOUT(MYDEBUG_INPUT, "Piezo matrix I:" << std::endl << PiezoMat[0][0] << PiezoMat[1][0]);
		DEBUGLCOUT(MYDEBUG_INPUT, "Piezo matrix II:" << std::endl << PiezoMat[0][1] << PiezoMat[1][1]);
#endif /* 0 */
	}

	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	unsigned uFlags = Beam::OUTPUT_NONE;
	ReadOptionalBeamCustomOutput(pDM, HP, uLabel, Type, uFlags, od);

	flag fOut = pDM->fReadOutput(HP, Elem::BEAM);
	if (fOut) {
		fOut |= uFlags;
	}

	/* Se necessario, interpola i parametri di rotazione delle sezioni */
	if (b_I || bII) {
		Mat3x3 R(R2*Rn2);
		Vec3 g1(CGR_Rot::Param, R.MulTM(R1*Rn1));
		Vec3 g3(CGR_Rot::Param, R.MulTM(R3*Rn3));
		if (b_I) {
			R_I = R*Mat3x3(CGR_Rot::MatR, Beam::InterpState(g1, Zero3, g3, Beam::S_I));
		}
		if (bII) {
			RII = R*Mat3x3(CGR_Rot::MatR, Beam::InterpState(g1, Zero3, g3, Beam::SII));
		}
	}

	std::ostream& out = pDM->GetLogFile();
	out << "beam3: " << uLabel
		<< " " << pNode1->GetLabel()
		<< " ", f1.Write(out, " ")
		<< " " << pNode2->GetLabel()
		<< " ", f2.Write(out, " ")
		<< " " << pNode3->GetLabel()
		<< " ", f3.Write(out, " ")
		<< std::endl;

	Elem* pEl = NULL;

	if ((CLType_I == ConstLawType::ELASTIC)
		&& (CLTypeII == ConstLawType::ELASTIC))
	{
		if (fPiezo == 0) {
			SAFENEWWITHCONSTRUCTOR(pEl,
				Beam,
				Beam(uLabel,
					pNode1, pNode2, pNode3,
					f1, f2, f3,
					Rn1, Rn2, Rn3,
					R_I, RII,
					pD_I, pDII,
					od, fOut));
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
				PiezoActuatorBeam,
				PiezoActuatorBeam(uLabel,
					pNode1, pNode2, pNode3,
					f1, f2, f3,
					Rn1, Rn2, Rn3,
					R_I, RII,
					pD_I, pDII,
					iNumElec,
					pvElecDofs,
					PiezoMat[0][0], PiezoMat[1][0],
					PiezoMat[0][1], PiezoMat[1][1],
					od, fOut));
		}

	} else /* At least one is VISCOUS or VISCOELASTIC */ {
		if (fPiezo == 0) {
			SAFENEWWITHCONSTRUCTOR(pEl,
				ViscoElasticBeam,
				ViscoElasticBeam(uLabel,
					pNode1, pNode2, pNode3,
					f1, f2, f3,
					Rn1, Rn2, Rn3,
					R_I, RII,
					pD_I, pDII,
					od, fOut));
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
				PiezoActuatorVEBeam,
				PiezoActuatorVEBeam(uLabel,
					pNode1, pNode2, pNode3,
					f1, f2, f3,
					Rn1, Rn2, Rn3,
					R_I, RII,
					pD_I, pDII,
					iNumElec,
					pvElecDofs,
					PiezoMat[0][0], PiezoMat[1][0],
					PiezoMat[0][1], PiezoMat[1][1],
					od, fOut));
		}
	}

	/* Costruttore normale
	 * Beam(unsigned int uL,
	 *      const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
	 *	   const Vec3& X1, const Vec3& X2, const Vec3& X3,
	 *	   const Vec3& F1, const Vec3& F2, const Vec3& F3,
	 *	   const Mat3x3& r_I, const Mat3x3& rII,
	 *	   const Mat3x3& d11_I, const Mat3x3& d12_I,
	 *	   const Mat3x3& d21_I, const Mat3x3& d22_I,
	 *	   const Mat3x3& d11II, const Mat3x3& d12II,
	 *	   const Mat3x3& d21II, const Mat3x3& d22II,
	 *	   const Vec3& eps0_I, const Vec3& k0_I,
	 *	   const Vec3& eps0II, const Vec3& k0II);
	 */


	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* End of ReadBeam() */

