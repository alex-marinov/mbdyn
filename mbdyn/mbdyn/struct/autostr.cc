/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/math.h>
#include <ac/float.h>

#include <autostr.h>

/* Costruttore */
AutomaticStructElem::AutomaticStructElem(const DynamicStructNode* pN)
: Elem(pN->GetLabel(), Elem::AUTOMATICSTRUCTURAL, pN->fToBeOutput()), 
pNode((DynamicStructNode *)pN), Q(0.), G(0.), QP(0.), GP(0.),
m(0.), S(0.), J(0.)
{ 
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
	/* FIXME: we export the test because we don't want Inv() to fail
	 * or issue error messages */
	Mat3x3 Jcg = J + Mat3x3(Xcg, S);
	doublereal dDet = Jcg.dDet();
	if (fabs(dDet) > DBL_EPSILON) {
		WP = Jcg.Inv(dDet, GP - Xcg.Cross(QP)
			- pNode->GetWCurr().Cross(G));
	} else {
		WP = Zero3;
	}
	XPP = QP/m + Xcg.Cross(WP)
		- pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(Xcg));
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
AutomaticStructElem::Init(const Vec3& q, const Vec3& g, 
			  const Vec3& qp, const Vec3& gp)
{
   Q = q;
   G = g;
   QP = qp;
   GP = gp;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream& 
AutomaticStructElem::Restart(std::ostream& out) const
{
   out << "    automatic structural: " << GetLabel() << ", "
     "reference, global, ", Q.Write(out, ", ") << ", "
     "reference, global, ", G.Write(out, ", ") << ", "
     "reference, global, ", QP.Write(out, ", ") << ", "
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
   WM.ResizeReset(24, 0);
      
   /* Setta gli indici della matrice - le incognite sono ordinate come:
    *   - posizione (3)
    *   - parametri di rotazione (3)
    *   - quantita' di moto (3)
    *   - momento della quantita' di moto 
    * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
    * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
    * e' dato da iGetFirstPositionIndex()+i */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
  
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutItem(iCnt, iFirstPositionIndex+iCnt,
		  iFirstMomentumIndex+iCnt, -dCoef);
      WM.PutItem(6+iCnt, iFirstMomentumIndex+iCnt,
		  iFirstMomentumIndex+iCnt, 1.);    
   }

   WM.PutCross(13, iFirstMomentumIndex+3, iFirstMomentumIndex, 
		   pNode->GetVCurr()*dCoef);
   WM.PutCross(19, iFirstMomentumIndex+3, iFirstPositionIndex, -Q);
   
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
   WMA.ResizeReset(6, 0);
   WMB.ResizeReset(6, 0);
      
   /* Setta gli indici della matrice - le incognite sono ordinate come:
    *   - posizione (3)
    *   - parametri di rotazione (3)
    *   - quantita' di moto (3)
    *   - momento della quantita' di moto 
    * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
    * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
    * e' dato da iGetFirstPositionIndex()+i */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
  
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WMA.PutItem(iCnt, iFirstPositionIndex+iCnt,
		   iFirstMomentumIndex+iCnt, -1.);
      WMB.PutItem(iCnt, iFirstMomentumIndex+iCnt,
		   iFirstMomentumIndex+iCnt, 1.);    
   }   
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
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
   }   
   
   /* Collects data */
   Q = Vec3(XCurr, iFirstMomentumIndex+1);
   G = Vec3(XCurr, iFirstMomentumIndex+4);
   QP = Vec3(XPrimeCurr, iFirstMomentumIndex+1);
   GP = Vec3(XPrimeCurr, iFirstMomentumIndex+4);
   
   /*
    * Momentum and momenta moment (about node):
    *
    * Q = m V + W /\ S
    *
    * G = S /\ V + J W
    *
    * Qp = F
    *
    * Gp + V /\ Q = M
    */
   WorkVec.Add(1, Q);
   WorkVec.Add(4, G);
   WorkVec.Sub(7, QP);
   WorkVec.Sub(10, GP + pNode->GetVCurr().Cross(Q));

   m = 0.;
   S = Zero3;
   J = Zero3x3;

   return WorkVec;
}


void 
AutomaticStructElem::Output(OutputHandler& OH) const
{
   ASSERT(pNode != NULL);
   if(pNode->fToBeOutput()) {
      OH.Inertia() << std::setw(8) << GetLabel() << " " 
	<< Q << " " << G << " " << QP << " " << GP << std::endl;
   }
}


/* Setta i valori iniziali delle variabili (e fa altre cose) 
 * prima di iniziare l'integrazione */
void 
AutomaticStructElem::SetValue(VectorHandler& /* X */ , VectorHandler& XP) const
{
   integer iIndex = pNode->iGetFirstMomentumIndex();
   
   XP.Put(iIndex + 1, QP);
   XP.Put(iIndex + 3 + 1, GP);
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
	if (strncmp(s, "beta", sizeof("beta") - 1) == 0) {
		s += sizeof("beta") - 1;
	} else if (strncmp(s, "gamma", sizeof("gamma") - 1) == 0) {
		s += sizeof("gamma") - 1;
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
		return QP(i);

	} else {
		if (type) {
			return G(i);
		}
		return Q(i);
	}
}

