/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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

#include <autostr.h>

/* Costruttore */
AutomaticStructElem::AutomaticStructElem(const DynamicStructNode* pN)
: Elem(pN->GetLabel(), Elem::AUTOMATICSTRUCTURAL, pN->fToBeOutput()), 
pNode(pN), Q(0.), G(0.), QP(0.), GP(0.)
{ 
   NO_OP;
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
   WM.ResizeInit(24, 0, 0.);
      
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
      WM.fPutItem(iCnt, iFirstPositionIndex+iCnt,
		  iFirstMomentumIndex+iCnt, -dCoef);
      WM.fPutItem(6+iCnt, iFirstMomentumIndex+iCnt,
		  iFirstMomentumIndex+iCnt, 1.);    
   }

   WM.fPutCross(13, iFirstMomentumIndex+3, iFirstMomentumIndex, 
		   pNode->GetVCurr()*dCoef);
   WM.fPutCross(19, iFirstMomentumIndex+3, iFirstPositionIndex, -Q);
   
   return WorkMat;
}


/* assemblaggio autoval */
void 
AutomaticStructElem::AssEig(VariableSubMatrixHandler& WorkMatA,
			    VariableSubMatrixHandler& WorkMatB,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("AutomaticStructElem::AssEig");

   /* Casting di WorkMat */
   SparseSubMatrixHandler& WMA = WorkMatA.SetSparse();
   SparseSubMatrixHandler& WMB = WorkMatB.SetSparse();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WMA.ResizeInit(6, 0, 0.);
   WMB.ResizeInit(6, 0, 0.);
      
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
      WMA.fPutItem(iCnt, iFirstPositionIndex+iCnt,
		   iFirstMomentumIndex+iCnt, -1.);
      WMB.fPutItem(iCnt, iFirstMomentumIndex+iCnt,
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

   WorkVec.Resize(12);
   WorkVec.Reset(0.);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 12; iCnt++) {      
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
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
   WorkVec.Sub(10, GP+pNode->GetVCurr().Cross(Q));
   
   return WorkVec;
}


void 
AutomaticStructElem::Output(OutputHandler& OH) const
{
   ASSERT(pNode != NULL);
   if(pNode->fToBeOutput()) {
#ifdef DEBUG   
      OH.Output() << "Automatic structural element " << GetLabel() << std::endl
	<< "Momentum: " << std::endl << Q <<std::endl << G << std::endl
	<< "Momentum derivative: " << std::endl << QP << std::endl << GP << std::endl;
#endif   
   
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
   
   XP.Put(iIndex+1, QP);
   XP.Put(iIndex+4, GP);
}

