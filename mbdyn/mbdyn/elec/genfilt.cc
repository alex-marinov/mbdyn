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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_ELECTRIC_NODES

#include <genfilt.h>
#ifdef USE_MESCHACH
#include <mschwrap.h>
#undef catch
#endif /* USE_MESCHACH */

#if 0
/* GenelFilter - begin */

GenelFilter::GenelFilter(unsigned int uLabel, const DofOwner* pDO,
			 const ScalarDof& y, const ScalarDof& u,
			 unsigned int na, unsigned int nb,
			 doublereal* p, doublereal* tau,
			 flag fOutput)
: Elem(uLabel, Elem::GENEL, fOutput),
Genel(uLabel, Genel::SCALARFILTER, pDO, fOutput),
SD_y(y), SD_u(u),
Na(na), Nb(nb),
pdP(p), pdTau(tau) 
{
   ASSERT(na > 0 ? pdP != NULL : 1);
   ASSERT(nb > 0 ? pdTau != NULL : 1);
   ASSERT(SD_y.iOrder == 0);
   ASSERT(SD_u.iOrder == 0);
   iNumDofs = max(integer(na-1), 0)+max(integer(nb-1), 0);
   DEBUGCOUT("GenelFilter " << uLabel << ", NumDofs: " << iNumDofs << endl);
}


GenelFilter::~GenelFilter(void) 
{
   if (pdTau != NULL) {
      SAFEDELETEARR(pdTau);
   }
   if (pdP != NULL) {
      SAFEDELETEARR(pdP);
   }
}
 
   
unsigned int GenelFilter::iGetNumDof(void) const
{
   return iNumDofs;
}

   
/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order GenelFilter::SetDof(unsigned int i) const
{
   ASSERT(i < iNumDofs);
   return DofOrder::DIFFERENTIAL;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& GenelFilter::Restart(ostream& out) const
{
   return out << "GenelFilter: not implemented yet!" << endl;
}


/* Dimensioni del workspace */
void GenelFilter::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = max(integer(Na)-1, 0)+max(integer(Nb)-1, 0)+1;
   *piNumCols = max(integer(Na)-1, 0)+max(integer(Nb)-1, 0)+2;
}


void GenelFilter::Output(OutputHandler& /* OH */ ) const
{
   NO_OP;
}
 

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
GenelFilter::AssJac(VariableSubMatrixHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& /* XCurr */ ,
		    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering GenelFilter::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
   integer iColIndex_y = SD_y.pNode->iGetFirstColIndex()+1;
   integer iColIndex_u = SD_u.pNode->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   unsigned int row_u = 1;
   unsigned int col_u = 2;
   WM.fPutRowIndex(1, iRowIndex_y);
   WM.fPutColIndex(1, iColIndex_y);
   WM.fPutCoef(1, 1, dCoef); /* sarebbe *(pdP), ma e' == 1. per def. */
   if (Na >= 1) {
      if (Na > 1) {
	 row_u = Na;
	 col_u = Na+1;
	 for (unsigned int iCnt = 1; iCnt <= Na-1; iCnt++) {
	    WM.fPutRowIndex(iCnt+1, iFirstIndex+iCnt);
	    WM.fPutColIndex(iCnt+1, iFirstIndex+iCnt);
	    WM.fPutCoef(iCnt+1, iCnt, 1.);
	    WM.fPutCoef(iCnt+1, iCnt+1, -dCoef);
	    WM.fPutCoef(1, iCnt+1, pdP[iCnt]*dCoef);
	 }
	 iFirstIndex += Na-1;
      }
      WM.fIncCoef(1, Na, pdP[Na]);
   }	          
   
   WM.fPutColIndex(col_u, iColIndex_u);
   WM.fPutCoef(1, col_u, -pdTau[0]*dCoef);
   if (Nb >= 1) {
      if (Nb > 1) {
	 for (unsigned int iCnt = 1; iCnt <= Nb-1; iCnt++) {
	    WM.fPutRowIndex(row_u+iCnt, iFirstIndex+iCnt);
	    WM.fPutColIndex(col_u+iCnt, iFirstIndex+iCnt);
	    WM.fPutCoef(row_u+iCnt, col_u-1+iCnt, 1.);
	    WM.fPutCoef(row_u+iCnt, col_u+iCnt, -dCoef);
	    WM.fPutCoef(1, col_u+iCnt, -pdTau[iCnt]*dCoef);
	 }	      
      }
      WM.fIncCoef(1, col_u+Nb-1, -pdTau[Nb]);
   }
   
   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
GenelFilter::AssRes(SubVectorHandler& WorkVec,
		    doublereal /* dCoef */ ,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering GenelFilter::AssRes()" << endl);
   
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   integer row_u = 1;
   
   doublereal y = SD_y.pNode->dGetX();
   doublereal d = -y;
   if (Na >= 1) {
      doublereal yp = SD_y.pNode->dGetXPrime();
      if (Na > 1) {
	 row_u = Na;
	 for (unsigned int iCnt = 1; iCnt <= Na-1; iCnt++) {
	    WorkVec.fPutRowIndex(iCnt+1, iFirstIndex+iCnt);
	    y = XCurr.dGetCoef(iFirstIndex+iCnt);	       
	    d -= pdP[iCnt]*y;
	    WorkVec.fPutCoef(iCnt+1, y-yp);	       
	    yp = XPrimeCurr.dGetCoef(iFirstIndex+iCnt);
	 }
	 iFirstIndex += Na-1;
      } 
      d -= pdP[Na]*yp;	 
   }
   
   doublereal u = SD_u.pNode->dGetX();
   d += pdTau[0]*u;
   if (Nb >= 1) {
      doublereal up = SD_u.pNode->dGetXPrime();
      if (Nb > 1) {
	 for (unsigned int iCnt = 1; iCnt <= Nb-1; iCnt++) {
	    WorkVec.fPutRowIndex(row_u+iCnt, iFirstIndex+iCnt);
	    u = XCurr.dGetCoef(iFirstIndex+iCnt);
	    d += pdTau[iCnt]*u;
	    WorkVec.fPutCoef(row_u+iCnt, u-up);	       
	    up = XPrimeCurr.dGetCoef(iFirstIndex+iCnt);
	 }
      }
      d += pdTau[Nb]*up;
   }      
   
   WorkVec.fPutItem(1, iRowIndex_y, d);
   
   return WorkVec;
}

/* GenelFilter - end */
#endif

/* GenelFilterEq - begin */

GenelFilterEq::GenelFilterEq(unsigned int uLabel, const DofOwner* pDO,
			     const ScalarDof& y, const ScalarDof& u,
			     unsigned int na, unsigned int nb,
			     doublereal* pa, doublereal* pb,
			     flag fSt, flag fOutput)
: Elem(uLabel, Elem::GENEL, fOutput),
Genel(uLabel, Genel::SCALARFILTER, pDO, fOutput),
SD_y(y), SD_u(u),
Na(na), Nb(nb),
pdA(pa), pdB(pb),
pdAlpha(NULL), pdBeta(NULL),
fSteady(fSt)
{
   ASSERT(Na > 0 ? pdA != NULL : 1);
   ASSERT(Nb > 0 ? pdB != NULL : 1);
   ASSERT(Na >= Nb);
   ASSERT(SD_y.iOrder == 0);
   
   DEBUGCOUT("GenelFilterEq " << uLabel << ", NumDofs: " << Na << endl);
   
   SAFENEWARR(pdAlpha, doublereal, Na+Nb+1);
      
   for (unsigned long i = 0; i < Na; i++) {
      pdAlpha[i] = -pdA[i]/pdA[Na];
   }
   
   pdBeta = pdAlpha+Na;
   if (Nb == Na) {
      pdBeta[Nb] = pdB[Nb]/pdA[Na];
      for (unsigned long i = 0; i < Nb; i++) {
	 pdBeta[i] = pdB[i]/pdA[Na]-pdA[i]/pdA[Na]*pdBeta[Nb];
      }
   } else {
      for (unsigned long i = 0; i <= Nb; i++) {
	 pdBeta[i] = pdB[i]/pdA[Na];
      }
   }
}


GenelFilterEq::~GenelFilterEq(void) 
{
   if (pdAlpha != NULL) {
      SAFEDELETEARR(pdAlpha);
   }
   if (pdB != NULL) {
      SAFEDELETEARR(pdB);
   }
   if (pdA != NULL) {
      SAFEDELETEARR(pdA);
   }
}
 
   
unsigned int GenelFilterEq::iGetNumDof(void) const
{
   return Na;
}

   
/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order GenelFilterEq::SetDof(unsigned int i) const
{
   ASSERT(i < Na);
   return DofOrder::DIFFERENTIAL;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& GenelFilterEq::Restart(ostream& out) const
{
   return out << "GenelFilterEq: not implemented yet!" << endl;
}


/* Dimensioni del workspace */
void GenelFilterEq::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = Na+1;
   *piNumCols = Na+1;
}


void GenelFilterEq::Output(OutputHandler& /* OH */ ) const
{
   NO_OP;
}
 

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
GenelFilterEq::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering GenelFilterEq::AssJac()" << endl);
   
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.Resize(3*Na, 0);
   /* 
    * Na coefficienti per la diagonale,
    * Na-1 per la subdiagonale,
    * Na-1 per la colonna degli alpha
    * 2 per la riga di y = x_n
    */
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
   integer iColIndex_y = SD_y.pNode->iGetFirstColIndex()+1;  
   integer iFirstIndex = iGetFirstIndex();

   for (unsigned long i = 1; i < Na; i++) {
      /* diagonale con 1. */
      WM.fPutItem(i, iFirstIndex+i, iFirstIndex+i, 1.);
      /* Subdiagonale */
      WM.fPutItem(Na+i, iFirstIndex+1+i, iFirstIndex+i, -dCoef);
      /* colonna con Alpha */
      WM.fPutItem(2*Na+i, iFirstIndex+i, iFirstIndex+Na, -dCoef*pdAlpha[i-1]);
   }
   /* Ultimo termine della diagonale + colonna */
   WM.fPutItem(Na, iFirstIndex+Na, iFirstIndex+Na, 1.-dCoef*pdAlpha[Na-1]);
   
   /* Ultimo termine subdiagonale */
   WM.fPutItem(2*Na, iRowIndex_y, iFirstIndex+Na, -dCoef);
   
   /* coeff di y */
   WM.fPutItem(3*Na, iRowIndex_y, iColIndex_y, dCoef);
   
   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
GenelFilterEq::AssRes(SubVectorHandler& WorkVec,
		      doublereal /* dCoef */ ,
		      const VectorHandler& XCurr,
		      const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering GenelFilterEq::AssRes()" << endl);
   
   WorkVec.Resize(Na+1); 
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal y = SD_y.pNode->dGetX();
   doublereal u = SD_u.pNode->dGetX();
   
   doublereal dXn = XCurr.dGetCoef(iFirstIndex+Na);
   doublereal dX = XCurr.dGetCoef(iFirstIndex+1);
   doublereal dXP = XPrimeCurr.dGetCoef(iFirstIndex+1);

   if (Na > 0) {
      WorkVec.fPutItem(1, iFirstIndex+1, pdAlpha[0]*dXn-dXP+pdBeta[0]*u);
   }
  
   if (Na == Nb) {
      for (unsigned long i = 2; i <= Na; i++) {
	 dXP = XPrimeCurr.dGetCoef(iFirstIndex+i);
	 WorkVec.fPutItem(i, iFirstIndex+i, 
			  pdAlpha[i-1]*dXn-dXP+dX+pdBeta[i-1]*u);
	 dX = XCurr.dGetCoef(iFirstIndex+i);
      }
      WorkVec.fPutItem(Na+1, iRowIndex_y, dXn-y+pdBeta[Nb]*u);
   } else {
      for (unsigned long i = 2; i <= Na; i++) {
	 dXP = XPrimeCurr.dGetCoef(iFirstIndex+i);
	 WorkVec.fPutItem(i, iFirstIndex+i, 
			  pdAlpha[i-1]*dXn-dXP+dX);
	 dX = XCurr.dGetCoef(iFirstIndex+i);
      }
      for (unsigned long i = 2; i <= Nb; i++) {       
	 WorkVec.fIncCoef(i, pdBeta[i-1]*u);
      }
      WorkVec.fPutItem(Na+1, iRowIndex_y, dXn-y);
   }

   return WorkVec;
}


/* Setta i valori iniziali delle variabili (e fa altre cose) 
 * prima di iniziare l'integrazione */
void 
GenelFilterEq::SetValue(VectorHandler& X, VectorHandler& XP) const
{
#ifdef USE_MESCHACH
   if (fSteady) {
      DEBUGCOUT("Finding initial conditions for scalar filter " << GetLabel() << endl);
      
      if (Na == 0) {
	 integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
	 doublereal u = SD_u.pNode->dGetX();
	 X.fPutCoef(iRowIndex_y, pdBeta[Nb]*u);
# ifdef DEBUG
	    cout << "Y = " << pdBeta[Nb]*u << endl;
# endif	 
	 return;
      } 
            
      MeschachSparseLUSolutionManager sm(Na);
      
      /* preparo matrice */
      DEBUGCOUT("Preparing matrix ..." << endl);
      sm.MatrInit();
      
      MatrixHandler* mh = sm.pMatHdl();
      for (unsigned long i = 1; i < Na; i++) {
	 /* Subdiagonale */
	 mh->fPutCoef(i+1, i, -1.);
	 /* colonna con Alpha */
	 mh->fPutCoef(i, Na, -pdAlpha[i-1]);
      }
      /* Ultimo termine della diagonale + colonna */
      mh->fPutCoef(Na, Na, -pdAlpha[Na-1]);
      
      /* preparo termine noto */
      DEBUGCOUT("Preparing rhs ..." << endl);
      VectorHandler* rh = sm.pResHdl();
      doublereal u = SD_u.pNode->dGetX();

      for (unsigned long i = 1; i <= Nb; i++) {
	 rh->fPutCoef(i, pdBeta[i-1]*u);
      }
      
      /* risolvo */
# ifdef USE_EXCEPTIONS
      try {
# endif 
	 DEBUGCOUT("Solving ..." << endl);
	 sm.Solve();
# ifdef USE_EXCEPTIONS
      } 
      catch (...) {
	 cerr << "Matrix might be singular; skipping state initialization" << endl;
	 return;
      }
# endif 
      
      /* setto coefficienti */
      DEBUGCOUT("Solution:" << endl);
      integer iFirstIndex = iGetFirstIndex();
      VectorHandler* sh = sm.pSolHdl();
      for (unsigned long i = 1; i <= Na; i++) {
# ifdef DEBUG
	 cout << "X[" << i << "] = " << sh->dGetCoef(i) << endl;
# endif
	 X.fPutCoef(iFirstIndex+i, sh->dGetCoef(i));
      }
      
      integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
      doublereal dXn = sh->dGetCoef(Na);
      if (Na == Nb) {
	 X.fPutCoef(iRowIndex_y, dXn+pdBeta[Nb]*u);
# ifdef DEBUG
	 cout << "Y = " << dXn+pdBeta[Nb]*u << endl;
# endif
      } else {
# ifdef DEBUG
	 cout << "Y = " << dXn << endl;
# endif	 
	 X.fPutCoef(iRowIndex_y, dXn);
      }
   }

#endif /* USE_MESCHACH */
}
/* GenelFilterEq - end */


/* GenelStateSpaceSISO - begin */

GenelStateSpaceSISO::GenelStateSpaceSISO(unsigned int uLabel,
					 const DofOwner* pDO, 
					 const ScalarDof& y, 
					 const ScalarDof& u,
					 unsigned int Order,
					 doublereal* pA,
					 doublereal* pB,
					 doublereal* pC, 
					 doublereal D,
					 flag fOutput)
: Elem(uLabel, Elem::GENEL, fOutput),
Genel(uLabel, Genel::STATESPACESISO, pDO, fOutput),
SD_y(y), SD_u(u),
iNumDofs(Order),
pdA(pA), pdB(pB), pdC(pC), dD(D),
pdX(NULL), pdXP(NULL) 
{
   ASSERT(Order > 0);
   ASSERT(pdA != NULL);
   ASSERT(pdB != NULL);
   ASSERT(pdC != NULL);	
   ASSERT(SD_y.iOrder == 0);
   ASSERT(SD_u.iOrder == 0);	
   DEBUGCOUT("GenelStateSpaceSISO " << uLabel 
	     << ", NumDofs: " << iNumDofs << endl);
   
   SAFENEWARR(pdX, doublereal, 2*Order);
   pdXP = pdX+Order;
}


GenelStateSpaceSISO::~GenelStateSpaceSISO(void) 
{
   if (pdX != NULL) {
      SAFEDELETEARR(pdX);
   }
   if (pdC != NULL) {
      SAFEDELETEARR(pdC);
   }
   if (pdB != NULL) {
      SAFEDELETEARR(pdB);
   }
   if (pdA != NULL) {
      SAFEDELETEARR(pdA);
   }
}


unsigned int GenelStateSpaceSISO::iGetNumDof(void) const 
{
   return iNumDofs;
}


/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order GenelStateSpaceSISO::SetDof(unsigned int i) const 
{
   ASSERT(i < iNumDofs);
   return DofOrder::DIFFERENTIAL;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& GenelStateSpaceSISO::Restart(ostream& out) const 
{
   return out << "GenelStateSpaceSISO: not implemented yet!" << endl; 
}


/* Dimensioni del workspace */
void GenelStateSpaceSISO::WorkSpaceDim(integer* piNumRows, 
				       integer* piNumCols) const 
{
   *piNumRows = iNumDofs+1;
   *piNumCols = iNumDofs+1;
}


void GenelStateSpaceSISO::Output(OutputHandler& /* OH */ ) const 
{
   NO_OP;
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
GenelStateSpaceSISO::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ ) 
{
   DEBUGCOUT("Entering GenelStateSpaceSISO::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;
   integer iColIndex_y = SD_y.pNode->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   
   WM.fPutRowIndex(iNumRows, iRowIndex_y);
   WM.fPutColIndex(iNumCols, iColIndex_y);
   
   WM.fPutCoef(iNumRows, iNumCols, dCoef);
   
   doublereal* pda = pdA+iNumDofs*iNumDofs-1;
   doublereal* pdc = pdC-1;
   for (unsigned int i = iNumDofs; i > 0; i--) {
      WM.fPutRowIndex(i, iFirstIndex+i);
      WM.fPutColIndex(i, iFirstIndex+i);
      WM.fPutCoef(iNumRows, i, -pdc[i]*dCoef);
      pda -= iNumDofs;
      for (unsigned int j = iNumDofs; j > 0; j--) {
	 /* Attenzione: si assume A orientata per righe:
	  * a_11, a_12, ..., a_1n, a_21, ..., a_2n, ..., a_nn */
	 WM.fPutCoef(i, j, -pda[j]*dCoef);
      }
      WM.fIncCoef(i, i, 1.);
   }
   
   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
GenelStateSpaceSISO::AssRes(SubVectorHandler& WorkVec,
			    doublereal /* dCoef */,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr) 
{
   DEBUGCOUT("Entering GenelStateSpaceSISO::AssRes()" << endl);
   
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   integer iRowIndex_y = SD_y.pNode->iGetFirstRowIndex()+1;    
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal y = SD_y.pNode->dGetX();
   doublereal u = SD_u.pNode->dGetX();
   
   WorkVec.fPutRowIndex(iNumRows, iRowIndex_y);
   
   doublereal* pdx = pdX-1;
   doublereal* pdxp = pdXP-1;            
   doublereal* pdc = pdC-1;
   doublereal d = dD*u-y;
   for (unsigned int i = iNumDofs; i > 0; i--) {
      WorkVec.fPutRowIndex(i, iFirstIndex+i);
      pdx[i] = XCurr.dGetCoef(iFirstIndex+i);
      pdxp[i] = XPrimeCurr.dGetCoef(iFirstIndex+i);
      d += pdc[i]*pdx[i];
   }
   WorkVec.fPutCoef(iNumRows, d);
   
   doublereal* pda = pdA+iNumDofs*iNumDofs;
   doublereal* pdb = pdB;
   pdxp = pdXP;     
   pdx = pdX;       
   for (unsigned int i = iNumDofs; i-- > 0; ) {
      d = pdb[i]*u-pdxp[i];      
      pda -= iNumDofs;
      for (unsigned int j = iNumDofs; j-- > 0; ) {
	 d += pda[j]*pdx[j];
      }
      WorkVec.fPutCoef(i+1, d);
   }
   
   return WorkVec;
}

/* GenelStateSpaceSISO - end */


/* GenelStateSpaceMIMO - begin */

GenelStateSpaceMIMO::GenelStateSpaceMIMO(unsigned int uLabel,
					 const DofOwner* pDO, 
					 unsigned int iNumOut,
					 const ScalarDof* y, 
					 unsigned int iNumIn,
					 const ScalarDof* u,
					 unsigned int Order,
					 doublereal* pA,
					 doublereal* pB,
					 doublereal* pC, 
					 doublereal* pD,
					 flag fOutput)
: Elem(uLabel, Elem::GENEL, fOutput),
Genel(uLabel, Genel::STATESPACEMIMO, pDO, fOutput),
iNumOutputs(iNumOut), iNumInputs(iNumIn),
pvSD_y((ScalarDof*)y), pvSD_u((ScalarDof*)u),
iNumDofs(Order),
pdA(pA), pdB(pB), pdC(pC), pdD(pD),
pdX(NULL), pdXP(NULL) 
{
#ifdef DEBUG
   ASSERT(iNumDofs > 0);
   ASSERT(iNumOutputs > 0);
   ASSERT(pvSD_y != NULL);
   for (int i = iNumOutputs; i-- > 0; ) {
      ASSERT(pvSD_y[i].iOrder == 0);
   }   
   ASSERT(iNumInputs > 0);
   ASSERT(pvSD_u != NULL);
   for (int i = iNumInputs; i-- > 0; ) {
      ASSERT(pvSD_u[i].iOrder == 0);
   }   
   ASSERT(pdA != NULL);
   ASSERT(pdB != NULL);
   ASSERT(pdC != NULL);	
#if 0
   ASSERT(pdD != NULL);	
#endif /* 0 */
   DEBUGCOUT("GenelStateSpaceMIMO " << uLabel 
	     << ", NumDofs: " << iNumDofs << endl);
#endif /* DEBUG */
   
   SAFENEWARR(pdX, doublereal, 2*Order);
   pdXP = pdX+Order;
}


GenelStateSpaceMIMO::~GenelStateSpaceMIMO(void) 
{
   if (pdX != NULL) {
      SAFEDELETEARR(pdX);
   }
   if (pdD != NULL) {
      SAFEDELETEARR(pdD);
   }
   if (pdD != NULL) {
      SAFEDELETEARR(pdC);
   }
   if (pdB != NULL) {
      SAFEDELETEARR(pdB);
   }
   if (pdA != NULL) {
      SAFEDELETEARR(pdA);
   }
   if (pvSD_u != NULL) {
      SAFEDELETEARR(pvSD_u);
   }   
   if (pvSD_y != NULL) {
      SAFEDELETEARR(pvSD_y);
   }   
}


unsigned int GenelStateSpaceMIMO::iGetNumDof(void) const 
{
   return iNumDofs;
}


/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order GenelStateSpaceMIMO::SetDof(unsigned int i) const 
{
   ASSERT(i < iNumDofs);
   return DofOrder::DIFFERENTIAL;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& GenelStateSpaceMIMO::Restart(ostream& out) const 
{
   return out << "GenelStateSpaceMIMO: not implemented yet!" << endl; 
}


/* Dimensioni del workspace */
void GenelStateSpaceMIMO::WorkSpaceDim(integer* piNumRows, 
				       integer* piNumCols) const 
{
   *piNumRows = iNumDofs+iNumOutputs;
   *piNumCols = iNumDofs+iNumOutputs;
}


void GenelStateSpaceMIMO::Output(OutputHandler& /* OH */ ) const 
{
   NO_OP;
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
GenelStateSpaceMIMO::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ ) 
{
   DEBUGCOUT("Entering GenelStateSpaceMIMO::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstIndex = iGetFirstIndex();
   doublereal* pdc = pdC+iNumOutputs*iNumDofs-1;
   for (unsigned int i = iNumOutputs; i > 0; i--) {
      integer iRowIndex_y = pvSD_y[i-1].pNode->iGetFirstRowIndex()+1;
      integer iColIndex_y = pvSD_y[i-1].pNode->iGetFirstColIndex()+1;

      WM.fPutRowIndex(iNumDofs+i, iRowIndex_y);
      WM.fPutColIndex(iNumDofs+i, iColIndex_y);
      WM.fPutCoef(iNumDofs+i, iNumDofs+i, dCoef); // 1 sulla diagonale

      pdc -= iNumDofs;
      for (unsigned int j = iNumDofs; j > 0; j--) {
	 /* Attenzione: si assume C orientata per righe:
	  * c_11, c_12, ..., c_1n, c_21, ..., c_2n, ..., c_nn */
	 WM.fPutCoef(iNumDofs+i, j, -pdc[j]*dCoef);
      }
   }

   doublereal* pda = pdA+iNumDofs*iNumDofs-1;
   for (unsigned int i = iNumDofs; i > 0; i--) {
      WM.fPutRowIndex(i, iFirstIndex+i);
      WM.fPutColIndex(i, iFirstIndex+i);
      pda -= iNumDofs;
      for (unsigned int j = iNumDofs; j > 0; j--) {
	 /* Attenzione: si assume A orientata per righe:
	  * a_11, a_12, ..., a_1n, a_21, ..., a_2n, ..., a_nn */
	 WM.fPutCoef(i, j, -pda[j]*dCoef);
      }
      WM.fIncCoef(i, i, 1.);
   }
   
   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
GenelStateSpaceMIMO::AssRes(SubVectorHandler& WorkVec,
			    doublereal /* dCoef */,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr) 
{
   DEBUGCOUT("Entering GenelStateSpaceMIMO::AssRes()" << endl);
   
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal* pdx = pdX-1;
   doublereal* pdxp = pdXP-1;
   for (unsigned int i = iNumDofs; i > 0; i--) {
      WorkVec.fPutRowIndex(i, iFirstIndex+i);
      pdx[i] = XCurr.dGetCoef(iFirstIndex+i);
      pdxp[i] = XPrimeCurr.dGetCoef(iFirstIndex+i);
   }   

   doublereal* pdc = pdC+iNumOutputs*iNumDofs-1;
   if (pdD != NULL) {      
      doublereal* pdd = pdD+iNumOutputs*iNumInputs-1;
      for (int i = iNumOutputs; i > 0; i--) {      
	 integer iRowIndex_y = pvSD_y[i-1].pNode->iGetFirstRowIndex()+1;    
	 WorkVec.fPutRowIndex(iNumDofs+i, iRowIndex_y);      
	 doublereal y = pvSD_y[i-1].pNode->dGetX();
	 doublereal d = -y;
	 pdc -= iNumDofs;
	 for (unsigned int j = iNumDofs; j > 0; j--) {
	    d += pdc[j]*pdx[j];
	 }
	 pdd -= iNumInputs;
	 for (unsigned int j = iNumInputs; j > 0; j--) {
	    d += pdd[j]*pvSD_u[j-1].pNode->dGetX();
	 }      
	 WorkVec.fPutCoef(iNumDofs+i, d);
      }
   } else {
      for (int i = iNumOutputs; i > 0; i--) {      
	 integer iRowIndex_y = pvSD_y[i-1].pNode->iGetFirstRowIndex()+1;    
	 WorkVec.fPutRowIndex(iNumDofs+i, iRowIndex_y);      
	 doublereal y = pvSD_y[i-1].pNode->dGetX();
	 doublereal d = -y;
	 pdc -= iNumDofs;
	 for (unsigned int j = iNumDofs; j > 0; j--) {
	    d += pdc[j]*pdx[j];
	 }
	 WorkVec.fPutCoef(iNumDofs+i, d);
      }
   }   
      
   doublereal* pda = pdA+iNumDofs*iNumDofs;
   doublereal* pdb = pdB+iNumDofs*iNumInputs;
   pdxp = pdXP;
   pdx = pdX;
   for (unsigned int i = iNumDofs; i-- > 0; ) {     
      doublereal d = -pdxp[i];
      pdb -= iNumInputs;
      for (unsigned int j = iNumInputs; j-- > 0; ) {
	 d += pdb[j]*pvSD_u[j].pNode->dGetX();
      }
      pda -= iNumDofs;
      for (unsigned int j = iNumDofs; j-- > 0; ) {
	 d += pda[j]*pdx[j];
      }
      WorkVec.fPutCoef(i+1, d);
   }

   return WorkVec;
}

/* GenelStateSpaceMIMO - end */

#endif /* USE_ELECTRIC_NODES */

