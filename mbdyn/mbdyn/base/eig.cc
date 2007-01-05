/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
 * Eigenvalue analysis by means of POD 
 * Copyright 2003-2007 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */ 
    

#include<ac/lapack.h>
#include "eig.h"

PODEig::PODEig(doublereal Period, 
	       integer S, 
	       integer N, 
	       doublereal Start,
	       doublereal Tr,
	       integer  NEig,
	       doublereal pres)
:T(Period),
iSSize(S),
iSnapN(N),
iSnapCount(0),
dStartTime(Start),
dTimeCount(Start),
dFinalTime(Start+Period*N)
dTresh(Tr),
iNumEig(NEig),
dPrecision(pres),
bEigDone(false)
{
	DEBUGCOUTFNAME("PODEig::PODEig");
	ASSERT(iSnapN > 0);
	pXMat = new std::vector<doublereal> [iSSize];
	for(int i=0; i<iSSize; i++) {
		pXMat[i].resize(iSnapN);
	}
}

		
PODEig::~PODEig(void)
{
	if (pXMat != NULL) {	
		delete pXMat;
	}
}
	
void PODEig::LogData(doublereal t, VectorHandler& AV, VectorHandler& PV)
{
	if (t > dFinalTime) return;
	else {
		if (t >= dTimeCount) {
			for(int i=0; i < iSSize; i++) {
				pXMat[i][iCountSnap] = AV.dGetCoef(i+1) - ( ( AV.dGetCoef(i+1) - PV.dGetCoef(i+1) ) * (t - dTimeCount)/T );   
			}	
    			dTimeCount += T;
			iSnapCount++;
		}
		return;
	}
			
}

struct less_mag : public binary_function<doublereal, doublereal, bool> {
	bool operator()(doublereal x, doublereal y) { return abs(x) < abs(y); }
};


void PODEig::ComputeEigenvalues(doublereal t)
{
	if (t <= dFinalTime) {
		silent_cerr("Warning: you cannot require to compute the eigenvalues " 
			 << "before getting all the required snapshot!! No Eig computed!!"
			 << std::endl);
		 return;
	}
	
	if (bEigDone) {
		silent_cerr("Warning: eigenvalues already obtained. No new eigenvalues computed"
				<< std::endl); 
		return;
	}
	
	/* depura ogni canale del valore medio */
	for (int i=0; i < iSnapN; i++) {
		std::vector<doublereal>::iterator iv;
		std::vector<doublereal>::const_iterator iv_end = pXMat[i].end();
		std::vector<doublereal>::const_iterator iM = std::max(pXMat[i].begin(), pXMat[i].end(), less_mag());
		for(iv = pXMat[i].begin(); iv != iv_end; iv++) { 
			*iv /= abs(*iM);	
		}
	}
	
	/* crea la matrice X^T X utilizzata per il calcolo degli SVD e dei relativi vettori 
	   la matrice e' simmetrica e viene salvata in formato lapack packed up triangualar */
	sdt::vector:<doublereal> RMat(int(iSnapN*(iSnapN+1/2)), 0);
	
	int j = 0;
	for (int iCol=0; iCol < iSnapN; iCol++) { 
		for(int iRow=0; iRow <= iCol; iRow++) {
			for(int k=0; k < iSSize; k++) {		
				RMat[j++] += pXMat[k][iCol] * pXMat[k][iRow];
			}
		}
	}
	
	char sJ[2] = "V"; /* calcola autoval e autovett */
   	char sR[2];
	
	if (iNumEig) sR[0] = "I";
	else sR[0] = "A";
	
	char sUP[2] ="U";
	
	doublereal zero = 0.;
	
	integer ilb, iub = ilb = 1;
	if (iNumEig) iub = iNumEig;
	
	doublereal dTol = 1e-6;
	if (dPrecision != 0) dTol = dPrecision;
	
	integer NEig;    /* numero di autovalori restituiti come output */	
	doublereal* pEig = new doublereal [iSnapN];			/*vettore degli autovalori */
	
	doublereal* pEigVec = new doublereal [iSnapN*iSnapN];
	 	
	integer* piSupp = new integer [2*iSnapN];
	
	doublereal* pdW = new doublereal [2];
	
	integer iWDim = -1;
	
	integer* piW = new integer [2];
	
	integer iWIDim = -1;
	
	integer iInfo;
			 
	__FC_DECL__(dsyevr)(sJ, SR, sUP, 
			    &iSnapN,
			    &RMat[0], &iSnapN,
			    &zero, &zero,
			    &ilb, &iub,
			    &dTol,
			    &NEig, pEig, pEigVec,
			    &iSnapN, piSupp,
			    pdW,
			    &iWDim,
			    piW,
			    &iWIDim,
			    &iInfo);
			    
	if (iInfo != 0) {
      		/* = 0:  successful exit */
      		silent_cerr("Eigenvalues error computing lapack work spaces dimensions !!"
				<< std::endl);
		return;
	}
	 
	delete pdW;
	delete piW;
	
	doublereal* pdW = new doublereal [iWDim];
	integer*    piW = new integer [iWIDim];
	
	__FC_DECL__(dsyevr)(sJ, SR, sUP, 
			    &iSnapN,
			    &RMat[0], &iSnapN,
			    &zero, &zero,
			    &ilb, &iub,
			    &dTol,
			    &NEig, pEig, pEigVec,
			    &iSnapN, piSupp,
			    pdW,
			    &iWDim,
			    piW,
			    &iWIDim,
			    &iInfo);
			    
   	if (iInfo < 0) {
      		char *th = "th";

      		/* Aaaaah, English! :) */
      		if (-iInfo/10 != 10) {
         		switch ((-iInfo+20)%10) {
         		case 1:
	    			th = "st";
	    			break;
         		case 2:
	    			th = "nd";
	    			break;
	 		case 3:
	    			th = "rd";
	    			break;
	 		}
      		}
      		/* < 0:  if INFO = -i, the i-th argument had an illegal value. */
      		silent_cerr("the " << -iInfo << "-" << th 
	      		  << " argument had an illegal value" << std::endl);
   	}
   
   /* ottenuti i valori singolari e i modi
        1) calcolo U^T*X^T;
   	2)  calcolo gli svd come norma delle righe di U^T*X^T;
   	3) li stampo in output (quelli al di sotto del Tresh) aasieme agli autovettori (righe di U^T*X^T normalizzate); 
	4) calcolo le storie temporali dei pom
    */
   
	FullMatrixHandler dPOM(iSState, iSnapN, 0);
	for (int iCol=0; iCol < NEig; iCol++) { 
		for(int iRow=0; iRow <= iSSize; iRow++) {
			for(int k=0; k < iSnapN; k++) {		
				dPOM.IncCoef(iCol+1, iRow+1, pEigVec[piSupp[2*NEig+1]+k] * pXMat[iRow][k]);
			}
		}
	}
	
	std::vector<doublereal> SVD(iSnapN, 0);
	doublereal dd;
	for (int iS=0; iS < iSnapN; iS++) { 
		for (int iRow=0; iRow < iSSize; iRow++) { 
			dd = dPOM.dGetCoef(iSSize+1, iS+1);
			SVD[iS] += dd*dd;
		}
		SVD[iS] = sqrt(SVD[iS]);
	}
	
	doublereal svTot = 0;
	for (int iS=0; iS < iSnapN; iS++) { 
		svTot += SVN[iS];
		for (int iRow=0; iRow < iSSize; iRow++) { 
			dd = dPOM.dGetCoef(iSSize+1, iS+1);
			dPOM.PutCoef(iSSize+1, iS+1, dd/SVD[iS]);
		}
	}
		
	integer iNumValidPOM = 0;
	doublereal tmpSum;	
	do  {
		tmpSum = SVD[iNumValidPOM++]/svTot;
	} while	(tmpSum < dTresh);
	
	
	/* OUTPUT */
	
	FullMatrixHandler dPOMTimeHs(iNumValidPOM, iSnapN, 0);
	
	/* dPOM^T * dXMat */
	for (int iRow=0; iRow < iNumValidPOM; iRow++) { 
		for(int iCol=0; iCol <= iSnapN; iCol++) {
			for(int k=0; k < iSSize; k++) {		
				dPOMTimeHs.IncCoef(iRow+1, iCol+1, dPOM.dGetCoef(iRow+1, k+1) * pXMat[k][iCol]);
			}
		}
	}
	
	
	
		
			
    /* identificazione delle storie temporali dei pom con l'algoritmo AR */
			
    /* calcolo degli autovalori della matrice A con le lapack ancora */			 

}


		
void PODEig::OutputEigenvalues(void)
		
#endif //PODEIG_H
