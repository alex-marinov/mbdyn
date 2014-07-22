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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "forgfact.h"


ForgettingFactor::ForgettingFactor(integer i) 
: iNumErr(i) 
{
   ASSERT(iNumErr >= 0);
}


ForgettingFactor::~ForgettingFactor(void) 
{
   NO_OP;
}


ConstForgettingFactor::ConstForgettingFactor(doublereal d) 
: ForgettingFactor(0), dk(d) 
{
   NO_OP;
}


ConstForgettingFactor::~ConstForgettingFactor(void) 
{
   NO_OP;
}
 

void ConstForgettingFactor::Update(const doublereal* /* pErr */ )
{
   NO_OP;
}


DynamicForgettingFactor::DynamicForgettingFactor(integer n1,
						 integer n2, 
						 integer i,
						 doublereal r, 
						 doublereal f,
						 doublereal kr, 
						 doublereal kl)
: ForgettingFactor(i), 
N1(n1), N2(n2), dRho(r), dFact(f), dkRef(kr), dkLim(kl),
dk(kr), pdErrM(NULL), pdErrS(NULL), iRef1(n1-1), iRef2(n1-n2-1),
dErr1M(0.), dErr1S(0.), dErr2M(0.), dErr2S(0.) {
   ASSERT(N1 > 0);
   ASSERT(N2 > 0);
   ASSERT(N2 < N1);
   ASSERT(dRho > 0. && dRho < 1.);
   ASSERT(dFact > 0. && dFact < 1.);
   ASSERT(dkRef > 0. && dkRef < 1.);
   ASSERT(dkLim > dkRef && dkLim <= 1.);
   
   SAFENEWARR(pdErrM, doublereal, 2*N1);
   
   pdErrS = pdErrM+N1;
   for (integer i = 0; i < 2*N1; i++) {
      pdErrM[i] = 0.;
   }	
}

DynamicForgettingFactor::~DynamicForgettingFactor(void) {
   /* pdErr must be non-null */
   SAFEDELETEARR(pdErrM);
}

void DynamicForgettingFactor::Update(const doublereal* pErr) {
   /* setta i riferimenti */
   if (++iRef1 == N1) {
      iRef1 = 0;
   }
   if (++iRef2 == N1) {
      iRef2 = 0;
   }
   
   /* calcola la media e la norma dell'errore */
   doublereal dM = 0.;
   doublereal dS = 0.;
   for (integer i = iNumErr; i--;) {
      dS += pErr[i]*pErr[i];
      dM += pErr[i];
   }
   
   /* somma agli errori la norma dell'errore corrente, poi sottrae
    * l'errore al limite posteriore delle rispettive finestre */
   dErr1M += dM;
   dErr1M -= pdErrM[iRef1];
   dErr1S += dS;
   dErr1S -= pdErrS[iRef1];
   
   dErr2M += dM;
   dErr2M -= pdErrM[iRef2];
   dErr2S += dS;
   dErr2S -= pdErrS[iRef2];
   
   /* aggiorna il limite superiore della finestra "lunga" */
   pdErrM[iRef1] = dM;
   pdErrS[iRef1] = dS;
   
   /* differenza tra la varianza dell'errore "corto"
    * rispetto a quello "lungo" */
   doublereal r = (dErr1S/N1-pow(dErr1M/N1,2));
   if (fabs(r) < 1.e-6) {
      return;
   }
   r = 1.-(dErr2S/N2-pow(dErr2M/N2,2))/r;	 
   
   /* se viene superato il limite, il forgetting factor scatta al valore 
    * minimo (massima dimenticanza) */
   if (fabs(r) > dFact) {
      dk = dkRef;
   } else {
      
      /* il forgetting factor viene fatto variare con una sua dinamica */
      dk = dRho*dk+(1-dRho)*dkLim;
   }
}


DynamicForgettingFactor2::DynamicForgettingFactor2(integer n1,
						   integer n2, 
						   integer i,
						   doublereal r,
						   doublereal f,
						   doublereal kr,
						   doublereal kl)
: ForgettingFactor(i), 
N1(n1), N2(n2), dRho(r), dFact(f), dkRef(kr), dkLim(kl),
dk(kr), pdErr(NULL), ppdErr(NULL), iRef1(n1-1), iRef2(n1-n2-1),
pdErr1M(NULL), pdErr1S(NULL), pdErr2M(NULL), pdErr2S(NULL) {
   ASSERT(N1 > 0);
   ASSERT(N2 > 0);
   ASSERT(N2 < N1);
   ASSERT(dRho > 0. && dRho < 1.);
   ASSERT(dFact > 0. && dFact < 1.);
   ASSERT(dkRef > 0. && dkRef < 1.);
   ASSERT(dkLim > dkRef && dkLim <= 1.);
   
   integer sz = (N1+4)*iNumErr;
   
   SAFENEWARR(pdErr, doublereal, sz);
   SAFENEWARR(ppdErr, doublereal*, N1);
   
   pdErr1M = pdErr+N1*iNumErr;
   pdErr1S = pdErr1M+iNumErr;
   pdErr2M = pdErr1S+iNumErr;
   pdErr2S = pdErr2M+iNumErr;
   for (integer i = sz; i-- > 0; ) {
      pdErr[i] = 0.;
   }
   for (integer i = N1; i-- > 0; ) {
      ppdErr[i] = pdErr+i*iNumErr;
   }	
}

DynamicForgettingFactor2::~DynamicForgettingFactor2(void) {
   /* pdErr must be non-null */
   SAFEDELETEARR(pdErr);
   SAFEDELETEARR(ppdErr);
}

void DynamicForgettingFactor2::Update(const doublereal* pE) {
   /* setta i riferimenti */
   if (++iRef1 == N1) {
      iRef1 = 0;
   }
   if (++iRef2 == N1) {
      iRef2 = 0;
   }
   
   /* somma agli errori la norma dell'errore corrente, poi sottrae
    * l'errore al limite posteriore delle rispettive finestre */    
   doublereal dS1 = 0.;  
   doublereal dS2 = 0.;
   for (integer i = iNumErr; i-- > 0; ) {
      doublereal d = pE[i]*pE[i];
      pdErr1S[i] += d;
      pdErr1S[i] -= ppdErr[iRef1][i]*ppdErr[iRef1][i];
      pdErr1M[i] += pE[i];
      pdErr1M[i] -= ppdErr[iRef1][i];
      dS1 += pdErr1S[i]/N1-pow(pdErr1M[i]/N1,2);
      pdErr2S[i] += d;
      pdErr2S[i] -= ppdErr[iRef2][i]*ppdErr[iRef2][i];
      pdErr2M[i] += pE[i];
      pdErr2M[i] -= ppdErr[iRef2][i];
      dS2 += pdErr2S[i]/N2-pow(pdErr2M[i]/N2,2);
      
      /* aggiorna il limite superiore della finestra "lunga" */
      ppdErr[iRef1][i] = pE[i];       
   }

   
   
   /* per debug
   cout << setw(16) << dS1 << setw(16) << dS2 << setw(16) 
     << (dS1 != 0. ? 1.-dS2/dS1 : 0.) << endl;
    */
   
   
   
   /* differenza tra la varianza dell'errore "corto"
    * rispetto a quello "lungo"
    * se viene superato il limite, il forgetting factor scatta al valore 
    * minimo (massima dimenticanza) */
   doublereal d = 0.;
   if (fabs(dS1) > std::numeric_limits<doublereal>::epsilon() && (d = fabs(1.-dS2/dS1)) > dFact) {
#if 0
      dk = (dkLim-dkRef)*exp(-d)+dkRef;
#else 
      dk = dkRef;
#endif
   } else {
      /* il forgetting factor viene fatto variare con una sua dinamica */
      dk = dRho*dk+(1-dRho)*dkLim;
   }
}

