/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* identificatore generico */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "id.h"

/* Ident - begin */

Ident::Ident(integer size, integer nout, ForgettingFactor* pf, 
	     const doublereal& ldl_init)
: size(size), nout(nout), pdBase(NULL), ppdBase(NULL),
ldl(NULL), vldl(NULL), z(0), vz(0), theta(NULL), vtheta(NULL),
phi(NULL), y(NULL), err(NULL), pF(pf), k(0.), w(0.)
{
   ASSERT(size > 0);
   ASSERT(nout > 0);
   
   integer i = size*size               // ldl
     +size*nout                        // z
     +size*nout                        // theta
     +size                             // phi
     +nout                             // y
     +nout;                            // err
     
   integer iv = size                   // vldl
     +nout                             // vz
     +nout;                            // vtheta

   SAFENEWARR(pdBase, doublereal, i);
   SAFENEWARR(ppdBase, doublereal*, iv);
   
   ldl = pdBase;
   z = ldl+size*size;
   theta = z+size*nout;
   phi = theta+size*nout;
   y = phi+size;
   err = y+nout;

   for (integer j = i; j-- > 0; ) {
      pdBase[j] = 0.;
   }

   vldl = ppdBase;
   vz = vldl+size;
   vtheta = vz+nout;

   for (integer j = size; j-- > 0; ) {
      vldl[j] = ldl+j*size;
      vldl[j][j] = ldl_init;
   }
   for (integer j = nout; j-- > 0; ) {
      vz[j] = z+j*size;
      vtheta[j] = theta+j*size;
   }   
}


Ident::~Ident(void) {
   SAFEDELETEARR(ppdBase);
   SAFEDELETEARR(pdBase);
   SAFEDELETE(pF);   
}

doublereal Ident::dGetForgettingFactor(void) const {
   return 1./(k*k);
}
   
void Ident::SetForgettingFactor(const doublereal& kk) {
   k = kk;
   w *= kk;
   if (w > 1.e3) {
      for (integer i = size; i-- > 0; ) {       
	 for (integer j = nout; j-- > 0; ) {
	    vz[j][i] /= w;
	 }
	 vldl[i][i] *= w;
      }
      w = 1.;
   } else if (w < 1.e-20) {
      w = 1.e-20;
   }
}

doublereal* Ident::pdGetTheta(void) {
   // copia il termine noto z in theta (richiesto da ldlsol_())
   for (integer i = size*nout; i-- > 0; ) {
      theta[i] = z[i];
   }

   __FC_DECL__(ldlsol) (ldl, &size, theta, &size, &size, &nout);
   return theta;
}

/* in realta' si va a correggere theta senza che Ident lo sappia! */
void Ident::UpdateTheta(const doublereal* pd)
{
   for (integer i = size*nout; i-- > 0; ) {
      theta[i] = pd[i];
   }
}

doublereal* Ident::pdGetErr(void) {
   return err;
}

void Ident::Update(const doublereal* pphi, const doublereal* yy) {
   // calcolo l'errore
   for (integer j = nout; j-- > 0; ) {
      err[j] = yy[j];
      for (integer i = size; i-- > 0; ) {
	 err[j] -= vtheta[j][i]*pphi[i];
      }
   }
   
   /* calcolo il nuovo forgetting factor */
   pF->Update(err);
   doublereal d = pF->dGet();
   SetForgettingFactor(k = 1./sqrt(d));

   for (integer i = size; i-- > 0; ) {
      phi[i] = pphi[i]*k;
   }
   for (integer j = nout; j-- > 0; ) {
      y[j] = yy[j]*k;
   }
   
   /* aggiorno la matrice fattorizzata */
   __FC_DECL__(uldlad) (ldl, &size, &size, phi, z, &size, &nout, y);
}

/* Ident - end */


/* IdentProcess - begin */

IdentProcess::IdentProcess(unsigned int iOut, unsigned int iIn,
			   unsigned int iA, unsigned int iB)
: iNumOutput(iOut), iNumInput(iIn), iOrdA(iA), iOrdB(iB), pIdent(NULL) 
{
   ASSERT(iOut > 0);
   ASSERT(iIn > 0);
   // ASSERT(iOrdA >= 0);  // = 0 in FIR
   // ASSERT(iOrdB >= 0);  // = 0 in AR
}
   
IdentProcess::~IdentProcess(void) 
{
   SAFEDELETE(pIdent);
}
   
void IdentProcess::CreateIdent(integer size, integer nout, ForgettingFactor* pf) 
{
   ASSERT(size > 0);
   ASSERT(nout > 0);
   ASSERT(pf != NULL);     
   
   SAFENEWWITHCONSTRUCTOR(pIdent, Ident, Ident(size, nout, pf));
}

void IdentProcess::GetErr(doublereal* pdE) 
{
   NO_OP;
}

/* IdentProcess - end */


/* IdentARXProcess - begin */

IdentARXProcess::IdentARXProcess(unsigned int iOut, unsigned int iIn,
				 unsigned int iA, unsigned int iB,
				 ForgettingFactor* pf)
: IdentProcess(iOut, iIn, iA, iB), 
size(iOut*iA+iIn*(iB+1)), pdBase(NULL), pdPhi(NULL), pdY(NULL)
{
   // size e' il lato lungo di Theta
   
   CreateIdent(size, integer(iOut), pf);
   
   integer i = size        // Phi
     +iOut;                // Y
     
   SAFENEWARR(pdBase, doublereal, i);
   
   pdPhi = pdBase;
   pdY = pdPhi+size;
   
   for (integer j = i; j-- > 0; ) {
      pdBase[j] = 0.;
   }	
}
 
IdentARXProcess::~IdentARXProcess(void)
{
   SAFEDELETEARR(pdBase);
}  
   
void IdentARXProcess::Update(const doublereal* pdYTmp, 
			     const doublereal* pdUTmp)
{
   if (iOrdA > 0) {
      for (integer i = iNumOutput*(iOrdA-1); i-- > 0; ) {       
	 pdPhi[i+iNumOutput] = pdPhi[i];
      }      
      for (integer i = iNumOutput; i-- > 0; ) {       
	 pdPhi[i] = pdY[i];
	 pdY[i] = pdYTmp[i];
      }
   } else {
      // FIR
      for (integer i = iNumOutput; i-- > 0; ) {
	 pdY[i] = pdYTmp[i];
      }
   }

   doublereal* pdTmp = pdPhi+iNumOutput*iOrdA;
   // if (iOrdB == 0) // AR is implicitly satisfied
   for (integer i = iNumInput*iOrdB; i-- > 0; ) {
      pdTmp[i+iNumInput] = pdTmp[i];
   }
   for (integer i = iNumInput; i-- > 0; ) {
      pdTmp[i] = pdUTmp[i];
   }

   pIdent->Update(pdPhi, pdY);
}

void IdentARXProcess::GetTheta(doublereal* pdTheta)
{
   doublereal* p = pIdent->pdGetTheta();
   for (integer i = iNumOutput*size; i-- > 0; ) {
      pdTheta[i] = p[i];
   }
}

/* IdentARXProcess - end */


/* IdentARMAXProcess - begin */

IdentARMAXProcess::IdentARMAXProcess(unsigned int iOut, unsigned int iIn,
				     unsigned int iA, unsigned int iB,
				     ForgettingFactor* pf)
: IdentProcess(iOut, iIn, iA, iB), 
size(2*iOut*iA+iIn*(iB+1)), pdBase(NULL), pdPhi(NULL), pdY(NULL), pdErr(NULL) 
{
   // size e' il lato lungo di Theta.
   
   CreateIdent(size, integer(iOut), pf);
   
   integer i = size        // Phi
     +iOut                 // Y
     +iOut;                // Err
     
   SAFENEWARR(pdBase, doublereal, i);
   
   pdPhi = pdBase;
   pdY = pdPhi+size;
   pdErr = pdY+iOut;
 
   
   for (integer j = i; j-- > 0; ) {
      pdBase[j] = 0.;
   }	
}
 
IdentARMAXProcess::~IdentARMAXProcess(void) {    
   SAFEDELETEARR(pdBase);
}  


// questa operazione deve essere svolta ad ogni passo (del controllo)
void IdentARMAXProcess::Update(const doublereal* pdYTmp, 
			       const doublereal* pdUTmp) 
{
   // qui iniziano gli errori
   doublereal* pdTmp = pdPhi+iNumOutput*iOrdA+iNumInput*(iOrdB+1);
   
   // sposta di un passo (iNumOutput) le uscite e gli errori
   for (integer i = iNumOutput*(iOrdA-1); i-- > 0; ) {
      pdPhi[i+iNumOutput] = pdPhi[i]; // Y vecchi
      pdTmp[i+iNumOutput] = pdTmp[i]; // Err vecchi
   }
   // aggiunge in testa le uscite e gli errori al passo precedente
   for (integer i = iNumOutput; i-- > 0; ) {       
      pdPhi[i] = pdY[i];    // Y all'ultima iterazione
      pdY[i] = pdYTmp[i];   // Y nuovo
      pdTmp[i] = pdErr[i];  // Err all'ultima iterazione 
   }

   // qui iniziano gli ingressi
   pdTmp = pdPhi+iNumOutput*iOrdA;
   // sposta di un passo (iNumInput) gli ingressi
   // if (iOrdB == 0) // ARMA is implicitly satisfied 
   for (integer i = iNumInput*iOrdB; i-- > 0; ) {      
      pdTmp[i+iNumInput] = pdTmp[i];
   }
   // aggiunge in testa gli ingressi al passo corrente
   for (integer i = iNumInput; i-- > 0; ) {     
      pdTmp[i] = pdUTmp[i];
   }
   
   // fa aggiornare la stima di Theta
   pIdent->Update(pdPhi, pdY);
 
   // err al passo corrente e' stato calcolato da Ident; lo ricopia
   doublereal* pd = pIdent->pdGetErr();
   for (integer i = iNumOutput; i-- > 0; ) {
      pdErr[i] = pd[i];
   }
}


// questa operazione puo' essere svolta in modo asincrono, ovvero ogni tanto
void IdentARMAXProcess::GetTheta(doublereal* pdTheta) 
{
   // fa calcolare il nuovo theta
   doublereal* p = pIdent->pdGetTheta();

   // se e' il caso, corregge la parte di MA di Theta
   if (fCheckMA(p)) {
      // se la parte MA di theta e' stata corretta, e' stata aggiornata in loco
      NO_OP;
   }   
   
   // restituisce il nuovo theta
   for (integer i = iNumOutput*size; i-- > 0; ) {     
      pdTheta[i] = p[i];
   }   
}

void IdentARMAXProcess::GetErr(doublereal* pdE) 
{
   for (integer i = iNumOutput; i-- > 0; ) {
      pdE[i] = pdErr[i];
   }          
}
 
flag IdentARMAXProcess::fCheckMA(doublereal* /* pdTheta */ )
{
   // fa la correzione di Theta in base agli autovalori ed autovettori
   return flag(0);
}

/* IdentARMAXProcess - end */

