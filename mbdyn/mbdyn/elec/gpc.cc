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

/* funzioni di base per le operazioni con le matrici */

#include <mbconfig.h>

#include <myassert.h>
#include <mynewmem.h>

extern "C" {
#include <mymath.h>
}

#include <gpc.h>

/* gpc utility functions - begin */

/*
 * esegue d = m1*m2+m3 
 */

int gpc_addmul(integer ndim1, integer nrow1, integer ncol1, const doublereal* m1,
	       integer ndim2, integer nrow2, integer ncol2, const doublereal* m2,
	       integer ndim3, integer nrow3, integer ncol3, const doublereal* m3,
	       integer ndimd, integer nrowd, integer ncold, doublereal* d)
{
   ASSERT(m1 != NULL);
   ASSERT(m2 != NULL);
   ASSERT(m3 != NULL);
   ASSERT(d != NULL);
   
   ASSERT(nrow2 == ncol1);
   ASSERT(nrow3 == nrow1);
   ASSERT(ncol3 == ncol2);
   ASSERT(nrowd == nrow1);
   ASSERT(ncold == ncol2);
   
   ASSERT(ndim1 >= nrow1);
   ASSERT(ndim2 >= ncol1);
   ASSERT(ndim3 >= nrow1);
   ASSERT(ndimd >= nrow1);
   
   doublereal* mm1 = (doublereal*)m1;
   doublereal* mm2 = (doublereal*)m2+ncol2*ndim2;
   doublereal* mm3 = (doublereal*)m3+ncol3*ndim3;
   doublereal* dd = d+ncold*ndimd;
   
   for (integer c = ncold; c-- > 0; ) {
      dd -= ndimd-nrowd;
      mm2 -= ndim2;
      mm3 -= ndim3-nrow3;
      mm1 += nrow1;
      
      for (integer r = nrowd; r-- > 0; ) {
	 dd--;
	 mm1--;
	 mm3--;
	 
	 dd[0] = mm3[0];	
	 for (integer k = ncol1; k-- > 0; ) {
	    dd[0] += mm1[k*ndim1]*mm2[k];
	 }
      }
   }
      
   return 0;
}


/*
 * esegue d = m1*m2
 */

int gpc_mul(integer ndim1, integer nrow1, integer ncol1, const doublereal* m1,
	    integer ndim2, integer nrow2, integer ncol2, const doublereal* m2,
	    integer ndimd, integer nrowd, integer ncold, doublereal* d)
{
   ASSERT(m1 != NULL);
   ASSERT(m2 != NULL);
   ASSERT(d != NULL);
   
   ASSERT(nrow2 == ncol1);
   ASSERT(nrowd == nrow1);
   ASSERT(ncold == ncol2);
   
   ASSERT(ndim1 >= nrow1);
   ASSERT(ndim2 >= ncol1);
   ASSERT(ndimd >= nrow1);

   doublereal* mm1 = (doublereal*)m1;
   doublereal* mm2 = (doublereal*)m2+ncol2*ndim2;
   doublereal* dd = d+ncold*ndimd;

   for (integer c = ncold; c-- > 0; ) {
      dd -= ndimd-nrowd;
      mm2 -= ndim2;	
      mm1 += nrow1;
      
      for (integer r = nrowd; r-- > 0; ) {
	 dd--;
	 mm1--;	 
	 
	 dd[0] = 0.;
	 for (integer k = ncol1; k-- > 0; ) {
	    dd[0] += mm1[k*ndim1]*mm2[k];
	 }
      }
   }
      
   return 0;
}


/*
 * esegue d = m1*m2; se m3 e' definita, d = m1*m2+m3
 */

int gpc_add_mul(integer ndim1, integer nrow1, integer ncol1, const doublereal* m1,
		integer ndim2, integer nrow2, integer ncol2, const doublereal* m2,
		integer ndim3, integer nrow3, integer ncol3, const doublereal* m3,
		integer ndimd, integer nrowd, integer ncold, doublereal* d)
{
   if (m3) {
      return gpc_addmul(ndim1, nrow1, ncol1, m1,
			ndim2, nrow2, ncol2, m2,
			ndim3, nrow3, ncol3, m3,
			ndimd, nrowd, ncold, d);
   }
   return gpc_mul(ndim1, nrow1, ncol1, m1,
		  ndim2, nrow2, ncol2, m2,
		  ndimd, nrowd, ncold, d);   
}


/*
 * esegue d = m1*m2+m3 senza check sulle dimensioni
 */
 
int gpc_addmul(integer ndim1, integer nrow1, integer ncol1, const doublereal* m1,
	       integer ndim2, integer ncol2, const doublereal* m2,
	       integer ndim3, const doublereal* m3,
	       integer ndimd, doublereal* d)
{
   ASSERT(m1 != NULL);
   ASSERT(m2 != NULL);
   ASSERT(m3 != NULL);
   ASSERT(d != NULL);

   ASSERT(ndim1 >= nrow1);
   ASSERT(ndim2 >= ncol1);
   ASSERT(ndim3 >= nrow1);
   ASSERT(ndimd >= nrow1);
   
   doublereal* dd = d+ncol2*ndimd;
   doublereal* mm1 = (doublereal*)m1;
   doublereal* mm2 = (doublereal*)m2+ncol2*ndim2;
   doublereal* mm3 = (doublereal*)m3+ncol2*ndim3;
   
   for (integer c = ncol2; c-- > 0; ) {
      dd -= ndimd-nrow1;
      mm2 -= ndim2;
      mm3 -= ndim3-nrow1;
      mm1 += nrow1;
      
      for (integer r = nrow1; r-- > 0; ) {
	 dd--;
	 mm1--;
	 mm3--;
	 
	 dd[0] = mm3[0];	
	 for (integer k = ncol1; k-- > 0; ) {
	    dd[0] += mm1[k*ndim1]*mm2[k];
	 }
      }
   }
   
   return 0;
}


/*
 * copia una matrice su un'altra
 */

int gpc_mcopy(integer ndims, integer nrows, integer ncols, const doublereal* s,
	      integer ndimd, doublereal* d)
{
   doublereal* ss = (doublereal*)s+ndims*ncols;
   doublereal* dd = d+ndimd*ncols;
   
   for (integer j = ncols; j-- > 0; ) {
      ss -= ndims;
      dd -= ndimd;
      for (integer i = nrows; i-- > 0; ) {
	 dd[i] = ss[i];
      }
   }

   return 0;
}


/*
 * copia una matrice orientata per righe su una orientata per colonne
 */

int gpc_mcopy_t(integer ndims, integer nrows, integer ncols, const doublereal* s,
		integer ndimd, doublereal* d)
{
   doublereal* ss = (doublereal*)s+ncols;
   doublereal* dd = d+ndimd*ncols;
   
   for (integer j = ncols; j-- > 0; ) {   
      dd -= ndimd;
      ss--;
      for (integer i = nrows; i-- > 0; ) {
	 dd[i] = ss[ndims*i];
      }
   }

   return 0;
}


/*
 * azzera un vettore
 */

int gpc_zero(integer size, doublereal* v)
{
   for (integer i = size; i-- > 0; ) {
      v[i] = 0.;
   }
   return 0;
}


/*
 * setta un vettore
 */

int gpc_set(integer size, doublereal* v, const doublereal d)
{
   for (integer i = size; i-- > 0; ) {
      v[i] = d;
   }
   return 0;
}


/*
 * costruisce le matrici grezze di predizione
 */

int gpc_build_matrices(integer ndimA, doublereal* A,
		       integer ndimB, doublereal* B,
		       integer ndimP, doublereal* P,
		       integer ndimC, doublereal* C,
		       integer nout, integer nin,
		       integer s, integer pa, integer pb,
		       doublereal* Theta)
{
   /*
    * matrici: A, B, P (, C se definita) 
    * dimensioni: dimA, dimB, dimP (, dimC se definita) = s*nout
    *             ncolA = pa*nout
    *             ncolB = pb*nin
    *             ncolP = s*nin
    *              (ncolC = pa*nout se definita)
    */
   
   /* inizializzo le matrici (attenzione che theta e' organizzata per righe) */
   integer size = nout*pa+nin*(pb+1);
   if (C != NULL && pa > 0) {
      size += nout*pa;
   }

   if (pa > 0) {
      /* inizializza A con a_1 .. a_pa */
      doublereal* source = Theta;
      doublereal* dest = A+nout*(s-1);
      gpc_mcopy_t(size, nout, nout*pa, source, ndimA, dest);

      if (C != NULL) {
	 /* inizializza c con c_1 .. c_pa */
	 source = Theta+nout*pa+nin*(pb+1);
	 dest = C+nout*(s-1);
	 gpc_mcopy_t(size, nout, nin*(pb+1), source, ndimB, dest);
      }      
   }
 
   /* inizializza B con b_1 .. b_pb */
   doublereal* source = Theta+nout*pa+nin;
   doublereal* dest = B+nout*(s-1);
   gpc_mcopy_t(size, nout, nin*pb, source, ndimB, dest);
      
   /* inizializza P con b_0 */
   gpc_zero(ndimP*(nin*s), P);
   source = Theta+nout*pa;
   dest = P+nout*(s-1)+ndimP*(nin*(s-1));
   gpc_mcopy_t(size, nout, nin, source, ndimP, dest);
 
   /* le matrici si assumono gia' inizializzate;
    * significa che all'ultimo blocco-riga contengono:
    *     A le matrici a_1 .. a_pa,
    *     B le matrici b_1 .. b_pb,
    *     P zeri terminati dalla matrice b_0,
    *     ( C le matrici c_1 .. c_pa se definita)
    */

   /* ciclo principale */
   for (integer i = s-1; i > 0; i--) {
      doublereal* m1 = A+nout*i; /* a_1^(l-1) */

      if (pa > 0) {
	 for (integer j = 1; j <= pa-1; j++) {
	    /* a_j^l = a_1^(l-1)*a_j^0+a_(j+1)^(l-1) */
	    doublereal* dest = A+nout*(i-1)+ndimA*(nout*(j-1)); /* a_j^l */
	    doublereal* m2 = A+nout*(s-1)+ndimA*(nout*(j-1));   /* a_j^0 */
	    doublereal* m3 = A+nout*i+ndimA*(nout*j);           /* a_(j+1)^(l-1) */
	    gpc_addmul(ndimA, nout, nout, m1,
		       ndimA, nout, nout, m2,
		       ndimA, nout, nout, m3,
		       ndimA, nout, nout, dest);
	 }

	 /* a_pa^l = a_1^(l-1)*a_pa^0 */
	 doublereal* dest = A+nout*(i-1)+ndimA*(nout*(pa-1)); /* a_pa^l */
	 doublereal* m2 = A+nout*(s-1)+ndimA*(nout*(pa-1));   /* a_pa^0 */
	 gpc_mul(ndimA, nout, nout, m1,
		 ndimA, nout, nout, m2,
		 ndimA, nout, nout, dest);
      }

      if (pb > 0) {	 
	 for (integer j = 1; j <= pb-1; j++) {
	    /* b_j^l = a_1^(l-1)*b_j^0+b_(j+1)^(l-1) */
	    doublereal* dest = B+nout*(i-1)+ndimB*(nin*(j-1)); /* b_j^l */
	    doublereal* m2 = B+nout*(s-1)+ndimB*(nin*(j-1));   /* b_j^0 */
	    doublereal* m3 = B+nout*i+ndimB*(nin*j);           /* b_(j+1)^(l-1) */
	    gpc_addmul(ndimA, nout, nout, m1,
		       ndimB, nout, nin, m2,
		       ndimB, nout, nin, m3,
		       ndimB, nout, nin, dest);
	 }
	 
	 /* b_pb^l = a_1^(l-1)*b_pb^0 */
	 doublereal* dest = B+nout*(i-1)+ndimB*(nin*(pb-1)); /* b_pb^l */
	 doublereal* m2 = B+nout*(s-1)+ndimB*(nin*(pb-1));   /* b_pb^0 */
	 gpc_mul(ndimA, nout, nout, m1,
		 ndimB, nout, nin, m2,
		 ndimB, nout, nin, dest);
      }
	 
      if (C != NULL && pa > 0) {
	 for (integer j = 1; j <= pa-1; j++) {
	    /* c_j^l = a_1^(l-1)*c_j^0+c_(j+1)^(l-1) */
	    doublereal* dest = C+nout*(i-1)+ndimC*(nout*(j-1)); /* c_j^l */
	    doublereal* m2 = C+nout*(s-1)+ndimC*(nout*(j-1));   /* c_j^0 */
	    doublereal* m3 = C+nout*i+ndimC*(nout*j);           /* c_(j+1)^(l-1) */
	    gpc_addmul(ndimA, nout, nout, m1,
		       ndimC, nout, nout, m2,
		       ndimC, nout, nout, m3,
		       ndimC, nout, nout, dest);
	 }

	 /* c_pa^l = a_1^(l-1)*c_pa^0 */
	 doublereal* dest = C+nout*(i-1)+ndimC*(nout*(pa-1)); /* c_pa^l */
	 doublereal* m2 = C+nout*(s-1)+ndimC*(nout*(pa-1));   /* c_pa^0 */
	 gpc_mul(ndimA, nout, nout, m1,
		 ndimC, nout, nout, m2,
		 ndimC, nout, nout, dest);
      }

      do {
	 do {
	    /* p_l^l .. p_(s-1)^l = p_(l+1)^(l-1) .. p_s^(l-1) */
	    doublereal* source = P+nout*i+ndimP*(nin*i);       /* p_(l+1)^(l-1) .. p_s^(l-1) */
	    doublereal* dest = P+nout*(i-1)+ndimP*(nin*(i-1)); /* p_l^l .. p_(s-1)^l */
	    gpc_mcopy(ndimP, nout, nin*(s-i), source,
		      ndimP, dest);
	 } while (0);

	 /* p_s^l = a_1^(l-1)*p_s^0+b_1^(l-1) */
	 doublereal* dest = P+nout*(i-1)+ndimP*(nin*(s-1)); /* p_s^l */
	 doublereal* m2 = P+nout*(s-1)+ndimP*(nin*(s-1));   /* p_s^0 */
	 doublereal* m3 = B+nout*i;                         /* b_1^(l-1) */
	 gpc_addmul(ndimA, nout, nout, m1,
		    ndimP, nout, nin, m2,
		    ndimB, nout, nin, m3,
		    ndimP, nout, nin, dest);
      } while (0);
   }

   return 0;
}


extern "C" {
   extern int __FC_DECL__(dgesvd)(char* jobu, char* jobvt, 
				  integer* m, integer* n, 
				  doublereal* a, integer* lda, 
				  doublereal* s, 
				  doublereal* u, integer* ldu, 
				  doublereal* vt, integer* ldvt,
				  doublereal* work, integer* lwork, 
				  integer* info);
}


/*
 * esegue la pseudoinversa di una matrice, 
 * sovrascrivendola orientata per righe
 */
 
integer gpc_pinv(integer lda, integer m, integer n, doublereal* a,
		 doublereal* s, 
		 integer ldu, doublereal* u,
		 integer ldvt, doublereal* vt,
		 integer lwork, doublereal* work)
{
#ifdef USE_LAPACK
   /* esegue la SVD di a, m x n, con routines di LAPACK;
    * ritorna la pseudo-inversa in a, con dimensioni invertite: n x m
    * la matrice pinv(a) viene organizzata per righe */
   static char jobu = 'S';
   static char jobvt = 'S';

   __FC_DECL__(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
		       work, &lwork, &info);
#else
   cerr << "warning: LAPACK libraries are not available," << endl
     << "so the pseudo-inversion will not be performed" << endl;
#endif

   integer info = 0;

   if (info != 0) {
      return info;
   }   
   
   /* corregge con gli inversi dei SV */
   const doublereal THR = s[0]*1.e-12; /* mettere qualcosa tipo RCOND */
   integer x = min(m, n);
   for (integer i = x; i-- > 0; ) {   
      doublereal d = s[i];
      if (fabs(d) < THR) {
	 d = 0.;
      } else {
	 d = 1./d;
      }

      doublereal* p = u+i*ldu;
      for (integer k = m; k-- > 0; ) {
	 p[k] *= d;
      }
   }

   /* mette l'inversa in a (organizzata per righe) */
   for (integer i = n; i-- > 0; ) {
      doublereal* aa = a+m*i;
      doublereal* vvt = vt+i*ldvt;
      for (integer j = m; j-- > 0; ) {
	 doublereal* uu = u+j;
	 aa[j] = 0.;
	 for (integer k = x; k-- > 0; ) {
	    aa[j] += vvt[k]*uu[k*ldu];
	 }
      }
   }

   return info;
}


/* GPCInv - begin */

GPCInv::GPCInv(void) : pdBase(NULL) 
{
   NO_OP;
}
   
GPCInv::~GPCInv(void) 
{
   if (pdBase != NULL) {
      SAFEDELETEARR(pdBase, DMmm);
   }
}
   
/* GPCInv - end */



/* GPC_LAPACK_pinv - begin */

GPC_LAPACK_pinv::GPC_LAPACK_pinv(integer m, integer n)
: GPCInv(), m(m), n(n), iMin(0), iMax(0), iWork(0),
pdS(NULL), pdU(NULL), pdVt(NULL), pdWork(NULL)
{
   ASSERT(m > 0);
   ASSERT(n > 0);   
   
   iMin = min(m, n);                /* usati per vari parametri */
   iMax = max(m, n);                /* serve solo qui, non vale la pena tenerlo */
   
   iWork = max(3*iMin+iMax, 5*iMin-4); /* aumentare per maggiore efficienza */
   
   /* aree di lavoro per SVD:
    * s:    (min(m,n))
    * u:    (ldu*min(m,n))
    *       ldu >= m
    * vt:   (ldvt*n)
    *       ldvt >= min(m,n)
    * work: lwork
    *       lwork >= max(3*min(m,n)+max(m,n),5*min(m,n)-4)
    */
     
   integer i = iMin                           /* s */
     +m*iMin                                  /* u */
     +iMin*n                                  /* vt */
#ifndef OPTIMIZE_WORK
     +iWork;                                  /* work */
#else
     ;

   /* per ottimizzare le dimensioni di work: viene allocato a parte */
   SAFENEWARR(pdWork, doublereal, iWork, DMmm);
#endif
     
   SAFENEWARR(pdBase, doublereal, i, DMmm);
   
   pdS = pdBase;
   pdU = pdS+iMin;
   pdVt = pdU+m*iMin;
   
#ifndef OPTIMIZE_WORK
   pdWork = pdVt+iMin*n;
#endif
}


GPC_LAPACK_pinv::~GPC_LAPACK_pinv(void)
{
#ifndef OPTIMIZE_WORK
   NO_OP;                 /* pdBase viene distrutto da GPCInv */
#else
   SAFEDELETEARR(pdWork, DMmm);
#endif
}


integer GPC_LAPACK_pinv::Inv(integer ndima, integer nrowa, 
			  integer ncola, doublereal* a)
{
   ASSERT(nrowa == m);
   ASSERT(ncola == n);
   
   integer info = gpc_pinv(ndima, nrowa, ncola, a,
			   pdS, m, pdU, iMin, pdVt, iWork, pdWork);
   
   if (info != 0) {
      return info;
   }

#ifdef OPTIMIZE_WORK
   if (((integer*)pdWork)[0] > iWork) {
      iWork = ((integer*)pdWork)[0];
      SAFEDELETEARR(pdWork, DMmm);
      SAFENEWARR(pdWork, doublereal, iWork, DMmm);
   }
#endif
   
   return integer(0);
}

/* GPC_LAPACK_pinv - end */


#ifdef USE_MESCHACH 
/* GPC_Meschach_QRinv - begin */
GPC_Meschach_QRinv::GPC_Meschach_QRinv(integer m, integer n)
: m(m), n(n), A(MNULL), diag(VNULL), x(VNULL), b(VNULL), pivot(PNULL)
{
   ASSERT(m > 0);
   ASSERT(n > 0);
   if ((A = m_get(m, n)) == MNULL) {
      cerr << "A = m_get(" << m << ',' << n << ") failed" << endl;
      THROW(ErrGeneric());
   }
   if ((diag = v_get(min(m, n))) == VNULL) {
      cerr << "diag = v_get(" << min(m, n) << ") failed" << endl;
      THROW(ErrGeneric());
   }
   if ((x = v_get(n)) == VNULL) {
      cerr << "x = v_get(" << n << ") failed" << endl;
      THROW(ErrGeneric());
   }
   if ((b = v_get(m)) == VNULL) {
      cerr << "b = v_get(" << m << ") failed" << endl;
      THROW(ErrGeneric());
   }
   if ((pivot = px_get(n)) == PNULL) {
      cerr << "pivot = px_get(" << n << ") failed" << endl;
      THROW(ErrGeneric());
   }
}


GPC_Meschach_QRinv::~GPC_Meschach_QRinv(void)
{
   PX_FREE(pivot);
   V_FREE(b);
   V_FREE(x);
   V_FREE(diag);
   M_FREE(A);
}


integer GPC_Meschach_QRinv::Inv(integer ndima, 
				integer nrowa, 
				integer ncola, 
				doublereal* a)
{
   ASSERT(nrowa == m);
   ASSERT(ncola == n);
   ASSERT(ndima >= m);
   ASSERT(a != NULL);
   
   /* copia a in A */
   for (int j = n; j-- > 0; ) {
      doublereal* p = a+ndima*j;
      for (int i = m; i-- > 0; ) {
	 A->me[i][j] = p[i];
      }
   }
   
   QRCPfactor(A, diag, pivot);
   
   /* puo' essere ottimizzata con accesso diretto alla fattorizzata */
   for (int j = m; j-- > 0; ) {
      v_zero(b);
      b->ve[j] = 1.;
      QRCPsolve(A, diag, pivot, b, x);
      doublereal* p = a+j;
      for (int i = n; i-- > 0; ) {
	 p[i*m] = x->ve[i];
      }
   }
   
   return 0;
}
 
/* GPC_Meschach_QRinv - end */
#endif /* USE_MESCHACH */
 

/* GPCDesigner - begin */

GPCDesigner::GPCDesigner(integer iNumOut, integer iNumIn, 
			 integer iOrdA, integer iOrdB,
			 integer iPredS, integer iContrS,
			 integer iPredH, integer iContrH)
: iNumOutputs(iNumOut), iNumInputs(iNumIn),
iOrderA(iOrdA), iOrderB(iOrdB),
iPredStep(iPredS), iContrStep(iContrS),
iPredHor(iPredH), iContrHor(iContrH),
pdBase(NULL), pdA(NULL), pdB(NULL), pdP(NULL),pdC(NULL),
pdac(NULL), pdbc(NULL), pdmd(NULL), pdcc(NULL)
{
   ASSERT(iNumOutputs > 0);
   ASSERT(iNumInputs > 0);
   ASSERT(iPredStep > 0);
   ASSERT(iContrStep > 0);
}


GPCDesigner::~GPCDesigner(void)
{
   SAFEDELETEARR(pdBase, DMmm);
}


void GPCDesigner::DesignControl(const doublereal* /* pdTheta */ ,
				doublereal** ppda, 
				doublereal** ppdb,
				doublereal** ppdm, 
				doublereal** ppdc) 
{
   if (ppda != NULL) {
      *ppda = pdac;
   }
   if (ppdb != NULL) {
      *ppdb = pdbc;
   }
   if (ppdm != NULL) {
      *ppdm = pdmd;
   }
   if (ppdc != NULL) {
      *ppdc = pdcc;
   }
}

/* GPCDesigner - end */


/* DeadBeat - begin */

DeadBeat::DeadBeat(integer iNumOut, integer iNumIn,
		   integer iOrdA, integer iOrdB,
		   integer iPredS, integer iContrS, flag f)
: GPCDesigner(iNumOut, iNumIn, iOrdA, iOrdB, iPredS, iContrS, iContrS, 0),
iDim(iNumOutputs*iPredStep), 
iTmpRows(iNumOutputs*(iPredStep-iPredHor)),
iTmpCols(iNumInputs*(iContrStep-0)),
f_armax(f),
pInv(NULL)
{
#if !defined(USE_MESCHACH) && !defined(USE_LAPACK)
#error "need USE_LAPACK or USE_MESCHACH"
#endif
   switch (1) { /* mettere un flag */
#if defined(USE_MESCHACH)
    case 1:
      SAFENEWWITHCONSTRUCTOR(pInv, 
			     GPC_Meschach_QRinv, 
			     GPC_Meschach_QRinv(iTmpRows, iTmpCols),
			     DMmm);
      break;
#endif // USE_MESCHACH
#if defined(USE_LAPACK)
    case 2:
      SAFENEWWITHCONSTRUCTOR(pInv, 
			     GPC_LAPACK_pinv, 
			     GPC_LAPACK_pinv(iTmpRows, iTmpCols),
			     DMmm);
#endif // USE_LAPACK
   }
   
   /* note that iPredHor == iContrStep */
     
   integer i = iDim*(iNumOutputs*iOrderA)            /* A */
     +iDim*(iNumInputs*iOrderB)                      /* B */
     +iDim*(iNumInputs*iPredStep)                    /* P */
     +iNumInputs*(iNumOutputs*iOrderA)               /* ac */
     +iNumInputs*(iNumInputs*iOrderB)                /* bc */
     +iNumInputs*iTmpRows;                           /* md */

   if (f_armax) {
      i += iDim*(iNumOutputs*iOrderA)            /* C */
	+iNumInputs*(iNumOutputs*iOrderA);       /* cc */
   }

   SAFENEWARR(pdBase, doublereal, i, DMmm);
   
   for (integer j = i; j-- > 0; ) {
      pdBase[j] = 0.;
   }
   
   pdA = pdBase;
   pdB = pdA+iDim*(iNumOutputs*iOrderA);
   pdP = pdB+iDim*(iNumInputs*iOrderB);
   pdac = pdP+iDim*(iNumInputs*iPredStep);
   pdbc = pdac+iNumInputs*(iNumOutputs*iOrderA);
   pdmd = pdbc+iNumInputs*(iNumInputs*iOrderB);

   pdPTmp = pdP+iDim*(iNumInputs*(iPredStep-iContrStep));
   
   if (f_armax) {
      pdC = pdmd+iNumInputs*iTmpRows;
      pdcc = pdC+iDim*(iNumOutputs*iOrderA);
   }   
}


DeadBeat::~DeadBeat(void)
{
   NO_OP;
}


void DeadBeat::DesignControl(const doublereal* pdTheta,
			     doublereal** ppda,
			     doublereal** ppdb,
			     doublereal** ppdm,
			     doublereal** ppdc)
{
   ASSERT(pdTheta != NULL);
   
   /* si assume che pdTheta contenga le stime delle matrici a_1 .. a_pa,
    * b_0 .. b_pb disposte per righe, ovvero:
    * a_1(1,1), a_1(1,2), ..., a_1(1,m), a_2(1,1), ..., a_pa(1,m),b_0(1,1),...  */
   
   gpc_build_matrices(iDim, pdA,
		      iDim, pdB,
		      iDim, pdP,
		      iDim, pdC,  /* se f_armax == 0, pdC = NULL */
		      iNumOutputs, iNumInputs,
		      iPredStep, iOrderA, iOrderB,
		      (doublereal*)pdTheta);

   /* l'inversore sa come invertire, e mette l'inversa in pdPTmp */
   integer info = pInv->Inv(iDim, iTmpRows, iTmpCols, pdPTmp);
   
   if (info < 0) {
      cerr << "DeadBeat::DesignControl(): illegal value in " 
	<< -info << "-th argument of dgesvd()" << endl;
      THROW(ErrGeneric());
   } else if (info > 0) {
      cerr << "DeadBeat::DesignControl(): error in dgesvd()" << endl;
      THROW(ErrGeneric());
   } /* else: OK */
     
   /* recupero delle matrici (pinv(P) e' organizzata per righe) */

   /* copia le ultime righe dell'inversa di p, che sta in pdPTmp, in md */
   doublereal* p = pdPTmp+iTmpRows*(iNumInputs*(iContrStep-1));
   for (integer i = iTmpRows*iNumInputs; i-- > 0; ) {
      pdmd[i] = p[i];
   }

   /* anche ac, bc (e cc) vengono organizzate per righe */
   for (integer i = iNumInputs; i-- > 0; ) {
      doublereal* pm = pdmd+iTmpRows*i;
      
      /* ac = -md*A */
      doublereal* p = pdac+i*iNumOutputs*iOrderA;
      for (integer j = iNumOutputs*iOrderA; j-- > 0; ) {
	 doublereal* pa = pdA+iDim*j;
	 p[j] = 0.;
	 for (integer k = iTmpRows; k-- > 0; ) {
	    /* p[j] -= pdmd[iTmpRows*i+k]*pdA[k+iDim*j]; */
	    p[j] -= pm[k]*pa[k];
	 }	 
      }
      
      /* bc = -md*B */
      p = pdbc+i*iNumInputs*iOrderB;     
      for (integer j = iNumInputs*iOrderB; j-- > 0; ) {
	 doublereal* pb = pdB+iDim*j;
	 p[j] = 0.;
	 for (integer k = iTmpRows; k-- > 0; ) {
	    /* p[j] -= pdmd[iTmpRows*i+k]*pdB[k+iDim*j]; */
	    p[j] -= pm[k]*pb[k];
	 }	 
      }

      if (f_armax) {
	 /* cc = -md*C */
	 p = pdcc+i*iNumOutputs*iOrderA;     
	 for (integer j = iNumOutputs*iOrderA; j-- > 0; ) {
	    doublereal* pc = pdC+iDim*j;
	    p[j] = 0.;
	    for (integer k = iTmpRows; k-- > 0; ) {
	       /* p[j] -= pdmd[iTmpRows*i+k]*pdC[k+iDim*j]; */
	       p[j] -= pm[k]*pc[k];
	    }
	 }
      }
   }    
   
   /* linka i puntatori alle matrici, se richiesto */
   GPCDesigner::DesignControl(NULL, ppda, ppdb, ppdm, ppdc);
}

/* DeadBeat - end */
 

/* GPC - begin */

GPC::GPC(integer iNumOut, integer iNumIn,
	 integer iOrdA, integer iOrdB,
	 integer iPredS, integer iContrS, integer iPredH, integer iContrH,
	 doublereal* pW, doublereal* pR, DriveCaller* pDC,
	 flag f)
: GPCDesigner(iNumOut, iNumIn, iOrdA, iOrdB, iPredS, iContrS, iPredH, iContrH),
iDim(iNumOutputs*iPredStep), 
iTmpRows(iNumOutputs*(iPredStep-iPredHor)),
iTmpCols(iNumInputs*iContrStep),
pdW(pW),
pdR(pR),
Weight(pDC),
pdM(NULL),
pdInvP(NULL),
f_armax(f),
pInv(NULL)
{
   ASSERT(pW != NULL);
   ASSERT(pW != NULL);
   ASSERT(pDC != NULL);
   
#if !defined(USE_MESCHACH) && !defined(USE_LAPACK)
#error "need USE_LAPACK or USE_MESCHACH"
#endif
   switch (1) { /* mettere un flag */
#if defined(USE_MESCHACH)
    case 1:
      SAFENEWWITHCONSTRUCTOR(pInv, 
			     GPC_Meschach_QRinv, 
			     GPC_Meschach_QRinv(iTmpCols, iTmpCols),
			     DMmm);
      break;
#endif // USE_MESCHACH
#if defined(USE_LAPACK)
    case 2:
      SAFENEWWITHCONSTRUCTOR(pInv, 
			     GPC_LAPACK_pinv, 
			     GPC_LAPACK_pinv(iTmpCols, iTmpCols),
			     DMmm);
#endif // USE_LAPACK
   }
      
   if (pInv == NULL) {
      cerr << "unable to create GPCInv" << endl;
      THROW(ErrGeneric());
   }   
   
   /* note that iPredHor == iContrStep */
     
   integer i = iDim*(iNumOutputs*iOrderA)            /* A */
     +iDim*(iNumInputs*iOrderB)                      /* B */
     +iDim*(iNumInputs*iPredStep)                    /* P */
     +iTmpCols*iTmpCols                              /* M */
     +iTmpCols*iTmpRows                              /* InvP */
     +iNumInputs*(iNumOutputs*iOrderA)               /* ac */
     +iNumInputs*(iNumInputs*iOrderB)                /* bc */
     +iNumInputs*iTmpRows;                           /* md */

   if (f_armax) {
      i += iDim*(iNumOutputs*iOrderA)            /* C */
	+iNumInputs*(iNumOutputs*iOrderA);       /* cc */
   }

   SAFENEWARR(pdBase, doublereal, i, DMmm);
   
   for (integer j = i; j-- > 0; ) {
      pdBase[j] = 0.;
   }
   
   pdA = pdBase;
   pdB = pdA+iDim*(iNumOutputs*iOrderA);
   pdP = pdB+iDim*(iNumInputs*iOrderB);   
   pdM = pdP+iDim*(iNumInputs*iPredStep);
   pdInvP = pdM+iTmpCols*iTmpCols;
   pdac = pdInvP+iTmpCols*iTmpRows;
   pdbc = pdac+iNumInputs*(iNumOutputs*iOrderA);
   pdmd = pdbc+iNumInputs*(iNumInputs*iOrderB);

   pdPTmp = pdP+iDim*(iNumInputs*(iPredStep-iContrStep));
   
   if (f_armax) {
      pdC = pdmd+iNumInputs*iTmpRows;
      pdcc = pdC+iDim*(iNumOutputs*iOrderA);
   }   
}


GPC::~GPC(void)
{
   SAFEDELETEARR(pdW, DMmm);
   SAFEDELETEARR(pdR, DMmm);
}


int make_m(integer ndimp, integer nrowp, integer ncolp, doublereal* P,
	   integer nout, integer nin, integer s,
	   doublereal* W, doublereal* R, doublereal lambda, doublereal* M)
{  
   
   for (integer j = ncolp; j-- > 0; ) {
      doublereal* m = M+ncolp*j;
      doublereal* pj = P+ndimp*j;
      for (integer i = ncolp; i-- > 0; ) {
	 doublereal* mm = m+i;
	 doublereal* pi = P+ndimp*i;	 
	 mm[0] = 0.;
	 for (integer k = s; k-- > 0; ) {
	    doublereal* pik = pi+nout*k;
	    doublereal* pjk = pj+nout*k;
	    for (integer l = nout; l-- > 0; ) {
	       /* 
		* M[i+ncolp*j] = P[nout*k+l+ndimp*i]*P[nout*k+l+ndimp*j]*W[k];
		*/
	       mm[0] += pik[l]*pjk[l]*W[k];
	    }
	 }
      }
      /*
       * M[j+ncolp*j] += lambda*R[j/nin];
       */
      m[j] += lambda*R[j/nin];
   }
   
   return 0;
}


/*
 * moltiplica una matrice L, organizzata per righe, per una R, organizzata
 * per colonne, e trasposta, e mette il risultato in dest, organizzata per
 * righe 
 */

int gpc_mult_t(integer ndiml, integer nrowl, integer ncoll, doublereal* L,
	       integer ndimr, integer nrowr, integer ncolr, doublereal* R,
	       integer ndimd, doublereal* dest)
{
   for (integer i = nrowl; i-- > 0; ) {
      doublereal* di = dest+ndimd*i;
      doublereal* Li = L+ndiml*i;
      for (integer j = ncolr; j-- > 0; ) {
	 doublereal* dij = di+j;
	 doublereal* Rj = R+j;
	 /*
	  * dest[j+ndimd*i] = 0.;
	  */
	 dij[0] = 0.;
	 for (integer k = ncoll; k-- > 0; ) {
	    /*
	     * dest[j+ndimd*i] += L[k+ndiml*i]*R[j+ndimr*k];
	     */
	    dij[0] += Li[k]*Rj[ndimr*k];
	 }
      }
   }
   
   return 0;
}



void GPC::DesignControl(const doublereal* pdTheta,
			doublereal** ppda,
			doublereal** ppdb,
			doublereal** ppdm,
			doublereal** ppdc)
{
   ASSERT(pdTheta != NULL);
   
   /* si assume che pdTheta contenga le stime delle matrici a_1 .. a_pa,
    * b_0 .. b_pb disposte per righe, ovvero:
    * a_1(1,1), a_1(1,2), ..., a_1(1,m), a_2(1,1), ..., a_pa(1,m),b_0(1,1),... */
   
   gpc_build_matrices(iDim, pdA,
		      iDim, pdB,
		      iDim, pdP,
		      iDim, pdC,  /* se f_armax == 0, pdC = NULL */
		      iNumOutputs, iNumInputs,
		      iPredStep, iOrderA, iOrderB,
		      (doublereal*)pdTheta);

   
   /* M = P^T*W*P+R */
   make_m(iDim, iTmpRows, iTmpCols, pdPTmp,
	  iNumOutputs, iNumInputs, iPredStep-iPredHor,
	  pdW, pdR, Weight.dGet(), pdM);
   
   /* l'inversore sa come invertire, e mette l'inversa in pdM */
   integer info = pInv->Inv(iTmpCols, iTmpCols, iTmpCols, pdM);
   
   /* M^-1*P^T */
   gpc_mult_t(iTmpCols, iTmpCols, iTmpCols, pdM,
	      iDim, iTmpRows, iTmpCols, pdPTmp, 
	      iTmpRows, pdInvP);
   
   if (info < 0) {
      cerr << "GPC::DesignControl(): illegal value in " 
	<< -info << "-th argument of dgesvd()" << endl;
      THROW(ErrGeneric());
   } else if (info > 0) {
      cerr << "GPC::DesignControl(): error in dgesvd()" << endl;
      THROW(ErrGeneric());
   } /* else: OK */
     
   /* recupero delle matrici (pinv(P) e' organizzata per righe) */

   /* copia le ultime righe dell'inversa di p, che sta in pdPTmp, in md */
   doublereal* p = pdInvP+iTmpRows*(iNumInputs*(iContrStep-1));
   
#if 0
   for (integer i = iTmpRows*iNumInputs; i-- > 0; ) {
      pdmd[i] = p[i];
   }
#else
   for (integer i = iNumInputs; i-- > 0; ) {
      doublereal* pi = p+iTmpRows*i;
      doublereal* pmi = pdmd+iTmpRows*i;
      for (integer k = iPredStep-iPredHor; k-- > 0; ) {
	 doublereal* pik = pi+iNumOutputs*k;
	 doublereal* pmik = pmi+iNumOutputs*k;
	 for (integer l = iNumOutputs; l-- > 0; ) {
	    /*
	     * pdmd[iTmpRows*i+iNumOutputs*k+l] =
	     *   p[iTmpRows*i+iNumOutputs*k+l]*pdW[k];
	     */
	    pmik[l] = pik[l]*pdW[k];
	 }
      }
   }
#endif
   
   /* anche ac, bc (e cc) vengono organizzate per righe */
   for (integer i = iNumInputs; i-- > 0; ) {
      doublereal* pm = pdmd+iTmpRows*i;
      
      /* ac = -md*A */
      doublereal* p = pdac+i*iNumOutputs*iOrderA;
      for (integer j = iNumOutputs*iOrderA; j-- > 0; ) {
	 doublereal* pa = pdA+iDim*j;
	 p[j] = 0.;
	 for (integer k = iTmpRows; k-- > 0; ) {
	    /* p[j] -= pdmd[iTmpRows*i+k]*pdA[k+iDim*j]; */
	    p[j] -= pm[k]*pa[k];
	 }	 
      }
      
      /* bc = -md*B */
      p = pdbc+i*iNumInputs*iOrderB;     
      for (integer j = iNumInputs*iOrderB; j-- > 0; ) {
	 doublereal* pb = pdB+iDim*j;
	 p[j] = 0.;
	 for (integer k = iTmpRows; k-- > 0; ) {
	    /* p[j] -= pdmd[iTmpRows*i+k]*pdB[k+iDim*j]; */
	    p[j] -= pm[k]*pb[k];
	 }	 
      }

      if (f_armax) {
	 /* cc = -md*C */
	 p = pdcc+i*iNumOutputs*iOrderA;     
	 for (integer j = iNumOutputs*iOrderA; j-- > 0; ) {
	    doublereal* pc = pdC+iDim*j;
	    p[j] = 0.;
	    for (integer k = iTmpRows; k-- > 0; ) {
	       /* p[j] -= pdmd[iTmpRows*i+k]*pdC[k+iDim*j]; */
	       p[j] -= pm[k]*pc[k];
	    }
	 }
      }
   }    
   
   /* linka i puntatori alle matrici, se richiesto */
   GPCDesigner::DesignControl(NULL, ppda, ppdb, ppdm, ppdc);
}

/* GPC - end */
 




#ifdef TEST_DPC 
int main(void)
{
   const integer nout = 2;
   const integer nin = 2;
   const integer s = 4;
   const integer pa = 2;
   const integer pb = 2;
   
   doublereal d[nout*nout*pa*s+nout*nin*pb*s+nout*s*nin*s+nout*nout*pa*s];
   doublereal theta[nout*nout*pa+nout*nin*(pb+1)+nout*nout*pa] = {
      2., 0.,-1., 0.,     1., 0.,-1., 0., 1., 0.,    -2., 0., 1., 0.,
      0., 2., 0.,-1.,     0., 1., 0.,-1., 0., 1.,     0.,-2., 0., 1.
      };
   
   /*
   for (integer i = 0; i < nout*nout*pa+nout*nin*(pb+1)+nout*nout*pa; i++) {
      theta[i] = 1.;
   }
    */
   
   gpc_build_matrices(nout*s, d,
		      nout*s, d+nout*nout*pa*s,
		      nout*s, d+nout*nout*pa*s+nout*nin*pb*s,
		      nout*s, d+nout*nout*pa*s+nout*nin*pb*s+nout*s*nin*s,
		      nout, nin,
		      s, pa, pb,
		      theta);
   
   doublereal* pA = d;
   cout << "A:" << endl;
   for (integer i = 0; i < nout*s; i++) {     
      for (integer j = 0; j < nout*pa; j++) {
	 cout << setw(8) << pA[i+nout*s*j];
      }
      cout << endl;
   }
   doublereal* pB = d+nout*nout*pa*s;
   cout << "B:" << endl;
   for (integer i = 0; i < nout*s; i++) {     
      for (integer j = 0; j < nin*pb; j++) {
	 cout << setw(8) << pB[i+nout*s*j];
      }
      cout << endl;
   }
   doublereal* pP = d+nout*nout*pa*s+nout*s*nin*pb;
   cout << "P:" << endl;
   for (integer i = 0; i < nout*s; i++) {
      for (integer j = 0; j < nin*s; j++) {
	 cout << setw(8) << pP[i+nout*s*j];
      }
      cout << endl;
   }
   doublereal* pC = d+nout*nout*pa*s+nout*s*nin*pb+nout*s*nin*s;
   cout << "C:" << endl;
   for (integer i = 0; i < nout*s; i++) {     
      for (integer j = 0; j < nout*pa; j++) {
	 cout << setw(8) << pC[i+nout*s*j];
      }
      cout << endl;
   }      

   integer q = s/2;
   integer tmprow = (s-q)*nout;
   integer tmpcol = q*nin;
   GPC_Meschach_QRinv inv(tmprow, tmpcol);
   doublereal* pT = pP+s*nout*(s-q)*nin;
   inv.Inv(s*nout, tmprow, tmpcol, pT);
   
   cout << "Inv(P'):" << endl;
   for (integer i = 0; i < tmpcol; i++) {
      for (integer j = 0; j < tmprow; j++) {
	 cout << setw(12) << pT[tmprow*i+j];
      }
      cout << endl;
   }
   
   DeadBeat db(nout, nin,
	       pa, pb,
	       s, q, 1);
   
   doublereal* pac;
   doublereal* pbc;
   doublereal* pcc;
   doublereal* pmd;
   db.DesignControl(theta, &pac, &pbc, &pmd, &pcc);
   
   cout << "Ac:" << endl;
   for (integer i = 0; i < nin; i++) {
      for (integer j = 0; j < nout*pa; j++) {
	 cout << setw(12) << pac[i*nout*pa+j];
      }
      cout << endl;
   }
   cout << "Bc:" << endl;
   for (integer i = 0; i < nin; i++) {
      for (integer j = 0; j < nin*pb; j++) {
	 cout << setw(12) << pbc[i*nin*pb+j];
      }
      cout << endl;
   }
   cout << "Cc:" << endl;
   for (integer i = 0; i < nin; i++) {
      for (integer j = 0; j < nout*pa; j++) {
	 cout << setw(12) << pcc[i*nout*pa+j];
      }
      cout << endl;
   }
   cout << "Md:" << endl;
   for (integer i = 0; i < nin; i++) {
      for (integer j = 0; j < tmprow; j++) {
	 cout << setw(12) << pmd[i*tmprow+j];
      }
      cout << endl;
   }

   return 0;
}
#endif
