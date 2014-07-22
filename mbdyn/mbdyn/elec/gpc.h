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

#ifndef GPC_H
#define GPC_H

#include <ac/f2c.h>
#include <drive.h>

#ifdef USE_MESCHACH
extern "C" {
#include <meschach/matrix2.h>
}
#endif /* USE_MESCHACH */


/*
 * classe di routines che consentono di invertire matrici;
 * allocano autonomamente lo spazio di lavoro lo gestiscono e lo distruggono;
 * usate per l'inversione delle matrici nel progetto dei controllori
 */

/* GPCInv - begin */

class GPCInv {
protected:
   	doublereal* pdBase;   

public:
   	GPCInv(void);
   	virtual ~GPCInv(void);
   
   	virtual integer Inv(integer ndima, integer nrowa,
		            integer ncola, doublereal* a) = 0;
};

/* GPCInv - end */


/* 
 * si appoggia sulla routine LAPACK dgesvd() che esegue la decomposizione SVD;
 * quindi la pseudoinversa viene ricostruita in loco, orientata per righe 
 */

/* GPC_LAPACK_pinv - begin */

class GPC_LAPACK_pinv : public GPCInv {
protected:
   	integer m; /* in realta' non servono, e' solo per sicurezza */
   	integer n;
   
   	integer iMin;
   	integer iMax;

   	integer iWork; /* per ridimensionare pdWork */

   	doublereal* pdS;
   	doublereal* pdU;
   	doublereal* pdVt;
   	doublereal* pdWork;
   
public:
   	GPC_LAPACK_pinv(integer n, integer m);
   	~GPC_LAPACK_pinv(void);
   	integer Inv(integer ndima, integer nrowa, integer ncola, doublereal* a);
};

/* GPC_LAPACK_pinv - end */


#ifdef USE_MESCHACH
/* GPC_Meschach_QRinv - begin */

class GPC_Meschach_QRinv : public GPCInv {
protected:
   	integer m; /* in realta' non servono, e' solo per sicurezza */
   	integer n;
   
   	MAT* A;
   	VEC* diag;
   	VEC* x;
   	VEC* b;
   	PERM* pivot;
   
public:
   	GPC_Meschach_QRinv(integer m, integer n);
   	~GPC_Meschach_QRinv(void);
   	integer Inv(integer ndima, integer nrowa, integer ncola, doublereal* a);
};

/* GPC_Meschach_QRinv - end */
#endif /* USE_MESCHACH */


/*
 * Progetta le matrici di controllo per un controllore discreto MIMO
 * a partire dal sistema identificato
 */

/* GPCDesigner - begin */

class GPCDesigner {
protected:
   	integer iNumOutputs;
   	integer iNumInputs;
   	integer iOrderA;
   	integer iOrderB;

   	integer iPredStep;      /* passo della predizione piu' in avanti */
   	integer iContrStep;     /* passo del controllo piu' avanti */
   	integer iPredHor;       /* passo iniziale della predizione */
   	integer iContrHor;      /* passo iniziale del controllo (sempre 0) */

	doublereal dPeriodicFactor;

   	doublereal* pdBase;

   	doublereal* pdA;
   	doublereal* pdB;
   	doublereal* pdP;
   	doublereal* pdC;
   
   	doublereal* pdac;
   	doublereal* pdbc;
   	doublereal* pdmd;
   	doublereal* pdcc;

public:
   	GPCDesigner(integer iNumOut, integer iNumIn, 
		    integer iOrdA, integer iOrdB,
		    integer iPredS, integer iContrS, 
		    integer iPredH, integer iContrH,
		    doublereal dPF);
   
   	virtual ~GPCDesigner(void);
   
   	virtual void DesignControl(const doublereal* /* pdTheta */ ,
			      	   doublereal** ppda = NULL, 
				   doublereal** ppdb = NULL,
				   doublereal** ppdm = NULL, 
				   doublereal** ppdc = NULL);
   
   	inline doublereal* pdGetAc(void) const;
   	inline doublereal* pdGetBc(void) const;
   	inline doublereal* pdGetMd(void) const;
   	inline doublereal* pdGetCc(void) const;

   	inline integer iGetPredStep(void) const;
   	inline integer iGetContrStep(void) const;
   	inline integer iGetPredHor(void) const;
   	inline integer iGetContrHor(void) const;
};

inline doublereal* 
GPCDesigner::pdGetAc(void) const 
{
   	return pdac;
}

inline doublereal* 
GPCDesigner::pdGetBc(void) const 
{
   	return pdbc;
}

inline doublereal* 
GPCDesigner::pdGetMd(void) const 
{
   	return pdmd;
}

inline doublereal* 
GPCDesigner::pdGetCc(void) const 
{
   	return pdcc;
}

inline integer 
GPCDesigner::iGetPredStep(void) const 
{
   	return iPredStep;
}

inline integer 
GPCDesigner::iGetContrStep(void) const 
{
   	return iContrStep;
}

inline integer 
GPCDesigner::iGetPredHor(void) const
{
   	return iPredHor;
}

inline integer 
GPCDesigner::iGetContrHor(void) const
{
   	return iContrHor;
}

/* GPCControlDesigner - end */


/*
 * progetta un controllore deadbeat basato su un sistema ARX o ARMAX
 */
 
/* DeadBeat - begin */

#ifdef USE_DBC
class DeadBeat : public GPCDesigner{
protected:
   	integer iDim;
   	integer iTmpRows;
   	integer iTmpCols;
   
   	doublereal* pdPTmp;
   
   	flag f_armax;

   	GPCInv* pInv;
   
public:
   	DeadBeat(integer iNumOut, integer iNumIn, integer iOrdA, integer iOrdB,
	    	 integer iPredS, integer iContrS, doublereal dPF, flag f);
   	virtual ~DeadBeat(void);
   
   	void DesignControl(const doublereal* pdTheta,
		           doublereal** ppda = NULL, doublereal** ppdb = NULL,
		      	   doublereal** ppdm = NULL, doublereal** ppdc = NULL);
};
#endif /* USE_DBC */

/* DeadBeat - end */


/*
 * progetta un controllore GPC basato su un sistema ARX o ARMAX
 */
 
/* GPC - begin */

class GPC : public GPCDesigner{
protected:
   	integer iDim;
   	integer iTmpRows;
   	integer iTmpCols;
   
   	doublereal* pdPTmp;
   
   	doublereal* pdW;
   	doublereal* pdR;
   	DriveOwner Weight;
   
   	doublereal* pdM;
   	doublereal* pdInvP;
   
   	flag f_armax;

   	GPCInv* pInv;
   
public:
   	GPC(integer iNumOut, integer iNumIn, integer iOrdA, integer iOrdB,
       	    integer iPredS, integer iContrS, integer iPredH, integer iContrH,
	    doublereal* pW, doublereal* pR, DriveCaller* pDC, doublereal dPF,
	    flag f);
   
   	virtual ~GPC(void);
   
   	void DesignControl(const doublereal* pdTheta,
		           doublereal** ppda = NULL, doublereal** ppdb = NULL,
			   doublereal** ppdm = NULL, doublereal** ppdc = NULL);
};

/* GPC - end */

#endif /* GPC_H */

