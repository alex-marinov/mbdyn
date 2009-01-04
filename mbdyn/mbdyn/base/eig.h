/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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
 * Copyright 2003-2009 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 */
 
#ifndef PODEIG_H
#define PODEIG_H 

#include <myassert.h>
#include <mynewmem.h>
#include <solman.h>
#include <fullmh.h>

#include <vector>
#include <algorithm>
#include <cmath>

class PODEig
{

   public:
   	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

   private:

  	doublereal  T;			/* periodo calcolo */	  		
	integer     iSSize;		/* dimensione del vettore degli stati */
	
	
	integer     iSnapN;             /* numero di snapshot */      
	integer     iSnapCount; 
	
	doublereal  dTimeCount;
	doublereal  dStartTime;		/* istante iniziale per l'estazione degli snapshot */
	doublereal  dFinalTime;
	doublereal  dTresh;		/* soglia per la scelta del numero di modi da identificare */
	integer     iNumEig;		/* numero di autovalori da estrarre con il metodo degli snapdshot
					   0 indica default = tutti */
	doublereal  dPecision; 		/* precisione da utilizzare nell'estrazione degli autovalori 
					   0 indica default = 10^-6 */				   
	
	std::vector<doublereal>*  pXMat;      /* matrice che contiene gli 
					  [snapshot numero di stati] x [numero di snapshot]
					  ordinati per righe in modo da poter normalizzare facilmente le storie 
					  temporali */

	bool bEigDone;		        /* bool per indicare che il calcolo degli autovalori e' 
					   avvenuto */
    public:
    
    	 PODEig(doublereal Period, 
	        integer S, 
		integer N, 
		doublereal Start,
		doublereal Tr);
		
	~PODEig();
	
	void LogData(doublereal t, VectorHandler& SV);		
    
        void ComputeEigenvalues(doublereal t);
		 
        void OutputEigenvalues(void);
};

#endif /* PODEIG_H */
