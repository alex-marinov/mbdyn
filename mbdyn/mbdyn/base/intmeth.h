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

#ifndef INTMETH_H
#define INTMETH_H

#include "myassert.h"
#include "dofown.h"
#include "drive.h"

/* DerivativePrediction - begin

class DerivativePrediction {
 public:
   virtual DerivativePrediction(void) {
      NO_OP;
   };
   
   virtual doublereal 
     dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
   
   virtual inline doublereal 
     dPredDer(const doublereal& dXm1,
	      const doublereal& dXm2,
	      const doublereal& dXPm1,
	      const doublereal& dXPm2) const = 0;
   
   virtual inline doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const = 0;  
};

 DerivativePrediction - end */


/* MultiStepIntegrationMethod - begin */

class MultiStepIntegrationMethod {
 public:
   enum { DIFFERENTIAL = 0, ALGEBRAIC = 1 };
   enum StepChange { NEWSTEP, REPEATSTEP };   
   
 public:
   virtual ~MultiStepIntegrationMethod(void) { 
      NO_OP; 
   };
   
   virtual void SetCoef(doublereal dT,
			doublereal dAlpha,
			enum StepChange NewStep,
			doublereal& db0Differential,
			doublereal& db0Algebraic) = 0;
   
   virtual void SetDriveHandler(const DriveHandler* pDH) = 0;
   
   // Overridden by dedicated inline functions
   virtual doublereal 
     dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
   
   // Overridden by dedicated inline functions
   virtual doublereal 
     dPredictState(const doublereal& dXm1,
		   const doublereal& dXm2,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
 
   virtual inline doublereal 
     dPredDer(const doublereal& dXm1,
	      const doublereal& dXm2,
	      const doublereal& dXPm1,
	      const doublereal& dXPm2) const = 0;
   
   virtual inline doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const = 0;
   
   virtual inline doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const = 0;
   
   virtual inline doublereal 
     dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2) const = 0;
};

/* MultiStepIntegrationMethod - end */


/* CrankNicholson - begin */

class CrankNicholson : public MultiStepIntegrationMethod {
 protected:
   doublereal db0;
   
 public:
   CrankNicholson(void);
   
   virtual ~CrankNicholson(void) { 
      NO_OP; 
   };
   
   virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep,
			doublereal& db0Differential,
			doublereal& db0Algebraic);
   
   virtual void SetDriveHandler(const DriveHandler* /* pDH */ ) {
      NO_OP;
   };
   
   virtual doublereal 
     dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
   virtual doublereal 
     dPredictState(const doublereal& dXm1,
		   const doublereal& dXm2,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
   // Nota: usa predizione lineare per le derivate (massimo ordine possibile)
   virtual inline doublereal 
     dPredDer(const doublereal& /* dXm1 */ ,
	      const doublereal& /* dXm2 */ ,
	      const doublereal& dXPm1,
	      const doublereal& /* dXPm2 */ ) const {
		 return dXPm1;
	      };
   
   virtual inline doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& /* dXm2 */ ,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& /* dXPm2 */ ) const {
		   return dXm1+db0*(dXP+dXPm1);
		};
   
   virtual inline doublereal 
     dPredDerAlg(const doublereal& /* dXm1 */ ,
		 const doublereal& dXPm1,
		 const doublereal& /* dXPm2 */ ) const {
		    return dXPm1;
		 };
   
   virtual inline doublereal 
     dPredStateAlg(const doublereal& /* dXm1 */ ,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& /* dXPm2 */ ) const {
		      return db0*(dXP+dXPm1);
		   };      
};

/* CrankNicholson - end */


/* NostroMetodo - begin */

class NostroMetodo : public MultiStepIntegrationMethod {
 protected:
   DriveOwner Rho;
   DriveOwner AlgebraicRho;
   
   doublereal a[2][2];
   doublereal b[3][2];
   
   doublereal mp[2];
   doublereal np[2];
   
 public:
   NostroMetodo(const DriveCaller* pRho, const DriveCaller* pAlgRho);
   
   virtual ~NostroMetodo(void) { 
      NO_OP; 
   };
   
   virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep,
			doublereal& db0Differential,
			doublereal& db0Algebraic);
   
   virtual void SetDriveHandler(const DriveHandler* pDH) {
      Rho.pGetDriveCaller()->SetDrvHdl(pDH);
      AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
   };
   
   
   virtual doublereal 
     dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
   virtual doublereal 
     dPredictState(const doublereal& dXm1,
		   const doublereal& dXm2,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
   // Nota: usa predizione cubica per le derivate (massimo ordine possibile)
   virtual inline doublereal 
     dPredDer(const doublereal& dXm1,
	      const doublereal& dXm2,
	      const doublereal& dXPm1,
	      const doublereal& dXPm2) const {
		 return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
	      };
   
   virtual inline doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const {
		   return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
		     +b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1+b[2][DIFFERENTIAL]*dXPm2;
		};
   
   virtual inline doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const {
		    return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		 };
   
   virtual inline doublereal 
     dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2) const {
		      return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2
			-a[1][ALGEBRAIC]*dXm1;
		   };
};

/* NostroMetodo - end */


/* Hope - begin */

class Hope : public MultiStepIntegrationMethod {
 protected:
   DriveOwner Rho;
   DriveOwner AlgebraicRho;
   
   flag fStep;
   
   doublereal a[2][2];
   doublereal b[2][2];
   
   doublereal mp[2];
   doublereal np[2];
   
 public:
   Hope(const DriveCaller* pRho, const DriveCaller* pAlgRho);
   
   virtual ~Hope(void) { 
      NO_OP; 
   };
   
   virtual void SetCoef(doublereal dT,
			doublereal dAlpha,
			enum StepChange NewStep,
			doublereal& db0Differential,
			doublereal& db0Algebraic);
   
   virtual void SetDriveHandler(const DriveHandler* pDH) {
      Rho.pGetDriveCaller()->SetDrvHdl(pDH);
      AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
   };

   virtual doublereal 
     dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
   virtual doublereal 
     dPredictState(const doublereal& dXm1,
		   const doublereal& dXm2,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const;      

   // Nota: usa predizione cubica per le derivate (massimo ordine possibile)
   virtual inline doublereal 
     dPredDer(const doublereal& dXm1,
	      const doublereal& dXm2,
	      const doublereal& dXPm1,
	      const doublereal& dXPm2) const {
		 return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
	      };
   
   virtual inline doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& /* dXPm2 */ ) const {
		   if (fStep) {
		      return dXm1+b[0][ALGEBRAIC]*(dXP+dXPm1);
		   } else {
		      return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
			+b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1;
		   }
		};
   
   virtual inline doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const {
		    return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		 };
   
   virtual inline doublereal 
     dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& /* dXPm2 */ ) const {
		      if (fStep) {
			 return b[0][ALGEBRAIC]*(dXP+dXPm1);
		      } else {
			 return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1
						   -a[1][ALGEBRAIC]*dXm1;
		      }
		   };
};

/* Hope - end */

#endif // INTMETH_H
