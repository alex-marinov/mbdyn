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

#include <mbconfig.h>

#include <intmeth.h>

/* MultiStepIntegrationMethods - begin */

/* CrankNicholson - begin */

CrankNicholson::CrankNicholson(void)
: db0(0.)
{
   NO_OP;
}


void CrankNicholson::SetCoef(doublereal dT,
			     doublereal dAlpha,
			     enum StepChange /* NewStep */ ,
			     doublereal& db0Differential,
			     doublereal& db0Algebraic)
{
   /* valori di ritorno */
   db0 = db0Differential = db0Algebraic = dT*dAlpha/2.;
}


doublereal CrankNicholson::dPredictDerivative(const doublereal& /* dXm1 */ ,
					      const doublereal& /* dXm2 */ ,
					      const doublereal& dXPm1,
					      const doublereal& /* dXPm2 */ ,
					      DofOrder::Order /* o */ ) const
{
   return dXPm1;
}


doublereal CrankNicholson::dPredictState(const doublereal& dXm1,
					 const doublereal& /* dXm2 */ ,
					 const doublereal& dXP,
					 const doublereal& dXPm1,
					 const doublereal& /* dXPm2 */ ,
					 DofOrder::Order o) const
{
   if (o == DofOrder::ALGEBRAIC) {
      return db0*(dXP+dXPm1);     
   } /* else if (o == DofOrder::DIFFERENTIAL) */   
   return dXm1+db0*(dXP+dXPm1);
}

/* CrankNicholson - end */


/* NostroMetodo - begin */

NostroMetodo::NostroMetodo(const DriveCaller* pRho,
			   const DriveCaller* pAlgRho)
: Rho(pRho), AlgebraicRho(pAlgRho)
{
   ASSERT(pRho != NULL);
   ASSERT(pAlgRho != NULL);
}


void NostroMetodo::SetCoef(doublereal dT,
			   doublereal dAlpha,
			   enum StepChange /* NewStep */ ,
			   doublereal& db0Differential,
			   doublereal& db0Algebraic)
{
   doublereal dRho = Rho.dGet();
   doublereal dAlgebraicRho = AlgebraicRho.dGet();
   
   doublereal dDen = 2.*(1.+dAlpha)-(1.-dRho)*(1.-dRho);
   doublereal dBeta = 
     dAlpha*((1.-dRho)*(1.-dRho)*(2.+dAlpha)
	     +2.*(2.*dRho-1.)*(1.+dAlpha))/dDen;
   doublereal dDelta = 
     .5*dAlpha*dAlpha*(1.-dRho)*(1.-dRho)/dDen;

   mp[0] = -6.*dAlpha*(1.+dAlpha);
   mp[1] = -mp[0];
   np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
   np[1] = dAlpha*(2.+3.*dAlpha);
   
   a[0][DIFFERENTIAL] = 1.-dBeta;
   a[1][DIFFERENTIAL] = dBeta;
   b[0][DIFFERENTIAL] = dT*(dDelta/dAlpha+dAlpha/2);
   b[1][DIFFERENTIAL] = dT*(dBeta/2.+dAlpha/2.-dDelta/dAlpha*(1.+dAlpha));
   b[2][DIFFERENTIAL] = dT*(dBeta/2.+dDelta);
   
   DEBUGCOUT("Predict()" << endl
	     << "Alpha = " << dAlpha << endl
	     << "Differential coefficients:" << endl
	     << "beta  = " << dBeta << endl
	     << "delta = " << dDelta << endl
	     << "a1    = " << a[0][DIFFERENTIAL] << endl
	     << "a2    = " << a[1][DIFFERENTIAL] << endl
	     << "b0    = " << b[0][DIFFERENTIAL] << endl
	     << "b1    = " << b[1][DIFFERENTIAL] << endl
	     << "b2    = " << b[2][DIFFERENTIAL] << endl);
   
   /* Coefficienti del metodo - variabili algebriche */
   if (dAlgebraicRho != dRho) {
      dDen = 2.*(1.+dAlpha)-(1.-dAlgebraicRho)*(1.-dAlgebraicRho);
      dBeta = dAlpha*((1.-dAlgebraicRho)*(1.-dAlgebraicRho)*(2.+dAlpha)
		      +2.*(2.*dAlgebraicRho-1.)*(1.+dAlpha))/dDen;      
      dDelta = .5*dAlpha*dAlpha*(1.-dAlgebraicRho)*(1.-dAlgebraicRho)/dDen;
            
      a[1][ALGEBRAIC] = dBeta;
      b[0][ALGEBRAIC] = dT*(dDelta/dAlpha+dAlpha/2.);
      b[1][ALGEBRAIC] = dT*(dBeta/2.+dAlpha/2.-dDelta/dAlpha*(1.+dAlpha));
      b[2][ALGEBRAIC] = dT*(dBeta/2.+dDelta);

   } else {
      a[1][ALGEBRAIC] = a[1][DIFFERENTIAL];
      b[0][ALGEBRAIC] = b[0][DIFFERENTIAL];
      b[1][ALGEBRAIC] = b[1][DIFFERENTIAL];
      b[2][ALGEBRAIC] = b[2][DIFFERENTIAL];
   }
   
   DEBUGCOUT("Algebraic coefficients:" << endl
	     << "beta  = " << dBeta << endl
	     << "delta = " << dDelta << endl
	     << "a2    = " << a[1][ALGEBRAIC] << endl
	     << "b0    = " << b[0][ALGEBRAIC] << endl
	     << "b1    = " << b[1][ALGEBRAIC] << endl
	     << "b2    = " << b[2][ALGEBRAIC] << endl);
      
   DEBUGCOUT("Asymptotic rho: " 
	     << -b[1][DIFFERENTIAL]/(2.*b[0][DIFFERENTIAL]) << endl
	     << "Discriminant: " 
	     << b[1][DIFFERENTIAL]*b[1][DIFFERENTIAL]-4.*b[2][DIFFERENTIAL]*b[0][DIFFERENTIAL] 
	     << endl
	     << "Asymptotic rho for algebraic variables: " 
	     << -b[1][ALGEBRAIC]/(2.*b[0][ALGEBRAIC]) << endl
	     << "Discriminant: " 
	     << b[1][ALGEBRAIC]*b[1][ALGEBRAIC]-4.*b[2][ALGEBRAIC]*b[0][ALGEBRAIC] 
	     << endl);
  
   /* Vengono modificati per la predizione, dopo che sono stati usati per
    * costruire gli altri coefficienti */
   mp[0] /= dT;
   mp[1] /= dT;
   
   /* valori di ritorno */
   db0Differential = b[0][DIFFERENTIAL];
   db0Algebraic = b[0][ALGEBRAIC];
}


doublereal NostroMetodo::dPredictDerivative(const doublereal& dXm1,
					    const doublereal& dXm2,
					    const doublereal& dXPm1,
					    const doublereal& dXPm2,
					    DofOrder::Order o) const
{
   if (o == DofOrder::ALGEBRAIC) {
      return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
   } /* else if (o == DofOrder::DIFFERENTIAL) */   
   return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
}


doublereal NostroMetodo::dPredictState(const doublereal& dXm1,
				       const doublereal& dXm2,
				       const doublereal& dXP,
				       const doublereal& dXPm1,
				       const doublereal& dXPm2,
				       DofOrder::Order o) const
{
   if (o == DofOrder::ALGEBRAIC) {
      return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2
	-a[1][ALGEBRAIC]*dXm1;
   } /* else if (o == DofOrder::DIFFERENTIAL) */   
   return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
     +b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1+b[2][DIFFERENTIAL]*dXPm2;
}

/* NostroMetodo - end */


/* Hope - begin */

Hope::Hope(const DriveCaller* pRho, const DriveCaller* pAlgRho)
: Rho(pRho), AlgebraicRho(pAlgRho), fStep(0)
{
   ASSERT(pRho != NULL);
   ASSERT(pAlgRho != NULL);
}


void Hope::SetCoef(doublereal dT,
		   doublereal dAlpha,
		   enum StepChange NewStep,
		   doublereal& db0Differential,
		   doublereal& db0Algebraic)
{
   /*
   if (dAlpha != 1.) {
      cerr << "HOPE time step integrator is not implemented yet in variable step form" << endl;
      THROW(ErrNotImplementedYet());
   }
    */
   
   if (NewStep == NEWSTEP) {
      ASSERT(fStep == flag(0) || fStep == flag(1));
      fStep = 1-fStep; // Commuta il valore di fStep
   }

   doublereal dTMod = dT*dAlpha;
   
   /* Differential coefficients */
   mp[0] = -6.*dAlpha*(1.+dAlpha);
   mp[1] = -mp[0];
   np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
   np[1] = dAlpha*(2.+3.*dAlpha);
      
   if (fStep) {
      b[0][DIFFERENTIAL] = b[1][DIFFERENTIAL] 
	= b[0][ALGEBRAIC] = b[1][ALGEBRAIC] 
	= db0Algebraic = db0Differential = dTMod/2.; // dT/4.;
           
   } else {
      doublereal dRho = Rho.dGet();
      doublereal dALPHA = 4.*dRho/(3.+dRho);      
      
      a[0][DIFFERENTIAL] = (4.-dALPHA)/3.;
      a[1][DIFFERENTIAL] = (dALPHA-1.)/3.;
      b[0][DIFFERENTIAL] = dTMod*(4.-dALPHA)/6.; // dT*(4.-dALPHA)/12.;
      b[1][DIFFERENTIAL] = dTMod*dALPHA/2.; // dT*dALPHA/4.;
      
      DEBUGCOUT("Predict()" << endl
		<< "Alpha = " << dAlpha << endl
		<< "Differential coefficients:" << endl
		<< "HOPE-Alpha = " << dALPHA << endl
		<< "a1    = " << a[0][DIFFERENTIAL] << endl
		<< "a2    = " << a[1][DIFFERENTIAL] << endl
		<< "b0    = " << b[0][DIFFERENTIAL] << endl
		<< "b1    = " << b[1][DIFFERENTIAL] << endl);
                  
      /* Coefficienti del metodo - variabili algebriche */
      doublereal dAlgebraicRho = AlgebraicRho.dGet();   
      doublereal dAlgebraicALPHA = 4.*dAlgebraicRho/(3.+dAlgebraicRho);     
            
      if (dAlgebraicRho != dRho) {                 
	 a[1][ALGEBRAIC] = (dAlgebraicALPHA-1.)/3.;
	 b[0][ALGEBRAIC] = dTMod*(4.-dAlgebraicALPHA)/6.; // dT*(4.-dAlgebraicALPHA)/12.;
	 b[1][ALGEBRAIC] = dTMod*dAlgebraicALPHA/2.; // dT*dAlgebraicALPHA/4.;
	 
      } else {
	 a[1][ALGEBRAIC] = a[1][DIFFERENTIAL];
	 b[0][ALGEBRAIC] = b[0][DIFFERENTIAL];
	 b[1][ALGEBRAIC] = b[1][DIFFERENTIAL];   
      }
      
      DEBUGCOUT("Algebraic coefficients:" << endl
		<< "HOPE-Alpha = " << dAlgebraicALPHA << endl
		<< "a2    = " << a[1][ALGEBRAIC] << endl
		<< "b0    = " << b[0][ALGEBRAIC] << endl
		<< "b1    = " << b[1][ALGEBRAIC] << endl);
            
      /* valori di ritorno */     
      db0Differential = b[0][DIFFERENTIAL];
      db0Algebraic = b[0][ALGEBRAIC];
   }
   
   /* Vengono modificati per la predizione, dopo che sono stati usati per
    * costruire gli altri coefficienti */
   mp[0] /= dT;
   mp[1] /= dT;
}


doublereal Hope::dPredictDerivative(const doublereal& dXm1,
				    const doublereal& dXm2,
				    const doublereal& dXPm1,
				    const doublereal& dXPm2,
				    DofOrder::Order o) const
{
   if (o == DofOrder::ALGEBRAIC) {
      return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
   } /* else if (o == DofOrder::DIFFERENTIAL) */   
   return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
}


doublereal Hope::dPredictState(const doublereal& dXm1,
			       const doublereal& dXm2,
			       const doublereal& dXP,
			       const doublereal& dXPm1,
			       const doublereal& /* dXPm2 */ ,
			       DofOrder::Order o) const
{
   if (fStep) {
      if (o == DofOrder::ALGEBRAIC) {
	 return b[0][ALGEBRAIC]*(dXP+dXPm1);
      } /* else if (o == DofOrder::DIFFERENTIAL) */
      return dXm1+b[0][ALGEBRAIC]*(dXP+dXPm1);
   } else {
      if (o == DofOrder::ALGEBRAIC) {
	 return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1
	   -a[1][ALGEBRAIC]*dXm1;
      } /* else if (o == DofOrder::DIFFERENTIAL) */
      return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
	+b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1;
   }
}

/* Hope - end */

/* MultiStepIntegrationMethods - end */
