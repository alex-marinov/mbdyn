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

/* Legami costitutivi */


#ifndef CONSTLTP__H
#define CONSTLTP__H

extern "C" {
#include <math.h>
/* #include <float.h> */
}

#include "constltp.h"

/* ElasticConstitutiveLaw - begin */

template <class T, class Tder>
class ElasticConstitutiveLaw
: public ConstitutiveLaw<T, Tder>, public TplDriveOwner<T> {
 protected:
   T PreStress;
   
   virtual std::ostream& Restart_(std::ostream& out) const {
      out << ", prestress, ", 
        Write(out, PreStress /* + GetF() */ , ", ");
      if (pGetDriveCaller()) {
	 out << ", prestrain, single, ",
	   Write(out, -Epsilon, ", ") << ", one /* ",
	   pGetDriveCaller()->Restart(out) << " */ ";
      }
      return out;
   };
   
 public:
   ElasticConstitutiveLaw(const TplDriveCaller<T>* pDC, const T& PStress) 
     : TplDriveOwner<T>(pDC), PreStress(PStress) {
      NO_OP;
   };
   
   virtual ~ElasticConstitutiveLaw(void) {
      NO_OP;
   };
};

typedef ElasticConstitutiveLaw<doublereal, doublereal> ElasticConstitutiveLaw1D;
typedef ElasticConstitutiveLaw<Vec3, Mat3x3> ElasticConstitutiveLaw3D;
typedef ElasticConstitutiveLaw<Vec6, Mat6x6> ElasticConstitutiveLaw6D;

/* ElasticConstitutiveLaw - end */


/* LinearElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffness;
   
 public:
   LinearElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
					 const T& PStress,
					 doublereal dStiff)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress), dStiffness(dStiff) { 
      ConstitutiveLaw<T, Tder>::FDE = dStiffness;   
   };
   
   virtual ~LinearElasticIsotropicConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;

      typedef LinearElasticIsotropicConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               dStiffness));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear elastic isotropic, " << dStiffness;
      return Restart_(out);
   };
   
   virtual void Update(const T& Eps, const T& /* EpsPrime */  = 0.) {
      ElasticConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+(ElasticConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get())*dStiffness;
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& /* EpsPrime */ = 0.) {
      Update(ElasticConstitutiveLaw<T, Tder>::Epsilon+DeltaEps);
   };
};

typedef LinearElasticIsotropicConstitutiveLaw<doublereal, doublereal> LinearElasticIsotropicConstitutiveLaw1D;
typedef LinearElasticIsotropicConstitutiveLaw<Vec3, Mat3x3> LinearElasticIsotropicConstitutiveLaw3D;
typedef LinearElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> LinearElasticIsotropicConstitutiveLaw6D;


/* LinearElasticIsotropicConstitutiveLaw - end */


/* LinearElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticGenericConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {   
 public:
   LinearElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
				const T& PStress, 
				const Tder& Stiff) 
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      ConstitutiveLaw<T, Tder>::FDE = Stiff;
   };
   
   virtual ~LinearElasticGenericConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef LinearElasticGenericConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               ConstitutiveLaw<T, Tder>::FDE));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear elastic generic, ",
        Write(out, ConstitutiveLaw<T, Tder>::FDE, ", ");
      return Restart_(out);
   };
   
   virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get());
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& /* EpsPrime */ = 0.) {
      Update(ConstitutiveLaw<T, Tder>::Epsilon+DeltaEps);
   };
};

typedef LinearElasticGenericConstitutiveLaw<doublereal, doublereal> LinearElasticGenericConstitutiveLaw1D;
typedef LinearElasticGenericConstitutiveLaw<Vec3, Mat3x3> LinearElasticGenericConstitutiveLaw3D;
typedef LinearElasticGenericConstitutiveLaw<Vec6, Mat6x6> LinearElasticGenericConstitutiveLaw6D;

/* LinearElasticGenericConstitutiveLaw - end */


/* LinearElasticGenericAxialTorsionCouplingConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticGenericAxialTorsionCouplingConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<T>* pDC,
							   const T& PStress,
							   const Tder& = 0.,
							   doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {     
      THROW(Err(std::cerr, "axial-torsion coupling constitutive law "
			      "is allowed only for beams (6x6)"));
   };
   
   virtual ~LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
      NO_OP;
   };

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
   
   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};
   
   
class LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6>
: public ElasticConstitutiveLaw6D {   
 protected:
   Mat6x6 Stiffness;
   doublereal dRefTorsion;
   doublereal dAxialTorsionCoupling;
   
 public:
   LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<Vec6>* pDC,
							   const Vec6& PStress, 
							   const Mat6x6& Stiff,
							   doublereal dAxTors)
     : ElasticConstitutiveLaw6D(pDC, PStress), 
     Stiffness(Stiff), 
     dRefTorsion(Stiffness.dGet(4, 4)),
     dAxialTorsionCoupling(dAxTors) {	
      FDE = Stiffness;    
   };
   
   virtual ~LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy(void) const {
      ConstitutiveLaw<Vec6, Mat6x6>* pCL = NULL;
      
      typedef LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress, 
                               Stiffness,
                               dAxialTorsionCoupling));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear eleastic generic axial torsion coupling, ",
        Write(out, FDE, ", ") << ", " << dAxialTorsionCoupling;
      return Restart_(out);
   };
   
   virtual void Update(const Vec6& Eps, const Vec6& /* EpsPrime */ = 0.) {
      Epsilon = Eps;
      doublereal d = Epsilon.dGet(1);    
      FDE.Put(4, 4, dRefTorsion+d*dAxialTorsionCoupling);
      F = PreStress+FDE*(Epsilon-Get());
   };
   
   virtual void IncrementalUpdate(const Vec6& DeltaEps, const Vec6& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};

/* LinearElasticGenericAxialTorsionCouplingConstitutiveLaw - end */


/* LogConstitutiveLaw - begin */

template <class T, class Tder>
class LogConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
   
 public:
   LogConstitutiveLaw(const TplDriveCaller<T>* pDC,
		      const T& PStress, 
		      doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {      
      THROW(Err(std::cerr, "log constitutive law is allowed only for rods"));
   };
   
   virtual ~LogConstitutiveLaw(void) { 
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
   
   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};


class LogConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;
   doublereal dCurrEps;
   
 public:
   LogConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
		      const doublereal& PStress, 
		      doublereal dStiff)
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff), dCurrEps(0.) {
      ASSERT(Get() < 1.);	      
   };
   
   virtual ~LogConstitutiveLaw(void) { 
      NO_OP;
   };
   
   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;
      
      typedef LogConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress, 
                               dStiffness));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "log elastic, " << dStiffness;
      return Restart_(out);
   };
   
   virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.) {
      Epsilon = Eps;
      
      doublereal dPreStrain = Get();
      dCurrEps = 1.+Epsilon-dPreStrain;
      ASSERT(dCurrEps > DBL_EPSILON);
      
      F = PreStress+dStiffness*log(dCurrEps);
      FDE = dStiffness/dCurrEps;        
   };
   
   virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};

/* LogConstitutiveLaw - end */


/* DoubleLinearElasticConstitutiveLaw - begin */

template <class T, class Tder>
class DoubleLinearElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {   
 public:
   DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
				      const T& PStress,
				      doublereal = 0.,
				      doublereal = 0.,
				      doublereal = 0.,
				      doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {    
      THROW(Err(std::cerr, "double linear elastic constitutive law "
			      "is allowed only for rods and 3D hinges"));
   };
   
   virtual ~DoubleLinearElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
   
   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};


class DoubleLinearElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   flag fSecondStiff;    
   
 public:
   DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
				      const doublereal& PStress,
				      doublereal dStiff,
				      doublereal dUppLimStrain,
				      doublereal dLowLimStrain,
				      doublereal dSecondStiff)
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff), 
     dUpperLimitStrain(dUppLimStrain),
     dLowerLimitStrain(dLowLimStrain),
     dSecondStiffness(dSecondStiff) {
      FDE = dStiffness;
   };
       
   virtual ~DoubleLinearElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;
      
      typedef DoubleLinearElasticConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress, 
                               dStiffness, 
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear elastic, "
        << dStiffness << ", "
        << dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness;
      return Restart_(out);
   };
   
   virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.) {
      Epsilon = Eps;   
      
      doublereal dPreStrain = Get();
      doublereal dCurrStrain = Epsilon-dPreStrain;
      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE = dStiffness;
	 F = PreStress+dStiffness*dCurrStrain;
      } else {      
	 FDE = dSecondStiffness;
	 
	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress+dStiffness*dUpperLimitStrain
	      +dSecondStiffness*(dCurrStrain-dUpperLimitStrain);
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {	    
	    F = PreStress+dStiffness*dLowerLimitStrain
	      +dSecondStiffness*(dCurrStrain-dLowerLimitStrain);
	 }
      }
   };
   
   virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};


class DoubleLinearElasticConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
 private:
   doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   
 public:
   DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
				      const Vec3& PStress,
				      doublereal dStiff,
				      doublereal dUppLimStrain,
				      doublereal dLowLimStrain,
				      doublereal dSecondStiff)
     : ElasticConstitutiveLaw3D(pDC, PStress),
     dStiffness(dStiff), 
     dUpperLimitStrain(dUppLimStrain),
     dLowerLimitStrain(dLowLimStrain),
     dSecondStiffness(dSecondStiff) {
      FDE = dStiffness;
   };
       
   virtual ~DoubleLinearElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
      ConstitutiveLaw<Vec3, Mat3x3>* pCL = NULL;
      
      typedef DoubleLinearElasticConstitutiveLaw<Vec3, Mat3x3> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress, 
                               dStiffness, 
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear elastic, "
        << dStiffness << ", "
	<< dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness;
      return Restart_(out);
   };
   
   virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */ = 0.) {
      Epsilon = Eps;   
      
      Vec3 PreStrain = Get();
      Vec3 CurrStrain = Epsilon-PreStrain;
      doublereal dCurrStrain = CurrStrain.dGet(3);
      
      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE.Put(3, 3, dStiffness);
	 F = PreStress+CurrStrain*dStiffness;
      } else {      
	 FDE.Put(3, 3, dSecondStiffness);
	 
	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress+Vec3(CurrStrain.dGet(1)*dStiffness,
			       CurrStrain.dGet(2)*dStiffness,
			       dUpperLimitStrain*dStiffness
			       +(dCurrStrain-dUpperLimitStrain)*dSecondStiffness);
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {	    
	    F = PreStress+Vec3(CurrStrain.dGet(1)*dStiffness,
			       CurrStrain.dGet(2)*dStiffness,
			       dLowerLimitStrain*dStiffness
			       +(dCurrStrain-dLowerLimitStrain)*dSecondStiffness);
	 }
      }
   };
   
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};

/* DoubleLinearElasticConstitutiveLaw - end */


/* IsotropicHardeningConstitutiveLaw - begin */

template <class T, class Tder>
class IsotropicHardeningConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffness;
   doublereal dAlpha;
   /* Legge costitutiva:
    * 
    *          k*Alpha*x^3
    *    f = ---------------
    *         (1+Alpha*x^2)
    */
   
 public:
   IsotropicHardeningConstitutiveLaw(const TplDriveCaller<T>* pDC,
				     const T& PStress, 
				     doublereal dStiff, 
				     doublereal dEpsHard) 
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
     dStiffness(dStiff), dAlpha(0.) {
      ASSERT(dEpsHard > DBL_EPSILON);
      dAlpha = 3./(dEpsHard*dEpsHard);
   };
   
   virtual ~IsotropicHardeningConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef IsotropicHardeningConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               dStiffness,
                               dAlpha));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "isoropic hardening elastic, " << dStiffness << ", " 
        << sqrt(3./dAlpha);
      return Restart_(out);
   };      
   
   virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;         
      T x = ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get();
      doublereal dx2 = x*x;
      doublereal dDen = 1.+dAlpha*dx2;
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+x*(dStiffness*dAlpha*dx2/dDen);
      ConstitutiveLaw<T, Tder>::FDE = dStiffness*dAlpha*(3.+dAlpha*dx2)*dx2/(dDen*dDen);      
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& /* EpsPrime */ = 0.) {
      Update(ConstitutiveLaw<T, Tder>::Epsilon+DeltaEps);
   };
};

typedef IsotropicHardeningConstitutiveLaw<doublereal, doublereal> IsotropicHardeningConstitutiveLaw1D;
typedef IsotropicHardeningConstitutiveLaw<Vec3, Mat3x3> IsotropicHardeningConstitutiveLaw3D;
typedef IsotropicHardeningConstitutiveLaw<Vec6, Mat6x6> IsotropicHardeningConstitutiveLaw6D;

/* IsotropicHardeningConstitutiveLaw - end */


/* ContactConstitutiveLaw - begin */

template <class T, class Tder>
class ContactConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
   
 public:
   ContactConstitutiveLaw(const TplDriveCaller<T>* pDC,
			  const T& PStress,
			  const doublereal = 0.,
			  const doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {      
      THROW(Err(std::cerr, "contact constitutive law "
			      "is allowed only for rods"));
   };
   
   virtual ~ContactConstitutiveLaw(void) { 
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
   
   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};



/*
 * 
 * F = k * ( 1-l_0 / l ) ^ gamma
 *     ==>
 * F = k * ( eps / ( 1 + eps ) ) ^ gamma
 * 
 * d F = gamma * k * ( eps / ( 1 + eps ) ) ^ ( gamma - 1 ) * 1 / ( 1 + eps ) ^ 2 d eps
 * 
 */

class ContactConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dKappa;
   doublereal dGamma;
   
 public:
   ContactConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			  const doublereal& PStress,
			  const doublereal& dKappa,			  
			  const doublereal& dGamma)
     : ElasticConstitutiveLaw1D(pDC, PStress), dKappa(dKappa), dGamma(dGamma) {
	NO_OP;
   };
   
   virtual ~ContactConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;
      
      typedef ContactConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(), 
                               PreStress, 
                               dKappa,
                               dGamma));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "contact elastic, "
        << dKappa << ", "
	<< dGamma;
      return Restart_(out);
   };
   
   virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */  = 0.) {
      doublereal dE;
      
      Epsilon = Eps;
      dE = Epsilon-Get();
      if ( dE >= 0. ) {
	 F = PreStress;
	 FDE = 0.;
      } else {
	 F = PreStress+dKappa*(1.-1./pow(1.+dE, dGamma));
	 FDE = dGamma*dKappa/pow(1.+dE, dGamma+1.);
      }
   };
   
   virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};


class ContactConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
 private:
   doublereal dKappa;
   doublereal dGamma;
   
 public:
   ContactConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
			  const Vec3& PStress,
			  const doublereal& dKappa,
			  const doublereal& dGamma)
     : ElasticConstitutiveLaw3D(pDC, PStress), dKappa(dKappa), dGamma(dGamma) {
	F = Zero3;
	FDE = Zero3x3;
   };
   
   virtual ~ContactConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
      ConstitutiveLaw<Vec3, Mat3x3>* pCL = NULL;
      
      typedef ContactConstitutiveLaw<Vec3, Mat3x3> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(), 
                               PreStress, 
                               dKappa,
                               dGamma));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "contact elastic, "
        << dKappa << ", "
	<< dGamma;
      return Restart_(out);
   };
   
   virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */  = 0.) {
      doublereal dE;
      
      Epsilon = Eps;
      dE = Epsilon.dGet(3)-Get().dGet(3);
      if ( dE >= 0. ) {
	 F.Put(3, 0.);
	 FDE.Put(3, 3, 0.);
      } else {
	 F.Put(3, dKappa*(1.-1./pow(1.+dE, dGamma)));
	 FDE.Put(3, 3, dGamma*dKappa/pow(1.+dE, dGamma+1.));
      }
   };
   
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& /* EpsPrime */ = 0.) {
      Update(Epsilon+DeltaEps);
   };
};

typedef ContactConstitutiveLaw<doublereal, doublereal> ContactConstitutiveLaw1D;
typedef ContactConstitutiveLaw<Vec3, Mat3x3> ContactConstitutiveLaw3D;
typedef ContactConstitutiveLaw<Vec6, Mat6x6> ContactConstitutiveLaw6D;

/* ContactConstitutiveLaw - end */


/* LinearViscousIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscousIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffnessPrime;  /* Isotropa: Eye*dStiffnessPrime */
   
 public:
   LinearViscousIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
					 const T& PStress,
					 doublereal dStiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
     dStiffnessPrime(dStiffPrime) {
      ConstitutiveLaw<T, Tder>::FDEPrime = dStiffnessPrime;
   };
   
   virtual ~LinearViscousIsotropicConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef LinearViscousIsotropicConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               dStiffnessPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscous isotropic, "
        << dStiffnessPrime;
      return Restart_(out);
   };
   
   virtual void Update(const T& /* Eps */ , const T& EpsPrime = 0.) {
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::EpsilonPrime*dStiffnessPrime;
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
      Update(DeltaEps, EpsPrime);
   };
};

/* LinearViscousIsotropicConstitutiveLaw - end */


/* LinearViscousGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscousGenericConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {   
 public:
   LinearViscousGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
				       const T& PStress,
				       const Tder& StiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      ConstitutiveLaw<T, Tder>::FDEPrime = StiffPrime;
   };
   
   virtual ~LinearViscousGenericConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef LinearViscousGenericConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               ConstitutiveLaw<T, Tder>::FDEPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscous generic, ",
        Write(out, ConstitutiveLaw<T, Tder>::FDEPrime, ", ");
      return Restart_(out);
   };
   
   virtual void Update(const T& /* Eps */ , const T& EpsPrime = 0.) {
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
      Update(DeltaEps, EpsPrime);
   };
};

/* LinearViscousGenericConstitutiveLaw - end */


/* LinearViscoElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffness;
   doublereal dStiffnessPrime;
   
 public:
   LinearViscoElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
					      const T& PStress,
					      doublereal dStiff,
					      doublereal dStiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress), 
     dStiffness(dStiff), dStiffnessPrime(dStiffPrime) { 
      ConstitutiveLaw<T, Tder>::FDE = dStiffness;
      ConstitutiveLaw<T, Tder>::FDEPrime = dStiffnessPrime;          
   };
   
   virtual ~LinearViscoElasticIsotropicConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef LinearViscoElasticIsotropicConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               dStiffness,
                               dStiffnessPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscoelastic isotropic, "
	<< dStiffness << ", "
	<< dStiffnessPrime;
      return Restart_(out);      
   };
   
   virtual void Update(const T& Eps, const T& EpsPrime = 0.) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
          
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get())*dStiffness+ConstitutiveLaw<T, Tder>::EpsilonPrime*dStiffnessPrime;
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
      Update(ConstitutiveLaw<T, Tder>::Epsilon+DeltaEps, EpsPrime);
   };
};

typedef LinearViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal> LinearViscoElasticIsotropicConstitutiveLaw1D;
typedef LinearViscoElasticIsotropicConstitutiveLaw<Vec3, Mat3x3> LinearViscoElasticIsotropicConstitutiveLaw3D;
typedef LinearViscoElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticIsotropicConstitutiveLaw6D;

/* LinearViscoElasticIsotropicConstitutiveLaw - end */


/* LinearViscoElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticGenericConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   LinearViscoElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
					    const T& PStress,
					    const Tder& Stiff,
					    const Tder& StiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      ConstitutiveLaw<T, Tder>::FDE = Stiff;
      ConstitutiveLaw<T, Tder>::FDEPrime = StiffPrime;
   };
   
   virtual ~LinearViscoElasticGenericConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = NULL;
      
      typedef LinearViscoElasticGenericConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
                               ElasticConstitutiveLaw<T, Tder>::PreStress, 
                               ConstitutiveLaw<T, Tder>::FDE,
                               ConstitutiveLaw<T, Tder>::FDEPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
     out << "linear viscoelastic generic, ",
       Write(out, ConstitutiveLaw<T, Tder>::FDE, ", ") << ", ",
       Write(out, ConstitutiveLaw<T, Tder>::FDEPrime, ", ");
       return Restart_(out);
   };
   
   virtual void Update(const T& Eps, const T& EpsPrime = 0.) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get())
	+ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
   };
   
   virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
      Update(ConstitutiveLaw<T, Tder>::Epsilon+DeltaEps, EpsPrime);
   };
};

typedef LinearViscoElasticGenericConstitutiveLaw<doublereal, doublereal> LinearViscoElasticGenericConstitutiveLaw1D;
typedef LinearViscoElasticGenericConstitutiveLaw<Vec3, Mat3x3> LinearViscoElasticGenericConstitutiveLaw3D;
typedef LinearViscoElasticGenericConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticGenericConstitutiveLaw6D;

/* LinearViscoElasticGenericConstitutiveLaw - end */


/* DoubleLinearViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class DoubleLinearViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
					   const T& PStress,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.) 
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {             
      THROW(Err(std::cerr, "doublelinear viscoelastic constitutive law "
			      "is allowed only for rods ad 3D hinges"));
   };
   
   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };

   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};


class DoubleLinearViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   doublereal dStiffnessPrime;
   
 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
					   const doublereal& PStress,
					   doublereal dStiff,
					   doublereal dUpp,
					   doublereal dLow,
					   doublereal dSecondS,
					   doublereal dStiffPrime) 
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff),
     dUpperLimitStrain(dUpp), dLowerLimitStrain(dLow),
     dSecondStiffness(dSecondS),
     dStiffnessPrime(dStiffPrime) { 
      FDEPrime = dStiffnessPrime;
   };
   
   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;
      
      typedef DoubleLinearViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(), 
                               PreStress, 
                               dStiffness,
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness,
                               dStiffnessPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear viscoelastic, "
	<< dStiffness << ", "
	<< dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness << ", "
	<< dStiffnessPrime << ", ";
      return Restart_(out);
   };
   
   virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
      Epsilon = Eps;
      EpsilonPrime = EpsPrime;
      
      doublereal dPreStrain = Get();
      doublereal dCurrStrain = Epsilon-dPreStrain;
      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE = dStiffness;
	 F = PreStress+dStiffness*dCurrStrain
	   +dStiffnessPrime*EpsilonPrime;
      } else {      
	 FDE = dSecondStiffness;
	 
	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress+dStiffness*dUpperLimitStrain
	      +dSecondStiffness*(dCurrStrain-dUpperLimitStrain)
		+dStiffnessPrime*EpsilonPrime;
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {	    
	    F = PreStress+dStiffness*dLowerLimitStrain
	      +dSecondStiffness*(dCurrStrain-dLowerLimitStrain)
		+dStiffnessPrime*EpsilonPrime;
	 }
      }
   };

   virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& EpsPrime = 0.) {
      Update(Epsilon+DeltaEps, EpsPrime);
   };
};


class DoubleLinearViscoElasticConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
 private:
   doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   doublereal dStiffnessPrime;
   
 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
					   const Vec3& PStress,
					   doublereal dStiff,
					   doublereal dUppLimStrain,
					   doublereal dLowLimStrain,
					   doublereal dSecondStiff,
					   doublereal dStiffPrime)
     : ElasticConstitutiveLaw3D(pDC, PStress),
     dStiffness(dStiff), 
     dUpperLimitStrain(dUppLimStrain),
     dLowerLimitStrain(dLowLimStrain),
     dSecondStiffness(dSecondStiff),
     dStiffnessPrime(dStiffPrime) {
      FDE = dStiffness;
      FDEPrime = dStiffnessPrime;
   };
       
   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
      ConstitutiveLaw<Vec3, Mat3x3>* pCL = NULL;
      
      typedef DoubleLinearViscoElasticConstitutiveLaw<Vec3, Mat3x3> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(), 
                               PreStress, 
                               dStiffness,
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness,
                               dStiffnessPrime));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear viscoelastic, "
        << dStiffness << ", "
	<< dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness << ", "
	<< dStiffnessPrime << ", ";
      return Restart_(out);
   };
   
   virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = 0.) {
      Epsilon = Eps;   
      EpsilonPrime = EpsPrime;
      
      Vec3 PreStrain = Get();
      Vec3 CurrStrain = Epsilon-PreStrain;
      doublereal dCurrStrain = CurrStrain.dGet(3);
      
      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE.Put(3, 3, dStiffness);
	 F = PreStress+CurrStrain*dStiffness+EpsilonPrime*dStiffnessPrime;
      } else {      
	 FDE.Put(3, 3, dSecondStiffness);
	 
	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress
	      +Vec3(CurrStrain.dGet(1)*dStiffness,
		    CurrStrain.dGet(2)*dStiffness,
		    dUpperLimitStrain*dStiffness
		    +(dCurrStrain-dUpperLimitStrain)*dSecondStiffness)
		+EpsilonPrime*dStiffnessPrime;
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {	    
	    F = PreStress
	      +Vec3(CurrStrain.dGet(1)*dStiffness,
		    CurrStrain.dGet(2)*dStiffness,
		    dLowerLimitStrain*dStiffness
		    +(dCurrStrain-dLowerLimitStrain)*dSecondStiffness)
		+EpsilonPrime*dStiffnessPrime;
	 }
      }
   };
   
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& EpsPrime = 0.) {
      Update(Epsilon+DeltaEps, EpsPrime);
   };
};

/* DoubleLinearViscoElasticConstitutiveLaw - end */


/* TurbulentViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class TurbulentViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   TurbulentViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
					const T& PStress,
					doublereal = 0.,
					doublereal = 0.,
					doublereal = 0.,
					doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {     
      THROW(Err(std::cerr, "Turbulent viscoelastic constitutive law "
			      "is allowed only for rods"));
   };
   
   virtual ~TurbulentViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return NULL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };
   
   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
   
   virtual void IncrementalUpdate(const T& /* DeltaEps */ , const T& /* EpsPrime */ = 0.) {
      NO_OP;
   };
};


class TurbulentViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;
   doublereal dStiffnessPrime;
   doublereal dTreshold;
   doublereal dParabolicStiffness;
   
 public:
   TurbulentViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
					const doublereal& PStress,
					doublereal dStiff,
					doublereal dStiffPrime,
					doublereal dTres,
					doublereal dParabStiff)
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff), dStiffnessPrime(dStiffPrime),
     dTreshold(dTres), dParabolicStiffness(dParabStiff) { 
      FDE = dStiffness;
   };
   
   virtual ~TurbulentViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;
      
      typedef TurbulentViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL, 
                            cl, 
                            cl(pGetDriveCaller()->pCopy(), 
                               PreStress, 
                               dStiffness,
                               dStiffnessPrime,
                               dTreshold, 
                               dParabolicStiffness));
      
      return pCL;
   };
   
   virtual std::ostream& Restart(std::ostream& out) const {
      out << "turbulent viscoelastic, " 
	<< dStiffness << ", "
	<< dStiffnessPrime << ", "
	<< dTreshold << ", "
	<< dParabolicStiffness << ", ";
      return Restart_(out);
   };
   
   virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
      Epsilon = Eps;
      EpsilonPrime = EpsPrime;
      
      doublereal dPreStrain = Get();
      
      doublereal d = fabs(EpsilonPrime);
      if (d < dTreshold) {
	 FDEPrime = dStiffnessPrime;
	 F = PreStress+dStiffness*(Epsilon-dPreStrain)
	   +dStiffnessPrime*EpsilonPrime;
      } else {
	 FDEPrime = 2.*dParabolicStiffness*d;
	 F = PreStress+dStiffness*(Epsilon-dPreStrain)
	   +dParabolicStiffness*d*EpsilonPrime;
      }	
   };
   
   virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& EpsPrime = 0.) {
      Update(Epsilon+DeltaEps, EpsPrime);
   };
};

/* TurbulentViscoElasticConstitutiveLaw - end */


/* LinearViscoElasticBiStopConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticBiStopConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	enum Status { INACTIVE, ACTIVE };
private:
	enum Status status;
	const DriveCaller *pActivatingCondition;
	const DriveCaller *pDeactivatingCondition;
	Tder FDECurr;
	Tder FDEPrimeCurr;
	T EpsRef;
public:
	LinearViscoElasticBiStopConstitutiveLaw(
			const TplDriveCaller<T>* pDC,
			const T& PStress,
			const Tder& Stiff,
			const Tder& StiffPrime,
			enum Status initialStatus,
			const DriveCaller *pA,
			const DriveCaller *pD
	) : ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
	status(initialStatus), 
	pActivatingCondition(pA), pDeactivatingCondition(pD), EpsRef(0.) {
		ASSERT(pActivatingCondition != NULL);
		ASSERT(pDeactivatingCondition != NULL);
		FDECurr = Stiff;
		FDEPrimeCurr = StiffPrime;
		if (status == ACTIVE) {
			ConstitutiveLaw<T, Tder>::FDE = Stiff;
			ConstitutiveLaw<T, Tder>::FDEPrime = StiffPrime;
		}
	};
   
	virtual 
	~LinearViscoElasticBiStopConstitutiveLaw(void) {
		NO_OP;
	};
	
	virtual 
	ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;
		
		typedef LinearViscoElasticBiStopConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, 
			cl, 
			cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<T, Tder>::PreStress, 
				ConstitutiveLaw<T, Tder>::FDE,
				ConstitutiveLaw<T, Tder>::FDEPrime,
				status, 
				pActivatingCondition->pCopy(), 
				pDeactivatingCondition->pCopy()));
		
		return pCL;
	};
   
	virtual std::ostream& 
	Restart(std::ostream& out) const {
		out << "linear viscoelastic bistop, ",
			Write(out, ConstitutiveLaw<T, Tder>::FDE, ", ") << ", ",
			Write(out, ConstitutiveLaw<T, Tder>::FDEPrime, ", ")
			<< ", initial status, ";
		if (status == INACTIVE) {
			out << "inactive";
		} else {
			out << "active";
		}
		out << ", ",
			pActivatingCondition->Restart(out) << ", ",
			pDeactivatingCondition->Restart(out);
		return Restart_(out);
	};
	
	virtual void 
	Update(const T& Eps, const T& EpsPrime = 0.) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		switch (status) {
		case INACTIVE:
			if (pActivatingCondition->dGet() == 0.) {
				break;
			}
			status = ACTIVE;
			EpsRef = ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get();
			ConstitutiveLaw<T, Tder>::FDE = FDECurr;
			ConstitutiveLaw<T, Tder>::FDEPrime = FDEPrimeCurr;

		case ACTIVE:
			if (pDeactivatingCondition->dGet() != 0.) {
				status = INACTIVE;
				ConstitutiveLaw<T, Tder>::FDE = 0.;
				ConstitutiveLaw<T, Tder>::FDEPrime = 0.;
				ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress;
				break;
			}
			ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
				+ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get()-EpsRef)
				+ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
			break;
		}
	};
	
	virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
		Update(ConstitutiveLaw<T, Tder>::Epsilon+DeltaEps, EpsPrime);
	};
};

typedef LinearViscoElasticBiStopConstitutiveLaw<doublereal, doublereal> LinearViscoElasticBiStopConstitutiveLaw1D;
typedef LinearViscoElasticBiStopConstitutiveLaw<Vec3, Mat3x3> LinearViscoElasticBiStopConstitutiveLaw3D;
typedef LinearViscoElasticBiStopConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticBiStopConstitutiveLaw6D;

/* LinearViscoElasticBiStopConstitutiveLaw - end */
















/* TurbulentViscousIsotropicConstitutiveLaw3D - begin */
#if 0
class TurbulentViscousIsotropicConstitutiveLaw3D : public ConstitutiveLaw3D {
 private:
   doublereal dParabolicStiffnessPrime;
   
 public:
   TurbulentViscousIsotropicConstitutiveLaw3D(doublereal dParabStiffPrime);
   virtual ~TurbulentViscousIsotropicConstitutiveLaw3D(void);
   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = 0.);
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& EpsPrime = 0.);
};
#endif
/* TurbulentViscousIsotropicConstitutiveLaw3D - end */


/* LinearTurbulentViscousIsotropicConstitutiveLaw3D - begin */
#if 0
class LinearTurbulentViscousIsotropicConstitutiveLaw3D : public ConstitutiveLaw3D {
 private:
   doublereal dStiffnessPrime;  /* Isotropa: Eye*dStiffnessPrime */
   doublereal dTreshold;
   doublereal dParabolicStiffnessPrime;
   
 public:
   LinearTurbulentViscousIsotropicConstitutiveLaw3D(doublereal dStiffPrime,
						    doublereal dTres,
						    doublereal dParabStiffPrime);
   virtual ~LinearTurbulentViscousIsotropicConstitutiveLaw3D(void);
   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = 0.);
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& EpsPrime = 0.);
};
#endif
/* LinearTurbulentViscousIsotropicConstitutiveLaw3D - end */


/* VelocityDamperConstitutiveLaw3D - begin */
#if 0
class VelocityDamperConstitutiveLaw3D : public ConstitutiveLaw3D {
 private:
   doublereal dRefValue;        /* Forza di riferimento */
   doublereal dStiffnessPrime;  /* Isotropa: Eye*dStiffnessPrime */
   doublereal dRefVelocity;     /* Velocita' di riferimento */
   
 public:
   VelocityDamperConstitutiveLaw3D(doublereal dRefVal, 
				   doublereal dStiffPrime, 
				   doublereal dRefVel);
   virtual ~VelocityDamperConstitutiveLaw3D(void);
   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = 0.);
   virtual void IncrementalUpdate(const Vec3& DeltaEps, const Vec3& EpsPrime = 0.);
};
#endif
/* VelocityDamperConstitutiveLaw3D - end */

#endif /* CONSTLTP__H */
