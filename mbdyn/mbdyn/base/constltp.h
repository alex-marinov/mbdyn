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

/* Legami costitutivi */

#ifndef CONSTLTP_H
#define CONSTLTP_H

#include <tpldrive.h>
#include <matvec3.h>
#include <matvec6.h>

/* Tipi di cerniere deformabili */
class DefHingeType {
 public:
   enum Type {
      UNKNOWN = -1,
	
	ELASTIC = 0,
	VISCOUS,
	VISCOELASTIC,
	
	LASTDEFHINGETYPE
   };
};

/* ConstitutiveLaw - begin */
template <class T, class Tder>
class ConstitutiveLaw {
 public:
   class ErrNotAvailable {
    public:
      ErrNotAvailable(void) {
	 cerr << "Constitutive law not available" << endl;
      };
      ErrNotAvailable(ostream& out) {
	 out << "Constitutive law not available" << endl;
      };
      ErrNotAvailable(ostream& out, const char* const s) {
	 out << s << endl;
      };
   };
   typedef ConstitutiveLaw<T, Tder>::ErrNotAvailable Err;   
   
 protected:
   T Epsilon;
   T EpsilonPrime;
   
   T F;
   Tder FDE;
   Tder FDEPrime;

   /*
    * To allow the duplication of the constitutive laws, 
    * a "virtual constructor-like" function must be implemented,
    * with type resolution for the constitutive laws and
    * for the tpldrives too.
    */
   
 public:
   ConstitutiveLaw(void)
     : Epsilon(0.), EpsilonPrime(0.),
     F(0.), FDE(0.), FDEPrime(0.) {
      NO_OP;
   };
   
   virtual ~ConstitutiveLaw(void) {
      NO_OP;
   };
   
   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const = 0;
   
   virtual ostream& Restart(ostream& out) const = 0;
   
   virtual void Update(const T& Eps, const T& EpsPrime = 0.) = 0;
   virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) = 0;
   
   virtual const T& GetF(void) const {
      return F;
   };
   
   virtual const Tder& GetFDE(void) const {
      return FDE;
   };
   
   virtual const Tder& GetFDEPrime(void) const {
      return FDEPrime;
   };
};

typedef ConstitutiveLaw<doublereal, doublereal> ConstitutiveLaw1D;
typedef ConstitutiveLaw<Vec3, Mat3x3> ConstitutiveLaw3D;
typedef ConstitutiveLaw<Vec6, Mat6x6> ConstitutiveLaw6D;

/* ConstitutiveLaw - end */


/* ConstitutiveLawOwner - begin */

template <class T, class Tder>
class ConstitutiveLawOwner {
 protected:
   const ConstitutiveLaw<T, Tder>* pConstLaw;
   
 public:
   ConstitutiveLawOwner(const ConstitutiveLaw<T, Tder>* pCL)
     : pConstLaw(pCL) { 
      ASSERT(pCL != NULL);
   };
   
   virtual ~ConstitutiveLawOwner(void) {
      ASSERT(pConstLaw != NULL);
      if (pConstLaw != NULL) {
	 SAFEDELETE(pConstLaw, DMmm);
      }	
   };
   
   inline ConstitutiveLaw<T, Tder>* pGetConstLaw(void) const {
      ASSERT(pConstLaw != NULL);
      return (ConstitutiveLaw<T, Tder>*)pConstLaw; 
   };
   
   inline void Update(const T& Eps, const T& EpsPrime = 0.) {
      ASSERT(pConstLaw != NULL);
      ((ConstitutiveLaw<T, Tder>*)pConstLaw)->Update(Eps, EpsPrime);
   };
   
   inline void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
      ASSERT(pConstLaw != NULL);
      ((ConstitutiveLaw<T, Tder>*)pConstLaw)->IncrementalUpdate(DeltaEps, EpsPrime);
   };
   
   inline const T& GetF(void) const {
      ASSERT(pConstLaw != NULL);
      return pConstLaw->GetF();
   };
   
   inline const Tder& GetFDE(void) const {	
      ASSERT(pConstLaw != NULL);
      return pConstLaw->GetFDE();
   };
   
   inline const Tder& GetFDEPrime(void) const {
      ASSERT(pConstLaw != NULL);
      return pConstLaw->GetFDEPrime();
   };
};

typedef ConstitutiveLawOwner<doublereal, doublereal> ConstitutiveLaw1DOwner;
typedef ConstitutiveLawOwner<Vec3, Mat3x3> ConstitutiveLaw3DOwner;
typedef ConstitutiveLawOwner<Vec6, Mat6x6> ConstitutiveLaw6DOwner;

/* ConstitutiveLawOwner - end */

#endif /* CONSTLTP_H */
