/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* Symbolic constitutive law */


#ifndef SYMCLTP_H
#define SYMCLTP_H

#include <constltp.h>

/* SymbolicElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicElasticIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const char *epsilon,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicElasticIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
			const T& /* EpsPrime */ = 0.);
};

typedef SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal> 
	SymbolicElasticIsotropicConstitutiveLaw1D;
typedef SymbolicElasticIsotropicConstitutiveLaw<Vec3, Mat3x3>
	SymbolicElasticIsotropicConstitutiveLaw3D;
typedef SymbolicElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicElasticIsotropicConstitutiveLaw6D;

template <class T, class Tder>
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::SymbolicElasticIsotropicConstitutiveLaw(
		const TplDriveCaller<T>* pDC, const T& PStress, 
		const char *epsilon, const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{ 
	typedef typename ConstitutiveLaw<T, Tder>::ErrNotAvailable E;
	THROW(E(std::cerr, "symbolic constitutive law "
				"is allowed only for scalar")); 
}
 
template <class T, class Tder>
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::~SymbolicElasticIsotropicConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return NULL;
}

template <class T, class Tder> std::ostream& 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

/* SymbolicElasticIsotropicConstitutiveLaw - end */

/* SymbolicViscousIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicViscousIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicViscousIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const char *epsilon,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicViscousIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
			const T& /* EpsPrime */ = 0.);
};

typedef SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal> 
	SymbolicViscousIsotropicConstitutiveLaw1D;
typedef SymbolicViscousIsotropicConstitutiveLaw<Vec3, Mat3x3>
	SymbolicViscousIsotropicConstitutiveLaw3D;
typedef SymbolicViscousIsotropicConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicViscousIsotropicConstitutiveLaw6D;

template <class T, class Tder>
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::SymbolicViscousIsotropicConstitutiveLaw(
		const TplDriveCaller<T>* pDC, const T& PStress, 
		const char *epsilon, const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{ 
	typedef typename ConstitutiveLaw<T, Tder>::ErrNotAvailable E;
	THROW(E(std::cerr, "symbolic constitutive law "
				"is allowed only for scalar")); 
}
 
template <class T, class Tder>
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::~SymbolicViscousIsotropicConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return NULL;
}

template <class T, class Tder> std::ostream& 
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
SymbolicViscousIsotropicConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

/* SymbolicViscousIsotropicConstitutiveLaw - end */

/* SymbolicViscoElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicViscoElasticIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicViscoElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			const char *epsilon, const char *epsilonPrime,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicViscoElasticIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
			const T& /* EpsPrime */ = 0.);
};

typedef SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal> 
	SymbolicViscoElasticIsotropicConstitutiveLaw1D;
typedef SymbolicViscoElasticIsotropicConstitutiveLaw<Vec3, Mat3x3>
	SymbolicViscoElasticIsotropicConstitutiveLaw3D;
typedef SymbolicViscoElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicViscoElasticIsotropicConstitutiveLaw6D;

template <class T, class Tder>
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::SymbolicViscoElasticIsotropicConstitutiveLaw(
		const TplDriveCaller<T>* pDC, const T& PStress, 
		const char *epsilon,  const char *epsilonPrime,
		const char *expression, const char *const *symbols = NULL)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{ 
	typedef typename ConstitutiveLaw<T, Tder>::ErrNotAvailable E;
	THROW(E(std::cerr, "symbolic constitutive law "
				"is allowed only for scalar")); 
}
 
template <class T, class Tder>
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::~SymbolicViscoElasticIsotropicConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return NULL;
}

template <class T, class Tder> std::ostream& 
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& EpsPrime )
{
	NO_OP;
}

/* SymbolicViscoElasticIsotropicConstitutiveLaw - end */

#endif /* SYMCLTP_H */

