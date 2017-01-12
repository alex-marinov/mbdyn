/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "constltp.h"

/* SymbolicElasticConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicElasticConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicElasticConstitutiveLaw(
		const TplDriveCaller<T>* pDC,
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& expression);
     	virtual ~SymbolicElasticConstitutiveLaw(void);
};

typedef SymbolicElasticConstitutiveLaw<doublereal, doublereal> 
	SymbolicElasticConstitutiveLaw1D;
typedef SymbolicElasticConstitutiveLaw<Vec3, Mat3x3>
	SymbolicElasticConstitutiveLaw3D;
typedef SymbolicElasticConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicElasticConstitutiveLaw6D;

template <class T, class Tder>
SymbolicElasticConstitutiveLaw<T, Tder>::SymbolicElasticConstitutiveLaw(
	const TplDriveCaller<T>* pDC,
	const T& PStress,
	std::vector<std::string>& epsilon,
	std::vector<std::string>& expression)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{
	NO_OP;
}

template <class T, class Tder>
SymbolicElasticConstitutiveLaw<T, Tder>::~SymbolicElasticConstitutiveLaw(void)
{
	NO_OP;
}

/* SymbolicElasticConstitutiveLaw - end */

/* SymbolicViscousConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicViscousConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicViscousConstitutiveLaw(
		const T& PStress,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~SymbolicViscousConstitutiveLaw(void);
	virtual ConstLawType::Type GetConstLawType(void) const;
};

typedef SymbolicViscousConstitutiveLaw<doublereal, doublereal> 
	SymbolicViscousConstitutiveLaw1D;
typedef SymbolicViscousConstitutiveLaw<Vec3, Mat3x3>
	SymbolicViscousConstitutiveLaw3D;
typedef SymbolicViscousConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicViscousConstitutiveLaw6D;

template <class T, class Tder>
SymbolicViscousConstitutiveLaw<T, Tder>::SymbolicViscousConstitutiveLaw(
	const T& PStress,
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: ElasticConstitutiveLaw<T, Tder>(0, PStress)
{
	NO_OP;
}

template <class T, class Tder>
SymbolicViscousConstitutiveLaw<T, Tder>::~SymbolicViscousConstitutiveLaw(void)
{
	NO_OP;
}

template <class T, class Tder>
ConstLawType::Type
SymbolicViscousConstitutiveLaw<T, Tder>::GetConstLawType(void) const
{
	return ConstLawType::VISCOUS;
}

/* SymbolicViscousConstitutiveLaw - end */

/* SymbolicViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicViscoElasticConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicViscoElasticConstitutiveLaw(
		const TplDriveCaller<T>* pDC,
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~SymbolicViscoElasticConstitutiveLaw(void);
	virtual ConstLawType::Type GetConstLawType(void) const;
};

typedef SymbolicViscoElasticConstitutiveLaw<doublereal, doublereal> 
	SymbolicViscoElasticConstitutiveLaw1D;
typedef SymbolicViscoElasticConstitutiveLaw<Vec3, Mat3x3>
	SymbolicViscoElasticConstitutiveLaw3D;
typedef SymbolicViscoElasticConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicViscoElasticConstitutiveLaw6D;

template <class T, class Tder>
SymbolicViscoElasticConstitutiveLaw<T, Tder>::SymbolicViscoElasticConstitutiveLaw(
	const TplDriveCaller<T>* pDC,
	const T& PStress,
	std::vector<std::string>& epsilon,
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{
	NO_OP;
}

template <class T, class Tder>
SymbolicViscoElasticConstitutiveLaw<T, Tder>::~SymbolicViscoElasticConstitutiveLaw(void)
{
	NO_OP;
}

template <class T, class Tder>
ConstLawType::Type
SymbolicViscoElasticConstitutiveLaw<T, Tder>::GetConstLawType(void) const
{
	return ConstLawType::VISCOELASTIC;
}

/* SymbolicViscoElasticConstitutiveLaw - end */

#endif /* SYMCLTP_H */

