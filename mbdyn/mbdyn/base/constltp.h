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

/* Legami costitutivi */

#ifndef CONSTLTP_H
#define CONSTLTP_H

#include "withlab.h"
#include "simentity.h"
#include "tpldrive.h"

#ifdef USE_AUTODIFF
#include "matvec.h"
#endif

/* Tipi di cerniere deformabili */
class ConstLawType {
public:
	enum Type {
		UNKNOWN = 0x0,
		
		ELASTIC = 0x1,
		VISCOUS = 0x2,
		VISCOELASTIC = (ELASTIC | VISCOUS),
	
		LASTCONSTLAWTYPE
	};
};

/* ConstitutiveLaw - begin */
template <class T, class Tder>
class ConstitutiveLaw : public WithLabel, public SimulationEntity {
public:
	class ErrNotAvailable : public MBDynErrBase {
	public:
		ErrNotAvailable(MBDYN_EXCEPT_ARGS_DECL)
		: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{
			silent_cerr("Constitutive law not available "
				"for this dimensionality"
				<< std::endl);
		};
		ErrNotAvailable(std::ostream& out, MBDYN_EXCEPT_ARGS_DECL)
		: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{
			out << "Constitutive law not available "
				"for this dimensionality"
				<< what()
				<< std::endl;
		};
	};

	typedef typename ConstitutiveLaw<T, Tder>::ErrNotAvailable Err;   
   
protected:
	T Epsilon;		/* strain */
	T EpsilonPrime;		/* strain rate */

	T F;			/* force */
	Tder FDE;		/* force strain derivative */
	Tder FDEPrime;		/* force strain rate derivative */

public:
	ConstitutiveLaw(void)
	: WithLabel(0),
	Epsilon(mb_zero<T>()), EpsilonPrime(mb_zero<T>()),
	F(mb_zero<T>()), FDE(mb_zero<Tder>()), FDEPrime(mb_zero<Tder>()) {
		NO_OP;
	};

	virtual ~ConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const = 0;
	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};
	
	virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) = 0;

	virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		NO_OP;
	};
   
	virtual const T& GetEpsilon(void) const {
		return Epsilon;
	};

	virtual const T& GetEpsilonPrime(void) const {
		return EpsilonPrime;
	};

	virtual const T& GetF(void) const {
		return F;
	};

	virtual const Tder& GetFDE(void) const {
		return FDE;
	};

	virtual const Tder& GetFDEPrime(void) const {
		return FDEPrime;
	};

	/* simentity */
	virtual unsigned int iGetNumDof(void) const {
		return 0;
	};

	virtual std::ostream& DescribeDof(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const
	{
		return out;
	};
	virtual void DescribeDof(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const
	{
		ASSERT(i <= 0);
		desc.resize(0);
	};

	virtual std::ostream& DescribeEq(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const
	{
		return out;
	};

	virtual void DescribeEq(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const
	{
		ASSERT(i <= 0);
		desc.resize(0);
	};

	virtual DofOrder::Order GetDofType(unsigned int i) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

protected:
#ifdef USE_AUTODIFF
	template <typename ConstLaw>
	static inline void UpdateViscoelastic(ConstLaw* pCl, const T& Eps, const T& EpsPrime);

	template <typename ConstLaw>
	static inline void UpdateElastic(ConstLaw* pCl, const T& Eps);
#endif
};

typedef ConstitutiveLaw<doublereal, doublereal> ConstitutiveLaw1D;
typedef ConstitutiveLaw<Vec3, Mat3x3> ConstitutiveLaw3D;
typedef ConstitutiveLaw<Vec6, Mat6x6> ConstitutiveLaw6D;

#ifdef USE_AUTODIFF
template <class T, class Tder>
template <typename ConstLaw>
inline void ConstitutiveLaw<T, Tder>::UpdateViscoelastic(ConstLaw* pCl, const T& Eps, const T& EpsPrime)
{


	using namespace grad;
	const index_type N = VectorSize<T>::N;
	Vector<Gradient<2 * N>, N> gEps, gEpsPrime, gF;

	pCl->Epsilon = Eps;
	pCl->EpsilonPrime = EpsPrime;

	for (int i = 0; i < N; ++i)
	{
		gEps(i + 1).SetValuePreserve(Eps(i + 1));
		gEps(i + 1).DerivativeResizeReset(0, i, MapVectorBase::LOCAL, 1.);
		gEpsPrime(i + 1).SetValuePreserve(EpsPrime(i + 1));
		gEpsPrime(i + 1).DerivativeResizeReset(0, i + N, MapVectorBase::LOCAL, 1.);
	}

	pCl->UpdateViscoelasticTpl(gEps, gEpsPrime, gF);

	for (int i = 1; i <= N; ++i)
	{
		pCl->F(i) = gF(i).dGetValue();

		for (int j = 0; j < N; ++j)
		{
			pCl->FDE(i, j + 1) = gF(i).dGetDerivativeLocal(j);
			pCl->FDEPrime(i, j + 1) = gF(i).dGetDerivativeLocal(j + N);
		}
	}
}

template <class T, class Tder>
template <typename ConstLaw>
inline void ConstitutiveLaw<T, Tder>::UpdateElastic(ConstLaw* pCl, const T& Eps)
{


	using namespace grad;
	const index_type N = VectorSize<T>::N;
	Vector<Gradient<N>, N> gEps, gF;

	pCl->Epsilon = Eps;

	for (int i = 0; i < N; ++i)
	{
		gEps(i + 1).SetValuePreserve(Eps(i + 1));
		gEps(i + 1).DerivativeResizeReset(0, i, MapVectorBase::LOCAL, 1.);
	}

	pCl->UpdateElasticTpl(gEps, gF);

	for (int i = 1; i <= N; ++i)
	{
		pCl->F(i) = gF(i).dGetValue();

		for (int j = 0; j < N; ++j)
		{
			pCl->FDE(i, j + 1) = gF(i).dGetDerivativeLocal(j);
		}
	}
}
#endif

/* ConstitutiveLaw - end */


/* ConstitutiveLawOwner - begin */

template <class T, class Tder>
class ConstitutiveLawOwner : public SimulationEntity {
protected:
	mutable ConstitutiveLaw<T, Tder>* pConstLaw;
   
public:
	ConstitutiveLawOwner(const ConstitutiveLaw<T, Tder>* pCL)
	: pConstLaw(const_cast<ConstitutiveLaw<T, Tder> *>(pCL)) { 
		ASSERT(pCL != NULL);
	};
   
	virtual ~ConstitutiveLawOwner(void) {
		ASSERT(pConstLaw != NULL);
		if (pConstLaw != NULL) {
			SAFEDELETE(pConstLaw);
		}	
	};

	inline ConstitutiveLaw<T, Tder>* pGetConstLaw(void) const {
		ASSERT(pConstLaw != NULL);
		return pConstLaw; 
	};

	inline void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ASSERT(pConstLaw != NULL);
		pConstLaw->Update(Eps, EpsPrime);
	};

	inline void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ASSERT(pConstLaw != NULL);
		pConstLaw->AfterConvergence(Eps, EpsPrime);
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

	/* simentity */
	virtual unsigned int iGetNumDof(void) const {
		ASSERT(pConstLaw != NULL);
		return pConstLaw->iGetNumDof();
	};

	virtual std::ostream& DescribeDof(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const
	{
		return out;
	};

	virtual void DescribeDof(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const
	{
		ASSERT(i <= 0);
		desc.resize(0);
	};

	virtual std::ostream& DescribeEq(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const
	{
		return out;
	};

	virtual void DescribeEq(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const
	{
		ASSERT(i <= 0);
		desc.resize(0);
	};

	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(pConstLaw != NULL);
		return pConstLaw->GetDofType(i);
	};

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const {
		return pConstLaw->iGetNumPrivData();
	};

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return pConstLaw->iGetPrivDataIdx(s);
	};

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const {
		return pConstLaw->dGetPrivData(i);
	};

	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return pConstLaw->OutputAppend(out);
	};
};

typedef ConstitutiveLawOwner<doublereal, doublereal> ConstitutiveLaw1DOwner;
typedef ConstitutiveLawOwner<Vec3, Mat3x3> ConstitutiveLaw3DOwner;
typedef ConstitutiveLawOwner<Vec6, Mat6x6> ConstitutiveLaw6DOwner;

/* ConstitutiveLawOwner - end */

/* functions that read a constitutive law */
extern ConstitutiveLaw<doublereal, doublereal> *
ReadCL1D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
extern ConstitutiveLaw<Vec3, Mat3x3> *
ReadCL3D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
extern ConstitutiveLaw<Vec6, Mat6x6> *
ReadCL6D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

/* prototype of the template functional object: reads a constitutive law */
template <class T, class Tder>
struct ConstitutiveLawRead {
	virtual ~ConstitutiveLawRead<T, Tder>( void ) { NO_OP; };
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) = 0;
};

/* constitutive law registration functions: call to register one */
extern bool
SetCL1D(const char *name, ConstitutiveLawRead<doublereal, doublereal> *rf);
extern bool
SetCL3D(const char *name, ConstitutiveLawRead<Vec3, Mat3x3> *rf);
extern bool
SetCL6D(const char *name, ConstitutiveLawRead<Vec6, Mat6x6> *rf);

/* create/destroy */
extern void InitCL(void);
extern void DestroyCL(void);

#endif /* CONSTLTP_H */

