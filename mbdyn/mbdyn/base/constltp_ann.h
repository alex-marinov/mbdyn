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

#ifndef CONSTLAW_ANN_H
#define CONSTLAW_ANN_H

#include "dataman.h"
#include "constltp.h"
#include "ann.h"

template <class T, class Tder>
class AnnElasticConstitutiveLaw
: public ConstitutiveLaw<T, Tder> {
protected:
	bool bUnit;
	unsigned dim;
	ANN *net;
	std::string fname;

	void AnnInit(void)
	{
		if (typeid(T) == typeid(Vec3)) {
			dim = 3;

		} else if (typeid(T) == typeid(Vec6)) {
			dim = 6;

		} else {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		net = new ANN;
		if (ANN_init(net, fname.c_str()) != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	virtual void
	AnnSanity(void)
	{
		/* add sanity checks */
		NO_OP;
	}

	virtual void Update(const T& Eps, const T& EpsPrime, int feedback) {

		ConstitutiveLaw<T, Tder>::Epsilon = Eps;

		T E = ConstitutiveLaw<T, Tder>::Epsilon;

		ConstitutiveLaw<T, Tder>::F = 0.;
		ConstitutiveLaw<T, Tder>::FDE = 0.;

		for (unsigned r = 0; r < dim; r++) {
                	net->input.vec[r] = E(r + 1);
		}
                if (ANN_sim(net, &net->input, &net->output, feedback)) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
		ANN_jacobian_matrix( net, &net->jacobian );

		for (unsigned r = 0; r < dim; r++) {
			ConstitutiveLaw<T, Tder>::F(r + 1) = -net->output.vec[r];
			for (unsigned c = 0; c < dim; c++) {
				ConstitutiveLaw<T, Tder>::FDE(r + 1, c + 1) = -net->jacobian.mat[r][c];
			}
		}
	};


public:
	AnnElasticConstitutiveLaw(const std::string& f, bool b = false)
	: bUnit(b), fname(f)
	{
		AnnInit();
		AnnSanity();
	}

	virtual ~AnnElasticConstitutiveLaw(void) {
		if (net != 0) {
			(void)ANN_destroy(net);
			delete net;
		}
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const T& Eps, const T& EpsPrime = 0.)
	{
		Update(Eps, EpsPrime, ANN_FEEDBACK_UPDATE);
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef AnnElasticConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(fname));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << " ann elastic, \"" << fname << "\",";
	};

	virtual void Update(const T& Eps, const T& EpsPrime = 0.) {
		Update(Eps, EpsPrime, ANN_FEEDBACK_NONE);
	};
};

template <>
class AnnElasticConstitutiveLaw<doublereal, doublereal>
: public ConstitutiveLaw<doublereal, doublereal> {
protected:
	bool bUnit;
	ANN *net;
	std::string fname;

	void AnnInit(void)
	{
		net = new ANN;
		if (ANN_init(net, fname.c_str()) != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	virtual void
	AnnSanity(void)
	{
		/* add sanity checks */
		NO_OP;
	}

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime, int feedback) {

		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

		doublereal E = ConstitutiveLaw<doublereal, doublereal>::Epsilon;

		ConstitutiveLaw<doublereal, doublereal>::F = 0.;
		ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;

               	net->input.vec[0] = E;
		if (net->N_input == 2) {
			if (bUnit) {
				net->input.vec[1] = 1.;
			} else {
				net->input.vec[1] = 0.;
			}
		}

                if (ANN_sim(net, &net->input, &net->output, feedback)) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
		ANN_jacobian_matrix(net, &net->jacobian);

		ConstitutiveLaw<doublereal, doublereal>::F = -net->output.vec[0];
		ConstitutiveLaw<doublereal, doublereal>::FDE = -net->jacobian.mat[0][0];
	};

public:
	AnnElasticConstitutiveLaw(const std::string& f, bool b = false)
	: bUnit(b), fname(f)
	{
		AnnInit();
		AnnSanity();
	}

	virtual ~AnnElasticConstitutiveLaw(void) {
		if (net != 0) {
			(void)ANN_destroy(net);
			delete net;
		}
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const doublereal& Eps, const doublereal& EpsPrime = 0.)
	{
		Update(Eps, EpsPrime, ANN_FEEDBACK_UPDATE);
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef AnnElasticConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(fname));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << " ann elastic, \"" << fname << "\",";
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		Update(Eps, EpsPrime, ANN_FEEDBACK_NONE);
	};
};

template <class T, class Tder>
class AnnViscoElasticConstitutiveLaw
: public AnnElasticConstitutiveLaw<T, Tder> {
protected:
	/* overrides AnnElasticConstitutiveLaw::AnnSanity() */
	virtual void
	AnnSanity(void)
	{
		/* add sanity checks */
		NO_OP;
	}

	virtual void Update(const T& Eps, const T& EpsPrime, int feedback) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		T E = ConstitutiveLaw<T, Tder>::Epsilon;

		ConstitutiveLaw<T, Tder>::F = 0.;
		ConstitutiveLaw<T, Tder>::FDE = 0.;

		ANN *net = AnnElasticConstitutiveLaw<T, Tder>::net;

		/* loop according to dimensions... */
		for (unsigned r = 0; r < AnnViscoElasticConstitutiveLaw<T, Tder>::dim; r++) {
			net->input.vec[r] = E(r + 1);
			net->input.vec[AnnViscoElasticConstitutiveLaw<T, Tder>::dim + r] = EpsPrime(r + 1);
		}
                if (ANN_sim(net, &net->input, &net->output, feedback))
		{
                        silent_cout("AnnViscoElasticConstitutiveLaw: Network simulation error" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
		ANN_jacobian_matrix(net, &net->jacobian);

		for (unsigned r = 0; r < AnnViscoElasticConstitutiveLaw<T, Tder>::dim; r++) {
			ConstitutiveLaw<T, Tder>::F(r + 1) = -net->output.vec[r];
			for (unsigned c = 0; c < AnnViscoElasticConstitutiveLaw<T, Tder>::dim; c++) {
				ConstitutiveLaw<T, Tder>::FDE(r + 1, c + 1) = -net->jacobian.mat[r][c];
			}
		}
	};

public:
	AnnViscoElasticConstitutiveLaw(const std::string& f, bool b = false)
	: AnnElasticConstitutiveLaw<T, Tder>(f, b)
	{
		NO_OP;
	};

	virtual ~AnnViscoElasticConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef AnnViscoElasticConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(cl::fname));
		return pCL;
	};
};

template <>
class AnnViscoElasticConstitutiveLaw<doublereal, doublereal>
: public AnnElasticConstitutiveLaw<doublereal, doublereal> {
protected:
	/* overrides AnnElasticConstitutiveLaw::AnnSanity() */
	virtual void
	AnnSanity(void)
	{
		/* add sanity checks */
		NO_OP;
	}

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime, int feedback ) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal E = ConstitutiveLaw<doublereal, doublereal>::Epsilon;

		ANN *net = AnnElasticConstitutiveLaw<doublereal, doublereal>::net;
		net->input.vec[0] = E;
		net->input.vec[1] = EpsPrime;
		if (net->N_input == 3) {
			if (bUnit) {
				net->input.vec[2] = 1.;
			} else {
				net->input.vec[2] = 0.;
			}
		}
                if (ANN_sim(net, &net->input, &net->output, feedback))
		{
                        silent_cout("AnnViscoElasticConstitutiveLaw: Network simulation error" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
		ANN_jacobian_matrix(net, &net->jacobian );
	
		ConstitutiveLaw<doublereal, doublereal>::F = -net->output.vec[0];
		ConstitutiveLaw<doublereal, doublereal>::FDE = -net->jacobian.mat[0][0];
	};

public:
	AnnViscoElasticConstitutiveLaw(const std::string& f, bool b = false)
	: AnnElasticConstitutiveLaw<doublereal, doublereal>(f, b)
	{
		NO_OP;
	};

	virtual ~AnnViscoElasticConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef AnnViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(cl::fname));
		return pCL;
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct AnnElasticCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::ELASTIC;

		bool bUnit(false);
		if (HP.IsKeyWord("unit" "input")) {
			bUnit = HP.GetYesNoOrBool();
		}

		const char *s = HP.GetFileName();
		if (s == 0) {
			silent_cerr("AnnElasticCLR: "
				"unable to get ann file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		typedef AnnElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(s, bUnit));

		return pCL;
	};
};

template <class T, class Tder>
struct AnnViscoElasticCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		bool bUnit(false);
		if (HP.IsKeyWord("unit" "input")) {
			bUnit = HP.GetYesNoOrBool();
		}

		const char *s = HP.GetFileName();
		if (s == 0) {
			silent_cerr("AnnViscoElasticCLR: "
				"unable to get ann file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		typedef AnnViscoElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(s, bUnit));

		return pCL;
	};
};

#endif // CONSTLAW_ANN_H
