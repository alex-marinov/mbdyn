/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Giuseppe Quaranta <quaranta@aero.polimi.it>
 *
 * Based on:	Gandhi, F. and Chopra, I., "An analytical model for
 * 		a nonlinear elastomeric lag damper and its effect on
 *		aeromechanical stability in hover"
 *		Journal of the American Helicopter Society,
 *		Vol. 39, 1994, pp. 59-69.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

#include "dataman.h"
#include "constltp.h"

class DamperGandhiConstitutiveLaw;

class DamperGandhi
: virtual public Elem, public UserDefinedElem {
	friend class DamperGandhiConstitutiveLaw;

private:
	doublereal m_K_1;
	doublereal m_C_1;

	doublereal m_c_1;
	doublereal m_c_2;
	doublereal m_c_3;
	doublereal m_c_4;

	doublereal m_f;
	doublereal m_dot_f;

	doublereal m_g;
	doublereal m_dg_df;

	const DamperGandhiConstitutiveLaw *m_pCL;

	void SetDGCL(const DamperGandhiConstitutiveLaw *pCL);
	void Update_int(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
	void Get(doublereal& f, doublereal& fde, doublereal& fdeprime) const;

public:
	DamperGandhi(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~DamperGandhi(void);

	virtual unsigned int iGetNumDof(void) const;
	virtual std::ostream& DescribeDof(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const;
	virtual void DescribeDof(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const;
	virtual std::ostream& DescribeEq(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const;
	virtual void DescribeEq(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	virtual void
	Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
	virtual void Output(OutputHandler& OH) const;
	unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
};

class DamperGandhiConstitutiveLaw
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	DamperGandhi *m_pElem;

public:
	DamperGandhiConstitutiveLaw(DamperGandhi *pElem)
	: m_pElem(pElem)
	{
		ASSERT(m_pElem != 0);
		m_pElem->SetDGCL(this);
		Update(0., 0.);
	};

	virtual ~DamperGandhiConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef DamperGandhiConstitutiveLaw cl;
		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_pElem));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "damper Gandhi, " << m_pElem->GetLabel();
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		m_pElem->Get(ConstitutiveLaw<doublereal, doublereal>::F,
				ConstitutiveLaw<doublereal, doublereal>::FDE,
				ConstitutiveLaw<doublereal, doublereal>::FDEPrime);
	};
};

/* specific functional object(s) */
struct DamperGandhiCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cout("DamperGandhiConstitutiveLaw\n        \"damper gandhi\" , <damper_gandhi_elem_label>" << std::endl);
		}

		DamperGandhi *pElem = dynamic_cast<DamperGandhi *>(const_cast<DataManager *>(pDM)->ReadElem(HP, Elem::LOADABLE));

		typedef DamperGandhiConstitutiveLaw L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pElem));

		return pCL;
	};
};

void
DamperGandhi::SetDGCL(const DamperGandhiConstitutiveLaw *pCL)
{
	m_pCL = pCL;
}

void
DamperGandhi::Update_int(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	integer iIndex = iGetFirstIndex() + 1;
	m_f = XCurr(iIndex);
	m_dot_f = XPrimeCurr(iIndex);

	doublereal abs_f = std::abs(m_f);

	m_g = m_f*(m_c_1 + abs_f*(m_c_2 + abs_f*(m_c_3 + abs_f*m_c_4)));
	m_dg_df = m_c_1 + abs_f*(2*m_c_2 + abs_f*(3*m_c_3 + abs_f*4*m_c_4));
}

void
DamperGandhi::Get(doublereal& f, doublereal& fde, doublereal& fdeprime) const
{
	f = m_f;

	// rough estimate...
	doublereal d = 1. + m_K_1*m_dg_df;
	fde = m_K_1/d;
	fdeprime = m_C_1/d;
}

DamperGandhi::DamperGandhi(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_f(0.),
m_dot_f(0.),
m_g(0.),
m_dg_df(0.),
m_pCL(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        damper Gandhi\n"
"Author:        Pierangelo Masarati <masarati@aero.polimi.it>\n"
"		Giuseppe Quaranta <quaranta@aero.polimi.it>\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               <http://www.aero.polimi.it/>\n"
"\n"
"               All rights reserved\n"
"\n"
"               user defined : <label> , <K1> , <C1> , <c1> , <c2> , <c3> , <c4> ;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_K_1 = HP.GetReal();
	m_C_1 = HP.GetReal();

	m_c_1 = HP.GetReal();
	m_c_2 = HP.GetReal();
	m_c_3 = HP.GetReal();
	m_c_4 = HP.GetReal();

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

DamperGandhi::~DamperGandhi(void)
{
	NO_OP;
}

void
DamperGandhi::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << GetLabel() 
			<< " " << m_f
			<< " " << m_dot_f
			<< std::endl;
	}
}

void
DamperGandhi::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Update_int(X, XP);
}

unsigned int
DamperGandhi::iGetNumDof(void) const
{
	return 1;
}

std::ostream&
DamperGandhi::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix << iGetFirstIndex() + 1 << ": damper force" << std::endl;
}

void
DamperGandhi::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT((i == -1) || (i == 1));

	std::ostringstream os;
	os << "DamperGandhi(" << GetLabel() << "): damper force";
	desc.resize(1);
	desc[0] = os.str();
}

std::ostream&
DamperGandhi::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix << iGetFirstIndex() + 1 << ": damper constitutive equation" << std::endl;
}

void
DamperGandhi::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT((i == -1) || (i == 1));

	std::ostringstream os;
	os << "DamperGandhi(" << GetLabel() << "): damper constitutive equation";
	desc.resize(1);
	desc[0] = os.str();
}

DofOrder::Order
DamperGandhi::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

DofOrder::Order
DamperGandhi::GetEqType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
DamperGandhi::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
DamperGandhi::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.ResizeReset(1, 1);

	integer iIndex = iGetFirstIndex() + 1;
	WM.PutRowIndex(1, iIndex);
	WM.PutColIndex(1, iIndex);

	WM.PutCoef(1, 1, m_C_1*m_dg_df + (1. + m_K_1*m_dg_df)*dCoef);

	return WorkMat;
}

SubVectorHandler& 
DamperGandhi::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(1);

	integer iIndex = iGetFirstIndex() + 1;

	doublereal xi = m_pCL->GetEpsilon();
	doublereal dot_xi = m_pCL->GetEpsilonPrime();

	WorkVec.PutItem(1, iIndex, m_K_1*(xi - m_g) + m_C_1*(dot_xi - m_dg_df*m_dot_f) - m_f);

	return WorkVec;
}

void
DamperGandhi::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	Update_int(XCurr, XPrimeCurr);
}

unsigned int
DamperGandhi::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
DamperGandhi::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "f") == 0) {
		return 1;
	}

	return 0;
}

doublereal
DamperGandhi::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
		return m_f;
	}

	return 0.;
}

int
DamperGandhi::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
DamperGandhi::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
DamperGandhi::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iIndex = iGetFirstIndex() + 1;
	X(iIndex) = m_f;
	XP(iIndex) = m_dot_f;
}

std::ostream&
DamperGandhi::Restart(std::ostream& out) const
{
	// don't worry about "soft" restart by now
	return out << "# DamperGandhi: not implemented" << std::endl;
}

unsigned int
DamperGandhi::iGetInitialNumDof(void) const
{
	return 0;
}

void 
DamperGandhi::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
DamperGandhi::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
DamperGandhi::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<DamperGandhi>;

	if (!SetUDE("damper" "gandhi", rf)) {
		delete rf;

		silent_cerr("module-damper-gandhi: "
			"module_init(" << module_name << ") "
			"unable to register \"damper gandhi\" user-defined element" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new DamperGandhiCLR;
	if (!SetCL1D("damper" "gandhi", rf1D)) {
		delete rf1D;

		silent_cerr("module-damper-gandhi: "
			"module_init(" << module_name << ") "
			"unable to register \"damper gandhi\" 1D constitutive law" << std::endl);

		return -1;
	}

	return 0;
}

