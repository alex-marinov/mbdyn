/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati <pierangelo.masarati@polimi.it>
 *
 * Sponsored by Hutchinson CdR, 2018-2019
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

#include "hfelem.h"

class HarmonicForcingElem
: virtual public Elem, public UserDefinedElem {
private:
	// add private data
	enum Priv {
		DT = 1,
		F,
		OMEGA,
		COUNT
	};

	const DataManager* m_pDM;

	integer m_iN;	// must be even

	doublereal m_dTInit;
	doublereal m_dT0;
	doublereal m_dDeltaT;
	doublereal m_dOmega;
	doublereal m_dOmega0;
	doublereal m_dOmegaMax;
	doublereal m_dF;

	doublereal m_dTol;
	integer m_iMinPeriods;
	doublereal m_dOmegaAddInc;
	doublereal m_dOmegaMulInc;

	integer m_iOmegaCnt;
	integer m_iPeriod;
	integer m_iPeriodCnt;
	bool m_bStarted;
	bool m_bConverged;

	bool m_bOut;
	integer m_iPeriodOut;
	doublereal m_dOmegaOut;

	std::vector<doublereal> m_cos_psi;
	std::vector<doublereal> m_sin_psi;
	std::vector<std::vector<doublereal> > m_X;
	std::vector<doublereal> m_Xcos;
	std::vector<doublereal> m_Xsin;
	std::vector<doublereal> m_XcosPrev;
	std::vector<doublereal> m_XsinPrev;

	std::vector<DriveCaller *> m_Input;

public:
	HarmonicForcingElem(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~HarmonicForcingElem(void);

	virtual void Output(OutputHandler& OH) const;
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
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
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
};

HarmonicForcingElem::HarmonicForcingElem(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pDM(pDM),
m_dTInit(-std::numeric_limits<doublereal>::max()),
m_dF(0.),
m_dOmegaAddInc(0.),
m_dOmegaMulInc(1.),
m_iOmegaCnt(0),
m_iPeriod(0),
m_iPeriodCnt(0),
m_bStarted(false),
m_bConverged(false),
m_bOut(false),
m_iPeriodOut(0),
m_dOmegaOut(0.)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	HarmonicExcitationElem					\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// do something useful

	// inputs
	if (!HP.IsKeyWord("inputs" "number")) {
		silent_cerr("HarmonicExcitationElem: \"inputs number\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNInput;
	try {
		iNInput = HP.GetInt(1, HighParser::range_gt<integer>(0));
	
	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("HarmonicExcitationElem: number of input must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	m_Input.resize(iNInput);
	for (integer i = 0; i < iNInput; ++i) {
		m_Input[i] = HP.GetDriveCaller();
	}

	// block size
	if (!HP.IsKeyWord("steps" "number")) {
		silent_cerr("HarmonicExcitationElem: \"steps number\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_iN = HP.GetInt(2, HighParser::range_gt<integer>(0));

	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("HarmonicExcitationElem: N must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	if ((m_iN/2)*2 != m_iN) {
		silent_cerr("HarmonicExcitationElem: N must be even, at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// frequency
	if (!HP.IsKeyWord("initial" "angular" "frequency")) {
		silent_cerr("HarmonicExcitationElem: \"initial angular frequency\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_dOmega0 = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

	} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
		silent_cerr("HarmonicExcitationElem: initial frequency must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	if (HP.IsKeyWord("max" "angular" "frequency")) {
		if (HP.IsKeyWord("forever")) {
			m_dOmegaMax = std::numeric_limits<doublereal>::max();

		} else {
			try {
				m_dOmegaMax = HP.GetReal(m_dOmega0, HighParser::range_gt<doublereal>(m_dOmega0));

			} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
				silent_cerr("HarmonicExcitationElem: maximum frequency must be greater than initial frequency " << m_dOmega0 << ", at line " << HP.GetLineData() << std::endl);
				throw e;
			}
		}

	} else {
		m_dOmegaMax = std::numeric_limits<doublereal>::max();
	}

	m_dOmega = m_dOmega0;

	m_dDeltaT = (2*M_PI/m_dOmega)/m_iN;

	// frequency increment
	if (!HP.IsKeyWord("angular" "frequency" "increment")) {
		silent_cerr("HarmonicExcitationElem: \"angular frequency increment\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("additive")) {
		try {
			m_dOmegaAddInc = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

		} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
			silent_cerr("HarmonicExcitationElem: frequency additive increment must be positive, at line " << HP.GetLineData() << std::endl);
			throw e;
		}

	} else if (HP.IsKeyWord("multiplicative")) {
		try {
			m_dOmegaMulInc = HP.GetReal(1., HighParser::range_gt<doublereal>(1.));

		} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
			silent_cerr("HarmonicExcitationElem: frequency multiplicative increment must be greater than 1, at line " << HP.GetLineData() << std::endl);
			throw e;
		}

	} else {
		silent_cerr("HarmonicExcitationElem: unknown or missing frequency increment method at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// initial time
	if (HP.IsKeyWord("initial" "time")) {
		m_dTInit = HP.GetReal();
	}

	// tolerance
	if (!HP.IsKeyWord("tolerance")) {
		silent_cerr("HarmonicExcitationElem: \"tolerance\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_dTol = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

	} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
		silent_cerr("HarmonicExcitationElem: tolerance must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	// minimum number of blocks
	if (!HP.IsKeyWord("min" "periods")) {
		silent_cerr("HarmonicExcitationElem: \"min periods\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_iMinPeriods = HP.GetInt(2, HighParser::range_gt<integer>(1));

	} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
		silent_cerr("HarmonicExcitationElem: minimum number of periods must be greater than 1, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	m_cos_psi.resize(m_iN);
	m_sin_psi.resize(m_iN);
	m_X.resize(m_iN);

	m_Xcos.resize(iNInput);
	m_Xsin.resize(iNInput);
	m_XcosPrev.resize(iNInput);
	m_XsinPrev.resize(iNInput);

	for (integer i = 0; i < m_iN; i++) {
		doublereal dPsi = (2*M_PI*i)/m_iN;
		m_cos_psi[i] = cos(dPsi)/(m_iN/2);
		m_sin_psi[i] = sin(dPsi)/(m_iN/2);

		m_X[i].resize(iNInput);
	}

	// reset
	for (unsigned i = 0; i < m_Xcos.size(); ++i) {
		m_Xcos[i] = 0.;
		m_Xsin[i] = 0.;
	}
}

HarmonicForcingElem::~HarmonicForcingElem(void)
{
	// destroy private data
	for (std::vector<DriveCaller *>::iterator i = m_Input.begin(); i != m_Input.end(); ++i) {
		delete *i;
	}
}

void
HarmonicForcingElem::Output(OutputHandler& OH) const
{
	// NOTE: output occurs only at convergence,
	// not with fixed periodicity
	if (bToBeOutput()) {
		if (m_bOut && m_iPeriod == 0 && m_iPeriodCnt == 0) {
			std::ostream& out = OH.Loadable();

			out << std::setw(8) << GetLabel()	// 1:	label
				<< " " << m_pDM->dGetTime()	// 2:	when
				<< " " << m_dOmegaOut		// 3:	what frequency
				<< " " << m_iPeriodOut;		// 4:	how many periods required

			for (unsigned i = 0; i < m_Input.size(); ++i) {
				out << " " << m_Xsin[i] << " " << m_Xcos[i];
			}

			out << std::endl;
		}
	}

	// TODO: NetCDF support
}

void
HarmonicForcingElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
HarmonicForcingElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
HarmonicForcingElem::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkVec.ResizeReset(0);

	return WorkVec;
}

void
HarmonicForcingElem::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	// initialize if needed
	if (!m_bStarted) {
		if (m_pDM->dGetTime() >= m_dTInit) {
			m_bStarted = true;
			m_dT0 = m_pDM->dGetTime();

		}

	} else {

		// if just initialized
		if (m_iPeriod == 0 && m_iPeriodCnt == 0) {
			m_dT0 = m_pDM->dGetTime();

			// reset vectors
			for (unsigned i = 0; i < m_Input.size(); ++i) {
				m_Xcos[i] = 0.;
				m_Xsin[i] = 0.;
			}

			m_bOut = true;
		}
	}

	if (m_bStarted) {
		m_dF = sin(m_dOmega*(m_pDM->dGetTime() - m_dT0));
	}
}

void
HarmonicForcingElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	if (!m_bStarted) {
		return;
	}

	// collect input
	doublereal dErr = 0.;

	for (unsigned i = 0; i < m_Input.size(); ++i) {
		if (m_iPeriod > 0) {
			m_XcosPrev[i] = m_Xcos[i];
			m_XsinPrev[i] = m_Xsin[i];

			m_Xcos[i] -= m_X[m_iPeriodCnt][i]*m_cos_psi[m_iPeriodCnt];
			m_Xsin[i] -= m_X[m_iPeriodCnt][i]*m_sin_psi[m_iPeriodCnt];
		}

		m_X[m_iPeriodCnt][i] = m_Input[i]->dGet();

		m_Xcos[i] += m_X[m_iPeriodCnt][i]*m_cos_psi[m_iPeriodCnt];
		m_Xsin[i] += m_X[m_iPeriodCnt][i]*m_sin_psi[m_iPeriodCnt];

		if (m_iPeriod > 0) {
			doublereal d = m_Xcos[i] - m_XcosPrev[i];
			dErr += d*d;
			d = m_Xsin[i] - m_XsinPrev[i];
			dErr += d*d;
		}
	}

	// check for convergence
	if (m_iPeriod > 0) {
		dErr = sqrt(dErr);
		if (dErr <= m_dTol) {
			m_bConverged = true;

		} else {
			m_bConverged = false;
		}
	}

	++m_iPeriodCnt;
	if (m_iPeriodCnt == m_iN) {
		m_iPeriodCnt = 0;
		++m_iPeriod;
		if (m_bConverged && m_iPeriod >= m_iMinPeriods) {
			m_iPeriodOut = m_iPeriod;
			m_dOmegaOut = m_dOmega;
			m_iOmegaCnt++;

			m_bConverged = false;
			m_iPeriod = 0;

			if (m_dOmegaAddInc > 0.) {
				m_dOmega += m_dOmegaAddInc;

			} else {
				m_dOmega *= m_dOmegaMulInc;
			}

			if (m_dOmega >= m_dOmegaMax) {
				// schedule for termination!
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}

			m_dDeltaT = (2*M_PI/m_dOmega)/m_iN;
		}
	}
}

unsigned int
HarmonicForcingElem::iGetNumPrivData(void) const
{
	return 4;
}

unsigned int
HarmonicForcingElem::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "timestep") == 0) {
		return DT;

	} else if (strcmp(s, "excitation") == 0) {
		return F;

	} else if (strcmp(s, "omega") == 0) {
		return OMEGA;

	} else if (strcmp(s, "count") == 0) {
		return COUNT;
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
HarmonicForcingElem::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case DT:
		return m_dDeltaT;

	case F:
		return m_dF;

	case OMEGA:
		return m_dOmega;

	case COUNT:
		return doublereal(m_iOmegaCnt);

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

int
HarmonicForcingElem::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
HarmonicForcingElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
HarmonicForcingElem::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
HarmonicForcingElem::Restart(std::ostream& out) const
{
	return out << "# HarmonicForcingElem: not implemented" << std::endl;
}

unsigned int
HarmonicForcingElem::iGetInitialNumDof(void) const
{
	return 0;
}

void 
HarmonicForcingElem::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
HarmonicForcingElem::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
HarmonicForcingElem::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool
hfelem_set(void)
{
	UserDefinedElemRead *rf = new UDERead<HarmonicForcingElem>;

	if (!SetUDE("harmonic" "excitation", rf)) {
		delete rf;

		silent_cerr("hfdrive: hfelem_set failed" << std::endl);

		return false;
	}

	return true;
}

