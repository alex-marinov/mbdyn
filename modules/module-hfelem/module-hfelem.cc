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

#include "mbpar.h"
#include "privdrive.h"
#include "drive_.h"

#include "module-hfelem.h"

static bool bHFElem(false);

class HarmonicForcingElem
: virtual public Elem, public UserDefinedElem, public DriveOwner {
private:
	// add private data
	enum Priv {
		Priv_DT = 1,
		Priv_F,
		Priv_PSI,
		Priv_OMEGA,
		Priv_AMPLITUDE,
		Priv_COUNT,

		Priv_LAST
	};

	const DataManager* m_pDM;

	integer m_iN;

	doublereal m_dTInit;
	doublereal m_dT0;
	doublereal m_dDeltaT;
	doublereal m_dMaxDeltaT;
	doublereal m_dOmega;
	doublereal m_dOmega0;
	doublereal m_dOmegaMax;
	doublereal m_dPsi;
	doublereal m_dF;
	doublereal m_dAmplitude;

	bool bRMSTest;
	bool bRMSTestTarget;
	doublereal m_dTol;
	integer m_iMinPeriods;
	integer m_iMaxPeriods;
	enum NoConv {
		NoConv_CONTINUE,
		NoConv_ABORT
	} m_NoConvStrategy;
	enum OmegaInc {
		Inc_ADDITIVE,
		Inc_MULTIPLICATIVE,
		Inc_CUSTOM
	} m_OmegaInc;
	doublereal m_dOmegaAddInc;
	doublereal m_dOmegaMulInc;
	std::vector<doublereal> m_Omega;

	integer m_iOmegaCnt;
	integer m_iPeriod;
	integer m_iPeriodCnt;
	integer m_iNumSubSteps;
	integer m_iCurrentStep;
	bool m_bStarted;
	bool m_bConverged;
	bool m_bDone;

	bool m_bPrintAllPeriods;

	bool m_bOut;
	integer m_iPeriodOut;
	doublereal m_dOmegaOut;

	enum OutputFormat {
		Out_COMPLEX,
		Out_MAGNITUDE_PHASE
	} m_OutputFormat;

	std::vector<doublereal> m_cos_psi;
	std::vector<doublereal> m_sin_psi;
	std::vector<std::vector<doublereal> > m_X;
	std::vector<doublereal> m_XPrev;
	std::vector<doublereal> m_Xcos;		// imaginary part (input is assumed to be "sin")
	std::vector<doublereal> m_Xsin;		// real part (input is assumed to be "sin")
	std::vector<doublereal> m_XcosPrev;
	std::vector<doublereal> m_XsinPrev;
	std::vector<doublereal> SUM;
	std::vector<std::vector<doublereal>> RMS;
	std::vector<doublereal> RMSOld;
	std::vector<doublereal> dAmplitude;
	std::vector<doublereal> RMSTarget;
	integer RMS_pol_order;
	integer RMS_cur_idx;
	
	struct HFInput {
		enum {
			HF_TEST = 0x1U,
			HF_OUTPUT = 0x2U
		};

		DriveCaller *m_pDC;
		unsigned m_Flag;
	};
	std::vector<HFInput> m_Input;

	doublereal compute_amplitude(
		const std::vector<std::vector<doublereal>>& RMS, 
		const std::vector<doublereal>& am,
		const doublereal& RMSt,
		const integer i_test,
		const integer RMS_cur_idx,
		const integer deg) {
		doublereal an;
		doublereal a0 = am[RMS_cur_idx];		// current
		doublereal a1 = am[(RMS_cur_idx+2)%3];	// minus 1
		doublereal a2 = am[(RMS_cur_idx+1)%3];	// minus 2
		doublereal R0 = RMS[RMS_cur_idx][i_test];	// current
		doublereal R1 = RMS[(RMS_cur_idx+2)%3][i_test];	// minus 1
		doublereal R2 = RMS[(RMS_cur_idx+1)%3][i_test];	// minus 2
		switch (deg) {
			case 0 : {
				an = RMSt / R0 * a0;
				break;
			}
			case 1 : {
				an = a0 + (RMSt - R0) * (a0 - a1) / (R0 - R1);
				if (an <= 0.) an = compute_amplitude(RMS, am, RMSt, i_test, RMS_cur_idx, 0);
				break;
			}
			case 2 : {
				doublereal R0d0 = R0 / (a0 - a1) / (a0 - a2);
				doublereal R1d1 = R1 / (a1 - a0) / (a1 - a2);
				doublereal R2d2 = R2 / (a2 - a0) / (a2 - a1);
				doublereal a = R0d0 + R1d1 + R2d2;
				doublereal b = -(R0d0 * (a1 + a2) + R1d1 * (a0 + a2) + R2d2 * (a0 + a1));
				doublereal c = (R0d0 * a1 * a2 + R1d1 * a0 * a2 + R2d2 * a0 * a1) - RMSt;
				doublereal delta = b * b - 4 * a * c;
				doublereal sqr_delta;
				if (delta >= 0) {
					sqr_delta = std::sqrt(b * b - 4 * a * c);
				} else {
					return compute_amplitude(RMS, am, RMSt, i_test, RMS_cur_idx, 1);
				}
				doublereal an1 = (-b + sqr_delta) / 2. / a;
				doublereal an2 = (-b - sqr_delta) / 2. / a;
				if (an1 > 0.)
					if (std::abs(a0 - an1) > std::abs(a0 - an2) && an2 > 0.) {
						an = an2;
					} else {
						an = an1;
				} else if (an2 > 0.) {
					an = an2;
				} else {
					an = compute_amplitude(RMS, am, RMSt, i_test, RMS_cur_idx, 1);;
				}
				
				break;
			}
			default : {
				silent_cerr("HarmonicExcitationElem::compute_amplitude unhandled deg = " << deg << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				break;
			}
		}
// 		if (a0 / an > 10) {
// 			an = a0 / 10;
// 		} else if (an / a0 > 10) {
// 			an = a0 * 10;
// 		}
		return an;
	};


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
DriveOwner(0),
m_pDM(pDM),
m_dTInit(-std::numeric_limits<doublereal>::max()),
m_dMaxDeltaT(std::numeric_limits<doublereal>::max()),
m_dPsi(0.),
m_dF(0.),
m_dAmplitude(1.),
bRMSTest(false),
bRMSTestTarget(false),
m_iMaxPeriods(0),
m_NoConvStrategy(NoConv_CONTINUE),
m_dOmegaAddInc(0.),
m_dOmegaMulInc(1.),
m_iOmegaCnt(0),
m_iPeriod(0),
m_iPeriodCnt(0),
m_iNumSubSteps(1),
m_iCurrentStep(0),
m_bStarted(false),
m_bConverged(false),
m_bDone(false),
m_bPrintAllPeriods(false),
m_bOut(false),
m_iPeriodOut(0),
m_dOmegaOut(0.),
m_OutputFormat(Out_COMPLEX),
RMS_pol_order(0),
RMS_cur_idx(0)
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
	if (bHFElem) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): warning, another HarmonicExcitationElem has already been defined, which might conflict with the one defined here at line " << HP.GetLineData() << std::endl);
		// do nothing by now
	}

	bHFElem = true;

	// inputs
	if (!HP.IsKeyWord("inputs" "number")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"inputs number\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNInput;
	try {
		iNInput = HP.GetInt(1, HighParser::range_gt<integer>(0));
	
	} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): number of input must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	m_Input.resize(iNInput);
	RMSTarget.resize(iNInput);
	bool bTest(false);
	for (integer i = 0; i < iNInput; ++i) {
		m_Input[i].m_pDC = HP.GetDriveCaller();

		m_Input[i].m_Flag = HFInput::HF_TEST | HFInput::HF_OUTPUT;
		while (true) {
			if (HP.IsKeyWord("test")) {
				bool b = HP.GetYesNoOrBool(false);
				if (b) {
					m_Input[i].m_Flag |= HFInput::HF_TEST;
				} else {
					m_Input[i].m_Flag &= ~HFInput::HF_TEST;
				}

			} else if (HP.IsKeyWord("output")) {
				bool b = HP.GetYesNoOrBool(false);
				if (b) {
					m_Input[i].m_Flag |= HFInput::HF_OUTPUT;
				} else {
					m_Input[i].m_Flag &= ~HFInput::HF_OUTPUT;
				}

			} else {
				break;
			}
		}

		if (m_Input[i].m_Flag == 0) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): warning, input #" << i << " unused, at line " << HP.GetLineData() << std::endl);
		}

		if (m_Input[i].m_Flag & HFInput::HF_TEST) {
			bTest = true;
		}
	}

	if (!bTest) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): no input is used for convergence test, at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("output" "format")) {
		if (HP.IsKeyWord("complex")) {
			m_OutputFormat = Out_COMPLEX;

		} else if (HP.IsKeyWord("magnitude" "phase")) {
			m_OutputFormat = Out_MAGNITUDE_PHASE;

		} else {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): unknown output format at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	// block size
	if (!HP.IsKeyWord("steps" "number")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"steps number\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_iN = HP.GetInt(2, HighParser::range_gt<integer>(0));

	} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): N must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	// maximum delta t
	if (HP.IsKeyWord("max" "delta" "t")) {
		try {
			m_dMaxDeltaT = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

		} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): maximum Delta t must be positive, at line " << HP.GetLineData() << std::endl);
			throw e;
		}
	}

	// frequency
	if (!HP.IsKeyWord("initial" "angular" "frequency")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"initial angular frequency\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_dOmega0 = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

	} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): initial frequency must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	if (HP.IsKeyWord("max" "angular" "frequency")) {
		if (HP.IsKeyWord("forever")) {
			m_dOmegaMax = std::numeric_limits<doublereal>::max();

		} else {
			try {
				m_dOmegaMax = HP.GetReal(m_dOmega0, HighParser::range_gt<doublereal>(m_dOmega0));

			} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): maximum frequency must be greater than initial frequency " << m_dOmega0 << ", at line " << HP.GetLineData() << std::endl);
				throw e;
			}
		}

	} else {
		m_dOmegaMax = std::numeric_limits<doublereal>::max();
	}

	m_dOmega = m_dOmega0;

	m_dDeltaT = (2*M_PI/m_dOmega)/m_iN;
	if (m_dDeltaT > m_dMaxDeltaT) {
		m_iNumSubSteps = std::ceil(m_dDeltaT / m_dMaxDeltaT);
		m_dDeltaT /= m_iNumSubSteps;
	} else {
		m_iNumSubSteps = 1;
	}
	//std::cerr << m_iNumSubSteps << " " << m_dDeltaT << std::endl;

	// frequency increment
	if (!HP.IsKeyWord("angular" "frequency" "increment")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"angular frequency increment\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("additive")) {
		m_OmegaInc = Inc_ADDITIVE;
		try {
			m_dOmegaAddInc = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

		} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): frequency additive increment must be positive, at line " << HP.GetLineData() << std::endl);
			throw e;
		}

	} else if (HP.IsKeyWord("multiplicative")) {
		m_OmegaInc = Inc_MULTIPLICATIVE;
		try {
			m_dOmegaMulInc = HP.GetReal(1., HighParser::range_gt<doublereal>(1.));

		} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): frequency multiplicative increment must be greater than 1, at line " << HP.GetLineData() << std::endl);
			throw e;
		}

	} else if (HP.IsKeyWord("custom")) {
		m_OmegaInc = Inc_CUSTOM;

		doublereal dOmega = m_dOmega0;
		m_Omega.push_back(dOmega);

		for (integer i = 1; i++; ) {
			if (HP.IsKeyWord("last")) {
				break;
			}

			try {
				dOmega = HP.GetReal(1., HighParser::range_gt<doublereal>(dOmega));

			} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): Omega[" << i << "] must be greater than Omega[" << i - 1 << "] = " << m_Omega[i-1] << ", at line " << HP.GetLineData() << std::endl);
				throw e;
			}

			m_Omega.push_back(dOmega);
		}

	} else {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): unknown or missing frequency increment method at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// initial time
	if (HP.IsKeyWord("initial" "time")) {
		m_dTInit = HP.GetReal();
		if (HP.IsKeyWord("prescribed" "time" "step")) {
			DriveCaller *pDC = HP.GetDriveCaller();
			DriveOwner::Set(pDC);
			if (pDC->dGet(m_dTInit) != m_dDeltaT) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): Prescribed initial time step "
					<< pDC->dGet(m_dTInit) <<
					"\nevaluated at the initial time t = " << m_dTInit <<
					"is different from initial time step " << m_dDeltaT << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	} else {
		m_dTInit = m_pDM->dGetTime();
	}

	// tolerance
	if (HP.IsKeyWord("RMS" "test")) {
		bRMSTest = true;
		if (HP.IsKeyWord("RMS" "target"))  {
			integer count_tests = 0;
			for (integer i = 0; i < iNInput; ++i) {
				if (m_Input[i].m_Flag & HFInput::HF_TEST) {
					count_tests++;
				}
			}
			if (count_tests != 1)  {
				silent_cerr("\nHarmonicExcitationElem(" << GetLabel() << "): \"RMS target\" required but"
					"\nthere are " << count_tests <<
					" test variables. There should be a single test variable" 
					"\nfor this options to make sense.\n"<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bRMSTestTarget = true;
			for (integer i = 0; i < iNInput; ++i) {
				if (m_Input[i].m_Flag & HFInput::HF_TEST) {
					try {
						RMSTarget[i] = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));
					} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
						silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): "
						"RMS target value must be positive, at line " << HP.GetLineData() << std::endl);
						throw e;
					}
				}
			}
		}
	}
	if (!HP.IsKeyWord("tolerance")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"tolerance\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_dTol = HP.GetReal(0., HighParser::range_gt<doublereal>(0.));

	} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): tolerance must be positive, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	// minimum number of blocks
	if (!HP.IsKeyWord("min" "periods")) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): \"min periods\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	try {
		m_iMinPeriods = HP.GetInt(2, HighParser::range_gt<integer>(1));

	} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
		silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): minimum number of periods must be greater than 1, at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	// maximum number of blocks
	if (HP.IsKeyWord("max" "periods")) {
		try {
			m_iMaxPeriods = HP.GetInt(m_iMinPeriods + 1, HighParser::range_gt<integer>(m_iMinPeriods));

		} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): maximum number of periods must be greater than minimum number of periods " << m_iMinPeriods << ", at line " << HP.GetLineData() << std::endl);
			throw e;
		}

		m_NoConvStrategy = NoConv_CONTINUE;
		if (HP.IsKeyWord("no" "convergence" "strategy")) {
			if (HP.IsKeyWord("continue")) {
				m_NoConvStrategy = NoConv_CONTINUE;

			} else if (HP.IsKeyWord("abort")) {
				m_NoConvStrategy = NoConv_ABORT;

			} else {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): "
					"unknown no convergence strategy at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	// instantiate timestep drive caller
	if (HP.IsKeyWord("timestep" "drive" "label")) {
		unsigned uTSLabel;
		try {
			uTSLabel = HP.GetInt(0, HighParser::range_ge<integer>(0));

		} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): timestep label must be non-negative, at line " << HP.GetLineData() << std::endl);
			throw e;
		}

		// check if driver already exists
		if (HP.GetDrive(uTSLabel) != 0) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): "
				"drive caller " << uLabel << " already defined "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		DriveCaller *pDirectDC = 0;
		SAFENEWWITHCONSTRUCTOR(pDirectDC,
			DirectDriveCaller,
			DirectDriveCaller(pDM->pGetDrvHdl()));

		DriveCaller *pDC = 0;
		SAFENEWWITHCONSTRUCTOR(pDC,
			PrivDriveCaller,
			PrivDriveCaller(pDM->pGetDrvHdl(),	// pDrvHdl,
				pDirectDC, 			// pTmp; direct
				this,				// pSE
				Priv_DT,			// iIndex,
				"timestep"));			// sIndexName.c_str()

		try {
			(void)HP.SetDrive(uTSLabel, pDC);
		} catch (const MBDynParser::ErrGeneric& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): unable to create timestep driver: " << e.what() << " at line " << HP.GetLineData() << std::endl);
			throw e;
		}
	}

	if (HP.IsKeyWord("print" "all" "periods"))
		m_bPrintAllPeriods = true;

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	m_cos_psi.resize(m_iN);
	m_sin_psi.resize(m_iN);
	m_X.resize(m_iN);

	m_XPrev.resize(iNInput);
	m_Xcos.resize(iNInput);
	m_Xsin.resize(iNInput);
	m_XcosPrev.resize(iNInput);
	m_XsinPrev.resize(iNInput);
	SUM.resize(iNInput);
	RMS.resize(3);
	dAmplitude.resize(3);
	for (integer i = 0; i < 3; i++) {
		RMS[i].resize(iNInput);
		std::fill(RMS[i].begin(), RMS[i].end(), 0.);
	}
	RMSOld.resize(iNInput);
	std::fill(RMSOld.begin(), RMSOld.end(), 0.);

	for (integer i = 0; i < m_iN; i++) {
		doublereal dPsi = (2*M_PI*i)/m_iN;
		m_cos_psi[i] = cos(dPsi)/(m_iN/2);
		m_sin_psi[i] = sin(dPsi)/(m_iN/2);

		m_X[i].resize(iNInput);
	}

	// reset
	std::fill(dAmplitude.begin(), dAmplitude.end(), m_dAmplitude);
	std::fill(m_XPrev.begin(), m_XPrev.end(), 0.);
	std::fill(m_Xcos.begin(), m_Xcos.end(), 0.);
	std::fill(m_Xsin.begin(), m_Xsin.end(), 0.);
	std::fill(SUM.begin(), SUM.end(), 0.);
}

HarmonicForcingElem::~HarmonicForcingElem(void)
{
	// destroy private data
	for (std::vector<HFInput>::iterator i = m_Input.begin(); i != m_Input.end(); ++i) {
		delete i->m_pDC;
	}
}

void
HarmonicForcingElem::Output(OutputHandler& OH) const
{
	// NOTE: output occurs only at convergence,
	// not with fixed periodicity
	if (bToBeOutput()) {
		if (((m_bOut && m_iPeriod == 0) || (m_bPrintAllPeriods && m_pDM->dGetTime()>m_dTInit))
			&& m_iPeriodCnt == 0 && m_iCurrentStep == 0) {
			std::ostream& out = OH.Loadable();

			out << std::setw(8) << GetLabel()	// 1:	label
				<< " " << m_pDM->dGetTime();	// 2:	when
			if(m_iPeriod==0) // just converged
				out << " " << m_dOmegaOut		// 3:	what frequency
					<< " " << m_iPeriodOut;		// 4:	how many periods required (-1 if not converged?)
			else
				out << " " << m_dOmega		// 3:	what frequency
					<< " " << m_iPeriod;		// 4:	how many periods required (-1 if not converged?)

			switch (m_OutputFormat) {
			case Out_COMPLEX:
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_OUTPUT) {
						// real imag
						out << " " << m_Xsin[i] / m_dAmplitude << " " << m_Xcos[i] / m_dAmplitude;
					}
				}
				break;

			case Out_MAGNITUDE_PHASE:
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_OUTPUT) {
						// magnitude phase
						out << " " << std::sqrt(m_Xsin[i]*m_Xsin[i] + m_Xcos[i]*m_Xcos[i])  / m_dAmplitude << " " << std::atan2(m_Xcos[i], m_Xsin[i]);
					}
				}
				break;
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
		if (m_bDone) {
			// throw it here, otherwise output might not be performed
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}

		// if just initialized
		if (m_iPeriod == 0 && m_iPeriodCnt == 0) {
			m_dT0 = m_pDM->dGetTime();

			// reset vectors
			std::fill(m_Xcos.begin(), m_Xcos.end(), 0.);
			std::fill(m_Xsin.begin(), m_Xsin.end(), 0.);
		}
	}

	if (m_bStarted) {
		m_dPsi = m_dOmega*(m_pDM->dGetTime() - m_dT0);
		m_dF = m_dAmplitude * sin(m_dPsi);
	}
}

void
HarmonicForcingElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	//std::cerr << "X" << std::endl;
	if (!m_bStarted) {
		if (DriveOwner::pGetDriveCaller())
			m_dDeltaT = DriveOwner::dGet(m_pDM->dGetTime());
		return;
	}

	++m_iCurrentStep;
	if (m_iCurrentStep != m_iNumSubSteps) {
		return;
	}
	
	m_iCurrentStep = 0;

	// collect input
	doublereal dErr = 0.;	

	for (unsigned i = 0; i < m_Input.size(); ++i) {
		if (m_iPeriod > 0) {
			m_XcosPrev[i] = m_Xcos[i];
			m_XsinPrev[i] = m_Xsin[i];

			m_Xcos[i] -= m_X[m_iPeriodCnt][i]*m_cos_psi[m_iPeriodCnt];
			m_Xsin[i] -= m_X[m_iPeriodCnt][i]*m_sin_psi[m_iPeriodCnt];
		}

		m_X[m_iPeriodCnt][i] = m_Input[i].m_pDC->dGet();

		m_Xcos[i] += m_X[m_iPeriodCnt][i]*m_cos_psi[m_iPeriodCnt];
		m_Xsin[i] += m_X[m_iPeriodCnt][i]*m_sin_psi[m_iPeriodCnt];

		if (m_iPeriod > 0 && (m_Input[i].m_Flag & HFInput::HF_TEST)) {
			doublereal d = m_Xcos[i] - m_XcosPrev[i];
			dErr += d*d;
			d = m_Xsin[i] - m_XsinPrev[i];
			dErr += d*d;
			RMS[RMS_cur_idx][i] += (m_X[m_iPeriodCnt][i] * m_X[m_iPeriodCnt][i] + m_XPrev[i] *  m_XPrev[i]) / 2.;
			SUM[i] += (m_X[m_iPeriodCnt][i] + m_XPrev[i]) / 2.;
		}
		m_XPrev[i] = m_X[m_iPeriodCnt][i];
	}

	// check for convergence
	++m_iPeriodCnt;
	bool bRMSConverged = false;
	if (m_iPeriod > 0) {
		if (bRMSTest) {
			if (m_iPeriodCnt == m_iN && bRMSTest) {
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					RMS[RMS_cur_idx][i] = std::sqrt(RMS[RMS_cur_idx][i] / m_iN - SUM[i] * SUM[i] / m_iN / m_iN);
				}
				m_bConverged = true;
				bRMSConverged = true;
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_TEST) {
						bRMSConverged = bRMSConverged && (std::abs(RMS[RMS_cur_idx][i]-RMSOld[i]) /
							std::max(RMS[RMS_cur_idx][i], RMSOld[i]) <= m_dTol);
					}
				}
				RMSOld = RMS[RMS_cur_idx];
				m_bConverged = bRMSConverged;
				if (bRMSConverged && bRMSTestTarget) {
					for (unsigned i = 0; i < m_Input.size(); ++i) {
						if (m_Input[i].m_Flag & HFInput::HF_TEST) {
							m_bConverged = m_bConverged && (std::abs(RMS[RMS_cur_idx][i]-RMSTarget[i]) /
								RMSTarget[i] <= m_dTol);
						}
					}
					if (!m_bConverged) {
						for (unsigned i = 0; i < m_Input.size(); ++i) {
							if (m_Input[i].m_Flag & HFInput::HF_TEST) {
								m_dAmplitude = compute_amplitude(RMS, dAmplitude, RMSTarget[i], i, RMS_cur_idx, RMS_pol_order);

// 								std::cerr << 
// 									"Input " << i << " " << RMS_cur_idx << " " << RMS_pol_order << " "
// 									" " << RMSTarget[i] << 
// 									" " << RMS[RMS_cur_idx][i] <<  
// 									" " << RMS[(RMS_cur_idx+2)%3][i] << 
// 									" " << RMS[(RMS_cur_idx+1)%3][i] << 
// 									" " << m_dAmplitude << 
// 									" " << dAmplitude[(RMS_cur_idx)%3] << 
// 									" " << dAmplitude[(RMS_cur_idx+2)%3] << 
// 									" " << dAmplitude[(RMS_cur_idx+1)%3] << 
// 									std::endl;
							}
						}
						dAmplitude[(RMS_cur_idx+1)%3] = m_dAmplitude;
						RMS_pol_order = std::min(RMS_pol_order+1, 2);
					} else {
						for (unsigned i = 0; i < m_Input.size(); ++i) {
							if (m_Input[i].m_Flag & HFInput::HF_TEST) {
// 								std::cerr << 
// 									"Ouput " << i << " " << RMS_cur_idx << " " << RMS_pol_order << " "
// 									" " << RMSTarget[i] << 
// 									" " << RMS[RMS_cur_idx][i] <<  
// 									" " << RMS[(RMS_cur_idx+2)%3][i] << 
// 									" " << RMS[(RMS_cur_idx+1)%3][i] << 
// 									" " << m_dAmplitude << 
// 									" " << dAmplitude[(RMS_cur_idx)%3] << 
// 									" " << dAmplitude[(RMS_cur_idx+2)%3] << 
// 									" " << dAmplitude[(RMS_cur_idx+1)%3] << 
// 									std::endl;
							}
						}
						std::fill(dAmplitude.begin(), dAmplitude.end(), m_dAmplitude);
						RMS_pol_order = 0;
					}
				}
			}
		} else {
			dErr = sqrt(dErr);
			if (dErr <= m_dTol) {
				m_bConverged = true;
			} else {
				m_bConverged = false;
			}
		}
	}

	if (m_iPeriodCnt == m_iN) {
		m_iPeriodCnt = 0;
		if (bRMSConverged) {
			RMS_cur_idx = (RMS_cur_idx+1)%3;//RMSOld = RMS;
		}
		std::fill(RMS[RMS_cur_idx].begin(), RMS[RMS_cur_idx].end(), 0.);
		std::fill(SUM.begin(), SUM.end(), 0.);
		++m_iPeriod;
		if (m_iPeriod >= m_iMinPeriods) {
			if (m_bConverged || ((m_iMaxPeriods > 0) && (m_iPeriod >= m_iMaxPeriods) && (m_NoConvStrategy == NoConv_CONTINUE))) {
				m_bOut = true;
				if (m_bConverged) {
					m_iPeriodOut = m_iPeriod;
				} else {
					m_iPeriodOut = -m_iPeriod;
				}
				m_dOmegaOut = m_dOmega;
				m_iOmegaCnt++;

				m_bConverged = false;
				m_iPeriod = 0;

				switch (m_OmegaInc) {
				case Inc_ADDITIVE:
					m_dOmega += m_dOmegaAddInc;
					break;

				case Inc_MULTIPLICATIVE:
					m_dOmega *= m_dOmegaMulInc;
					break;

				case Inc_CUSTOM:
					if ((unsigned)m_iOmegaCnt >= m_Omega.size()) {
						// schedule for termination!
						m_bDone = true;
					}
					m_dOmega = m_Omega[m_iOmegaCnt];
					break;

				default:
					// impossible
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (m_dOmega > m_dOmegaMax) {
					// schedule for termination!
					m_bDone = true;
				}

				m_dDeltaT = (2*M_PI/m_dOmega)/m_iN;
				if (m_dDeltaT > m_dMaxDeltaT) {
					m_iNumSubSteps = std::ceil(m_dDeltaT / m_dMaxDeltaT);
					m_dDeltaT /= m_iNumSubSteps;
				} else {
					m_iNumSubSteps = 1;
				}
				std::cerr << m_iNumSubSteps << " " << m_dDeltaT << std::endl;

			} else if ((m_iMaxPeriods > 0) && (m_iPeriod > m_iMaxPeriods) && (m_NoConvStrategy == NoConv_ABORT)) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): max periods " << m_iMaxPeriods << " exceeded for frequency " << m_dOmega << "; aborting..." << std::endl);

				// schedule for termination!
				m_bDone = true;
			}
		}
	}
}

unsigned int
HarmonicForcingElem::iGetNumPrivData(void) const
{
	return (Priv_LAST - 1);
}

unsigned int
HarmonicForcingElem::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "timestep") == 0) {
		return Priv_DT;

	} else if (strcmp(s, "excitation") == 0) {
		return Priv_F;

	} else if (strcmp(s, "psi") == 0) {
		return Priv_PSI;

	} else if (strcmp(s, "omega") == 0) {
		return Priv_OMEGA;

	} else if (strcmp(s, "amplitude") == 0) {
		return Priv_AMPLITUDE;

	} else if (strcmp(s, "count") == 0) {
		return Priv_COUNT;
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
HarmonicForcingElem::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case Priv_DT:
		return m_dDeltaT;

	case Priv_F:
		return m_dF;

	case Priv_PSI:
		return m_dPsi;

	case Priv_OMEGA:
		return m_dOmega;

	case Priv_AMPLITUDE:
		return m_dAmplitude;

	case Priv_COUNT:
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

#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!hfelem_set()) {
		silent_cerr("HFElem: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}
	return 0;
}
#endif // ! STATIC_MODULES
