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

#define HFELEM_COMMENTS 0

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
		Priv_OUTPUT,
		Priv_LAST
	};

	const DataManager* m_pDM;

	integer m_iN;

	doublereal m_dTInit;
    doublereal m_dTimestepsCompareTol; // hard coded variable to check if time steps from drive callers are different
	doublereal m_dT0;
	doublereal m_dDeltaT;
	doublereal m_dMaxDeltaT;
	DriveCaller *m_pDC_MaxDeltaT;
	doublereal m_dOmega;
	doublereal m_dOmega0;
	doublereal m_dOmegaMax;
	doublereal m_dPsi;
	doublereal m_dF;

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
	integer m_iPeriodRMS;
	integer m_iPeriodCnt;
	integer m_iNumSubSteps;
	integer m_iCurrentStep;
	bool m_bStarted;
	bool m_bConverged;
	bool m_bDone;

	bool m_bPrintAllPeriods;

	integer m_iPeriodOut;
	doublereal m_dOmegaOut;

	enum OutputFormat {
		Out_COMPLEX,
		Out_MAGNITUDE_PHASE
	} m_OutputFormat;
	bool b_OutputNormalized;

	std::vector<doublereal> m_cos_psi;
	std::vector<doublereal> m_sin_psi;
	std::vector<std::vector<doublereal> > m_X;
	std::vector<doublereal> m_XPrev;
	std::vector<doublereal> m_Xcos;		// imaginary part (input is assumed to be "sin")
	std::vector<doublereal> m_Xsin;		// real part (input is assumed to be "sin")
	std::vector<doublereal> m_XcosPrev;
	std::vector<doublereal> m_XsinPrev;
	std::vector<doublereal> SUM;
	std::vector<doublereal> RMSNow;
	std::vector<doublereal> RMSUpp;
	std::vector<doublereal> RMSLow;
	std::vector<std::vector<doublereal>> RMSPrev;
	std::vector<doublereal> RMSBest;
	std::vector<DriveCaller*> m_pDC_RMSTarget;
	std::vector<doublereal> RMSTarget;
	doublereal m_dAmplitude;
	doublereal m_dAmplitudeOut;
	doublereal m_dAmplitudeUpp;
	doublereal m_dAmplitudeLow;
	doublereal m_dAmplitudeBest;
	doublereal m_dRMSovershoot;
	integer m_iAmpliConvergences;
	integer m_iMinAmpliConvergences;
	integer m_iRMSrestarts;
	integer m_iMaxRMSrestarts;
	integer m_iAfterMaxRMSrestartsPeriods;
	integer m_iRMSprevs;
	bool bRMSUpp;
	bool bRMSLow;

	integer m_bShouldBeOutput;
	integer m_iWritePeriodsAfterConvergence;
	integer m_iWriteAfterConvergenceEvery;
	integer m_iWriteAfterConvergenceCount;
	integer m_iConvergedPeriod;

	struct HFInput {
		enum {
			HF_TEST = 0x1U,
			HF_OUTPUT = 0x2U,
			HF_TARGET = 0x4U
		};

		DriveCaller *m_pDC;
		unsigned m_Flag;
	};
	std::vector<HFInput> m_Input;

	doublereal compute_amplitude(
		const doublereal& R0,
		const doublereal& R1,
		const doublereal& a0,
		const doublereal& a1,
		const doublereal& RMSt,
		const integer operational) {
		doublereal an;
		switch (operational) {
			case 0 : {
				an = RMSt / R0 * a0;
				break;
			}
			case 1 : {
				an = a0 + (RMSt - R0) * (a0 - a1) / (R0 - R1);
				break;
			}
			default : {
				silent_cerr("HarmonicExcitationElem::compute_amplitude unhandled operational = " << operational << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				break;
			}
		}
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
m_dTimestepsCompareTol(1.0e-5),
m_pDC_MaxDeltaT(0),
m_dPsi(0.),
m_dF(0.),
bRMSTest(false),
bRMSTestTarget(false),
m_iMaxPeriods(0),
m_NoConvStrategy(NoConv_CONTINUE),
m_dOmegaAddInc(0.),
m_dOmegaMulInc(1.),
m_iOmegaCnt(0),
m_iPeriod(0),
m_iPeriodRMS(0),
m_iPeriodCnt(0),
m_iNumSubSteps(1),
m_iCurrentStep(0),
m_bStarted(false),
m_bConverged(false),
m_bDone(false),
m_bPrintAllPeriods(false),
m_iPeriodOut(0),
m_dOmegaOut(0.),
m_OutputFormat(Out_COMPLEX),
b_OutputNormalized(false),
m_dAmplitude(1.),
m_dAmplitudeOut(1.),
m_dAmplitudeBest(1.),
m_dRMSovershoot(0.),
m_iAmpliConvergences(0),
m_iMinAmpliConvergences(1),
m_iRMSrestarts(0),
m_iMaxRMSrestarts(5),
bRMSUpp(false),
bRMSLow(false),
m_bShouldBeOutput(false),
m_iWritePeriodsAfterConvergence(0),
m_iWriteAfterConvergenceEvery(1),
m_iWriteAfterConvergenceCount(0),
m_iConvergedPeriod(0)
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
	m_pDC_RMSTarget.resize(iNInput);
	bool bTest(false);
	for (integer i = 0; i < iNInput; ++i) {
		m_Input[i].m_pDC = HP.GetDriveCaller();

		m_Input[i].m_Flag = HFInput::HF_TEST | HFInput::HF_OUTPUT;
		while (true) {
            if (HP.IsKeyWord("test")) { // variable to be tested for convergence - default: true
				bool b = HP.GetYesNoOrBool(false);
				if (b) {
					m_Input[i].m_Flag |= HFInput::HF_TEST;
				} else {
					m_Input[i].m_Flag &= ~HFInput::HF_TEST;
				}

            } else if (HP.IsKeyWord("output")) { // variable to be written in output - default: true
				bool b = HP.GetYesNoOrBool(false);
				if (b) {
					m_Input[i].m_Flag |= HFInput::HF_OUTPUT;
				} else {
					m_Input[i].m_Flag &= ~HFInput::HF_OUTPUT;
				}

            } else if (HP.IsKeyWord("target")) { // variable to be controlled - default: false
				bool b = HP.GetYesNoOrBool(false);
				if (b) {
					m_Input[i].m_Flag |= HFInput::HF_TARGET;
				} else {
					m_Input[i].m_Flag &= ~HFInput::HF_TARGET;
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
		if (HP.IsKeyWord("normalized")) {
			b_OutputNormalized = true;
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

    // maximum delta t - drive caller function of omega
	if (HP.IsKeyWord("max" "delta" "t")) {
		m_pDC_MaxDeltaT = HP.GetDriveCaller();

		if (m_pDC_MaxDeltaT->dGet(1.0)<=0) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): max delta t drive caller for test omega=1 "
				<< m_pDC_MaxDeltaT->dGet(1.0) <<
				"is not positive" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	if(m_pDC_MaxDeltaT)
		m_dMaxDeltaT = m_pDC_MaxDeltaT->dGet(m_dOmega);
	else
		m_dMaxDeltaT = m_dDeltaT;
	if (m_dDeltaT > m_dMaxDeltaT) {
		m_iNumSubSteps = std::ceil(m_dDeltaT / m_dMaxDeltaT);
        if(((m_dDeltaT/m_dMaxDeltaT)/std::floor(m_dDeltaT/m_dMaxDeltaT)-1.0)<m_dTimestepsCompareTol && ((m_dDeltaT/m_dMaxDeltaT)/std::floor(m_dDeltaT/m_dMaxDeltaT)-1.0)>0.0)
			m_iNumSubSteps -= 1;
		m_dDeltaT /= m_iNumSubSteps;
	} else {
		m_iNumSubSteps = 1;
	}
#if HFELEM_COMMENTS == 1
	std::cerr << "omega: " << m_dOmega << " | DeltaT0: " << (2*M_PI/m_dOmega)/m_iN << " | MaxDeltaT: " << m_dMaxDeltaT << " | substeps: " << m_iNumSubSteps << " | DeltaT: " << m_dDeltaT << std::endl;
	std::cerr << "m_dDeltaT / m_dMaxDeltaT: " << m_dDeltaT / m_dMaxDeltaT << " | err: " << (m_dDeltaT/m_dMaxDeltaT)/std::floor(m_dDeltaT/m_dMaxDeltaT)-1.0 << std::endl;
#endif

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
			//if (pDC->dGet(m_dTInit) != m_dDeltaT) {
			//Biondani
            if (std::abs(pDC->dGet(m_dTInit) - m_dDeltaT)/m_dDeltaT > m_dTimestepsCompareTol ) {
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
		integer count_targets = 0;
		for (integer i = 0; i < iNInput; ++i) {
			if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
				count_targets++;
			}
		}
		if (HP.IsKeyWord("RMS" "target"))  {
			if (count_targets > 1)  {
				silent_cerr("\nHarmonicExcitationElem(" << GetLabel() << "): \nthere are " << count_targets << " target variables." 
					"\nThere should be at most one target variable for this options to make sense.\n"<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			} else if (count_targets < 1){
				silent_cerr("\nHarmonicExcitationElem(" << GetLabel() << "): \nthere are " << count_targets << " target variables." 
					"\nThere should be at least one target variable for this options to make sense.\n"<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bRMSTestTarget = true;
			for (integer i = 0; i < iNInput; ++i) {
				if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
					m_pDC_RMSTarget[i] = HP.GetDriveCaller();
					RMSTarget[i]=m_pDC_RMSTarget[i]->dGet(m_dOmega);
					if (RMSTarget[i]<=0) {
						silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): RMS target drive caller for test omega=omega0 "
							<< RMSTarget[i] <<
							"is not positive" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}
			if (HP.IsKeyWord("initial" "amplitude"))  {
				try {
					m_dAmplitude = HP.GetReal(0., HighParser::range_ge<doublereal>(0.0));
					m_dAmplitudeOut = m_dAmplitude;
				} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
					silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): RMS target, initial amplitude must be greater than 0, at line " << HP.GetLineData() << std::endl);
					throw e;
				}
			}
			if (HP.IsKeyWord("overshoot"))  {
				try {
					m_dRMSovershoot = HP.GetReal(0., HighParser::range_ge<doublereal>(0.0));
				} catch (const HighParser::ErrValueOutOfRange<doublereal>& e) {
					silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): RMS target, tolerance overshoot must be greater than 0, at line " << HP.GetLineData() << std::endl);
					throw e;
				}
			}
		}else{
			if(count_targets>0){
				silent_cerr("\nHarmonicExcitationElem(" << GetLabel() << "): \nthere are " << count_targets << " target variables," 
					"\nbut no target value specified.\n"<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		if (HP.IsKeyWord("RMS" "prev" "periods"))  {
			try {
				m_iRMSprevs = HP.GetInt(1, HighParser::range_ge<integer>(1));
			} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): RMS prev periods must be greater than 0, at line " << HP.GetLineData() << std::endl);
				throw e;
			}
		}else{
			m_iRMSprevs = 2;
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

	if (HP.IsKeyWord("write" "after" "convergence" "periods")){
		try {
			m_iWritePeriodsAfterConvergence = HP.GetInt(0, HighParser::range_ge<integer>(0));

		} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): periods to be written after convergence must be non-negative, at line " << HP.GetLineData() << std::endl);
			throw e;
		}
		if (HP.IsKeyWord("frequency")){
			try {
				m_iWriteAfterConvergenceEvery = HP.GetInt(0, HighParser::range_gt<integer>(0));

			} catch (const HighParser::ErrValueOutOfRange<integer>& e) {
				silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): writing frequency after convergence must be positive, at line " << HP.GetLineData() << std::endl);
				throw e;
			}
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	m_cos_psi.resize(m_iN);
	m_sin_psi.resize(m_iN);
	m_X.resize(m_iN);

	m_XPrev.resize(iNInput);
	m_Xcos.resize(iNInput);
	m_Xsin.resize(iNInput);
	m_XcosPrev.resize(iNInput);
	m_XsinPrev.resize(iNInput);
	if(bRMSTest){
		SUM.resize(iNInput);
		RMSNow.resize(iNInput);
		std::fill(RMSNow.begin(), RMSNow.end(), 0.);
		RMSPrev.resize(std::max(m_iRMSprevs,3));
		for (unsigned i = 0; i < RMSPrev.size(); i++) {
			RMSPrev[i].resize(iNInput);
			std::fill(RMSPrev[i].begin(), RMSPrev[i].end(), -1.);
		}
		RMSBest.resize(iNInput);
		std::fill(RMSBest.begin(), RMSBest.end(), 0.);
		std::fill(SUM.begin(), SUM.end(), 0.);
		m_iAfterMaxRMSrestartsPeriods=(m_iMaxPeriods-1)/10 + 1;
    }

	for (integer i = 0; i < m_iN; i++) {
		doublereal dPsi = (2*M_PI*i)/m_iN;
		m_cos_psi[i] = cos(dPsi)/(m_iN/2);
		m_sin_psi[i] = sin(dPsi)/(m_iN/2);

		m_X[i].resize(iNInput);
	}

	// reset
	std::fill(m_XPrev.begin(), m_XPrev.end(), 0.);
	std::fill(m_Xcos.begin(), m_Xcos.end(), 0.);
	std::fill(m_Xsin.begin(), m_Xsin.end(), 0.);
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
		if (((m_iPeriodOut!=0 && m_iPeriod == 0) || (m_bPrintAllPeriods && m_pDM->dGetTime()>m_dTInit))
			&& m_iPeriodCnt == 0 && m_iCurrentStep == 0) {
			std::ostream& out = OH.Loadable();

			out << std::setprecision(16) << GetLabel()	// 1:	label
				<< " " << m_pDM->dGetTime()	// 2:	when
                << std::setprecision(8);
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
						if (b_OutputNormalized) {
							out << " " << m_Xsin[i] / m_dAmplitude << " " << m_Xcos[i] / m_dAmplitude;
						} else {
							out << " " << m_Xsin[i] << " " << m_Xcos[i];
						}
					}
				}
				break;

			case Out_MAGNITUDE_PHASE:
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_OUTPUT) {
						// magnitude phase
						if (b_OutputNormalized) {
							out << " " << std::sqrt(m_Xsin[i]*m_Xsin[i] + m_Xcos[i]*m_Xcos[i]) / m_dAmplitude << " " << std::atan2(m_Xcos[i], m_Xsin[i]);
						} else {
							out << " " << std::sqrt(m_Xsin[i]*m_Xsin[i] + m_Xcos[i]*m_Xcos[i]) << " " << std::atan2(m_Xcos[i], m_Xsin[i]);
						}
					}
				}
				break;
			}
			if(bRMSTestTarget){
				out << " " << m_dAmplitudeOut;
				for (unsigned i = 0; i < m_Input.size(); ++i)
					if (m_Input[i].m_Flag & HFInput::HF_TARGET)
						out << " " << RMSPrev[0][i];
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
#if HFELEM_COMMENTS == 1
			std::cerr << "STARTED" << std::endl;
			std::cerr << "omega: " << m_dOmega << " | DeltaT0: " << (2*M_PI/m_dOmega)/m_iN << " | MaxDeltaT: " << m_dMaxDeltaT << " | substeps: " << m_iNumSubSteps << " | DeltaT: " << m_dDeltaT << std::endl;
#endif
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
	else{
        m_dPsi = 0;
        m_dF = 0;
    }
}

void
HarmonicForcingElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	if (!m_bStarted) {
		if(m_pDM->dGetTime()==0)
			m_bShouldBeOutput=true;
		else
			m_bShouldBeOutput=false;
		if (DriveOwner::pGetDriveCaller())
			m_dDeltaT = DriveOwner::dGet(m_pDM->dGetTime());
		return;
	}

	if(m_bConverged && m_iConvergedPeriod>0){
		m_bShouldBeOutput=true;
	}else{
		m_bShouldBeOutput=false;
	}

	if(m_bShouldBeOutput){
		if (m_iWriteAfterConvergenceCount%m_iWriteAfterConvergenceEvery>0)
			m_bShouldBeOutput=false;
		m_iWriteAfterConvergenceCount++;
	}else
		m_iWriteAfterConvergenceCount=0;

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
			if (bRMSTest) {
				RMSNow[i] += (m_X[m_iPeriodCnt][i] * m_X[m_iPeriodCnt][i] + m_XPrev[i] *  m_XPrev[i]) / 2.;
				SUM[i] += (m_X[m_iPeriodCnt][i] + m_XPrev[i]) / 2.;
			}
		}
		m_XPrev[i] = m_X[m_iPeriodCnt][i];
	}
	// check for convergence
	++m_iPeriodCnt;
	bool bRMSConverged = false;
	bool bRMSConvergedDouble = true;
	bool bRMSTooLow = false;
	bool bRMSTooHigh = false;
	if (m_iPeriod > 0) {
		if (bRMSTest) {
			if (m_iPeriodCnt == m_iN) { // the period is over
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					RMSNow[i] = std::sqrt(RMSNow[i] / m_iN - SUM[i] * SUM[i] / m_iN / m_iN);
				}

				// check if RMS converged
				m_iPeriodRMS++;
				bRMSConverged = m_iPeriodRMS >= m_iMinPeriods && m_iPeriod > m_iConvergedPeriod+m_iWritePeriodsAfterConvergence;
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_TEST) {
						bool bRMSConvergedSingle = true; // true when there are no alternating amplitudes
						for (int j = 0; j < m_iRMSprevs; ++j){
							bRMSConvergedSingle = bRMSConvergedSingle && (std::abs(RMSNow[i]-RMSPrev[j][i]) /
								std::max(RMSNow[i], RMSPrev[j][i]) <= m_dTol);
						}
						if(!bRMSConvergedSingle && (std::abs(RMSNow[i]-RMSPrev[0][i]) /
								std::max(RMSNow[i], RMSPrev[0][i]))>= 4*m_dTol  ){ // check for double alternating amplitude
							for (unsigned j = 0; j < (RMSPrev.size()-1)/2; ++j){
								bRMSConvergedDouble = bRMSConvergedDouble && (std::abs(RMSNow[i]-RMSPrev[2*j+1][i]) /
									std::max(RMSNow[i], RMSPrev[2*j+1][i]) <= m_dTol) &&
								(std::abs(RMSPrev[0][i]-RMSPrev[2*j+2][i]) / std::max(RMSPrev[0][i],RMSPrev[2*j+2][i]) <= m_dTol);
							}
							bRMSConverged = bRMSConverged && bRMSConvergedDouble;
						}else{
							bRMSConvergedDouble = false;
							bRMSConverged = bRMSConverged && bRMSConvergedSingle;
						}
					}
				}

				// check if response would be anyway too low or too high
				if(!bRMSConverged && bRMSTestTarget && m_iPeriodRMS >= m_iMinPeriods && m_iRMSrestarts==0){
					for (unsigned i = 0; i < m_Input.size(); ++i) {
						if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
							bRMSTooLow=true;
							bRMSTooHigh=true;
							doublereal RMSmin = RMSNow[i];
							integer iRMSmin = -1;
							doublereal RMSmax = RMSNow[i];
							integer iRMSmax = -1;
							for (unsigned j = 0; j < RMSPrev.size(); ++j){
								bRMSTooLow  = bRMSTooLow  && ((1+m_dTol)*RMSPrev[j][i]<RMSTarget[i]);
								bRMSTooHigh = bRMSTooHigh && ((1-m_dTol)*RMSPrev[j][i]>RMSTarget[i]);
								if(RMSPrev[j][i]<RMSmin){
									RMSmin=RMSPrev[j][i];
									iRMSmin=j;
								}
								if(RMSPrev[j][i]>RMSmax){
									RMSmax=RMSPrev[j][i];
									iRMSmax=j;
								}
							}
							if((iRMSmin==-1 && iRMSmax==int(RMSPrev.size())-1)||(iRMSmin==int(RMSPrev.size())-1 && iRMSmax==-1)){
								bRMSTooLow  = false;
								bRMSTooHigh = false;
							}
							bRMSTooLow  = bRMSTooLow  && ((1+m_dTol)*RMSNow[i]<RMSTarget[i]);
							bRMSTooHigh = bRMSTooHigh && ((1-m_dTol)*RMSNow[i]>RMSTarget[i]);
							if(bRMSTooLow){
								RMSNow[i]=RMSmin;
								bRMSConverged=true;
							}
							else if(bRMSTooHigh){
								RMSNow[i]=RMSmax;
								bRMSConverged=true;
							}
						}
					}
				}

				if(m_iPeriod > m_iConvergedPeriod+m_iWritePeriodsAfterConvergence)
					m_bConverged = bRMSConverged;
				if (bRMSConverged) {
					if (bRMSConvergedDouble){
						for (unsigned i = 0; i < m_Input.size(); ++i) {
							RMSNow[i] = (RMSNow[i]+RMSPrev[0][i])/2;
						}
					}
					m_iPeriodRMS=0;
					if (bRMSTestTarget) {
						// check if RMS converged to target
						if(m_iAmpliConvergences<m_iMinAmpliConvergences)
							m_bConverged=false;
						for (unsigned i = 0; i < m_Input.size(); ++i) {
							if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
								m_bConverged = m_bConverged && (std::abs(RMSNow[i]-RMSTarget[i]) /
									RMSTarget[i] <= m_dTol);
								if(std::abs(RMSNow[i]-RMSTarget[i])<std::abs(RMSBest[i]-RMSTarget[i])){
									RMSBest=RMSNow;
									m_dAmplitudeBest=m_dAmplitude;
								}
							}
						}
						if (!m_bConverged) { // compute a new amplitude to get closer to target
#if HFELEM_COMMENTS == 1
							std::cerr << "\nAMPLITUDE CONVERGED - iteration " << m_iPeriod << std::endl;
							if(bRMSConvergedDouble){
								std::cerr << "   **** DOUBLE AMPLITUDE ****\n";
								for (unsigned i = 0; i < m_Input.size(); ++i) {
									if (m_Input[i].m_Flag & HFInput::HF_TEST) {
										std::cerr << "     RMS1:" << RMSPrev[1][i] << "  RMS2:" << RMSPrev[0][i] << std::endl;
									}
								}
							}
							else if(bRMSTooHigh)
								std::cerr << "   **** RMS PROBABLY TOO HIGH ****\n";
							else if(bRMSTooLow)
								std::cerr << "   **** RMS PROBABLY TOO LOW  ****\n";
#endif
							m_iAmpliConvergences++;
							doublereal dAmplitudeNew = 0.;
							for (unsigned i = 0; i < m_Input.size(); ++i) {
								if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
									if(m_iRMSrestarts>m_iMaxRMSrestarts){
										dAmplitudeNew = compute_amplitude(RMSNow[i], 0.0,
										                  m_dAmplitude, 0.0, RMSTarget[i], 0);
#if HFELEM_COMMENTS == 1
										std::cerr << "  too many restarts" << std::endl <<
										             "now - RMS:" << RMSNow[i]    << "  ampli:" << m_dAmplitude << std::endl <<
										             "new - target:" << RMSTarget[i] << "  ampli:" << dAmplitudeNew <<std::endl;
#endif
										continue;
									}
#if HFELEM_COMMENTS == 1
									if(bRMSLow)
										std::cerr << "low - RMS:" << RMSLow[i] << "  ampli:" << m_dAmplitudeLow << std::endl;
									if(bRMSUpp)
										std::cerr << "upp - RMS:" << RMSUpp[i] << "  ampli:" << m_dAmplitudeUpp << std::endl;
#endif
									
									if(RMSNow[i]<RMSTarget[i] && (!bRMSLow || (bRMSLow && RMSNow[i]>RMSLow[i]))){
										bRMSLow = true;
										m_dAmplitudeLow = m_dAmplitude;
										RMSLow = RMSNow;
										RMSLow[i] = (1.0-m_dRMSovershoot*m_dTol)*RMSNow[i];
#if HFELEM_COMMENTS == 1
										std::cerr << "new lower point - RMS:" << RMSLow[i] << "  ampli:" << m_dAmplitudeLow << std::endl;
#endif
									}else if(RMSNow[i]>RMSTarget[i] && (!bRMSUpp || (bRMSUpp && RMSNow[i]<RMSUpp[i]))){
										bRMSUpp = true;
										m_dAmplitudeUpp = m_dAmplitude;
										RMSUpp = RMSNow;
										RMSUpp[i] = (1.0+m_dRMSovershoot*m_dTol)*RMSNow[i];
#if HFELEM_COMMENTS == 1
										std::cerr << "new upper point - RMS: " << RMSUpp[i] << "  ampli:" << m_dAmplitudeUpp << std::endl;
#endif
									}else{
										bRMSLow=false;
										bRMSUpp=false;
#if HFELEM_COMMENTS == 1
										std::cerr << "reset lower and upper points" << std::endl;
#endif
									}
									if(bRMSLow && bRMSUpp && ((RMSUpp[i]-RMSLow[i])/std::max(RMSUpp[i],RMSLow[i])<2*m_dRMSovershoot*m_dTol || m_dAmplitudeUpp-m_dAmplitudeLow<m_dTol*m_dTol)){
										bRMSLow=false;
										bRMSUpp=false;
#if HFELEM_COMMENTS == 1
										std::cerr << "lower and upper points too close: reset" << std::endl;
#endif
									}
									if (bRMSLow && bRMSUpp){
#if HFELEM_COMMENTS == 1
										std::cerr << "  true true" << std::endl;
#endif
										dAmplitudeNew = compute_amplitude(RMSUpp[i], RMSLow[i],
										                  m_dAmplitudeUpp, m_dAmplitudeLow, RMSTarget[i], 1);
									}
									else{
										if (bRMSLow){
#if HFELEM_COMMENTS == 1
											std::cerr << "  only low" << std::endl;
#endif
										}else if(bRMSUpp){
#if HFELEM_COMMENTS == 1
											std::cerr << "  only upp" << std::endl;
#endif
										}else{
											m_iRMSrestarts++;
#if HFELEM_COMMENTS == 1
											std::cerr << "  restarting" << std::endl;
#endif
											if(m_iRMSrestarts>m_iMaxRMSrestarts)
												m_iPeriod = m_iMaxPeriods-m_iAfterMaxRMSrestartsPeriods;
										}
										dAmplitudeNew = compute_amplitude(RMSNow[i], 0.0,
										                  m_dAmplitude, 0.0, RMSTarget[i], 0);
									}
#if HFELEM_COMMENTS == 1
									std::cerr << "now - RMS:" << RMSNow[i] << "  ampli:" << m_dAmplitude << std::endl <<
									             "new - target:" << RMSTarget[i] << "  ampli:" << dAmplitudeNew <<std::endl;
#endif
								}
							}
							m_dAmplitude    = dAmplitudeNew;
							m_dAmplitudeOut = m_dAmplitude;
						}
					}
				}

				// update RMS for previous periods
				for (unsigned j = RMSPrev.size()-1; j > 0; --j){
					RMSPrev[j] = RMSPrev[j-1];
				}
				RMSPrev[0] = RMSNow;
			}
		} else {
			if(m_iPeriod > m_iConvergedPeriod+m_iWritePeriodsAfterConvergence){
				dErr = sqrt(dErr);
				if (dErr <= m_dTol) {
					m_bConverged = true;
				} else {
					m_bConverged = false;
				}
			}
		}
	}else{
		m_iPeriodRMS=0;
	}

	if (m_iPeriodCnt == m_iN) { // the period is over
		m_iPeriodCnt = 0;
		std::fill(RMSNow.begin(), RMSNow.end(), 0.);
		std::fill(SUM.begin(), SUM.end(), 0.);
		++m_iPeriod;
		if ( m_iPeriod >= m_iMinPeriods && (m_bConverged || (m_iMaxPeriods > 0 && m_iPeriod >= m_iMaxPeriods && m_NoConvStrategy == NoConv_CONTINUE )) && (m_iPeriod > m_iConvergedPeriod+m_iWritePeriodsAfterConvergence) ) {
			if (m_bConverged) {
				m_iPeriodOut =  m_iPeriod+m_iWritePeriodsAfterConvergence;
#if HFELEM_COMMENTS == 1
				std::cerr << "\n--- FREQUENCY CONVERGED " << std::setprecision(4) << m_dOmega/(2*M_PI) << " ---  iteration " << m_iPeriod;
				if(bRMSTest)
					for (unsigned i = 0; i < m_Input.size(); ++i)
						if (m_Input[i].m_Flag & HFInput::HF_OUTPUT)
							std::cerr << " ,  RMSvar" << i << " : " << RMSPrev[0][i];
				std::cerr << std::endl;
#endif
			} else {
				m_iPeriodOut = -m_iPeriod-m_iWritePeriodsAfterConvergence;
#if HFELEM_COMMENTS == 1
				std::cerr << "\n--- FREQUENCY NOT CONVERGED " << std::setprecision(4) <<m_dOmega/(2*M_PI) << " ---  iteration " << m_iPeriod;
				if(bRMSTest)
					for (unsigned i = 0; i < m_Input.size(); ++i)
						if (m_Input[i].m_Flag & HFInput::HF_OUTPUT)
							std::cerr << " ,  RMSvar" << i << " : " << RMSPrev[0][i];
				std::cerr << std::endl;
#endif
			}
			if(bRMSTest){
#if HFELEM_COMMENTS == 1
				if(bRMSConvergedDouble){
					std::cerr << "   **** DOUBLE AMPLITUDE ****\n";
					for (unsigned i = 0; i < m_Input.size(); ++i) {
						if (m_Input[i].m_Flag & HFInput::HF_TEST) {
							std::cerr << "     RMS1:" << RMSPrev[1][i] << "  RMS2:" << RMSPrev[2][i] << std::endl;
						}
					}
				}
#endif
				if(bRMSTestTarget){
					if (m_bConverged) {
						m_dAmplitudeOut = m_dAmplitude;
					} else {
						m_dAmplitudeOut = m_dAmplitudeBest;
						m_dAmplitude = m_dAmplitudeBest;
					}
#if HFELEM_COMMENTS == 1
					std::cerr << "Converged amplitudes: " << m_iAmpliConvergences << "/" << m_iMinAmpliConvergences << std::endl;
					std::cerr << "Upp Low Restarts: " << m_iRMSrestarts << "/" << m_iMaxRMSrestarts << std::endl;
					for (unsigned i = 0; i < m_Input.size(); ++i) {
						if (m_Input[i].m_Flag & HFInput::HF_TEST) {
							if(bRMSLow)
								std::cerr << "low - RMS:" << RMSLow[i] << "  ampli:" << m_dAmplitudeLow << std::endl;
							if(bRMSUpp)
								std::cerr << "upp - RMS:" << RMSUpp[i] << "  ampli:" << m_dAmplitudeUpp << std::endl;
							std::cerr << "now - RMS:"  << RMSPrev[0][i] << "  ampli:" << m_dAmplitude << std::endl;
							std::cerr << "best - RMS:" << RMSBest[i]   << "  ampli:" << m_dAmplitudeBest << std::endl;
							std::cerr << "target RMS:" << RMSTarget[i] << "  out ampli:" << m_dAmplitudeOut << std::endl;
						}
					}
#endif
					m_iRMSrestarts=0;
					m_iAmpliConvergences=0;
					std::fill(RMSBest.begin(), RMSBest.end(), 0.);
					bRMSUpp = false;
					bRMSLow = false;
				}
			}
#if HFELEM_COMMENTS == 1
			std::cerr << "----------" << std::endl << std::endl;
#endif
			m_bConverged = true;
			m_iConvergedPeriod = m_iPeriod;
		} else if ((m_iPeriod >= m_iMinPeriods) && (m_iMaxPeriods > 0) && (m_iPeriod > m_iMaxPeriods) && (m_NoConvStrategy == NoConv_ABORT)) {
			silent_cerr("HarmonicExcitationElem(" << GetLabel() << "): max periods " << m_iMaxPeriods << " exceeded for frequency " << m_dOmega << "; aborting..." << std::endl);

			// schedule for termination!
			m_bDone = true;
		}

		if(m_bConverged && m_iPeriod == m_iConvergedPeriod+m_iWritePeriodsAfterConvergence && !m_bDone){ // update omega
			m_dOmegaOut = m_dOmega;
			m_iOmegaCnt++;

			m_iPeriod = 0;
			m_iConvergedPeriod = 0;

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
			if (bRMSTestTarget) {
				for (unsigned i = 0; i < m_Input.size(); ++i) {
					if (m_Input[i].m_Flag & HFInput::HF_TARGET) {
						RMSTarget[i]=m_pDC_RMSTarget[i]->dGet(m_dOmega);
					}
				}
			}

			m_dDeltaT = (2*M_PI/m_dOmega)/m_iN;
			if(m_pDC_MaxDeltaT)
				m_dMaxDeltaT = m_pDC_MaxDeltaT->dGet(m_dOmega);
			else
				m_dMaxDeltaT = m_dDeltaT;
			if (m_dDeltaT > m_dMaxDeltaT) {
				m_iNumSubSteps = std::ceil(m_dDeltaT / m_dMaxDeltaT);
                if(((m_dDeltaT/m_dMaxDeltaT)-std::floor(m_dDeltaT/m_dMaxDeltaT))<m_dTimestepsCompareTol && ((m_dDeltaT/m_dMaxDeltaT)/std::floor(m_dDeltaT/m_dMaxDeltaT)-1.0)>0.0)
					m_iNumSubSteps -= 1;
				m_dDeltaT /= m_iNumSubSteps;
			} else {
				m_iNumSubSteps = 1;
			}
#if HFELEM_COMMENTS == 1
			if(m_dOmega>0 && !m_bDone)
				std::cerr << "omega: " << m_dOmega << " | DeltaT0: " << (2*M_PI/m_dOmega)/m_iN << " | MaxDeltaT: " << m_dMaxDeltaT << " | substeps: " << m_iNumSubSteps << " | DeltaT: " << m_dDeltaT << std::endl;
#endif
			m_bConverged = false;
		}

		if ((m_iPeriodOut!=0 && m_iPeriod == 0) || (m_bPrintAllPeriods && m_pDM->dGetTime()>m_dTInit)){
			m_bShouldBeOutput=true;
			m_iWriteAfterConvergenceCount=1;
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

	} else if (strcmp(s, "output") == 0) {
		return Priv_OUTPUT;
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

	case Priv_OUTPUT:
		return integer(m_bShouldBeOutput);

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
