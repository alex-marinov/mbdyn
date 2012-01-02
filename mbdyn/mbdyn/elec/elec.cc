/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Elementi elettrici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "elec.h"
#include "elecnode.h"
#include "drive.h"
#include "strnode.h"
#include "accelerometer.h"
#include "displacement.h"
#include "motor.h"
#include "dataman.h"
#include "discctrl.h"

/* Electric - begin */

Electric::Electric(unsigned int uL,
	const DofOwner* pDO, flag fOut)
: Elem(uL, fOut),
ElemWithDofs(uL, pDO, fOut)
{
	NO_OP;
}

Electric::~Electric(void)
{
	NO_OP;
}

/* Contributo al file di restart
 * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
 * relativa alla classe derivata */
std::ostream&
Electric::Restart(std::ostream& out) const
{
	return out << "  electric: " << GetLabel();
}


/* Tipo dell'elemento (usato solo per debug ecc.) */
Elem::Type
Electric::GetElemType(void) const
{
	return Elem::ELECTRIC;
}

/* Electric - end */

/* Legge un forgetting factor */

ForgettingFactor *
ReadFF(MBDynParser& HP, integer iNumOutputs)
{
	ForgettingFactor* pFF = 0;

	if (HP.IsKeyWord("forgetting" "factor")) {
		if (HP.IsKeyWord("const")) {
			doublereal d = HP.GetReal();

			SAFENEWWITHCONSTRUCTOR(pFF, ConstForgettingFactor,
				ConstForgettingFactor(d));

		} else if (HP.IsKeyWord("dynamic")) {
			integer n1 = HP.GetInt();
			integer n2 = HP.GetInt();
			doublereal dRho = HP.GetReal();
			doublereal dFact = HP.GetReal();
			doublereal dKRef = HP.GetReal();
			doublereal dKLim = HP.GetReal();

			SAFENEWWITHCONSTRUCTOR(pFF, DynamicForgettingFactor2,
				DynamicForgettingFactor2(n1, n2, iNumOutputs,
					dRho, dFact, dKRef, dKLim));

		} else {
			silent_cerr("Unknown forgetting factor type "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		/* default */
		SAFENEWWITHCONSTRUCTOR(pFF, ConstForgettingFactor,
			ConstForgettingFactor(1.));
	}

	ASSERT(pFF != 0);
	return pFF;
}

/* Legge un eccitatore persistente */

PersistentExcitation *
ReadPX(DataManager* pDM, MBDynParser& HP, integer iNumInputs)
{
	PersistentExcitation* pPX = 0;

	if (HP.IsKeyWord("excitation")) {
		if (iNumInputs == 1) {
			DriveCaller* pDC = HP.GetDriveCaller();
			SAFENEWWITHCONSTRUCTOR(pPX, ScalarPX, ScalarPX(pDC));

		} else {
			DriveCaller** ppDC = 0;
			SAFENEWARR(ppDC, DriveCaller*, iNumInputs);

			for (integer i = iNumInputs; i-- > 0; ) {
				ppDC[i] = HP.GetDriveCaller();
			}
			SAFENEWWITHCONSTRUCTOR(pPX, VectorPX, VectorPX(iNumInputs, ppDC));
		}

	} else {
		/* Null excitation */
		SAFENEW(pPX, NullPX);
	}

	ASSERT(pPX != 0);
	return pPX;
}

/* Legge un elemento elettrico */

Elem *
ReadElectric(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadElectric()");

	const char* sKeyWords[] = {
		"accelerometer",
		"displacement",
		"motor",
		"discrete" "control",
			"identification",
			"control",
			"adaptive" "control"
	};

	// enum delle parole chiave
	enum KeyWords {
		UNKNOWN = -1,

		ACCELEROMETER = 0,
		DISPLACEMENT,

		MOTOR,

		DISCRETECONTROL,
			IDENTIFICATION,
			CONTROL,
			ADAPTIVECONTROL,

		LASTKEYWORD
	};

	// tabella delle parole chiave
	KeyTable K(HP, sKeyWords);

	// lettura del tipo di elemento elettrico
	KeyWords CurrKeyWord = KeyWords(HP.GetWord());

	Elem* pEl = 0;

	switch (CurrKeyWord) {
	case ACCELEROMETER: {
		int f = 0;
		if (HP.IsKeyWord("translational")) {
			f = 1;
		} else if(HP.IsKeyWord("rotational")) {
			f = 2;
		}

		if (f) {
			// TODO: check if downgradable to ScalarNode

			// nodo strutturale collegato
			const StructNode* pStrNode
				= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

			// nodo astratto collegato
			const ScalarDifferentialNode* pAbsNode
				= pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);

			// Direzione
			Vec3 Dir;
			try {
				Dir = HP.GetUnitVecRel(ReferenceFrame(pStrNode));
			} catch (ErrNullNorm) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"null direction "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
			}
			DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);

			// offset
			Vec3 Tmpf(Zero3);
			if (HP.IsKeyWord("position") || HP.IsKeyWord("offset") || f == 1) {
				Tmpf = HP.GetPosRel(ReferenceFrame(pStrNode));
			}
			flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

			switch (f) {
			case 1:
				SAFENEWWITHCONSTRUCTOR(pEl, TranslAccel,
					TranslAccel(uLabel, pDO,
						pStrNode, pAbsNode,
						Dir, Tmpf, fOut));
				break;

			case 2:
				SAFENEWWITHCONSTRUCTOR(pEl, RotAccel,
					RotAccel(uLabel, pDO,
						pStrNode, pAbsNode,
						Dir, fOut));
				break;

			default:
				silent_cerr("you shouldn't be here!" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else {

			// nodo strutturale collegato
			const StructNode* pStrNode
				= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

			// nodo astratto collegato
			const ScalarDifferentialNode* pAbsNode
				= pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);

			// Direzione
			Vec3 Dir;
			try {
				Dir = HP.GetUnitVecRel(ReferenceFrame(pStrNode));
			} catch (ErrNullNorm) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"null direction "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
			}
			DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);

			// Parametri
			doublereal dOmega = HP.GetReal();
			if (dOmega <= 0.) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"illegal Omega "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			doublereal dTau = HP.GetReal();
			if (dTau <= 0.) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"illegal Tau "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			doublereal dCsi = HP.GetReal();
			if (dCsi <= 0. || dCsi > 1.) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"illegal Csi "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			doublereal dKappa = HP.GetReal();
			if (dKappa == 0.) {
				silent_cerr("Accelerometer(" << uLabel << "): "
					"illegal Kappa "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			DEBUGCOUT("Omega: " << dOmega
				<< ", Tau: " << dTau
				<< ", Csi: " << dCsi
				<< ", Kappa: " << dKappa
				<< std::endl);

			flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

			SAFENEWWITHCONSTRUCTOR(pEl, Accelerometer,
				Accelerometer(uLabel, pDO, pStrNode, pAbsNode,
					Dir, dOmega, dTau, dCsi, dKappa,
					fOut));
		}

		} break;

	case DISPLACEMENT: {
		// nodo strutturale collegato 1
		const StructNode* pStrNode1
			= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		// offset 1
		Vec3 Tmpf1(Zero3);
		if (HP.IsKeyWord("position")) {
			NO_OP;
		}
		Tmpf1 = HP.GetPosRel(ReferenceFrame(pStrNode1));

		// nodo strutturale collegato 2
		const StructNode* pStrNode2
			= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		// offset 2
		Vec3 Tmpf2(Zero3);
		if (HP.IsKeyWord("position")) {
			NO_OP;
		}
		Tmpf2 = HP.GetPosRel(ReferenceFrame(pStrNode2));

		// nodo astratto collegato
		const ScalarDifferentialNode* pAbsNode
			= pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);

		flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

		SAFENEWWITHCONSTRUCTOR(pEl, DispMeasure,
			DispMeasure(uLabel, pDO,
				pStrNode1, pStrNode2, pAbsNode,
				Tmpf1, Tmpf2, fOut));
		} break;

	case MOTOR: {
		// nodo strutturale collegato 1
		const StructNode* pStrNode1
			= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		// nodo strutturale collegato 2
		const StructNode* pStrNode2
			= pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		// direzione
		Vec3 TmpDir;
		try {
			TmpDir = HP.GetUnitVecRel(ReferenceFrame(pStrNode1));
		} catch (ErrNullNorm) {
			silent_cerr("Motor(" << uLabel << "): "
				"illegal motor direction "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
		}

		// nodo elettrico1 collegato
		ElectricNode* pVoltage1
			= dynamic_cast<ElectricNode *>(pDM->ReadNode(HP, Node::ELECTRIC));

		// nodo elettrico2 collegato
		ElectricNode* pVoltage2
			= dynamic_cast<ElectricNode *>(pDM->ReadNode(HP, Node::ELECTRIC));

		doublereal dG = HP.GetReal();
		doublereal dl = HP.GetReal();
		doublereal dr = HP.GetReal();

		flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

		SAFENEWWITHCONSTRUCTOR(pEl, Motor,
			Motor(uLabel, pDO,
				pStrNode1, pStrNode2, pVoltage1, pVoltage2,
				TmpDir, dG, dl, dr, fOut));
		} break;

	case DISCRETECONTROL: {
		// Dati generali del controllore
		integer iNumOutputs = HP.GetInt();
		integer iNumInputs = HP.GetInt();
		integer iOrderA = HP.GetInt();
		integer iOrderB = iOrderA;
		if (HP.IsKeyWord("fir")) {
			iOrderB = HP.GetInt();
		}

		integer iNumIter = HP.GetInt();

		DEBUGCOUT("Discrete controller of order " << iOrderA);
		if (iOrderB != iOrderA) {
			DEBUGCOUT(" (fir order " << iOrderB << ')' << std::endl);
		}
		DEBUGCOUT(": " << iNumOutputs << " output(s) and "
			<< iNumInputs << " input(s)" << std::endl
			<< "Update every " << iNumIter << " iterations" << std::endl);

		// Tipo di controllo
		DiscreteControlProcess* pDCP = 0;
		switch (HP.GetWord()) {
		case CONTROL: {
			// Add the file with control data
			const char* s = HP.GetFileName();
			std::string infile = s;

			DEBUGCOUT("Getting control matrices "
				"from file \"" << infile << "\""
				<< std::endl);

			/* Construction of controller */
			SAFENEWWITHCONSTRUCTOR(pDCP,
				DiscreteControlARXProcess_Debug,
				DiscreteControlARXProcess_Debug(iNumOutputs,
					iNumInputs,
					iOrderA,
					iOrderB,
					infile));
			} break;

		case IDENTIFICATION: {
			unsigned f_proc = DiscreteControlProcess::DISCPROC_UNKNOWN;
			if (HP.IsKeyWord("arx")) {
				f_proc = DiscreteControlProcess::DISCPROC_ARX;

			} else if (HP.IsKeyWord("armax")) {
				f_proc = DiscreteControlProcess::DISCPROC_ARMAX;

			} else {
				silent_cerr("DiscreteControl(" << uLabel << "): "
					"unknown identification type "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			// Forgetting factor
			ForgettingFactor* pFF = ReadFF(HP, iNumOutputs);

			// Persistent excitation
			PersistentExcitation* pPX = ReadPX(pDM, HP, iNumInputs);

			std::string outfile;
			if (HP.IsKeyWord("file")) {
				const char *s = HP.GetFileName();
				outfile = s;
			}

			/* Construction of controller */
			SAFENEWWITHCONSTRUCTOR(pDCP,
				DiscreteIdentProcess_Debug,
				DiscreteIdentProcess_Debug(iNumOutputs,
					iNumInputs,
					iOrderA,
					iOrderB,
					pFF, pPX,
					f_proc, outfile));
	   		} break;

		case ADAPTIVECONTROL: {
#ifdef USE_DBC
			unsigned f_proc = DiscreteControlProcess::DISCPROC_ARX;
			doublereal dPeriodicFactor(0.);

			if (HP.IsKeyWord("arx")) {
				DEBUGCOUT("ARX adaptive control" << std::endl);
				f_proc = DiscreteControlProcess::DISCPROC_ARX;

			} else if (HP.IsKeyWord("armax")) {
				DEBUGCOUT("ARMAX adaptive control" << std::endl);
				f_proc = DiscreteControlProcess::DISCPROC_ARMAX;
			}

			if (HP.IsKeyWord("periodic")) {
				dPeriodicFactor = HP.GetReal();
			}

			GPCDesigner* pCD = 0;
			if (HP.IsKeyWord("gpc")) {
				DEBUGCOUT("GPC adaptive control" << std::endl);

				integer iPredS = HP.GetInt();
				integer iContrS = HP.GetInt();
				integer iPredH = HP.GetInt();
				integer iContrH = 0;

				DEBUGCOUT("prediction advancing horizon: " << iPredS << std::endl
					<< "prediction receding horizon: " << iPredH << std::endl
					<< "control advancing horizon: " << iContrS << std::endl
					<< "control receding horizon: " << iContrH << std::endl);

				if (iPredS < 0) {
					silent_cerr("DiscreteControl(" << uLabel << "): "
						"prediction advancing horizon "
						"must be positive "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (iPredH < 0) {
					silent_cerr("DiscreteControl(" << uLabel << "): "
						"prediction receding horizon "
						"must be positive "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (iPredH >= iPredS) {
					silent_cerr("DiscreteControl(" << uLabel << "): "
						"prediction receding horizon "
						"must be less than "
						"prediction advancing horizon "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (iContrS < 0) {
					silent_cerr("DiscreteControl(" << uLabel << "): "
						"control advancing horizon "
						"must be positive "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				doublereal* pW = 0;
				doublereal* pR = 0;
				SAFENEWARR(pW, doublereal, iPredS - iPredH);
				SAFENEWARR(pR, doublereal, iContrS - iContrH);

				if (HP.IsKeyWord("prediction" "weights")) {
					DEBUGCOUT("prediction weights:" << std::endl);
					for (integer i = iPredS - iPredH; i-- > 0; ) {
						pW[i] = HP.GetReal();
						DEBUGCOUT("W[" << i + 1 << "] = " << pW[i] << std::endl);
					}

				} else {
					for (integer i = 0; i < iPredS - iPredH; i++) {
						pW[i] = 1.;
					}
				}

				if (HP.IsKeyWord("control" "weights")) {
					DEBUGCOUT("control weights:" << std::endl);
					for (integer i = iContrS - iContrH; i-- > 0; ) {
						pR[i] = HP.GetReal();
						DEBUGCOUT("R[" << i + 1 << "] = " << pR[i] << std::endl);
					}
				} else {
					for (integer i = 0; i < iContrS-iContrH; i++) {
						pR[i] = 1.;
					}
				}

				DEBUGCOUT("Weight Drive:" << std::endl);
				DriveCaller* pLambda = HP.GetDriveCaller();

				SAFENEWWITHCONSTRUCTOR(pCD, GPC,
					GPC(iNumOutputs, iNumInputs,
						iOrderA, iOrderB,
						iPredS, iContrS,
						iPredH, iContrH,
						pW, pR, pLambda,
						dPeriodicFactor, f_proc));

			} else if (HP.IsKeyWord("deadbeat")) {
				DEBUGCOUT("DeadBeat adaptive control" << std::endl);

				int iPredS = HP.GetInt();
				int iContrS = HP.GetInt();
				SAFENEWWITHCONSTRUCTOR(pCD, DeadBeat,
					DeadBeat(iNumOutputs, iNumInputs,
						iOrderA, iOrderB,
						iPredS, iContrS,
						dPeriodicFactor, f_proc));

			} else {
				silent_cerr("DiscreteControl(" << uLabel << "): "
					"unknown adaptive control type at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			// Forgetting factor
			DEBUGCOUT("Forgetting Factor:" << std::endl);
			ForgettingFactor* pFF = ReadFF(HP, iNumOutputs);

			// Persistent excitation
			DEBUGCOUT("Persistent Excitation:" << std::endl);
			PersistentExcitation* pPX = ReadPX(pDM, HP, iNumInputs);

			DriveCaller* pTrig = 0;
			if (HP.IsKeyWord("trigger")) {
				DEBUGCOUT("Trigger:" << std::endl);
				pTrig = HP.GetDriveCaller();
			} else {
				SAFENEW(pTrig, OneDriveCaller);
			}

			// desired output
			std::vector<DriveCaller*> vDesiredOut;
			if (HP.IsKeyWord("desired" "output")) {
				DEBUGCOUT("Desired output:" << std::endl);
				vDesiredOut.resize(iNumOutputs);

				for (integer i = 0; i < iNumOutputs; i++) {
					DEBUGCOUT("output[" << i + 1 << "]:" << std::endl);
					vDesiredOut[i] = HP.GetDriveCaller();
				}
			}

			std::string outfile;
			if (HP.IsKeyWord("file") || HP.IsKeyWord("output" "file")) {
				const char *s = HP.GetFileName();
				outfile = s;
				DEBUGCOUT("Identified matrices will be output in file \"" << s << "\"" << std::endl);
			}

			/* Construction of controller */
			ASSERT(f_proc == DiscreteControlProcess::DISCPROC_ARX
				|| f_proc == DiscreteControlProcess::DISCPROC_ARMAX);
			SAFENEWWITHCONSTRUCTOR(pDCP, DAC_Process_Debug,
				DAC_Process_Debug(iNumOutputs, iNumInputs,
					iOrderA, iOrderB,
					pFF, pCD, pPX, pTrig, vDesiredOut,
					outfile, f_proc));
#else /* !USE_DBC */
			silent_cerr("DiscreteControl(" << uLabel << "): "
				"GPC/deadbeat control is not available "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_DBC */
			} break;

		default:
			silent_cerr("Sorry, not implemented yed" << std::endl);
			throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
		}



		if (!HP.IsKeyWord("outputs")) {
			silent_cerr("DiscreteControl(" << uLabel << "): "
				"\"outputs\" expected "
				"at line " << HP.GetLineData()
				<< std::endl);

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::vector<ScalarValue *> vOutputs(iNumOutputs);
		std::vector<DriveCaller *> vOutScaleFact(iNumOutputs);

		// Allocazione nodi e connessioni
		for (int i = 0; i < iNumOutputs; i++) {
			vOutputs[i] = ReadScalarValue(pDM, HP);
			if (HP.IsKeyWord("scale")) {
				vOutScaleFact[i] = HP.GetDriveCaller();

			} else {
				vOutScaleFact[i] = 0;
				SAFENEW(vOutScaleFact[i], OneDriveCaller);
			}
		}

		if (!HP.IsKeyWord("inputs")) {
			silent_cerr("DiscreteControl(" << uLabel << "): "
				"\"inputs\" expected "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		// Same thing for input nodes
		ScalarDof* pInputs = 0;
		SAFENEWARRNOFILL(pInputs, ScalarDof, iNumInputs);

		// Allocazione nodi e connessioni
		for (int i = 0; i < iNumInputs; i++) {
			pInputs[i] = ReadScalarDof(pDM, HP, false, true);
		}

		flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

		SAFENEWWITHCONSTRUCTOR(pEl, DiscreteControlElem,
			DiscreteControlElem(uLabel, pDO,
				iNumOutputs, vOutputs, vOutScaleFact,
				iNumInputs, pInputs,
				pDCP, iNumIter, fOut));
		} break;

	// Aggiungere altri elementi elettrici

	default:
		silent_cerr("Electric(" << uLabel << "): "
			"unknown type "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// Se non c'e' il punto e virgola finale
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* ReadElectric() */

