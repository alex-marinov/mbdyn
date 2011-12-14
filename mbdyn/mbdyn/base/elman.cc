/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

/* ElementManager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "dataman.h"
#include "search.h"
#include "gravity.h"
#include "aerodyn.h"
#include "driven.h"

#include "aeroelem.h"
#include "beam.h"

/* DataManager - begin */

/* costruttore: resetta i dati */
void
DataManager::ElemManager(void)
{
	/* Reset della struttura ElemData */
	for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
		ElemData[i].iExpectedNum = 0;
		ElemData[i].DofOwnerType = DofOwner::UNKNOWN;
		ElemData[i].uFlags = 0U;
		ElemData[i].DefaultOut(::fDefaultOut == 1); /* Da "output.h" */
		ElemData[i].OutFile = OutputHandler::UNKNOWN;	/* "output.h" */
		ElemData[i].uOutputFlags = 0U;
		ElemData[i].od = EULER_123;
	}

	/* Add dof type */
	ElemData[Elem::JOINT].DofOwnerType = DofOwner::JOINT;
	ElemData[Elem::PLATE].DofOwnerType = DofOwner::PLATE;
	ElemData[Elem::ELECTRIC].DofOwnerType = DofOwner::ELECTRIC;
	ElemData[Elem::THERMAL].DofOwnerType = DofOwner::THERMAL;
	ElemData[Elem::ELECTRICBULK].DofOwnerType = DofOwner::ELECTRICBULK;
	ElemData[Elem::GENEL].DofOwnerType = DofOwner::GENEL;
	ElemData[Elem::INDUCEDVELOCITY].DofOwnerType = DofOwner::INDUCEDVELOCITY;
	ElemData[Elem::AEROMODAL].DofOwnerType = DofOwner::AEROMODAL;
	ElemData[Elem::AERODYNAMIC].DofOwnerType = DofOwner::AERODYNAMIC;
	ElemData[Elem::HYDRAULIC].DofOwnerType = DofOwner::HYDRAULIC;
	ElemData[Elem::LOADABLE].DofOwnerType = DofOwner::LOADABLE;

	/* Add file type */
	ElemData[Elem::AUTOMATICSTRUCTURAL].OutFile = OutputHandler::INERTIA;
	ElemData[Elem::AUTOMATICSTRUCTURAL].Desc = "AutomaticStructural";
	ElemData[Elem::AUTOMATICSTRUCTURAL].ShortDesc = "autostruct";

	ElemData[Elem::JOINT].OutFile = OutputHandler::JOINTS;
	ElemData[Elem::JOINT].Desc = "Joint";
	ElemData[Elem::JOINT].ShortDesc = "joint";

	ElemData[Elem::FORCE].OutFile = OutputHandler::FORCES;
	ElemData[Elem::FORCE].Desc = "Force";
	ElemData[Elem::FORCE].ShortDesc = "force";

	ElemData[Elem::BEAM].OutFile = OutputHandler::BEAMS;
	ElemData[Elem::BEAM].Desc = "Beam";
	ElemData[Elem::BEAM].ShortDesc = "beam";
	ElemData[Elem::BEAM].uOutputFlags = Beam::OUTPUT_DEFAULT;

	ElemData[Elem::PLATE].OutFile = OutputHandler::PLATES;
	ElemData[Elem::PLATE].Desc = "Plate";
	ElemData[Elem::PLATE].ShortDesc = "plate";

	ElemData[Elem::AIRPROPERTIES].OutFile = OutputHandler::AIRPROPS;
	ElemData[Elem::AIRPROPERTIES].Desc = "AirProperties";
	ElemData[Elem::AIRPROPERTIES].ShortDesc = "airprops";

	ElemData[Elem::GRAVITY].OutFile = OutputHandler::GRAVITY;
	ElemData[Elem::GRAVITY].Desc = "Gravity";
	ElemData[Elem::GRAVITY].ShortDesc = "gravity";
	ElemData[Elem::GRAVITY].DefaultOut(false);

	ElemData[Elem::INDUCEDVELOCITY].OutFile = OutputHandler::ROTORS;
	ElemData[Elem::INDUCEDVELOCITY].Desc = "InducedVelocity";
	ElemData[Elem::INDUCEDVELOCITY].ShortDesc = "indvel";

	ElemData[Elem::AEROMODAL].OutFile = OutputHandler::AEROMODALS;
	ElemData[Elem::AEROMODAL].Desc = "AerodynamicModal";
	ElemData[Elem::AEROMODAL].ShortDesc = "aeromodal";

	ElemData[Elem::AERODYNAMIC].OutFile = OutputHandler::AERODYNAMIC;
	ElemData[Elem::AERODYNAMIC].Desc = "Aerodynamic";
	ElemData[Elem::AERODYNAMIC].ShortDesc = "aero";
	ElemData[Elem::BEAM].uOutputFlags = AerodynamicOutput::OUTPUT_DEFAULT;

	ElemData[Elem::HYDRAULIC].OutFile = OutputHandler::HYDRAULIC;
	ElemData[Elem::HYDRAULIC].Desc = "Hydraulic";
	ElemData[Elem::HYDRAULIC].ShortDesc = "hydr";

	ElemData[Elem::LOADABLE].OutFile = OutputHandler::LOADABLE;
	ElemData[Elem::LOADABLE].Desc = "Loadable";
	ElemData[Elem::LOADABLE].ShortDesc = "loadable";

	ElemData[Elem::GENEL].OutFile = OutputHandler::GENELS;
	ElemData[Elem::GENEL].Desc = "Genel";
	ElemData[Elem::GENEL].ShortDesc = "genel";

	ElemData[Elem::EXTERNAL].OutFile = OutputHandler::EXTERNALS;
	ElemData[Elem::EXTERNAL].Desc = "External";
	ElemData[Elem::EXTERNAL].ShortDesc = "ext";

	ElemData[Elem::THERMAL].OutFile = OutputHandler::THERMALELEMENTS;
	ElemData[Elem::THERMAL].Desc = "Thermal";
	ElemData[Elem::THERMAL].ShortDesc = "thermal";

	/* Inheritance scheme */
	ElemData[Elem::AUTOMATICSTRUCTURAL].iDerivation = ELEM;
	ElemData[Elem::GRAVITY].iDerivation = ELEM;
	ElemData[Elem::BODY].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
	ElemData[Elem::JOINT].iDerivation = ELEM | GRAVITYOWNER | DOFOWNER | INITIALASSEMBLY;
	ElemData[Elem::GENEL].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::FORCE].iDerivation = ELEM | INITIALASSEMBLY;
	ElemData[Elem::BEAM].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
	ElemData[Elem::PLATE].iDerivation = ELEM | GRAVITYOWNER | DOFOWNER | INITIALASSEMBLY;
	ElemData[Elem::AIRPROPERTIES].iDerivation = ELEM | INITIALASSEMBLY;
	ElemData[Elem::INDUCEDVELOCITY].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
	ElemData[Elem::AEROMODAL].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
	ElemData[Elem::AERODYNAMIC].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER | INITIALASSEMBLY;
	ElemData[Elem::ELECTRIC].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::THERMAL].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::ELECTRICBULK].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::BULK].iDerivation = ELEM;
	ElemData[Elem::LOADABLE].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY | DOFOWNER;
	ElemData[Elem::DRIVEN].iDerivation = ELEM;
	ElemData[Elem::INERTIA].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;

	/* Types that must be unique */
	ElemData[Elem::GRAVITY].IsUnique(true);
	ElemData[Elem::AIRPROPERTIES].IsUnique(true);

	/* Types that are used by default in initial assembly */
	ElemData[Elem::JOINT].ToBeUsedInAssembly(true);
	ElemData[Elem::JOINT_REGULARIZATION].ToBeUsedInAssembly(true);
	ElemData[Elem::BEAM].ToBeUsedInAssembly(true);
	ElemData[Elem::PLATE].ToBeUsedInAssembly(true);

	/* Aggiungere qui se un tipo genera forze d'inerzia e quindi
	 * deve essere collegato all'elemento accelerazione di gravita' */
	ElemData[Elem::BODY].GeneratesInertiaForces(true);
	ElemData[Elem::JOINT].GeneratesInertiaForces(true);
	ElemData[Elem::LOADABLE].GeneratesInertiaForces(true);
	ElemData[Elem::PLATE].GeneratesInertiaForces(true);

	/* Aggiungere qui se un tipo usa le proprieta' dell'aria e quindi
	 * deve essere collegato all'elemento proprieta' dell'aria */
	ElemData[Elem::INDUCEDVELOCITY].UsesAirProperties(true);
	ElemData[Elem::AEROMODAL].UsesAirProperties(true);
	ElemData[Elem::AERODYNAMIC].UsesAirProperties(true);
	ElemData[Elem::LOADABLE].UsesAirProperties(true);
	ElemData[Elem::EXTERNAL].UsesAirProperties(true);

	/* Reset della struttura DriveData */
	for (int i = 0; i < Drive::LASTDRIVETYPE; i++) {
		DriveData[i].ppFirstDrive = NULL;
		DriveData[i].iNum = 0;
	}
}

/* distruttore */
void
DataManager::ElemManagerDestructor(void)
{
	DEBUGCOUT("Entering DataManager::ElemManagerDestructor()" << std::endl);

	/* Distruzione matrici di lavoro per assemblaggio */
	if (pWorkMatB != NULL) {
		DEBUGCOUT("deleting assembly structure, SubMatrix B" << std::endl);
		SAFEDELETE(pWorkMatB);
	}

	if (pWorkMatA != NULL) {
		DEBUGCOUT("deleting assembly structure, SubMatrix A" << std::endl);
		SAFEDELETE(pWorkMatA);
	}

	if (pWorkVec != NULL) {
		DEBUGCOUT("deleting assembly structure, SubVector" << std::endl);
		SAFEDELETE(pWorkVec);
	}

	for (ElemVecType::iterator p = Elems.begin(); p != Elems.end(); ++p) {
		DEBUGCOUT("deleting element "
			<< psElemNames[(*p)->GetElemType()]
			<< "(" << (*p)->GetLabel() << ")"
			<< std::endl);
		SAFEDELETE(*p);
	}

	/* Distruzione drivers */
	if (ppDrive != NULL) {
		Drive** pp = ppDrive;
		while (pp < ppDrive+iTotDrive) {
			if (*pp != NULL) {
				DEBUGCOUT("deleting driver "
						<< (*pp)->GetLabel()
						<< ", type "
						<< psDriveNames[(*pp)->GetDriveType()]
						<< std::endl);
				SAFEDELETE(*pp);
			}
			pp++;
		}

		DEBUGCOUT("deleting drivers structure" << std::endl);
		SAFEDELETEARR(ppDrive);
	}
}


/* Inizializza la struttura dei dati degli elementi ed alloca
 * l'array degli elementi */
void
DataManager::ElemDataInit(void)
{
	unsigned iTotElem = 0;

	/* struttura degli elementi */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		iTotElem += ElemData[iCnt].iExpectedNum;
	}

	DEBUGCOUT("iTotElem = " << iTotElem << std::endl);

	/* FIXME: reverse this:
	 * - read and populate ElemData[iCnt].ElemContainer first
	 * - then create Elems and fill it with Elem*
	 */
	if (iTotElem > 0) {
		Elems.resize(iTotElem);

		/* Inizializza l'iteratore degli elementi usato all'interno
		 * dell'ElemManager */
		ElemIter.Init(&Elems[0], iTotElem);

		for (unsigned iCnt = 0; iCnt < iTotElem; iCnt++) {
			Elems[iCnt] = 0;
		}

	} else {
		silent_cerr("warning, no elements are defined" << std::endl);
	}

	/* struttura dei drivers */
	for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE; iCnt++) {
		iTotDrive += DriveData[iCnt].iNum;
	}

	DEBUGCOUT("iTotDrive = " << iTotDrive << std::endl);

	if (iTotDrive > 0) {
		SAFENEWARR(ppDrive, Drive*, iTotDrive);

		/* Puo' essere sostituito con un opportuno iteratore: */
#if 0
		VecIter<Drive*> DriveIter(ppDrive, iTotDrive);
		Drive* pTmp = NULL;
		if (DriveIter.bGetFirst(pTmp)) {
			do {
				pTmp = NULL;
			} while (DriveIter.bGetNext(pTmp));
		}
#endif

		Drive** ppTmp = ppDrive;
		while (ppTmp < ppDrive + iTotDrive) {
			*ppTmp++ = NULL;
		}

		DriveData[0].ppFirstDrive = ppDrive;
		for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE-1; iCnt++) {
			DriveData[iCnt+1].ppFirstDrive =
				DriveData[iCnt].ppFirstDrive +
				DriveData[iCnt].iNum;
		}

	} else {
		DEBUGCERR("warning, no drivers are defined" << std::endl);
	}
}

/* Prepara per l'assemblaggio */
void
DataManager::ElemAssInit(void)
{
	Elem* pE = 0;
	if (ElemIter.bGetFirst(pE)) {
		do {
			/* Aggiorna le dimensioni massime degli spazi di lavoro */
			integer iNumRows = 0;
			integer iNumCols = 0;
			pE->WorkSpaceDim(&iNumRows, &iNumCols);
			if (iNumRows > iMaxWorkNumRows) {
				iMaxWorkNumRows = iNumRows;
				DEBUGLCOUT(MYDEBUG_INIT, "Current max work rows number: "
					<< iMaxWorkNumRows << std::endl);
			}
			if (iNumCols > iMaxWorkNumCols) {
				iMaxWorkNumCols = iNumCols;
				DEBUGLCOUT(MYDEBUG_INIT, "Current max work cols number: "
					<< iMaxWorkNumCols << std::endl);
			}
		} while (ElemIter.bGetNext(pE));
	}

	ASSERT(iMaxWorkNumRows > 0);
	ASSERT(iMaxWorkNumCols > 0);

	/* Per evitare problemi, alloco tanto spazio quanto necessario per
	 * scrivere in modo sparso la matrice piu' grande
	 *
	 * Nota: le dimensioni sono state moltiplicate per due per
	 * poter creare due matrici (in quanto la seconda e'
	 * richiesta per gli autovalori)
	 */
	iWorkIntSize = 2*iMaxWorkNumRows*iMaxWorkNumCols;
	iWorkDoubleSize = iMaxWorkNumRows*iMaxWorkNumCols;

	if (iWorkIntSize > 0) {

		/* SubMatrixHandlers */
		SAFENEWWITHCONSTRUCTOR(pWorkMatA,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iMaxWorkNumRows,
					iMaxWorkNumCols));

		SAFENEWWITHCONSTRUCTOR(pWorkMatB,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iMaxWorkNumRows,
					iMaxWorkNumCols));

		pWorkMat = pWorkMatA;

		SAFENEWWITHCONSTRUCTOR(pWorkVec,
				MySubVectorHandler,
				MySubVectorHandler(iMaxWorkNumRows));

		DEBUGCOUT("Creating working matrices:" << iMaxWorkNumRows
				<< " x " << iMaxWorkNumCols << std::endl);

	} else {
		silent_cerr("warning, null size of working matrix"
				<< std::endl);
	}
}

/* Assemblaggio dello jacobiano.
 * Di questa routine e' molto importante l'efficienza, quindi vanno valutate
 * correttamente le opzioni. */
void
DataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	ASSERT(pWorkMat != NULL);
	ASSERT(Elems.begin() != Elems.end());

	AssJac(JacHdl, dCoef, ElemIter, *pWorkMat);
}

void
DataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef,
		VecIter<Elem *> &Iter,
		VariableSubMatrixHandler& WorkMat)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	JacHdl.Reset();

#if 0
	for (int r = 1; r <= JacHdl.iGetNumRows(); r++) {
		for (int c = 1; c <= JacHdl.iGetNumCols(); c++) {
			if (JacHdl(r, c) != 0.) {
				silent_cerr("JacHdl(" << r << "," << c << ")="
						<< JacHdl(r, c) << std::endl);
			}
			JacHdl(r, c) = std::numeric_limits<doublereal>::epsilon();
		}
	}
#endif

	Elem* pTmpEl = NULL;
	if (Iter.bGetFirst(pTmpEl)) {
		do {
			try {
				JacHdl += pTmpEl->AssJac(WorkMat, dCoef,
						*pXCurr, *pXPrimeCurr);
			}
			catch (ErrDivideByZero) {
				silent_cerr("AssJac: divide by zero "
					"in " << psElemNames[pTmpEl->GetElemType()]
					<< "(" << pTmpEl->GetLabel() << ")"
					<< std::endl);
				throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
			}
		} while (Iter.bGetNext(pTmpEl));
	}
}

void
DataManager::AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl)
{
	DEBUGCOUT("Entering DataManager::AssMats()" << std::endl);

	ASSERT(pWorkMatA != NULL);
	ASSERT(pWorkMatB != NULL);
	ASSERT(Elems.begin() != Elems.end());

	AssMats(A_Hdl, B_Hdl, ElemIter, *pWorkMatA, *pWorkMatB);
}

void
DataManager::AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl,
		VecIter<Elem *> &Iter,
		VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB)
{
	DEBUGCOUT("Entering DataManager::AssMats()" << std::endl);

	/* Versione con iteratore: */
	Elem* pTmpEl = NULL;
	if (Iter.bGetFirst(pTmpEl)) {
		/* Nuova versione, piu' compatta.
		 * La funzione propria AssJac, comune a tutti gli elementi,
		 * scrive nella WorkMat (passata come reference) il contributo
		 * dell'elemento allo jacobiano e restituisce un reference
		 * alla workmat stessa, che viene quindi sommata allo jacobiano.
		 * Ogni elemento deve provvedere al resizing della WorkMat e al
		 * suo reset ove occorra */

		/* il SubMatrixHandler e' stato modificato in modo da essere
		 * in grado di trasformarsi agevolmente da Full a Sparse
		 * e quindi viene gestito in modo automatico, e del tutto
		 * trasparente, il passaggio da un modo all'altro.
		 * L'elemento switcha la matrice nel modo che ritiene
		 * opportuno; l'operatore += capisce di quale matrice
		 * si sta occupando ed agisce di conseguenza.
		 */

		/* Con VariableSubMatrixHandler */
		do {
			pTmpEl->AssMats(WorkMatA, WorkMatB,
					*pXCurr, *pXPrimeCurr);
			A_Hdl += WorkMatA;
			B_Hdl += WorkMatB;
		} while (Iter.bGetNext(pTmpEl));
	}
}

/* Assemblaggio del residuo */
void
DataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef) 
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	AssRes(ResHdl, dCoef, ElemIter, *pWorkVec);
}

void
DataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef,
		VecIter<Elem *> &Iter,
		SubVectorHandler& WorkVec)
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	Elem* pTmpEl = NULL;
	bool ChangedEqStructure(false);
	if (Iter.bGetFirst(pTmpEl)) {
		do {
			try {
				ResHdl += pTmpEl->AssRes(WorkVec, dCoef,
					*pXCurr, *pXPrimeCurr);
			}
			catch(Elem::ChangedEquationStructure) {
				ResHdl += WorkVec;
				ChangedEqStructure = true;
			}
			catch (ErrDivideByZero) {
				silent_cerr("AssRes: divide by zero "
					"in " << psElemNames[pTmpEl->GetElemType()]
					<< "(" << pTmpEl->GetLabel() << ")"
					<< std::endl);
				throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
			}

#if 0
			silent_cerr("### " << psElemNames[pTmpEl->GetElemType()] << "(" << pTmpEl->GetLabel() <<"):" << std::endl);
			PrintResidual(ResHdl, -1);
#endif
		} while (Iter.bGetNext(pTmpEl));
	}
	if (ChangedEqStructure) {
		throw ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

void
DataManager::ElemOutputPrepare(OutputHandler& OH)
{
#ifdef USE_NETCDF
	for (unsigned et = 0; et < Elem::LASTELEMTYPE; et++) {
		if (ElemData[et].ElemContainer.size() && OH.UseNetCDF(ElemData[et].OutFile)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));
			ASSERT(ElemData[et].Desc != 0);
			ASSERT(ElemData[et].ShortDesc != 0);

			integer iNumElems = ElemData[et].ElemContainer.size();

			OutputHandler::AttrValVec attrs(1);
			attrs[0] = OutputHandler::AttrVal("description", std::string(ElemData[et].Desc) + " elements labels");

			OutputHandler::NcDimVec dim(1);
			dim[0] = OH.CreateDim(std::string(ElemData[et].ShortDesc) + "_elem_labels_dim", iNumElems);

			NcVar *VarLabels = OH.CreateVar(std::string("elem.") + ElemData[et].ShortDesc, ncInt, attrs, dim);

			ElemContainerType::const_iterator p = ElemData[et].ElemContainer.begin();
			for (unsigned i = 0; i < unsigned(iNumElems); i++, p++) {
				VarLabels->set_cur(i);
				const long l = p->second->GetLabel();
				VarLabels->put(&l, 1);
			}
		}
	}
#endif // USE_NETCDF

	Elem* pTmpEl = 0;
	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->OutputPrepare(OH);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

void
DataManager::ElemOutput(OutputHandler& OH) const
{
	Elem* pTmpEl = NULL;

	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->Output(OH);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

void
DataManager::ElemOutput( OutputHandler& OH, const VectorHandler& X,
		const VectorHandler& XP) const
{
	Elem* pTmpEl = NULL;

	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->Output(OH, X, XP);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

/* cerca un elemento qualsiasi */
Elem *
DataManager::pFindElem(Elem::Type Typ, unsigned int uL) const
{
	if (ElemData[Typ].bIsUnique() && uL == (unsigned)(-1)) {
		if (!ElemData[Typ].ElemMapToList.empty()) {
			return ElemData[Typ].ElemMapToList.begin()->second->second;
		}
		return 0;
	}

	ElemMapToListType::const_iterator p = ElemData[Typ].ElemMapToList.find(uL);
	if (p == ElemData[Typ].ElemMapToList.end()) {
		return 0;
	}

	return p->second->second;
}


/* cerca un elemento qualsiasi */
Elem **
DataManager::ppFindElem(Elem::Type Typ, unsigned int uL) const
{
	if (ElemData[Typ].bIsUnique() && uL == (unsigned)(-1)) {
		if (!ElemData[Typ].ElemMapToList.empty()) {
			return &ElemData[Typ].ElemMapToList.begin()->second->second;
		}
		return 0;
	}

	ASSERT(uL > 0);

	ElemMapToListType::const_iterator p = ElemData[Typ].ElemMapToList.find(uL);
	if (p == ElemData[Typ].ElemMapToList.end()) {
		return 0;
	}

	return const_cast<Elem **>(&p->second->second);
}

/* cerca un elemento qualsiasi */
Elem *
DataManager::pFindElem(Elem::Type Typ, unsigned int uL,
		unsigned int iDeriv) const
{
	ASSERT(iDeriv == int(ELEM) || ElemData[Typ].iDerivation & iDeriv);

	if (ElemData[Typ].bIsUnique() && uL == (unsigned)(-1)) {
		if (!ElemData[Typ].ElemMapToList.empty()) {
			return pChooseElem(ElemData[Typ].ElemMapToList.begin()->second->second, iDeriv);
		}
		return 0;
	}

	ASSERT(uL > 0);

	ElemMapToListType::const_iterator p = ElemData[Typ].ElemMapToList.find(uL);
	if (p == ElemData[Typ].ElemMapToList.end()) {
		return 0;
	}

	return pChooseElem(p->second->second, iDeriv);
}

/* Usata dalle due funzioni precedenti */
Elem *
DataManager::pChooseElem(Elem* p, unsigned int iDeriv) const
{
	ASSERT(p != NULL);

	switch (iDeriv) {
	case ELEM:
		return p;

	case DOFOWNER:
		ASSERT(dynamic_cast<ElemWithDofs *>(p) != NULL);
		return p;

	case GRAVITYOWNER:
		ASSERT(dynamic_cast<ElemGravityOwner *>(p) != NULL);
		return p;

	case AIRPROPOWNER:
		ASSERT(dynamic_cast<AerodynamicElem *>(p) != NULL);
		return p;

	case INITIALASSEMBLY:
		ASSERT(dynamic_cast<InitialAssemblyElem *>(p) != NULL);
		return p;
	}

	/* default */
	return NULL;
}

/* cerca un drive qualsiasi */
Drive *
DataManager::pFindDrive(Drive::Type Typ, unsigned int uL) const
{
	ASSERT(uL > 0);

	if (DriveData[Typ].iNum == 0) {
		silent_cerr("FileDrive(" << uL << "): "
			"no file drivers defined" << std::endl);
		return 0;
	}

	if (DriveData[Typ].ppFirstDrive == 0) {
		silent_cerr("FileDrive(" << uL << "): "
			"file drivers can only be dereferenced "
			"after the \"drivers\" block" << std::endl);
		return 0;
	}

	return pLabelSearch(DriveData[Typ].ppFirstDrive, DriveData[Typ].iNum, uL);
}


flag
DataManager::fGetDefaultOutputFlag(const Elem::Type& t) const
{
	return ElemData[t].bDefaultOut();
}

/* DataManager - end */


/* InitialAssemblyIterator - begin */

InitialAssemblyIterator::InitialAssemblyIterator(
	const DataManager::ElemDataStructure (*pED)[Elem::LASTELEMTYPE])
: pElemData(pED),
m_FirstType(Elem::UNKNOWN),
m_CurrType(Elem::UNKNOWN)
{
	int iCnt = 0;

	while (!(*pElemData)[iCnt].bToBeUsedInAssembly()
		|| (*pElemData)[iCnt].ElemContainer.empty())
	{
		if (++iCnt >= Elem::LASTELEMTYPE) {
			break;
		}
	}

	// NOTE: this fails if by chance there's no element
	// that participates in initial assembly
	ASSERT(iCnt < Elem::LASTELEMTYPE);
	ASSERT((*pElemData)[iCnt].ElemContainer.begin() != (*pElemData)[iCnt].ElemContainer.end());

	m_FirstType = m_CurrType = Elem::Type(iCnt);
	m_CurrElem = (*pElemData)[m_FirstType].ElemContainer.begin();
}

InitialAssemblyElem *
InitialAssemblyIterator::GetFirst(void) const
{
	// FIXME: this should not happen: if there are no elements
	// that need initial assembly, we shouldn't get here in the
	// first place
	if (m_FirstType == Elem::LASTELEMTYPE) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_CurrType = m_FirstType;
	m_CurrElem = (*pElemData)[m_FirstType].ElemContainer.begin();

	/* La variabile temporanea e' necessaria per il debug. */
	InitialAssemblyElem* p = dynamic_cast<InitialAssemblyElem *>(m_CurrElem->second);
	if (p == 0) {
		DrivenElem *pDE = dynamic_cast<DrivenElem *>(m_CurrElem->second);
		if (pDE != 0) {
			p = dynamic_cast<InitialAssemblyElem *>(pDE->pGetElem());
		}
		if (p == 0) {
			p = GetNext();
		}
	}

	return p;
}

InitialAssemblyElem*
InitialAssemblyIterator::GetNext(void) const
{
	InitialAssemblyElem* p = 0;
	do {
		++m_CurrElem;
		if (m_CurrElem == (*pElemData)[m_CurrType].ElemContainer.end()) {
			int iCnt = int(m_CurrType);

			do {
				if (++iCnt >= Elem::LASTELEMTYPE) {
					return 0;
				}
			} while (!(*pElemData)[iCnt].bToBeUsedInAssembly()
					|| (*pElemData)[iCnt].ElemContainer.empty());

			ASSERT((*pElemData)[iCnt].ElemContainer.begin() != (*pElemData)[iCnt].ElemContainer.end());
			m_CurrType = Elem::Type(iCnt);
			m_CurrElem = (*pElemData)[iCnt].ElemContainer.begin();
		}

		/* La variabile temporanea e' necessaria per il debug. */
		p = dynamic_cast<InitialAssemblyElem *>(m_CurrElem->second);

		if (p == 0) {
			DrivenElem *pDE = dynamic_cast<DrivenElem *>(m_CurrElem->second);
			if (pDE != 0) {
				p = dynamic_cast<InitialAssemblyElem *>(pDE->pGetElem());
			}
		}

#ifdef DEBUG
		if (p == 0) {
			silent_cerr(psElemNames[m_CurrElem->second->GetElemType()]
				<< "(" << m_CurrElem->second->GetLabel() << ")"
				" is not subjected to initial assembly"
				<< std::endl);
		}
#endif
	} while (p == 0);

	return p;
}

/* InitialAssemblyIterator - end */

