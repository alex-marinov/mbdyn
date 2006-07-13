/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "search.h"

/* DataManager - begin */

/* costruttore: resetta i dati */
void
DataManager::ElemManager(void)
{
	/* Reset della struttura ElemData */
	for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
		ElemData[i].ppFirstElem = NULL;
		ElemData[i].iNum = 0;
		ElemData[i].DofOwnerType = DofOwner::UNKNOWN;
		ElemData[i].uFlags = 0U;
		ElemData[i].DefaultOut(::fDefaultOut == 1); /* Da "output.h" */
		ElemData[i].OutFile = OutputHandler::UNKNOWN;	/* "output.h" */
	}

	/* Add dof type */
	ElemData[Elem::JOINT].DofOwnerType = DofOwner::JOINT;
	ElemData[Elem::ELECTRIC].DofOwnerType = DofOwner::ELECTRIC;
	ElemData[Elem::ELECTRICBULK].DofOwnerType = DofOwner::ELECTRICBULK;
	ElemData[Elem::GENEL].DofOwnerType = DofOwner::GENEL;
	ElemData[Elem::ROTOR].DofOwnerType = DofOwner::ROTOR;
	ElemData[Elem::AEROMODAL].DofOwnerType = DofOwner::AEROMODAL;
	ElemData[Elem::HYDRAULIC].DofOwnerType = DofOwner::HYDRAULIC;
	ElemData[Elem::LOADABLE].DofOwnerType = DofOwner::LOADABLE;

	/* Add file type */
	ElemData[Elem::AUTOMATICSTRUCTURAL].OutFile = OutputHandler::INERTIA;
	ElemData[Elem::JOINT].OutFile = OutputHandler::JOINTS;
	ElemData[Elem::FORCE].OutFile = OutputHandler::FORCES;
	ElemData[Elem::BEAM].OutFile = OutputHandler::BEAMS;
	ElemData[Elem::ROTOR].OutFile = OutputHandler::ROTORS;
	ElemData[Elem::AEROMODAL].OutFile = OutputHandler::AEROMODALS;
	ElemData[Elem::AERODYNAMIC].OutFile = OutputHandler::AERODYNAMIC;
	ElemData[Elem::HYDRAULIC].OutFile = OutputHandler::HYDRAULIC;
	ElemData[Elem::LOADABLE].OutFile = OutputHandler::LOADABLE;
	ElemData[Elem::GENEL].OutFile = OutputHandler::GENELS;
	ElemData[Elem::EXTERNAL].OutFile = OutputHandler::EXTERNALS;

	/* Inheritance scheme */
	ElemData[Elem::AUTOMATICSTRUCTURAL].iDerivation = ELEM;
	ElemData[Elem::GRAVITY].iDerivation = ELEM;
	ElemData[Elem::BODY].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
	ElemData[Elem::JOINT].iDerivation = ELEM | GRAVITYOWNER | DOFOWNER | INITIALASSEMBLY;
	ElemData[Elem::GENEL].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::FORCE].iDerivation = ELEM | INITIALASSEMBLY;
	ElemData[Elem::BEAM].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
	ElemData[Elem::PLATE].iDerivation = ELEM | INITIALASSEMBLY;
	ElemData[Elem::AIRPROPERTIES].iDerivation = ELEM | INITIALASSEMBLY;
	ElemData[Elem::ROTOR].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
	ElemData[Elem::AEROMODAL].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
	ElemData[Elem::AERODYNAMIC].iDerivation = ELEM | AIRPROPOWNER |INITIALASSEMBLY;
	ElemData[Elem::ELECTRIC].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::ELECTRICBULK].iDerivation = ELEM | DOFOWNER;
	ElemData[Elem::BULK].iDerivation = ELEM;
	ElemData[Elem::LOADABLE].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY | DOFOWNER;
	ElemData[Elem::DRIVEN].iDerivation = ELEM;
	ElemData[Elem::AEROMODAL].iDerivation = ELEM | AIRPROPOWNER;

	/* Types that must be unique */
	ElemData[Elem::GRAVITY].IsUnique(true);
	ElemData[Elem::AIRPROPERTIES].IsUnique(true);

	/* Types that are used by default in initial assembly */
	ElemData[Elem::JOINT].ToBeUsedInAssembly(true);
	ElemData[Elem::BEAM].ToBeUsedInAssembly(true);

	/* Aggiungere qui se un tipo genera forze d'inerzia e quindi
	 * deve essere collegato all'elemento accelerazione di gravita' */
	ElemData[Elem::BODY].GeneratesInertiaForces(true);
	ElemData[Elem::JOINT].GeneratesInertiaForces(true);
	ElemData[Elem::LOADABLE].GeneratesInertiaForces(true);

	/* Aggiungere qui se un tipo usa le proprieta' dell'aria e quindi
	 * deve essere collegato all'elemento proprieta' dell'aria */
	ElemData[Elem::ROTOR].UsesAirProperties(true);
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

	/* Distruzione elementi */
	ASSERT(ppElems != NULL);

	if (ppElems != NULL) {
		Elem* p = NULL;
		if (ElemIter.bGetFirst(p)) {
			do {
				ASSERT(p != NULL);
				if (p != NULL) {
					DEBUGCOUT("deleting element "
							<< p->GetLabel()
							<< ", type "
							<< psElemNames[p->GetElemType()]
							<< std::endl);
					SAFEDELETE(p);
				}
			} while (ElemIter.bGetNext(p));
		}

		DEBUGCOUT("deleting elements structure" << std::endl);
		SAFEDELETEARR(ppElems);
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
	/* struttura degli elementi */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		iTotElem += ElemData[iCnt].iNum;
	}

	DEBUGCOUT("iTotElem = " << iTotElem << std::endl);

	if (iTotElem > 0) {
		SAFENEWARR(ppElems, Elem*, iTotElem);

		/* Inizializza l'iteratore degli elementi usato all'interno
		 * dell'ElemManager */
		ElemIter.Init(ppElems, iTotElem);

		Elem** ppTmp = ppElems;
		while (ppTmp < ppElems + iTotElem) {
			*ppTmp++ = NULL;
		}

		ElemData[0].ppFirstElem = ppElems;
		for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE-1; iCnt++) {
			ElemData[iCnt+1].ppFirstElem =
				ElemData[iCnt].ppFirstElem + ElemData[iCnt].iNum;
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
	ASSERT(ppElems != NULL);

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
			JacHdl(r, c) = 1e-16;
		}
	}
#endif

	Elem* pTmpEl = NULL;
	if (Iter.bGetFirst(pTmpEl)) {
		do {
			JacHdl += pTmpEl->AssJac(WorkMat, dCoef,
					*pXCurr, *pXPrimeCurr);

		} while (Iter.bGetNext(pTmpEl));
	}
}

void
DataManager::AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl)
{
	DEBUGCOUT("Entering DataManager::AssMats()" << std::endl);

	ASSERT(pWorkMatA != NULL);
	ASSERT(pWorkMatB != NULL);
	ASSERT(ppElems != NULL);

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

		} while (Iter.bGetNext(pTmpEl));
	}
	if (ChangedEqStructure) {
		throw ChangedEquationStructure();
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

#if 0
void
DataManager::ElemOutput_pch( std::ostream& pch) const
{
	Elem* pTmpEl = NULL;

	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->Output_pch(pch);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

void
DataManager::ElemOutput_f06( std::ostream& f06, const VectorHandler& X) const
{
	Elem* pTmpEl = NULL;

	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->Output_f06(f06, X);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

void
DataManager::ElemOutput_f06( std::ostream& f06, const VectorHandler& Xr,
		const VectorHandler& Xi) const
{
	Elem* pTmpEl = NULL;

	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->Output_f06(f06, Xr, Xi);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}
#endif


/* cerca un elemento qualsiasi */
void *
DataManager::pFindElem(Elem::Type Typ, unsigned int uL) const
{
	ASSERT(ElemData[Typ].ppFirstElem != NULL);
	ASSERT(ElemData[Typ].iNum > 0);
	ASSERT(uL > 0);

	Elem* p = pLabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);

	if (p == NULL) {
		return NULL;
	}

	ASSERT(p->pGetElem() != NULL);
	return p->pGetElem();
}


/* cerca un elemento qualsiasi */
Elem **
DataManager::ppFindElem(Elem::Type Typ, unsigned int uL) const
{
	ASSERT(ElemData[Typ].ppFirstElem != NULL);
	ASSERT(ElemData[Typ].iNum > 0);
	ASSERT(uL > 0);

	int i = LabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);

	if (i < 0) {
		return NULL;
	}

	return ElemData[Typ].ppFirstElem+i;
}

/* cerca un elemento qualsiasi */
void *
DataManager::pFindElem(Elem::Type Typ, unsigned int uL,
		unsigned int iDeriv) const
{
	ASSERT(ElemData[Typ].ppFirstElem != NULL);
	ASSERT(ElemData[Typ].iNum > 0);
	ASSERT(uL > 0);
	ASSERT(iDeriv == int(ELEM) || ElemData[Typ].iDerivation & iDeriv);

	Elem* p = pLabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);

	if (p == NULL) {
		return NULL;
	}

	return pChooseElem(p, iDeriv);
}

/* Usata dalle due funzioni precedenti */
void *
DataManager::pChooseElem(Elem* p, unsigned int iDeriv) const
{
	ASSERT(p != NULL);

	switch (iDeriv) {
	case ELEM:
		ASSERT(p->pGetElem() != NULL);
		return p->pGetElem();

	case DOFOWNER:
		ASSERT(p->pGetElemWithDofs() != NULL);
		return p->pGetElemWithDofs();

	case GRAVITYOWNER:
		ASSERT(p->pGetElemGravityOwner() != NULL);
		return p->pGetElemGravityOwner();

	case AIRPROPOWNER:
		ASSERT(p->pGetAerodynamicElem() != NULL);
		return p->pGetAerodynamicElem();

	case INITIALASSEMBLY:
		ASSERT(p->pGetInitialAssemblyElem() != NULL);
		return p->pGetInitialAssemblyElem();
	}

	/* default */
	return NULL;
}

/* cerca un drive qualsiasi */
void *
DataManager::pFindDrive(Drive::Type Typ, unsigned int uL) const
{
	ASSERT(DriveData[Typ].ppFirstDrive != NULL);
	ASSERT(DriveData[Typ].iNum > 0);
	ASSERT(uL > 0);

	Drive* p = pLabelSearch(DriveData[Typ].ppFirstDrive, DriveData[Typ].iNum, uL);

	if (p == NULL) {
		return NULL;
	}

	return p;
}


flag
DataManager::fGetDefaultOutputFlag(const Elem::Type& t) const
{
	return ElemData[t].bDefaultOut();
}

/* DataManager - end */


/* InitialAssemblyIterator - begin */

InitialAssemblyIterator::InitialAssemblyIterator(
		const DataManager::ElemDataStructure (*pED)[Elem::LASTELEMTYPE]
		)
: pElemData(pED),
FirstType(Elem::UNKNOWN), ppFirst(NULL),
CurrType(Elem::UNKNOWN), ppCurr(NULL)
{
	int iCnt = 0;

	while (!(*pElemData)[iCnt].bToBeUsedInAssembly()
			|| (*pElemData)[iCnt].iNum == 0)
	{
		if (++iCnt >= Elem::LASTELEMTYPE) {
			break;
		}
	}

	ASSERT(iCnt < Elem::LASTELEMTYPE);
	ASSERT((*pElemData)[iCnt].ppFirstElem != NULL);

	(Elem**&)ppFirst = (*pElemData)[iCnt].ppFirstElem;
	ppCurr = (Elem**)ppFirst;
	(Elem::Type&)FirstType = CurrType = Elem::Type(iCnt);
}

InitialAssemblyElem *
InitialAssemblyIterator::GetFirst(void) const
{
	CurrType = FirstType;
	ppCurr = (Elem**)ppFirst;

	/* La variabile temporanea e' necessaria per il debug. */
	InitialAssemblyElem* p;
	for (p = (*ppCurr)->pGetInitialAssemblyElem();
			p == 0;
			p = GetNext())
	{
#ifdef DEBUG
		if (p == NULL) {
			silent_cerr(psElemNames[(*ppCurr)->GetElemType()]
				<< "(" << (*ppCurr)->GetLabel() << ")"
				" is not subject to initial assembly" << std::endl);
		}
#endif
	}

	return p;
}

InitialAssemblyElem* InitialAssemblyIterator::GetNext(void) const
{
	InitialAssemblyElem* p = 0;
	do {
		ppCurr++;
		if (ppCurr >= (*pElemData)[CurrType].ppFirstElem
				+ (*pElemData)[CurrType].iNum)
		{
			int iCnt = int(CurrType);

			do {
				if (++iCnt >= Elem::LASTELEMTYPE) {
					return NULL;
				}
			} while (!(*pElemData)[iCnt].bToBeUsedInAssembly()
					|| (*pElemData)[iCnt].iNum == 0);

			ASSERT((*pElemData)[iCnt].ppFirstElem != NULL);
			CurrType = Elem::Type(iCnt);
			ppCurr = (*pElemData)[iCnt].ppFirstElem;
		}

		/* La variabile temporanea e' necessaria per il debug. */
		p = (*ppCurr)->pGetInitialAssemblyElem();

#ifdef DEBUG
		if (p == 0) {
			silent_cerr(psElemNames[(*ppCurr)->GetELemType()]
				<< "(" << (*ppCurr)->GetLabel() << ")"
				" is not subjected to initial assembly"
				<< std::endl);
		}
#endif
	} while (p == 0);

	return p;
}

/* InitialAssemblyIterator - end */

