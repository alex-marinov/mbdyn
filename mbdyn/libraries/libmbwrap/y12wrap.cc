/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

/*****************************************************************************
 *                                                                           *
 *                            Y12 C++ WRAPPER                                *
 *                                                                           *
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_Y12

#include <y12wrap.h>
#include <y12lib.h>

/* Y12LUSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
Y12LUSolver::Y12LUSolver(integer iMatOrd, integer iSize,
			 integer** ppiTmpRow, integer** ppiTmpCol,
			 doublereal** ppdTmpMat,
			 doublereal* pdTmpRhs, integer iPivotParam)
: iMatSize(iSize),
ppiRow(ppiTmpRow),
ppiCol(ppiTmpCol),
ppdMat(ppdTmpMat),
iN(iMatOrd),
iNonZeroes(0),
pdRhs(pdTmpRhs),
piHA(NULL),
pdPIVOT(NULL),
iFirstSol(-1)
{
	ASSERT(iMatSize > 0);
	ASSERT(ppiRow != NULL);
	ASSERT(ppiCol != NULL);
	ASSERT(ppdMat != NULL);
	ASSERT(*ppiRow != NULL);
	ASSERT(*ppiCol != NULL);
	ASSERT(*ppdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsValid(*ppiRow, 
				iMatSize*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsValid(*ppiCol, 
				iMatSize*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsValid(*ppdMat, 
				iMatSize*sizeof(doublereal)));
	ASSERT(defaultMemoryManager.fIsValid(pdRhs, 
				iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

	SAFENEWARR(piHA, integer, 11*iN);
	SAFENEWARR(pdPIVOT, doublereal, iN);
	
#ifdef DEBUG
	for (int iCnt = 0; iCnt < 11*iN; iCnt++) {
		piHA[iCnt] = 0;
	}
	for (int iCnt = 0; iCnt < iN; iCnt++) {
		pdPIVOT[iCnt] = 0;
	}
#endif /* DEBUG */

        iIFLAG[0] = 0;
	iIFLAG[1] = 3;
	iIFLAG[2] = iPivotParam;
	iIFLAG[3] = 0;
	iIFLAG[4] = 2;
					
        dAFLAG[0] = 8.;         /* Should be 4.<dAFLAG[0]<16. for stability */
	dAFLAG[1] = 0.;         /* Should be 0.<dAFLAG[1]<1.e-12 */
	dAFLAG[2] = 1.e6;       /* Should be dAFLAG[2]>1.e5 */
	dAFLAG[3] = 0.;  	/* FIXME: Should be 0<dAFLAG[3]<1.e-12 */
}

/* Distruttore */
Y12LUSolver::~Y12LUSolver(void)
{
	if (pdPIVOT != NULL) {
		SAFEDELETEARR(pdPIVOT);
	}
	if (piHA != NULL) {
		SAFEDELETEARR(piHA);
	}
}

void 
Y12LUSolver::IsValid(void) const
{
	ASSERT(iMatSize > 0);
	ASSERT(iMatSize > 0);
	ASSERT(ppiRow != NULL);
	ASSERT(ppiCol != NULL);
	ASSERT(ppdMat != NULL);
	ASSERT(*ppiRow != NULL);
	ASSERT(*ppiCol != NULL);
	ASSERT(*ppdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0); 
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsValid(*ppiRow, 
				iMatSize*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsValid(*ppiCol, 
				iMatSize*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsValid(*ppdMat, 
				iMatSize*sizeof(doublereal)));
	
	ASSERT(defaultMemoryManager.fIsValid(pdRhs, 
				iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

	ASSERT(piHA != NULL);
	ASSERT(pdPIVOT != NULL);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsBlock(piHA, 11*iN*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsBlock(pdPIVOT, 1*iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}

/* Fattorizza la matrice */
flag
Y12LUSolver::fLUFactor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	/*
	 * FIXME: This is set by Y12SparseLUSolutionManager in PacVec;
	 * better move such info to the matrix handler!
	 */
	ASSERT(iNonZeroes > 0);
	
	/* Sets parameters */
	integer iIFAIL = 0;
	iIFLAG[4] = 2;
	iFirstSol = 1;
	
	/* Prepares for factorization */
	__FC_DECL__(y12mbf)(&iN, &iNonZeroes, *ppdMat,
			    *ppiCol, &iMatSize,
			    *ppiRow, &iMatSize,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);
			    
	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during pre-factorization, code "
			<< iIFAIL << ":" << endl;
		PutError(cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	/* actual factorization */
	__FC_DECL__(y12mcf)(&iN, &iNonZeroes, *ppdMat,
			    *ppiCol, &iMatSize,
			    *ppiRow, &iMatSize,
			    pdPIVOT, pdRhs,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);

	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during factorization, code "
			<< iIFAIL << ":" << endl;
		PutError(cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}

	if (dAFLAG[7] < 1.e-12) {
		cerr << "Y12LUSolver:"
			" warning, possible bad conditioning of matrix" << endl;
	}
	
	return iIFAIL;
}

/* Risolve */
void
Y12LUSolver::Solve(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	integer iIFAIL = 0;
	
	__FC_DECL__(y12mdf)(&iN, *ppdMat, &iMatSize, pdRhs,
			    pdPIVOT, *ppiCol,
			    piHA, &iN,
			    iIFLAG, &iIFAIL);
	
	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during back substitution, code "
			<< iIFAIL << ":" << endl;
		PutError(cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	if (iFirstSol == 1) {
		iIFLAG[4] = 3;
		iFirstSol = 0;
	}
}

void 
Y12LUSolver::PutError(ostream& out, int rc) const
{
	out << endl;

	switch (rc ) {
	case 1:
		out 
			<< "\tThe    coefficient   matrix   A   is   not" << endl
			<< "\tfactorized, i.e. the  call  of  subroutine" << endl
			<< "\tY12MD  was not preceded by a call of Y12MC" << endl
			<< "\tduring the solution of   Ax=b   or  during" << endl
			<< "\tthe  solution  of  the  first  system in a" << endl
			<< "\tsequence ( Ax1 = b1 , Ax2 = b2,.....,Axp =" << endl
			<< "\tbp)  of  systems with the same coefficient" << endl
			<< "\tmatrix. This will work in all  cases  only" << endl
			<< "\tif  the  user  sets IFLAG(1) .ge. 0 before" << endl
			<< "\tthe call of package Y12M (i.e. before  the" << endl
			<< "\tfirst   call   of  a  subroutine  of  this" << endl
			<< "\tpackage)." << endl;
		break;

	case 2:
		out
			<< "\tThe coefficient matrix A is  not  ordered," << endl
			<< "\ti.e.  the call of subroutine Y12MC was not" << endl
			<< "\tpreceded by a call  of  Y12MB.  This  will" << endl
			<< "\twork  in  all  cases only if the user sets" << endl
			<< "\tIFLAG(1) .ge. 0 before the call of package" << endl
			<< "\tY12M  (i.e.  before  the  first  call of a" << endl
			<< "\tsubroutine of this package)." << endl;
		break;

	case 3:
		out
			<< "\tA pivotal element abs(a(i,i;j)) < AFLAG(4)" << endl
			<< "\t*  AFLAG(6) is selected.  When AFLAG(4) is" << endl
			<< "\tsufficiently small this is  an  indication" << endl
			<< "\tthat the coefficient matrix is numerically" << endl
			<< "\tsingular." << endl;
		break;
		
	case 4:
		out
			<< "\tAFLAG(5), the  growth  factor,  is  larger" << endl
			<< "\tthan    AFLAG(3).    When    AFLAG(3)   is" << endl
			<< "\tsufficiently large this indicates that the" << endl
			<< "\telements  of the coefficient matrix A grow" << endl
			<< "\tso quickly during the  factorization  that" << endl
			<< "\tthe continuation of the computation is not" << endl
			<< "\tjustified.  The  choice   of   a   smaller" << endl
			<< "\tstability   factor,   AFLAG(1),  may  give" << endl
			<< "\tbetter results in this case." << endl;
		break;

	case 5:
		out
			<< "\tThe length NN of arrays A and SNR  is  not" << endl
			<< "\tsufficient.   Larger  values  of  NN  (and" << endl
			<< "\tpossibly of NN1) should be used." << endl;
		break;

	case 6:
		out
			<< "\tThe  length  NN1  of  array  RNR  is   not" << endl
			<< "\tsufficient.   Larger  values  of  NN1 (and" << endl
			<< "\tpossibly of NN) should be used." << endl;
		break;

	case 7:
		out
			<< "\tA row without  non-zero  elements  in  its" << endl
			<< "\tactive    part   is   found   during   the" << endl
			<< "\tdecomposition.  If   the   drop-tolerance," << endl
			<< "\tAFLAG(2),   is  sufficiently  small,  then" << endl
			<< "\tIFAIL = 7 indicates  that  the  matrix  is" << endl
			<< "\tnumerically  singular. If a large value of" << endl
			<< "\tthe drop-tolerance AFLAG(2) is used and if" << endl
			<< "\tIFAIL = 7  on exit, this is not certain. A" << endl
			<< "\trun  with  a  smaller  value  of  AFLAG(2)" << endl
			<< "\tand/or  a  careful check of the parameters" << endl
			<< "\tAFLAG(8) and AFLAG(5)  is  recommended  in" << endl
			<< "\tthe latter case." << endl;
		break;

	case 8:
		out
			<< "\tA  column without non-zero elements in its" << endl
			<< "\tactive   part   is   found   during    the" << endl
			<< "\tdecomposition.   If   the  drop-tolerance," << endl
			<< "\tAFLAG(2),  is  sufficiently  small,   then" << endl
			<< "\tIFAIL  =  8  indicates  that the matrix is" << endl
			<< "\tnumerically singular. If a large value  of" << endl
			<< "\tthe drop-tolerance AFLAG(2) is used and if" << endl
			<< "\tIFAIL = 8  on exit, this is not certain. A" << endl
			<< "\trun  with  a  smaller  value  of  AFLAG(2)" << endl
			<< "\tand/or a careful check of  the  parameters" << endl
			<< "\tAFLAG(8)  and  AFLAG(5)  is recommended in" << endl
			<< "\tthe latter case." << endl;
		break;
	
	case 9:
		out
			<< "\tA pivotal element  is  missing.  This  may" << endl
			<< "\toccur  if  AFLAG(2)  >  0 and IFLAG(4) = 2" << endl
			<< "\t(i.e. some system after the first one in a" << endl
			<< "\tsequence   of   systems   with   the  same" << endl
			<< "\tstructure is solved using a positive value" << endl
			<< "\tfor  the drop-tolerance). The value of the" << endl
			<< "\tdrop-tolerance   AFLAG(2),    should    be" << endl
			<< "\tdecreased  and  the  coefficient matrix of" << endl
			<< "\tthe system refactorized.  This  error  may" << endl
			<< "\talso occur when one of the special pivotal" << endl
			<< "\tstrategies (IFLAG(3)=0 or  IFLAG(3)=2)  is" << endl
			<< "\tused  and  the  matrix is not suitable for" << endl
			<< "\tsuch a strategy." << endl;
		break;

	case 10:
		out
			<< "\tSubroutine Y12MF is called with IFLAG(5) =" << endl
			<< "\t1  (i.e.  with  a  request  to  remove the" << endl
			<< "\tnon-zero elements of the lower  triangular" << endl
			<< "\tmatrix    L).     IFLAG(5)=2     must   be" << endl
			<< "\tinitialized instead of IFLAG(5)=1." << endl;
		break;

	case 11:
		out
			<< "\tThe coefficient matrix A contains at least" << endl
			<< "\ttwo  elements  in the same position (i,j)." << endl
			<< "\tThe  input   data   should   be   examined" << endl
			<< "\tcarefully in this case." << endl;
		break;

	case 12:
		out
			<< "\tThe number of equations in the system Ax=b" << endl
			<< "\tis smaller than 2 (i.e.  N<2).  The  value" << endl
			<< "\tof N should be checked." << endl;
		break;
		
	case 13:
		out
			<< "\tThe  number  of  non-zero  elements of the" << endl
			<< "\tcoefficient matrix is  non-positive  (i.e." << endl
			<< "\tZ.le.0  ).   The  value of the parameter Z" << endl
			<< "\t(renamed NZ in Y12MF) should be checked." << endl;
		break;

	case 14:
		out
			<< "\tThe number of  non-zero  elements  in  the" << endl
			<< "\tcoefficient  matrix  is  smaller  than the" << endl
			<< "\tnumber of equations (i.e. Z  <  N  ).   If" << endl
			<< "\tthere  is no mistake (i.e. if parameter Z," << endl
			<< "\trenamed NZ in Y12MF, is correctly assigned" << endl
			<< "\ton  entry)  then the coefficient matrix is" << endl
			<< "\tstructurally singular in this case." << endl;
		break;

	case 15:
		out
			<< "\tThe length IHA of the first  dimension  of" << endl
			<< "\tarray  HA  is  smaller  than  N.  IHA.ge.N" << endl
			<< "\tshould be assigned." << endl;
		break;

	case 16:
		out
			<< "\tThe value of  parameter  IFLAG(4)  is  not" << endl
			<< "\tassigned  correctly.  IFLAG(4)  should  be" << endl
			<< "\tequal to 0, 1 or 2. See  more  details  in" << endl
			<< "\tthe description of this parameter." << endl;
		break;
		
	case 17:
		out
			<< "\tA  row  without non-zero elements has been" << endl
			<< "\tfound in the coefficient matrix A  of  the" << endl
			<< "\tsystem  before the Gaussian elimination is" << endl
			<< "\tinitiated.  Matrix   A   is   structurally" << endl
			<< "\tsingular." << endl;
		break;

	case 18:
		out
			<< "\tA  column  without  non-zero  elements has" << endl
			<< "\tbeen found in the coefficient matrix A  of" << endl
			<< "\tthe system before the Gaussian elimination" << endl
			<< "\tis initiated.  Matrix  A  is  structurally" << endl
			<< "\tsingular." << endl;
		break;

	case 19:
		out
			<< "\tParameter  IFLAG(2) is smaller than 1. The" << endl
			<< "\tvalue of IFLAG(2)  should  be  a  positive" << endl
			<< "\tinteger (IFLAG(2) = 3 is recommended)." << endl;
		break;

	case 20:
		out
			<< "\tParameter   IFLAG(3)   is  out  of  range." << endl
			<< "\tIFLAG(3) should be equal to 0, 1 or 2." << endl;
		break;

	case 21:
		out
			<< "\tParameter  IFLAG(5)  is  out   of   range." << endl
			<< "\tIFLAG(5) should be equal to 1, 2 or 3 (but" << endl
			<< "\twhen IFLAG(5) = 3 Y12MB and  Y12MC  should" << endl
			<< "\tnot  be  called;  see also the message for" << endl
			<< "\tIFAIL = 22 below)." << endl;
		break;

	case 22:
		out
			<< "\tEither  subroutine  Y12MB  or   subroutine" << endl
			<< "\tY12MC is called with IFLAG(5) = 3. Each of" << endl
			<< "\tthese subroutines should  be  called  with" << endl
			<< "\tIFLAG(5) equal to 1 or 2." << endl;
		break;

	case 23:
		out
			<< "\tThe    number    of   allowed   iterations" << endl
			<< "\t(parameter IFLAG(11) when Y12MF  is  used)" << endl
			<< "\tis  smaller  than  2.   IFLAG(11)  .ge.  2" << endl
			<< "\tshould be assigned." << endl;
		break;

	case 24:
		out
			<< "\tAt least one element whose  column  number" << endl
			<< "\tis  either larger than N or smaller than 1" << endl
			<< "\tis found." << endl;
		break;

	case 25:
		out
			<< "\tAt least one element whose row  number  is" << endl
			<< "\teither  larger than N or smaller than 1 is" << endl
			<< "\tfound." << endl;
		break;

	default:
		out
			<<"\t Unhandled code." << endl;
		break;
	}

	out << endl;
}

/* Y12LUSolver - end */


/* Y12SparseLUSolutionManager - begin: code */

/* Costruttore */
Y12SparseLUSolutionManager::Y12SparseLUSolutionManager(integer iSize, 
						       integer iWorkSpaceSize,
						       const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
piRow(NULL),
piCol(NULL), 
pdMat(NULL),
pdVec(NULL),
pMH(NULL),
pVH(NULL),
pLU(NULL),
fHasBeenReset(1)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
      		iWorkSpaceSize = 2*iSize*iSize;
   	}

	integer iPivot;
	if (dPivotFactor == 0.) {
		iPivot = 0;
	} else {
		iPivot = 1;
	}
      
   	/* Alloca arrays */
   	SAFENEWARR(piRow, integer, iWorkSpaceSize);
   	SAFENEWARR(piCol, integer, iWorkSpaceSize);
   	SAFENEWARR(pdMat, doublereal, iWorkSpaceSize);
   	SAFENEWARR(pdVec, doublereal, iMatSize);
   
   	/* Alloca handlers ecc. */
   	SAFENEWWITHCONSTRUCTOR(pMH, 
			       SparseMatrixHandler,
			       SparseMatrixHandler(iMatSize, &piRow, 
			       			   &piCol, &pdMat,
			       			   iWorkSpaceSize));
   	SAFENEWWITHCONSTRUCTOR(pVH,
			       MyVectorHandler,
			       MyVectorHandler(iMatSize, pdVec));
   	SAFENEWWITHCONSTRUCTOR(pLU, 
			       Y12LUSolver,
			       Y12LUSolver(iMatSize, iWorkSpaceSize,
			       		   &piRow, &piCol,
					   &pdMat, pdVec, iPivot));
   
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
Y12SparseLUSolutionManager::~Y12SparseLUSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	/* Dealloca oggetti strani */
   	if (pLU != NULL) {	
      		SAFEDELETE(pLU);
   	}
   	if (pVH != NULL) {      
      		SAFEDELETE(pVH);
   	}
   	if (pMH != NULL) {
      		SAFEDELETE(pMH);
   	}
   
   	/* Dealloca arrays */
   	if (pdVec != NULL) {	
      		SAFEDELETEARR(pdVec);
   	}
   	if (pdMat != NULL) {	
      		SAFEDELETEARR(pdMat);
   	}
   	if (piCol != NULL) {	
      		SAFEDELETEARR(piCol);
   	}
   	if (piRow != NULL) {	
      		SAFEDELETEARR(piRow);
   	}
}

/* Test di validita' del manager */
void 
Y12SparseLUSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatMaxSize > 0);
   	ASSERT(iMatSize > 0);
   	ASSERT(pMH != NULL);
   	ASSERT(pdMat != NULL);
   	ASSERT(piRow != NULL);
   	ASSERT(piCol != NULL);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(piRow));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(piCol));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pdMat));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pMH));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pVH));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLU));
#endif /* DEBUG_MEMMANAGER */
   
   	ASSERT((pMH->IsValid(), 1));
   	ASSERT((pVH->IsValid(), 1));
   	ASSERT((pLU->IsValid(), 1));
}

/* Prepara i vettori e la matrice per il solutore */
void
Y12SparseLUSolutionManager::PacVec(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	ASSERT(fHasBeenReset == 1);
	
	/* FIXME: move this to the matrix handler! */
   	pLU->iNonZeroes = pMH->iPacVec();
}

/* Inizializza il gestore delle matrici */
void
Y12SparseLUSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	pMH->Init(dResetVal);
   	fHasBeenReset = flag(1);
}

/* Risolve il problema */
void
Y12SparseLUSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	if (fHasBeenReset == 1) {
      		this->PacVec();
      		fHasBeenReset = 0;
      		if (pLU->fLUFactor() < 0) {	 
	 		THROW(Y12SparseLUSolutionManager::ErrGeneric());
      		}
   	}
   	pLU->Solve();
}

/* Y12SparseLUSolutionManager - end */

#endif /* USE_Y12 */
