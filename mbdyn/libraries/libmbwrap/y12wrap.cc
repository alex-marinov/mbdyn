/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_Y12

#include <cstring>
#include "spmh.h"
#include "spmapmh.h"
#include "dirccmh.h"
#include "ccmh.h"
#include "y12lib.h"
#include "y12wrap.h"

/* Y12Solver - begin */

/* Costruttore: si limita ad allocare la memoria */
Y12Solver::Y12Solver(integer iMatOrd, integer iWorkSpaceSize,
			 doublereal* pdTmpRhs, 
			 integer iPivotParam,
			 bool bDupInd)
: LinearSolver(0),
iMaxSize(iWorkSpaceSize),
iCurSize(iWorkSpaceSize),
piRow(0),
piCol(0),
pdMat(0),
pir(0),
pic(0),
bDuplicateIndices(bDupInd),
iN(iMatOrd),
iNonZeroes(0),
piHA(NULL),
pdPIVOT(NULL)
{
	(void)pdSetResVec(pdTmpRhs);
	(void)pdSetSolVec(pdTmpRhs);

	ASSERT(pdTmpRhs != NULL);
	ASSERT(iN > 0);

	if (bDuplicateIndices) {
		/*
		 * NOTE: Y12 alters the index arrays :(
		 */
		iRow.reserve(iCurSize);
		iCol.reserve(iCurSize);
	}
	
	SAFENEWARR(piHA, integer, 11*iN);
	SAFENEWARR(pdPIVOT, doublereal, iN);
	
#ifdef DEBUG
	for (integer iCnt = 0; iCnt < 11*iN; iCnt++) {
		piHA[iCnt] = 0;
	}

	for (integer iCnt = 0; iCnt < iN; iCnt++) {
		pdPIVOT[iCnt] = 0;
	}
#endif /* DEBUG */

        iIFLAG[I_1] = 0;
	iIFLAG[I_2] = 3;	/* recommended row number for pivoting */
	iIFLAG[I_3] = iPivotParam;
	iIFLAG[I_4] = 0;
	iIFLAG[I_5] = 2;	/* store non-zero elements of L */
					
        dAFLAG[I_1] = 8.;	/* Should be 4.<dAFLAG[0]<16. for stability */
	dAFLAG[I_2] = 0.;	/* Should be 0.<dAFLAG[1]<1.e-12 */
	dAFLAG[I_3] = 1.e6;	/* Should be dAFLAG[2]>1.e5 */
	dAFLAG[I_4] = 0.;   	/* FIXME: Should be 0 < dAFLAG[3]<1.e-12 */
}

/* Distruttore */
Y12Solver::~Y12Solver(void)
{
	if (pdPIVOT != NULL) {
		SAFEDELETEARR(pdPIVOT);
	}
	if (piHA != NULL) {
		SAFEDELETEARR(piHA);
	}
}

#ifdef DEBUG
void 
Y12Solver::IsValid(void) const
{
	ASSERT(iCurSize > 0);
	ASSERT(piRow != NULL);
	ASSERT(piCol != NULL);
	ASSERT(pdMat != NULL);
	ASSERT(iN > 0); 
	
	ASSERT(piHA != NULL);
	ASSERT(pdPIVOT != NULL);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsBlock(piHA, 11*iN*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsBlock(pdPIVOT, 1*iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */

/* Fattorizza la matrice */
void
Y12Solver::Factor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	ASSERT(iNonZeroes > 0);
	
	/* Sets parameters */
	integer iIFAIL = 0;

	/*
	 * must be set to 0 before the first call to a routine 
	 * of the y12m package
	 */
	iIFLAG[I_1] = 0;
	iIFLAG[I_5] = 2;

	if (bDuplicateIndices) {
		/*
		 * NOTE: Y12 alters the index arrays :(
		 *
		 * FIXME: make it stl-ish
		 */
		iRow.resize(iCurSize);
		iCol.resize(iCurSize);

		pir = &iRow[0];
		pic = &iCol[0];

#ifdef HAVE_MEMMOVE
		memmove(pir, piRow, sizeof(integer)*iNonZeroes);
		memmove(pic, piCol, sizeof(integer)*iNonZeroes);
#else /* ! HAVE_MEMMOVE */
		for (unsigned i = 0; i < iNonZeroes; i++) {
			pir[i] = piRow[i];
			pic[i] = piCol[i];
		}
#endif /* ! HAVE_MEMMOVE */

	} else {
		pir = piRow;
		pic = piCol;
	}

	y12prefactor(&iN, &iNonZeroes, pdMat,
			    pic, &iCurSize,
			    pir, &iCurSize,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);
			    
	if (iIFAIL != 0) {
		silent_cerr("Y12Solver (y12prefactor): "
			"error during pre-factorization, code " 
			<< iIFAIL << ":" << std::endl);
		PutError(iIFAIL);
		throw Y12Solver::ErrFactorization(iIFAIL, MBDYN_EXCEPT_ARGS);
	}

	/* actual factorization */
	y12factor(&iN, &iNonZeroes, pdMat,
			    pic, &iCurSize,
			    pir, &iCurSize,
			    pdPIVOT, LinearSolver::pdRhs,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);

	if (iIFAIL != 0) {
		silent_cerr("Y12Solver (y12factor): "
			"error during factorization, code " 
			<< iIFAIL << ":" << std::endl);
		PutError(iIFAIL);
		throw Y12Solver::ErrFactorization(iIFAIL, MBDYN_EXCEPT_ARGS);
	}

	if (dAFLAG[7] < 1.e-12) {
		silent_cerr("Y12Solver (y12factor):"
			" warning, possible bad conditioning of matrix" 
			<< std::endl);
	}
}

/* Risolve */
void
Y12Solver::Solve(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	
	if (bHasBeenReset) {
      		((Y12Solver *)this)->Factor();
	}
		
	integer iIFAIL = 0;

	y12solve(&iN, pdMat, &iCurSize, LinearSolver::pdRhs,
			    pdPIVOT, pic,
			    piHA, &iN,
			    iIFLAG, &iIFAIL);
	
	if (iIFAIL != 0) {
		silent_cerr("Y12Solver (y12solve): "
			"error during back substitution, code "
			<< iIFAIL << ":" << std::endl);
		PutError(iIFAIL);
		throw Y12Solver::ErrFactorization(iIFAIL, MBDYN_EXCEPT_ARGS);
	}
	
	if (bHasBeenReset) {
		iIFLAG[I_5] = 3;
      		bHasBeenReset = false;
	}
}

/* Index Form */
void
Y12Solver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Ar, std::vector<integer>& Ac,
		std::vector<integer>& Ap) const
{
	if (!bHasBeenReset) {
		return;
	}
	
	iNonZeroes = mh.MakeIndexForm(Ax, Ar, Ac, Ap, 1);
	ASSERT(iNonZeroes > 0);

	pdMat = &Ax[0];
	piRow = &Ar[0];
	piCol = &Ac[0];

	/* iCurSize should be between 3 and 5 times iNonZeroes ... */
	if (iCurSize > 5*iNonZeroes) {
		iCurSize = 5*iNonZeroes;
		
	} else if (iCurSize < 3*iNonZeroes) {
		if (iMaxSize < 5*iNonZeroes) {
			iCurSize = iMaxSize;
		} else {
			iCurSize = 5*iNonZeroes;
		}
	}
}

void 
Y12Solver::PutError(integer rc) const
{
	silent_cerr(std::endl);

	switch (rc) {
	case 1:
		silent_cerr("\tThe    coefficient   matrix   A   is   not" << std::endl
			<< "\tfactorized, i.e. the  call  of  subroutine" << std::endl
			<< "\tY12MD  was not preceded by a call of Y12MC" << std::endl
			<< "\tduring the solution of   Ax=b   or  during" << std::endl
			<< "\tthe  solution  of  the  first  system in a" << std::endl
			<< "\tsequence ( Ax1 = b1 , Ax2 = b2,.....,Axp =" << std::endl
			<< "\tbp)  of  systems with the same coefficient" << std::endl
			<< "\tmatrix. This will work in all  cases  only" << std::endl
			<< "\tif  the  user  sets IFLAG(1) .ge. 0 before" << std::endl
			<< "\tthe call of package Y12M (i.e. before  the" << std::endl
			<< "\tfirst   call   of  a  subroutine  of  this" << std::endl
			<< "\tpackage)." << std::endl);
		break;

	case 2:
		silent_cerr("\tThe coefficient matrix A is  not  ordered," << std::endl
			<< "\ti.e.  the call of subroutine Y12MC was not" << std::endl
			<< "\tpreceded by a call  of  Y12MB.  This  will" << std::endl
			<< "\twork  in  all  cases only if the user sets" << std::endl
			<< "\tIFLAG(1) .ge. 0 before the call of package" << std::endl
			<< "\tY12M  (i.e.  before  the  first  call of a" << std::endl
			<< "\tsubroutine of this package)." << std::endl);
		break;

	case 3:
		silent_cerr("\tA pivotal element abs(a(i,i;j)) < AFLAG(4)" << std::endl
			<< "\t*  AFLAG(6) is selected.  When AFLAG(4) is" << std::endl
			<< "\tsufficiently small this is  an  indication" << std::endl
			<< "\tthat the coefficient matrix is numerically" << std::endl
			<< "\tsingular." << std::endl);
		break;
		
	case 4:
		silent_cerr("\tAFLAG(5), the  growth  factor,  is  larger" << std::endl
			<< "\tthan    AFLAG(3).    When    AFLAG(3)   is" << std::endl
			<< "\tsufficiently large this indicates that the" << std::endl
			<< "\telements  of the coefficient matrix A grow" << std::endl
			<< "\tso quickly during the  factorization  that" << std::endl
			<< "\tthe continuation of the computation is not" << std::endl
			<< "\tjustified.  The  choice   of   a   smaller" << std::endl
			<< "\tstability   factor,   AFLAG(1),  may  give" << std::endl
			<< "\tbetter results in this case." << std::endl);
		break;

	case 5:
		silent_cerr("\tThe length NN of arrays A and SNR  is  not" << std::endl
			<< "\tsufficient.   Larger  values  of  NN  (and" << std::endl
			<< "\tpossibly of NN1) should be used." << std::endl);
		break;

	case 6:
		silent_cerr("\tThe  length  NN1  of  array  RNR  is   not" << std::endl
			<< "\tsufficient.   Larger  values  of  NN1 (and" << std::endl
			<< "\tpossibly of NN) should be used." << std::endl);
		break;

	case 7:
		silent_cerr("\tA row without  non-zero  elements  in  its" << std::endl
			<< "\tactive    part   is   found   during   the" << std::endl
			<< "\tdecomposition.  If   the   drop-tolerance," << std::endl
			<< "\tAFLAG(2),   is  sufficiently  small,  then" << std::endl
			<< "\tIFAIL = 7 indicates  that  the  matrix  is" << std::endl
			<< "\tnumerically  singular. If a large value of" << std::endl
			<< "\tthe drop-tolerance AFLAG(2) is used and if" << std::endl
			<< "\tIFAIL = 7  on exit, this is not certain. A" << std::endl
			<< "\trun  with  a  smaller  value  of  AFLAG(2)" << std::endl
			<< "\tand/or  a  careful check of the parameters" << std::endl
			<< "\tAFLAG(8) and AFLAG(5)  is  recommended  in" << std::endl
			<< "\tthe latter case." << std::endl);
		break;

	case 8:
		silent_cerr("\tA  column without non-zero elements in its" << std::endl
			<< "\tactive   part   is   found   during    the" << std::endl
			<< "\tdecomposition.   If   the  drop-tolerance," << std::endl
			<< "\tAFLAG(2),  is  sufficiently  small,   then" << std::endl
			<< "\tIFAIL  =  8  indicates  that the matrix is" << std::endl
			<< "\tnumerically singular. If a large value  of" << std::endl
			<< "\tthe drop-tolerance AFLAG(2) is used and if" << std::endl
			<< "\tIFAIL = 8  on exit, this is not certain. A" << std::endl
			<< "\trun  with  a  smaller  value  of  AFLAG(2)" << std::endl
			<< "\tand/or a careful check of  the  parameters" << std::endl
			<< "\tAFLAG(8)  and  AFLAG(5)  is recommended in" << std::endl
			<< "\tthe latter case." << std::endl);
		break;
	
	case 9:
		silent_cerr("\tA pivotal element  is  missing.  This  may" << std::endl
			<< "\toccur  if  AFLAG(2)  >  0 and IFLAG(4) = 2" << std::endl
			<< "\t(i.e. some system after the first one in a" << std::endl
			<< "\tsequence   of   systems   with   the  same" << std::endl
			<< "\tstructure is solved using a positive value" << std::endl
			<< "\tfor  the drop-tolerance). The value of the" << std::endl
			<< "\tdrop-tolerance   AFLAG(2),    should    be" << std::endl
			<< "\tdecreased  and  the  coefficient matrix of" << std::endl
			<< "\tthe system refactorized.  This  error  may" << std::endl
			<< "\talso occur when one of the special pivotal" << std::endl
			<< "\tstrategies (IFLAG(3)=0 or  IFLAG(3)=2)  is" << std::endl
			<< "\tused  and  the  matrix is not suitable for" << std::endl
			<< "\tsuch a strategy." << std::endl);
		break;

	case 10:
		silent_cerr("\tSubroutine Y12MF is called with IFLAG(5) =" << std::endl
			<< "\t1  (i.e.  with  a  request  to  remove the" << std::endl
			<< "\tnon-zero elements of the lower  triangular" << std::endl
			<< "\tmatrix    L).     IFLAG(5)=2     must   be" << std::endl
			<< "\tinitialized instead of IFLAG(5)=1." << std::endl);
		break;

	case 11:
		silent_cerr("\tThe coefficient matrix A contains at least" << std::endl
			<< "\ttwo  elements  in the same position (i,j)." << std::endl
			<< "\tThe  input   data   should   be   examined" << std::endl
			<< "\tcarefully in this case." << std::endl);
		break;

	case 12:
		silent_cerr("\tThe number of equations in the system Ax=b" << std::endl
			<< "\tis smaller than 2 (i.e.  N<2).  The  value" << std::endl
			<< "\tof N should be checked." << std::endl);
		break;
		
	case 13:
		silent_cerr("\tThe  number  of  non-zero  elements of the" << std::endl
			<< "\tcoefficient matrix is  non-positive  (i.e." << std::endl
			<< "\tZ.le.0  ).   The  value of the parameter Z" << std::endl
			<< "\t(renamed NZ in Y12MF) should be checked." << std::endl);
		break;

	case 14:
		silent_cerr("\tThe number of  non-zero  elements  in  the" << std::endl
			<< "\tcoefficient  matrix  is  smaller  than the" << std::endl
			<< "\tnumber of equations (i.e. Z  <  N  ).   If" << std::endl
			<< "\tthere  is no mistake (i.e. if parameter Z," << std::endl
			<< "\trenamed NZ in Y12MF, is correctly assigned" << std::endl
			<< "\ton  entry)  then the coefficient matrix is" << std::endl
			<< "\tstructurally singular in this case." << std::endl);
		break;

	case 15:
		silent_cerr("\tThe length IHA of the first  dimension  of" << std::endl
			<< "\tarray  HA  is  smaller  than  N.  IHA.ge.N" << std::endl
			<< "\tshould be assigned." << std::endl);
		break;

	case 16:
		silent_cerr("\tThe value of  parameter  IFLAG(4)  is  not" << std::endl
			<< "\tassigned  correctly.  IFLAG(4)  should  be" << std::endl
			<< "\tequal to 0, 1 or 2. See  more  details  in" << std::endl
			<< "\tthe description of this parameter." << std::endl);
		break;
		
	case 17:
		silent_cerr("\tA  row  without non-zero elements has been" << std::endl
			<< "\tfound in the coefficient matrix A  of  the" << std::endl
			<< "\tsystem  before the Gaussian elimination is" << std::endl
			<< "\tinitiated.  Matrix   A   is   structurally" << std::endl
			<< "\tsingular." << std::endl);
		break;

	case 18:
		silent_cerr("\tA  column  without  non-zero  elements has" << std::endl
			<< "\tbeen found in the coefficient matrix A  of" << std::endl
			<< "\tthe system before the Gaussian elimination" << std::endl
			<< "\tis initiated.  Matrix  A  is  structurally" << std::endl
			<< "\tsingular." << std::endl);
		break;

	case 19:
		silent_cerr("\tParameter  IFLAG(2) is smaller than 1. The" << std::endl
			<< "\tvalue of IFLAG(2)  should  be  a  positive" << std::endl
			<< "\tinteger (IFLAG(2) = 3 is recommended)." << std::endl);
		break;

	case 20:
		silent_cerr("\tParameter   IFLAG(3)   is  out  of  range." << std::endl
			<< "\tIFLAG(3) should be equal to 0, 1 or 2." << std::endl);
		break;

	case 21:
		silent_cerr("\tParameter  IFLAG(5)  is  out   of   range." << std::endl
			<< "\tIFLAG(5) should be equal to 1, 2 or 3 (but" << std::endl
			<< "\twhen IFLAG(5) = 3 Y12MB and  Y12MC  should" << std::endl
			<< "\tnot  be  called;  see also the message for" << std::endl
			<< "\tIFAIL = 22 below)." << std::endl);
		break;

	case 22:
		silent_cerr("\tEither  subroutine  Y12MB  or   subroutine" << std::endl
			<< "\tY12MC is called with IFLAG(5) = 3. Each of" << std::endl
			<< "\tthese subroutines should  be  called  with" << std::endl
			<< "\tIFLAG(5) equal to 1 or 2." << std::endl);
		break;

	case 23:
		silent_cerr("\tThe    number    of   allowed   iterations" << std::endl
			<< "\t(parameter IFLAG(11) when Y12MF  is  used)" << std::endl
			<< "\tis  smaller  than  2.   IFLAG(11)  .ge.  2" << std::endl
			<< "\tshould be assigned." << std::endl);
		break;

	case 24:
		silent_cerr("\tAt least one element whose  column  number" << std::endl
			<< "\tis  either larger than N or smaller than 1" << std::endl
			<< "\tis found." << std::endl);
		break;

	case 25:
		silent_cerr("\tAt least one element whose row  number  is" << std::endl
			<< "\teither  larger than N or smaller than 1 is" << std::endl
			<< "\tfound." << std::endl);
		break;

	default:
		silent_cerr("\t Unhandled code." << std::endl);
		break;
	}

	silent_cerr(std::endl);
}

/* Y12Solver - end */


/* Y12SparseSolutionManager - begin: code */

/* Costruttore */
Y12SparseSolutionManager::Y12SparseSolutionManager(integer iSize, 
		integer iWorkSpaceSize,
		const doublereal& dPivotFactor, bool bDupInd)
: iMatSize(iSize), 
iColStart(iSize + 1),
dVec(iSize),
MH(iSize),
VH(iSize, &dVec[0])
{
   	ASSERT(iSize > 0);
   	ASSERT(((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0)) || dPivotFactor == -1.0);


   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
		/*
		 * y12 requires at least 3*numzeros to store factors
		 * for multiple backsubs
		 */
      		iWorkSpaceSize = 3*iSize*iSize;
   	}
	
	integer iPivot;
	if (dPivotFactor == 0.) {
		iPivot = 0;

	} else {
		iPivot = 1;
	}

	iRow.reserve(iWorkSpaceSize);
	iCol.reserve(iWorkSpaceSize);
	dMat.reserve(iWorkSpaceSize);

   	SAFENEWWITHCONSTRUCTOR(SolutionManager::pLS, 
			       Y12Solver,
			       Y12Solver(iMatSize, iWorkSpaceSize,
					   &dVec[0], iPivot, bDupInd));
   
	pLS->SetSolutionManager(this);

#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
Y12SparseSolutionManager::~Y12SparseSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	/* Dealloca arrays */
}

#ifdef DEBUG
/* Test di validita' del manager */
void 
Y12SparseSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatSize > 0);

#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLS));
#endif /* DEBUG_MEMMANAGER */
   	pLS->IsValid();
}
#endif /* DEBUG */

void
Y12SparseSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

void
Y12SparseSolutionManager::MakeIndexForm(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
	pLS->MakeCompactForm(MH, dMat, iRow, iCol, iColStart);
}

/* Risolve il problema */
void
Y12SparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
   	MakeIndexForm();

#if 0
	std::cerr << "### after MakeIndexForm:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < iColStart.size(); i++) {
		std::cerr << i << " " << iColStart[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	unsigned c = 0;
	for (unsigned i = 0; i < dMat.size(); i++) {
		std::cerr << i << " " << iRow[i] << " " << c << " " << dMat[i] << std::endl;
		if (i == iColStart[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif

   	pLS->Solve();

#if 0
	std::cerr << "### after Solve:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < iColStart.size(); i++) {
		std::cerr << i << " " << iColStart[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	c = 0;
	for (unsigned i = 0; i < dMat.size(); i++) {
		std::cerr << i << " " << iRow[i] << " " << c << " " << dMat[i] << std::endl;
		if (i == iColStart[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif
}

/* Y12SparseSolutionManager - end */

/* Y12SparseCCSolutionManager - begin */

template <class CC>
Y12SparseCCSolutionManager<CC>::Y12SparseCCSolutionManager(integer Dim,
		integer dummy, doublereal dPivot)
: Y12SparseSolutionManager(Dim, dummy, dPivot, true),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
Y12SparseCCSolutionManager<CC>::~Y12SparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
Y12SparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
template <class CC>
void
Y12SparseCCSolutionManager<CC>::MakeIndexForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(MH, dMat, iRow, iCol, iColStart);

		if (Ac == 0) {
			SAFENEWWITHCONSTRUCTOR(Ac, CC,
					CC(dMat, iRow, iColStart));
		}

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
template <class CC>
void
Y12SparseCCSolutionManager<CC>::MatrInitialize(void)
{
	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
Y12SparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &MH;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class Y12SparseCCSolutionManager<CColMatrixHandler<1> >;
template class Y12SparseCCSolutionManager<DirCColMatrixHandler<1> >;

/* Y12SparseCCSolutionManager - end */

#endif /* USE_Y12 */

