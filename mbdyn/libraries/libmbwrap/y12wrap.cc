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
			 std::vector<integer>*const piTmpRow, 
			 std::vector<integer>*const piTmpCol,
			 std::vector<doublereal>*const pdTmpMat,
			 doublereal* pdTmpRhs, 
			 integer iPivotParam)
: iMatSize(iSize),
iCurSize(iSize),
piRow(piTmpRow),
piCol(piTmpCol),
pdMat(pdTmpMat),
iN(iMatOrd),
iNonZeroes(0),
pdRhs(pdTmpRhs),
piHA(NULL),
pdPIVOT(NULL),
iFirstSol(-1)
{
	ASSERT(iMatSize > 0);
	ASSERT(piRow != NULL);
	ASSERT(piCol != NULL);
	ASSERT(pdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0);
	
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

        iIFLAG[I_1] = 0;
	iIFLAG[I_2] = 3;	/* recommended row number for pivoting */
	iIFLAG[I_3] = iPivotParam;
	iIFLAG[I_4] = 0;
	iIFLAG[I_5] = 2;	/* store non-zero elements of L */
					
        dAFLAG[I_1] = 8.;	/* Should be 4.<dAFLAG[0]<16. for stability */
	dAFLAG[I_2] = 0.;	/* Should be 0.<dAFLAG[1]<1.e-12 */
	dAFLAG[I_3] = 1.e6;	/* Should be dAFLAG[2]>1.e5 */
	dAFLAG[I_4] = 0.;   	/* FIXME: Should be 0<dAFLAG[3]<1.e-12 */
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
	ASSERT(iCurSize > 0 && iCurSize <= iMatSize);
	ASSERT(piRow != NULL);
	ASSERT(piCol != NULL);
	ASSERT(pdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0); 
	
	ASSERT(piHA != NULL);
	ASSERT(pdPIVOT != NULL);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsBlock(piHA, 11*iN*sizeof(integer)));
	ASSERT(defaultMemoryManager.fIsBlock(pdPIVOT, 1*iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}

bool
Y12LUSolver::SetCurSize(integer i)
{
	if (i < 1 || i > iMatSize) {
		return false;
	}

	iCurSize = i;

	return true;
}

integer
Y12LUSolver::iGetCurSize(void) const
{
	return iCurSize;
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
	iFirstSol = 1;

	/*
	 * must be set to 0 before the first call to a routine 
	 * of the y12m package
	 */
	iIFLAG[I_1] = 0;
	iIFLAG[I_5] = 2;
	
	__FC_DECL__(y12mbf)(&iN, &iNonZeroes, &((*pdMat)[0]),
			    &((*piCol)[0]), &iCurSize,
			    &((*piRow)[0]), &iCurSize,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);
			    
	if (iIFAIL != 0) {
		std::cerr << "Y12LUSolver (y12mbf): "
			"error during pre-factorization, code " 
			<< iIFAIL << ":" << std::endl;
		PutError(std::cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}

	/* actual factorization */
	__FC_DECL__(y12mcf)(&iN, &iNonZeroes, &((*pdMat)[0]),
			    &((*piCol)[0]), &iCurSize,
			    &((*piRow)[0]), &iCurSize,
			    pdPIVOT, pdRhs,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);

	if (iIFAIL != 0) {
		std::cerr << "Y12LUSolver (y12mcf): "
			"error during factorization, code " 
			<< iIFAIL << ":" << std::endl;
		PutError(std::cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}

	if (dAFLAG[7] < 1.e-12) {
		std::cerr << "Y12LUSolver (y12mcf):"
			" warning, possible bad conditioning of matrix" 
			<< std::endl;
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
	
	__FC_DECL__(y12mdf)(&iN, &((*pdMat)[0]), &iCurSize, pdRhs,
			    pdPIVOT, &((*piCol)[0]),
			    piHA, &iN,
			    iIFLAG, &iIFAIL);
	
	if (iIFAIL != 0) {
		std::cerr << "Y12LUSolver (y12mdf): "
			"error during back substitution, code "
			<< iIFAIL << ":" << std::endl;
		PutError(std::cerr, iIFAIL);
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	if (iFirstSol == 1) {
		iIFLAG[I_5] = 3;
		iFirstSol = 0;
	}
}

void 
Y12LUSolver::PutError(std::ostream& out, int rc) const
{
	out << std::endl;

	switch (rc ) {
	case 1:
		out 
			<< "\tThe    coefficient   matrix   A   is   not" << std::endl
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
			<< "\tpackage)." << std::endl;
		break;

	case 2:
		out
			<< "\tThe coefficient matrix A is  not  ordered," << std::endl
			<< "\ti.e.  the call of subroutine Y12MC was not" << std::endl
			<< "\tpreceded by a call  of  Y12MB.  This  will" << std::endl
			<< "\twork  in  all  cases only if the user sets" << std::endl
			<< "\tIFLAG(1) .ge. 0 before the call of package" << std::endl
			<< "\tY12M  (i.e.  before  the  first  call of a" << std::endl
			<< "\tsubroutine of this package)." << std::endl;
		break;

	case 3:
		out
			<< "\tA pivotal element abs(a(i,i;j)) < AFLAG(4)" << std::endl
			<< "\t*  AFLAG(6) is selected.  When AFLAG(4) is" << std::endl
			<< "\tsufficiently small this is  an  indication" << std::endl
			<< "\tthat the coefficient matrix is numerically" << std::endl
			<< "\tsingular." << std::endl;
		break;
		
	case 4:
		out
			<< "\tAFLAG(5), the  growth  factor,  is  larger" << std::endl
			<< "\tthan    AFLAG(3).    When    AFLAG(3)   is" << std::endl
			<< "\tsufficiently large this indicates that the" << std::endl
			<< "\telements  of the coefficient matrix A grow" << std::endl
			<< "\tso quickly during the  factorization  that" << std::endl
			<< "\tthe continuation of the computation is not" << std::endl
			<< "\tjustified.  The  choice   of   a   smaller" << std::endl
			<< "\tstability   factor,   AFLAG(1),  may  give" << std::endl
			<< "\tbetter results in this case." << std::endl;
		break;

	case 5:
		out
			<< "\tThe length NN of arrays A and SNR  is  not" << std::endl
			<< "\tsufficient.   Larger  values  of  NN  (and" << std::endl
			<< "\tpossibly of NN1) should be used." << std::endl;
		break;

	case 6:
		out
			<< "\tThe  length  NN1  of  array  RNR  is   not" << std::endl
			<< "\tsufficient.   Larger  values  of  NN1 (and" << std::endl
			<< "\tpossibly of NN) should be used." << std::endl;
		break;

	case 7:
		out
			<< "\tA row without  non-zero  elements  in  its" << std::endl
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
			<< "\tthe latter case." << std::endl;
		break;

	case 8:
		out
			<< "\tA  column without non-zero elements in its" << std::endl
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
			<< "\tthe latter case." << std::endl;
		break;
	
	case 9:
		out
			<< "\tA pivotal element  is  missing.  This  may" << std::endl
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
			<< "\tsuch a strategy." << std::endl;
		break;

	case 10:
		out
			<< "\tSubroutine Y12MF is called with IFLAG(5) =" << std::endl
			<< "\t1  (i.e.  with  a  request  to  remove the" << std::endl
			<< "\tnon-zero elements of the lower  triangular" << std::endl
			<< "\tmatrix    L).     IFLAG(5)=2     must   be" << std::endl
			<< "\tinitialized instead of IFLAG(5)=1." << std::endl;
		break;

	case 11:
		out
			<< "\tThe coefficient matrix A contains at least" << std::endl
			<< "\ttwo  elements  in the same position (i,j)." << std::endl
			<< "\tThe  input   data   should   be   examined" << std::endl
			<< "\tcarefully in this case." << std::endl;
		break;

	case 12:
		out
			<< "\tThe number of equations in the system Ax=b" << std::endl
			<< "\tis smaller than 2 (i.e.  N<2).  The  value" << std::endl
			<< "\tof N should be checked." << std::endl;
		break;
		
	case 13:
		out
			<< "\tThe  number  of  non-zero  elements of the" << std::endl
			<< "\tcoefficient matrix is  non-positive  (i.e." << std::endl
			<< "\tZ.le.0  ).   The  value of the parameter Z" << std::endl
			<< "\t(renamed NZ in Y12MF) should be checked." << std::endl;
		break;

	case 14:
		out
			<< "\tThe number of  non-zero  elements  in  the" << std::endl
			<< "\tcoefficient  matrix  is  smaller  than the" << std::endl
			<< "\tnumber of equations (i.e. Z  <  N  ).   If" << std::endl
			<< "\tthere  is no mistake (i.e. if parameter Z," << std::endl
			<< "\trenamed NZ in Y12MF, is correctly assigned" << std::endl
			<< "\ton  entry)  then the coefficient matrix is" << std::endl
			<< "\tstructurally singular in this case." << std::endl;
		break;

	case 15:
		out
			<< "\tThe length IHA of the first  dimension  of" << std::endl
			<< "\tarray  HA  is  smaller  than  N.  IHA.ge.N" << std::endl
			<< "\tshould be assigned." << std::endl;
		break;

	case 16:
		out
			<< "\tThe value of  parameter  IFLAG(4)  is  not" << std::endl
			<< "\tassigned  correctly.  IFLAG(4)  should  be" << std::endl
			<< "\tequal to 0, 1 or 2. See  more  details  in" << std::endl
			<< "\tthe description of this parameter." << std::endl;
		break;
		
	case 17:
		out
			<< "\tA  row  without non-zero elements has been" << std::endl
			<< "\tfound in the coefficient matrix A  of  the" << std::endl
			<< "\tsystem  before the Gaussian elimination is" << std::endl
			<< "\tinitiated.  Matrix   A   is   structurally" << std::endl
			<< "\tsingular." << std::endl;
		break;

	case 18:
		out
			<< "\tA  column  without  non-zero  elements has" << std::endl
			<< "\tbeen found in the coefficient matrix A  of" << std::endl
			<< "\tthe system before the Gaussian elimination" << std::endl
			<< "\tis initiated.  Matrix  A  is  structurally" << std::endl
			<< "\tsingular." << std::endl;
		break;

	case 19:
		out
			<< "\tParameter  IFLAG(2) is smaller than 1. The" << std::endl
			<< "\tvalue of IFLAG(2)  should  be  a  positive" << std::endl
			<< "\tinteger (IFLAG(2) = 3 is recommended)." << std::endl;
		break;

	case 20:
		out
			<< "\tParameter   IFLAG(3)   is  out  of  range." << std::endl
			<< "\tIFLAG(3) should be equal to 0, 1 or 2." << std::endl;
		break;

	case 21:
		out
			<< "\tParameter  IFLAG(5)  is  out   of   range." << std::endl
			<< "\tIFLAG(5) should be equal to 1, 2 or 3 (but" << std::endl
			<< "\twhen IFLAG(5) = 3 Y12MB and  Y12MC  should" << std::endl
			<< "\tnot  be  called;  see also the message for" << std::endl
			<< "\tIFAIL = 22 below)." << std::endl;
		break;

	case 22:
		out
			<< "\tEither  subroutine  Y12MB  or   subroutine" << std::endl
			<< "\tY12MC is called with IFLAG(5) = 3. Each of" << std::endl
			<< "\tthese subroutines should  be  called  with" << std::endl
			<< "\tIFLAG(5) equal to 1 or 2." << std::endl;
		break;

	case 23:
		out
			<< "\tThe    number    of   allowed   iterations" << std::endl
			<< "\t(parameter IFLAG(11) when Y12MF  is  used)" << std::endl
			<< "\tis  smaller  than  2.   IFLAG(11)  .ge.  2" << std::endl
			<< "\tshould be assigned." << std::endl;
		break;

	case 24:
		out
			<< "\tAt least one element whose  column  number" << std::endl
			<< "\tis  either larger than N or smaller than 1" << std::endl
			<< "\tis found." << std::endl;
		break;

	case 25:
		out
			<< "\tAt least one element whose row  number  is" << std::endl
			<< "\teither  larger than N or smaller than 1 is" << std::endl
			<< "\tfound." << std::endl;
		break;

	default:
		out
			<<"\t Unhandled code." << std::endl;
		break;
	}

	out << std::endl;
}

/* Y12LUSolver - end */


/* Y12SparseLUSolutionManager - begin: code */

/* Costruttore */
Y12SparseLUSolutionManager::Y12SparseLUSolutionManager(integer iSize, 
						       integer iWorkSpaceSize,
						       const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
// iRow(iWorkSpaceSize,0),
// iCol(iWorkSpaceSize,0), 
// dMat(iWorkSpaceSize,0.),
MH(iSize),
pVH(NULL),
pLU(NULL),
fHasBeenReset(1),
optimizeWorkSize(false)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));


   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
		/*
		 * y12 requires at least 3*numzeros to store factors
		 * for multiple backsubs
		 */
      		iWorkSpaceSize = 3*iSize*iSize;

		/*
		 * work size will be optimized
		 */
		optimizeWorkSize = true;

   	}
	
	integer iPivot;
	if (dPivotFactor == 0.) {
		iPivot = 0;
	} else {
		iPivot = 1;
	}

   	/* Alloca arrays */
	dVec.resize(iMatSize,0.);
   
   	/* Alloca handlers ecc. */
   	SAFENEWWITHCONSTRUCTOR(pVH,
			       MyVectorHandler,
			       MyVectorHandler(iMatSize, &(dVec[0])));
	iRow.reserve(iWorkSpaceSize);
	iCol.reserve(iWorkSpaceSize);
	dMat.reserve(iWorkSpaceSize);
   	SAFENEWWITHCONSTRUCTOR(pLU, 
			       Y12LUSolver,
			       Y12LUSolver(iMatSize, iWorkSpaceSize,
			       		   &iRow, &iCol,
					   &dMat, &(dVec[0]), iPivot));
   
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
   
   	/* Dealloca arrays */
}

/* Test di validita' del manager */
void 
Y12SparseLUSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatMaxSize > 0);
   	ASSERT(iMatSize > 0);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pVH));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLU));
#endif /* DEBUG_MEMMANAGER */
   
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
   	pLU->iNonZeroes = MH.MakeIndexForm(dMat,iRow,iCol,1);

#if 0
	/*
	 * This is not much important because we actually have
	 * all the space we need; the real optimization is in reducing
	 * the space for the sparse matrix.
	 */

	if (optimizeWorkSize) {
		integer cs = pLU->iGetCurSize();
		integer nz = pLU->iNonZeroes;
		integer ns = -1;

		/*
		 * y12m (to save the factored LU matrices) requires 
		 * 3 * nonzeros <= NN, NN1 <= 5 * nonzeros
		 */
		if (cs < 3*nz) {
			/* compromise (tune 3 => 5 ?) */
			ns = 5*nz;

		} else {
			/* asyntotically go to 4*nz */
			ns = (5*nz+cs)/2;
		}

		if (ns != -1 && ns != cs) {
			// pMH->SetCurSize(ns);	/* same as solver */
			pLU->SetCurSize(ns);
		}
		pMH->SetCurSize(2*nz);	/* optimal fill? */
	}
#endif
}

/* Inizializza il gestore delle matrici */
void
Y12SparseLUSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	MH.Init(dResetVal);
   	fHasBeenReset = flag(1);
}

/* Risolve il problema */
void
Y12SparseLUSolutionManager::Solve(const doublereal /* dCoef */)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	if (fHasBeenReset == 1) {
// 		std::fill(iRow.begin(),iRow.end(),0);
// 		std::fill(iCol.begin(),iCol.end(),0);
// 		std::fill(dMat.begin(),dMat.end(),0.);
// 		iRow.resize(iRow.capacity(),0);
// 		iCol.resize(iCol.capacity(),0);
// 		dMat.resize(dMat.capacity(),0.);
// 		std::fill(iRow.begin(),iRow.end(),0);
// 		std::fill(iCol.begin(),iCol.end(),0);
// 		std::fill(dMat.begin(),dMat.end(),0.);
      		PacVec();
      		if (pLU->fLUFactor() < 0) {	 
	 		THROW(Y12SparseLUSolutionManager::ErrGeneric());
      		}
	
      		fHasBeenReset = 0;
   	}

   	pLU->Solve();
}

/* Y12SparseLUSolutionManager - end */

#endif /* USE_Y12 */

