/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <random>

#include "dgeequ.h"

#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "naivemh.h"

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient_spmh.h"
#include "cscmhtpl.h"
#endif

#ifdef USE_TRILINOS
#undef HAVE_BLAS
#include "epetraspmh.h"
#include <Epetra_SerialComm.h>
#endif

static doublereal mat[5][5] = {
     { 11.,  0., 13.,  0., 15. },
     {  0., 22.,  0., 24.,  0. },
     { 31.,  0., 33.,  0., 35. },
     {  0., 42.,  0., 44.,  0. },
     { 51.,  0., 53.,  0., 55. }
};

static doublereal sol[5] = {1., 2., 3., 4., 5.};

static doublereal res[5] = {125., 140., 305., 260., 485};

void ReportMatScale(const char* title, bool fOK, const MatrixScaleBase& matScale, const MatrixHandler& mh, const double cond[2]);

template<typename T>
void ScaleMatrix(const char* title, MatrixScale<T>& matScale, T& mh, const VectorHandler& x, const VectorHandler& b, bool bTrans = false)
{
     const doublereal dTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.5);

     MyVectorHandler xo(mh.iGetNumCols()), xs(mh.iGetNumCols());
     MyVectorHandler bo(mh.iGetNumRows()), bs(mh.iGetNumRows());

     for (integer i = 1; i <= x.iGetSize(); ++i) {
          xo.PutCoef(i, x.dGetCoef(i));
     }

     if (bTrans) {
          mh.MatTVecMul(bo, xo);
     } else {
          mh.MatVecMul(bo, xo);
     }

     bo -= b;

     std::cout << "norm(bo) = " << bo.Norm() << std::endl;

     assert(bo.Norm() < dTol);

     doublereal cond[2];

     try {
          cond[0] = mh.ConditionNumber();
     } catch (const ErrNotImplementedYet&) {
          cond[0] = -1.;
     }

     const bool fOK = matScale.ComputeScaleFactors(mh);

     matScale.ScaleMatrix(mh);

     try {
          cond[1] = mh.ConditionNumber();
     } catch (const ErrNotImplementedYet&) {
          cond[1] = -1.;
     }

     ReportMatScale(title, fOK, matScale, mh, cond);

     const auto& rs = matScale.GetRowScale();
     const auto& cs = matScale.GetColScale();

     for (integer i = 1; i <= mh.iGetNumCols(); ++i) {
          xs.PutCoef(i, x.dGetCoef(i));
     }

     if (bTrans) {
          if (!rs.empty()) {
               for (integer i = 1; i <= mh.iGetNumCols(); ++i) {
                    xs(i) /= rs[i - 1];
               }
          }

          mh.MatTVecMul(bs, xs);

          if (!cs.empty()) {
               for (integer i = 1; i <= mh.iGetNumRows(); ++i) {
                    bs(i) /= cs[i - 1];
               }
          }
     } else {
          if (!cs.empty()) {
               for (integer i = 1; i <= mh.iGetNumCols(); ++i) {
                    xs(i) /= cs[i - 1];
               }
          }

          mh.MatVecMul(bs, xs);

          if (!rs.empty()) {
               for (integer i = 1; i <= mh.iGetNumRows(); ++i) {
                    bs(i) /= rs[i - 1];
               }
          }
     }

     bs -= b;

     std::cout << "norm(bs) = " << bs.Norm() << std::endl;
     assert(bs.Norm() < dTol);
}

int
main(int argc, char* argv[])
{
     SolutionManager::ScaleOpt scale;
     scale.uFlags |= SolutionManager::SCALEF_VERBOSE | SolutionManager::SCALEF_WARN;
     scale.iMaxIter = argc >= 2 ? atoi(argv[1]) : 100;
     scale.dTol = argc >= 3 ? atof(argv[2]) : sqrt(std::numeric_limits<doublereal>::epsilon());
     bool bRand = false;
     int iCount = 1;

     for (int i = 1; i < argc; ++i) {
          if (0 == strcmp(argv[i], "-r")) {
               bRand = true;
          }

          if (0 == strcmp(argv[i], "-c") && argc > i + 1) {
               iCount = atoi(argv[i + 1]);
          }
     }

     std::mt19937 e1;
     std::uniform_real_distribution<doublereal> uniform_dist(-1e3, 1e3);

     for (int i = 0; i < iCount; ++i) {
          if (bRand) {
               for (integer iRow = 0; iRow < 5; ++iRow) {
                    for (integer iCol = 0; iCol < 5; ++iCol) {
                         mat[iRow][iCol] = uniform_dist(e1);
                    }
               }

               for (integer iRow = 0; iRow < 5; ++iRow) {
                    sol[iRow] = uniform_dist(e1);
               }

               for (integer iRow = 0; iRow < 5; ++iRow) {
                    double d = 0.;
                    for (integer iCol = 0; iCol < 5; ++iCol) {
                         d += mat[iRow][iCol] * sol[iCol];
                    }
                    res[iRow] = d;
               }
          }

          struct {
               MatrixScale<NaiveMatrixHandler>* pNaive;
               MatrixScale<NaivePermMatrixHandler>* pNaivePerm;
               MatrixScale<FullMatrixHandler>* pFull;
               MatrixScale<CColMatrixHandler<0> >* pCCol0;
               MatrixScale<CColMatrixHandler<1> >* pCCol1;
               MatrixScale<DirCColMatrixHandler<0> >* pDirCCol0;
               MatrixScale<DirCColMatrixHandler<1> >* pDirCCol1;
#ifdef USE_SPARSE_AUTODIFF
               MatrixScale<SpGradientSparseMatrixHandler>* pGrad;
               MatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >* pCSC0;
               MatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >* pCSC1;
#ifdef USE_TRILINOS
               MatrixScale<EpetraSparseMatrixHandler>* pEpetra;
#endif
#endif
               MatrixScale<SpMapMatrixHandler>* pMap;
          } matScale[] = {
               { new RowMaxMatrixScale<NaiveMatrixHandler>(scale),
                 new RowMaxMatrixScale<NaivePermMatrixHandler>(scale),
                 new RowMaxMatrixScale<FullMatrixHandler>(scale),
                 new RowMaxMatrixScale<CColMatrixHandler<0> >(scale),
                 new RowMaxMatrixScale<CColMatrixHandler<1> >(scale),
                 new RowMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new RowMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new RowMaxMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new RowMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new RowMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new RowMaxMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new RowMaxMatrixScale<SpMapMatrixHandler>(scale)},
               { new RowSumMatrixScale<NaiveMatrixHandler>(scale),
                 new RowSumMatrixScale<NaivePermMatrixHandler>(scale),
                 new RowSumMatrixScale<FullMatrixHandler>(scale),
                 new RowSumMatrixScale<CColMatrixHandler<0> >(scale),
                 new RowSumMatrixScale<CColMatrixHandler<1> >(scale),
                 new RowSumMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new RowSumMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new RowSumMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new RowSumMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new RowSumMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new RowSumMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new RowSumMatrixScale<SpMapMatrixHandler>(scale)},
               { new ColMaxMatrixScale<NaiveMatrixHandler>(scale),
                 new ColMaxMatrixScale<NaivePermMatrixHandler>(scale),
                 new ColMaxMatrixScale<FullMatrixHandler>(scale),
                 new ColMaxMatrixScale<CColMatrixHandler<0> >(scale),
                 new ColMaxMatrixScale<CColMatrixHandler<1> >(scale),
                 new ColMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new ColMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new ColMaxMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new ColMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new ColMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new ColMaxMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new ColMaxMatrixScale<SpMapMatrixHandler>(scale)},
               { new ColSumMatrixScale<NaiveMatrixHandler>(scale),
                 new ColSumMatrixScale<NaivePermMatrixHandler>(scale),
                 new ColSumMatrixScale<FullMatrixHandler>(scale),
                 new ColSumMatrixScale<CColMatrixHandler<0> >(scale),
                 new ColSumMatrixScale<CColMatrixHandler<1> >(scale),
                 new ColSumMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new ColSumMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new ColSumMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new ColSumMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new ColSumMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new ColSumMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new ColSumMatrixScale<SpMapMatrixHandler>(scale)},
               { new LapackMatrixScale<NaiveMatrixHandler>(scale),
                 new LapackMatrixScale<NaivePermMatrixHandler>(scale),
                 new LapackMatrixScale<FullMatrixHandler>(scale),
                 new LapackMatrixScale<CColMatrixHandler<0> >(scale),
                 new LapackMatrixScale<CColMatrixHandler<1> >(scale),
                 new LapackMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new LapackMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new LapackMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new LapackMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new LapackMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new LapackMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new LapackMatrixScale<SpMapMatrixHandler>(scale)},
               { new IterativeMatrixScale<NaiveMatrixHandler>(scale),
                 new IterativeMatrixScale<NaivePermMatrixHandler>(scale),
                 new IterativeMatrixScale<FullMatrixHandler>(scale),
                 new IterativeMatrixScale<CColMatrixHandler<0> >(scale),
                 new IterativeMatrixScale<CColMatrixHandler<1> >(scale),
                 new IterativeMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new IterativeMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new IterativeMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new IterativeMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new IterativeMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new IterativeMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new IterativeMatrixScale<SpMapMatrixHandler>(scale)},
               { new RowMaxColMaxMatrixScale<NaiveMatrixHandler>(scale),
                 new RowMaxColMaxMatrixScale<NaivePermMatrixHandler>(scale),
                 new RowMaxColMaxMatrixScale<FullMatrixHandler>(scale),
                 new RowMaxColMaxMatrixScale<CColMatrixHandler<0> >(scale),
                 new RowMaxColMaxMatrixScale<CColMatrixHandler<1> >(scale),
                 new RowMaxColMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
                 new RowMaxColMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
#ifdef USE_SPARSE_AUTODIFF
                 new RowMaxColMaxMatrixScale<SpGradientSparseMatrixHandler>(scale),
                 new RowMaxColMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 0> >(scale),
                 new RowMaxColMaxMatrixScale<CSCMatrixHandlerTpl<doublereal, integer, 1> >(scale),
#ifdef USE_TRILINOS
                 new RowMaxColMaxMatrixScale<EpetraSparseMatrixHandler>(scale),
#endif
#endif
                 new RowMaxColMaxMatrixScale<SpMapMatrixHandler>(scale)}
          };

          const int N = sizeof(matScale)/sizeof(matScale[0]);

          for (int iMatScale = 0; iMatScale < N; ++iMatScale) {
               std::vector<integer> perm(5), invperm(5);
               perm[0] = 4;
               perm[1] = 3;
               perm[2] = 2;
               perm[3] = 1;
               perm[4] = 0;
               for (int i = 0; i < 5; i++) {
                    invperm[perm[i]] = i;
               }

               NaiveMatrixHandler nm(5);
               NaivePermMatrixHandler npm(5, perm, invperm);
               FullMatrixHandler fm(5);
               SpMapMatrixHandler spm(5, 5);
               MyVectorHandler x(5), b(5);

#ifdef USE_SPARSE_AUTODIFF
               SpGradientSparseMatrixHandler spgmh(5, 5);
#endif

#ifdef USE_TRILINOS
               Epetra_SerialComm oComm;
               EpetraSparseMatrixHandler epmh(5, 5, 5, oComm);
#endif
               nm.Reset();
               npm.Reset();
               fm.Reset();
               spm.Reset();
#ifdef USE_TRILINOS
               epmh.Reset();
#endif
               for (unsigned ir = 0; ir < 5; ir++) {
                    for (unsigned ic = 0; ic < 5; ic++) {
                         if (mat[ir][ic] != 0.) {
                              nm(ir + 1, ic + 1) = mat[ir][ic];
                              npm(ir + 1, ic + 1) = mat[ir][ic];
                              fm(ir + 1, ic + 1) = mat[ir][ic];
                              spm(ir + 1, ic + 1) = mat[ir][ic];
#ifdef USE_SPARSE_AUTODIFF
                              sp_grad::SpGradient g;
                              g.Reset(0., ic + 1, mat[ir][ic]);
                              spgmh.AddItem(ir + 1, g);
#ifdef USE_TRILINOS
                              epmh.AddItem(ir + 1, g);
#endif
#endif
                         }
                    }
               }

               for (integer i = 0; i < 5; ++i) {
                    x.PutCoef(i + 1, sol[i]);
                    b.PutCoef(i + 1, res[i]);
               }

               std::vector<doublereal> Ax0, Ax1, Axd0, Axd1;
               std::vector<integer> Ai0, Ai1, Ap0, Ap1, Aid0, Apd0, Aid1, Apd1;

               spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);
               spm.MakeCompressedColumnForm(Ax1, Ai1, Ap1, 1);
               spm.MakeCompressedColumnForm(Axd0, Aid0, Apd0, 0);
               spm.MakeCompressedColumnForm(Axd1, Aid1, Apd1, 1);

               CColMatrixHandler<0> ccm0(Ax0, Ai0, Ap0);
               CColMatrixHandler<1> ccm1(Ax1, Ai1, Ap1);
               DirCColMatrixHandler<0> dirccm0(Axd0, Aid0, Apd0);
               DirCColMatrixHandler<1> dirccm1(Axd1, Aid1, Apd1);

#ifdef USE_SPARSE_AUTODIFF
               std::vector<doublereal> Ax0g, Ax1g, Ax0gT, Ax1gT;
               std::vector<integer> Ai0g, Ap0g, Ai1g, Ap1g, Ai0gT, Ap0gT, Ai1gT, Ap1gT;
               spgmh.MakeCompressedColumnForm(Ax0g, Ai0g, Ap0g, 0);
               spgmh.MakeCompressedColumnForm(Ax1g, Ai1g, Ap1g, 1);
               spgmh.MakeCompressedRowForm(Ax0gT, Ai0gT, Ap0gT, 0);
               spgmh.MakeCompressedRowForm(Ax1gT, Ai1gT, Ap1gT, 1);
               CSCMatrixHandlerTpl<doublereal, integer, 0> csc0(&Ax0g.front(), &Ai0g.front(), &Ap0g.front(), spgmh.iGetNumCols(), Ai0g.size());
               CSCMatrixHandlerTpl<doublereal, integer, 1> csc1(&Ax1g.front(), &Ai1g.front(), &Ap1g.front(), spgmh.iGetNumCols(), Ai1g.size());
               CSCMatrixHandlerTpl<doublereal, integer, 0> csc0T(&Ax0gT.front(), &Ai0gT.front(), &Ap0gT.front(), spgmh.iGetNumRows(), Ai0gT.size());
               CSCMatrixHandlerTpl<doublereal, integer, 1> csc1T(&Ax1gT.front(), &Ai1gT.front(), &Ap1gT.front(), spgmh.iGetNumRows(), Ai1gT.size());
#ifdef USE_TRILINOS
               std::vector<doublereal> Axepc0, Axepc1, Axepr0, Axepr1;
               std::vector<integer> Aiepc0, Apepc0, Aiepc1, Apepc1, Aiepr0, Apepr0, Aiepr1, Apepr1;
               epmh.PacMat();
               epmh.MakeCompressedColumnForm(Axepc0, Aiepc0, Apepc0, 0);
               epmh.MakeCompressedColumnForm(Axepc1, Aiepc1, Apepc1, 1);
               epmh.MakeCompressedRowForm(Axepr0, Aiepr0, Apepr0, 0);
               epmh.MakeCompressedRowForm(Axepr1, Aiepr1, Apepr1, 1);
               CSCMatrixHandlerTpl<doublereal, integer, 0> epcsc0(&Axepc0.front(), &Aiepc0.front(), &Apepc0.front(), epmh.iGetNumCols(), Aiepc0.size());
               CSCMatrixHandlerTpl<doublereal, integer, 1> epcsc1(&Axepc1.front(), &Aiepc1.front(), &Apepc1.front(), epmh.iGetNumCols(), Aiepc1.size());
               CSCMatrixHandlerTpl<doublereal, integer, 0> epcsc0T(&Axepr0.front(), &Aiepr0.front(), &Apepr0.front(), epmh.iGetNumRows(), Aiepr0.size());
               CSCMatrixHandlerTpl<doublereal, integer, 1> epcsc1T(&Axepr1.front(), &Aiepr1.front(), &Apepr1.front(), epmh.iGetNumRows(), Aiepr1.size());
#endif
#endif
               ScaleMatrix("Naive", *matScale[iMatScale].pNaive, nm, x, b);
               ScaleMatrix("NaivePerm", *matScale[iMatScale].pNaivePerm, npm, x, b);
               ScaleMatrix("Full", *matScale[iMatScale].pFull, fm, x, b);
               ScaleMatrix("Dir0", *matScale[iMatScale].pDirCCol0, dirccm0, x, b);
               ScaleMatrix("Dir1", *matScale[iMatScale].pDirCCol1, dirccm1, x, b);
               ScaleMatrix("CC0", *matScale[iMatScale].pCCol0, ccm0, x, b);
               ScaleMatrix("CC1", *matScale[iMatScale].pCCol1, ccm1, x, b);
               ScaleMatrix("Map", *matScale[iMatScale].pMap, spm, x, b);
#ifdef USE_SPARSE_AUTODIFF
               ScaleMatrix("Grad", *matScale[iMatScale].pGrad, spgmh, x, b);
               ScaleMatrix("csc0", *matScale[iMatScale].pCSC0, csc0, x, b);
               ScaleMatrix("csc1", *matScale[iMatScale].pCSC1, csc1, x, b);
               ScaleMatrix("csc0^T", *matScale[iMatScale].pCSC0, csc0T, x, b, true);
               ScaleMatrix("csc1^T", *matScale[iMatScale].pCSC1, csc1T, x, b, true);
#ifdef USE_TRILINOS
               ScaleMatrix("epmh", *matScale[iMatScale].pEpetra, epmh, x, b);
               ScaleMatrix("epcsc0", *matScale[iMatScale].pCSC0, epcsc0, x, b);
               ScaleMatrix("epcsc0", *matScale[iMatScale].pCSC1, epcsc1, x, b);
               ScaleMatrix("epcsc^T", *matScale[iMatScale].pCSC0, epcsc0T, x, b, true);
               ScaleMatrix("epcsc^T", *matScale[iMatScale].pCSC1, epcsc1T, x, b, true);
#endif
#endif
          }

          for (unsigned i = 0; i < sizeof(matScale)/sizeof(matScale[0]); ++i) {
               delete matScale[i].pNaive;
               delete matScale[i].pNaivePerm;
               delete matScale[i].pFull;
               delete matScale[i].pCCol0;
               delete matScale[i].pCCol1;
               delete matScale[i].pDirCCol0;
               delete matScale[i].pDirCCol1;
               delete matScale[i].pMap;
#ifdef USE_SPARSE_AUTODIFF
               delete matScale[i].pGrad;
               delete matScale[i].pCSC0;
               delete matScale[i].pCSC1;
#ifdef USE_TRILINOS
               delete matScale[i].pEpetra;
#endif
#endif
          }
     }

     return 0;
}

void ReportMatScale(const char* title, bool fOK, const MatrixScaleBase& matScale, const MatrixHandler& mh, const double cond[2]) {
     std::cout << "------------------------------------------------------------" << std::endl;
     std::cout << title << ": " << (fOK ? "OK" : "NOK") << " : " << typeid(matScale).name() << std::endl;

     std::cout << "condition number before scaling:" << cond[0] << std::endl;
     std::cout << "condition number after scaling:" << cond[1] << std::endl;

     matScale.Report(std::cout);
     const std::vector<doublereal>& r = matScale.GetRowScale();
     const std::vector<doublereal>& c = matScale.GetColScale();

     const int N = std::min(mh.iGetNumRows(), mh.iGetNumCols());

     for (int i = 0; i < N; ++i) {
          std::cout
               << "   r[" << i << "]=" << std::setw(12) << (r.empty() ? 1. : r[i])
               << "       c[" << i << "]=" << std::setw(12) << (c.empty() ? 1. : c[i])
               << std::endl;
     }

     std::cout << "matrix after scaling:\n";

     mh.Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);

     std::cout << "\n------------------------------------------------------------\n";
}

template
class RowMaxMatrixScale<SpMapMatrixHandler>;

template
class RowMaxMatrixScale<DirCColMatrixHandler<0> >;

template
class RowMaxMatrixScale<CColMatrixHandler<0> >;

template
class RowMaxMatrixScale<NaiveMatrixHandler>;

template
class RowMaxMatrixScale<NaivePermMatrixHandler>;

template
class RowMaxMatrixScale<FullMatrixHandler>;

template
class RowSumMatrixScale<SpMapMatrixHandler>;

template
class RowSumMatrixScale<DirCColMatrixHandler<0> >;

template
class RowSumMatrixScale<CColMatrixHandler<0> >;

template
class RowSumMatrixScale<NaiveMatrixHandler>;

template
class RowSumMatrixScale<NaivePermMatrixHandler>;

template
class RowSumMatrixScale<FullMatrixHandler>;

template
class ColMaxMatrixScale<SpMapMatrixHandler>;

template
class ColMaxMatrixScale<DirCColMatrixHandler<0> >;

template
class ColMaxMatrixScale<CColMatrixHandler<0> >;

template
class ColMaxMatrixScale<NaiveMatrixHandler>;

template
class ColMaxMatrixScale<NaivePermMatrixHandler>;

template
class ColMaxMatrixScale<FullMatrixHandler>;

template
class ColSumMatrixScale<SpMapMatrixHandler>;

template
class ColSumMatrixScale<DirCColMatrixHandler<0> >;

template
class ColSumMatrixScale<CColMatrixHandler<0> >;

template
class ColSumMatrixScale<NaiveMatrixHandler>;

template
class ColSumMatrixScale<NaivePermMatrixHandler>;

template
class ColSumMatrixScale<FullMatrixHandler>;

template
class LapackMatrixScale<SpMapMatrixHandler>;

template
class LapackMatrixScale<DirCColMatrixHandler<0> >;

template
class LapackMatrixScale<CColMatrixHandler<0> >;

template
class LapackMatrixScale<NaiveMatrixHandler>;

template
class LapackMatrixScale<NaivePermMatrixHandler>;

template
class LapackMatrixScale<FullMatrixHandler>;

template
class IterativeMatrixScale<SpMapMatrixHandler>;

template
class IterativeMatrixScale<DirCColMatrixHandler<0> >;

template
class IterativeMatrixScale<CColMatrixHandler<0> >;

template
class IterativeMatrixScale<NaiveMatrixHandler>;

template
class IterativeMatrixScale<NaivePermMatrixHandler>;

template
class IterativeMatrixScale<FullMatrixHandler>;

template
class RowMaxColMaxMatrixScale<SpMapMatrixHandler>;

template
class RowMaxColMaxMatrixScale<DirCColMatrixHandler<0> >;

template
class RowMaxColMaxMatrixScale<CColMatrixHandler<0> >;

template
class RowMaxColMaxMatrixScale<NaiveMatrixHandler>;

template
class RowMaxColMaxMatrixScale<NaivePermMatrixHandler>;

template
class RowMaxColMaxMatrixScale<FullMatrixHandler>;
