/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

static doublereal mat[5][5] = {
	{ 11.,  0., 13.,  0., 15. },
	{  0., 22.,  0., 24.,  0. },
	{ 31.,  0., 33.,  0., 35. },
	{  0., 42.,  0., 44.,  0. },
	{ 51.,  0., 53.,  0., 55. }
};

void ReportMatScale(const char* title, bool fOK, const MatrixScaleBase& matScale, const MatrixHandler& mh, const double cond[2]);

template<typename T>
void ScaleMatrix(const char* title, MatrixScale<T>& matScale, T& mh)
{
	doublereal cond[2];

	//std::cout << "matrix before scaling:" << std::endl << mh << std::endl;

	cond[0] = mh.ConditionNumber();

	const bool fOK = matScale.ComputeScaleFactors(mh);

	matScale.ScaleMatrix(mh);

	cond[1] = mh.ConditionNumber();

	ReportMatScale(title, fOK, matScale, mh, cond);
}

int
main(int argc, char* argv[])
{
	SolutionManager::ScaleOpt scale;
	scale.uFlags |= SolutionManager::SCALEF_VERBOSE | SolutionManager::SCALEF_WARN;
	scale.iMaxIter = argc >= 2 ? atoi(argv[1]) : 100;
	scale.dTol = argc >= 3 ? atof(argv[2]) : sqrt(std::numeric_limits<doublereal>::epsilon());

	for (int i = 1; i < argc; ++i) {
	     if (0 == strcmp(argv[i], "-r")) {
		  std::mt19937 e1;
		  std::uniform_real_distribution<doublereal> uniform_dist(-1e16, 1e16);
		 
		  for (integer iRow = 0; iRow < 5; ++iRow) {
		       for (integer iCol = 0; iCol < 5; ++iCol) {
			    mat[iRow][iCol] = uniform_dist(e1);
		       }
		  }
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
		MatrixScale<SpMapMatrixHandler>* pMap;
	} matScale[] = {
			{ new RowMaxMatrixScale<NaiveMatrixHandler>(scale),
			  new RowMaxMatrixScale<NaivePermMatrixHandler>(scale),
			  new RowMaxMatrixScale<FullMatrixHandler>(scale),
			  new RowMaxMatrixScale<CColMatrixHandler<0> >(scale),
			  new RowMaxMatrixScale<CColMatrixHandler<1> >(scale),
			  new RowMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new RowMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new RowMaxMatrixScale<SpMapMatrixHandler>(scale)},
			{ new RowSumMatrixScale<NaiveMatrixHandler>(scale),
			  new RowSumMatrixScale<NaivePermMatrixHandler>(scale),
			  new RowSumMatrixScale<FullMatrixHandler>(scale),
			  new RowSumMatrixScale<CColMatrixHandler<0> >(scale),
			  new RowSumMatrixScale<CColMatrixHandler<1> >(scale),
			  new RowSumMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new RowSumMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new RowSumMatrixScale<SpMapMatrixHandler>(scale)},
			{ new ColMaxMatrixScale<NaiveMatrixHandler>(scale),
			  new ColMaxMatrixScale<NaivePermMatrixHandler>(scale),
			  new ColMaxMatrixScale<FullMatrixHandler>(scale),
			  new ColMaxMatrixScale<CColMatrixHandler<0> >(scale),
			  new ColMaxMatrixScale<CColMatrixHandler<1> >(scale),
			  new ColMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new ColMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new ColMaxMatrixScale<SpMapMatrixHandler>(scale)},
			{ new ColSumMatrixScale<NaiveMatrixHandler>(scale),
			  new ColSumMatrixScale<NaivePermMatrixHandler>(scale),
			  new ColSumMatrixScale<FullMatrixHandler>(scale),
			  new ColSumMatrixScale<CColMatrixHandler<0> >(scale),
			  new ColSumMatrixScale<CColMatrixHandler<1> >(scale),
			  new ColSumMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new ColSumMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new ColSumMatrixScale<SpMapMatrixHandler>(scale)},
			{ new LapackMatrixScale<NaiveMatrixHandler>(scale),
			  new LapackMatrixScale<NaivePermMatrixHandler>(scale),
			  new LapackMatrixScale<FullMatrixHandler>(scale),
			  new LapackMatrixScale<CColMatrixHandler<0> >(scale),
			  new LapackMatrixScale<CColMatrixHandler<1> >(scale),
			  new LapackMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new LapackMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new LapackMatrixScale<SpMapMatrixHandler>(scale)},
			{ new IterativeMatrixScale<NaiveMatrixHandler>(scale),
			  new IterativeMatrixScale<NaivePermMatrixHandler>(scale),
			  new IterativeMatrixScale<FullMatrixHandler>(scale),
			  new IterativeMatrixScale<CColMatrixHandler<0> >(scale),
			  new IterativeMatrixScale<CColMatrixHandler<1> >(scale),
			  new IterativeMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new IterativeMatrixScale<DirCColMatrixHandler<1> >(scale),
			  new IterativeMatrixScale<SpMapMatrixHandler>(scale)},
			{ new RowMaxColMaxMatrixScale<NaiveMatrixHandler>(scale),
			  new RowMaxColMaxMatrixScale<NaivePermMatrixHandler>(scale),
			  new RowMaxColMaxMatrixScale<FullMatrixHandler>(scale),
			  new RowMaxColMaxMatrixScale<CColMatrixHandler<0> >(scale),
			  new RowMaxColMaxMatrixScale<CColMatrixHandler<1> >(scale),
			  new RowMaxColMaxMatrixScale<DirCColMatrixHandler<0> >(scale),
			  new RowMaxColMaxMatrixScale<DirCColMatrixHandler<1> >(scale),
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

		nm.Reset();
		npm.Reset();
		fm.Reset();
		spm.Reset();
		for (unsigned ir = 0; ir < 5; ir++) {
			for (unsigned ic = 0; ic < 5; ic++) {
				if (mat[ir][ic] != 0.) {
					nm(ir + 1, ic + 1) = mat[ir][ic];
					npm(ir + 1, ic + 1) = mat[ir][ic];
					fm(ir + 1, ic + 1) = mat[ir][ic];
					spm(ir + 1, ic + 1) = mat[ir][ic];
				}
			}
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

		ScaleMatrix("Naive", *matScale[iMatScale].pNaive, nm);
		ScaleMatrix("NaivePerm", *matScale[iMatScale].pNaivePerm, npm);
		ScaleMatrix("Full", *matScale[iMatScale].pFull, fm);
		ScaleMatrix("Dir0", *matScale[iMatScale].pDirCCol0, dirccm0);
		ScaleMatrix("Dir1", *matScale[iMatScale].pDirCCol1, dirccm1);
		ScaleMatrix("CC0", *matScale[iMatScale].pCCol0, ccm0);
		ScaleMatrix("CC1", *matScale[iMatScale].pCCol1, ccm1);
		ScaleMatrix("Map", *matScale[iMatScale].pMap, spm);
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

	std::cout << "matrix after scaling:" << std::endl << mh << std::endl;
	std::cout << "------------------------------------------------------------" << std::endl;
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
