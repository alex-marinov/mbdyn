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

#include <iostream>
#include <iomanip>

#include "fullmh.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "naivemh.h"

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient_spmh.h"
#include "cscmhtpl.h"
#ifdef USE_TRILINOS
#undef HAVE_BLAS
#include "epetraspmh.h"
#include "Epetra_SerialComm.h"
#endif
#endif

static constexpr doublereal mat[5][5] = {
        { 11.,  0., 13.,  0., 15. },
        {  0., 22.,  0., 24.,  0. },
        { 31.,  0., 33.,  0., 35. },
        {  0., 42.,  0., 44.,  0. },
        { 51.,  0., 53.,  0., 55. }
};

int
main(void)
{
        std::vector<integer> perm(5), invperm(5);
        perm[0] = 4;
        perm[1] = 3;
        perm[2] = 2;
        perm[3] = 1;
        perm[4] = 0;
        for (int i = 0; i < 5; i++) {
                invperm[perm[i]] = i;
        }

        FullMatrixHandler fm(5);
        SpMapMatrixHandler spm(5, 5);
        NaiveMatrixHandler nm(5);
        NaivePermMatrixHandler npm(5, perm, invperm);

#ifdef USE_SPARSE_AUTODIFF
        SpGradientSparseMatrixHandler spgmh(5, 5);
#ifdef USE_TRILINOS
        Epetra_SerialComm comm;
        EpetraSparseMatrixHandler epmh(5, 5, 5, comm);
#endif
#endif

        for (int r = 0; r < 5; r++) {
                for (int c = 0; c < 5; c++) {
                        if (mat[r][c] != 0.) {
                                fm(r + 1, c + 1) = mat[r][c];
                                spm(r + 1, c + 1) = mat[r][c];
                                nm(r + 1, c + 1) = mat[r][c];
                                npm(r + 1, c + 1) = mat[r][c];
#ifdef USE_SPARSE_AUTODIFF
                                sp_grad::SpGradient g;
                                g.Reset(0., c + 1, mat[r][c]);
                                spgmh.AddItem(r + 1, g);
#ifdef USE_TRILINOS
                                epmh.AddItem(r + 1, g);
#endif
#endif
                        }
                }
        }

        std::cout << "matrix in full form: " << std::endl
                << fm << std::endl;

        std::cout << "matrix in sparse form: " << std::endl
                << spm << std::endl;

        std::cout << "matrix in naive form: " << std::endl
                << nm << std::endl;

        std::cout << "matrix in naive permuted form: " << std::endl
                << npm << std::endl;

#ifdef USE_SPARSE_AUTODIFF
        std::cout << "matrix in sparse gradient form: " << std::endl
                  << spgmh << std::endl;
#ifdef USE_TRILINOS
        epmh.FillComplete();

        std::cout << "matrix in Epetra sparse form:\n"
                  << epmh << std::endl;
#endif
#endif

        std::vector<doublereal> Ax0;
        std::vector<integer> Ai0, Ap0;
        spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);

#ifdef USE_SPARSE_AUTODIFF
        std::vector<doublereal> Ax0g, Ax1g, Ax0gT, Ax1gT;
        std::vector<integer> Ai0g, Ap0g, Ai1g, Ap1g, Ai0gT, Ap0gT, Ai1gT, Ap1gT;

        spgmh.MakeCompressedColumnForm(Ax0g, Ai0g, Ap0g, 0);
        spgmh.MakeCompressedColumnForm(Ax1g, Ai1g, Ap1g, 1);
        spgmh.MakeCompressedRowForm(Ax0gT, Ai0gT, Ap0gT, 0);
        spgmh.MakeCompressedRowForm(Ax1gT, Ai1gT, Ap1gT, 1);

        const CSCMatrixHandlerTpl<doublereal, integer, 0> csc0(&Ax0g.front(), &Ai0g.front(), &Ap0g.front(), spgmh.iGetNumCols(), Ai0g.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 1> csc1(&Ax1g.front(), &Ai1g.front(), &Ap1g.front(), spgmh.iGetNumCols(), Ai1g.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 0> csc0T(&Ax0gT.front(), &Ai0gT.front(), &Ap0gT.front(), spgmh.iGetNumRows(), Ai0gT.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 1> csc1T(&Ax1gT.front(), &Ai1gT.front(), &Ap1gT.front(), spgmh.iGetNumRows(), Ai1gT.size());
#ifdef USE_TRILINOS
        std::vector<doublereal> Ax0e, Ax1e, Ax0eT, Ax1eT;
        std::vector<integer> Ai0e, Ap0e, Ai1e, Ap1e, Ai0eT, Ap0eT, Ai1eT, Ap1eT;

        epmh.MakeCompressedColumnForm(Ax0e, Ai0e, Ap0e, 0);
        epmh.MakeCompressedColumnForm(Ax1e, Ai1e, Ap1e, 1);
        epmh.MakeCompressedRowForm(Ax0eT, Ai0eT, Ap0eT, 0);
        epmh.MakeCompressedRowForm(Ax1eT, Ai1eT, Ap1eT, 1);

        const CSCMatrixHandlerTpl<doublereal, integer, 0> epc0(&Ax0e.front(), &Ai0e.front(), &Ap0e.front(), epmh.iGetNumCols(), Ai0e.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 1> epc1(&Ax1e.front(), &Ai1e.front(), &Ap1e.front(), epmh.iGetNumCols(), Ai1e.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 0> epc0T(&Ax0eT.front(), &Ai0eT.front(), &Ap0eT.front(), epmh.iGetNumRows(), Ai0eT.size());
        const CSCMatrixHandlerTpl<doublereal, integer, 1> epc1T(&Ax1eT.front(), &Ai1eT.front(), &Ap1eT.front(), epmh.iGetNumRows(), Ai1eT.size());
#endif
#endif

        CColMatrixHandler<0> ccm0(Ax0, Ai0, Ap0);

        std::cout << "matrix in cc<0> form: " << std::endl
                << ccm0 << std::endl;

        std::cout << "matrix in cc<0> form again: " << std::endl;
        const CColMatrixHandler<0>& const_ccm0 = ccm0;
        for (int ir = 1; ir <= 5; ir++) {
                for (int ic = 1; ic <= 5; ic++) {
                        std::cout << std::setw(16) << const_ccm0(ir, ic);
                }
                std::cout << std::endl;
        }

        std::vector<doublereal> Ax1;
        std::vector<integer> Ai1, Ap1;
        spm.MakeCompressedColumnForm(Ax1, Ai1, Ap1, 1);

        CColMatrixHandler<1> ccm1(Ax1, Ai1, Ap1);

        std::cout << "matrix in cc<1> form: " << std::endl
                << ccm1 << std::endl;

        std::cout << "matrix in cc<1> form again: " << std::endl;
        const CColMatrixHandler<1>& const_ccm1 = ccm1;
        for (int ir = 1; ir <= 5; ir++) {
                for (int ic = 1; ic <= 5; ic++) {
                        std::cout << std::setw(16) << const_ccm1(ir, ic);
                }
                std::cout << std::endl;
        }

        DirCColMatrixHandler<0> dirm0(Ax0, Ai0, Ap0);
        std::cout << "matrix in dir<0> form: " << std::endl
                << dirm0 << std::endl;

        std::cout << "matrix in dir<0> form again: " << std::endl;
        const DirCColMatrixHandler<0>& const_dirm0 = dirm0;
        for (int ir = 1; ir <= 5; ir++) {
                for (int ic = 1; ic <= 5; ic++) {
                        std::cout << std::setw(16) << const_dirm0(ir, ic);
                }
                std::cout << std::endl;
        }

        DirCColMatrixHandler<1> dirm1(Ax1, Ai1, Ap1);
        std::cout << "matrix in dir<1> form: " << std::endl
                << dirm1 << std::endl;

        std::cout << "matrix in dir<1> form again: " << std::endl;
        const DirCColMatrixHandler<1>& const_dirm1 = dirm1;
        for (int ir = 1; ir <= 5; ir++) {
                for (int ic = 1; ic <= 5; ic++) {
                        std::cout << std::setw(16) << const_dirm1(ir, ic);
                }
                std::cout << std::endl;
        }

        std::cout << "***************************" << std::endl
                << "full matrix handler:" << std::endl;

        for (FullMatrixHandler::const_iterator i = fm.begin();
                i != fm.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "sparse map matrix handler:" << std::endl;

        for (SpMapMatrixHandler::const_iterator i = spm.begin();
                i != spm.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "naive sparse matrix handler:" << std::endl;
        for (NaiveMatrixHandler::const_iterator i = nm.begin();
                i != nm.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "naive permuted sparse matrix handler:" << std::endl;
        for (NaivePermMatrixHandler::const_iterator i = npm.begin();
                i != npm.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "column compressed <0> sparse matrix handler:" << std::endl;
        for (CColMatrixHandler<0>::const_iterator i = ccm0.begin();
                i != ccm0.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "column compressed <1> sparse matrix handler:" << std::endl;
        for (CColMatrixHandler<1>::const_iterator i = ccm1.begin();
                i != ccm1.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "dir column compressed <0> sparse matrix handler:" << std::endl;
        for (DirCColMatrixHandler<0>::const_iterator i = dirm0.begin();
                i != dirm0.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                << "dir column compressed <1> sparse matrix handler:" << std::endl;
        for (DirCColMatrixHandler<1>::const_iterator i = dirm1.begin();
                i != dirm1.end(); ++i)
        {
                std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
                if (mat[i->iRow][i->iCol] != i->dCoef) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

#ifdef USE_SPARSE_AUTODIFF
        const SpGradientSparseMatrixHandler& cspgmh(spgmh);
        
        std::cout << "***************************" << std::endl
                  << "gradient sparse matrix handler:" << std::endl;
        for (const auto& i: spgmh)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || cspgmh(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                  << "csc0 matrix handler:" << std::endl;
        for (const auto& i: csc0)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || csc0(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                  << "csc1 matrix handler:" << std::endl;
        for (const auto& i: csc1)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || csc1(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                  << "csc0T matrix handler:" << std::endl;
        for (const auto& i: csc0T)
        {
                std::cout << "(" << i.iCol << ", " << i.iRow << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iCol][i.iRow] != i.dCoef || csc0T(i.iCol + 1, i.iRow + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }


        std::cout << "***************************" << std::endl
                  << "csc1T matrix handler:" << std::endl;
        for (const auto& i: csc1T)
        {
                std::cout << "(" << i.iCol << ", " << i.iRow << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iCol][i.iRow] != i.dCoef || csc1T(i.iCol + 1, i.iRow + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }
#ifdef USE_TRILINOS
        std::cout << "***************************" << std::endl
                  << "Epetra sparse matrix handler:" << std::endl;

        const EpetraSparseMatrixHandler& cepmh(epmh);
        
        for (const auto& i: epmh)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || cepmh(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);
                }
        }

        std::cout << "***************************" << std::endl
                  << "Epetra csc0 matrix handler:" << std::endl;
        for (const auto& i: epc0)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || epc0(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                  << "Epetra csc1 matrix handler:" << std::endl;
        for (const auto& i: epc1)
        {
                std::cout << "(" << i.iRow << ", " << i.iCol << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iRow][i.iCol] != i.dCoef || epc1(i.iRow + 1, i.iCol + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }

        std::cout << "***************************" << std::endl
                  << "Epetra csc0^T matrix handler:" << std::endl;
        for (const auto& i: epc0T)
        {
                std::cout << "(" << i.iCol << ", " << i.iRow << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iCol][i.iRow] != i.dCoef || epc0T(i.iCol + 1, i.iRow + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }


        std::cout << "***************************" << std::endl
                  << "Epetra csc1^T matrix handler:" << std::endl;
        for (const auto& i: epc1T)
        {
                std::cout << "(" << i.iCol << ", " << i.iRow << ", " << i.dCoef << ")" << std::endl;
                if (mat[i.iCol][i.iRow] != i.dCoef || epc1T(i.iCol + 1, i.iRow + 1) != mat[i.iRow][i.iCol]) {
                        std::cout << "==> failed!" << std::endl;
                        ASSERT(0);                        
                }
        }
#endif
#endif
        return 0;
}
