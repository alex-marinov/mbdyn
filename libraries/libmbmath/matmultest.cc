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

#include "sp_gradient_spmh.h"
#include "cscmhtpl.h"
#ifdef USE_TRILINOS
#undef HAVE_BLAS
#include "epetraspmh.h"
#include "epetravh.h"
#include <Epetra_SerialComm.h>
#endif

static constexpr doublereal mat[5][5] = {
        { 11.,  0., 13.,  0., 15. },
        {  0., 22.,  0., 24.,  0. },
        { 31.,  0., 33.,  0., 35. },
        {  0., 42.,  0., 44.,  0. },
        { 51.,  0., 53.,  0., 55. }
};

static int
check_mat(MatrixHandler& mh)
{
        for (int r = 0; r < 5; r++) {
                for (int c = 0; c < 5; c++) {
                        if (mat[r][c] != mh(r + 1, c + 1)) {
                                ASSERT(0);
                                return -1;
                        }
                }
        }

        return 0;
}

static int
check_mat_transpose(MatrixHandler& mh)
{
        for (int r = 0; r < 5; r++) {
                for (int c = 0; c < 5; c++) {
                        if (mat[c][r] != mh(r + 1, c + 1)) {
                                ASSERT(0);
                                return -1;
                        }
                }
        }

        return 0;
}

static int
check_vec(VectorHandler& vh, unsigned c)
{
        c--;
        for (int r = 0; r < 5; r++) {
                if (mat[r][c] != vh(r + 1)) {
                        ASSERT(0);
                        return -1;
                }
        }

        return 0;
}

static int
check_vec_transpose(VectorHandler& vh, unsigned r)
{
        r--;
        for (int c = 0; c < 5; c++) {
                if (mat[r][c] != vh(c + 1)) {
                        ASSERT(0);
                        return -1;
                }
        }

        return 0;
}

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

        SpMapMatrixHandler spm(5, 5);
        NaiveMatrixHandler nm(5);
        NaivePermMatrixHandler npm(5, perm, invperm);
        FullMatrixHandler fm(5, 5);

        SpGradientSparseMatrixHandler spgmh(5, 5);
#ifdef USE_TRILINOS
        Epetra_SerialComm comm;
        EpetraSparseMatrixHandler epmh(5, 5, 5, comm);
#endif
        for (int r = 0; r < 5; r++) {
                for (int c = 0; c < 5; c++) {
                        if (mat[r][c] != 0.) {
                                spm(r + 1, c + 1) = mat[r][c];
                                nm(r + 1, c + 1) = mat[r][c];
                                npm(r + 1, c + 1) = mat[r][c];
                                fm(r + 1, c + 1) = mat[r][c];
                                sp_grad::SpGradient g;
                                g.Reset(0., c + 1, mat[r][c]);
                                spgmh.AddItem(r + 1, g);
#ifdef USE_TRILINOS
                                epmh.AddItem(r + 1, g);
#endif
                        }
                }
        }

        std::cout << "matrix in sparse form: " << std::endl
                << spm << std::endl;

        std::cout << "matrix in naive form: " << std::endl
                << nm << std::endl;

        std::cout << "matrix in naive permuted form: " << std::endl
                << npm << std::endl;

        std::cout << "matrix in full form: " << std::endl
                << fm << std::endl;

        std::cout << "matrix in sparse gradient form: " << std::endl
                << spgmh << std::endl;

#ifdef USE_TRILINOS
        epmh.PacMat(); // Call before using operator()
        
        std::cout << "matrix in Epetra sparse form: " << std::endl
                  << epmh << std::endl;
#endif

        std::vector<doublereal> Ax0;
        std::vector<integer> Ai0, Ap0;
        spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);

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
        std::vector<doublereal> Ax0e, Ax1e, Ax0eT, Ax1eT;
        std::vector<integer> Ai0e, Ap0e, Ai1e, Ap1e, Ai0eT, Ap0eT, Ai1eT, Ap1eT;
        epmh.MakeCompressedColumnForm(Ax0e, Ai0e, Ap0e, 0);
        epmh.MakeCompressedColumnForm(Ax1e, Ai1e, Ap1e, 1);
        epmh.MakeCompressedRowForm(Ax0eT, Ai0eT, Ap0eT, 0);
        epmh.MakeCompressedRowForm(Ax1eT, Ai1eT, Ap1eT, 1);
        CSCMatrixHandlerTpl<doublereal, integer, 0> epcsc0(&Ax0e.front(), &Ai0e.front(), &Ap0e.front(), epmh.iGetNumCols(), Ai0e.size());
        CSCMatrixHandlerTpl<doublereal, integer, 1> epcsc1(&Ax1e.front(), &Ai1e.front(), &Ap1e.front(), epmh.iGetNumCols(), Ai1e.size());
        CSCMatrixHandlerTpl<doublereal, integer, 0> epcsc0T(&Ax0eT.front(), &Ai0eT.front(), &Ap0eT.front(), epmh.iGetNumRows(), Ai0eT.size());
        CSCMatrixHandlerTpl<doublereal, integer, 1> epcsc1T(&Ax1eT.front(), &Ai1eT.front(), &Ap1eT.front(), epmh.iGetNumRows(), Ai1eT.size());
#endif

        MyVectorHandler v(5), out(5);
        v.Reset();
        for (int i = 1; i <= 5; i++) {
                v(i) = 1.;

                spm.MatVecMul(out, v);
                std::cout << "sp*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                spm.MatTVecMul(out, v);
                std::cout << "sp^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                nm.MatVecMul(out, v);
                std::cout << "naive*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                nm.MatTVecMul(out, v);
                std::cout << "naive^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                npm.MatVecMul(out, v);
                std::cout << "naiveperm*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                npm.MatTVecMul(out, v);
                std::cout << "naiveperm^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                ccm0.MatVecMul(out, v);
                std::cout << "cc<0>*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                spgmh.MatVecMul(out, v);
                std::cout << "spgmh*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc0.MatVecMul(out, v);
                std::cout << "csc0*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc1.MatVecMul(out, v);
                std::cout << "csc1*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc0T.MatTVecMul(out, v);
                std::cout << "csc0T^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc1T.MatTVecMul(out, v);
                std::cout << "csc1T^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

#ifdef USE_TRILINOS
                epmh.MatVecMul(out, v);
                std::cout << "epmh^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc0.MatVecMul(out, v);
                std::cout << "epcsc0*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc1.MatVecMul(out, v);
                std::cout << "epcsc1*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc0T.MatTVecMul(out, v);
                std::cout << "epcsc0T^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc1T.MatTVecMul(out, v);
                std::cout << "epcsc1T^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }
#endif
                ccm0.MatTVecMul(out, v);
                std::cout << "cc<0>^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                ccm1.MatVecMul(out, v);
                std::cout << "cc<1>*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                ccm1.MatTVecMul(out, v);
                std::cout << "cc<1>^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                dirm0.MatVecMul(out, v);
                std::cout << "dir<0>*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                dirm0.MatTVecMul(out, v);
                std::cout << "dir<0>^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                dirm1.MatVecMul(out, v);
                std::cout << "dir<1>*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                dirm1.MatTVecMul(out, v);
                std::cout << "dir<1>^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                fm.MatVecMul(out, v);
                std::cout << "fm*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                fm.MatTVecMul(out, v);
                std::cout << "fm^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                spgmh.MatTVecMul(out, v);
                std::cout << "spgmh^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc0.MatTVecMul(out, v);
                std::cout << "csc0^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc1.MatTVecMul(out, v);
                std::cout << "csc1^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc0T.MatVecMul(out, v);
                std::cout << "csc0T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                csc1T.MatVecMul(out, v);
                std::cout << "csc1T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

#ifdef USE_TRILINOS
                epmh.MatTVecMul(out, v);
                std::cout << "epmh*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc0.MatTVecMul(out, v);
                std::cout << "epcsc0^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc1.MatTVecMul(out, v);
                std::cout << "epcsc1^T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc0T.MatVecMul(out, v);
                std::cout << "epcsc0T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }

                epcsc1T.MatVecMul(out, v);
                std::cout << "epcsc1T*v(" << i << ")=" << std::endl
                        << out << std::endl;
                if (check_vec_transpose(out, i)) {
                        std::cerr << "*** failed!" << std::endl;
                }
#endif

                v(i) = 0.;
        }

        FullMatrixHandler fmin(5, 5), fmout(5, 5);
        fmin.Reset();

        fmin(1, 1) = 1.;
        fmin(2, 2) = 1.;
        fmin(3, 3) = 1.;
        fmin(4, 4) = 1.;
        fmin(5, 5) = 1.;

        spm.MatMatMul(fmout, fmin);
        std::cout << "sp*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        spm.MatTMatMul(fmout, fmin);
        std::cout << "sp^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        nm.MatMatMul(fmout, fmin);
        std::cout << "naive*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        nm.MatTMatMul(fmout, fmin);
        std::cout << "naive^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        npm.MatMatMul(fmout, fmin);
        std::cout << "naiveperm*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        npm.MatTMatMul(fmout, fmin);
        std::cout << "naiveperm^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        ccm0.MatMatMul(fmout, fmin);
        std::cout << "cc<0>*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        ccm0.MatTMatMul(fmout, fmin);
        std::cout << "cc<0>^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        ccm1.MatMatMul(fmout, fmin);
        std::cout << "cc<1>*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        ccm1.MatTMatMul(fmout, fmin);
        std::cout << "cc<1>^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        dirm0.MatMatMul(fmout, fmin);
        std::cout << "dir<0>*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        dirm0.MatTMatMul(fmout, fmin);
        std::cout << "dir<0>^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        dirm1.MatMatMul(fmout, fmin);
        std::cout << "dir<1>*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        dirm1.MatTMatMul(fmout, fmin);
        std::cout << "dir<1>^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        fm.MatMatMul(fmout, fmin);
        std::cout << "fm*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        fm.MatTMatMul(fmout, fmin);
        std::cout << "fm^T*eye=" << std::endl
                << fmout << std::endl;
        if (check_mat_transpose(fmout)) {
                std::cerr << "*** failed!" << std::endl;
        }

        return 0;
}
