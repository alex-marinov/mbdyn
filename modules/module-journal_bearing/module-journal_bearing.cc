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

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2015(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <cstring>
#include <ctime>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-journal_bearing.h"

#ifdef USE_SPARSE_AUTODIFF
#include <sp_gradient.h>
#include <sp_matrix_base.h>
#include <sp_matvecass.h>

using namespace sp_grad;


class JournalBearing: virtual public Elem, public UserDefinedElem
{
public:
        JournalBearing(unsigned uLabel, const DofOwner *pDO,
                       DataManager* pDM, MBDynParser& HP);
        virtual ~JournalBearing(void);
        virtual unsigned int iGetNumDof(void) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
        virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
        virtual unsigned int iGetNumPrivData(void) const;
        virtual unsigned int iGetPrivDataIdx(const char *s) const;
        virtual doublereal dGetPrivData(unsigned int i) const;
        virtual void Output(OutputHandler& OH) const;
        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        virtual VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr);
        virtual void
        AssJac(VectorHandler& JacY,
               const VectorHandler& Y,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               VariableSubMatrixHandler& WorkMat) override;
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr);
        template <typename T>
        inline void
        AssRes(SpGradientAssVec<T>& WorkVec,
               doublereal dCoef,
               const SpGradientVectorHandler<T>& XCurr,
               const SpGradientVectorHandler<T>& XPrimeCurr,
               enum SpFunctionCall func);
        int iGetNumConnectedNodes(void) const;
        void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
        void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                      SimulationEntity::Hints *ph);
        std::ostream& Restart(std::ostream& out) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr);
        SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
        template <typename T>
        inline void
        InitialAssRes(SpGradientAssVec<T>& WorkVec,
                      const SpGradientVectorHandler<T>& XCurr,
                      enum SpFunctionCall func);

private:
        inline void SaveLambda(const sp_grad::SpColVector<doublereal, 2>& lambda);
        inline void SaveLambda(const sp_grad::SpColVector<sp_grad::SpGradient, 2>&) {}
        inline void SaveLambda(const sp_grad::SpColVector<sp_grad::GpGradProd, 2>&) {}
        inline void SaveFriction(doublereal omega, doublereal mf);
        inline void SaveFriction(const sp_grad::SpGradient&, const sp_grad::SpGradient&) {}
        inline void SaveFriction(const sp_grad::GpGradProd&, const sp_grad::GpGradProd&) {}

        StructNode* pNode1;
        SpColVectorA<doublereal, 3> o1;
        SpMatrixA<doublereal, 3, 3> e;
        StructNode* pNode2;
        SpColVectorA<doublereal, 3> o2;

        sp_grad::SpColVectorA<doublereal, 2> lambda;
        doublereal z, zP, omega, mf;
        doublereal d;
        doublereal sigma0, sigma1;
        doublereal muc, mus, vs, a, kv;

        static const index_type iWorkSpace = 14 + 1;
        static const index_type iInitialWorkSpace = 28;
};

JournalBearing::JournalBearing(
        unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP)
:       Elem(uLabel, flag(0)),
        UserDefinedElem(uLabel, pDO),
        pNode1(0),
        o1(Zero3),
        e(Eye3),
        pNode2(0),
        o2(Zero3),
        z(0.),
        zP(0.),
        omega(0.),
        mf(0.),
        d(0.),
        sigma0(0.),
        sigma1(0.),
        muc(0.),
        mus(0),
        vs(1.),
        a(1.),
        kv(0.)
{
        // help
        if (HP.IsKeyWord("help")) {
                silent_cout(
                        "\n"
                        "Module:        InlineAD\n"
                        "\n"
                        "	This element implements a journal bearing with lugre friction\n"
                        "\n"
                        "	journal bearing,\n"
                        "		node1, (label) <node1>,\n"
                        "			[ offset, (Vec3) <offset>, ]\n"
                        "			[ hinge, (Mat3x3) <orientation>, ]\n"
                        "		node2, (label) <node2>,\n"
                        "			[ offset, (Vec3) <offset>, ]\n"
                        "       [friction, "
                        "		 diameter, (real) <d>, \n"
                        "		 coulomb friction coefficient, (real) <muc>,\n"
                        "		 [static friction coefficient, (real) <mus>,]\n"
                        "		 micro stick stiffness, (real) <sigma0>\n"
                        "		 [,micro stick damping, (real) <sigma1>]\n"
                        "		 [,sliding velocity coefficient, (real) <vs>]\n"
                        "		 [,sliding velocity exponent, (real) <a>]\n"
                        "		 [,viscous friction coefficient, (real) <kv>]\n"
                        "\n"
                        << std::endl);

                if (!HP.IsArg()) {
                        /*
                         * Exit quietly if nothing else is provided
                         */
                        throw NoErr(MBDYN_EXCEPT_ARGS);
                }
        }

        if ( !HP.IsKeyWord("node1") ) {
                silent_cerr("journal bearing(" << GetLabel() << "): keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        pNode1 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

        if (!pNode1) {
                silent_cerr("journal bearing(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const ReferenceFrame refNode1(pNode1);

        if (HP.IsKeyWord("offset")) {
                o1 = HP.GetPosRel(refNode1);
        }

        if (HP.IsKeyWord("hinge") || HP.IsKeyWord("orientation")) {
                e = HP.GetRotRel(refNode1);
        }

        if (!HP.IsKeyWord("node2")) {
                silent_cerr("journal bearing(" << GetLabel() << "): keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        pNode2 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

        if (!pNode2) {
                silent_cerr("journal bearing(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("offset")) {
                const ReferenceFrame refNode2(pNode2);

                o2 = HP.GetPosRel(refNode2);
        }

        if (HP.IsKeyWord("friction")) {
                if (!HP.IsKeyWord("diameter")) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): keyword \"diameter\" expected at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                d = HP.GetReal();

                if (d <= 0.) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"diameter\" must be greater than zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!HP.IsKeyWord("coulomb" "friction" "coefficient")) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): keyword \"coulomb friction coefficient\" expected at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                muc = HP.GetReal();

                if (muc < 0) {
                        silent_cerr("\"coulomb friction coefficient\" "
                                    "must be greater than zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (HP.IsKeyWord("static" "friction" "coefficient")) {
                        mus = HP.GetReal();
                } else {
                        mus = muc;
                }

                if (mus < muc) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"static friction coefficient\" must be greater "
                                    "than or equal to \"coulomb friction coefficient\" "
                                    "at line " << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!HP.IsKeyWord("micro" "stick" "stiffness")) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): keyword \"micro stick stiffness\" expected at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                sigma0 = HP.GetReal();

                if (sigma0 <= 0.) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"micro stick stiffness\" must be greater than zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (HP.IsKeyWord("micro" "stick" "damping")) {
                        sigma1 = HP.GetReal();
                }

                if (sigma1 < 0.) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"micro stick damping\" must be greater than or equal to zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (HP.IsKeyWord("sliding" "velocity" "coefficient")) {
                        vs = HP.GetReal();
                }

                if (vs <= 0) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"sliding velocity coefficient "
                                    "must be greater than zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (HP.IsKeyWord("sliding" "velocity" "exponent")) {
                        a = HP.GetReal();
                }

                if (HP.IsKeyWord("viscous" "friction" "coefficient")) {
                        kv = HP.GetReal();
                }

                if (kv < 0.) {
                        silent_cerr("journal bearing(" << GetLabel()
                                    << "): \"viscous friction coefficient\" "
                                    "must be greater than or equal to zero at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (HP.IsKeyWord("initialize") && HP.GetYesNoOrBool()) {
                        using namespace sp_grad;
                        const Mat3x3& R1 = pNode1->GetRCurr();
                        const Vec3& omega1 = pNode1->GetWCurr();
                        const Vec3& omega2 = pNode2->GetWCurr();
                        const SpColVector<doublereal, 3> domega(R1.MulTV(omega2 - omega1));
                        const double domega_proj = Dot(e.GetCol(1), domega);
                        const double v = (0.5 * d) * domega_proj;

                        if (fabs(v) > 0.) {
                                const double g = muc + (mus - muc) * exp(-pow(fabs(v / vs), a));
                                z = v * g / (sigma0 * fabs(v));
                        }
                }
        }

        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

        std::ostream& out = pDM->GetLogFile();

        out << "journal bearing: " << GetLabel() << " "
                << pNode1->GetLabel() << " "
                << o1 << " "
                << e << " "
                << pNode2->GetLabel() << " " << o2 << " "
                << std::endl;
}

JournalBearing::~JournalBearing(void)
{
        // destroy private data
}

unsigned int JournalBearing::iGetNumDof(void) const
{
        return 2u + (muc > 0 ? 1u : 0u);
}

DofOrder::Order JournalBearing::GetDofType(unsigned int i) const
{
        switch (i) {
                case 0:
                case 1:
                        return DofOrder::ALGEBRAIC;

                case 2:
                        return DofOrder::DIFFERENTIAL;

                default:
                        ASSERT(0);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

DofOrder::Order JournalBearing::GetEqType(unsigned int i) const
{
        switch (i) {
                case 0:
                case 1:
                        return DofOrder::ALGEBRAIC;

                case 2:
                        return DofOrder::DIFFERENTIAL;

                default:
                        ASSERT(0);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

std::ostream& JournalBearing::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
        const integer iFirstIndex = iGetFirstIndex();

        out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": reaction forces [lambda1, lambda2]" << std::endl;

        if (bInitial) {
                out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": reaction force derivatives [lambdaP1, lambdaP2]" << std::endl;
        } else if (muc > 0) {
                out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 3 << ": sticktion state [z]" << std::endl;
        }

        return out;
}

std::ostream& JournalBearing::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
        const integer iFirstIndex = iGetFirstIndex();

        out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": position constraints [c1, c2]" << std::endl;

        if (bInitial) {
                out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": velocity constraints [cP1, cP2]" << std::endl;
        } else if (muc > 0) {
                out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 3 << ": sticktion state equation [f(z, zP)]" << std::endl;
        }

        return out;
}

unsigned int JournalBearing::iGetNumPrivData(void) const
{
        return 7u;
}

unsigned int JournalBearing::iGetPrivDataIdx(const char *s) const
{
        static const struct {
                unsigned int index;
                char name[8];
        } data[] = {
                { 1u, "lambda1" },
                { 2u, "lambda2" },
                { 3u, "z" },
                { 4u, "zP" },
                { 5u, "omega" },
                { 6u, "mf" },
                { 7u, "Pf" }
        };

        const unsigned N = iGetNumPrivData();

        ASSERT(N <= sizeof(data) / sizeof(data[0]));

        for (unsigned i = 0; i < N; ++i) {
                if (0 == strcmp(data[i].name, s)) {
                        return data[i].index;
                }
        }

        return 0;
}

doublereal JournalBearing::dGetPrivData(unsigned int i) const
{
        switch (i) {
                case 1:
                case 2:
                        return lambda(i);

                case 3:
                        return z;

                case 4:
                        return zP;

                case 5:
                        return omega;

                case 6:
                        return mf;

                case 7:
                        return -omega * mf;

                default:
                        silent_cerr("journal bearing(" << GetLabel() << "): invalid private data index " << i << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

void
JournalBearing::Output(OutputHandler& OH) const
{
        if ( bToBeOutput() )
        {
                if ( OH.UseText(OutputHandler::LOADABLE) )
                {
                        std::ostream& os = OH.Loadable();

                        os << std::setw(8) << GetLabel();

                        for (int i = 1; i <= lambda.iGetNumRows(); ++i)
                        {
                                os << " " << lambda(i);
                        }

                        if (muc > 0)
                        {
                                os << " " << z << " " << zP << " " << omega << " " << mf << " " << -omega * mf;
                        }

                        os << std::endl;
                }
        }
}

void
JournalBearing::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
        *piNumRows = iWorkSpace;
        *piNumCols = 0;
}

VariableSubMatrixHandler&
JournalBearing::AssJac(VariableSubMatrixHandler& WorkMat,
                       doublereal dCoef,
                       const VectorHandler& XCurr,
                       const VectorHandler& XPrimeCurr)
{
        SpGradientAssVec<SpGradient>::AssJac(this,
                                             WorkMat.SetSparseGradient(),
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC);

        return WorkMat;
}

void
JournalBearing::AssJac(VectorHandler& JacY,
                       const VectorHandler& Y,
                       doublereal dCoef,
                       const VectorHandler& XCurr,
                       const VectorHandler& XPrimeCurr,
                       VariableSubMatrixHandler& WorkMat)
{
        using namespace sp_grad;
     
        SpGradientAssVec<GpGradProd>::AssJac(this,
                                             JacY,
                                             Y,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
JournalBearing::AssRes(SubVectorHandler& WorkVec,
                       doublereal dCoef,
                       const VectorHandler& XCurr,
                       const VectorHandler& XPrimeCurr)
{
        SpGradientAssVec<doublereal>::AssRes(this,
                                             WorkVec,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_RES);

        return WorkVec;
}

template <typename T>
inline void
JournalBearing::AssRes(SpGradientAssVec<T>& WorkVec,
                       doublereal dCoef,
                       const SpGradientVectorHandler<T>& XCurr,
                       const SpGradientVectorHandler<T>& XPrimeCurr,
                       enum SpFunctionCall func) {

        const integer iFirstIndex = iGetFirstIndex();
        const integer iFirstMomentumIndexNode1 = pNode1->iGetFirstMomentumIndex();
        const integer iFirstMomentumIndexNode2 = pNode2->iGetFirstMomentumIndex();

        SpColVector<T, 3> X1(3, 1), X2(3, 1);
        SpColVector<T, 3> omega1(3, 3), omega2(3, 3);
        SpMatrix<T, 3, 3> R1(3, 3, 3), R2(3, 3, 3);

        pNode1->GetXCurr(X1, dCoef, func);
        pNode1->GetRCurr(R1, dCoef, func);
        pNode2->GetXCurr(X2, dCoef, func);
        pNode2->GetRCurr(R2, dCoef, func);

        SpColVector<T, 2> lambda(2, 1);

        XCurr.GetVec(iFirstIndex + 1, lambda, 1.); // Note: for algebraic variables dCoef is always one

        SaveLambda(lambda);

        const SpColVector<T, 3> R2o2 = R2 * o2;
        const SpColVector<T, 3> l1 = X2 + R2o2 - X1;

        const SpColVector<T, 3> F1 = R1 * (e.GetCol(2) * lambda(1) + e.GetCol(3) * lambda(2));
        SpColVector<T, 3> M1 = Cross(l1, F1);
        const SpColVector<T, 3> F2 = -F1;
        SpColVector<T, 3> M2 = Cross(R2o2, F2);

        if (muc > 0 || kv > 0)
        {
                pNode1->GetWCurr(omega1, dCoef, func);
                pNode2->GetWCurr(omega2, dCoef, func);

                const T domega = Dot(e.GetCol(1), Transpose(R1) * (omega2 - omega1));
                T mf = kv * domega;

                if (muc > 0)
                {
                        T z, zP;

                        XCurr.dGetCoef(iFirstIndex + 3, z, dCoef);
                        XPrimeCurr.dGetCoef(iFirstIndex + 3, zP, 1.);

                        T g;

                        const T v = (0.5 * d) * domega;

                        if (v != 0.) {
                                g = muc + (mus - muc) * exp(-pow(fabs(v / vs), a));
                        } else  {
                             SpGradientTraits<T>::ResizeReset(g, mus, 0);
                        }

                        const T f = v - sigma0 * fabs(v) / g * z - zP;
                        const T mu = sigma0 * z + sigma1 * zP;

                        if (lambda(1) != 0. || lambda(2) != 0) {
                                mf += (0.5 * d) * mu * sqrt(lambda(1) * lambda(1) + lambda(2) * lambda(2));
                        }

                        WorkVec.AddItem(iFirstIndex + 3, f);
                }

                SaveFriction(domega, mf);

                const SpColVector<T, 3> Mf = (R1 * e.GetCol(1)) * mf;

                M1 += Mf;
                M2 -= Mf;
        }

        WorkVec.AddItem(iFirstMomentumIndexNode1 + 1, F1);
        WorkVec.AddItem(iFirstMomentumIndexNode1 + 4, M1);
        WorkVec.AddItem(iFirstMomentumIndexNode2 + 1, F2);
        WorkVec.AddItem(iFirstMomentumIndexNode2 + 4, M2);

        const SpColVector<T, 3> a1 = Transpose(R1) * l1 - o1;

        for (integer i = 1; i <= 2; ++i) {
                WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a1) / dCoef);
        }
}

void JournalBearing::SaveLambda(const sp_grad::SpColVector<doublereal, 2>& lambda)
{
        this->lambda = lambda;
}

void JournalBearing::SaveFriction(doublereal omega, doublereal mf)
{
        this->omega = omega;
        this->mf = mf;
}

int
JournalBearing::iGetNumConnectedNodes(void) const
{
        return 2;
}

void
JournalBearing::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
        connectedNodes.resize(iGetNumConnectedNodes());
        connectedNodes[0] = pNode1;
        connectedNodes[1] = pNode2;
}

void
JournalBearing::SetValue(DataManager *pDM,
                         VectorHandler& X, VectorHandler& XP,
                         SimulationEntity::Hints *ph)
{
        const integer iFirstIndex = iGetFirstIndex();

        integer i;

        for (i = 1; i <= 2; ++i) {
                X.PutCoef(iFirstIndex + i, lambda(i));
        }

        if (muc > 0.) {
                X.PutCoef(iFirstIndex + i, z);
                XP.PutCoef(iFirstIndex + i, zP);
        }
}

std::ostream&
JournalBearing::Restart(std::ostream& out) const
{
        return out;
}

unsigned int
JournalBearing::iGetInitialNumDof(void) const
{
        return 4u;
}

void
JournalBearing::InitialWorkSpaceDim(
        integer* piNumRows,
        integer* piNumCols) const
{
        *piNumRows = iInitialWorkSpace;
        *piNumCols = 0;
}

VariableSubMatrixHandler&
JournalBearing::InitialAssJac(
        VariableSubMatrixHandler& WorkMat,
        const VectorHandler& XCurr)
{

        SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                    WorkMat.SetSparseGradient(),
                                                    XCurr,
                                                    INITIAL_ASS_JAC);

        return WorkMat;
}

SubVectorHandler&
JournalBearing::InitialAssRes(
        SubVectorHandler& WorkVec,
        const VectorHandler& XCurr)
{
        SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                    WorkVec,
                                                    XCurr,
                                                    INITIAL_ASS_RES);

        return WorkVec;
}

template <typename T>
inline void
JournalBearing::InitialAssRes(SpGradientAssVec<T>& WorkVec,
                              const SpGradientVectorHandler<T>& XCurr,
                              enum SpFunctionCall func) {
        SpColVector<T, 3> X1(3, 1), XP1(3, 1), X2(3, 1), XP2(3, 1), omega1(3, 3), omega2(3, 3);
        SpMatrix<T, 3, 3> R1(3, 3, 3), R2(3, 3, 3);

        pNode1->GetXCurr(X1, 1., func);	// Note: during initial assembly dCoef is always one
        pNode1->GetRCurr(R1, 1., func);
        pNode1->GetVCurr(XP1, 1., func);
        pNode1->GetWCurr(omega1, 1., func);

        pNode2->GetXCurr(X2, 1., func);
        pNode2->GetRCurr(R2, 1., func);
        pNode2->GetVCurr(XP2, 1., func);
        pNode2->GetWCurr(omega2, 1., func);

        const integer iFirstIndexNode1 = pNode1->iGetFirstIndex();
        const integer iFirstIndexNode2 = pNode2->iGetFirstIndex();
        const integer iFirstIndex = iGetFirstIndex();

        SpColVector<T, 2> lambda(2, 1), lambdaP(2, 1);

        for (integer i = 1; i <= 2; ++i) {
                XCurr.dGetCoef(iFirstIndex + i, lambda(i), 1.);
                XCurr.dGetCoef(iFirstIndex + i + 2, lambdaP(i), 1.);
        }

        SaveLambda(lambda);

        const SpColVector<T, 3> R2o2 = R2 * o2;
        const SpColVector<T, 3> l1 = X2 + R2o2 - X1;

        const SpColVector<T, 3> F1 = R1 * (e.GetCol(2) * lambda(1) + e.GetCol(3) * lambda(2));
        const SpColVector<T, 3> M1 = Cross(l1, F1);
        const SpColVector<T, 3> FP1 = Cross(omega1, R1 * (e.GetCol(2) * lambda(1) + e.GetCol(3) * lambda(2)))
                + R1 * (e.GetCol(2) * lambdaP(1) + e.GetCol(3) * lambdaP(2));
        const SpColVector<T, 3> MP1 = -Cross(F1, XP2 + Cross(omega2, R2o2) - XP1) + Cross(l1, FP1);
        const SpColVector<T, 3> F2 = -F1;
        const SpColVector<T, 3> M2 = Cross(R2o2, F2);
        const SpColVector<T, 3> FP2 = -FP1;
        const SpColVector<T, 3> MP2 = Cross(Cross(omega2, R2o2), F2) + Cross(R2o2, FP2);

        const SpColVector<T, 3> a = Transpose(R1) * l1 - o1;
        const SpColVector<T, 3> aP = Transpose(R1) * (Cross(l1, omega1) + XP2 + Cross(omega2, R2o2) - XP1);

        WorkVec.AddItem(iFirstIndexNode1 + 1, F1);
        WorkVec.AddItem(iFirstIndexNode1 + 4, M1);
        WorkVec.AddItem(iFirstIndexNode1 + 7, FP1);
        WorkVec.AddItem(iFirstIndexNode1 + 10, MP1);

        WorkVec.AddItem(iFirstIndexNode2 + 1, F2);
        WorkVec.AddItem(iFirstIndexNode2 + 4, M2);
        WorkVec.AddItem(iFirstIndexNode2 + 7, FP2);
        WorkVec.AddItem(iFirstIndexNode2 + 10, MP2);

        for (integer i = 1; i <= 2; ++i) {
                WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a));
                WorkVec.AddItem(iFirstIndex + i + 2, Dot(e.GetCol(i + 1), aP));
        }
}
#endif

bool journal_bearing_set(void)
{
#ifdef USE_SPARSE_AUTODIFF
        UserDefinedElemRead *rf = new UDERead<JournalBearing>;

        if (!SetUDE("journal" "bearing", rf))
        {
                delete rf;
                return false;
        }

        return true;
#else
        return false;
#endif
}

#ifndef STATIC_MODULES

extern "C"
{

int module_init(const char *module_name, void *pdm, void *php)
{
        if (!journal_bearing_set())
        {
                silent_cerr("journal_bearing: "
                            "module_init(" << module_name << ") "
                            "failed" << std::endl);

                return -1;
        }

        return 0;
}

}

#endif // ! STATIC_MODULE
