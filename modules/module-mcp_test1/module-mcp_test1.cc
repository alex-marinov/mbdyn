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
  Copyright (C) 2022(-2022) all rights reserved.

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

#include "strnodead.h"
#include <sp_matvecass.h>

using namespace sp_grad;


class MCPTest1: virtual public Elem, public UserDefinedElem
{
public:
     MCPTest1(unsigned uLabel, const DofOwner *pDO,
              DataManager* pDM, MBDynParser& HP);
     virtual ~MCPTest1();
     virtual unsigned int iGetNumDof() const override;
     virtual DofOrder::Order GetDofType(unsigned int i) const override;
     virtual DofOrder::Order GetEqType(unsigned int i) const override;
     virtual DofOrder::Complementarity GetCompType(unsigned int i) const override;
     virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const override;
     virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const override;
     virtual unsigned int iGetNumPrivData(void) const override;
     virtual unsigned int iGetPrivDataIdx(const char *s) const override;
     virtual doublereal dGetPrivData(unsigned int i) const override;
     virtual void Output(OutputHandler& OH) const override;
     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;
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
            const VectorHandler& XPrimeCurr) override;
     template <typename T>
     inline void
     AssRes(SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const SpGradientVectorHandler<T>& XCurr,
            const SpGradientVectorHandler<T>& XPrimeCurr,
            enum SpFunctionCall func);
     int iGetNumConnectedNodes(void) const;
     void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override;
     void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                   SimulationEntity::Hints *ph) override;
     std::ostream& Restart(std::ostream& out) const override;
     virtual unsigned int iGetInitialNumDof(void) const override;
     virtual void
     InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
                   const VectorHandler& XCurr) override;
     SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr) override;
     template <typename T>
     inline void
     InitialAssRes(SpGradientAssVec<T>& WorkVec,
                   const SpGradientVectorHandler<T>& XCurr,
                   enum SpFunctionCall func);

private:
     void SaveVar(sp_grad::SpGradient& vg, doublereal& v) {}
     void SaveVar(sp_grad::GpGradProd& vg, doublereal& v) {}
     void SaveVar(doublereal vg, doublereal& v) {
          v = vg;
     }
     doublereal m[2], k[2], q[2], qdot[2], qddot[2], lambda;
};

MCPTest1::MCPTest1(unsigned uLabel, const DofOwner *pDO,
                   DataManager* pDM, MBDynParser& HP)
     :Elem(uLabel, flag(0)),
      UserDefinedElem(uLabel, pDO)
{
     // help
     if (HP.IsKeyWord("help")) {
          silent_cout(
               "\n"
               "Module:        mcp_test1\n"
               "\n"
               "	This element implements a mixed complementarity problem\n"
               " (real) m1\n"
               " (real) k1\n"
               " (real) q1\n"
               " (real) qdot1\n"
               " (real) m2\n"
               " (real) k2\n"
               " (real) q2\n"
               " (real) qdot2\n"
               "\n"
               << std::endl);

          if (!HP.IsArg()) {
               /*
                * Exit quietly if nothing else is provided
                */
               throw NoErr(MBDYN_EXCEPT_ARGS);
          }
     }

     for (integer i = 0; i < 2; ++i) {
          m[i] = HP.GetReal();
          k[i] = HP.GetReal();
          q[i] = HP.GetReal();
          qdot[i] = HP.GetReal();
     }

     SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

     lambda = 0.;
     qddot[0] = 0.;
     qddot[1] = 0.;
}

MCPTest1::~MCPTest1(void)
{
     // destroy private data
}

unsigned int MCPTest1::iGetNumDof(void) const
{
     return 5u;
}

DofOrder::Order MCPTest1::GetDofType(unsigned int i) const
{
     switch (i) {
     case 0:
     case 1:
     case 2:
     case 3:
          return DofOrder::DIFFERENTIAL;

     case 4:
          return DofOrder::ALGEBRAIC;

     default:
          ASSERT(0);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

DofOrder::Order MCPTest1::GetEqType(unsigned int i) const
{
     switch (i) {
     case 0:
     case 1:
     case 2:
     case 3:
          return DofOrder::DIFFERENTIAL;

     case 4:
          return DofOrder::ALGEBRAIC;

     default:
          ASSERT(0);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

DofOrder::Complementarity MCPTest1::GetCompType(unsigned int i) const
{
     switch (i) {
     case 0:
     case 1:
     case 2:
     case 3:
          return DofOrder::EQUALITY;

     case 4:
          return DofOrder::COMPLEMENTARY;

     default:
          ASSERT(0);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

std::ostream& MCPTest1::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": displacement [q1, q2]" << std::endl;
     out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": velocity [q1P, q2P]" << std::endl;
     out << prefix << iFirstIndex + 5 << "->" << iFirstIndex + 5 << ": reaction force [lambda]" << std::endl;
     return out;
}

std::ostream& MCPTest1::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": equation of motion [f1, f2]" << std::endl;
     out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": definition of acceleration [f3, f4]" << std::endl;
     out << prefix << iFirstIndex + 5 << "->" << iFirstIndex + 5 << ": definition of reaction force [f5]" << std::endl;

     return out;
}

unsigned int MCPTest1::iGetNumPrivData(void) const
{
     return 7u;
}

unsigned int MCPTest1::iGetPrivDataIdx(const char *s) const
{
     static const struct {
          unsigned int index;
          char name[8];
     } data[] = {
          { 1u, "q1" },
          { 2u, "q2" },
          { 3u, "q1P" },
          { 4u, "q2P" },
          { 5u, "q1PP" },
          { 6u, "q2PP" },
          { 7u, "lambda" }
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

doublereal MCPTest1::dGetPrivData(unsigned int i) const
{
     switch (i) {
     case 1:
     case 2:
          return q[i - 1];
     case 3:
     case 4:
          return qdot[i - 3];
     case 5:
     case 6:
          return qddot[i - 5];
     case 7:
          return lambda;

     default:
          silent_cerr("mcp_test1(" << GetLabel() << "): invalid private data index " << i << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

void
MCPTest1::Output(OutputHandler& OH) const
{
     if ( bToBeOutput() )
     {
          if ( OH.UseText(OutputHandler::LOADABLE) )
          {
               std::ostream& os = OH.Loadable();

               os << std::setw(8) << GetLabel();

               for (integer i = 0; i < 2; ++i)
               {
                    os << ' ' << q[i] << ' ' << qdot[i] << ' ' << qddot[i] << ' ';
               }

               os << lambda << '\n';
          }
     }
}

void
MCPTest1::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 5;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
MCPTest1::AssJac(VariableSubMatrixHandler& WorkMat,
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
MCPTest1::AssJac(VectorHandler& JacY,
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
MCPTest1::AssRes(SubVectorHandler& WorkVec,
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
MCPTest1::AssRes(SpGradientAssVec<T>& WorkVec,
                 doublereal dCoef,
                 const SpGradientVectorHandler<T>& XCurr,
                 const SpGradientVectorHandler<T>& XPrimeCurr,
                 enum SpFunctionCall func) {

     const integer iFirstIndex = iGetFirstIndex();

     T q[2], qdot[2], v[2], vdot[2], lambda;

     for (integer i = 0; i < 2; ++i) {
          XCurr.dGetCoef(iFirstIndex + i + 1, q[i], dCoef);
          XPrimeCurr.dGetCoef(iFirstIndex + i + 1, qdot[i], 1.);
          XCurr.dGetCoef(iFirstIndex + i + 3, v[i], dCoef);
          XPrimeCurr.dGetCoef(iFirstIndex + i + 3, vdot[i], 1.);
     }

     XCurr.dGetCoef(iFirstIndex + 5, lambda, 1.);

     T f1 = m[0] * vdot[0] + k[0] * q[0] + k[1] * (q[0] - q[1]);
     T f2 = m[1] * vdot[1] + k[1] * (q[1] - q[0]) - lambda;
     T f3 = qdot[0] - v[0];
     T f4 = qdot[1] - v[1];
     T f5 = -q[1];

     WorkVec.AddItem(iFirstIndex + 1, f1);
     WorkVec.AddItem(iFirstIndex + 2, f2);
     WorkVec.AddItem(iFirstIndex + 3, f3);
     WorkVec.AddItem(iFirstIndex + 4, f4);
     WorkVec.AddItem(iFirstIndex + 5, f5);

     for (integer i = 0; i < 2; ++i) {
          SaveVar(q[i], this->q[i]);
          SaveVar(qdot[i], this->qdot[i]);
          SaveVar(vdot[i], this->qddot[i]);
     }

     SaveVar(lambda, this->lambda);
}

int
MCPTest1::iGetNumConnectedNodes(void) const
{
     return 0;
}

void
MCPTest1::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(iGetNumConnectedNodes());
}

void
MCPTest1::SetValue(DataManager *pDM,
                   VectorHandler& X, VectorHandler& XP,
                   SimulationEntity::Hints *ph)
{
     const integer iFirstIndex = iGetFirstIndex();

     for (integer i = 1; i <= 2; ++i) {
          X.PutCoef(iFirstIndex + i, q[i - 1]);
          XP.PutCoef(iFirstIndex + i, qdot[i - 1]);
          X.PutCoef(iFirstIndex + i + 2, qdot[i - 1]);
          XP.PutCoef(iFirstIndex + i + 2, qddot[i - 1]);
     }
}

std::ostream&
MCPTest1::Restart(std::ostream& out) const
{
     return out;
}

unsigned int
MCPTest1::iGetInitialNumDof(void) const
{
     return 0u;
}

void
MCPTest1::InitialWorkSpaceDim(
     integer* piNumRows,
     integer* piNumCols) const
{
     *piNumRows = 0;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
MCPTest1::InitialAssJac(
     VariableSubMatrixHandler& WorkMat,
     const VectorHandler& XCurr)
{

     WorkMat.SetNullMatrix();

     return WorkMat;
}

SubVectorHandler&
MCPTest1::InitialAssRes(
     SubVectorHandler& WorkVec,
     const VectorHandler& XCurr)
{
     WorkVec.ResizeReset(0);

     return WorkVec;
}

bool mcp_test1_set(void)
{
     UserDefinedElemRead *rf = new UDERead<MCPTest1>;

     if (!SetUDE("mcp" "test1", rf))
     {
          delete rf;
          return false;
     }

     return true;
}

extern "C"
{

     int module_init(const char *module_name, void *pdm, void *php)
     {
          if (!mcp_test1_set())
          {
               silent_cerr("journal_bearing: "
                           "module_init(" << module_name << ") "
                           "failed" << std::endl);

               return -1;
          }

          return 0;
     }

}
