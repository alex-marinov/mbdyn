/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

/* FIXME: gravity in modal elements is eXperimental; undefine to disable */
#define MODAL_USE_GRAVITY

#include "modalad.h"

ModalAd::ModalAd(unsigned int uL,
                 const ModalNodeAd* pR,
                 const Vec3& x0,
                 const Mat3x3& R0,
                 const DofOwner* pDO,
                 unsigned int NM,         /* numero modi */
                 unsigned int NI,         /* numero nodi d'interfaccia */
                 unsigned int NF,         /* numero nodi FEM */
                 doublereal dMassTmp,     /* inv. inerzia (m, m.stat., d'in.) */
                 const Vec3& STmp,
                 const Mat3x3& JTmp,
                 const std::vector<unsigned int>& uModeNumber,
                 MatNxN *pGenMass,
                 MatNxN *pGenStiff,
                 MatNxN *pGenDamp,
                 const std::vector<std::string>& IdFEMNodes,    /* label nodi FEM */
                 Mat3xN *pN,               /* posizione dei nodi FEM */
                 const std::vector<Modal::StrNodeData>& snd,
                 Mat3xN *pPHItStrNode,     /* forme modali nodi d'interfaccia */
                 Mat3xN *pPHIrStrNode,
                 Mat3xN *pModeShapest,     /* autovettori: servono a aeromodal */
                 Mat3xN *pModeShapesr,
                 Mat3xN *pInv3,            /* invarianti d'inerzia I3...I11 */
                 Mat3xN *pInv4,
                 Mat3xN *pInv5,
                 Mat3xN *pInv8,
                 Mat3xN *pInv9,
                 Mat3xN *pInv10,
                 Mat3xN *pInv11,
                 VecN *aa,
                 VecN *bb,
                 flag fOut)
:Elem(uL, fOut),
 Modal(uL, pR, x0, R0, pDO, NM, NI, NF, dMassTmp, STmp, JTmp, uModeNumber, pGenMass, pGenStiff, pGenDamp, IdFEMNodes, pN, snd, pPHItStrNode, pPHIrStrNode, pModeShapest, pModeShapesr, pInv3, pInv4, pInv5, pInv8, pInv9, pInv10, pInv11, aa, bb, fOut),
 pModalNode(pR)
{
     ASSERT(pModalNode == nullptr || pModalNode->GetStructNodeType() == StructNode::MODAL);
}

ModalAd::~ModalAd(void)
{
}

void
ModalAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = (pModalNode ? 12 : 0) + 2 * NModes  + 2 * NModes * NStrNodes + 12 * NStrNodes;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
ModalAd::AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr)
{
     DEBUGCOUT("Entering ModalAd::AssJac()" << std::endl);

     SpGradientSubMatrixHandler& WorkMatSp = WorkMat.SetSparseGradient();

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMatSp,
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
ModalAd::AssJac(VectorHandler& JacY,
                const VectorHandler& Y,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr,
                VariableSubMatrixHandler& WorkMat)
{
     DEBUGCOUT("Entering ModalAd::AssJac()" << std::endl);
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
ModalAd::AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr)
{
     DEBUGCOUT("Entering ModalAd::AssRes()" << std::endl);

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

void
ModalAd::UpdateStrNodeData(ModalAd::StrNodeData& oNode,
                           const sp_grad::SpColVector<doublereal, 3>& d1tot,
                           const sp_grad::SpMatrix<doublereal, 3, 3>& R1tot,
                           const sp_grad::SpColVector<doublereal, 3>& F,
                           const sp_grad::SpColVector<doublereal, 3>& M,
                           const sp_grad::SpMatrix<doublereal, 3, 3>& R2)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= 3; ++i) {
          oNode.d1tot(i) = d1tot(i);
          oNode.F(i) = F(i);
          oNode.M(i) = M(i);
     }

     for (index_type j = 1; j <= 3; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               oNode.R1tot(i, j) = R1tot(i, j);
               oNode.R2(i, j) = R2(i, j);
          }
     }
}

void
ModalAd::UpdateModalNode(const sp_grad::SpColVector<doublereal, 3>& xTmp,
                         const sp_grad::SpMatrix<doublereal, 3, 3>& RTmp)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= 3; ++i) {
          x(i) = xTmp(i);
     }

     for (index_type j = 1; j <= 3; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               R(i, j) = RTmp(i, j);
          }
     }
}

void
ModalAd::UpdateState(const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& aTmp,
                     const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& aPrimeTmp,
                     const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& bTmp,
                     const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& bPrimeTmp)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= a.iGetNumRows(); ++i) {
          a(i) = aTmp(i);
     }

     for (index_type i = 1; i <= aPrime.iGetNumRows(); ++i) {
          aPrime(i) = aPrimeTmp(i);
     }

     for (index_type i = 1; i <= b.iGetNumRows(); ++i) {
          b(i) = bTmp(i);
     }

     for (index_type i = 1; i <= bPrime.iGetNumRows(); ++i) {
          bPrime(i) = bPrimeTmp(i);
     }
}

void
ModalAd::UpdateInvariants(const sp_grad::SpColVector<doublereal, 3>& Inv3jajTmp,
                          const sp_grad::SpMatrix<doublereal, 3, 3>& Inv8jajTmp,
                          const sp_grad::SpMatrix<doublereal, 3, 3>& Inv9jkajakTmp)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= Inv3jajTmp.iGetNumRows(); ++i) {
          Inv3jaj(i) = Inv3jajTmp(i);
     }

     for (index_type j = 1; j <= Inv8jajTmp.iGetNumCols(); ++j) {
          for (index_type i = 1; i <= Inv8jajTmp.iGetNumRows(); ++i) {
               Inv8jaj(i, j) = Inv8jajTmp(i, j);
          }
     }

     for (index_type j = 1; j <= Inv9jkajakTmp.iGetNumCols(); ++j) {
          for (index_type i = 1; i <= Inv9jkajakTmp.iGetNumRows(); ++i) {
               Inv9jkajak(i, j) = Inv9jkajakTmp(i, j);
          }
     }
}

template <typename T>
void
ModalAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                const doublereal dCoef,
                const sp_grad::SpGradientVectorHandler<T>& XCurr,
                const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                const sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iModalIndex = iGetFirstIndex();

     SpColVector<T> a(NModes, 1), aPrime(NModes, 1), b(NModes, 1), bPrime(NModes, 1);

     for (unsigned int i = 1; i <= NModes; ++i) {
          XCurr.dGetCoef(iModalIndex + i, a(i), dCoef);
          XPrimeCurr.dGetCoef(iModalIndex + i, aPrime(i), 1.);
     }

     for (unsigned int i = 1; i <= NModes; ++i) {
          XCurr.dGetCoef(iModalIndex + NModes + i, b(i), dCoef);
          XPrimeCurr.dGetCoef(iModalIndex + NModes + i, bPrime(i), 1.);
     }

     const SpColVector<T> MaPP_CaP_Ka = -*pModalMass * bPrime - *pModalDamp * b - *pModalStiff * a;

     SpColVector<T, 3> Inv3jaj(3, NModes), Inv3jaPj(3, NModes), Inv3jaPPj(3, NModes);

     if (pInv3) {
          Inv3jaj = *pInv3 * a;
          Inv3jaPj = *pInv3 * b;
          Inv3jaPPj = *pInv3 * bPrime;
     }

     SpMatrix<T, 3, 3> Inv8jaj(3, 3, pInv8 ? NModes : 0), Inv8jaPj(3, 3, pInv8 ? NModes : 0);
     SpMatrix<T, 3> Inv5jaj(3, NModes, pInv5 ? NModes : 0), Inv5jaPj(3, NModes, pInv5 ? NModes: 0);
     SpMatrix<T, 3, 3> Inv9jkajak(3, 3, (pInv8 && pInv9) ? 2 * NModes * NModes : 0);
     SpMatrix<T, 3, 3> Inv9jkajaPk(3, 3, (pInv8 && pInv9) ? 2 * NModes * NModes: 0);
     SpMatrix<T, 3, 3> Inv10jaPj(3, 3, pInv10 ? NModes : 0);

     if (pInv5 || pInv8 || pInv9 || pInv10) {
          for (unsigned int jMode = 1; jMode <= NModes; jMode++)  {
               const T& a_jMode = a(jMode);
               const T& aP_jMode = b(jMode);

               if (pInv5) {
                    Inv5jaj += SubMatrix<1, 1, 3>(*pInv5, (jMode - 1) * NModes + 1, 1, NModes) * a_jMode;
                    Inv5jaPj += SubMatrix<1, 1, 3>(*pInv5, (jMode - 1) * NModes + 1, 1, NModes) * aP_jMode;
               }

               if (pInv8) {
                    Inv8jaj += SubMatrix<3, 3>(*pInv8, 1, 1, (jMode - 1) * 3 + 1, 1) * a_jMode;
                    Inv8jaPj += SubMatrix<3, 3>(*pInv8, 1, 1, (jMode - 1) * 3 + 1, 1) * aP_jMode;

                    if (pInv9) {
                         for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
                              const index_type iOffset = (jMode - 1) * 3 * NModes + (kMode - 1) * 3 + 1;
                              Inv9jkajak += SubMatrix<3, 3>(*pInv9, 1, 1, iOffset, 1) * a_jMode * a(kMode);
                              Inv9jkajaPk += SubMatrix<3, 3>(*pInv9, 1, 1, iOffset, 1) * a_jMode * b(kMode);
                         }
                    }
               }

               if (pInv10) {
                    Inv10jaPj += SubMatrix<3, 3>(*pInv10, 1, 1, (jMode - 1) * 3 + 1, 1) * aP_jMode;
               }
          }
     }

#ifdef MODAL_USE_GRAVITY
     /* forza di gravita' (decidere come inserire g) */
     /* FIXME: use a reasonable reference point where compute gravity */
     ::Vec3 GravityAcceleration(::Zero3);
     const bool bGravity = GravityOwner::bGetGravity(this->x, GravityAcceleration);
#endif /* MODAL_USE_GRAVITY */

     const integer iRigidIndex = pModalNode ? pModalNode->iGetFirstIndex() : -1;

     SpColVector<T, 3> x{this->x}, vP(3, 1), wP(3, 1), w(3, 1), RTw(3, 6);
     SpMatrix<T, 3, 3> R{this->R};
     SpColVector<T, 3> FTmp(3, 0), MTmp(3, 0);

     if (pModalNode) {
          SpColVector<T, 3> xP(3, 1), g(3, 1), gP(3, 1), v(3, 1);

          pModalNode->GetXCurr(x, dCoef, func);
          pModalNode->GetVCurr(xP, dCoef, func);
          pModalNode->GetgCurr(g, dCoef, func);
          pModalNode->GetgPCurr(gP, dCoef, func);
          pModalNode->GetXPPCurr(vP, dCoef, func);
          const ::Vec3& wr = pModalNode->GetWRef();
          pModalNode->GetWPCurr(wP, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               XCurr.dGetCoef(iRigidIndex + 6 + i, v(i), dCoef);
               XCurr.dGetCoef(iRigidIndex + 9 + i, w(i), dCoef);
          }

          pModalNode->GetRCurr(R, dCoef, func);

          RTw = Transpose(R) * w;

          SpMatrix<T, 3, 3> J(3, 3, (NModes + 1) * NModes);

          J = Inv7;

          if (pInv8) {
               J += Inv8jaj + Transpose(Inv8jaj);
               if (pInv9 != 0) {
                    J -= Inv9jkajak;
               }
          }

          SpMatrix<T, 3, 3> Jtmp = EvalUnique(J);

          J = (R * Jtmp) * Transpose(R);

          SpColVector<T, 3> STmp(3, NModes);

          STmp = Inv2;

          if (pInv3) {
               STmp += Inv3jaj;
          }

          const SpColVector<T, 3> S = R * STmp;

          FTmp = vP * -dMass + Cross(S, wP) - Cross(w, Cross(w, S));

          if (pInv3 != 0) {
               FTmp -= R * Inv3jaPPj + Cross(w, R * Inv3jaPj) * 2.;
          }

#ifdef MODAL_USE_GRAVITY
          if (bGravity) {
               FTmp += GravityAcceleration * dMass;
          }
#endif /* MODAL_USE_GRAVITY */

          MTmp = -Cross(S, vP) - J * wP - Cross(w, J * w);

          if (pInv4) {
               MTmp -= R * ((*pInv4) * bPrime);
          }

          if (pInv5) {
               MTmp -= R * (Inv5jaj * bPrime);
          }

          if (pInv8) {
               SpMatrix<T, 3, 3> Tmp = Inv8jaPj;
               if (pInv9 != 0) {
                    Tmp -= Inv9jkajaPk;
               }
               MTmp -= R * (Tmp * RTw * 2.);
          }

          if (pInv10) {
               SpColVector<T, 3> VTmp = (Inv10jaPj + Transpose(Inv10jaPj)) * RTw;
               if (pInv11) {
                    VTmp += Cross(w, (R * (*pInv11 * b)));
               }
               MTmp -= R * VTmp;
          }

#ifdef MODAL_USE_GRAVITY
          if (bGravity) {
               MTmp += Cross(S, GravityAcceleration);
          }
#endif
          const SpColVector<T, 3> f1 = v - xP;
          const SpColVector<T, 3> f2 = w - MatGVec(g) * gP - MatRVec(g) * wr;

          WorkVec.AddItem(iRigidIndex + 1, f1);
          WorkVec.AddItem(iRigidIndex + 4, f2);
     }

     for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
          const index_type jOffset = (jMode - 1) * 3 + 1;
          T d = MaPP_CaP_Ka(jMode);

          if (pInv3) {
               const SpColVector<T, 3> RInv3j = R * SubMatrix<3, 1>(*pInv3, 1, 1, jMode, 1);

               d -= Dot(RInv3j, vP);

#ifdef MODAL_USE_GRAVITY
               if (bGravity) {
                    d += Dot(RInv3j, GravityAcceleration);
               }
#endif /* MODAL_USE_GRAVITY */
          }

          if (pInv4) {
               SpColVector<T, 3> Inv4j(3, NModes);

               Inv4j = pInv4->GetVec(jMode);

               if (pInv5) {
                    Inv4j += Inv5jaj.GetCol(jMode);
                    d -= Dot(R * Inv5jaPj.GetCol(jMode), w) * 2.;
               }

               d -= Dot(R * Inv4j, wP);
          }

          if (pInv8 || pInv9 || pInv10) {
               SpMatrix<T, 3, 3> MatTmp(3, 3, NModes);

               if (pInv8) {
                    MatTmp += Transpose(SubMatrix<3, 3>(*pInv8, 1, 1, jOffset, 1));

                    if (pInv9) {
                         for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++) {
                              const index_type kOffset = (jMode - 1) * 3 * NModes + kModem1 * 3 + 1;

                              MatTmp -= SubMatrix<3, 3>(*pInv9, 1, 1, kOffset, 1) * a(kModem1 + 1);
                         }
                    }
               }

               if (pInv10) {
                    MatTmp += SubMatrix<3, 3>(*pInv10, 1, 1, jOffset, 1);
               }

               d += Dot(w, R * (MatTmp * RTw));
          }

          WorkVec.AddItem(iModalIndex + NModes + jMode, d);
     }

     for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
          WorkVec.AddItem(iModalIndex + iCnt, b(iCnt) - aPrime(iCnt));
     }

     for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
          const index_type iStrNodem1 = iStrNode - 1;
          const integer iNodeFirstMomIndex = SND[iStrNodem1].pNode->iGetFirstMomentumIndex();

          SpColVector<T, 3> PHIta(3, NModes), PHIra(3, NModes);

          for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
               const index_type iOffset = (jMode - 1) * NStrNodes + iStrNode;
               PHIta += SubMatrix<3, 1>(*pPHIt, 1, 1, iOffset, 1) * a(jMode);
               PHIra += SubMatrix<3, 1>(*pPHIr, 1, 1, iOffset, 1) * a(jMode);
          }

          const SpColVector<T, 3> d1tot = R * (PHIta + SND[iStrNodem1].OffsetFEM);
          const SpMatrix<T, 3, 3> R1tot = R * MatCrossVec(PHIra, 1.);

          SpColVector<T, 3> F(3, 1);

          for (index_type i = 1; i <= 3; ++i) {
               XCurr.dGetCoef(iModalIndex + 2 * NModes + 6 * iStrNodem1 + i, F(i), 1.);
          }

          SpColVector<T, 3> x2(3, 1);

          auto pStrNodem1 = dynamic_cast<const StructNodeAd*>(SND[iStrNodem1].pNode);
          
          pStrNodem1->GetXCurr(x2, dCoef, func);

          SpMatrix<T, 3, 3> R2(3, 3, 3);

          pStrNodem1->GetRCurr(R2, dCoef, func);

          const SpColVector<T, 3> dTmp2(R2 * SND[iStrNodem1].OffsetMB);

          if (pModalNode) {
               FTmp -= F;
               MTmp -= Cross(d1tot, F);
          }

          SpColVector<T, 3> vtemp = Transpose(R) * F;

          for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
               const index_type iOffset = (jMode - 1) * NStrNodes + iStrNode;
               const T d = Dot(-vtemp, pPHIt->GetVec(iOffset));
               WorkVec.AddItem(iModalIndex + NModes + jMode, d);
          }

          SpColVector<T, 3> M(3, 1);

          for (index_type i = 1; i <= 3; ++i) {
               XCurr.dGetCoef(iModalIndex + 2 * NModes + 6 * iStrNodem1 + 3 + i, M(i), 1.);
          }

          const SpMatrix<T, 3, 3> DeltaR = Transpose(R2) * R1tot;
          const SpColVector<T, 3> ThetaCurr = VecRotMat(DeltaR);
          const SpColVector<T, 3> R2_M = R2 * M;

          if (pModalNode) {
               MTmp -= R2_M;
          }

          vtemp = Transpose(R) * R2_M;

          for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
               const index_type iOffset = (jMode - 1) * NStrNodes + iStrNode;
               const T d = Dot(-vtemp, pPHIr->GetVec(iOffset));
               WorkVec.AddItem(iModalIndex + NModes + jMode, d);
          }

          ASSERT(dCoef != 0.);

          const SpColVector<T, 3> f1 = (x2 + dTmp2 - x - d1tot) / dCoef;
          const SpColVector<T, 3> f2 = ThetaCurr / -dCoef;

          WorkVec.AddItem(iModalIndex + 2 * NModes + 6 * iStrNodem1 + 1, f1);
          WorkVec.AddItem(iModalIndex + 2 * NModes + 6 * iStrNodem1 + 4, f2);

          const SpColVector<T, 3> MTmp2 = Cross(dTmp2, F) + R2_M;

          WorkVec.AddItem(iNodeFirstMomIndex + 1, F);
          WorkVec.AddItem(iNodeFirstMomIndex + 4, MTmp2);

          UpdateStrNodeData(SND[iStrNodem1], d1tot, R1tot, F, M, R2);
     }

     if (pModalNode) {
          WorkVec.AddItem(iRigidIndex + 7, FTmp);
          WorkVec.AddItem(iRigidIndex + 10, MTmp);
     }

     if (pModalNode) {
          UpdateModalNode(x, R);
     }

     UpdateState(a, aPrime, b, bPrime);
     UpdateInvariants(Inv3jaj, Inv8jaj, Inv9jkajak);
}
