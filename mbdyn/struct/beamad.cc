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

#include "beamad.h"
#include "shapefnc.h"

BeamAd::BeamAd(unsigned int uL,
               const StructNodeAd* pN1,
               const StructNodeAd* pN2,
               const StructNodeAd* pN3,
               const Vec3& F1,
               const Vec3& F2,
               const Vec3& F3,
               const Mat3x3& R1,
               const Mat3x3& R2,
               const Mat3x3& R3,
               const Mat3x3& r_I,
               const Mat3x3& rII,
               const ConstitutiveLaw6D* pD_I,
               const ConstitutiveLaw6D* pDII,
               OrientationDescription ood,
               flag fOut)
:Elem(uL, fOut),
 Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut),
 pNode{pN1, pN2, pN3}
{

}

BeamAd::~BeamAd()
{
}

void BeamAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 18;
     *piNumCols = 0;
}

void
BeamAd::AddInternalForces(sp_grad::SpColVector<doublereal, 6>& AzLoc, unsigned int iSez)
{
}

void
BeamAd::AddInternalForces(sp_grad::SpColVector<sp_grad::SpGradient, 6>& AzLoc, unsigned int iSez)
{
}

void
BeamAd::AddInternalForces(sp_grad::SpColVector<sp_grad::GpGradProd, 6>& AzLoc, unsigned int iSez)
{
}

template <typename T>
sp_grad::SpColVector<T, 3>
BeamAd::InterpState(const sp_grad::SpColVector<T, 3>& v1,
                    const sp_grad::SpColVector<T, 3>& v2,
                    const sp_grad::SpColVector<T, 3>& v3,
                    Section Sec)
{
     using namespace sp_grad;

     SpColVector<T, 3> p(3, 0);

     for (index_type i = 1; i <= 3; ++i) {
          p(i) = v1(i) * dN3[Sec][0] + v2(i) * dN3[Sec][1] + v3(i) * dN3[Sec][2];
     }

     return p;
}

template <typename T>
sp_grad::SpColVector<T, 3>
BeamAd::InterpDeriv(const sp_grad::SpColVector<T, 3>& v1,
                    const sp_grad::SpColVector<T, 3>& v2,
                    const sp_grad::SpColVector<T, 3>& v3,
                    Section Sec)
{
     using namespace sp_grad;
     SpColVector<T, 3> g(3, 0);

     for (index_type i = 1; i <= 3; ++i) {
          g(i) = (v1(i) * dN3P[Sec][0] + v2(i) * dN3P[Sec][1] + v3(i) * dN3P[Sec][2]) * dsdxi[Sec];
     }

     return g;
}

void BeamAd::UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& RTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& pTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLocTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzTmp,
                         const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLocTmp)
{
     using namespace sp_grad;

     for (index_type i = 0; i < NUMSEZ; ++i) {
          R[i] = RTmp[i];
          p[i] = pTmp[i];
          g[i] = gTmp[i];
          L[i] = LTmp[i];
          DefLoc[i] = DefLocTmp[i];
          Az[i] = AzTmp[i];
          AzLoc[i] = AzLocTmp[i];
     }
}

template <typename T>
void
BeamAd::AssReactionForce(sp_grad::SpGradientAssVec<T>& WorkVec,
                         const std::array<sp_grad::SpColVectorA<T, 3>, NUMSEZ>& p,
                         const std::array<sp_grad::SpColVectorA<T, 6>, NUMSEZ>& Az,
                         const std::array<sp_grad::SpColVectorA<T, 3>, NUMNODES>& X) const
{
     using namespace sp_grad;

     const index_type iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
     const index_type iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
     const index_type iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();

     const SpColVector<T, 3> F_I = SubColVector<1, 1, 3>(Az[S_I]);

     WorkVec.AddItem(iNode1FirstMomIndex + 1, F_I);

     const SpColVector<T, 3> M_I = Cross(p[S_I] - X[NODE1], SubColVector<1, 1, 3>(Az[S_I])) + SubColVector<4, 1, 3>(Az[S_I]);

     WorkVec.AddItem(iNode1FirstMomIndex + 4, M_I);

     const SpColVector<T, 3> F_II = SubColVector<1, 1, 3>(Az[SII]) - SubColVector<1, 1, 3>(Az[S_I]);
     const SpColVector<T, 3> M_II = SubColVector<4, 1, 3>(Az[SII]) - SubColVector<4, 1, 3>(Az[S_I])
          + Cross(p[SII] - X[NODE2], SubColVector<1, 1, 3>(Az[SII]))
          - Cross(p[S_I] - X[NODE2], SubColVector<1, 1, 3>(Az[S_I]));

     WorkVec.AddItem(iNode2FirstMomIndex + 1, F_II);
     WorkVec.AddItem(iNode2FirstMomIndex + 4, M_II);

     const SpColVector<T, 3> F_III = -SubColVector<1, 1, 3>(Az[SII]);
     const SpColVector<T, 3> M_III = Cross(SubColVector<1, 1, 3>(Az[SII]), p[SII] - X[NODE3]) - SubColVector<4, 1, 3>(Az[SII]);

     WorkVec.AddItem(iNode3FirstMomIndex + 1, F_III);
     WorkVec.AddItem(iNode3FirstMomIndex + 4, M_III);
}

template <typename T>
inline void
BeamAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
               doublereal dCoef,
               const sp_grad::SpGradientVectorHandler<T>& XCurr,
               const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
               enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     DEBUGCOUT("dCoef=" << dCoef << "\n");

     std::array<SpColVectorA<T, 3>, NUMNODES> X, xTmp, gNod;
     SpMatrixA<T, 3, 3> RNod;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetgCurr(gNod[i], dCoef, func);
          pNode[i]->GetRCurr(RNod, dCoef, func);
          pNode[i]->GetXCurr(X[i], dCoef, func);

          xTmp[i] = X[i] + RNod * f[i];

          DEBUGCOUT("f[" << i << "]=" << f[i] << "\n");
          DEBUGCOUT("gNod[" << i << "]=" << gNod[i] << "\n");
          DEBUGCOUT("RNod[" << i << "]=" << RNod << "\n");
          DEBUGCOUT("X[" << i << "]=" << X[i] << "\n");
          DEBUGCOUT("xTmp[" << i << "]=" << xTmp[i] << "\n");
     }

     std::array<SpMatrixA<T, 3, 3>, NUMSEZ> RDelta, R;
     std::array<SpColVectorA<T, 3>, NUMSEZ> gGrad, p, g, L;
     std::array<SpColVectorA<T, 6>, NUMSEZ> DefLoc, Az, AzLoc;

     DEBUGCOUT("beam3(" << GetLabel() << ")\n");
     DEBUGCOUT("Beam::AssRes bFirstRes = " << bFirstRes << std::endl);

     /* Aggiorna le grandezze della trave nei punti di valutazione */
     for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

          /* Posizione */
          p[iSez] = InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

          /* Matrici di rotazione */
          g[iSez] = InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
          RDelta[iSez] = MatRVec(g[iSez]);
          R[iSez] = RDelta[iSez] * RRef[iSez];

          /* Derivate della posizione */
          L[iSez] = InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));

          /* Derivate dei parametri di rotazione */
          gGrad[iSez] = InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));

          /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> GgGrad = MatGVec(g[iSez]) * gGrad[iSez];

          for (index_type i = 1; i <= 3; ++i) {
               DefLoc[iSez](i) = Dot(R[iSez].GetCol(i), L[iSez]) - L0[iSez](i);
               DefLoc[iSez](i + 3) = Dot(R[iSez].GetCol(i), GgGrad) + DefLocRef[iSez](i + 3);
          }

          /* Calcola le azioni interne */
          AzLoc[iSez] = pD[iSez]->pGetConstLaw()->Update(DefLoc[iSez]);

          /* corregge le azioni interne locali (piezo, ecc) */
          AddInternalForces(AzLoc[iSez], iSez);

          /* Porta le azioni interne nel sistema globale */
          for (integer i = 1; i <= 3; ++i) {
               Az[iSez](i) = Dot(Transpose(R[iSez].GetRow(i)), SubColVector<1, 1, 3>(AzLoc[iSez]));
               Az[iSez](i + 3) = Dot(Transpose(R[iSez].GetRow(i)), SubColVector<4, 1, 3>(AzLoc[iSez]));
          }

          DEBUGCOUT("p[" << iSez << "]=" << p[iSez] << std::endl);
          DEBUGCOUT("g[" << iSez << "]=" << g[iSez] << std::endl);
          DEBUGCOUT("RDelta[" << iSez << "]=" << RDelta[iSez] << std::endl);
          DEBUGCOUT("RPrev[" << iSez << "]=" << RPrev[iSez] << std::endl);
          DEBUGCOUT("RRef[" << iSez << "]=" << RRef[iSez] << std::endl);
          DEBUGCOUT("R[" << iSez << "]=" << R[iSez] << std::endl);
          DEBUGCOUT("L[" << iSez << "]=" << L[iSez] << std::endl);
          DEBUGCOUT("DefLoc[" << iSez << "]=" << DefLoc[iSez] << std::endl);
          DEBUGCOUT("Az[" << iSez << "]=" << Az[iSez] << std::endl);
     }

     AssReactionForce(WorkVec, p, Az, X);

     UpdateState(R, p, g, L, DefLoc, Az, AzLoc);

     bFirstRes = false;
}

VariableSubMatrixHandler&
BeamAd::AssJac(VariableSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("Beam::AssJac => AssStiffnessMat");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
BeamAd::AssJac(VectorHandler& JacY,
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
BeamAd::AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("BeamAd::AssRes => AssStiffnessVec");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

ViscoElasticBeamAd::ViscoElasticBeamAd(unsigned int uL,
                                       const StructNodeAd* pN1,
                                       const StructNodeAd* pN2,
                                       const StructNodeAd* pN3,
                                       const Vec3& F1,
                                       const Vec3& F2,
                                       const Vec3& F3,
                                       const Mat3x3& R1,
                                       const Mat3x3& R2,
                                       const Mat3x3& R3,
                                       const Mat3x3& r_I,
                                       const Mat3x3& rII,
                                       const ConstitutiveLaw6D* pD_I,
                                       const ConstitutiveLaw6D* pDII,
                                       OrientationDescription ood,
                                       flag fOut)
:Elem(uL, fOut),
 Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut),
 ViscoElasticBeam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut),
 BeamAd(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut) 
{
}

ViscoElasticBeamAd::~ViscoElasticBeamAd()
{
}

inline void
ViscoElasticBeamAd::UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& RTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& pTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gPrimeTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& OmegaTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LPrimeTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLocTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefPrimeLocTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzTmp,
                                const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLocTmp)
{
     using namespace sp_grad;

     for (index_type i = 0; i < NUMSEZ; ++i) {
          R[i] = RTmp[i];
          p[i] = pTmp[i];
          g[i] = gTmp[i];
          gPrime[i] = gPrimeTmp[i];
          L[i] = LTmp[i];
          LPrime[i] = LPrimeTmp[i];
          DefLoc[i] = DefLocTmp[i];
          DefPrimeLoc[i] = DefPrimeLocTmp[i];
          Az[i] = AzTmp[i];
          AzLoc[i] = AzLocTmp[i];
     }
}

template <typename T>
inline void
ViscoElasticBeamAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                           doublereal dCoef,
                           const sp_grad::SpGradientVectorHandler<T>& XCurr,
                           const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                           enum sp_grad::SpFunctionCall func)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssRes");

     /* Riceve il vettore gia' dimensionato e con gli indici a posto
      * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */

     /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex
      * (ora e' superato) */

     using namespace sp_grad;

     std::array<SpColVectorA<T, 3>, NUMNODES> gNod, xTmp, gPrimeNod, xPrimeTmp;
     SpMatrixA<T, 3, 3> RNod;
     std::array<SpColVectorA<T, 3>, NUMNODES> XNod, XPrimeNod;
     SpColVectorA<T, 3> WNod;

     for (unsigned int i = 0; i < NUMNODES; i++) {
          pNode[i]->GetgCurr(gNod[i], dCoef, func);
          pNode[i]->GetRCurr(RNod, dCoef, func);
          pNode[i]->GetgPCurr(gPrimeNod[i], dCoef, func);
          pNode[i]->GetWCurr(WNod, dCoef, func);
          pNode[i]->GetXCurr(XNod[i], dCoef, func);
          pNode[i]->GetVCurr(XPrimeNod[i], dCoef, func);

          SpColVector<T, 3> fTmp = RNod * f[i];

          xTmp[i] = XNod[i] + fTmp;
          xPrimeTmp[i] = XPrimeNod[i] + Cross(WNod, fTmp);
     }

     std::array<SpMatrixA<T, 3, 3>, NUMSEZ> R, RDelta;
     std::array<SpColVectorA<T, 3>, NUMSEZ> p, g, gGrad, gPrime, gPrimeGrad, Omega, L, LPrime;
     std::array<SpColVectorA<T, 6>, NUMSEZ> DefLoc, DefPrimeLoc, Az, AzLoc;

     /* Aggiorna le grandezze della trave nei punti di valutazione */
     for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {

          /* Posizione */
          p[iSez] = InterpState(xTmp[NODE1],
                                xTmp[NODE2],
                                xTmp[NODE3], Beam::Section(iSez));

          /* Matrici di rotazione */
          g[iSez] = InterpState(gNod[NODE1],
                                gNod[NODE2],
                                gNod[NODE3], Beam::Section(iSez));

          const SpMatrix<T, 3, 3> G = MatGVec(g[iSez]);

          RDelta[iSez] = MatRVec(g[iSez]);
          R[iSez] = RDelta[iSez] * RRef[iSez];

          /* Velocita' angolare della sezione */
          gPrime[iSez] = InterpState(gPrimeNod[NODE1],
                                     gPrimeNod[NODE2],
                                     gPrimeNod[NODE3], Beam::Section(iSez));

          Omega[iSez] = G * gPrime[iSez]
               + RDelta[iSez] * OmegaRef[iSez];

          /* rate of MatG */
          const T dtmp0 = Dot(g[iSez], g[iSez]);
          const T dtmp1 = 4. + dtmp0;
          const T dtmp2 = -4. / (dtmp1 * dtmp1);
          const T dtmp3 = 2. / dtmp1;

          const SpColVector<T, 3> GPrimeg =(gPrime[iSez] * dtmp0 + g[iSez] * Dot(gPrime[iSez], g[iSez])) * dtmp2
               + Cross(gPrime[iSez], g[iSez]) * dtmp3;

          /* Derivate della posizione */
          L[iSez] = InterpDeriv(xTmp[NODE1],
                                xTmp[NODE2],
                                xTmp[NODE3], Beam::Section(iSez));

          /* Derivate della velocita' */
          LPrime[iSez] = InterpDeriv(xPrimeTmp[NODE1],
                                     xPrimeTmp[NODE2],
                                     xPrimeTmp[NODE3], Beam::Section(iSez));

          /* Derivate dei parametri di rotazione */
          gGrad[iSez] = InterpDeriv(gNod[NODE1],
                                    gNod[NODE2],
                                    gNod[NODE3], Beam::Section(iSez));

          /* Derivate delle derivate spaziali dei parametri di rotazione */
          gPrimeGrad[iSez] = InterpDeriv(gPrimeNod[NODE1],
                                         gPrimeNod[NODE2],
                                         gPrimeNod[NODE3], Beam::Section(iSez));

          /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> GgGrad = G * gGrad[iSez];

          for (index_type i = 1; i <= 3; ++i) {
               DefLoc[iSez](i) = Dot(R[iSez].GetCol(i), L[iSez]) - L0[iSez](i);
               DefLoc[iSez](i + 3) = Dot(R[iSez].GetCol(i), GgGrad) + DefLocRef[iSez](i + 3);
          }

          /* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
          const SpColVector<T, 3> DL1 = EvalUnique(LPrime[iSez] + Cross(L[iSez], Omega[iSez]));
          const SpColVector<T, 3> DL2 = EvalUnique(G * gPrimeGrad[iSez] + GPrimeg + Cross(GgGrad, Omega[iSez]));

          for (index_type i = 1; i <= 3; ++i) {
               DefPrimeLoc[iSez](i) = Dot(R[iSez].GetCol(i), DL1);
               DefPrimeLoc[iSez](i + 3) = Dot(R[iSez].GetCol(i), DL2) + DefPrimeLocRef[iSez](i + 3);
          }

          /* Calcola le azioni interne */
          AzLoc[iSez] = pD[iSez]->pGetConstLaw()->Update(DefLoc[iSez], DefPrimeLoc[iSez]);

          /* corregge le azioni interne locali (piezo, ecc) */
          AddInternalForces(AzLoc[iSez], iSez);

          /* Porta le azioni interne nel sistema globale */
          for (index_type i = 1; i <= 3; ++i) {
               Az[iSez](i) = Dot(Transpose(R[iSez].GetRow(i)), SubColVector<1, 1, 3>(AzLoc[iSez]));
               Az[iSez](i + 3) = Dot(Transpose(R[iSez].GetRow(i)), SubColVector<4, 1, 3>(AzLoc[iSez]));
          }
     }

     AssReactionForce(WorkVec, p, Az, XNod);

     UpdateState(R, p, g, gPrime, Omega, L, LPrime, DefLoc, DefPrimeLoc, Az, AzLoc);

     bFirstRes = false;
}

SubVectorHandler&
ViscoElasticBeamAd::AssRes(SubVectorHandler& WorkVec,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

VariableSubMatrixHandler&
ViscoElasticBeamAd::AssJac(VariableSubMatrixHandler& WorkMat,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
ViscoElasticBeamAd::AssJac(VectorHandler& JacY,
                           const VectorHandler& Y,
                           doublereal dCoef,
                           const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr,
                           VariableSubMatrixHandler& WorkMat)
{
     DEBUGCOUTFNAME("ViscoElasticBeamAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}
