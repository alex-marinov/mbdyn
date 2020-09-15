/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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
  Copyright (C) 2013(-2020) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#define CREATE_PROFILE 1

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

#include "module-autodiff_test.h"
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)

#include <clock_time.h>

#if defined(USE_AUTODIFF)
#include <gradient.h>
#include <matvec.h>
#include <matvecass.h>
#endif

#if defined(USE_SPARSE_AUTODIFF)
#include <sp_gradient.h>
#include <sp_matrix_base.h>
#include <sp_matvecass.h>
#endif

#include <ac/f2c.h>

extern "C" int __FC_DECL__(def_joint_ad_ass_res)(const doublereal X1[3],
						 const doublereal g1[3],
						 const doublereal XP1[3],
						 const doublereal gP1[3],
						 const doublereal R1_0[3*3],
						 const doublereal W1_0[3],
						 const doublereal X2[2],
						 const doublereal g2[3],
						 const doublereal XP2[3],
						 const doublereal gP2[3],
						 const doublereal R2_0[3*3],
						 const doublereal W2_0[3],
						 const doublereal S1[3*3],
						 const doublereal D1[3*3],
						 const doublereal o1[3],
						 const doublereal o2[3],
						 doublereal F1[3],
						 doublereal M1[3],
						 doublereal F2[3],
						 doublereal M2[3]);

extern "C" int __FC_DECL__(def_joint_ad_ass_res_dv)(const doublereal X1[3],
						    const doublereal X1d[3][12],
						    const doublereal g1[3],
						    const doublereal g1d[3][12],
						    const doublereal XP1[3],
						    const doublereal XP1d[3][12],
						    const doublereal gP1[3],
						    const doublereal gP1d[3][12],
						    const doublereal R1_0[3*3],
						    const doublereal W1_0[3],
						    const doublereal X2[3],
						    const doublereal X2d[3][12],
						    const doublereal g2[3],
						    const doublereal g2d[3][12],
						    const doublereal XP2[3],
						    const doublereal XP2d[3][12],
						    const doublereal gP2[3],
						    const doublereal gp2d[3][12],
						    const doublereal R2_0[3*3],
						    const doublereal W2_0[3*3],
						    const doublereal S1[3*3],
						    const doublereal D1[3*3],
						    const doublereal o1[3],
						    const doublereal o2[3],
						    doublereal F1[3],
						    doublereal F1d[3][12],
						    doublereal M1[3],
						    doublereal M1d[3][12],
						    doublereal F2[3],
						    doublereal F2d[3][12],
						    doublereal M2[3],
						    doublereal M2d[3][12],
						    const integer &nbdirs);

class DeformableJointAD: virtual public Elem, public UserDefinedElem
{
public:
     DeformableJointAD(unsigned uLabel, const DofOwner *pDO,
		       DataManager* pDM, MBDynParser& HP);
     virtual ~DeformableJointAD(void);
     virtual void Output(OutputHandler& OH) const;
     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
     VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
     SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
     VariableSubMatrixHandler&
     AssJacTrad(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
     SubVectorHandler&
     AssResTrad(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
     VariableSubMatrixHandler&
     AssJacF77(VariableSubMatrixHandler& WorkMat,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
     SubVectorHandler&
     AssResF77(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
#ifdef USE_AUTODIFF
     template <typename T>
     inline void
     AssRes(grad::GradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const grad::GradientVectorHandler<T>& XCurr,
	    const grad::GradientVectorHandler<T>& XPrimeCurr,
	    grad::FunctionCall func);
#endif

#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    sp_grad::SpFunctionCall func);
#endif
     
     unsigned int iGetNumPrivData(void) const;
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

private:
#ifdef USE_AUTODIFF
     static const grad::index_type iNumADVars = 12;
#endif
     const StructNode* pNode1;
     const StructNode* pNode2;
     Mat3x3 S1;
     Mat3x3 D1;
     Vec3 o1;
     Vec3 o2;
     
#ifdef USE_AUTODIFF
     grad::LocalDofMap dof;   
#endif
     
     enum AssemblyFlag {
	  TRADITIONAL = 0,
#ifdef USE_AUTODIFF
	  TEMPLATE_META_PROG,
#endif
#ifdef USE_SPARSE_AUTODIFF
	  SPARSE_TEMPLATE_META_PROG,
#endif
	  F77
     } fRes, fJac;

#if CREATE_PROFILE == 1
     enum {
	  RESIDUAL = 0,
	  JACOBIAN = 1
     };
     struct {
	  doublereal dtRes;
	  doublereal dtJac;
	  doublereal dtInit[2];
	  doublereal dtGet[2];
	  doublereal dtCalc[2];
	  doublereal dtAss[2];
     } profile;

     static const char* AssemblyFuncName(enum AssemblyFlag flag);
#endif
};

DeformableJointAD::DeformableJointAD(
     unsigned uLabel, const DofOwner *pDO,
     DataManager* pDM, MBDynParser& HP)
     :	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	pNode1(nullptr),
	pNode2(nullptr),
	fRes(TRADITIONAL),
	fJac(TRADITIONAL)
{
#if CREATE_PROFILE == 1
     memset(&profile, 0, sizeof(profile));
#endif

     // help
     if (HP.IsKeyWord("help"))
     {
	  silent_cout(
	       "\n"
	       "Module:		DeformableJointAD\n"
	       "\n"
	       "	deformable joint ad,\n"
	       "		node1, (label) <node1>,\n"
	       "		[offset, (Vec3) <o1>,]\n"
	       "		node2, (label) <node2>,\n"
	       "		[offset, (Vec3) <o2>,]\n"
	       "		stiffness, (Mat3x3) S1,\n"
	       "		damping, (Mat3x3) D1,\n"
	       << std::endl);

	  if (!HP.IsArg())
	  {
	       /*
		* Exit quietly if nothing else is provided
		*/
	       throw NoErr(MBDYN_EXCEPT_ARGS);
	  }
     }

     if (!HP.IsKeyWord("node1"))
     {
	  silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

     if (!pNode1) {
	  silent_cerr("deformable joint ad" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if (HP.IsKeyWord("offset")) {
	  o1 = HP.GetPosRel(ReferenceFrame(pNode1));
     }

     if (!HP.IsKeyWord("node2"))
     {
	  silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

     if (!pNode2) {
	  silent_cerr("deformable joint ad" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if ( HP.IsKeyWord("offset") ) {
	  o2 = HP.GetPosRel(ReferenceFrame(pNode2));
     }

     if (!HP.IsKeyWord("stiffness")) {
	  silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"stiffness\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     S1 = HP.GetMat3x3();

     if (!HP.IsKeyWord("damping")) {
	  silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"damping\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     D1 = HP.GetMat3x3();

     if (HP.IsKeyWord("residual")) {
	  if (HP.IsKeyWord("traditional")) {
	       fRes = TRADITIONAL;
#ifdef USE_AUTODIFF
	  } else if (HP.IsKeyWord("template" "meta" "program")) {
	       fRes = TEMPLATE_META_PROG;
#endif
#ifdef USE_SPARSE_AUTODIFF
	  } else if (HP.IsKeyWord("sparse" "template" "meta" "program")) {
	       fRes = SPARSE_TEMPLATE_META_PROG;
#endif
	  } else if (HP.IsKeyWord("fortran77")) {
	       fRes = F77;
	  } else {
	       silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"traditional\""
#ifdef USE_AUTODIFF
			   ", \"template meta program\", "
#endif
#ifdef USE_SPARSE_AUTODIFF
			   ", \"sparse template meta program\""
#endif
			   " or \"fortran77\" expected at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
     }

     if (HP.IsKeyWord("jacobian")) {
	  if (HP.IsKeyWord("traditional")) {
	       fJac = TRADITIONAL;
#ifdef USE_AUTODIFF
	  } else if (HP.IsKeyWord("template" "meta" "program")) {
	       fJac = TEMPLATE_META_PROG;
#endif
#ifdef USE_SPARSE_AUTODIFF
	  } else if (HP.IsKeyWord("sparse" "template" "meta" "program")) {
	       fJac = SPARSE_TEMPLATE_META_PROG;
#endif
	  } else if (HP.IsKeyWord("fortran77")) {
	       fJac = F77;
	  } else {
	       silent_cerr("deformable joint ad" << GetLabel() << "): keyword \"traditional\" or \"template meta program\" expected at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
     }

     SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

DeformableJointAD::~DeformableJointAD(void)
{
#if CREATE_PROFILE == 1
     std::cerr << "deformable joint ad(" << GetLabel() << ")" << std::endl;

     std::cerr << "dtRes["
	       << AssemblyFuncName(fRes)
	       << "]=" << profile.dtRes << std::endl;

     std::cerr << "dtJac["
	       << AssemblyFuncName(fJac)
	       << "]=" << profile.dtJac << std::endl;


     std::cerr << "dtInit[RESIDUAL]=" << profile.dtInit[RESIDUAL] << std::endl;
     std::cerr << "dtGet[RESIDUAL]=" << profile.dtGet[RESIDUAL] << std::endl;
     std::cerr << "dtCalc[RESIDUAL]=" << profile.dtCalc[RESIDUAL] << std::endl;
     std::cerr << "dtAss[RESIDUAL]=" << profile.dtAss[RESIDUAL] << std::endl;

     std::cerr << "dtInit[JACOBIAN]=" << profile.dtInit[JACOBIAN] << std::endl;
     std::cerr << "dtGet[JACOBIAN]=" << profile.dtGet[JACOBIAN] << std::endl;
     std::cerr << "dtCalc[JACOBIAN]=" << profile.dtCalc[JACOBIAN] << std::endl;
     std::cerr << "dtAss[JACOBIAN]=" << profile.dtAss[JACOBIAN] << std::endl;

     std::cerr << std::endl;
#endif
}

void
DeformableJointAD::Output(OutputHandler& OH) const
{

}

void
DeformableJointAD::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     switch (fJac) {
     case SPARSE_TEMPLATE_META_PROG:
	  *piNumRows = -12;
	  *piNumCols = 24;
	  break;
     default:
	  *piNumRows = 12;
	  *piNumCols = 12;
     }
}

VariableSubMatrixHandler&
DeformableJointAD::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     const double start = mbdyn_clock_time();
#endif

     switch (fJac) {
     case TRADITIONAL:
	  AssJacTrad(WorkMat, dCoef, XCurr, XPrimeCurr);
	  break;
#ifdef USE_AUTODIFF
     case TEMPLATE_META_PROG:
	  grad::GradientAssVec<grad::Gradient<iNumADVars> >::AssJac(this,
								    WorkMat.SetSparse(),
								    dCoef,
								    XCurr,
								    XPrimeCurr,
								    grad::REGULAR_JAC,
								    &dof);
	  break;
#endif
#ifdef USE_SPARSE_AUTODIFF
     case SPARSE_TEMPLATE_META_PROG:
	  sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
								 WorkMat.SetSparse(),
								 dCoef,
								 XCurr,
								 XPrimeCurr,
								 sp_grad::REGULAR_JAC);
	  break;
#endif
     case F77:
	  AssJacF77(WorkMat, dCoef, XCurr, XPrimeCurr);
	  break;

     default:
	  WorkMat.SetNullMatrix();
     }

#if CREATE_PROFILE == 1
     profile.dtJac += mbdyn_clock_time() - start;
#endif

     return WorkMat;
}

SubVectorHandler&
DeformableJointAD::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     const double start = mbdyn_clock_time();
#endif

     switch (fRes) {
     case TRADITIONAL:
	  AssResTrad(WorkVec, dCoef, XCurr, XPrimeCurr);
	  break;
#ifdef USE_AUTODIFF
     case TEMPLATE_META_PROG:
	  grad::GradientAssVec<doublereal>::AssRes(this,
						   WorkVec,
						   dCoef,
						   XCurr,
						   XPrimeCurr,
						   grad::REGULAR_RES);
	  break;
#endif
#ifdef USE_SPARSE_AUTODIFF
     case SPARSE_TEMPLATE_META_PROG:
	  sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
							WorkVec,
							dCoef,
							XCurr,
							XPrimeCurr,
							sp_grad::REGULAR_RES);
	  break;
#endif
     case F77:
	  AssResF77(WorkVec, dCoef, XCurr, XPrimeCurr);
	  break;

     default:
	  WorkVec.ResizeReset(0);
     }

#if CREATE_PROFILE == 1
     profile.dtRes += mbdyn_clock_time() - start;
#endif

     return WorkVec;
}

VariableSubMatrixHandler&
DeformableJointAD::AssJacTrad(VariableSubMatrixHandler& WorkMatVar,
			      doublereal dCoef,
			      const VectorHandler& XCurr,
			      const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     const Vec3& X1 = pNode1->GetXCurr();
     const Vec3& V1 = pNode1->GetVCurr();
     const Mat3x3& R1 = pNode1->GetRCurr();
     const Mat3x3& R1_0 = pNode1->GetRRef();
     const Vec3& W1 = pNode1->GetWCurr();
     const Vec3& W1_0 = pNode1->GetWRef();

     const Vec3& X2 = pNode2->GetXCurr();
     const Vec3& V2 = pNode2->GetVCurr();
     const Mat3x3& R2 = pNode2->GetRCurr();
     const Mat3x3& R2_0 = pNode2->GetRRef();
     const Vec3& W2 = pNode2->GetWCurr();
     const Vec3& W2_0 = pNode2->GetWRef();

#if CREATE_PROFILE == 1
     profile.dtGet[JACOBIAN] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const Vec3 R1o1 = R1 * o1;
     const Mat3x3 skew_R1o1(MatCross, R1o1);
     const Vec3 R1_0o1 = R1_0 * o1;
     const Mat3x3 skew_R1_0o1(MatCross, R1_0o1);
     const Vec3 R2o2 = R2 * o2;
     const Mat3x3 skew_R2o2(MatCross, R2o2);
     const Vec3 R2_0o2 = R2_0 * o2;
     const Mat3x3 skew_R2_0o2(MatCross, R2_0o2);
     const Vec3 dX = R1.MulTV(X1 + R1o1 - X2 - R2o2);
     const Vec3 dV = R1.MulTV(V1 + W1.Cross(R1o1) - V2 - W2.Cross(R2o2));
     const Vec3 F1_R1 = -(S1 * dX + D1 * dV);
     const Vec3 F1 = R1 * F1_R1;
     const Vec3 F2 = -F1;

     const Mat3x3 dF1_dX1 = -R1 * S1.MulMT(R1);

     const Mat3x3 ddX_dg1 = R1_0.MulTM(Mat3x3(MatCross, X1 + R1o1 - X2 - R2o2)) - R1.MulTM(skew_R1_0o1);
     const Mat3x3 ddV_dg1 = R1_0.MulTM(Mat3x3(MatCross, V1 + W1.Cross(R1o1) - V2 - W2.Cross(R2o2)))
	  + R1.MulTM(Mat3x3(MatCrossCross, R1o1, W1_0) - Mat3x3(MatCrossCross, W1, R1_0o1));

     const Mat3x3 dF1_dg1 = Mat3x3(MatCross, R1_0 * (-F1_R1)) - R1 * (S1 * ddX_dg1 + D1 * ddV_dg1);

     const Mat3x3 dF1_dX2 = R1 * S1.MulMT(R1);
     const Mat3x3 ddX_dg2 = R1.MulTM(skew_R2_0o2);
     const Mat3x3 ddV_dg2 = R1.MulTM(Mat3x3(MatCrossCross, R2o2, -W2_0) + Mat3x3(MatCrossCross, W2, R2_0o2));
     const Mat3x3 dF1_dg2 = -R1 * (S1 * ddX_dg2 + D1 * ddV_dg2);

     const Mat3x3 dF2_dX1 = -dF1_dX1;
     const Mat3x3 dF2_dg1 = -dF1_dg1;
     const Mat3x3 dF2_dX2 = -dF1_dX2;
     const Mat3x3 dF2_dg2 = -dF1_dg2;

     const Mat3x3 dM1_dX1 = skew_R1o1 * dF1_dX1;
     const Mat3x3 dM1_dg1 = Mat3x3(MatCrossCross, F1, R1_0o1) + skew_R1o1 * dF1_dg1;
     const Mat3x3 dM1_dX2 = skew_R1o1 * dF1_dX2;
     const Mat3x3 dM1_dg2 = skew_R1o1 * dF1_dg2;

     const Mat3x3 dM2_dX1 = skew_R2o2 * dF2_dX1;
     const Mat3x3 dM2_dg1 = skew_R2o2 * dF2_dg1;
     const Mat3x3 dM2_dX2 = skew_R2o2 * dF2_dX2;
     const Mat3x3 dM2_dg2 = Mat3x3(MatCrossCross, F2, R2_0o2) + skew_R2o2 * dF2_dg2;

     const Mat3x3 dF1_dV1 = -R1 * D1.MulMT(R1);
     const Mat3x3 ddV_dgP1 = -R1.MulTM(skew_R1o1);
     const Mat3x3 dF1_dgP1 = -R1 * D1 * ddV_dgP1;

     const Mat3x3 dF1_dV2 = R1 * D1.MulMT(R1);
     const Mat3x3 ddV_dgP2 = R1.MulTM(skew_R2o2);
     const Mat3x3 dF1_dgP2 = -R1 * D1 * ddV_dgP2;

     const Mat3x3 dM1_dV1 = skew_R1o1 * dF1_dV1;
     const Mat3x3 dM1_dgP1 = skew_R1o1 * dF1_dgP1;
     const Mat3x3 dM1_dV2 = skew_R1o1 * dF1_dV2;
     const Mat3x3 dM1_dgP2 = skew_R1o1 * dF1_dgP2;

     const Mat3x3 dF2_dV1 = -dF1_dV1;
     const Mat3x3 dF2_dgP1 = -dF1_dgP1;
     const Mat3x3 dF2_dV2 = -dF1_dV2;
     const Mat3x3 dF2_dgP2 = -dF1_dgP2;

     const Mat3x3 dM2_dV1 = skew_R2o2 * dF2_dV1;
     const Mat3x3 dM2_dgP1 = skew_R2o2 * dF2_dgP1;
     const Mat3x3 dM2_dV2 = skew_R2o2 * dF2_dV2;
     const Mat3x3 dM2_dgP2 = skew_R2o2 * dF2_dgP2;

#if CREATE_PROFILE == 1
     profile.dtCalc[JACOBIAN] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

     integer iNumRows, iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);

     WorkMat.ResizeReset(iNumRows, iNumCols);

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();
     const integer iFirstPosIndexNode1 = pNode1->iGetFirstPositionIndex();
     const integer iFirstPosIndexNode2 = pNode2->iGetFirstPositionIndex();

     for (integer i = 1; i <= 6; ++i) {
	  WorkMat.PutRowIndex(i,     iFirstMomIndexNode1 + i);
	  WorkMat.PutRowIndex(i + 6, iFirstMomIndexNode2 + i);
	  WorkMat.PutColIndex(i,     iFirstPosIndexNode1 + i);
	  WorkMat.PutColIndex(i + 6, iFirstPosIndexNode2 + i);
     }

     WorkMat.Put(1,  1, -dF1_dV1  - dF1_dX1 * dCoef);
     WorkMat.Put(1,  4, -dF1_dgP1 - dF1_dg1 * dCoef);
     WorkMat.Put(1,  7, -dF1_dV2  - dF1_dX2 * dCoef);
     WorkMat.Put(1, 10, -dF1_dgP2 - dF1_dg2 * dCoef);

     WorkMat.Put(4,  1, -dM1_dV1  - dM1_dX1 * dCoef);
     WorkMat.Put(4,  4, -dM1_dgP1 - dM1_dg1 * dCoef);
     WorkMat.Put(4,  7, -dM1_dV2  - dM1_dX2 * dCoef);
     WorkMat.Put(4, 10, -dM1_dgP2 - dM1_dg2 * dCoef);

     WorkMat.Put(7,  1, -dF2_dV1  - dF2_dX1 * dCoef);
     WorkMat.Put(7,  4, -dF2_dgP1 - dF2_dg1 * dCoef);
     WorkMat.Put(7,  7, -dF2_dV2  - dF2_dX2 * dCoef);
     WorkMat.Put(7, 10, -dF2_dgP2 - dF2_dg2 * dCoef);

     WorkMat.Put(10,  1, -dM2_dV1  - dM2_dX1 * dCoef);
     WorkMat.Put(10,  4, -dM2_dgP1 - dM2_dg1 * dCoef);
     WorkMat.Put(10,  7, -dM2_dV2  - dM2_dX2 * dCoef);
     WorkMat.Put(10, 10, -dM2_dgP2 - dM2_dg2 * dCoef);

#if CREATE_PROFILE == 1
     profile.dtAss[JACOBIAN] += mbdyn_clock_time() - start;
#endif

     return WorkMatVar;
}

SubVectorHandler&
DeformableJointAD::AssResTrad(SubVectorHandler& WorkVec,
			      doublereal dCoef,
			      const VectorHandler& XCurr,
			      const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     const Vec3& X1 = pNode1->GetXCurr();
     const Vec3& XP1 = pNode1->GetVCurr();
     const Mat3x3& R1 = pNode1->GetRCurr();
     const Vec3& W1 = pNode1->GetWCurr();

     const Vec3& X2 = pNode2->GetXCurr();
     const Vec3& XP2 = pNode2->GetVCurr();
     const Mat3x3& R2 = pNode2->GetRCurr();
     const Vec3& W2 = pNode2->GetWCurr();

#if CREATE_PROFILE == 1
     profile.dtGet[RESIDUAL] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const Vec3 R1o1 = R1 * o1;
     const Vec3 R2o2 = R2 * o2;
     const Vec3 dX = R1.MulTV(X1 + R1o1 - X2 - R2o2);
     const Vec3 dXP = R1.MulTV(XP1 + W1.Cross(R1o1) - XP2 - W2.Cross(R2o2));
     const Vec3 F1 = -(R1 * (S1 * dX + D1 * dXP));
     const Vec3 F2 = -F1;
     const Vec3 M1 = R1o1.Cross(F1);
     const Vec3 M2 = R2o2.Cross(F2);

#if CREATE_PROFILE == 1
     profile.dtCalc[RESIDUAL] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();

     integer iNumRows, iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);

     WorkVec.ResizeReset(iNumRows);

     for (integer i = 1; i <= 6; ++i) {
	  WorkVec.PutRowIndex(i,     iFirstMomIndexNode1 + i);
	  WorkVec.PutRowIndex(i + 6, iFirstMomIndexNode2 + i);
     }

     WorkVec.Put(1, F1);
     WorkVec.Put(4, M1);
     WorkVec.Put(7, F2);
     WorkVec.Put(10, M2);

#if CREATE_PROFILE == 1
     profile.dtAss[RESIDUAL] += mbdyn_clock_time() - start;
#endif

     return WorkVec;
}

VariableSubMatrixHandler&
DeformableJointAD::AssJacF77(VariableSubMatrixHandler& WorkMatVar,
			     doublereal dCoef,
			     const VectorHandler& XCurr,
			     const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     const doublereal* X1 = pNode1->GetXCurr().pGetVec();
     const doublereal* g1 = pNode1->GetgCurr().pGetVec();
     const doublereal* XP1 = pNode1->GetVCurr().pGetVec();
     const doublereal* gP1 = pNode1->GetgPCurr().pGetVec();
     const doublereal* R1_0 = pNode1->GetRRef().pGetMat();
     const doublereal* W1_0 = pNode1->GetWRef().pGetVec();
     const doublereal* X2 = pNode2->GetXCurr().pGetVec();
     const doublereal* g2 = pNode2->GetgCurr().pGetVec();
     const doublereal* XP2 = pNode2->GetVCurr().pGetVec();
     const doublereal* gP2 = pNode2->GetgPCurr().pGetVec();
     const doublereal* R2_0 = pNode2->GetRRef().pGetMat();
     const doublereal* W2_0 = pNode2->GetWRef().pGetVec();

     doublereal X1d[3][12] = {{0.}}, g1d[3][12] = {{0.}};
     doublereal XP1d[3][12] = {{0.}}, gP1d[3][12] = {{0.}};
     doublereal X2d[3][12] = {{0.}}, g2d[3][12] = {{0.}};
     doublereal XP2d[3][12] = {{0.}}, gP2d[3][12] = {{0.}};

     for (integer i = 0; i < 3; ++i) {
	  X1d[i][i] = -dCoef;
	  XP1d[i][i] = -1.;
	  g1d[i][i + 3] = -dCoef;
	  gP1d[i][i + 3] = -1.;
	  X2d[i][i + 6] = -dCoef;
	  XP2d[i][i + 6] = -1.;
	  g2d[i][i + 9] = -dCoef;
	  gP2d[i][i + 9] = -1.;
     }

#if CREATE_PROFILE == 1
     profile.dtGet[JACOBIAN] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     doublereal F1[3], M1[3], F2[3], M2[3];
     doublereal F1d[3][12], M1d[3][12], F2d[3][12], M2d[3][12];

     __FC_DECL__(def_joint_ad_ass_res_dv)(X1,
					  X1d,
					  g1,
					  g1d,
					  XP1,
					  XP1d,
					  gP1,
					  gP1d,
					  R1_0,
					  W1_0,
					  X2,
					  X2d,
					  g2,
					  g2d,
					  XP2,
					  XP2d,
					  gP2,
					  gP2d,
					  R2_0,
					  W2_0,
					  S1.pGetMat(),
					  D1.pGetMat(),
					  o1.pGetVec(),
					  o2.pGetVec(),
					  F1,
					  F1d,
					  M1,
					  M1d,
					  F2,
					  F2d,
					  M2,
					  M2d,
					  12);

#if CREATE_PROFILE == 1
     profile.dtCalc[JACOBIAN] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

     integer iNumRows, iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);

     WorkMat.ResizeReset(iNumRows, iNumCols);

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();
     const integer iFirstPosIndexNode1 = pNode1->iGetFirstPositionIndex();
     const integer iFirstPosIndexNode2 = pNode2->iGetFirstPositionIndex();

     for (integer i = 1; i <= 6; ++i) {
	  WorkMat.PutRowIndex(i,     iFirstMomIndexNode1 + i);
	  WorkMat.PutRowIndex(i + 6, iFirstMomIndexNode2 + i);
	  WorkMat.PutColIndex(i,     iFirstPosIndexNode1 + i);
	  WorkMat.PutColIndex(i + 6, iFirstPosIndexNode2 + i);
     }

     for (integer i = 0; i < 3; ++i) {
	  for (integer j = 0; j < 12; ++j) {
	       WorkMat.PutCoef(i + 1, j + 1, F1d[i][j]);
	       WorkMat.PutCoef(i + 4, j + 1, M1d[i][j]);
	       WorkMat.PutCoef(i + 7, j + 1, F2d[i][j]);
	       WorkMat.PutCoef(i + 10, j + 1, M2d[i][j]);
	  }
     }

#if CREATE_PROFILE == 1
     profile.dtAss[JACOBIAN] += mbdyn_clock_time() - start;
#endif

     return WorkMatVar;
}

SubVectorHandler&
DeformableJointAD::AssResF77(SubVectorHandler& WorkVec,
			     doublereal dCoef,
			     const VectorHandler& XCurr,
			     const VectorHandler& XPrimeCurr)
{
#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     const doublereal* X1 = pNode1->GetXCurr().pGetVec();
     const doublereal* g1 = pNode1->GetgCurr().pGetVec();
     const doublereal* XP1 = pNode1->GetVCurr().pGetVec();
     const doublereal* gP1 = pNode1->GetgPCurr().pGetVec();
     const doublereal* R1_0 = pNode1->GetRRef().pGetMat();
     const doublereal* W1_0 = pNode1->GetWRef().pGetVec();
     const doublereal* X2 = pNode2->GetXCurr().pGetVec();
     const doublereal* g2 = pNode2->GetgCurr().pGetVec();
     const doublereal* XP2 = pNode2->GetVCurr().pGetVec();
     const doublereal* gP2 = pNode2->GetgPCurr().pGetVec();
     const doublereal* R2_0 = pNode2->GetRRef().pGetMat();
     const doublereal* W2_0 = pNode2->GetWRef().pGetVec();

#if CREATE_PROFILE == 1
     profile.dtGet[RESIDUAL] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     doublereal F1[3], M1[3], F2[3], M2[3];

     __FC_DECL__(def_joint_ad_ass_res)(X1,
				       g1,
				       XP1,
				       gP1,
				       R1_0,
				       W1_0,
				       X2,
				       g2,
				       XP2,
				       gP2,
				       R2_0,
				       W2_0,
				       S1.pGetMat(),
				       D1.pGetMat(),
				       o1.pGetVec(),
				       o2.pGetVec(),
				       F1,
				       M1,
				       F2,
				       M2);

#if CREATE_PROFILE == 1
     profile.dtCalc[RESIDUAL] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();

     integer iNumRows, iNumCols;

     WorkSpaceDim(&iNumRows, &iNumCols);

     WorkVec.ResizeReset(iNumRows);

     for (integer i = 1; i <= 6; ++i) {
	  WorkVec.PutRowIndex(i,     iFirstMomIndexNode1 + i);
	  WorkVec.PutRowIndex(i + 6, iFirstMomIndexNode2 + i);
     }

     for (integer i = 0; i < 3; ++i) {
	  WorkVec.PutCoef(1 + i, F1[i]);
	  WorkVec.PutCoef(4 + i, M1[i]);
	  WorkVec.PutCoef(7 + i, F2[i]);
	  WorkVec.PutCoef(10 + i, M2[i]);
     }

#if CREATE_PROFILE == 1
     profile.dtAss[RESIDUAL] += mbdyn_clock_time() - start;
#endif

     return WorkVec;
}

#ifdef USE_AUTODIFF
template <typename T>
inline void DeformableJointAD::AssRes(grad::GradientAssVec<T>& WorkVec,
				      doublereal dCoef,
				      const grad::GradientVectorHandler<T>& XCurr,
				      const grad::GradientVectorHandler<T>& XPrimeCurr,
				      enum grad::FunctionCall func)
{
     using namespace grad;
     
     typedef Matrix<T, 3, 3> Mat3x3;
     typedef Vector<T, 3> Vec3;

#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     Mat3x3 R1, R2;
     Vec3 X1, X2, XP1, XP2, W1, W2;

#if CREATE_PROFILE == 1
     const index_type iFunc = (func == grad::REGULAR_RES) ? RESIDUAL : JACOBIAN;
     profile.dtInit[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     pNode1->GetXCurr(X1, dCoef, func, &dof);
     pNode1->GetVCurr(XP1, dCoef, func, &dof);
     pNode1->GetRCurr(R1, dCoef, func, &dof);
     pNode1->GetWCurr(W1, dCoef, func, &dof);

     pNode2->GetXCurr(X2, dCoef, func, &dof);
     pNode2->GetVCurr(XP2, dCoef, func, &dof);
     pNode2->GetRCurr(R2, dCoef, func, &dof);
     pNode2->GetWCurr(W2, dCoef, func, &dof);

#if CREATE_PROFILE == 1
     profile.dtGet[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const Vec3 R1o1 = R1 * o1;
     const Vec3 R2o2 = R2 * o2;
     const Vec3 dX = Transpose(R1) * Vec3(X1 + R1o1 - X2 - R2o2);

     const Vec3 dXP = Transpose(R1) * Vec3(XP1 + Cross(W1, R1o1) - XP2 - Cross(W2, R2o2));
     const Vec3 F1 = -(R1 * Vec3(S1 * dX + D1 * dXP));

     const Vec3 F2 = -F1;
     const Vec3 M1 = Cross(R1o1, F1);
     const Vec3 M2 = Cross(R2o2, F2);

#if CREATE_PROFILE == 1
     profile.dtCalc[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();

     WorkVec.AddItem(iFirstMomIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstMomIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstMomIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstMomIndexNode2 + 4, M2);

#if CREATE_PROFILE == 1
     profile.dtAss[iFunc] += mbdyn_clock_time() - start;
#endif
}
#endif

#ifdef USE_SPARSE_AUTODIFF
template <typename T>
inline void DeformableJointAD::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
				      doublereal dCoef,
				      const sp_grad::SpGradientVectorHandler<T>& XCurr,
				      const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
				      sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;
     typedef SpMatrix<T, 3, 3> Mat3x3;
     typedef SpColVector<T, 3> Vec3;

#if CREATE_PROFILE == 1
     doublereal start = mbdyn_clock_time();
#endif

     Mat3x3 R1, R2;
     Vec3 X1, X2, XP1, XP2, W1, W2;
     
#if CREATE_PROFILE == 1
     const index_type iFunc = (func == REGULAR_RES) ? RESIDUAL : JACOBIAN;
     profile.dtInit[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     pNode1->GetXCurr(X1, dCoef, func);
     pNode1->GetVCurr(XP1, dCoef, func);
     pNode1->GetRCurr(R1, dCoef, func);
     pNode1->GetWCurr(W1, dCoef, func);

     pNode2->GetXCurr(X2, dCoef, func);
     pNode2->GetVCurr(XP2, dCoef, func);
     pNode2->GetRCurr(R2, dCoef, func);
     pNode2->GetWCurr(W2, dCoef, func);

#if CREATE_PROFILE == 1
     profile.dtGet[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const Vec3 R1o1 = R1 * o1;
     const Vec3 R2o2 = R2 * o2;
     const Vec3 dX = Transpose(R1) * (X1 + R1o1 - X2 - R2o2);

     const Vec3 dXP = Transpose(R1) * (XP1 + Cross(W1, R1o1) - XP2 - Cross(W2, R2o2));
     const Vec3 F1 = -(R1 * (S1 * dX + D1 * dXP));

     const Vec3 F2 = -F1;
     const Vec3 M1 = Cross(R1o1, F1);
     const Vec3 M2 = Cross(R2o2, F2);

#if CREATE_PROFILE == 1
     profile.dtCalc[iFunc] += mbdyn_clock_time() - start;
     start = mbdyn_clock_time();
#endif

     const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomIndexNode2 = pNode2->iGetFirstMomentumIndex();

     WorkVec.AddItem(iFirstMomIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstMomIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstMomIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstMomIndexNode2 + 4, M2);

#if CREATE_PROFILE == 1
     profile.dtAss[iFunc] += mbdyn_clock_time() - start;
#endif
}
#endif

unsigned int
DeformableJointAD::iGetNumPrivData(void) const
{
     return 0;
}

int
DeformableJointAD::iGetNumConnectedNodes(void) const
{
     return 2; // 1x shaft + 1x bearing
}

void
DeformableJointAD::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(iGetNumConnectedNodes());
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
}

void
DeformableJointAD::SetValue(DataManager *pDM,
			    VectorHandler& X, VectorHandler& XP,
			    SimulationEntity::Hints *ph)
{

}

std::ostream&
DeformableJointAD::Restart(std::ostream& out) const
{
     return out;
}

unsigned int
DeformableJointAD::iGetInitialNumDof(void) const
{
     return 0;
}

void
DeformableJointAD::InitialWorkSpaceDim(
     integer* piNumRows,
     integer* piNumCols) const
{
     *piNumRows = 0;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
DeformableJointAD::InitialAssJac(
     VariableSubMatrixHandler& WorkMat,
     const VectorHandler& XCurr)
{
     WorkMat.SetNullMatrix();

     return WorkMat;
}

SubVectorHandler&
DeformableJointAD::InitialAssRes(
     SubVectorHandler& WorkVec,
     const VectorHandler& XCurr)
{
     WorkVec.ResizeReset(0);

     return WorkVec;
}

#if CREATE_PROFILE == 1
const char* DeformableJointAD::AssemblyFuncName(enum AssemblyFlag flag)
{
     switch (flag) {
     case TRADITIONAL:
	  return "traditional";
#ifdef USE_AUTODIFF
     case TEMPLATE_META_PROG:
	  return "template meta program";
#endif
#ifdef USE_SPARSE_AUTODIFF
     case SPARSE_TEMPLATE_META_PROG:
	  return "sparse template meta prog";
#endif
     case F77:
	  return "Fortran 77";

     default:
	  return "?";
     }
}
#endif

class InlineJointAD: virtual public Elem, public UserDefinedElem
{
public:
     InlineJointAD(unsigned uLabel, const DofOwner *pDO,
		   DataManager* pDM, MBDynParser& HP);
     virtual ~InlineJointAD(void);
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
     VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
     SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     template <typename T>
     inline void
     AssRes(grad::GradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const grad::GradientVectorHandler<T>& XCurr,
	    const grad::GradientVectorHandler<T>& XPrimeCurr,
	    enum grad::FunctionCall func);
#endif
#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    sp_grad::SpFunctionCall func);
#endif
     int iGetNumConnectedNodes(void) const;
     void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
     void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph);
     virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
     virtual void Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
     std::ostream& Restart(std::ostream& out) const;
     virtual unsigned int iGetInitialNumDof(void) const;
     virtual void
     InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
     VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);
     SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     template <typename T>
     inline void
     InitialAssRes(grad::GradientAssVec<T>& WorkVec,
		   const grad::GradientVectorHandler<T>& XCurr,
		   enum grad::FunctionCall func);
#endif
#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     inline void
     InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
		   const sp_grad::SpGradientVectorHandler<T>& XCurr,
		   sp_grad::SpFunctionCall func);
#endif
private:
     StructNode* pNode1;
     StructNode* pNode2;

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     grad::Vector<doublereal, 3> o1;
     grad::Matrix<doublereal, 3, 3> e;
     grad::Vector<doublereal, 3> o2;
#endif
#ifdef USE_SPARSE_AUTODIFF     
     sp_grad::SpColVector<doublereal, 3> o1;
     sp_grad::SpMatrix<doublereal, 3, 3> e;
     sp_grad::SpColVector<doublereal, 3> o2;
#endif
     enum LagrangeMultiplierIndex
     {
	  L1 = 0,
	  L2 = 1
     };

     doublereal lambda[2];

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     LocalDofMap dof;
#endif
     static constexpr integer iWorkSpace = 14;
     static constexpr integer iInitialWorkSpace = 28;
};

InlineJointAD::InlineJointAD(
     unsigned uLabel, const DofOwner *pDO,
     DataManager* pDM, MBDynParser& HP)
     :	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	pNode1(nullptr),
	pNode2(nullptr),
	o1(Zero3),
	e(Eye3),
	o2(Zero3)
{
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     using namespace grad;
#else
     using namespace sp_grad;
#endif
     
     memset(lambda, 0, sizeof(lambda));

     // help
     if (HP.IsKeyWord("help")) {
	  silent_cout(
	       "\n"
	       "Module:		InlineAD\n"
	       "\n"
	       "	This element implements a inline joint with friction\n"
	       "\n"
	       "	inline friction,\n"
	       "		node1, (label) <node1>,\n"
	       "			[ offset, (Vec3) <offset>, ]\n"
	       "			[ hinge, (Mat3x3) <orientation>, ]\n"
	       "		node2, (label) <node2>,\n"
	       "			[ offset, (Vec3) <offset>, ]\n"
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
	  silent_cerr("inline friction(" << GetLabel() << "): keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pNode1 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

     if (!pNode1) {
	  silent_cerr("inline friction ad(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
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
	  silent_cerr("inline friction ad(" << GetLabel() << "): keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pNode2 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

     if (!pNode2) {
	  silent_cerr("inline friction ad(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if (HP.IsKeyWord("offset")) {
	  const ReferenceFrame refNode2(pNode2);

	  o2 = HP.GetPosRel(refNode2);
     }

     SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

     std::ostream& out = pDM->GetLogFile();

     out << "inline joint ad: " << GetLabel() << " "
	 << pNode1->GetLabel() << " "
	 << o1 << " "
	 << e << " "
	 << pNode2->GetLabel() << " " << o2 << " "
	 << std::endl;
}

InlineJointAD::~InlineJointAD(void)
{
     // destroy private data
}

unsigned int InlineJointAD::iGetNumDof(void) const
{
     return 2u;
}

DofOrder::Order InlineJointAD::GetDofType(unsigned int i) const
{
     switch (i) {
     case 0:
     case 1:
	  return DofOrder::ALGEBRAIC;

     default:
	  ASSERT(0);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

DofOrder::Order InlineJointAD::GetEqType(unsigned int i) const
{
     switch (i) {
     case 0:
     case 1:
	  return DofOrder::ALGEBRAIC;

     default:
	  ASSERT(0);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

std::ostream& InlineJointAD::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": reaction forces [lambda1, lambda2]" << std::endl;

     if (bInitial) {
	  out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": reaction force derivatives [lambdaP1, lambdaP2]" << std::endl;
     }

     return out;
}

std::ostream& InlineJointAD::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": position constraints [c1, c2]" << std::endl;

     if (bInitial) {
	  out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": velocity constraints [cP1, cP2]" << std::endl;
     }

     return out;
}

unsigned int InlineJointAD::iGetNumPrivData(void) const
{
     return 2;
}

unsigned int InlineJointAD::iGetPrivDataIdx(const char *s) const
{
     static const struct {
	  int index;
	  char name[8];
     } data[] = {
	  { 1, "lambda1" },
	  { 2, "lambda2" }
     };

     const int N = sizeof(data) / sizeof(data[0]);

     for (int i = 0; i < N; ++i) {
	  if (0 == strcmp(data[i].name, s)) {
	       return data[i].index;
	  }
     }

     return 0;
}

doublereal InlineJointAD::dGetPrivData(unsigned int i) const
{
     switch (i) {
     case 1:
     case 2:
	  return lambda[i - 1];

     default:
	  silent_cerr("inline ad(" << GetLabel() << "): invalid private data index " << i << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

void
InlineJointAD::Output(OutputHandler& OH) const
{
     if ( bToBeOutput() )
     {
	  if ( OH.UseText(OutputHandler::LOADABLE) )
	  {
	       std::ostream& os = OH.Loadable();

	       os << std::setw(8) << GetLabel() << " " << lambda[L1] << " " << lambda[L2] << std::endl;
	  }
     }
}

void
InlineJointAD::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = *piNumCols = iWorkSpace;
}

VariableSubMatrixHandler&
InlineJointAD::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& XCurr,
		      const VectorHandler& XPrimeCurr)
{

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     grad::GradientAssVec<grad::Gradient<iWorkSpace> >::AssJac(this,
							       WorkMat.SetSparse(),
							       dCoef,
							       XCurr,
							       XPrimeCurr,
							       grad::REGULAR_JAC,
							       &dof);
#else
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
							    WorkMat.SetSparse(),
							    dCoef,
							    XCurr,
							    XPrimeCurr,
							    sp_grad::REGULAR_JAC);
#endif     
     return WorkMat;
}


SubVectorHandler&
InlineJointAD::AssRes(SubVectorHandler& WorkVec,
		      doublereal dCoef,
		      const VectorHandler& XCurr,
		      const VectorHandler& XPrimeCurr)
{
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     grad::GradientAssVec<doublereal>::AssRes(this,
					      WorkVec,
					      dCoef,
					      XCurr,
					      XPrimeCurr,
					      grad::REGULAR_RES);
#else
     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
						   WorkVec,
						   dCoef,
						   XCurr,
						   XPrimeCurr,
						   sp_grad::REGULAR_RES);
#endif
     
     return WorkVec;
}

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
template <typename T>
inline void
InlineJointAD::AssRes(GradientAssVec<T>& WorkVec,
		      doublereal dCoef,
		      const GradientVectorHandler<T>& XCurr,
		      const GradientVectorHandler<T>& XPrimeCurr,
		      enum FunctionCall func) {
     using namespace grad;
     
     typedef Vector<T, 3> Vec3;
     typedef Matrix<T, 3, 3> Mat3x3;

     const integer iFirstIndex = iGetFirstIndex();
     const integer iFirstMomentumIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomentumIndexNode2 = pNode2->iGetFirstMomentumIndex();

     Vec3 X1, X2;
     Mat3x3 R1, R2;

     pNode1->GetXCurr(X1, dCoef, func, &dof);
     pNode1->GetRCurr(R1, dCoef, func, &dof);
     pNode2->GetXCurr(X2, dCoef, func, &dof);
     pNode2->GetRCurr(R2, dCoef, func, &dof);

     Vector<T, 2> lambda;

     XCurr.GetVec(iFirstIndex + 1, lambda, 1., &dof); // Note: for algebraic variables dCoef is always one

     const Vec3 R2o2 = R2 * o2;
     const Vec3 l1 = X2 + R2o2 - X1;

     const Vec3 F1 = R1 * Vec3(e.GetCol(2) * lambda(1) + e.GetCol(3) * lambda(2));
     const Vec3 M1 = Cross(l1, F1);
     const Vec3 F2 = -F1;
     const Vec3 M2 = Cross(R2o2, F2);

     const Vec3 a = Transpose(R1) * l1 - o1;

     WorkVec.AddItem(iFirstMomentumIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstMomentumIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstMomentumIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstMomentumIndexNode2 + 4, M2);

     for (integer i = 1; i <= 2; ++i) {
	  WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a) / dCoef);
     }
}
#endif

#ifdef USE_SPARSE_AUTODIFF
template <typename T>
inline void
InlineJointAD::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
		      doublereal dCoef,
		      const sp_grad::SpGradientVectorHandler<T>& XCurr,
		      const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
		      sp_grad::SpFunctionCall func) {
     using namespace sp_grad;
     
     typedef SpColVector<T, 3> Vec3;
     typedef SpMatrix<T, 3, 3> Mat3x3;

     const integer iFirstIndex = iGetFirstIndex();
     const integer iFirstMomentumIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstMomentumIndexNode2 = pNode2->iGetFirstMomentumIndex();

     Vec3 X1, X2;
     Mat3x3 R1, R2;

     pNode1->GetXCurr(X1, dCoef, func);
     pNode1->GetRCurr(R1, dCoef, func);
     pNode2->GetXCurr(X2, dCoef, func);
     pNode2->GetRCurr(R2, dCoef, func);

     SpColVector<T, 2> lambda(2, 1);

     XCurr.GetVec(iFirstIndex + 1, lambda, 1.); // Note: for algebraic variables dCoef is always one

     const Vec3 R2o2 = R2 * o2;
     const Vec3 l1 = X2 + R2o2 - X1;

     const Vec3 F1 = R1 * (e.GetCol(2) * lambda(1) + e.GetCol(3) * lambda(2));
     const Vec3 M1 = Cross(l1, F1);
     const Vec3 F2 = -F1;
     const Vec3 M2 = Cross(R2o2, F2);

     const Vec3 a = Transpose(R1) * l1 - o1;

     WorkVec.AddItem(iFirstMomentumIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstMomentumIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstMomentumIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstMomentumIndexNode2 + 4, M2);

     for (integer i = 1; i <= 2; ++i) {
	  WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a) / dCoef);
     }
}
#endif

int
InlineJointAD::iGetNumConnectedNodes(void) const
{
     return 2;
}

void
InlineJointAD::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(iGetNumConnectedNodes());
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
}

void
InlineJointAD::SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph)
{
     const integer iFirstIndex = iGetFirstIndex();

     for (int i = 1; i <= 2; ++i) {
	  X.PutCoef(iFirstIndex + i, lambda[i - 1]);
     }
}

void InlineJointAD::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
     Update(X, XP);
}

void InlineJointAD::Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr)
{
     const integer iFirstIndex = iGetFirstIndex();

     for (int i = 1; i <= 2; ++i) {
	  lambda[i - 1] = XCurr(iFirstIndex + i);
     }
}

std::ostream&
InlineJointAD::Restart(std::ostream& out) const
{
     return out;
}

unsigned int
InlineJointAD::iGetInitialNumDof(void) const
{
     return 4u;
}

void
InlineJointAD::InitialWorkSpaceDim(
     integer* piNumRows,
     integer* piNumCols) const
{
     *piNumRows = *piNumCols = iInitialWorkSpace;
}

VariableSubMatrixHandler&
InlineJointAD::InitialAssJac(
     VariableSubMatrixHandler& WorkMat,
     const VectorHandler& XCurr)
{
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     grad::GradientAssVec<grad::Gradient<iInitialWorkSpace> >::InitialAssJac(this,
									     WorkMat.SetSparse(),
									     XCurr,
									     grad::INITIAL_ASS_JAC,
									     &dof);
#else
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::InitialAssJac(this,
								   WorkMat.SetSparse(),
								   XCurr,
								   sp_grad::INITIAL_ASS_JAC);     
#endif
     return WorkMat;
}

SubVectorHandler&
InlineJointAD::InitialAssRes(
     SubVectorHandler& WorkVec,
     const VectorHandler& XCurr)
{
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
     grad::GradientAssVec<doublereal>::InitialAssRes(this,
						     WorkVec,
						     XCurr,
						     grad::INITIAL_ASS_RES);
#else
     sp_grad::SpGradientAssVec<doublereal>::InitialAssRes(this,
							  WorkVec,
							  XCurr,
							  sp_grad::INITIAL_ASS_RES);
#endif
     return WorkVec;
}

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
template <typename T>
inline void
InlineJointAD::InitialAssRes(GradientAssVec<T>& WorkVec,
			     const GradientVectorHandler<T>& XCurr,
			     enum FunctionCall func) {
     using namespace grad;
     typedef Vector<T, 3> Vec3;
     typedef Matrix<T, 3, 3> Mat3x3;

     Vec3 X1, XP1, X2, XP2, omega1, omega2;
     Mat3x3 R1, R2;

     pNode1->GetXCurr(X1, 1., func, &dof);	// Note: during initial assembly dCoef is always one
     pNode1->GetRCurr(R1, 1., func, &dof);
     pNode1->GetVCurr(XP1, 1., func, &dof);
     pNode1->GetWCurr(omega1, 1., func, &dof);

     pNode2->GetXCurr(X2, 1., func, &dof);
     pNode2->GetRCurr(R2, 1., func, &dof);
     pNode2->GetVCurr(XP2, 1., func, &dof);
     pNode2->GetWCurr(omega2, 1., func, &dof);

     const integer iFirstIndexNode1 = pNode1->iGetFirstIndex();
     const integer iFirstIndexNode2 = pNode2->iGetFirstIndex();
     const integer iFirstIndex = iGetFirstIndex();

     T lambda[2], lambdaP[2];

     for (integer i = 1; i <= 2; ++i) {
	  XCurr.dGetCoef(iFirstIndex + i, lambda[i - 1], 1., &dof);
	  XCurr.dGetCoef(iFirstIndex + i + 2, lambdaP[i - 1], 1., &dof);
     }

     const Vec3 R2o2 = R2 * o2;
     const Vec3 l1 = X2 + R2o2 - X1;

     const Vec3 F1 = R1 * Vec3(e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]);
     const Vec3 M1 = Cross(l1, F1);
     const Vec3 FP1 = Cross(omega1, R1 * Vec3(e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]))
	  + R1 * Vec3(e.GetCol(2) * lambdaP[L1] + e.GetCol(3) * lambdaP[L2]);
     const Vec3 MP1 = -Cross(F1, XP2 + Cross(omega2, R2o2) - XP1) + Cross(l1, FP1);
     const Vec3 F2 = -F1;
     const Vec3 M2 = Cross(R2o2, F2);
     const Vec3 FP2 = -FP1;
     const Vec3 MP2 = Cross(Cross(omega2, R2o2), F2) + Cross(R2o2, FP2);

     const Vec3 a = Transpose(R1) * l1 - o1;
     const Vec3 aP = Transpose(R1) * Vec3(Cross(l1, omega1) + XP2 + Cross(omega2, R2o2) - XP1);

     WorkVec.AddItem(iFirstIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstIndexNode1 + 7, FP1);
     WorkVec.AddItem(iFirstIndexNode1 + 10, MP1);

     WorkVec.AddItem(iFirstIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstIndexNode2 + 4, M2);
     WorkVec.AddItem(iFirstIndexNode2 + 7, FP2);
     WorkVec.AddItem(iFirstIndexNode2 + 10, MP2);

     for (int i = 1; i <= 2; ++i) {
	  WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a));
	  WorkVec.AddItem(iFirstIndex + i + 2, Dot(e.GetCol(i + 1), aP));
     }
}
#endif

#if defined(USE_SPARSE_AUTODIFF)
template <typename T>
inline void
InlineJointAD::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
			     const sp_grad::SpGradientVectorHandler<T>& XCurr,
			     sp_grad::SpFunctionCall func) {
     using namespace sp_grad;
     typedef SpColVector<T, 3> Vec3;
     typedef SpMatrix<T, 3, 3> Mat3x3;

     Vec3 X1, XP1, X2, XP2, omega1, omega2;
     Mat3x3 R1, R2;

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

     T lambda[2], lambdaP[2];

     for (integer i = 1; i <= 2; ++i) {
	  XCurr.dGetCoef(iFirstIndex + i, lambda[i - 1], 1.);
	  XCurr.dGetCoef(iFirstIndex + i + 2, lambdaP[i - 1], 1.);
     }

     const Vec3 R2o2 = R2 * o2;
     const Vec3 l1 = X2 + R2o2 - X1;

     const Vec3 F1 = R1 * Vec3(e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]);
     const Vec3 M1 = Cross(l1, F1);
     const Vec3 FP1 = Cross(omega1, R1 * Vec3(e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]))
	  + R1 * Vec3(e.GetCol(2) * lambdaP[L1] + e.GetCol(3) * lambdaP[L2]);
     const Vec3 MP1 = -Cross(F1, XP2 + Cross(omega2, R2o2) - XP1) + Cross(l1, FP1);
     const Vec3 F2 = -F1;
     const Vec3 M2 = Cross(R2o2, F2);
     const Vec3 FP2 = -FP1;
     const Vec3 MP2 = Cross(Cross(omega2, R2o2), F2) + Cross(R2o2, FP2);

     const Vec3 a = Transpose(R1) * l1 - o1;
     const Vec3 aP = Transpose(R1) * Vec3(Cross(l1, omega1) + XP2 + Cross(omega2, R2o2) - XP1);

     WorkVec.AddItem(iFirstIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstIndexNode1 + 7, FP1);
     WorkVec.AddItem(iFirstIndexNode1 + 10, MP1);

     WorkVec.AddItem(iFirstIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstIndexNode2 + 4, M2);
     WorkVec.AddItem(iFirstIndexNode2 + 7, FP2);
     WorkVec.AddItem(iFirstIndexNode2 + 10, MP2);

     for (int i = 1; i <= 2; ++i) {
	  WorkVec.AddItem(iFirstIndex + i, Dot(e.GetCol(i + 1), a));
	  WorkVec.AddItem(iFirstIndex + i + 2, Dot(e.GetCol(i + 1), aP));
     }
}
#endif
#endif

bool autodiff_test_set(void)
{
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
     UserDefinedElemRead *rf = new UDERead<DeformableJointAD>;

     if (!SetUDE("deformable" "joint" "ad", rf))
     {
	  delete rf;
	  return false;
     }

     rf = new UDERead<InlineJointAD>;

     if (!SetUDE("inline" "joint" "ad", rf))
     {
	  delete rf;
	  return false;
     }

     return true;
#else
     return false;
#endif
}

//#ifndef STATIC_MODULES

extern "C"
{

     int module_init(const char *module_name, void *pdm, void *php)
     {
	  if (!autodiff_test_set())
	  {
	       silent_cerr("autodiff_test: "
			   "module_init(" << module_name << ") "
			   "failed" << std::endl);

	       return -1;
	  }

	  return 0;
     }

}

//#endif // ! STATIC_MODULE
