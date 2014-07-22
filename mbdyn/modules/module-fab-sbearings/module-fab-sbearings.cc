/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/**
    Library of bearings for "digital fabrication" machines (alpha version) [2013]
    Eduardo Okabe (okabe@unicamp.br)
    Postdoc CNPq at Aero/Polimi
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "Rot.hh"
#include "module-fab-sbearings.h"
#include "constltp.h"

// Hydrodynamic Bearing Model 01 v0 - September 12th, 2013
// Warning: IDETC/CIE 2014 version uses ADOL-C instead of finite differences !!!!

HydrodynamicBearing01::HydrodynamicBearing01(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
pNode1(0), pNode2(0)
{
   DEBUGCOUT("Entering HydrodynamicBearing01 constructor" << std::endl);

	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Hydrodynamic Bearing			\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// Hydrodynamic Bearing processing .mbd file:

   // Registry all types of hydrodynamic bearing
   const char* sBModels[] = {
          "Ocvirk",
          "Capone" "1991",
          "Sommerfeld",
          NULL
          };
   KeyTable K(HP, sBModels);

   /* lettura del tipo di vincolo */
   bearing_model = bModels(HP.IsKeyWord());

	// Read the node of bearing from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
          silent_cerr("Hydrodynamic Bearing (" << GetLabel() << ") - bearing node: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   ReferenceFrame RF1(pNode1);

   X1tilde = Vec3(Zero3);

	// Read the position offset of bearing, if supplied:
   if (HP.IsKeyWord("position")) {
        X1tilde = HP.GetPosRel(RF1);
        DEBUGCOUT("Position offset of bearing node is supplied: " << X1tilde << std::endl);
   }

   R1tilde = Mat3x3(Eye3);

	// Read the relative reference frame of bearing, if supplied:
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of bearing node is supplied" << std::endl);
        R1tilde = HP.GetRotRel(RF1);
   }

	// Read the shaft node from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
      silent_cerr("Hydrodynamic Bearing  (" << GetLabel() << ") - shaft node : structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   ReferenceFrame RF2(pNode2);

   X2tilde = Vec3(Zero3);

	// Read the position offset of bearing, if supplied:
   if (HP.IsKeyWord("position")) {
        DEBUGCOUT("Position offset of bearing node is supplied" << std::endl);
        X2tilde = HP.GetPosRel(RF2);
   }

   R2tilde = Mat3x3(Eye3);

	// Read the relative reference frame of the shaft, if supplied:
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of shaft node is supplied" << std::endl);
        R2tilde = HP.GetRotRel(RF2);
   }

   // Read the radial clearance:
   cr = HP.GetReal();

   // Read the bearing diameter:
   d0 = HP.GetReal();

   // Read the bearing length:
   l0 = HP.GetReal();

   // Read the lubricant:
   hFluid = HP.GetHydraulicFluid();

   if (HP.IsKeyWord("minimum" "speed")) {
        DEBUGCOUT("Minimum rotation speed is supplied" << std::endl);
        wr_min = HP.GetReal();
   } else { wr_min = 0.;};

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

HydrodynamicBearing01::~HydrodynamicBearing01(void)
{
    NO_OP;
}

void
HydrodynamicBearing01::Output(OutputHandler& OH) const
{

   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << XRel
      << " " << XPRel
      << " " << HForce
      << " " << HMoment
      << std::endl;
   }
}

void
HydrodynamicBearing01::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler&
HydrodynamicBearing01::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering HydrodynamicBearing01::AssJac()" << std::endl);

    FullSubMatrixHandler& WM = WorkMat.SetFull();
    /* Change the dimension of the submatrix based on the constraint demand */
    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WM.ResizeReset(iNumRows, iNumCols);

    /* Recover the index of the nodes and reaction moment variables */
    integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
    integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
    integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

    /* Set the indexes of equation */
    for (int iCnt = 1; iCnt <= 6; iCnt++) {
          WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
          WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
          WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
    }

    // Get info from nodes:
    Mat3x3 R1(pNode1->GetRCurr());
    Mat3x3 R1r(R1*R1tilde);
    Mat3x3 R1rT(R1r.Transpose());
    Mat3x3 R2(pNode2->GetRCurr());
    Vec3 X1(pNode1->GetXCurr());
    Vec3 X2(pNode2->GetXCurr());
    Vec3 XP1(pNode1->GetVCurr());
    Vec3 XP2(pNode2->GetVCurr());
    Vec3 W1(pNode1->GetWCurr());
    Vec3 W2(pNode2->GetWCurr());

    // Calculate relative position and velocities:
    Vec3 vecB = R2*X2tilde + X2 - R1*X1tilde - X1;
    Vec3 vecBP = XP2 + W2.Cross(R2*X2tilde) - XP1 - W1.Cross(R1*X1tilde);
    Vec3 dX1r = R1rT*vecB;
    Vec3 dXP1r = R1rT*(vecBP - W1.Cross(vecB));
    Vec3 dTheta12 = RotManip::VecRot(R1rT*R2*R2tilde);
    Vec3 dW1r = R1r.MulTV(W2 - W1);
    Mat3x3 R1X1tCross = Mat3x3(MatCross, R1*X1tilde);
    Mat3x3 R2X2tCross = Mat3x3(MatCross, R2*X2tilde);

    // Prepare variables vector to the calculation of the Jacobian matrix:
    doublereal Xi[12];
    doublereal Fi[6];
    doublereal JacMat[6*12];

    for (int j = 0; j < 3; j++) {
      Xi[j] = dX1r[j];
      Xi[j+3] = dXP1r[j];
      Xi[j+6] = dTheta12[j];
      Xi[j+9] = dW1r[j];
    }

    // Call force calculation routine:
    FForce(Xi, Fi);

    Vec3 HydroForce(&Fi[0]);
    Vec3 HydroMoment(&Fi[3]);

    Mat3x3 R1rFhCross = Mat3x3(MatCross, R1r*HydroForce);

    // Calculation of the Jacobian matrix:
    JacNum(Xi, JacMat, 12, 6);

    // Calculation of the partial derivatives (force/global coordinates):
    Mat3x3 dFdXP1 = (-Mat3x3(&JacMat[0], 6)*R1rT + Mat3x3(&JacMat[3*6], 6)*R1rT*Mat3x3(MatCross, W1))*dCoef
        - Mat3x3(&JacMat[3*6], 6)*R1rT;

    Mat3x3 dMdXP1 = (-Mat3x3(&JacMat[0 + 3], 6)*R1rT + Mat3x3(&JacMat[3*6 + 3], 6)*R1rT*Mat3x3(MatCross, W1))*dCoef
        - Mat3x3(&JacMat[3*6 + 3], 6)*R1rT;

    Mat3x3 dXdGP1 = R1rT*(Mat3x3(MatCross, vecB) + R1X1tCross)*dCoef;

    Mat3x3 dXPdGP1 = R1rT*((Mat3x3(MatCross, vecBP) - R1X1tCross*Mat3x3(MatCross, W1) + Mat3x3(MatCross, W1.Cross(vecB))
        - Mat3x3(MatCrossCross, vecB, W1))*dCoef + R1X1tCross + Mat3x3(MatCross, vecB));

    Mat3x3 dThdGP1 = -RotManip::DRot_I(dTheta12)*R1rT*dCoef;

    Mat3x3 dWdGP1 = R1rT*(Mat3x3(MatCross, W2 - W1) + Mat3x3(MatCross, W1))*dCoef - R1rT;

    Mat3x3 dFdGP1 = Mat3x3(&JacMat[0], 6)*dXdGP1 + Mat3x3(&JacMat[3*6], 6)*dXPdGP1 +
      Mat3x3(&JacMat[6*6], 6)*dThdGP1 + Mat3x3(&JacMat[9*6], 6)*dWdGP1;

    Mat3x3 dMdGP1 = Mat3x3(&JacMat[0 + 3], 6)*dXdGP1 + Mat3x3(&JacMat[3*6 + 3], 6)*dXPdGP1 +
      Mat3x3(&JacMat[6*6 + 3], 6)*dThdGP1 + Mat3x3(&JacMat[9*6 + 3], 6)*dWdGP1;

    Mat3x3 dFdXP2 = -dFdXP1;

    Mat3x3 dMdXP2 = -dMdXP1;

    Mat3x3 dXdGP2 = -R1rT*R2X2tCross*dCoef;

    Mat3x3 dXPdGP2 = R1rT*((R2X2tCross*Mat3x3(MatCross, W2) - Mat3x3(MatCross, W2)*R2X2tCross +
        Mat3x3(MatCross, W1)*R2X2tCross)*dCoef - R2X2tCross);

    Mat3x3 dThdGP2 = -dThdGP1;

    Mat3x3 dWdGP2 = -R1rT*Mat3x3(MatCross, W2)*dCoef + R1rT;

    Mat3x3 dFdGP2 = Mat3x3(&JacMat[0], 6)*dXdGP2 + Mat3x3(&JacMat[3*6], 6)*dXPdGP2 +
      Mat3x3(&JacMat[6*6], 6)*dThdGP2 + Mat3x3(&JacMat[9*6], 6)*dWdGP2;

    Mat3x3 dMdGP2 = Mat3x3(&JacMat[0 + 3], 6)*dXdGP2 + Mat3x3(&JacMat[3*6 + 3], 6)*dXPdGP2 +
      Mat3x3(&JacMat[6*6 + 3], 6)*dThdGP2 + Mat3x3(&JacMat[9*6 + 3], 3)*dWdGP2;

    // Perturbation of force on node 1 (dX1):
    WM.Sub(1, 1, -R1r*dFdXP1);

    // Perturbation of moment on node 1 (dX1):
    WM.Sub(3 + 1, 1, -R1X1tCross*R1r*dFdXP1 - R1r*dMdXP1);


    // Perturbation of force on node 1 (dG1):
    WM.Sub(1, 3 + 1, R1rFhCross*dCoef - R1r*dFdGP1);

    // Perturbation of moment on node 1 (dG1):
    WM.Sub(3 + 1, 3 + 1, (-R1rFhCross*R1X1tCross + R1X1tCross*R1rFhCross + Mat3x3(MatCross, R1r*HydroMoment))*dCoef
      - R1X1tCross*R1r*dFdGP1 - R1r*dMdGP1);


    // Perturbation of force on node 1 (dX2):
    WM.Sub(1, 6 + 1, -R1r*dFdXP2);

    // Perturbation of moment on node 1 (dX2):
    WM.Sub(3 + 1, 6 + 1, -R1X1tCross*R1r*dFdXP2 - R1r*dMdXP2);

    // Perturbation of force on node 1 (dG2):
    WM.Sub(1, 9 + 1, -R1r*dFdGP2);

    // Perturbation of moment on node 1 (dG2):
    WM.Sub(3 + 1, 9 + 1, -R1X1tCross*R1r*dFdGP2 - R1r*dMdGP2);


    /* Perturbation on node 2: */

    // Perturbation of force on node 2 (dX1):
    WM.Sub(6 + 1, 1, R1r*dFdXP1);

    // Perturbation of moment on node 2 (dX1):
    WM.Sub(9 + 1, 1, R2X2tCross*R1r*dFdXP1 + R1r*dMdXP1);

    // Perturbation of force on node 2 (dG1):
    WM.Sub(6 + 1, 3 + 1, -R1rFhCross*dCoef + R1r*dFdGP1);

    // Perturbation of moment on node 2 (dG1):
    WM.Sub(9 + 1, 3 + 1, (-R2X2tCross*R1rFhCross - Mat3x3(MatCross, R1r*HydroMoment))*dCoef
      + R2X2tCross*R1r*dFdGP1 + R1r*dMdGP1);

    // Perturbation of force on node 2 (dX2):
    WM.Sub(6 + 1, 6 + 1, R1r*dFdXP2);

    // Perturbation of moment on node 2 (dX2):
    WM.Sub(9 + 1, 6 + 1, R2X2tCross*R1r*dFdXP2 + R1r*dMdXP2);

    // Perturbation of force on node 2 (dG2):
    WM.Sub(6 + 1, 9 + 1, R1r*dFdGP2);

    // Perturbation of moment on node 2 (dG2):
    WM.Sub(9 + 1, 9 + 1, R1rFhCross*R2X2tCross*dCoef + R2X2tCross*R1r*dFdGP2 + R1r*dMdGP2);

    return WorkMat;
}

void
HydrodynamicBearing01::JacNum(doublereal Xi[], doublereal JacMat[], integer n, integer m)
{
   /*
      Forward finite difference routine

      It evaluates n + 1 times the function force to calculate
      the Jacobian matrix numerically.

   */

   doublereal Xj[n];
   doublereal Fj[m];
   doublereal F0[m];
   doublereal deltaX0 = 1e-8;
   doublereal deltaX1 = deltaX0;

   // Copy variables to array Xj:
   memcpy(Xj, Xi, sizeof(doublereal)*n);

   // Calculate force in the central point:
   FForce(Xi, F0);
   DEBUGCOUT("JacMat(): " << std::endl);

   for (int j = 0; j < n; j++) {
      // Apply delta to jth variable:
      if (abs(Xj[j])>deltaX0) deltaX1 = deltaX0*abs(Xj[j]); else deltaX1 = deltaX0;
      if (deltaX1==0.) deltaX1 = deltaX0;

      Xj[j] += deltaX1;

      // Calculate force:
      FForce(Xj, Fj);
      //std::cout << j << ": " << Fj[0] << ", " << Fj[1] << ", " << Fj[2] << "\n";

      // Calculate forward finite difference:
      for (int k = 0; k < m; k++) JacMat[k+j*m] = (Fj[k]-F0[k])/deltaX1;

      // Return jth variable to its original value:
      Xj[j] = Xi[j];
   }

}


SubVectorHandler&
HydrodynamicBearing01::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering HydrodynamicBearing01::AssRes()" << std::endl);

    /* Change the dimension of the vector based on the constraint demand */
    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.ResizeReset(iNumRows);

    /* Get the index of the nodes and reaction moment variables */
    integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

    for (int iCnt = 1; iCnt <= 6; iCnt++) {
          WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
    }

    // Get info from nodes:
    Mat3x3 R1(pNode1->GetRCurr());
    Mat3x3 R1r(R1*R1tilde);
    Mat3x3 R2(pNode2->GetRCurr());
    Vec3 X1(pNode1->GetXCurr());
    Vec3 X2(pNode2->GetXCurr());
    Vec3 XP1(pNode1->GetVCurr());
    Vec3 XP2(pNode2->GetVCurr());
    Vec3 W1(pNode1->GetWCurr());
    Vec3 W2(pNode2->GetWCurr());

    // Calculate relative position and velocity:
    Vec3 vecB = R2*X2tilde + X2 - R1*X1tilde - X1;
    Vec3 dX1r = R1r.MulTV(vecB);
    Vec3 dXP1r = R1r.MulTV(XP2 + W2.Cross(R2*X2tilde) - XP1 - W1.Cross(R1*X1tilde + vecB));
    Vec3 dW1r = R1r.MulTV(W2 - W1);
    Vec3 dTheta12 = RotManip::VecRot(R1r.MulTM(R2)*R2tilde);

    // Calculate hydrodynamic force:
    doublereal Xi[12];
    doublereal Fi[6];
    for (int j = 0; j < 3; j++) {
      Xi[j] = dX1r[j];
      Xi[j+3] = dXP1r[j];
      Xi[j+6] = dTheta12[j];
      Xi[j+9] = dW1r[j];
    }

    // Call force calculation routine:
    FForce(Xi, Fi);

    Vec3 HydroForce(&Fi[0]);
    Vec3 HydroMoment(&Fi[3]);

    HForce = HydroForce;
    HMoment = HydroMoment;

    WorkVec.Sub(1, R1r*HydroForce);
    WorkVec.Sub(3 + 1,  (R1*X1tilde).Cross(R1r*HydroForce) + R1r*HydroMoment);
    WorkVec.Add(6 + 1, R1r*HydroForce);
    WorkVec.Add(9 + 1,  (R2*X2tilde).Cross(R1r*HydroForce) + R1r*HydroMoment);

    DEBUGCOUT("HydrodynamicBearing01::AssRes(), Force-Moment: " << Fi[1] << ", " <<
    Fi[2] <<  ", "  << Fi[3] <<   ", " << HydroForce <<  ", " << HydroMoment << std::endl);

    return WorkVec;
}


void
HydrodynamicBearing01::FForce(doublereal Xi[], doublereal Fi[]) const
{
   /*
      Vector Xi:
      Xi[0-2]: relative position;
      Xi[3-5]: relative linear velocity;
      Xi[6-8]: relative angle;
      Xi[9-11]: relative angular velocity;

      Vector Fi:
      Fi[0-2]: forces;
      Fi[3-5]: moments (couples);
   */

   for (int j=0; j<6; j++) Fi[j] = 0.;

   switch (bearing_model) {
      case OCVIRK:
         {
         /**
            Ocvirk model [1953] - Short Bearing - Cylindrical
            DuBois, G., Ocvirk, F., Analytical derivation and experimental evaluation of
            short-bearing approximation for full journal bearings.
            Report 1157 - NACA
         */

         doublereal r = d0/2. + cr; // Verify radius calculation!
         doublereal v = atan2(Xi[2], Xi[1]);
         doublereal e = sqrt(Xi[1]*Xi[1] + Xi[2]*Xi[2])/cr;

         if (e==0.) break;

         doublereal Py = 4.*e*e/pow(1.-e*e, 2.);
         doublereal Pz = M_PI*e/pow(1.-e*e, 3./2.);
         doublereal So = hFluid->dGetViscosity()*Xi[9]*l0*l0/(4.*cr*cr);

         Py = So*r*l0*Py;
         Pz = So*r*l0*Pz;

         Fi[1] = Py*cos(v) + Pz*sin(v);
         Fi[2] = Py*sin(v) - Pz*cos(v);
         break;
      }

      case CAPONE91:
      {
         /**
            Capone model [1991] - Short Bearing - Cylindrical
            Capone, G., Descrizione analitica del campo di forze fluidodinamica
            nei cuscinetti cilindrici lubrificati
            L'Energia Elettrica No. 3 - 1991
         */


         doublereal r0 = d0/2. + cr; // Verify radius calculation!
         doublereal wr = Xi[9];
         if (wr<wr_min) wr = wr_min;

         doublereal x = Xi[1]/cr;
         doublereal y = Xi[2]/cr;
         doublereal dx = Xi[1+3]/(cr*wr);
         doublereal dy = Xi[2+3]/(cr*wr);

         doublereal v = atan2(-y-2.*dx, -x+2.*dy);
         if (v != v) v = 0.;

         doublereal b = y*cos(v) - x*sin(v);
         doublereal b1 = x*cos(v) + y*sin(v);
         doublereal c = 1. - x*x - y*y;
         doublereal h = sqrt(c);

         if (c>0. && c<1.) {
             doublereal g = 2./h*(M_PI/2. + atan2(b,h));
             doublereal v1 = (2.+b*g)/c;
             doublereal f1 = b1/(1.-b1*b1);
             doublereal c1 = sqrt(pow(x-2.*dy, 2.) + pow(y+2.*dx, 2.))/c;
             doublereal fy = -c1*(3.*x*v1 - sin(v)*g - 2.*cos(v)*f1);
             doublereal fz = -c1*(3.*y*v1 + cos(v)*g - 2.*sin(v)*f1);
             doublereal So = hFluid->dGetViscosity()*wr*l0*l0/(4.*cr*cr);
             Fi[1] = So*r0*l0*fy;
             Fi[2] = So*r0*l0*fz;
        } else {
             Fi[1] = 0.;
             Fi[2] = 0.;
        };

         break;
      }

      case SOMMERFELD:
      {
         /**
            Half Sommerfeld model [1904] - Long Bearing - Cylindrical
            Sommerfeld, A., Zur Hydrodynamischen theorie der Schmiermittelreibung
            Z. Angew. Math. Phys., n. 50, pp. 97-155, 1904.
         */

         doublereal r0 = d0/2. + cr; // Verify radius calculation!
         doublereal v = atan2(Xi[2],Xi[1]);
         doublereal e = sqrt(Xi[2]*Xi[2] + Xi[1]*Xi[1])/cr;

         doublereal Py = M_PI*e/((2.+e*e)*sqrt(1.-e*e));
         doublereal Px = 2.*e*e/((2.+e*e)*(1.-e*e));

         doublereal So = 6.*hFluid->dGetViscosity()*Xi[9]*r0*pow(r0/cr,2.)*l0;
         Px = So*Px;
         Py = So*Py;

         Fi[1] = -Px*cos(v) - Py*sin(v);
         Fi[2] = -Px*sin(v) + Py*cos(v);

         break;
      }
   }
}


doublereal HydrodynamicBearing01::hi(doublereal cr, doublereal x, doublereal y, doublereal v) const
{
   return cr-x*cos(v)-y*sin(v);
};


doublereal HydrodynamicBearing01::dhi(doublereal dx, doublereal dy, doublereal v) const
{
   return -dx*cos(v)-dy*sin(v);
};


unsigned int
HydrodynamicBearing01::iGetNumPrivData(void) const
{
	return 0;
}

int
HydrodynamicBearing01::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
HydrodynamicBearing01::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
}

void
HydrodynamicBearing01::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
HydrodynamicBearing01::Restart(std::ostream& out) const
{
	return out << "# HydrodynamicBearing01: not implemented" << std::endl;
}

unsigned int
HydrodynamicBearing01::iGetInitialNumDof(void) const
{
	return 0;
}

void
HydrodynamicBearing01::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
HydrodynamicBearing01::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
HydrodynamicBearing01::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


unsigned int
HydrodynamicBearing01::iGetNumDof(void) const
{
   return 0;
}


DofOrder::Order
HydrodynamicBearing01::GetDofType(unsigned int i) const
{
   return DofOrder::DIFFERENTIAL;
}


void
HydrodynamicBearing01::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{

   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R1r(R1*R1tilde);
   Mat3x3 R1rT(R1r.Transpose());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 X1(pNode1->GetXCurr());
   Vec3 X2(pNode2->GetXCurr());
   Vec3 XP1(pNode1->GetVCurr());
   Vec3 XP2(pNode2->GetVCurr());
   Vec3 W1(pNode1->GetWCurr());
   Vec3 W2(pNode2->GetWCurr());

   // Calculate relative position and velocities:
   Vec3 vecB = R2*X2tilde + X2 - R1*X1tilde - X1;
   Vec3 vecBP = XP2 + W2.Cross(R2*X2tilde) - XP1 - W1.Cross(R1*X1tilde);
   XRel = R1rT*vecB;
   XPRel = R1rT*vecBP;

}


/**
    Rolling Bearing Constitutive Law (Gargiulo Model - Hambric, 2013) v0 - January 17th, 2014

    Gargiulo model presented in this article:
    Hambric, S. A., Shepherd, M. R., Campbell, R. L., and
    Hanford, A. D., 2013. “Simulations and measurements of the vibroacoustic effects of
    replacing rolling element bearings with journal bearings in a simple gearbox”.
    Journal of Vibration and Acoustics-Transactions of the Asme, 135(3).
*/

template <class T, class Tder>
class RollingBearingConstitutiveLaw
: public ConstitutiveLaw<T, Tder> {
private:
	doublereal DSphere;         // Sphere (or ball) diameter [m]
	doublereal ZRollElem;       // Number of rolling elements (balls or rollers)
	doublereal LRoller;         // Length of the rollers [m]
	doublereal AlphaContact;    // Contact angle (zero if contact is purely radial)
	doublereal cDamping;         // Proportional damping coefficient
    bool bIsBallBearing;        // Is it a ball bearing? If not, it is a roller bearing.

public:
	RollingBearingConstitutiveLaw(bool bBall, doublereal dDL, doublereal dNElms, doublereal dAlpha, doublereal dDamp)
	: DSphere(dDL), ZRollElem(dNElms), LRoller(dDL), AlphaContact(dAlpha), bIsBallBearing(bBall), cDamping(dDamp) {
		NO_OP;
	};

	virtual ~RollingBearingConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef RollingBearingConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(bIsBallBearing, DSphere, ZRollElem, AlphaContact, cDamping));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "rolling bearing, " << DSphere;
	};

	virtual void Update(const Vec3& Eps, const Vec3& EpsPrime) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
        ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		doublereal d_axial = fabs(Eps(1));
		doublereal d_radial = sqrt(Eps(2)*Eps(2) + Eps(3)*Eps(3));

		if (bIsBallBearing) {

            doublereal c_radial = 1.17732411925961e7;
            doublereal c_axial = 3.43778642823805e7;

            doublereal f_radial = pow(d_radial*c_radial,1.5)*sqrt(DSphere*ZRollElem*ZRollElem*
                pow(cos(AlphaContact),5.));

            doublereal f_axial = pow(d_axial*c_axial,1.5)*sqrt(DSphere*ZRollElem*ZRollElem*
                pow(sin(AlphaContact),5.));

            doublereal k_radial(0.);
            if (d_radial > 0.) k_radial = f_radial/d_radial;

            doublereal k_axial(0.);
            if (d_axial > 0.) k_axial = f_axial/d_axial;

            Mat3x3 FDETemp(Eye3);
            FDETemp(1,1) = k_axial;
            FDETemp(2,2) = k_radial;
            FDETemp(3,3) = k_radial;

            ConstitutiveLaw<T, Tder>::FDE = FDETemp;
            ConstitutiveLaw<T, Tder>::FDEPrime = FDETemp*cDamping;
            ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::FDE*ConstitutiveLaw<T, Tder>::Epsilon +
                ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;

            Vec3 dF = ConstitutiveLaw<T, Tder>::F;

            DEBUGCOUT("Rolling Bearing, f_radial, d_radial, f_axial, d_axial: " << f_radial << ", " <<
                d_radial << ", " << f_axial << ", " << d_axial << std::endl);
            DEBUGCOUT("Rolling Bearing, Eps: " << Eps << std::endl);
            DEBUGCOUT("Rolling Bearing, F: " << dF << std::endl);

        } else {

            doublereal c1 = 10./9.;

            doublereal c_radial = 8.54650149899497e8;
            doublereal c_axial = 3.64650730623785e9;

            doublereal f_radial = pow(d_radial*c_radial,c1)*ZRollElem*pow(LRoller, 0.8*c1)*
                pow(cos(AlphaContact),1.9*c1);

            doublereal f_axial = pow(d_axial*c_axial,c1)*ZRollElem*pow(LRoller, 0.8*c1)*
                pow(sin(AlphaContact),1.9*c1);

            doublereal k_radial(0.);
            if (d_radial > 0.) k_radial = f_radial/d_radial;

            doublereal k_axial(0.);
            if (d_axial > 0.) k_axial = f_axial/d_axial;


            Mat3x3 FDETemp(Eye3);
            FDETemp(1,1) = k_axial;
            FDETemp(2,2) = k_radial;
            FDETemp(3,3) = k_radial;

            ConstitutiveLaw<T, Tder>::FDE = FDETemp;
            ConstitutiveLaw<T, Tder>::FDEPrime = FDETemp*cDamping;
            ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::FDE*ConstitutiveLaw<T, Tder>::Epsilon +
                ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;

            Vec3 dF = ConstitutiveLaw<T, Tder>::F;

            DEBUGCOUT("Rolling Bearing, f_radial, d_radial, f_axial, d_axial: " << f_radial << ", " <<
                d_radial << ", " << f_axial << ", " << d_axial << std::endl);
            DEBUGCOUT("Rolling Bearing, Eps: " << Eps << std::endl);
            DEBUGCOUT("Rolling Bearing, EpsPrime: " << EpsPrime << std::endl);
            DEBUGCOUT("Rolling Bearing, FDEPrime: " << FDETemp*cDamping << std::endl);
            DEBUGCOUT("Rolling Bearing, F: " << dF << std::endl);

        }
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct RollingBearingCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

        bool bBall(false);

        if (HP.IsKeyWord("ball" "bearing")) {
            bBall = true;
        } else if (HP.IsKeyWord("roller" "bearing")) {
            bBall = false;
        }

        // Ball diameter or roller length
		doublereal dDL = HP.GetReal();
		/*if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness "
				"at line " << HP.GetLineData() << std::endl);
		}*/

        // Number of rolling elements
		doublereal dNElms = HP.GetReal();

        // Contact angle [radians]
		doublereal dAlpha = HP.GetReal();

        doublereal dDamp = 0.;

        if (HP.IsKeyWord("proportional")) {
            dDamp = HP.GetReal();
        }

		typedef RollingBearingConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(bBall, dDL, dNElms, dAlpha, dDamp));

		return pCL;
	};
};


extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<HydrodynamicBearing01>;

	if (!SetUDE("hydrodynamic" "bearing", rf1)) {
		delete rf1;

		silent_cerr("module-fabricate-bearings: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new RollingBearingCLR<Vec3, Mat3x3>;

	if (!SetCL3D("rolling" "bearing", rf3D)) {
		delete rf3D;

		silent_cerr("RollingBearingConstitutiveLaw3D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}
