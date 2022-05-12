/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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
 * Author:	Louis Gagnon <louis.gagnon.10@ulaval.ca>
 *		Departement de genie mecanique
 *		Universite Laval
 *		http://www.gmc.ulaval.ca
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <limits>
#include <cfloat>
#include <limits>


#include "module-wheel4.h"
#include "dataman.h"
#include "userelem.h"
#include "simentity.h"
#include "body.h"

Wheel4::Wheel4(unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
firstRes(true)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	wheel4							\n"
"Author: 	Louis Gagnon <louis.gagnon.10@ulaval.ca>		\n"
"Organization:	Departement de genie mecanique				\n"
"		Universite Laval					\n"
"		http://www.gmc.ulaval.ca				\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"		Marco Morandini <morandini@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"L. Gagnon, M. J. Richard, P. Masarati, M. Morandini, and G. Dore. Multibody simulation of tires operating on an uneven road. In Multibody Dynamics 2011, 4-7 July 2011"
				" And soon coming : L. Gagnon, M. J. Richard, P. Masarati, M. Morandini, and G. Dore. A Free Implicit Rigid Ring Tire Model"
"	All rights reserved						\n"
"	Wheel4 requires the ginac libraries to be installed		\n"
"									\n"
"Nodes:									\n"
"     -	Wheel								\n"
"     -	Ring								\n"
" 									\n"
"     -	Patch								\n"
"									\n"
"Note: 									\n"
"     -	The Ring and the Wheel structural nodes must be connected 	\n"
"	by a 6 DoF viscoelastic element	\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */	
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	i = Vec3(1.,0.,0.); // unit vector in x-dir
	j = Vec3(0.,1.,0.); // unit vector in y-dir
	k = Vec3(0.,0.,1.); // unit vector in z-dir

	// read wheel node
	pWheel = pDM->ReadNode<StructNode, StructDispNode, Node::STRUCTURAL>(HP); // input 1, wheel node (or ring node if using swift)
	pWheelB = pDM->ReadElem<Body, Elem::BODY>(HP); // input, wheel body

	// read wheel axle
	ReferenceFrame RF = ReferenceFrame(pWheel); //makes reference frame from wheel node
	WheelAxle = HP.GetVecRel(RF); //converts value to RF reference frame // input
	doublereal dWheelAxle = WheelAxle.Dot();
	if (dWheelAxle < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Wheel4(" << GetLabel() << "): "
			"wheel axle is too small "
			"for numeric limits" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	WheelAxle /= sqrt(dWheelAxle);

	// wheel geometry
	dR_0 = HP.GetReal(); // Rigid ring radius, (in MF calculations, a home-made effective radius is used)
	zZero = Vec3(1.,1.,0.); // this permits the multiplication of the contact patch displacement vectors by zero in the z direction, to allow for ground controlled position

	// SWIFT
	bSwift = false;
	if (HP.IsKeyWord("swift")) {
		bSwift = true;
	pRing = pDM->ReadNode<StructNode, StructDispNode, Node::STRUCTURAL>(HP); // get ring node if swift
	pRingB = pDM->ReadElem<Body, Elem::BODY>(HP); // get ring body

	Xpa = pRing->GetXCurr() - k*0.98*dR_0;
	XpaPrev=Xpa;

	XparPrev=XparPrevPrev=XparpPrev=XparpPrevPrev=FintPrev=FintPrevPrev=Vec3(0.,0.,0.); // simply for proper initialization


	dRoad = dRoadInitial = Xpa(3); //getting z component of initial elevation
//	dRoadPrev = XpaPrev(3); // not necessary
	Kpa = HP.GetVec3();
	pKpa = HP.GetDriveCaller(); //time variant driver for Kpa
	Cpa = HP.GetVec3();
	pCpa = HP.GetDriveCaller(); //time variant driver for Cpa
	Vpa = HP.GetVec3(); //velocity of patch in wheel ref frame
	VpaPrev = Vpa;
	VpaBC=Vpa;//
	XpaBC=Xpa; //
	Mpa = HP.GetReal(); // mass of patch
	dt = 0;
	curTime = 0;
	tdc.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));
	pRoad = HP.GetDriveCaller(); //road profile driver, makes a function "f(x)"
	dPls = HP.GetReal(); // patch contact-length to elleptical cam tandem base parameter (Zegelaar eq. 4.15) dPls = l_s/(2a)
	dR_a1 = HP.GetReal(); // r_a1 contact length parameter from Besselink eq. 4.85
	dR_a2 = HP.GetReal(); // r_a2 contact length parameter from Besselink eq. 4.85
	rRatio = HP.GetReal(); // ratio btwn portion of ring in contact with patch and the total ring
	Krz = HP.GetReal(); // vertical wheel-ring stiffness
	nPrev = k;
	dXxPrev = 0;
	deltaPrev = 0;
	bFirstAP=true;
	bFirstAC=false;
	}

	bLoadedRadius = false;
	if (HP.IsKeyWord("loadedRadius"))
	{
		bLoadedRadius = true;
	}
	// friction
	bSlip = false;
	if (HP.IsKeyWord("slip")) {
		bSlip = true;

		/*
		 * Parametri di attrito
		 */
		pMuX0 = HP.GetDriveCaller(); //makes a function "f(x)"
		pMuY0 = HP.GetDriveCaller();
		pTr = HP.GetDriveCaller();
		pMzr = HP.GetDriveCaller();
		S_ht = HP.GetReal();
		S_hf = HP.GetReal();
		q_sy1 = HP.GetReal(); // should be between 0.01 and 0.02 usually...
		q_sy3 = HP.GetReal(); //
		dvao  = HP.GetReal(); // velocity for rolling resistance velocity influence factor (reference velocity)
	
		TRH = 0.; //prevents division by zero at null x-velocity at the price of losing validity for velocities above -TRH*1.1 and below TRH*1.1
		TRHA = 0.; //buffer used to prevent division by zero
		TRHT = 0.; //prevents divison by zero for computation of the angle of the car (and wheels)
		TRHTA = 0.; //buffer used on angle zero division prevention
		TRHC = 0.; //cap on kappa

		if (HP.IsKeyWord("threshold")) {
			TRH = HP.GetReal();
			TRHA = HP.GetReal();
			TRHT = HP.GetReal();
			TRHTA = HP.GetReal();
			TRHC = HP.GetReal();
			TdLs = HP.GetReal();
			TdReDiv = HP.GetReal();
			TdR_e = HP.GetReal();
			dt_On = HP.GetReal();
			dt_maxH = HP.GetReal();
//			dt_numAhead = HP.GetInt();
			dt_Res = HP.GetReal();
			dt_maxstep = HP.GetReal();
			dt_minstep = HP.GetReal();
			dt_maxF = HP.GetReal(); // maximum divison of dt per timestep
			dt_minF = HP.GetReal(); // minimum division of dt per timestep
			dt_minStepsCycle = HP.GetInt(); // minimum number of steps wanted in a force cycle (approx., influences dt)
			dt_divF = HP.GetReal(); // factor by which to divide the timestep if the force oscillates more than wanted (as determined by TminS)
			dt_divF3 = HP.GetReal(); // 3 sign chg
			dt_divF4 = HP.GetReal(); // 4 sign chg or more
			RDA = HP.GetReal();
			RDB = HP.GetReal();
			RDL = HP.GetReal();
}
//		if (HP.IsKeyWord("master")) {
//			int numWheels = HP.GetInt();
//			for(int iCnt = 1; iCnt <= numWheels; iCnt++) {
//			//	pWheelE = (Elem *)pDM->ReadElem(HP, Elem::LOADABLE); // input, wheel elem
//
//				std::cout << "Wheel: " << HP.GetInt() << std::endl;
//			}
//		}

	}

	std::vector<int> row;
	for (int jCnt = 0; jCnt < 3; jCnt++) {
		row.push_back(0);
	}
	for(int iCnt = 0; iCnt < dt_minStepsCycle; iCnt++) {
		FintSignCk.push_back(row);
	}


	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
	
	std::ostream& out = pDM->GetLogFile();
	out << "wheel4: " << uLabel
		<< " " << pWheel->GetLabel()	//node label
		<< " " << WheelAxle		//wheel axle
//		<< " " << dRadius		//wheel radius
//		<< " " << dInternalRadius	//wheel internal radius
		<< std::endl;
}

Wheel4::~Wheel4(void)
{
	NO_OP;
}

void
Wheel4::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::LOADABLE)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

		        std::ostringstream os;
		        os << "loadable." << GetLabel();
		        (void)OH.CreateVar(os.str(), "wheel4");

		        // joint sub-data
		        os << '.';
		        std::string name = os.str();

			/* Add NetCDF (output) variables to the BinFile object
			 * and save the NcVar* pointer returned from add_var
			 * as handle for later write accesses.
			 * Define also variable attributes */
			OH.CreateVar<Vec3>(name + "Fint",
                                                OutputHandler::Dimensions::Force,
                                                "force btwn ring and patch acting on patch (x, y, z)");
			OH.CreateVar<Vec3>(name + "Xpar",
                                                OutputHandler::Dimensions::Length,
                                                "rel pos of patch (x, y, z)");
			OH.CreateVar<Vec3>(name + "Xparp",
                                                OutputHandler::Dimensions::Length,
                                                "rel pos of patch expressed in the ring non-rotating reference frame (x, y, z)");
			OH.CreateVar<doublereal>(name + "dXxProj",
                                                OutputHandler::Dimensions::Length,
                                                "patch center point x-value of the road profile (attention: this is not the same thing as the patch position)");
			OH.CreateVar<doublereal>(name + "dRoad",
                                                OutputHandler::Dimensions::Length,
                                                "road height");
			OH.CreateVar<Vec3>(name + "F",
                                                OutputHandler::Dimensions::Force,
                                                "Road-patch friction force expressed in absolute reference frame (x, y, z)");
			OH.CreateVar<doublereal>(name + "Fn",
                                                OutputHandler::Dimensions::Force,
                                                "Normal force at the patch");
			OH.CreateVar<doublereal>(name + "debug",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Debugging variable");
			OH.CreateVar<doublereal>(name + "dSr",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dSr");
			OH.CreateVar<doublereal>(name + "ddistM",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "ddistM");
			OH.CreateVar<doublereal>(name + "Fcent",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Fcent");
			OH.CreateVar<doublereal>(name + "dLs",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dLs");
			OH.CreateVar<doublereal>(name + "R_e",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "R_e");
			OH.CreateVar<doublereal>(name + "dSa",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dSa");
			OH.CreateVar<doublereal>(name + "dvax",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dvax");
			OH.CreateVar<doublereal>(name + "dvx",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dvx");
			OH.CreateVar<doublereal>(name + "dvay",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dvay");
			OH.CreateVar<doublereal>(name + "dMuY",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dMuY");
			OH.CreateVar<doublereal>(name + "dMuX",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dMuX");
			OH.CreateVar<doublereal>(name + "KE",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "KE");
			OH.CreateVar<doublereal>(name + "PE",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "PE");
			OH.CreateVar<doublereal>(name + "E",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "E");
			OH.CreateVar<doublereal>(name + "dRoadAhead",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dRoadAhead");
			OH.CreateVar<doublereal>(name + "dRoadBehind",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dRoadBehind");
			OH.CreateVar<doublereal>(name + "dCt",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "dCt");
			OH.CreateVar<Vec3>(name + "M",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "M");
			OH.CreateVar<Vec3>(name + "distM",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "distM");
			OH.CreateVar<Vec3>(name + "n",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "n");
			OH.CreateVar<Vec3>(name + "Xpa",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Xpa");
			OH.CreateVar<Vec3>(name + "Vpa",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Vpa");
			OH.CreateVar<Vec3>(name + "Vpar",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Vpar");
			OH.CreateVar<Vec3>(name + "fwd",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "fwd");
			OH.CreateVar<Vec3>(name + "fwdRing",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "fwdRing");
			OH.CreateVar<Vec3>(name + "fwdRingFlat",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "fwdRingFlat");
			OH.CreateVar<Vec3>(name + "pcRing",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "pcRing");
			OH.CreateVar<Vec3>(name + "VparWheel",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "VparWheel");
			OH.CreateVar<Vec3>(name + "Fr",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Fr");
			OH.CreateVar<Vec3>(name + "Mz",
                                                OutputHandler::Dimensions::UnknownDimension,
                                                "Mz");


		}
#endif // USE_NETCDF
	}

}

//void Wheel4::NetCDFPrepare(OutputHandler &OH, char &buf) const
//{
//	strcpy(&buf[l], "dRoad");
//	Var_dRoad = pBinFile->add_var(buf, ncDouble,
//		OH.DimTime(), OH.DimV1());
//	if (Var_dRoad == 0) {
//		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//	}
//
//	if (!Var_dRoad->add_att("units", "m")) {
//		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//	}
//
//	if (!Var_dRoad->add_att("type", "doublereal")) {
//		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//	}
//
//	if (!Var_dRoad->add_att("description",
//		"road height"))
//	{
//		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//	}
//}

void
Wheel4::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::LOADABLE)) {
			OH.WriteNcVar(Var_Fint, Fint);
			OH.WriteNcVar(Var_Xpar, Xpar);
			OH.WriteNcVar(Var_Xparp, Xparp);
			OH.WriteNcVar(Var_dXxProj, dXxProj);
			OH.WriteNcVar(Var_dRoad, dRoad);
			OH.WriteNcVar(Var_F, F);
			OH.WriteNcVar(Var_Fn, Fn);
			OH.WriteNcVar(Var_debug, dDebug);
			OH.WriteNcVar(Var_dSr, dSr);
			OH.WriteNcVar(Var_ddistM, ddistM);
			OH.WriteNcVar(Var_Fcent, Fcent);
			OH.WriteNcVar(Var_dLs, dLs);
			OH.WriteNcVar(Var_R_e, R_e);
			OH.WriteNcVar(Var_dSa, dSa);
			OH.WriteNcVar(Var_dvax, dvax);
			OH.WriteNcVar(Var_dvx, dvx);
			OH.WriteNcVar(Var_dvay, dvay);
			OH.WriteNcVar(Var_dMuY, dMuY);
			OH.WriteNcVar(Var_dMuX, dMuX);
			OH.WriteNcVar(Var_KE, KE);
			OH.WriteNcVar(Var_PE, PE);
			OH.WriteNcVar(Var_E, E);
			OH.WriteNcVar(Var_dRoadAhead, dRoadAhead);
			OH.WriteNcVar(Var_dRoadBehind, dRoadBehind);
			OH.WriteNcVar(Var_dCt, dCt);
			OH.WriteNcVar(Var_M, M);
			OH.WriteNcVar(Var_distM, distM);
			OH.WriteNcVar(Var_n, n);
			OH.WriteNcVar(Var_Xpa, Xpa);
			OH.WriteNcVar(Var_Vpa, Vpa);
			OH.WriteNcVar(Var_Vpar, Vpar);
			OH.WriteNcVar(Var_fwd, fwd);
			OH.WriteNcVar(Var_fwdRing, fwdRing);
			OH.WriteNcVar(Var_fwdRingFlat, fwdRingFlat);
			OH.WriteNcVar(Var_pcRing, pcRing);
			OH.WriteNcVar(Var_VparWheel, VparWheel);
			OH.WriteNcVar(Var_Fr, Fr);
			OH.WriteNcVar(Var_Mz, Mz);
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::LOADABLE)) {

			std::ostream& out = OH.Loadable();



			out << std::setw(8) << GetLabel()	// 1:	label
				<< " " << dvax		// 2: velocity of wheel in x-dir
				<< " " << dvay		// 3: velocity of wheel in y-dir
				<< " " << dvx		// 4: relative speed between center of wheel and contact point on tire in the forward direction
				<< " " << M			// 5-7:	moment
				<< " " << distM		// 8-10: moment arm on ring
				<< " " << dSr			// 11:	slip ratio
				<< " " << 180./M_PI*dAlpha	// 12:	slip angle
				<< " " << dMuX			// 13:	longitudinal friction coefficient
				<< " " << dMuY			// 14:	lateral friction coefficient
				<< " " << dRoad			// 15:	road height
				<< " " << n			// 16-18:	road normal
				<< " " << Xpa		// 19-21: pos of patch
				<< " " << Vpa		// 22-24: vel of patch
				<< " " << Xpar		// 25-27: rel pos of patch
				<< " " << Vpar		// 28-30: rel vel of patch
				<< " " << -Fint		// 31-33: force btwn ring and patch acting on patch
				<< " " << Fint_ring	// 34-36: force btwn ring and patch acting on ring
				<< " " << fwd		// 37-39: forward direction vector of the wheel
				<< " " << fwdRing	// 40-42: forward direction vector of the ring
				<< " " << fwdRingFlat	// 43-45: forward direction vector of the ring without the slope of the profile
				<< " " << pcRing	// 46-48: point of contact on ring between ring and springs
				<< " " << Fn		// 49: normal force for Pacejka's formulae
				<< " " << VparWheel // 50-52: relative vel btwn patch and wheel
				<< " " << KE   // 53: wheel+ring+patch kin energy
				<< " " << PE   // 54: wheel+ring pot energy (patch is not subject to gravity)
				<< " " << E   // 55: wheel+ring+patch tot energy
				<< " " << R_e // 56: virtually calculated effective rolling radius
				<< " " << dLs	// 57: half length of the tandem elliptical cam follower
				<< " " << ddistM // 58: home made loaded radius (distance between ring center and patch center)
				<< " " << Xparp // 59-61: distance between ring contact point and patch as seen from the ring in its own reference frame (not rotated with road slope)
				<< " " << Fcent // 62: Fcent, centrifugal force added to tire
				<< " " << Fr // 63-65: Rolling resistance force vector (only applied as a moment)
				<< " " << Mz // 66-68: Aligning moment
				<< " " << dXxProj // 69: center point x-value of the road profile (this is not actual position, but only position on the input road file
				<< " " << dRoadAhead // 70: front edge x-point of the tandem
				<< " " << dRoadBehind // 71: rear edge x-point of the tandem
				<< " " << dCt // 72: centrifugally induced virtual displacement of tire in the radial direction
				<< " " << (dRoad-dRoadPrev)/dt // 73 road-timestep calculated vertical velocity
				<< " " << dt; // 74 timestep

				out << std::endl;
			}
	}
}



void
Wheel4::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 16; // defines dimensions of work vector having X DOF
	*piNumCols = 16;
}

VariableSubMatrixHandler&
Wheel4::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(10, 16);

	// attributing rows and columns of the Jacobian matrix,

	integer iFirstIndex = iGetFirstIndex();
	integer iRingFirstPosIndex = pRing->iGetFirstPositionIndex();
	integer iWheelFirstPosIndex = pWheel->iGetFirstPositionIndex();
	integer iRingFirstMomIndex = pRing->iGetFirstMomentumIndex();
	// integer iWheelFirstMomIndex = pWheel->iGetFirstMomentumIndex();

	// patch x dir,
	WM.PutRowIndex(1, iFirstIndex + 1);
	WM.PutColIndex(1, iFirstIndex + 1);
	WM.PutRowIndex(2, iFirstIndex + 2);
	WM.PutColIndex(2, iFirstIndex + 2);
	// patch y dir,
	WM.PutRowIndex(3, iFirstIndex + 3);
	WM.PutColIndex(3, iFirstIndex + 3);
	WM.PutRowIndex(4, iFirstIndex + 4);
	WM.PutColIndex(4, iFirstIndex + 4);

	WM.PutRowIndex(5, iRingFirstMomIndex+1); // ring x
	WM.PutRowIndex(6, iRingFirstMomIndex+2); // ring y
	WM.PutRowIndex(7, iRingFirstMomIndex+3); // ring z
	WM.PutRowIndex(8, iRingFirstMomIndex+4); // ring Mom x
	WM.PutRowIndex(9, iRingFirstMomIndex+5); // ring Mom y
	WM.PutRowIndex(10, iRingFirstMomIndex+6); // ring Mom z

	WM.PutColIndex(5, iRingFirstPosIndex+1); // ring x
	WM.PutColIndex(6, iRingFirstPosIndex+2); // ring y
	WM.PutColIndex(7, iRingFirstPosIndex+3); // ring z
	WM.PutColIndex(8, iRingFirstPosIndex+4); // ring ang x
	WM.PutColIndex(9, iRingFirstPosIndex+5); // ring ang y
	WM.PutColIndex(10, iRingFirstPosIndex+6); // ring ang z
	WM.PutColIndex(11, iWheelFirstPosIndex+1); // wheel x
	WM.PutColIndex(12, iWheelFirstPosIndex+2); // wheel y
	WM.PutColIndex(13, iWheelFirstPosIndex+3); // wheel z
	WM.PutColIndex(14, iWheelFirstPosIndex+4); // wheel ang x
	WM.PutColIndex(15, iWheelFirstPosIndex+5); // wheel ang y
	WM.PutColIndex(16, iWheelFirstPosIndex+6); // wheel ang z

	// ,attributing rows and columns of the Jacobian matrix


	// contributions of Fint from patch pos to patch Jacobian,


	Vec3 dFnewdFintx = i*(fwdRingFlat(1)) + j*(fwdRingFlat(2));
	Vec3 dFnewdFinty = i*(latRingFlat(1)) + j*(latRingFlat(2));
	Vec3 dFnewdFintz = k;

	doublereal dFintxdVx = fwdRingFlat.Dot(i)*Cpatv(1);
	doublereal dFintxdXx = fwdRingFlat.Dot(i)*Kpatv(1);
	doublereal dFintxdVy = fwdRingFlat.Dot(j)*Cpatv(1);
	doublereal dFintxdXy = fwdRingFlat.Dot(j)*Kpatv(1);
	doublereal dFintxdVz = 0;
	doublereal dFintxdXz = 0;
	doublereal dFintydVy = latRingFlat.Dot(j)*Cpatv(2);
	doublereal dFintydXy = latRingFlat.Dot(j)*Kpatv(2);
	doublereal dFintydVx = latRingFlat.Dot(i)*Cpatv(2);
	doublereal dFintydXx = latRingFlat.Dot(i)*Kpatv(2);
	doublereal dFintydVz = 0;
	doublereal dFintydXz = 0;
	doublereal dFintzdVy = 0; // note that this one and the following 5 derivatives are defined to clarify the code since in reality the patch has no DoF in the z dir.
	doublereal dFintzdXy = 0;
	doublereal dFintzdVx = 0;
	doublereal dFintzdXx = 0;
	doublereal dFintzdVz = Cpatv(3);
	doublereal dFintzdXz = Kpatv(3);

	Vec3 dFnewdVx = dFnewdFintx*dFintxdVx+dFnewdFinty*dFintydVx+dFnewdFintz*dFintzdVx;
	Vec3 dFnewdXx = dFnewdFintx*dFintxdXx+dFnewdFinty*dFintydXx+dFnewdFintz*dFintzdXx;
	Vec3 dFnewdVy = dFnewdFintx*dFintxdVy+dFnewdFinty*dFintydVy+dFnewdFintz*dFintzdVy;
	Vec3 dFnewdXy = dFnewdFintx*dFintxdXy+dFnewdFinty*dFintydXy+dFnewdFintz*dFintzdXy;
	Vec3 dFnewdVz = dFnewdFintx*dFintxdVz+dFnewdFinty*dFintydVz+dFnewdFintz*dFintzdVz;
	Vec3 dFnewdXz = dFnewdFintx*dFintxdXz+dFnewdFinty*dFintydXz+dFnewdFintz*dFintzdXz;

	// TODO add threshold conditions if want close to zero velocities to work
	Vec3 dFdVx = -fwd*pMuX0->dGetP(dSr)*Fn/abs(dvax)*fwd.Dot(fwdRing); // TODO: CHECK THESE these four are negative of what they physically are... for consistency with previous derivs based on Fint seen from unrotated ring
	Vec3 dFdVy = -fwd*pMuX0->dGetP(dSr)*Fn/abs(dvax)*fwd.Dot(latRing);
	Vec3 dFdXx = fwd*pMuX0->dGetP(dSr)*Fn/abs(dvax)*fwd.Dot(pWheel->GetWCurr().Cross(n*sqrt((-n.Cross(WheelAxle)*fwdRingFlat.Dot(i)+j*fwdRingFlat.Dot(j)).Dot())));
	Vec3 dFdXy = fwd*pMuX0->dGetP(dSr)*Fn/abs(dvax)*fwd.Dot(pWheel->GetWCurr().Cross(n*sqrt((-n.Cross(WheelAxle)*latRingFlat.Dot(i)+j*latRingFlat.Dot(j)).Dot())));

	// x dir,
	WM.PutCoef(1, 1, Mpa - dCoef*(-dFnewdVx(1)+dFdVx(1)));	// -(dr_1/dVxprime+dCoef*dr_1/dVx) acting in x-dir, complete
	WM.PutCoef(1, 2, -dCoef*(-dFnewdXx(1)+dFdXx(1)));		// -(dr_1/dXxprime+dCoef*dr_1/dXx) acting in x-dir, complete
	WM.PutCoef(2, 1, -dCoef);					// -(dr_2/dVxprime+dCoef*dr_2/dVx) acting in x-dir, complete
	WM.PutCoef(2, 2, 1.);						// -(dr_2/dXxprime+dCoef*dr_2/dXx) acting in x-dir, complete
	WM.PutCoef(1, 3, dCoef*(dFnewdVy(1)-dFdVy(1)));		// -(-(dr_7/dVyprime+dCoef*dr_7/dVy) acting in y-dir, complete) taken from ring contrib
	WM.PutCoef(1, 4, dCoef*(dFnewdXy(1)-dFdXy(1)));		// -(-(dr_7/dXyprime+dCoef*dr_7/dXy) acting in y-dir, complete) taken from ring contrib
	// y dir,
	WM.PutCoef(3, 3, Mpa - dCoef*(-dFnewdVy(2)+dFdVy(2)));	// -(dr_3/dVyprime+dCoef*dr_3/dVy) acting in y-dir, complete
	WM.PutCoef(3, 4, dCoef*(dFnewdXy(2)-dFdXy(2)));		// -(dr_3/dXyprime+dCoef*dr_3/dXy) acting in y-dir, complete
	WM.PutCoef(4, 3, -dCoef);					// -(dr_4/dVyprime+dCoef*dr_4/dVy) acting in y-dir, complete
	WM.PutCoef(4, 4, 1.);						// -(dr_4/dXyprime+dCoef*dr_4/dXy) acting in y-dir, complete
	WM.PutCoef(3, 1, dCoef*(dFnewdVx(2)-dFdVx(2)));		// -(-(dr_8/dVxprime+dCoef*dr_8/dVx)) taken from ring contrib
	WM.PutCoef(3, 2, dCoef*(dFnewdXx(2)-dFdXx(2)));		// -(-(dr_8/dXxprime+dCoef*dr_8/dXx)) taken from ring contrib

	// ,contributions of Fint from patch pos to patch Jacobian

	// contributions of Fint from patch pos to ring Jacobian,
	Vec3 dFrintdFnewx = -(n.Cross(WheelAxle));
	Vec3 dFrintdFnewy = j;
	Vec3 dFrintdFnewz = n;
	Vec3 dFrintdVx = dFrintdFnewx*dFnewdVx(1)+dFrintdFnewy*dFnewdVx(2)+dFrintdFnewz*dFnewdVx(3);
	Vec3 dFrintdVy = dFrintdFnewx*dFnewdVy(1)+dFrintdFnewy*dFnewdVy(2)+dFrintdFnewz*dFnewdVy(3);
	Vec3 dFrintdVz = dFrintdFnewx*dFnewdVz(1)+dFrintdFnewy*dFnewdVz(2)+dFrintdFnewz*dFnewdVz(3);
	Vec3 dFrintdXx = dFrintdFnewx*dFnewdXx(1)+dFrintdFnewy*dFnewdXx(2)+dFrintdFnewz*dFnewdXx(3);
	Vec3 dFrintdXy = dFrintdFnewx*dFnewdXy(1)+dFrintdFnewy*dFnewdXy(2)+dFrintdFnewz*dFnewdXy(3);
	Vec3 dFrintdXz = dFrintdFnewx*dFnewdXz(1)+dFrintdFnewy*dFnewdXz(2)+dFrintdFnewz*dFnewdXz(3);



	WM.PutCoef(5, 1, -dCoef*(dFrintdVx(1)));		// -(dr_5/dVxprime+dCoef*dr_5/dVx) acting in x-dir, complete
	WM.PutCoef(5, 2, -dCoef*(dFrintdXx(1)));		// -(dr_5/dXxprime+dCoef*dr_5/dXx) acting in x-dir, complete
	WM.PutCoef(5, 3, -dCoef*(dFrintdVy(1)));		// -(dr_5/dVyprime+dCoef*dr_5/dVy) acting in y-dir, complete
	WM.PutCoef(5, 4, -dCoef*(dFrintdXy(1)));		// -(dr_5/dXyprime+dCoef*dr_5/dXy) acting in y-dir, complete
	WM.PutCoef(6, 1, -dCoef*(dFrintdVx(2)));		// -(dr_6/dVxprime+dCoef*dr_8/dVx)
	WM.PutCoef(6, 2, -dCoef*(dFrintdXx(2)));		// -(dr_6/dXxprime+dCoef*dr_8/dXx)
	WM.PutCoef(6, 3, -dCoef*(dFrintdVy(2)));		// -(dr_6/dVyprime+dCoef*dr_8/dVy)
	WM.PutCoef(6, 4, -dCoef*(dFrintdXy(2)));		// -(dr_6/dXyprime+dCoef*dr_8/dXy)
	WM.PutCoef(7, 1, -dCoef*(dFrintdVx(3)));		// -(dr_7/dVxprime+dCoef*dr_7/dVx) acting in x-dir, complete
	WM.PutCoef(7, 2, -dCoef*(dFrintdXx(3)));		// -(dr_7/dXxprime+dCoef*dr_7/dXx) acting in x-dir, complete
	WM.PutCoef(7, 3, -dCoef*(dFrintdVy(3)));		// -(dr_7/dVyprime+dCoef*dr_7/dVy) acting in y-dir, complete
	WM.PutCoef(7, 4, -dCoef*(dFrintdXy(3)));		// -(dr_7/dXyprime+dCoef*dr_7/dXy) acting in y-dir, complete




	// ,contributions of Fint from patch pos to the ring Jacobian





	// contributions of the ring position & rot mat to the jacobian of the patch,

	WM.PutCoef(1, 5, -(dFnewdVx(1)+dCoef*dFnewdXx(1)));	// -(dr_1/dXx,ring,prime+dCoef*dr_1/dXx,ring) acting in x-dir, complete
	WM.PutCoef(1, 6, -(dFnewdVy(1)+dCoef*dFnewdXy(1)));	// -(dr_1/dXy,ring,prime+dCoef*dr_1/dXy,ring) acting in x-dir, complete
	WM.PutCoef(3, 6, -(dFnewdVy(2)+dCoef*dFnewdXy(2)));	// -(dr_3/dXy,ring,prime+dCoef*dr_3/dXy,ring) acting in y-dir, complete
	WM.PutCoef(3, 5, -(dFnewdVx(2)+dCoef*dFnewdXx(2)));	// -(dr_3/dXx,ring,prime+dCoef*dr_3/dXx,ring) acting in y-dir, complete

	doublereal dFintz_dZ = -Kpatv(3);
	doublereal dFintz_dZp = -Cpatv(3);
	doublereal dFxpatch_dZ = -(fwd*dMuX*dFintz_dZ).Dot(i)-(lat*dMuY*dFintz_dZ).Dot(i);
	doublereal dFxpatch_dZp = -(fwd*dMuX*dFintz_dZp).Dot(i)-(lat*dMuY*dFintz_dZp).Dot(i);
	doublereal dFypatch_dZ = -(fwd*dMuX*dFintz_dZ).Dot(j)-(lat*dMuY*dFintz_dZ).Dot(j);
	doublereal dFypatch_dZp = -(fwd*dMuX*dFintz_dZp).Dot(j)-(lat*dMuY*dFintz_dZp).Dot(j);
	WM.PutCoef(1, 7, -dFxpatch_dZp - dCoef * dFxpatch_dZ);	// influence on the MF forces caused by the variation of the normal force
	WM.PutCoef(3, 7, -dFypatch_dZp - dCoef * dFypatch_dZ);	// influence on the MF forces caused by the variation of the normal force


// TODO move some of these vectors out of the loop to accelerate the calculation!!
	Vec3 RvrCrossn = (pRing->GetRCurr()*WheelAxle).Cross(n);
	Vec3 dfwdR_dRvx = (i.Cross(n))/sqrt(RvrCrossn.Dot())-RvrCrossn*pow(RvrCrossn.Dot(),-3/2)*((i.Cross(n)).Dot(RvrCrossn)); // d(fwd)/d(Rv)x
	Vec3 dfwdR_dRvy = (j.Cross(n))/sqrt(RvrCrossn.Dot())-RvrCrossn*pow(RvrCrossn.Dot(),-3/2)*((j.Cross(n)).Dot(RvrCrossn)); // d(fwd)/d(Rv)y
	Vec3 dfwdR_dRvz = (k.Cross(n))/sqrt(RvrCrossn.Dot())-RvrCrossn*pow(RvrCrossn.Dot(),-3/2)*((k.Cross(n)).Dot(RvrCrossn)); // d(fwd)/d(Rv)z
	Vec3 RvrCrossk = (pRing->GetRCurr()*WheelAxle).Cross(k);
	Vec3 dfwdRFlat_dRvx = (i.Cross(k))/sqrt(RvrCrossk.Dot())-RvrCrossk*pow(RvrCrossk.Dot(),-3/2)*((i.Cross(k)).Dot(RvrCrossk)); // d(fwd)/d(Rv)x
	Vec3 dfwdRFlat_dRvy = (j.Cross(k))/sqrt(RvrCrossk.Dot())-RvrCrossk*pow(RvrCrossk.Dot(),-3/2)*((j.Cross(k)).Dot(RvrCrossk)); // d(fwd)/d(Rv)y
	Vec3 dfwdRFlat_dRvz = (k.Cross(k))/sqrt(RvrCrossk.Dot())-RvrCrossk*pow(RvrCrossk.Dot(),-3/2)*((k.Cross(k)).Dot(RvrCrossk)); // d(fwd)/d(Rv)z
	Vec3 dR1_dRvr = Vec3((dfwdRFlat_dRvx.Dot(Xpar)+fwdRingFlat.Dot(-((n.Cross(dfwdR_dRvx)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvx))*dR_0))*Kpatv(1)+dfwdRFlat_dRvx.Dot(Vpar)*Cpatv(1),(dfwdRFlat_dRvy.Dot(Xpar)+fwdRingFlat.Dot(-((n.Cross(dfwdR_dRvy)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvy))*dR_0))*Kpatv(1)+dfwdRFlat_dRvy.Dot(Vpar)*Cpatv(1),(dfwdRFlat_dRvz.Dot(Xpar)+fwdRingFlat.Dot(-((n.Cross(dfwdR_dRvz)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvz))*dR_0))*Kpatv(1)+dfwdRFlat_dRvz.Dot(Vpar)*Cpatv(1));
	Vec3 dR3_dRvr = Vec3(((k.Cross(dfwdRFlat_dRvx)).Dot(Xpar)+(latRingFlat).Dot(-((n.Cross(dfwdR_dRvx)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvx))*dR_0))*Kpatv(2)+(k.Cross(dfwdRFlat_dRvx)).Dot(Vpar)*Cpatv(2),((k.Cross(dfwdRFlat_dRvy)).Dot(Xpar)+latRingFlat.Dot(-((n.Cross(dfwdR_dRvy)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvy))*dR_0))*Kpatv(2)+(k.Cross(dfwdRFlat_dRvy)).Dot(Vpar)*Cpatv(2),((k.Cross(dfwdRFlat_dRvz)).Dot(Xpar)+latRingFlat.Dot(-((n.Cross(dfwdR_dRvz)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvz))*dR_0))*Kpatv(2)+(k.Cross(dfwdRFlat_dRvz)).Dot(Vpar)*Cpatv(2));
	// NOTE: latRingFlat = k X fwdRingFlat = i(-fwdRingFlat(2)) + j(fwdRingFlat(1)), negative sign is there for clarity even though it is then canceled in the Gr equations
	Vec3 dR1new_dRvr = -Vec3(dfwdRFlat_dRvx(1)*Fint_old(1)+fwdRingFlat(1)*dR1_dRvr(1)-dfwdRFlat_dRvx(2)*Fint_old(2)+latRingFlat(1)*dR3_dRvr(1),dfwdRFlat_dRvy(1)*Fint_old(1)+fwdRingFlat(1)*dR1_dRvr(2)-dfwdRFlat_dRvy(2)*Fint_old(2)+latRingFlat(1)*dR3_dRvr(2),dfwdRFlat_dRvz(1)*Fint_old(1)+fwdRingFlat(1)*dR1_dRvr(3)-dfwdRFlat_dRvz(2)*Fint_old(2)+latRingFlat(1)*dR3_dRvr(3));
	Vec3 dR3new_dRvr = -Vec3(dfwdRFlat_dRvx(2)*Fint_old(1)+fwdRingFlat(2)*dR1_dRvr(1)+dfwdRFlat_dRvx(1)*Fint_old(2)+latRingFlat(2)*dR3_dRvr(1),dfwdRFlat_dRvy(2)*Fint_old(1)+fwdRingFlat(2)*dR1_dRvr(2)+dfwdRFlat_dRvy(1)*Fint_old(2)+latRingFlat(2)*dR3_dRvr(2),dfwdRFlat_dRvz(2)*Fint_old(1)+fwdRingFlat(2)*dR1_dRvr(3)+dfwdRFlat_dRvz(1)*Fint_old(2)+latRingFlat(2)*dR3_dRvr(3));
	Vec3 Rvr = (pRing->GetRCurr()*WheelAxle);
	Vec3 Gr1r = -dR1new_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr3r = -dR3new_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)


	WM.PutCoef(1, 8, -Gr1r(1)); // x-comp of  ring Rot Mat contrib on x dir of patch
	WM.PutCoef(1, 9, -Gr1r(2)); // y-comp of ring Rot Mat contrib on x dir of patch
	WM.PutCoef(1, 10, -Gr1r(3)); // z-comp of ring Rot Mat contrib on x dir of patch
	WM.PutCoef(3, 8, -Gr3r(1)); // x-comp of  ring Rot Mat contrib on y dir of patch
	WM.PutCoef(3, 9, -Gr3r(2)); // y-comp of ring Rot Mat contrib on y dir of patch
	WM.PutCoef(3, 10, -Gr3r(3)); // z-comp of ring Rot Mat contrib on y dir of patch
	// ,contributions of the ring position & rot mat to the jacobian of the patch


	// contributions of ring rot mat to the ring Jacobian,

	Vec3 dR5_dRvr = Vec3(-(((n.Cross(dfwdR_dRvx)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvx))*dR_0)(3)*Kpatv(3),-(((n.Cross(dfwdR_dRvy)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvy))*dR_0)(3)*Kpatv(3),-(((n.Cross(dfwdR_dRvz)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvz))*dR_0)(3)*Kpatv(3));
	Vec3 dR5new_dRvr = dR5_dRvr;
	Vec3 dRr5_dRvr = Vec3((n*(dR5new_dRvr(1))-(n.Cross(WheelAxle))*(dR1new_dRvr(1))).Dot(i),dR1new_dRvr(2),(n*(dR5new_dRvr(3))-(n.Cross(WheelAxle))*(dR1new_dRvr(3))).Dot(i)); // adjusted for inclination of slope in new equation (ignores projection correction, ie: only good for straqight line or very small steer)
//std::cout << "This should never be something other than 0: " << (n*(dR5new_dRvr(2))-(n.Cross(WheelAxle))*(dR1new_dRvr(2))).Dot(j) << std::endl;
	Vec3 dRr6_dRvr = dR3new_dRvr; // adjusted for inclination of slope in new equation (ignores steer projection correction, ie: only good for straight line or very small steer)
	Vec3 dRr7_dRvr = Vec3((n*(dR5new_dRvr(1))-(n.Cross(WheelAxle))*(dR1new_dRvr(1))).Dot(k),dR5new_dRvr(2),(n*(dR5new_dRvr(3))-(n.Cross(WheelAxle))*(dR1new_dRvr(3))).Dot(k)); // adjusted for inclination of slope in new equation (ignores projection correction, ie: only good for straqight line or very small steer)


	Vec3 Gr5r = -dRr5_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr6r = -dRr6_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr7r = -dRr7_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)

//	std::cout << "dR5_dRvr: " << dR5_dRvr << std::endl;
//	std::cout << "Gr7r: " << Gr7r << std::endl;

	WM.PutCoef(5, 8, Gr5r(1)); // x-comp of  ring Rot Mat contrib on x dir of ring
	WM.PutCoef(5, 9, Gr5r(2)); // y-comp of ring Rot Mat contrib on x dir of ring
	WM.PutCoef(5, 10, Gr5r(3)); // z-comp of ring Rot Mat contrib on x dir of ring
	WM.PutCoef(6, 8, Gr6r(1)); // x-comp of  ring Rot Mat contrib on y dir of ring
	WM.PutCoef(6, 9, Gr6r(2)); // y-comp of ring Rot Mat contrib on y dir of ring
	WM.PutCoef(6, 10, Gr6r(3)); // z-comp of ring Rot Mat contrib on y dir of ring
	WM.PutCoef(7, 8, Gr7r(1)); // x-comp of  ring Rot Mat contrib on z dir of ring
	WM.PutCoef(7, 9, Gr7r(2)); // y-comp of ring Rot Mat contrib on z dir of ring
	WM.PutCoef(7, 10, Gr7r(3)); // z-comp of ring Rot Mat contrib on z dir of ring

	Vec3 dcrcr_dRvx = ((n.Cross(dfwdR_dRvx)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvx))*dR_0;
	Vec3 dcrcr_dRvy = ((n.Cross(dfwdR_dRvy)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvy))*dR_0;
	Vec3 dcrcr_dRvz = ((n.Cross(dfwdR_dRvz)).Cross(fwdRing)+(n.Cross(fwdRing)).Cross(dfwdR_dRvz))*dR_0;
	Vec3 ddistM_dRvx = dcrcr_dRvx-n.Cross(WheelAxle)*(dfwdRFlat_dRvx(1)*Xpar(1)-fwdRingFlat(1)*dcrcr_dRvx(1)+(n.Cross(dcrcr_dRvx)).Dot(i)*Xpar(2)-latRingFlat(1)*(dfwdRFlat_dRvx(2)))-n*(dcrcr_dRvx(3))+j*(dfwdRFlat_dRvx(2)*Xpar(1)-fwdRingFlat(2)*dcrcr_dRvx(1)+(n.Cross(dcrcr_dRvx)).Dot(j)*Xpar(2)-latRingFlat(2)*(dfwdRFlat_dRvx(2)));
	Vec3 ddistM_dRvy = dcrcr_dRvy-n.Cross(WheelAxle)*(dfwdRFlat_dRvy(1)*Xpar(1)-fwdRingFlat(1)*dcrcr_dRvy(1)+(n.Cross(dcrcr_dRvy)).Dot(i)*Xpar(2)-latRingFlat(1)*(dfwdRFlat_dRvy(2)))-n*(dcrcr_dRvy(3))+j*(dfwdRFlat_dRvy(2)*Xpar(1)-fwdRingFlat(2)*dcrcr_dRvy(1)+(n.Cross(dcrcr_dRvy)).Dot(j)*Xpar(2)-latRingFlat(2)*(dfwdRFlat_dRvy(2)));
	Vec3 ddistM_dRvz = dcrcr_dRvz-n.Cross(WheelAxle)*(dfwdRFlat_dRvz(1)*Xpar(1)-fwdRingFlat(1)*dcrcr_dRvz(1)+(n.Cross(dcrcr_dRvz)).Dot(i)*Xpar(2)-latRingFlat(1)*(dfwdRFlat_dRvz(2)))-n*(dcrcr_dRvz(3))+j*(dfwdRFlat_dRvz(2)*Xpar(1)-fwdRingFlat(2)*dcrcr_dRvz(1)+(n.Cross(dcrcr_dRvz)).Dot(j)*Xpar(2)-latRingFlat(2)*(dfwdRFlat_dRvz(2)));



	Vec3 dMsa_dRvrx = n*boolFn*(-pTr->dGet(dAlpha_t)*(dR5_dRvr(1)*Fint_ring(2)+Fn*dRr6_dRvr(1))+pMzr->dGet(dAlpha_r)*dR5_dRvr(1));
	Vec3 dMsa_dRvry = n*boolFn*(-pTr->dGet(dAlpha_t)*(dR5_dRvr(2)*Fint_ring(2)+Fn*dRr6_dRvr(2))+pMzr->dGet(dAlpha_r)*dR5_dRvr(2));
	Vec3 dMsa_dRvrz = n*boolFn*(-pTr->dGet(dAlpha_t)*(dR5_dRvr(3)*Fint_ring(2)+Fn*dRr6_dRvr(3))+pMzr->dGet(dAlpha_r)*dR5_dRvr(3));

	Vec3 dRrM_dRvrx  = distM.Cross(Vec3(dRr5_dRvr(1), dRr6_dRvr(1), dRr7_dRvr(1)))+ddistM_dRvx.Cross(Fint_ring) + dMsa_dRvrx;
	Vec3 dRrM_dRvry  = distM.Cross(Vec3(dRr5_dRvr(2), dRr6_dRvr(2), dRr7_dRvr(2)))+ddistM_dRvy.Cross(Fint_ring) + dMsa_dRvry;
	Vec3 dRrM_dRvrz = distM.Cross(Vec3(dRr5_dRvr(3), dRr6_dRvr(3), dRr7_dRvr(3)))+ddistM_dRvz.Cross(Fint_ring)  + dMsa_dRvrz;

	Vec3 dRr8_dRvr  = Vec3(dRrM_dRvrx(1),dRrM_dRvry(1),dRrM_dRvrz(1));
	Vec3 dRr9_dRvr  = Vec3(dRrM_dRvrx(2),dRrM_dRvry(2),dRrM_dRvrz(2));
	Vec3 dRr10_dRvr = Vec3(dRrM_dRvrx(3),dRrM_dRvry(3),dRrM_dRvrz(3));

	Vec3 Gr8r = -dRr8_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr9r = -dRr9_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr10r = -dRr10_dRvr.Cross(Rvr)*dCoef; // dCoef*((-dr/dRv)((Rv)X)


	WM.PutCoef(8, 8, Gr8r(1) ); // x-comp of  ring Rot Mat contrib on x mom of ring
	WM.PutCoef(8, 9, Gr8r(2) ); // y-comp of ring Rot Mat contrib on x mom of ring
	WM.PutCoef(8, 10, Gr8r(3) ); // z-comp of ring Rot Mat contrib on x mom of ring
	WM.PutCoef(9, 8, Gr9r(1) ); // x-comp of  ring Rot Mat contrib on y mom of ring
	WM.PutCoef(9, 9, Gr9r(2) ); // y-comp of ring Rot Mat contrib on y mom of ring
	WM.PutCoef(9, 10, Gr9r(3) ); // z-comp of ring Rot Mat contrib on y mom of ring
	WM.PutCoef(10, 8, Gr10r(1) ); // x-comp of  ring Rot Mat contrib on z mom of ring
	WM.PutCoef(10, 9, Gr10r(2) ); // y-comp of ring Rot Mat contrib on z mom of ring
	WM.PutCoef(10, 10, Gr10r(3) ); // z-comp of ring Rot Mat contrib on z mom of ring

	// ,contributions of ring rot mat to the ring Jacobian


	// contributions of ring and patch positions to ring Momentum jacobian,

	Vec3 dM_dXx = distM.Cross(Vec3(dFrintdXx(1),dFrintdXx(2),dFrintdXx(3)))+(-n.Cross(WheelAxle)*fwdRingFlat(1)+j*fwdRingFlat(2)).Cross(Fint_ring)+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdXx(3)*Fint_ring(2)-tr*dFrintdXx(2)+pMzr->dGet(dAlpha_r)*dFnewdXx(3));
	Vec3 dM_dXy = distM.Cross(Vec3(dFrintdXy(1),dFrintdXy(2),dFrintdXy(3)))+(-n.Cross(WheelAxle)*latRingFlat(1)+j*latRingFlat(2)).Cross(Fint_ring)+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdXy(3)*Fint_ring(2)-tr*dFrintdXy(2)+pMzr->dGet(dAlpha_r)*dFnewdXy(3));
	Vec3 dM_dXz = distM.Cross(Vec3(dFrintdXz(1),dFrintdXz(2),dFrintdXz(3)))+n.Cross(Fint_ring)+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdXz(3)*Fint_ring(2)-tr*dFrintdXz(2)+pMzr->dGet(dAlpha_r)*dFnewdXz(3));
	Vec3 dM_dVx = distM.Cross(Vec3(dFrintdVx(1),dFrintdVx(2),dFrintdVx(3)))+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdVx(3)*Fint_ring(2)-tr*dFrintdVx(2)+pMzr->dGet(dAlpha_r)*dFnewdVx(3));
	Vec3 dM_dVy = distM.Cross(Vec3(dFrintdVy(1),dFrintdVy(2),dFrintdVy(3)))+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdVy(3)*Fint_ring(2)-tr*dFrintdVy(2)+pMzr->dGet(dAlpha_r)*dFnewdVy(3));
	Vec3 dM_dVz = distM.Cross(Vec3(dFrintdVz(1),dFrintdVz(2),dFrintdVz(3)))+n*boolFn*(-pTr->dGet(dAlpha_t)*dFnewdVz(3)*Fint_ring(2)-tr*dFrintdVz(2)+pMzr->dGet(dAlpha_r)*dFnewdVz(3));
	Vec3 dM_dXrx = -dM_dXx; // because ring pos is present in exactly the same way as patch pos but negative
	Vec3 dM_dXrxp = -dM_dVx; // because ring pos is present in exactly the same way as patch pos but negative
	Vec3 dM_dXry = -dM_dXy; // because ring pos is present in exactly the same way as patch pos but negative
	Vec3 dM_dXryp = -dM_dVy; // because ring pos is present in exactly the same way as patch pos but negative
	Vec3 dM_dXrz = -dM_dXz;
	Vec3 dM_dXrzp = -dM_dVz;

	WM.PutCoef(8, 1, -dCoef*dM_dVx(1));
	WM.PutCoef(8, 2, -dCoef*dM_dXx(1));
	WM.PutCoef(8, 3, -dCoef*dM_dVy(1));
	WM.PutCoef(8, 4, -dCoef*dM_dXy(1));
	WM.PutCoef(8, 5, -dCoef*dM_dXrx(1)-dM_dXrxp(1));
	WM.PutCoef(8, 6, -dCoef*dM_dXry(1)-dM_dXryp(1));
	WM.PutCoef(8, 7, -dCoef*dM_dXrz(1)-dM_dXrzp(1));
	WM.PutCoef(9, 1, -dCoef*dM_dVx(2));
	WM.PutCoef(9, 2, -dCoef*dM_dXx(2));
	WM.PutCoef(9, 3, -dCoef*dM_dVy(2));
	WM.PutCoef(9, 4, -dCoef*dM_dXy(2));
	WM.PutCoef(9, 5, -dCoef*dM_dXrx(2)-dM_dXrxp(2));
	WM.PutCoef(9, 6, -dCoef*dM_dXry(2)-dM_dXryp(2));
	WM.PutCoef(9, 7, -dCoef*dM_dXrz(2)-dM_dXrzp(2));
	WM.PutCoef(10, 1, -dCoef*dM_dVx(3));
	WM.PutCoef(10, 2, -dCoef*dM_dXx(3));
	WM.PutCoef(10, 3, -dCoef*dM_dVy(3));
	WM.PutCoef(10, 4, -dCoef*dM_dXy(3));
	WM.PutCoef(10, 5, -dCoef*dM_dXrx(3)-dM_dXrxp(3));
	WM.PutCoef(10, 6, -dCoef*dM_dXry(3)-dM_dXryp(3));
	WM.PutCoef(10, 7, -dCoef*dM_dXrz(3)-dM_dXrzp(3));


	// ,contributions of ring and patch positions to ring Momentum jacobian

	// contributions of ring pos to the ring Jacobian,

	WM.PutCoef(5, 5, dFrintdVx(1) + dCoef*dFrintdXx(1));
	WM.PutCoef(5, 6, dFrintdVy(1) + dCoef*dFrintdXy(1));
	WM.PutCoef(5, 7, dFrintdVz(1) + dCoef*dFrintdXz(1));
//	WM.PutCoef(5, 7, 0.); // contribution of z distance to x and y viscoelastic forces is null

	WM.PutCoef(6, 5, dFrintdVx(2) + dCoef*dFrintdXx(2));
	WM.PutCoef(6, 6, dFrintdVy(2) + dCoef*dFrintdXy(2));
	WM.PutCoef(6, 7, dFrintdVz(2) + dCoef*dFrintdXz(2));
//	WM.PutCoef(6, 7, 0.); // contribution of z distance to x and y viscoelastic forces is null

	WM.PutCoef(7, 5, dFrintdVx(3) + dCoef*dFrintdXx(3));
	WM.PutCoef(7, 6, dFrintdVy(3) + dCoef*dFrintdXy(3));
	WM.PutCoef(7, 7, dFrintdVz(3) + dCoef*dFrintdXz(3));

	// ,contributions of ring pos to the ring Jacobian



	// contributions of wheel position  velocity, angular velocity, and rotation matrix to jac of patch,


	// block to see where we are on the derivative and react accordingly
//		derivSign=-1  ensures the proper sign of the derivative if we are the the left side of the y-axis because we have an absolute value function to derive! This applies to both derivatives of lat and long slips
//		derivSign=1 ensures the proper sign of the derivative if we are the the right side of the y-axis because we have an absolute value function to derive! This applies to both derivatives of lat and long slips
derivSign = copysign(1, dvax);

if (std::abs(dvax) < TRH ) // max linear values of dAlpha and dSr estimated from the plots of the functions using data from Vincenzo
{
//	std::cout << "JacOFF, " << " dvax: " << dvax << " dAlpha: " << dAlpha << " dSr: " << dSr << std::endl;
	// returns a null Jacobian for the wheel contribution for velocities too close to zero or past uppermost point of the curve
	// TODO: check if lateral slip could be considered even when vel near zero...
	// TODO make condition for on/off distinguish between lateral contribs and long contribs?
}
else
{
//	std::cout << "JacON" << std::endl;



	doublereal Cmf = Fn*pMuX0->dGetP(dSr); // slope of MF long forces using deriv of force w.r.t. slip
	doublereal Cmf_a = Fn*pMuY0->dGetP(dAlpha); // slope of MF lat forces using deriv of force w.r.t. slip
	latBool = 1; // cancels the on/off switch because slope is calculated properly by ginac function
	fwdBool = 1;


	// TODO: add contribution of the rotation matrix? ie: dependence of Fint on Rring (so add 3 columns to Jacobian? Or is the rotation matrix assumed to be fixed over the iteration as pointed out by Marco if I understood correctly...?)
	Vec3 wheelW = pWheel->GetWCurr();
	Vec3 refW = pWheel->GetWRef();
	doublereal dAtan_dAlpha = 1/(1+pow((lat.Dot(va))/(fwd.Dot(va)),2)); // this deriv is actually d(atan(dvay/dvax)/d(dvay/dvax)
	doublereal dr1_dWx = derivSign*(-fwd.Dot(i.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dWx,wheel contrib of long slip MF on x-dir of patch
	doublereal dr1_dWy = derivSign*(-fwd.Dot(j.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dWy,wheel contrib of long slip MF on x-dir of patch
	doublereal dr1_dWz = derivSign*(-fwd.Dot(k.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dWz,wheel contrib of long slip MF on x-dir of patch
	doublereal dr1_dXxp = derivSign*((fwd.Dot(i))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(i))/pow((fwd.Dot(va)),2)*(-1))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dXx,prime,wheel contrib of long slip MF on x-dir of patch
	dr1_dXxp += derivSign*dAtan_dAlpha*((lat.Dot(i))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(fwd.Dot(i)))*(-Cmf_a)*lat(1)*latBool; // dr_1/dXx,prime,wheel contrib of lat slip MF on x-dir of patch
	doublereal dr1_dXyp = derivSign*((fwd.Dot(j))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(j))/pow((fwd.Dot(va)),2)*(-1))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dXy,prime,wheel contrib of long slip MF on x-dir of patch
	dr1_dXyp += derivSign*dAtan_dAlpha*((lat.Dot(j))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(fwd.Dot(j)))*(-Cmf_a)*lat(1)*latBool; // dr_1/dXy,prime,wheel contrib of lat slip MF on x-dir of patch
	doublereal dr1_dXzp = derivSign*((fwd.Dot(k))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(k))/pow((fwd.Dot(va)),2)*(-1))*fwd(1)*(-Cmf)*fwdBool; // dr_1/dXz,prime,wheel contrib of long slip MF on x-dir of patch
	dr1_dXzp += derivSign*dAtan_dAlpha*((lat.Dot(k))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(fwd.Dot(k)))*(-Cmf_a)*lat(1)*latBool; // dr_1/dXz,prime,wheel contrib of lat slip MF on x-dir of patch


	Vec3 dr1_dW(dr1_dWx,dr1_dWy,dr1_dWz);
	Vec3 dr1_dgdot(dr1_dW + refW.Cross(dr1_dW) * dCoef);
	Vec3 n_dr1_dXp = Vec3(dr1_dXxp,dr1_dXyp,dr1_dXzp)*(-1); // vec of Jac contribs from wheel local x,y,z forces to global x force on patch

	WM.PutCoef(1, 11, n_dr1_dXp(1)); // -dr_1/dXx,prime,wheel from J = -(dr_1/dXx,prime+dCoef*dr_1/dXx), acting in x-dir, complete
	WM.PutCoef(1, 12, n_dr1_dXp(2)); // -dr_1/dXy,prime,wheel from J = -(dr_1/dXy,prime+dCoef*dr_1/dXy), acting in x-dir, complete
	WM.PutCoef(1, 13, n_dr1_dXp(3)); // -dr_1/dXz,prime,wheel from J = -(dr_1/dXz,prime+dCoef*dr_1/dXz), acting in x-dir, complete


	Vec3 RvwCrossn = (pWheel->GetRCurr()*WheelAxle).Cross(n);
	Vec3 dfwdx_dRv = (i.Cross(n))/sqrt(RvwCrossn.Dot())-RvwCrossn*pow(RvwCrossn.Dot(),-3/2)*((i.Cross(n)).Dot(RvwCrossn)); // d(fwd)/d(Rv), x-comp
	Vec3 dfwdy_dRv = (j.Cross(n))/sqrt(RvwCrossn.Dot())-RvwCrossn*pow(RvwCrossn.Dot(),-3/2)*((j.Cross(n)).Dot(RvwCrossn)); // d(fwd)/d(Rv), y-comp
	Vec3 dfwdz_dRv = (k.Cross(n))/sqrt(RvwCrossn.Dot())-RvwCrossn*pow(RvwCrossn.Dot(),-3/2)*((k.Cross(n)).Dot(RvwCrossn)); // d(fwd)/d(Rv), z-comp
	Vec3 dR1_dRv = Vec3(derivSign*fwd(1)*(-Cmf)*((dfwdx_dRv.Dot(va)-dfwdx_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdx_dRv.Dot(va))), derivSign*fwd(1)*(-Cmf)*((dfwdy_dRv.Dot(va)-dfwdy_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdy_dRv.Dot(va))), derivSign*fwd(1)*(-Cmf)*((dfwdz_dRv.Dot(va)-dfwdz_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdz_dRv.Dot(va)))); //longitudinal slip forces contrib vector to x-dir of patch
//	dR1_dRv += Vec3(derivSign*lat(1)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdx_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdx_dRv.Dot(va))), 0, 0); // lateral slip contrib vector to x-dir of patch
	dR1_dRv += Vec3(derivSign*lat(1)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdx_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdx_dRv.Dot(va))), derivSign*lat(1)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdy_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdy_dRv.Dot(va))), derivSign*lat(1)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdz_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdz_dRv.Dot(va)))); // lateral slip contrib vector to x-dir of patch
	Vec3 dR3_dRv = Vec3(derivSign*fwd(2)*(-Cmf)*((dfwdx_dRv.Dot(va)-dfwdx_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdx_dRv.Dot(va))), derivSign*fwd(2)*(-Cmf)*((dfwdy_dRv.Dot(va)-dfwdy_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdy_dRv.Dot(va))), derivSign*fwd(2)*(-Cmf)*((dfwdz_dRv.Dot(va)-dfwdz_dRv.Dot(wheelW.Cross(n*dR_0)))/(fwd.Dot(va))+(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*dR_0)))/pow((fwd.Dot(va)),2)*(-1)*(dfwdz_dRv.Dot(va)))); //longitudinal slip forces contrib vector to y-dir of patch
	dR3_dRv += Vec3(derivSign*lat(2)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdx_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdx_dRv.Dot(va))), derivSign*lat(2)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdy_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdy_dRv.Dot(va))), derivSign*lat(2)*(-Cmf_a)*dAtan_dAlpha*(((n.Cross(dfwdz_dRv)).Dot(va))/(fwd.Dot(va))+(lat.Dot(va))/pow((fwd.Dot(va)),2)*(-1)*(dfwdz_dRv.Dot(va)))); // lateral slip contrib vector to y-dir of patch
	Vec3 Rv = (pWheel->GetRCurr()*WheelAxle);
	Vec3 Gr1 = -dR1_dRv.Cross(Rv)*dCoef; // dCoef*((-dr/dRv)((Rv)X) // 29.08.2012
	Vec3 Gr3 = -dR3_dRv.Cross(Rv)*dCoef; // dCoef*((-dr/dRv)((Rv)X) // simplified on 29.08.2012 as well as all Gr equations



	WM.PutCoef(1, 14, -dr1_dgdot(1)-Gr1(1)); // x-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on x dir of patch (formula changed)
	WM.PutCoef(1, 15, -dr1_dgdot(2)-Gr1(2)); // y-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on x dir of patch
	WM.PutCoef(1, 16, -dr1_dgdot(3)-Gr1(3)); // z-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on x dir of patch


	doublereal dr3_dWx = derivSign*(fwd.Dot(i.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(2)*Cmf*fwdBool; // dr_3/dWx,wheel on y-dir of patch
	doublereal dr3_dWy = derivSign*(fwd.Dot(j.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(2)*Cmf*fwdBool; // dr_3/dWy,wheel on y-dir of patch
	doublereal dr3_dWz = derivSign*(fwd.Dot(k.Cross(n*dR_0)))/(fwd.Dot(va))*fwd(2)*Cmf*fwdBool; // dr_3/dWz,wheel on y-dir of patch
	doublereal dr3_dXxp = derivSign*((fwd.Dot(i))/(fwd.Dot(va))-(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(i))/pow((fwd.Dot(va)),2))*fwd(2)*(-Cmf)*fwdBool; // dr_3/dXx,prime,wheel contrib of long slip MF on y-dir of patch
	dr3_dXxp += derivSign*dAtan_dAlpha*((lat.Dot(i))/(fwd.Dot(va))-(lat.Dot(va))/pow((fwd.Dot(va)),2)*(fwd.Dot(i)))*(-Cmf_a)*lat(2)*latBool; // dr_3/dXx,prime,wheel contrib of lat slip MF on y-dir of patch
	doublereal dr3_dXyp = derivSign*((fwd.Dot(j))/(fwd.Dot(va))-(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(j))/pow((fwd.Dot(va)),2))*fwd(2)*(-Cmf)*fwdBool; // dr_3/dXy,prime,wheel contrib of long slip MF on y-dir of patch
	dr3_dXyp += derivSign*dAtan_dAlpha*((lat.Dot(j))/(fwd.Dot(va))-(lat.Dot(va))/pow((fwd.Dot(va)),2)*(fwd.Dot(j)))*(-Cmf_a)*lat(2)*latBool; // dr_3/dXy,prime,wheel contrib of lat slip MF on y-dir of patch
	doublereal dr3_dXzp = derivSign*((fwd.Dot(k))/(fwd.Dot(va))-(fwd.Dot(va)-fwd.Dot(wheelW.Cross(n*(dR_0))))*(fwd.Dot(k))/pow((fwd.Dot(va)),2))*fwd(2)*(-Cmf)*fwdBool; // dr_3/dXz,prime,wheel contrib of long slip MF on y-dir of patch
	dr3_dXzp += derivSign*dAtan_dAlpha*((lat.Dot(k))/(fwd.Dot(va))-(lat.Dot(va))/pow((fwd.Dot(va)),2)*(fwd.Dot(k)))*(-Cmf_a)*lat(2)*latBool; // dr_3/dXz,prime,wheel contrib of lat slip MF on y-dir of patch


	Vec3 dr3_dW(dr3_dWx,dr3_dWy,dr3_dWz);
	Vec3 dr3_dgdot(dr3_dW + refW.Cross(dr3_dW) * dCoef);
	Vec3 n_dr3_dXp = Vec3(dr3_dXxp,dr3_dXyp,dr3_dXzp)*(-1); // vec of Jac contribs from wheel local x,y,z forces to global y force on patch

 	WM.PutCoef(3, 11, n_dr3_dXp(1)); // -dr_3/dXx,prime,wheel from J = -(dr_3/dXx,prime+dCoef*dr_3/dXx), acting in y-dir, complete
	WM.PutCoef(3, 12, n_dr3_dXp(2)); // -dr_3/dXy,prime,wheel from J = -(dr_3/dXy,prime+dCoef*dr_3/dXy), acting in y-dir, complete
	WM.PutCoef(3, 13, n_dr3_dXp(3)); // -dr_3/dXz,prime,wheel from J = -(dr_3/dXz,prime+dCoef*dr_3/dXz), acting in y-dir, complete
	WM.PutCoef(3, 14, -dr3_dgdot(1)-Gr3(1)); // x-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on y dir of patch
	WM.PutCoef(3, 15, -dr3_dgdot(2)-Gr3(2)); // y-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on y dir of patch
	WM.PutCoef(3, 16, -dr3_dgdot(3)-Gr3(3)); // z-comp of {refW X (dr_1/dW,wheel)} + Rot Mat contrib on y dir of patch


// NOTE: the wheel variables have no direct influence on the residual of the ring, thus r5-r10 are null for columns 11-16


// ,contributions of wheel position  velocity, angular velocity, and rotation matrix to jac of patch


	// contributions of wheel velocity and rot matrix to Jac of moment on ring,

	//TODO: sub the following in the code above where it occurs in the long form
	Vec3 dAlpha_dVw = (n.Cross(fwd))/(fwd.Dot(va))-fwd*(lat.Dot(va))/pow(fwd.Dot(va),2); // this deriv is d(dvay/dvax)/d(va x,y,z)
	Vec3 dM_dAlpha = n*boolFn*(-pTr->dGetP(dAlpha_t)*pow(sec(dAlpha),2)*dAtan_dAlpha*Fn*Fint_ring(2)+pMzr->dGetP(dAlpha_r)*pow(sec(dAlpha),2)*dAtan_dAlpha*Fn); // this deriv is d(Mz,vector)/d(dvay/dvax)
	Vec3 dM_dVwx = dM_dAlpha*dAlpha_dVw(1);
	Vec3 dM_dVwy = dM_dAlpha*dAlpha_dVw(2);
	Vec3 dM_dVwz = dM_dAlpha*dAlpha_dVw(3); // TODO could replace by a matrix...

	WM.PutCoef(8, 11, -dM_dVwx(1)); // -dr_8/dXx,prime,wheel from J = -(dr_8/dXx,prime+dCoef*dr_8/dXx)
	WM.PutCoef(8, 12, -dM_dVwy(1)); // -dr_8/dXy,prime,wheel from J = -(dr_8/dXy,prime+dCoef*dr_8/dXy)
	WM.PutCoef(8, 13, -dM_dVwz(1)); // -dr_8/dXz,prime,wheel from J = -(dr_8/dXz,prime+dCoef*dr_8/dXz)
	WM.PutCoef(9, 11, -dM_dVwx(2)); // -dr_9/dXx,prime,wheel from J = -(dr_9/dXx,prime+dCoef*dr_9/dXx)
	WM.PutCoef(9, 12, -dM_dVwy(2)); // -dr_9/dXy,prime,wheel from J = -(dr_9/dXy,prime+dCoef*dr_9/dXy)
	WM.PutCoef(9, 13, -dM_dVwz(2)); // -dr_9/dXz,prime,wheel from J = -(dr_9/dXz,prime+dCoef*dr_9/dXz)
	WM.PutCoef(10, 11, -dM_dVwx(3)); // -dr_10/dXx,prime,wheel from J = -(dr_10/dXx,prime+dCoef*dr_10/dXx)
	WM.PutCoef(10, 12, -dM_dVwy(3)); // -dr_10/dXy,prime,wheel from J = -(dr_10/dXy,prime+dCoef*dr_10/dXy)
	WM.PutCoef(10, 13, -dM_dVwz(3)); // -dr_10/dXz,prime,wheel from J = -(dr_10/dXz,prime+dCoef*dr_10/dXz)


	Vec3 dAlpha_dRv = Vec3((va.Dot(n.Cross(dfwdx_dRv)))/(fwd.Dot(va))-va.Dot(dfwdx_dRv)*(lat.Dot(va))/pow(va.Dot(fwd),2),(va.Dot(n.Cross(dfwdy_dRv)))/(fwd.Dot(va))-va.Dot(dfwdy_dRv)*(lat.Dot(va))/pow(va.Dot(fwd),2),(va.Dot(n.Cross(dfwdz_dRv)))/(fwd.Dot(va))-va.Dot(dfwdz_dRv)*(lat.Dot(va))/pow(va.Dot(fwd),2)); // this deriv is d(dvay/dvax)/d(Rv x,y,z)
	Vec3 dM_dRvx = dM_dAlpha*dAlpha_dRv(1);
	Vec3 dM_dRvy = dM_dAlpha*dAlpha_dRv(2);
	Vec3 dM_dRvz = dM_dAlpha*dAlpha_dRv(3); // TODO could replace by a matrix...
	Vec3 dR8_dRv  = Vec3(dM_dRvx(1),dM_dRvy(1),dM_dRvz(1));
	Vec3 dR9_dRv  = Vec3(dM_dRvx(2),dM_dRvy(2),dM_dRvz(2));
	Vec3 dR10_dRv = Vec3(dM_dRvx(3),dM_dRvy(3),dM_dRvz(3));


	Vec3 Gr8 = -dR8_dRv.Cross(Rv)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr9 = -dR9_dRv.Cross(Rv)*dCoef; // dCoef*((-dr/dRv)((Rv)X)
	Vec3 Gr10 = -dR10_dRv.Cross(Rv)*dCoef; // dCoef*((-dr/dRv)((Rv)X)


	WM.PutCoef(8, 14,  Gr8(1)); // Wheel rot mat contrib on x moment on ring
	WM.PutCoef(8, 15,  Gr8(2)); // Wheel rot mat contrib on x moment on ring
	WM.PutCoef(8, 16,  Gr8(3)); // Wheel rot mat contrib on x moment on ring
	WM.PutCoef(9, 14,  Gr9(1)); // Wheel rot mat contrib on y moment on ring
	WM.PutCoef(9, 15,  Gr9(2)); // Wheel rot mat contrib on y moment on ring
	WM.PutCoef(9, 16,  Gr9(3)); // Wheel rot mat contrib on y moment on ring
	WM.PutCoef(10, 14, Gr10(1)); // Wheel rot mat contrib on z moment on ring
	WM.PutCoef(10, 15, Gr10(2)); // Wheel rot mat contrib on z moment on ring
	WM.PutCoef(10, 16, Gr10(3)); // Wheel rot mat contrib on z moment on ring


	// ,contributions of wheel velocity and rot matrix to Jac of moment on ring


	// contributions of ring rot matrix to Jac of moment on ring,



	// ,contributions of ring rot matrix to Jac of moment on ring

}
	for (int i=1; i <= WM.iGetNumCols(); i++) {
		WM.PutCoef(1, i, WM.dGetCoef(1, i) / Kpatv(1));
		WM.PutCoef(3, i, WM.dGetCoef(3, i) / Kpatv(2));
	}
	return WorkMat;
}




SubVectorHandler& 
Wheel4::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{


		if (bFirstAC) {
			oldTime = curTime; // in if to avoid changing oldTime when only time step length changed
			dRoadPrevPrev = dRoadPrev;
			dRoadPrev = dRoad;
			dXxProjPrev = dXxProj;
			dLsProjPrev = dLsProj;
//			dtPrev = dt;
			pcRingPrev = pcRing;
			VpcRingPrevPrev = VpcRingPrev;
			VpcRingPrev = VpcRing;
			RingRadPrev = RingRad;
			nPrev = n;
			fwdRingPrev = fwdRing;
			FintPrevPrev = FintPrev;
			FintPrev = Fint;
			XparpPrevPrev = XparpPrev;
			XparpPrev = Xparp;
			XparPrevPrev = XparPrev;
			XparPrev = Xpar;
			XpaPrevPrev = XpaPrev;
			XpaPrev = Xpa;
			dXxPrev = dXx; // where dXx is the position degree of freedom of the patch in the x-dir
			deltaPrev = dR_0 - dEffRad;
		    if (deltaPrev < 0.)
		    	{
		    	deltaPrev = 0.;
//		    	std::cout << "Note for wheel: " << GetLabel() << ", the tire deflection is negative (ie: the road and tire should have lost contact in reality but the model forced it." << std::endl;
		    	}

//			std::cout << "bFirstAC" << std::endl;
	     }
	     bFirstAC = false;

	     if (bFirstAP) {
				curTime = tdc.dGet();
				dt = curTime - oldTime;
			 // variable normal direction,
		//  using previous timestep data
 		dLs = dR_a1*dR_0*(deltaPrev/dR_0 + dR_a2*sqrt(deltaPrev/dR_0)) * dPls; // multiplied by dPls as stated in Schmeitz. Remainder of equation from from Besslink eq. 4.85 (also shown in Pac2006 eq A3.10) NOTE: Besselink uses "vertical" deflection, but I use the actual deflection in the normal direction because Besselink only considers flat road while I consider unevennesses
	    if (dLs < TdLs) dLs = TdLs;
//	    RDLC = 0;
//	    if (RDL != 0) RDLC = (dXxPrev - RDB)/RDL; // this was changed to a more "proper" coding technique
		while (dXxProj >= RDL) dXxProj-=RDL;

		if (dXxPrev < RDA)
		{
		dXxProj = 0.;
		dRoadAhead = dRoadBehind = dRoadInitial;
		dLsProj = dLs;
		}
		else if (dXxPrev >= RDA && dXxPrev < RDB)
		{
//16.08.2012			dXxProj = pcRingPrev(1)-(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i))*XparPrev(1) - dXxPrev; // projected sloped x-pos of the patch
//16.08.2012			dXxProj = pcRingPrev(1)+(fwdRingPrev.Dot(i))*XparPrev(1) - dXxPrev; // projected sloped x-pos of the patch
			dXxProj = 0; // projected sloped x-pos of the patch
//			dLsProj = -(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i))*dLs;
			dLsProj = (fwdRingPrev.Dot(i))*dLs;
			dRoadAhead = dRoadInitial + ((dXxPrev-RDA)/(RDB-RDA))*(pRoad->dGet(dXxProj+dLsProj)-dRoadInitial); // projected sloped distance in the absolute x-plane
			dRoadBehind = dRoadInitial + ((dXxPrev-RDA)/(RDB-RDA))*(pRoad->dGet(dXxProj-dLsProj)-dRoadInitial); // more or less precise smoothing, only for initial stabilization of the vehicle
		}
		else if (dXxPrev >= RDB)
		{
//			dXxProj = pcRingPrev(1)-(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i))*XparPrev(1) - RDB - RDL*RDLC; // projected sloped x-pos of the patch
//16.08.2012			dXxProj = pcRingPrev(1)-(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i))*XparPrev(1) - RDB; // projected sloped x-pos of the patch
//16.08.2012			dXxProj = pcRingPrev(1)+(fwdRingPrev.Dot(i))*XparPrev(1) - RDB; // projected sloped x-pos of the patch
			dXxProj = dXxPrev - RDB; // projected sloped x-pos of the patch
//			dLsProj = -(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i))*dLs;
			dLsProj = (fwdRingPrev.Dot(i))*dLs;
			doublereal dXxProjAhead = dXxProj+dLsProj;
			doublereal dXxProjBehind = dXxProj-dLsProj;
			dXxProj = CapLoop(dXxProj);
			dXxProjAhead = CapLoop(dXxProjAhead);
			dXxProjBehind = CapLoop(dXxProjBehind);
			dRoadAhead = pRoad->dGet(dXxProjAhead);
			dRoadBehind = pRoad->dGet(dXxProjBehind);
		}
			dRoad = (dRoadAhead+dRoadBehind)/2.; // road height in the absolute ref frame

			n = (Vec3(dLsProj*2,0,dRoadAhead-dRoadBehind)).Cross(WheelAxle); // reads road file and applies a two-point follower in order to get the normal...
			dn = n.Dot();
//			std::cerr << dn << " " << dXxPrev << std::endl; //TODO: check dRoad influence on patch position/ normal / why patch moves forward when encountering a bump?
			n /= sqrt(dn);
			if (fabs(n(3)) < std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("Wheel4(" << GetLabel() << "): "
					"a road segment is vertical, "
					"wheel4 cannot properly deal with such a road" << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			// ,variable normal direction


			/////////////// calculate recommended maximum timestep for "next" step, ///////////////////
			if (dt_On)
			{
				doublereal dt_RDA = 0.5;
				if (dXxPrev >= RDB)
					{
	//				std::cout << curTime << "LOOP\n";

					doublereal dDebug2; //TEMP TODO DELETE
					dDebug2 = -std::numeric_limits<doublereal>::max();

					// are those two really needed? yes only if the vertical oscillations of the wheels and all the other members movements have lower frequencies...
					//doublereal dtMaxWheel = 1/(ang vel wheel *2);
					// doublereal dtMaxRing =

						doublereal dt_distAhead = 0;
						for (int iCnt = 0; iCnt <= ceil(log(dt_maxstep/dt_minstep)/log(dt_maxF)); iCnt++)
						{
							dt_distAhead += (dXxProj-dXxProjPrev)/pow(dt_maxF,iCnt);
						}
						dt_numAhead = ceil(dt_distAhead/dt_Res);
	//					std::cout << "dt: " << dt << " dt_distAhead: " << dt_distAhead << std::endl;

						dt_fNow = 1e-9; // initializing to the minimum value that the current timestep should be divided by
	//					dLsProjRatio = -(nPrev.Cross(pRing->GetRPrev()*WheelAxle).Dot(i));
						for (int iCnt = 1; iCnt <= dt_numAhead; iCnt++) {

							doublereal dt_maxHeight = -std::numeric_limits<doublereal>::max();
							doublereal dt_minHeight = std::numeric_limits<doublereal>::max();
							doublereal dt_dXxProj_dLsNow = dXxProj+(iCnt-1)*dt_Res - 0.35*dR_0*dPls ; // Besselink, p. 131, half contact patch length never really goes above 0.35*dR_0
							while (dt_dXxProj_dLsNow <= 2*dXxProj-dXxProjPrev+iCnt*dt_Res + 0.35*dR_0*dPls)
							{
								dt_maxHeight = fmax( dt_maxHeight, pRoad->dGet(CapLoop(dt_dXxProj_dLsNow)) );
								dt_minHeight = fmin( dt_minHeight, pRoad->dGet(CapLoop(dt_dXxProj_dLsNow)) );
								dt_dXxProj_dLsNow += dt_Res;
							}
	//						doublereal dt_stepHeight = dt_maxHeight - dt_minHeight;
							doublereal dt_stepHeight = fmax( dt_maxHeight - dt_minHeight ,1e-15);

							dt_adjFactor = dt_stepHeight/dt_maxH; // factor by which the current time step is too big for the bump to come in iCnt times the previous step distance

							dDebug2 = fmax( dDebug2, dt_adjFactor );

	//						int dt_nStepsNow = fmax(0,ceil(log(dt_adjFactor)/log(dt_maxF))); // this float max could be an integer function instead
							int dt_nStepsNow = fmax(0,ceil(log(dt_maxstep/dt_minstep)/log(dt_maxF))); // this float max could be an integer function instead
							doublereal dt_distReqd = 0;
	//						while (dt_distReqd < iCnt*dt_Res)
	//						{
								for (int jCnt = 1; jCnt <= dt_nStepsNow; jCnt++)
								{
									dt_distReqd += (dXxProj-dXxProjPrev)/pow(dt_maxF,jCnt);
								}
	//						}
							doublereal dt_fNowTmp; // "currently searched ahead position" adjustment factor (divider) wanted on timestep
							if (dt_distReqd <= dXxProjPrev+iCnt*dt_Res - 0.35*dR_0*dPls - dXxProj - (dXxProj-dXxProjPrev) && dt_adjFactor > 1.)
								{
									dt_fNowTmp = 1.;
	//								std::cout << "if 1" << std::endl;

								} else if ( dt_adjFactor > 1.)
								{
									doublereal dt_stepInflRatio = 0.;
									for (int iCnt = 1; iCnt <= dt_numAhead; iCnt++) {  // this loop checks if a smaller timestep guarantees a smaller bump
										dt_maxHeight = -std::numeric_limits<doublereal>::max();
										dt_minHeight = std::numeric_limits<doublereal>::max();
										dt_dXxProj_dLsNow = dXxProj+(iCnt-1)*dt_Res - 0.35*dR_0*dPls ; // Besselink, p. 131, half contact patch length never really goes above 0.35*dR_0
										while (dt_dXxProj_dLsNow <= dXxProj+(dXxProj-dXxProjPrev)/dt_adjFactor+iCnt*dt_Res + 0.35*dR_0*dPls)
										{
											dt_maxHeight = fmax( dt_maxHeight, pRoad->dGet(CapLoop(dt_dXxProj_dLsNow)) );
											dt_minHeight = fmin( dt_minHeight, pRoad->dGet(CapLoop(dt_dXxProj_dLsNow)) );
											dt_dXxProj_dLsNow += dt_Res;
										}
										dt_stepInflRatio = (dt_maxHeight - dt_minHeight)/dt_stepHeight; // this dt_stepHeight is the one calculated previously
									}
									if ( dt_stepInflRatio <= 1/dt_adjFactor )
									{
									dt_fNowTmp = fmin(pow(dt_adjFactor,1./dt_nStepsNow),dt_maxF);
									} else
									{
										dt_fNowTmp = dt_maxF;
									}
	//								std::cout << "if 2" << std::endl;
								} else
								{
									dt_fNowTmp = fmax(dt_minF,dt_adjFactor);

								}

	//						dt_fNowTmp = fmin(pow(dt_adjFactor,1./iCnt),dt_maxF); // allowing a maximum instant timestep change factor of dt_maxF
							dt_fNow = fmax(dt_fNow,dt_fNowTmp);

						}
						dDebug = dDebug2;
					}
				else if (dXxPrev <= dt_RDA) // dt_RDA: TOTO: make user input to control initially kept small timestep
				{
					doublereal dt_stepHeight = (pRoad->dGet(dXxProj)-dRoadInitial)*(XpaPrev(1)-XpaPrevPrev(1))/(RDB-RDA); // rough estimate of step height at each timestep for the initial buffer zone
					dt_adjFactor = fmax(dt_stepHeight/dt_maxH,1e-9); // factor by which the current time step is too big for the bump to come
					dt_fNow = fmin(dt_adjFactor,dt_maxF); // allowing a maximum instant timestep change factor of dt_maxF
					doublereal dt_RDAmin = 0.9999;
					dt_fNow = fmax(dt_fNow,dt_RDAmin); //
	//				dt_fNow = 1;
				}
				else
				{
					doublereal dt_stepHeight = (pRoad->dGet(dXxProj)-dRoadInitial)*(XpaPrev(1)-XpaPrevPrev(1))/(RDB-RDA); // rough estimate of step height at each timestep for the initial buffer zone
					dt_adjFactor = fmax(dt_stepHeight/dt_maxH,1e-9); // factor by which the current time step is too big for the bump to come
					dt_fNow = fmin(dt_adjFactor,dt_maxF); // allowing a maximum instant timestep change factor of dt_maxF
				}

				if (dt_fNow < dt_minF) dt_fNow = dt_minF; // allow minimum factor by which to divide the current timestep to be dt_minF


				bool b_dtInc = 1; // if true the timestep is allowed to increase (if road profile also allows it)
				doublereal dt_divFT = dt_divF; // initializes the dividing factor for the force induced timestep change

				doublereal dt_minMag = 0.001;
				for (int iCnt = 1; iCnt <= 3; iCnt++) {
					for (int jCnt = 1; jCnt <dt_minStepsCycle; jCnt++) {
						FintSignCk[jCnt-1][iCnt-1] = FintSignCk[jCnt][iCnt-1];
						FintSignCk[jCnt-1][iCnt-1+3] = FintSignCk[jCnt][iCnt-1+3];
					}
					FintSignCk[dt_minStepsCycle-1][iCnt-1] = sign(FintPrev(iCnt)-FintPrevPrev(iCnt)); // sign check
					FintSignCk[dt_minStepsCycle-1][iCnt-1+3] = fabs(FintPrev(iCnt)-FintPrevPrev(iCnt)); // magnitude check
					int signChg = 0;
					for (int jCnt = 2; jCnt <dt_minStepsCycle; jCnt++) {
					if ( FintSignCk[jCnt-2][iCnt-1] != FintSignCk[jCnt-1][iCnt-1] && fabs( FintSignCk[jCnt-2][iCnt-1+3] - FintSignCk[jCnt-1][iCnt-1+3] ) >= dt_minMag   ) {
						signChg++;
					}
					}
					if (signChg >=2) b_dtInc = 0;
					if (signChg ==3) dt_divFT = dt_divF3;
					if (signChg >=4) dt_divFT = dt_divF4;
				}



				if (b_dtInc)
					{
					dtMax = dt/dt_fNow;
	//				std::cout << "step allowed to increase, time: " << curTime << std::endl;
					} else
					{
						dtMax = dt/fmax(dt_fNow,dt_divFT); // reduce timestep if force cycle is underresolved





					}
			}
			else
			{
				dtMax = 0;
			}

			/////////////// ,calculate recommended maximum timestep for "next" step ///////////////////


	     }
	     bFirstAP = false;


//	     std::cerr << curTime << " " << dt;


			dKpa = pKpa->dGet();
			Kpatv = Kpa*dKpa; // time variant Kpa
			dCpa = pCpa->dGet();
			Cpatv = Cpa*dCpa;

//			std::cout << "dt: " << dt << " curTIme: " << curTime <<  " n: " << n << " oldTime: " << oldTime << "\n";

	// resize residual
		integer iNumRows = 0;
		integer iNumCols = 0;
		WorkSpaceDim(&iNumRows, &iNumCols);
		WorkVec.ResizeReset(iNumRows);
		// integer iWheelFirstMomIndex = pWheel->iGetFirstMomentumIndex();
		integer iRingFirstMomIndex = pRing->iGetFirstMomentumIndex();
		// equations indexes
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
//			WorkVec.PutRowIndex(10 + iCnt, iWheelFirstMomIndex + iCnt); //local, global indices
			WorkVec.PutRowIndex(4 + iCnt, iRingFirstMomIndex + iCnt);
		}


		integer iFirstIndex = iGetFirstIndex();

		for (int iCnt = 1; iCnt <= 4; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}


		if (firstRes)
	{
		firstRes =  false;
		SetInitialValue(const_cast <VectorHandler&>(XCurr));
	}

	// in x dir,
	dVPx = XPrimeCurr(iFirstIndex + 1);
	dXPx = XPrimeCurr(iFirstIndex + 2);
	dVx = XCurr(iFirstIndex + 1);
	dXx = XCurr(iFirstIndex + 2);
	// in y dir,
	dVPy = XPrimeCurr(iFirstIndex + 3);
	dXPy = XPrimeCurr(iFirstIndex + 4);
	dVy = XCurr(iFirstIndex + 3);
	dXy = XCurr(iFirstIndex + 4);


		// "forward" direction: axle cross normal to ground
		fwd = (pWheel->GetRCurr()*WheelAxle).Cross(n);
		doublereal d = fwd.Dot();
		if (d < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Wheel4(" << GetLabel() << "): "
				"wheel axle is (nearly) orthogonal "
				"to the ground" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		fwd /= sqrt(d);

		// "lateral" direction: normal to ground cross forward
		lat = n.Cross(fwd); // cross product


		fwdRing = (pRing->GetRCurr()*WheelAxle).Cross(n); // expresses local axle of the ring into global ref frame and then does cross prod with normal to get forward direction
		doublereal dRing = fwdRing.Dot();
		if (dRing < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Wheel4(" << GetLabel() << "): "
				"ring axle is (nearly) orthogonal "
				"to the ground" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		fwdRing /= sqrt(dRing);
		// "lateral" direction: normal to ground cross forward
		latRing = n.Cross(fwdRing); // cross product




	Xring = pRing->GetXCurr(); // absolute pos of ring axle
//29.08.2012	RingRad = latRing.Cross(fwdRing)*dR_0;// * (1 + fabs(VpcRing(1)*0.0001) ); // ring radius vector
	RingRad = -n*dR_0;// * (1 + fabs(VpcRing(1)*0.0001) ); // ring radius vector
//16.08.2012	VpcRing = (RingRad - RingRadPrev)/dt;

	pcRing = Xring + RingRad; // accounts for ring camber; position of ring "contact" point (NOT patch contact point!) absolute ref frame
	// we assume the contact point lies at the root of the normal
	// to the ground that passes through the center of the wheel


		if (curTime == 0.0) // initializes velocity and positions but gets position in z from the road driver
		{
			Vpa = VpaPrev;
			Xpa = XpaPrev;
		} else {
//			Vpa = Vec3(dVx,dVy,(dRoad-dRoadPrev)/dt); // this method calculated the patch speed explicitly and was changed for the "more implicit" method below
			Vpa = Vec3(dVx,dVy,(n(1))/(n(3))*-dVx); // TODO: jacobian for Vpa(3) ?
//			std::cout << curTime << "WheelMod " << GetLabel() << "using dt: " << (dRoad-dRoadPrev)/dt << " dVx: " <<  dVx << " --- " << Vpa(3) << std::endl;
			Xpa = Vec3(dXx,dXy,dRoad);
		}
	Vec3 XparOld = Xpa - pcRing; // patch rel. to ring, expressed in ABS ref frame
	Xpar = Xpa - pcRing; // patch rel. to ring, expressed in ABS ref frame
//	Xpar = Xpa - Xring + k*dR_0; // trying new nov 16 2011

	Vpar = Vpa - pRing->GetVCurr();// - VpcRing;// + VpcRing * (dR_0-ddistM) * 8*16*16 / ( fabs( (-( pRing->GetWCurr()).Cross(distM)).Dot(fwd) ) + 0.01 );// + VpcRing*15;// - VpcRing; // absolute ref frame
//	std::cout << VpcRing << std::endl;



// The following "Flat" direction vectors are used to properly model the ring rotation for the purpose of calculating spring elements forces btw ring and patch without actually rotating in the direction of slope which is accomplished by the projection of the forces acting on the ring, ie: Xpar(1) is the force in x of unrotated tire regardless of the slope
		fwdRingFlat = (pRing->GetRCurr()*WheelAxle).Cross(k); // expresses local axle of the ring into global ref frame and then does cross prod with normal to get forward direction
		doublereal dRingFlat = fwdRingFlat.Dot();
		if (dRingFlat < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Wheel4(" << GetLabel() << "): "
				"ringFlat axle is (nearly) orthogonal "
				"to the ground" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		fwdRingFlat /= sqrt(dRingFlat);
		// "lateral" direction: normal to ground cross forward
		latRingFlat = k.Cross(fwdRingFlat); // cross product

// here the absolute relative position and velocity of the patch are projected along the local components of the ring direction in order to calculate the forces using the proper spring and damper properties
		Xparp = Vec3(fwdRingFlat.Dot(Xpar),latRingFlat.Dot(Xpar),Xpar(3)); // projected relative displacement vector RING REF FRAME (but not with road slope)
//		Vec3 XparpOld = Vec3(fwdRingFlat.Dot(XparOld),latRingFlat.Dot(XparOld),XparOld(3)); // projected relative displacement vector RING REF FRAME (but not with road slope)
		Vparp = Vec3(fwdRingFlat.Dot(Vpar),latRingFlat.Dot(Vpar),Vpar(3)); // projected relative displacement vector RING REF FRAME (but not with road slope)


	// the "NEW" Xpar, gives the real dist in the absolute frame (because Xpar is absolute but considered in a rotated frame that rotates with the ring) (but not with road slope)



//	distM = RingRad+n*Xparp(3)-n.Cross(WheelAxle)*Xparp(1)+j*Xparp(2); // vector distance btwn ring center and patch // TODO: correct jacobian if needed
//	distM = (RingRad + Xparp ) *  sqrt(VpcRing.Dot()) * (dR_0-ddistM)*3200/(-( pRing->GetWCurr()).Cross(distM)); // vector distance btwn ring center and patch // TODO: correct jacobian if needed
//16.08.2012	distM = (RingRad + Xparp); // vector distance btwn ring center and patch // TODO: correct jacobian if needed
	distM = (RingRad + Xpar); // vector distance btwn ring center and patch // TODO: correct jacobian if needed
	ddistM = sqrt(distM.Dot());

	if (ddistM < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Wheel4(" << GetLabel() << "): "
			"distance between patch and ring center is (nearly) zero, adjust your elements or initial distance "
			 << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

// centripetal tire force due to velocity
//		std::cout << "Private data 3 :" << pRingB->dGetPrivData(3) << std::endl;
			doublereal Fcent_r = rRatio*(pRingB->dGetPrivData(3)+Mpa)*pow(fwdRing.Dot(pRing->GetWCurr().Cross(n)),2)*dR_0; // m*V^2/R using ring radius and V about ring center
//			doublereal Fcent_p = Mpa*pow(fwdRing.Dot(pRing->GetWCurr().Cross(-distM)+Vpar),2)/ddistM; // using patch to ring radius and V of patch // REMOVED because patch does not rotate, it only represents the local linear movement of the contact zone...
//			dCt = (Fcent_r+Fcent_p)/Krz + Fcent_p/Kpatv(3); // desired radial growth // adjusted after removing Fcent_p
			dCt = Fcent_r/Krz; // desired radial growth
			Fcent = Kpatv(3)*dCt; // force to apply on ring and patch in opposite dirs

	Fint = Vec3(Xparp.EBEMult(Kpatv))+Vec3(Vparp.EBEMult(Cpatv)) + k*Fcent; // rotated in the ring reference frame (but not with road slope) //TODO: jacobian for Fcent
	Fint_old = Fint; // Fint_old is the local forces expressed in the ring ref frame
	// the "NEW" Fint, AKA Fnew,
// simpler below //	Fint = i*( (fwdRingFlat*Fint(1) + latRingFlat*Fint(2)).Dot(i) ) + j*( (fwdRingFlat*Fint(1) + latRingFlat*Fint(2)).Dot(j) ) + k*(Fint(3)); // rotated back into absolute ref frame
	Fint = fwdRingFlat*Fint(1) + latRingFlat*Fint(2)  + k*(Fint(3)); // rotated back into absolute ref frame
	Fint_ring = n*(Fint(3))+j*(Fint(2))-(n.Cross(WheelAxle))*Fint(1); // artificial contact patch rotation due to the slope of the ground

		// relative speed between wheel axle and ground
		va = pWheel->GetVCurr(); // ground assumed fixed

	// disabled because wrongfully assumed that Vpa is in ring ref frame but it is in absolute...	VparWheel = fwdRing*Vpa(1)+latRing*Vpa(2)-va; // absolute ref frame.. should be in wheel ref frame
//		VparWheel = Vpa-va; // distorded absolute ref frame..
//		Vec3 veloM = n*Vparp(3)-n.Cross(WheelAxle)*Vparp(1)+j*Vparp(2);
//		Vec3 veloM = n*Vparp(3)-n.Cross(WheelAxle)*Vparp(1)+j*Vparp(2);
		VparWheel = pRing->GetVCurr()-va+Vpar;//veloM; changed...! // relative velocity between wheel and patch taking consideration of the x-y flat rotation of forces (ie: Vy will be the velocity in local y-dir of ring, etc.)



		// TODO: JACOBIAN,



		dvx = fwd.Dot( -(pRing->GetWCurr()).Cross(distM) - VparWheel ); // relative speed between center of wheel and patch in the forward direction; using fwd for projection because we are concerned with the angular velocity of the wheel and slip in the wheel's longitudinal direction...
//  * (1 + fabs(VpcRing(1)*0.01) )   ///   (1 + fabs(Vparp(1)*0.05)  ///  fmax(0,(1-fabs(VpcRing(1))/10))  ///   * (1 + 3.15 * (1 - pContrib) ) - VparWheel * pContrib
	dvax = fwd.Dot(va); // Calculated on wheel but force applied to contact patch. // axle velocity in the forward direction


	// Fn has been corrected to ensure it is either positive or null because without contact there are no forces! (Although once the truck suspension is modeled contact should always be present, or almost. This is a weakness of having the road profile force directly on the patch and not having the possibility of the wheel flying a little!)
	Fn = Fint(3); // not further influenced by slope because patch is never actually rotated
	if (Fn < 0)
	{
		Fn = 0;	// TODO: change derivatives in jacobian of Fn accordingly! (and what to do about discontinuity at zero?: should not be much of an issue because the MF forces will be almost null close to zero...)
		boolFn = 0; // using a different var name to clarify code in Jacobian
	}
	else
	{
		boolFn=1;
	}

	// lateral speed of wheel center
	dvay = lat.Dot(va);
	// wheel center drift angle
	// NOTE: getting the angle relative to the fwd axis to properly treat angles > pi/2
	dAlpha = atan2(dvay, fabs(dvax)); // assumes that camber is small (ie velocities taken at wheel center)



// slip ratio, attention, the definition can vary (ie: here divided by wheel center velocity but often uses patch velocity instead) slip ratio is dependent on contact patch velocity (ie: transient slips)
			if (std::abs(dvax) < TRH)
			{
				dSr = (dvax - dvx)/fabs(dvax+copysign(TRHA, dvax)); // copysign is used to either add or subtract the threshold added value depending on the sign of dvax in order to preserve its sign and thus have a more precise value of the ratio at that point
//				dSr = (dvax - dvxr + Vpar(1))/std::abs(dvax+copysign(TRHA, dvax));
			}
			else
			{
				dSr = (dvax -dvx)/fabs(dvax); // here and above, a more "SWIFT" definition of the slip should include the term  + Vpar(1) in the numerator BUT, it seems somewhat wrong because the slip measured by MF experimental tests reads wheel rotational velocity and the drum or moving belt velocity, which is actually vel between ground and wheel... The idea to have patch velocity is interesting but somewhat unclear... (Should it be patch vel w.r.t. wheel?
			}
			if (std::abs(dSr) >= TRHC)
					{
						dSr = copysign(TRHC, dSr);
					}
//			dSr = 1. / pow((fabs(VpcRing.Dot(fwd))+1),3) * dSr;
//			std::cout << VpcRing << std::endl;

					// longitudinal friction coefficient
					dMuX = pMuX0->dGet(dSr);
					F = -fwd*dMuX*Fn; // absolute ref frame, initializing vec. F, force is negative for a positive slip ratio dSr
						dMuY = pMuY0->dGet(dAlpha);
						// force correction
						F -= lat*dMuY*Fn;  // absolute ref frame
		Fpatch = F-Fint; // oct10.2011_pm, changed back after thorough check
		// TODO: jacobian


	// in x dir,
	WorkVec.PutCoef(1, (Fpatch(1)-Mpa*dVPx)/Kpatv(1));  // Fpatch(x,v) applied in absolute frame because patch IS in absolute frame
	WorkVec.PutCoef(2, dVx - dXPx); // V = Xprime
	// in y dir,
	WorkVec.PutCoef(3, (Fpatch(2)-Mpa*dVPy)/Kpatv(2));  // Fpatch(x,v)
	WorkVec.PutCoef(4, dVy - dXPy); // V = Xprime

 // using patch pos to simulate the eff. radius
	dAlpha_t = tan(dAlpha) + S_ht;
	dAlpha_r = tan(dAlpha) + S_hf;
	tr = pTr->dGet(dAlpha_t)*Fn; // pneumatic trail
	M_zr = pMzr->dGet(dAlpha_r)*Fn; // residual torque


	Fr = fwdRing*Fn*(q_sy1 + q_sy3*abs(dvax/dvao))*copysign(1,fwd.Dot(pWheel->GetWCurr().Cross(k))/dvax); // rolling resistance force, force change with velocity and standing waves portion are ignored, as in Pac2006, see swift params=0, note that Fr is considered to work in the direction of the ring because this is where most of the deformations occur, although there is pneumatic damping that will be in a direction between wheel and ring directions...
	//TODO: adjust JACOBIAN! for Fr

	Mz = n*( -tr*Fint_old(2) + M_zr ); // TODO update jacobian because (if was then updated) was using Fint_ring(2) (now using Fint_old because it is local)
	if (bLoadedRadius)
	{
		M = distM.Cross(Fint_ring+Fr) + Mz; // self-aligning torque is applied in the z-direction of the patch on the ring reference frame, thus in the n dir.
	}
	else
	{
		CalculateR_e();
//		std::cout << "Wheel: " << GetLabel() << "Re: " << R_e << " -latRing*(R_e - ddistM)*Fint_ring.Dot(fwdRing): " << -latRing*(R_e - ddistM)*Fint_ring.Dot(fwdRing) << std::endl;
		M = distM.Cross(Fint_ring+Fr) - latRing*(R_e - ddistM)*Fint_ring.Dot(fwdRing) + Mz; //	(with brake lever arm adjusted for Fx [see pac2006 p. 469, eq. 9.236])
	}
	// moment assuming that the patch is always at the proper angular position, in 3D, because it has no inertia.

		WorkVec.Add(5, Fint_ring); // on ring
		WorkVec.Add(8, M); // on ring




	return WorkVec;
}

void
Wheel4::AfterConvergence(
		const VectorHandler& X,
		const VectorHandler& XP)
{
	   bFirstAC = true;

		// calculations done before output,

		// loaded radius (between wheel center and contact center),
			EffRad = Xring-pWheel->GetXCurr()+distM; // vector dist btwn wheel center and patch
			dEffRad = sqrt(EffRad.Dot()); // loaded radius as defined in Pac2006 p. 464 TODO: adjust Jacobian accordingly
		//	Vec3 nEffRad = -EffRad/dEffRad;
		// ,loaded radius (between wheel center and contact center)

			KE = pWheelB->dGetPrivData(1) + pRingB->dGetPrivData(1) + (pow(Vpa(1),2)+pow(Vpa(2),2))*Mpa/2;
			PE = pWheelB->dGetPrivData(2) + pRingB->dGetPrivData(2);
			E = KE + PE;
			CalculateR_e();

			dSa = 180./M_PI*dAlpha;

		// ,calculations done before output
}

void
Wheel4::CalculateR_e()
{
	doublereal R_e_divider = fwd.Dot(pWheel->GetWCurr().Cross(-lat.Cross(fwd)));
	if (fabs(R_e_divider) < TdReDiv)
	{
		R_e = dvx/copysign(TdReDiv,R_e_divider); // calculates the virtual R_e (effective radius) by equating the slip ratio formula to a formula like dSr = (dvax - WheelW*R_e) / dvax
	}
	else
	{
		R_e = dvx/R_e_divider; // calculates the virtual R_e (effective radius) by equating the slip ratio formula to a formula like dSr = (dvax - WheelW*R_e) / dvax
	}
	if ( fabs(R_e) > TdR_e * dR_0 ) R_e = copysign(TdR_e,R_e);
}

doublereal
Wheel4::CapLoop(doublereal Xuncapped) const
{
	while (Xuncapped > RDL) Xuncapped-=RDL;
	return Xuncapped;
}

void 
Wheel4::SetInitialValue(VectorHandler& XCurr) 
{
integer iFirstIndex = iGetFirstIndex();
XCurr.PutCoef(iFirstIndex + 4, XpaPrev(2));
XCurr.PutCoef(iFirstIndex + 2, XpaPrev(1));
XCurr.PutCoef(iFirstIndex + 3, 0);
XCurr.PutCoef(iFirstIndex + 1, 0);
}


void
Wheel4::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	bFirstAP = true;
}


unsigned int
Wheel4::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
Wheel4::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

//	only using one option for now... to set timestep
	return 1;
}


doublereal Wheel4::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= iGetNumPrivData());
    	return dtMax; // this should return the maximum timestep that this wheel is able to take (to be fed into the strategy:change cirective in the MBDyn input file)
}


int
Wheel4::iGetNumConnectedNodes(void) const
{
	// wheel  + ring
	return 2;
}

void
Wheel4::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	// wheel  + ring
	connectedNodes.resize(2);

	connectedNodes[0] = pWheel;
	connectedNodes[1] = pRing;
}

void 
Wheel4::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
   	NO_OP;
}

std::ostream& 
Wheel4::Restart(std::ostream& out) const
{
   	return out << "# not implemented yet" << std::endl;
}

unsigned int
Wheel4::iGetNumDof(void) const
{
	return 4;
}
DofOrder::Order
Wheel4::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

DofOrder::Order
Wheel4::GetEqType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

unsigned
Wheel4::iGetInitialNumDof(void) const
{
	return 0;
}

void 
Wheel4::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Wheel4::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler& 
Wheel4::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{

	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.Resize(0);
	return WorkVec;
}


/* Wheel4 - end */

/* TimeStep - begin */

// FIXME: does not need to be an element, imho

TimeStep::TimeStep(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	timestep						\n"
"Author: 	Louis Gagnon < louis.gagnon.10@ulaval.ca		\n"
"Organization:	Universite Laval			\n"
"							\n"
"						\n"
"									\n"
"	All rights reserved						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// FIXME: why?
	for (DataManager::ElemContainerType::const_iterator e = pDM->begin(Elem::LOADABLE);
			e != pDM->end(Elem::LOADABLE); ++e)
	{
		pWheelE = e->second;
		pWheelsE.push_back(pWheelE);
	}
}

TimeStep::~TimeStep(void)
{
	// destroy private data
	NO_OP;
}

void
TimeStep::Output(OutputHandler& OH) const
{
	NO_OP;
}

void
TimeStep::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
TimeStep::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
TimeStep::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);
	return WorkVec;
}

unsigned int
TimeStep::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
TimeStep::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

//	only using one option for now... to set timestep
	return 1;
}

doublereal TimeStep::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= iGetNumPrivData());

	doublereal dtMax = pWheelsE[0]->dGetPrivData(1);
	for (unsigned iCnt = 1; iCnt < pWheelsE.size(); iCnt++) {
		dtMax = fmin(dtMax, pWheelsE[iCnt]->dGetPrivData(1));
	}
    	return dtMax; // this should return the maximum timestep that this wheel is able to take (to be fed into the strategy:change cirective in the MBDyn input file)
}

int
TimeStep::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
TimeStep::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
TimeStep::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
TimeStep::Restart(std::ostream& out) const
{
	return out << "# TimeStep: not implemented" << std::endl;
}

unsigned int
TimeStep::iGetInitialNumDof(void) const
{
	return 0;
}

void
TimeStep::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
TimeStep::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
TimeStep::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

/* TimeStep - end */


// #ifdef STATIC_MODULES, the function is registered by InitUDE()
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<Wheel4>;
	if (!SetUDE("rigid" "ring" "tire", rf)) {
		delete rf;
		silent_cerr("RigidRingTire: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	// DataManager *pDM = (DataManager *)pdm;
	MBDynParser *pHP = (MBDynParser *)php;

	while (pHP->IsArg()) {
		const char *s = pHP->GetString();
		if (s == 0) {
			silent_cerr("RigidRingTire: "
				"unable to get arg; "
				"module_init(" << module_name << ") "
				"failed" << std::endl);
			return -1;
		}

		rf = new UDERead<Wheel4>;
		if (!SetUDE(s, rf)) {
			delete rf;
			silent_cerr("Wheel4: "
				"unable to set parser for name \"" << s << "\"; "
				"module_init(" << module_name << ") "
				"failed" << std::endl);
			return -1;
		}
	}

	UserDefinedElemRead *rf2 = new UDERead<TimeStep>;
	if (!SetUDE("timestep", rf2)) {
		delete rf2;
		silent_cerr("module-timestep: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

