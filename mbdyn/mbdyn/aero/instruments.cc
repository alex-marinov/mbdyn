/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Aircraft Instruments */

#include "dataman.h"
#include "instruments.h"

AircraftInstruments::AircraftInstruments(unsigned int uLabel, 
	const DofOwner *pDO,
	const StructNode* pN, const Mat3x3 &R, flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
pNode(pN),
Rh(R)
{
	memset(&dMeasure[0], 0, sizeof(dMeasure));
}

AircraftInstruments::~AircraftInstruments(void)
{
	NO_OP;
}

Elem::Type
AircraftInstruments::GetElemType(void) const
{
	return Elem::AERODYNAMIC;
}

void
AircraftInstruments::Update(void)
{
	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr()*Rh);
	const Vec3& V(pNode->GetVCurr());
	const Vec3& Omega(pNode->GetWCurr());
	Vec3 VV = V;
	const Vec3 e1(R.GetVec(1));
	const Vec3 e2(R.GetVec(2));
	const Vec3 e3(R.GetVec(3));

	/*
	 * Assumptions:
	 * - world is "flat"
	 * - world "z" is positive upward
	 * (flight mechanics conventions?)
	 * - aircraft "x" is positive forward
	 * - aircraft "z" is positive downward
	 * - aircraft "y" is positive towards the right of the pilot
	 * - pitch is positive nose up
	 * - roll is positive right wing down
	 * - yaw is positive right wing backward
	 */

	Vec3 VTmp(Zero3);
      	if (fGetAirVelocity(VTmp, X)) {
		VV -= VTmp;
	}

	/* airspeed */
	dMeasure[AIRSPEED] = VV.Norm();

	/* groundspeed */
	VTmp = V;
	VTmp(3) = 0.;
	dMeasure[GROUNDSPEED] = VTmp.Norm();

	/* altitude */
	dMeasure[ALTITUDE] = X(3);

	/* attitude */
	/* FIXME: better asin(e1(3)) ? */
	dMeasure[ATTITUDE] = std::atan2(e1(3), e1(1));

	/* bank */
	/* FIXME: better asin(e2(3)) ? */
	dMeasure[BANK] = -std::atan2(e2(3), e2(2));

	/* turn */
	dMeasure[TURN] = 0.;	/* FIXME */

	/* slip */
	dMeasure[SLIP] = 0.;	/* FIXME */

	/* vertical speed */
	dMeasure[VERTICALSPEED] = VV(3);

	/* angle of attack */
	VTmp = R.MulTV(VV);
	dMeasure[AOA] = std::atan2(VTmp(3), VTmp(1));

	/* heading */
	/* FIXME: assumes a flat world! N={1,0,0}, W={0,1,0} */
	dMeasure[HEADING] = -std::atan2(e1(2), e1(1));

	/* longitude */
	/* FIXME: ??? */
	dMeasure[LONGITUDE] = -X(2);	// /EARTH_RADIUS?

	/* latitude */
	/* FIXME: ??? */
	dMeasure[LATITUDE] = X(1);	// /EARTH_RADIUS?

	VTmp = R.MulTV(Omega);
	dMeasure[ROLLRATE] = VTmp(1);
	dMeasure[PITCHRATE] = VTmp(2);
	dMeasure[YAWRATE] = VTmp(3);
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AircraftInstruments::Restart(std::ostream& out) const
{
	out << "aircraft instruments: " << GetLabel()
		<< ", " << pNode->GetLabel() 
		<< ", orientation, ", Rh.Write(out)
		<< ";" << std::endl;

	return out;
}
	
void
AircraftInstruments::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Aerodynamic()
			<< std::setw(8) << GetLabel();

		for (int iCnt = 1; iCnt < LASTMEASURE; iCnt++) {
			out << " " << dMeasure[iCnt];
		}

		out << std::endl;
	}
}

/* Dimensioni del workspace */
void
AircraftInstruments::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}
	
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
AircraftInstruments::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AircraftInstruments::AssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}
	
/* assemblaggio residuo */
SubVectorHandler&
AircraftInstruments::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	WorkVec.Resize(0);

	Update();

	return WorkVec;
}
	
/* Dati privati */
unsigned int
AircraftInstruments::iGetNumPrivData(void) const
{
	return LASTMEASURE-1;
}

unsigned int
AircraftInstruments::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	struct {
		const char *s;
		int i;
	} s2i[] = {
		{ "airspeed", AIRSPEED },
		{ "ground" "speed", GROUNDSPEED },
		{ "altitude", ALTITUDE },
		{ "attitude", ATTITUDE },
		{ "bank", BANK },
		{ "turn", TURN },
		{ "slip", SLIP },
		{ "vertical" "speed", VERTICALSPEED},
		{ "angle" "of" "attack", AOA },
		{ "aoa", AOA },
		{ "heading", HEADING },
		{ "longitude", LONGITUDE },
		{ "latitude", LATITUDE },
		{ "rollrate", ROLLRATE },
		{ "pitchrate", PITCHRATE },
		{ "yawrate", YAWRATE },
		{ 0 }
	};

	for (int i = 0; s2i[i].s != 0; i++) {
		if (strcasecmp(s, s2i[i].s) == 0) {
			return s2i[i].i;
		}
	}

	return 0;
}

doublereal
AircraftInstruments::dGetPrivData(unsigned int i) const
{
	if (i <= 0 || i >= LASTMEASURE) {
		silent_cerr("AircraftInstruments(" << GetLabel() << "): "
			"illegal measure " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return dMeasure[i];
}

void
AircraftInstruments::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pNode;
}

Elem *
ReadAircraftInstruments(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel)
{
	Elem *pEl = NULL;

	const StructNode* pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	Mat3x3 R = Eye3;
	if (HP.IsKeyWord("orientation")) {
		if (HP.IsKeyWord("flight" "mechanics")) {
			R = ::Eye3;

		} else if (HP.IsKeyWord("aeroelasticity")) {
			R = Mat3x3(-1., 0., 0., 0., 1., 0., 0., 0., -1.);

		} else {
			R = HP.GetRotRel(ReferenceFrame(pNode));
		}
	}
	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);

	SAFENEWWITHCONSTRUCTOR(pEl, AircraftInstruments,
		AircraftInstruments(uLabel, pDO, pNode, R, fOut));

	return pEl;
}

