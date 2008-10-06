/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

/* Aircraft Instruments */

#include <dataman.h>
#include <instruments.h>

AircraftInstruments::AircraftInstruments(unsigned int uLabel, 
		const StructNode* pN, const Mat3x3 &R, flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, fOut),
pNode(pN),
Rh(R)
{
	NO_OP;
}

AircraftInstruments::~AircraftInstruments(void)
{
	NO_OP;
}

void
AircraftInstruments::Update(void)
{
	Vec3 X(pNode->GetXCurr());
	Mat3x3 R(pNode->GetRCurr()*Rh);
	Vec3 V(pNode->GetVCurr());
	Vec3 VV = V;
	Vec3 e1(R.GetVec(1));
	Vec3 e2(R.GetVec(2));
	Vec3 e3(R.GetVec(3));

	Vec3 VTmp(0.);
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
	dMeasure[BANK] = std::atan2(e2(3), e2(2));

	/* turn */
	dMeasure[TURN] = 0.;	/* FIXME */

	/* slip */
	dMeasure[SLIP] = 0.;	/* FIXME */

	/* vertical speed */
	dMeasure[VERTICALSPEED] = VV(3);

	/* angle of attack */
	VTmp = R.Transpose()*VV;
	dMeasure[AOA] = std::atan2(VTmp(3), VTmp(1));

	/* heading */
	dMeasure[HEADING] = std::atan2(e1(2), e1(1));
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
	
	if (strcasecmp(s, "airspeed") == 0) {
		return AIRSPEED;
	}

	if (strcasecmp(s, "ground" "speed") == 0) {
		return GROUNDSPEED;
	}

	if (strcasecmp(s, "altitude") == 0) {
		return ALTITUDE;
	}

	if (strcasecmp(s, "attitude") == 0) {
		return ATTITUDE;
	}

	if (strcasecmp(s, "bank") == 0) {
		return BANK;
	}

	if (strcasecmp(s, "turn") == 0) {
		return TURN;
	}

	if (strcasecmp(s, "slip") == 0) {
		return SLIP;
	}

	if (strcasecmp(s, "vertical" "speed") == 0) {
		return VERTICALSPEED;
	}

	if (strcasecmp(s, "aoa") == 0
			|| strcasecmp(s, "angle" "of" "attack") == 0)
	{
		return AOA;
	}

	if (strcasecmp(s, "heading") == 0) {
		return HEADING;
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

Elem *
ReadAircraftInstruments(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	Elem *pEl = NULL;

	StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
	Mat3x3 R = Eye3;
	if (HP.IsKeyWord("orientation")) {
		R = HP.GetRotRel(ReferenceFrame(pNode));
	}
	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);

	SAFENEWWITHCONSTRUCTOR(pEl, AircraftInstruments,
			AircraftInstruments(uLabel, pNode, R, fOut));

	return pEl;
}

