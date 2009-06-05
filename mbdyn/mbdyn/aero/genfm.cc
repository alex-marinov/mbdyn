/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "aeroelem.h"
#include "genfm.h"
#include "bisec.h"

/* GenericAerodynamicForce - begin */

/* Assemblaggio residuo */
void
GenericAerodynamicForce::AssVec(SubVectorHandler& WorkVec)
{
	/* velocity at aerodynamic center */

	/* position of aerodynamic center */
	const Mat3x3& Rn(pNode->GetRCurr());
	Vec3 f(Rn*tilde_f);
	Mat3x3 R(Rn*tilde_Ra);
	Vec3 Xca(pNode->GetXCurr() + f);
	Vec3 Vca(pNode->GetVCurr() + f.Cross(pNode->GetWCurr()));

	Vec3 VTmp(0.);
      	if (fGetAirVelocity(VTmp, Xca)) {
		Vca -= VTmp;
	}

#if 0
	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * aggiunge alla velocita' la velocita' indotta
	 */
	if (pRotor != NULL) {
 		Vca += pRotor->GetInducedVelocity(Xca);
	}
#endif

	Vec3 V(R.Transpose()*Vca);

	/* FIXME: according to source, the rotation sequence
	 * is alpha, then beta */
	dAlpha = atan2(V(3), V(1));
	dBeta = -atan2(V(2), V(1));

	doublereal rho, c, p, T;
	GetAirProps(Xca, rho, c, p, T);	/* p, T no used yet */
	doublereal q = Vca.Dot()*rho/2.;

	doublereal dScaleForce = q*dRefSurface;
	doublereal dScaleMoment = dScaleForce*dRefSurface;

	int nAlpha = pData->Alpha.size() - 1;
	int nBeta = pData->Beta.size() - 1;

	ASSERT(nAlpha >= 0);
	ASSERT(nBeta >= 0);

	/* TODO: search for Alpha, Beta */
	if (dAlpha <= pData->Alpha[0]) {
		/* smooth out coefficients if Alpha does not span -180 => 180 */
		doublereal dAlphaX = (dAlpha - pData->Alpha[0])/(-M_PI - pData->Alpha[0]);
		doublereal dSmoothAlpha = (std::cos(M_PI*dAlphaX) + 1)/2.;

		if (dBeta <= pData->Beta[0]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[0])/(-M_PI - pData->Beta[0]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			F = Vec3(&pData->Data[0][0].dCoef[0])*(dScaleForce*dSmoothAlpha*dSmoothBeta);
			M = Vec3(&pData->Data[0][0].dCoef[3])*(dScaleMoment*dSmoothAlpha*dSmoothBeta);

		} else if (dBeta >= pData->Beta[nBeta]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[nBeta])/(M_PI - pData->Beta[nBeta]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			F = Vec3(&pData->Data[nBeta][0].dCoef[0])*(dScaleForce*dSmoothAlpha*dSmoothBeta);
			M = Vec3(&pData->Data[nBeta][0].dCoef[3])*(dScaleMoment*dSmoothAlpha*dSmoothBeta);

		} else {
			int iBeta = bisec<doublereal>(&pData->Beta[0], dBeta, 0, nBeta);

			ASSERT(iBeta >= 0);
			ASSERT(iBeta < nBeta);

			doublereal ddBeta = pData->Beta[iBeta + 1] - pData->Beta[iBeta];
			doublereal d1Beta = (pData->Beta[iBeta + 1] - dBeta)/ddBeta;
			doublereal d2Beta = (dBeta - pData->Beta[iBeta])/ddBeta;

			GenericAerodynamicData::GenericAerodynamicCoef c
				= pData->Data[iBeta][0]*d1Beta + pData->Data[iBeta + 1][0]*d2Beta;

			F = Vec3(&c.dCoef[0])*(dScaleForce*dSmoothAlpha);
			M = Vec3(&c.dCoef[3])*(dScaleMoment*dSmoothAlpha);
		}

	} else if (dAlpha >= pData->Alpha[nAlpha]) {
		/* smooth out coefficients if Alpha does not span -180 => 180 */
		doublereal dAlphaX = (dAlpha - pData->Alpha[nAlpha])/(-M_PI - pData->Alpha[nAlpha]);
		doublereal dSmoothAlpha = (std::cos(M_PI*dAlphaX) + 1)/2.;

		if (dBeta <= pData->Beta[0]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[0])/(-M_PI - pData->Beta[0]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			F = Vec3(&pData->Data[0][nAlpha].dCoef[0])*(dScaleForce*dSmoothAlpha*dSmoothBeta);
			M = Vec3(&pData->Data[0][nAlpha].dCoef[3])*(dScaleMoment*dSmoothAlpha*dSmoothBeta);

		} else if (dBeta >= pData->Beta[nBeta]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[nBeta])/(M_PI - pData->Beta[nBeta]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			F = Vec3(&pData->Data[nBeta][nAlpha].dCoef[0])*(dScaleForce*dSmoothAlpha*dSmoothBeta);
			M = Vec3(&pData->Data[nBeta][nAlpha].dCoef[3])*(dScaleMoment*dSmoothAlpha*dSmoothBeta);

		} else {
			int iBeta = bisec<doublereal>(&pData->Beta[0], dBeta, 0, nBeta);

			ASSERT(iBeta >= 0);
			ASSERT(iBeta < nBeta);

			doublereal ddBeta = pData->Beta[iBeta + 1] - pData->Beta[iBeta];
			doublereal d1Beta = (pData->Beta[iBeta + 1] - dBeta)/ddBeta;
			doublereal d2Beta = (dBeta - pData->Beta[iBeta])/ddBeta;

			GenericAerodynamicData::GenericAerodynamicCoef c
				= pData->Data[iBeta][nAlpha]*d1Beta + pData->Data[iBeta + 1][nAlpha]*d2Beta;

			F = Vec3(&c.dCoef[0])*(dScaleForce*dSmoothAlpha);
			M = Vec3(&c.dCoef[3])*(dScaleMoment*dSmoothAlpha);
		}

	} else {
		int iAlpha = bisec<doublereal>(&pData->Alpha[0], dAlpha, 0, nAlpha);

		ASSERT(iAlpha >= 0);
		ASSERT(iAlpha < nAlpha);

		doublereal ddAlpha = pData->Alpha[iAlpha + 1] - pData->Alpha[iAlpha];
		doublereal d1Alpha = (pData->Alpha[iAlpha + 1] - dAlpha)/ddAlpha;
		doublereal d2Alpha = (dAlpha - pData->Alpha[iAlpha])/ddAlpha;

		if (dBeta <= pData->Beta[0]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[0])/(-M_PI - pData->Beta[0]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			GenericAerodynamicData::GenericAerodynamicCoef c
				= pData->Data[0][iAlpha]*d1Alpha + pData->Data[0][iAlpha + 1]*d2Alpha;

			F = Vec3(&c.dCoef[0])*(dScaleForce*dSmoothBeta);
			M = Vec3(&c.dCoef[3])*(dScaleMoment*dSmoothBeta);

		} else if (dBeta >= pData->Beta[nBeta]) {
			/* smooth out coefficients if Beta does not span -180 => 180 */
			doublereal dBetaX = (dBeta - pData->Beta[nBeta])/(M_PI - pData->Beta[nBeta]);
			doublereal dSmoothBeta = (std::cos(M_PI*dBetaX) + 1)/2.;

			GenericAerodynamicData::GenericAerodynamicCoef c
				= pData->Data[nBeta][iAlpha]*d1Alpha + pData->Data[nBeta][iAlpha + 1]*d2Alpha;

			F = Vec3(&c.dCoef[0])*(dScaleForce*dSmoothBeta);
			M = Vec3(&c.dCoef[3])*(dScaleMoment*dSmoothBeta);

		} else {
			int iBeta = bisec<doublereal>(&pData->Beta[0], dBeta, 0, nBeta);

			ASSERT(iBeta >= 0);
			ASSERT(iBeta < nBeta);

			doublereal ddBeta = pData->Beta[iBeta + 1] - pData->Beta[iBeta];
			doublereal d1Beta = (pData->Beta[iBeta + 1] - dBeta)/ddBeta;
			doublereal d2Beta = (dBeta - pData->Beta[iBeta])/ddBeta;

			GenericAerodynamicData::GenericAerodynamicCoef c1
				= pData->Data[iBeta][iAlpha]*d1Alpha + pData->Data[iBeta][iAlpha + 1]*d2Alpha;
			GenericAerodynamicData::GenericAerodynamicCoef c2
				= pData->Data[iBeta + 1][iAlpha]*d1Alpha + pData->Data[iBeta + 1][iAlpha + 1]*d2Alpha;
			
			GenericAerodynamicData::GenericAerodynamicCoef c = c1*d1Beta + c2*d2Beta;

			F = Vec3(&c.dCoef[0])*dScaleForce;
			M = Vec3(&c.dCoef[3])*dScaleMoment;
		}
	}

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);
}

GenericAerodynamicForce::GenericAerodynamicForce(unsigned int uLabel,
	const DofOwner *pDO,
	const StructNode* pN,
	const Vec3& fTmp, const Mat3x3& RaTmp,
	doublereal dS, doublereal dL,
	GenericAerodynamicData *pD,
	flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut),
pNode(pN),
dRefSurface(dS),
dRefLength(dL),
tilde_f(fTmp),
tilde_Ra(RaTmp),
F(0.),
M(0.),
pData(pD)
{
	NO_OP;
}

GenericAerodynamicForce::~GenericAerodynamicForce(void)
{
	if (pData != 0) {
		delete pData;
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
GenericAerodynamicForce::Restart(std::ostream& out) const
{
	return out;
}

/* assemblaggio residuo */
SubVectorHandler&
GenericAerodynamicForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(6);

	integer iNodeFirstIndex = pNode->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
GenericAerodynamicForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
      		OH.Aerodynamic()
			<< std::setw(8) << GetLabel()
			<< " " << dAlpha*180./M_PI
			<< " " << dBeta*180./M_PI
			<< " " << F << " " << M << std::endl;
	}
}

/* assemblaggio residuo */
SubVectorHandler&
GenericAerodynamicForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	return WorkVec;
}

extern Elem *
ReadGenericAerodynamicForce(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel);

/* GenericAerodynamicForce - end */

/* GenericAerodynamicData - begin */

/* GenericAerodynamicData::GenericAerodynamicCoef - begin */

GenericAerodynamicData::GenericAerodynamicCoef::GenericAerodynamicCoef(void)
{
	NO_OP;
}

GenericAerodynamicData::GenericAerodynamicCoef::GenericAerodynamicCoef(
	const GenericAerodynamicData::GenericAerodynamicCoef& c)
{
	dCoef[0] = c.dCoef[0];
	dCoef[1] = c.dCoef[1];
	dCoef[2] = c.dCoef[2];
	dCoef[3] = c.dCoef[3];
	dCoef[4] = c.dCoef[4];
	dCoef[5] = c.dCoef[5];
}

GenericAerodynamicData::GenericAerodynamicCoef
GenericAerodynamicData::GenericAerodynamicCoef::operator + (
	const GenericAerodynamicData::GenericAerodynamicCoef& c) const
{
	GenericAerodynamicData::GenericAerodynamicCoef retval;

	retval.dCoef[0] = dCoef[0] + c.dCoef[0];
	retval.dCoef[1] = dCoef[1] + c.dCoef[1];
	retval.dCoef[2] = dCoef[2] + c.dCoef[2];
	retval.dCoef[3] = dCoef[3] + c.dCoef[3];
	retval.dCoef[4] = dCoef[4] + c.dCoef[4];
	retval.dCoef[5] = dCoef[5] + c.dCoef[5];

	return retval;
}

GenericAerodynamicData::GenericAerodynamicCoef
GenericAerodynamicData::GenericAerodynamicCoef::operator - (
	const GenericAerodynamicData::GenericAerodynamicCoef& c) const
{
	GenericAerodynamicData::GenericAerodynamicCoef retval;

	retval.dCoef[0] = dCoef[0] - c.dCoef[0];
	retval.dCoef[1] = dCoef[1] - c.dCoef[1];
	retval.dCoef[2] = dCoef[2] - c.dCoef[2];
	retval.dCoef[3] = dCoef[3] - c.dCoef[3];
	retval.dCoef[4] = dCoef[4] - c.dCoef[4];
	retval.dCoef[5] = dCoef[5] - c.dCoef[5];

	return retval;
}

GenericAerodynamicData::GenericAerodynamicCoef
GenericAerodynamicData::GenericAerodynamicCoef::operator * (const doublereal& d) const
{
	GenericAerodynamicData::GenericAerodynamicCoef retval;

	retval.dCoef[0] = dCoef[0]*d;
	retval.dCoef[1] = dCoef[1]*d;
	retval.dCoef[2] = dCoef[2]*d;
	retval.dCoef[3] = dCoef[3]*d;
	retval.dCoef[4] = dCoef[4]*d;
	retval.dCoef[5] = dCoef[5]*d;

	return retval;
}

GenericAerodynamicData::GenericAerodynamicCoef
GenericAerodynamicData::GenericAerodynamicCoef::operator / (const doublereal& d) const
{
	GenericAerodynamicData::GenericAerodynamicCoef retval;

	retval.dCoef[0] = dCoef[0]/d;
	retval.dCoef[1] = dCoef[1]/d;
	retval.dCoef[2] = dCoef[2]/d;
	retval.dCoef[3] = dCoef[3]/d;
	retval.dCoef[4] = dCoef[4]/d;
	retval.dCoef[5] = dCoef[5]/d;

	return retval;
}

/* GenericAerodynamicData::GenericAerodynamicCoef - end */


static GenericAerodynamicData *
ReadGenericAerodynamicData(const std::string& fname)
{
	std::ifstream in(fname.c_str());

	if (!in) {
		silent_cerr("ReadGenericAerodynamicData: "
			"unable to open file \"" << fname << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	char buf[LINE_MAX];
	int nAlpha, nBeta;
	int c;

	/* skip comments */
	for (c = in.get(); c == '%' || c == '#'; c = in.get()) {
		/* discard to end of line */
		in.getline(buf, sizeof(buf));
	}
	in.putback(c);

	/* get the size of the matrices */
	in >> nAlpha >> nBeta;
	/* discard to end of line */
	in.getline(buf, sizeof(buf));

	if (!in) {
		silent_cerr("ReadGenericAerodynamicData(" << fname << "): "
			"unable to read size of data matrix "
			"from file \"" << fname << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (nAlpha <= 0) {
		silent_cerr("ReadGenericAerodynamicData(" << fname << "): "
			"invalid number of angles of attack " << nAlpha << " "
			"from file \"" << fname << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (nBeta <= 0) {
		silent_cerr("ReadGenericAerodynamicData(" << fname << "): "
			"invalid number of sideslip angles " << nBeta << " "
			"from file \"" << fname << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* skip comments */
	for (c = in.get(); c == '%' || c == '#'; c = in.get()) {
		/* discard to end of line */
		in.getline(buf, sizeof(buf));
	}
	in.putback(c);

	if (!in) {
		silent_cerr("ReadGenericAerodynamicData(" << fname << "): "
			"unable to get to data "
			"in file \"" << fname << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	GenericAerodynamicData *pData = new GenericAerodynamicData;
	pData->name = fname;
	pData->nAlpha = nAlpha;
	pData->nBeta = nBeta;

	pData->Alpha.resize(nAlpha);
	pData->Beta.resize(nBeta);
	pData->Data.resize(nBeta);

	/* get the matrices */
	for (int iBeta = 0; iBeta < nBeta; iBeta++) {
		pData->Data[iBeta].resize(nAlpha);

		for (int iAlpha = 0; iAlpha < nAlpha; iAlpha++) {
			doublereal dCoef;

			/* read (and check) alpha */
			in >> dCoef;
			if (iBeta == 0) {
				pData->Alpha[iAlpha] = dCoef;

			} else if (dCoef != pData->Alpha[iAlpha]) {
				silent_cerr("ReadGenericAerodynamicData"
					"(" << fname << "): "
					"inconsistent data, "
					"Alpha[" << iAlpha << "]"
						"=" << dCoef << " "
					"for Beta[" << iBeta << "]"
						"=" << pData->Beta[iBeta] << " "
					"differs from previous, "
						<< pData->Alpha[iAlpha]
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			/* read (and check) beta */
			in >> dCoef;
			if (iAlpha == 0) {
				pData->Beta[iBeta] = dCoef;

			} else if (dCoef != pData->Beta[iBeta]) {
				silent_cerr("ReadGenericAerodynamicData"
					"(" << fname << "): "
					"inconsistent data, "
					"Beta[" << iBeta << "]"
						"=" << dCoef << " "
					"for Alpha[" << iAlpha << "]"
						"=" << pData->Alpha[iAlpha] << " "
					"differs from previous, "
						<< pData->Beta[iBeta]
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			for (int iCoef = 0; iCoef < 6; iCoef++) {
				in >> dCoef;

				pData->Data[iBeta][iAlpha].dCoef[iCoef] = dCoef;
			}

			/* discard to end of line */
			if (iAlpha < nAlpha - 1 && iBeta < nBeta - 1) {
				in.getline(buf, sizeof(buf));

				if (!in) {
					silent_cerr("ReadGenericAerodynamicData"
						"(" << fname << "): "
						"unable to read data past "
						"iAlpha=" << iAlpha << ", "
						"iBeta=" << iBeta << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
		}
	}

	/* deg => radian */
	for (int iAlpha = 0; iAlpha < nAlpha; iAlpha++) {
		pData->Alpha[iAlpha] *= M_PI/180.;

		if (iAlpha == 0) {
			continue;
		}

		if ( pData->Alpha[iAlpha] <= pData->Alpha[iAlpha - 1]) {
			silent_cerr("ReadGenericAerodynamicData"
				"(" << fname << "): "
				"strict ordering violated between "
				"Alpha #" << iAlpha - 1 << " and "
				"Alpha #" << iAlpha << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* deg => radian */
	for (int iBeta = 0; iBeta < nBeta; iBeta++) {
		pData->Beta[iBeta] *= M_PI/180.;

		if (iBeta == 0) {
			continue;
		}

		if ( pData->Beta[iBeta] <= pData->Beta[iBeta - 1]) {
			silent_cerr("ReadGenericAerodynamicData"
				"(" << fname << "): "
				"strict ordering violated between "
				"Beta #" << iBeta - 1 << " and "
				"Beta #" << iBeta << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	return pData;
}

Elem *
ReadGenericAerodynamicForce(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel)
{
   	/* Nodo */
	StructNode* pNode = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));

	/* The offset is in the reference frame of the node */
	ReferenceFrame RF(pNode);
	Vec3 f(HP.GetPosRel(RF));

	/* the orientation is in flight mechanics (FIXME?) reference frame:
	 * X forward, Y to the right, Z down */
	Mat3x3 Ra(HP.GetRotRel(RF));

	/* 1. by default, which means that coefficients are only normalized
	 * by the dynamic pressure */
	doublereal dRefSurface = 1.;
	doublereal dRefLength = 1.;

	if (HP.IsKeyWord("reference" "surface")) {
		dRefSurface = HP.GetReal();
	}

	if (HP.IsKeyWord("reference" "length")) {
		dRefLength = HP.GetReal();
	}

	/* TODO: allow to reference previously loaded data */
	std::string fname(HP.GetFileName());
	GenericAerodynamicData *pData = ReadGenericAerodynamicData(fname);

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, GenericAerodynamicForce,
		GenericAerodynamicForce(uLabel, pDO,
			pNode, f, Ra,
			dRefSurface, dRefLength,
			pData, fOut));

	return pEl;
}

/* GenericAerodynamicData - end */
