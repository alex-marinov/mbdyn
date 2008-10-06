/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2008
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

#include <dataman.h>
#include "extedge.h"
#include "modalext.h"
#include "modaledge.h"

#include <fstream>

ExtModalForceBase::~ExtModalForceBase(void)
{
	NO_OP;
}

/* ModalExt - begin */

/* Costruttore */
ModalExt::ModalExt(unsigned int uL,
	Modal *pmodal,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	ExtModalForceBase* pEMF,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pEFH, iCoupling, fOut), 
pModal(pmodal),
bOutputAccelerations(bOutputAccelerations),
pEMF(pEMF),
uFlags(ExtModalForceBase::EMF_NONE),
F(0.),
M(0.)
{
	ASSERT(pModal != 0);
	ASSERT(pEMF != 0);

	f.resize(pModal->uGetNModes());

	// Temporary?
	q.resize(pModal->uGetNModes());
	qP.resize(pModal->uGetNModes());

	if (pModal->uGetNModes() > 0) {
		uFlags |= ExtModalForceBase::EMF_MODAL;
	}

	if (pModal->pGetModalNode() != 0) {
		uFlags |= ExtModalForceBase::EMF_RIGID;
	}
}

ModalExt::~ModalExt(void)
{
	NO_OP;
}

/*
 * Send output to companion software
 */
void
ModalExt::Send(std::ostream& outf, bool bAfterConvergence)
{
	const StructNode *pNode = pModal->pGetModalNode();

	Vec3 x;
	Mat3x3 R;
	Vec3 v;
	Vec3 w;

	if (pNode) {
		x = pNode->GetXCurr();
		R = pNode->GetRCurr();
		v = pNode->GetVCurr();
		w = pNode->GetWCurr();

	} else {
		x = Zero3;
		R = Eye3;
		v = Zero3;
		w = Zero3;
	}

	// Temporary?
	const VecN& a = pModal->GetA();
	const VecN& b = pModal->GetB();
	for (unsigned i = 0; i < q.size(); i++) {
		q[i] = a(i + 1);
		qP[i] = b(i + 1);
	}

	pEMF->Send(outf, uFlags, GetLabel(), x, R, v, w, q, qP);

#if 0
	if (pNode) {
		Mat3x3 RT = pNode->GetRCurr().Transpose();
		outf << pNode->GetLabel()
			<< " " << pNode->GetXCurr()
			<< " " << pNode->GetRCurr()
			<< " " << RT*pNode->GetVCurr()
			<< " " << RT*pNode->GetWCurr();
		if (bOutputAccelerations) {
			outf
				<< " " << pNode->GetXPPCurr()
				<< " " << pNode->GetWPCurr();
		}
		outf << std::endl;

	} else {
		outf << 0
			<< " " << Zero3
			<< " " << Zero3x3
			<< " " << Zero3
			<< " " << Zero3;
		if (bOutputAccelerations) {
			outf
				<< " " << Zero3
				<< " " << Zero3;
		}
		outf << std::endl;
	}

	const VecN& a = pModal->GetA();
	const VecN& b = pModal->GetB();
	const VecN& bP = pModal->GetBP();
	for (integer i = 1; i < pModal->uGetNModes(); i++) {
		outf << a(i) << " " << b(i);
		if (bOutputAccelerations) {
			outf << " " << bP(i);
		}
		outf << std::endl;
	}
#endif
}

void
ModalExt::Recv(std::istream& inf)
{
	unsigned uLabel = 0;
	unsigned uOutFlags = pEMF->Recv(inf, uFlags, uLabel, F, M, f);

	if (uOutFlags & ExtModalForceBase::EMF_ERR) {
		silent_cerr("ModalExt(" << GetLabel() << "): "
			"error while reading modal force data" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (uFlags != uOutFlags) {
		silent_cerr("ModalExt(" << GetLabel() << "): "
			"error while reading modal force data "
			"(" << uOutFlags << "!=" << uFlags << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

SubVectorHandler&
ModalExt::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	integer iR, iC;
	WorkSpaceDim(&iR, &iC);
	WorkVec.ResizeReset(iR);

	integer iIdx = 1;
	const StructNode *pNode = pModal->pGetModalNode();
	if (pNode) {
		integer iFirstIndex = pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Put(1, pNode->GetRCurr()*F);
		WorkVec.Put(4, pNode->GetRCurr()*M);

		iIdx += 6;

	}

	integer iModalIndex = pModal->iGetModalIndex() + pModal->uGetNModes() + 1;
	for (integer iMode = 0; iMode < pModal->uGetNModes(); iMode++) {
		WorkVec.PutItem(iIdx + iMode, iModalIndex + iMode, f[iMode]);
	}

	return WorkVec;
}

void
ModalExt::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Forces();

		out << GetLabel() << " " << F << " " << M;
		for (std::vector<doublereal>::const_iterator i = f.begin(); i != f.end(); i++) {
			out << " " << *i;
		}
		out << std::endl;
	}
}
   
Elem*
ReadModalExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH;
	int iCoupling;
	
	ReadExtForce(pDM, HP, uLabel, pEFH, iCoupling);

	Modal *pModal = dynamic_cast<Modal *>(pDM->ReadElem(HP, Elem::JOINT));
	if (pModal == 0) {
		silent_cerr("ModalExt(" << uLabel << "): illegal Modal joint "
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	bool bOutputAccelerations(false);
	if (HP.IsKeyWord("accelerations")) {
		bOutputAccelerations = true;
	}

	ExtModalForceBase *pEMF = 0;
	if (dynamic_cast<ExtFileHandlerEDGE *>(pEFH) != 0) {
		SAFENEWWITHCONSTRUCTOR(pEMF, ExtModalForceEDGE,
			ExtModalForceEDGE(pDM));

	} else {
		silent_cerr("ModalExt(" << uLabel << "): "
			"unknown external force type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalExt,
		ModalExt(uLabel, pModal, bOutputAccelerations,
			pEFH, pEMF, iCoupling, fOut));

	return pEl;
}

