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
	DataManager *pDM,
	Modal *pmodal,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	ExtModalForceBase* pEMF,
	int iCoupling,
	ExtModalForceBase::BitMask bm,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pDM, pEFH, iCoupling, fOut), 
pModal(pmodal),
bOutputAccelerations(bOutputAccelerations),
pEMF(pEMF),
uFlags(ExtModalForceBase::EMF_NONE),
F(0.),
M(0.)
{
	ASSERT(pModal != 0);
	ASSERT(pEMF != 0);

	ASSERT((bm & ExtModalForceBase::EMF_ALL) != 0);
	ASSERT((bm & ~ExtModalForceBase::EMF_ALL) == 0);

	if (bm & ExtModalForceBase::EMF_MODAL) {
		f.resize(pModal->uGetNModes());

		// Temporary?
		q.resize(pModal->uGetNModes());
		qP.resize(pModal->uGetNModes());

		ASSERT(pModal->uGetNModes() > 0);
		uFlags |= ExtModalForceBase::EMF_MODAL;

	} else {
		f.resize(0);
		q.resize(0);
		qP.resize(0);
	}

	if (bm & ExtModalForceBase::EMF_RIGID) {
		ASSERT(pModal->pGetModalNode() != 0);
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
	Vec3 x;
	Mat3x3 R;
	Vec3 v;
	Vec3 w;

	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		const StructNode *pNode = pModal->pGetModalNode();

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
	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		const StructNode *pNode = pModal->pGetModalNode();

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
	}

	if (uFlags & ExtModalForceBase::EMF_MODAL) {
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

	const StructNode *pNode = pModal->pGetModalNode();
	integer iIdx = 1;
	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		integer iFirstIndex = pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Put(1, pNode->GetRCurr()*F);
		WorkVec.Put(4, pNode->GetRCurr()*M);
	}

	if (pNode) {
		iIdx += 6;
	}

	if (uFlags & ExtModalForceBase::EMF_MODAL) {
		integer iModalIndex = pModal->iGetModalIndex() + pModal->uGetNModes() + 1;
		for (integer iMode = 0; iMode < pModal->uGetNModes(); iMode++) {
			WorkVec.PutItem(iIdx + iMode, iModalIndex + iMode, f[iMode]);
		}
	}

	return WorkVec;
}

void
ModalExt::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Forces();

		out << GetLabel();

		if (uFlags & ExtModalForceBase::EMF_RIGID) {
			out << " " << F << " " << M;
		}

		if (uFlags & ExtModalForceBase::EMF_MODAL) {
			for (std::vector<doublereal>::const_iterator i = f.begin(); i != f.end(); i++) {
				out << " " << *i;
			}
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

	// safe default: both rigid and modal
	ExtModalForceBase::BitMask bm = ExtModalForceBase::EMF_ALL;
	if (HP.IsKeyWord("type")) {
		if (HP.IsKeyWord("rigid")) {
			bm = ExtModalForceBase::EMF_RIGID;
		} else if (HP.IsKeyWord("modal")) {
			bm = ExtModalForceBase::EMF_MODAL;
		} else if (HP.IsKeyWord("all")) {
			bm = ExtModalForceBase::EMF_ALL;
		} else {
			silent_cerr("ModalExt(" << uLabel << "): unknown ModalExt type "
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	ExtModalForceBase *pEMF = 0;
	if (dynamic_cast<ExtFileHandlerEDGE *>(pEFH) != 0) {
		// EDGE needs two separate ModalExt elements,
		// one for the rigid part and one for the modal part

		switch (bm & ExtModalForceBase::EMF_ALL) {
		case ExtModalForceBase::EMF_RIGID:
			SAFENEWWITHCONSTRUCTOR(pEMF, ExtRigidForceEDGE,
				ExtRigidForceEDGE(pDM));
			break;

		case ExtModalForceBase::EMF_MODAL:
			SAFENEWWITHCONSTRUCTOR(pEMF, ExtModalForceEDGE,
				ExtModalForceEDGE(pDM));
			break;

		case ExtModalForceBase::EMF_ALL:
			silent_cerr("ModalExt(" << uLabel << "): "
				"EDGE ExtFileHandler can only be used "
				"when ModalExt is either \"rigid\" "
				"or \"modal\" but not when it is \"all\" "
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
			ASSERT(0);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	// add more types

	} else {
		silent_cerr("ModalExt(" << uLabel << "): "
			"unknown external force type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);

	if (HP.IsArg()) {
		silent_cerr("ModalExt(" << uLabel << "): "
			"semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalExt,
		ModalExt(uLabel, pDM, pModal, bOutputAccelerations,
			pEFH, pEMF, iCoupling, bm, fOut));

	return pEl;
}

