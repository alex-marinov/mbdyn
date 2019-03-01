/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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

#include <fstream>
#include <algorithm>

#include "dataman.h"
#include "modalforce.h"

/* ModalForce - begin */

/* Costruttore */
ModalForce::ModalForce(unsigned int uL,
	const Modal *pmodal,
	const std::vector<unsigned int>& modeList,
	std::vector<DriveCaller *>& f,
	const Mat3xN *mt,
	const Mat3xN *mr,
	flag fOut)
: Elem(uL, fOut), 
Force(uL, fOut), 
pModal(pmodal),
modeList(modeList),
f(f),
Mt(mt),
Mr(mr),
F(Zero3),
M(Zero3)
{
	ASSERT(pModal != 0);
}

ModalForce::~ModalForce(void)
{
	if (!f.empty()) {
		for (unsigned i = 0; i < f.size(); i++) {
			delete f[i];
		}
	}

	if (Mt) {
		delete Mt;
	}

	if (Mr) {
		delete Mr;
	}
}

SubVectorHandler&
ModalForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iR, iC;
	WorkSpaceDim(&iR, &iC);
	WorkVec.ResizeReset(iR);

	integer iIdx = 1;
	const StructNode *pNode = pModal->pGetModalNode();
	integer iFirstIndex = 0;
	if (pNode) {
		iFirstIndex = pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		F = Zero3;
		M = Zero3;

		iIdx += 6;

	}

	integer iModalIndex = pModal->iGetModalIndex() + pModal->uGetNModes();
	for (unsigned iMode = 0; iMode < modeList.size(); iMode++) {
		doublereal d = f[iMode]->dGet();
		WorkVec.PutItem(iIdx + iMode, iModalIndex + modeList[iMode], d);
		if (pNode) {
			F += Mt->GetVec(iMode + 1)*d;
			M += Mr->GetVec(iMode + 1)*d;
		}
	}

	if (pNode) {
		WorkVec.Put(1, pNode->GetRCurr()*F);
		WorkVec.Put(4, pNode->GetRCurr()*M);
	}

	return WorkVec;
}

void
ModalForce::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Forces();

		out << GetLabel();

		if (pModal->pGetModalNode()) {
			out << " " << F << " " << M;
		}

		for (std::vector<DriveCaller *>::const_iterator i = f.begin(); i != f.end(); ++i) {
			out << " " << (*i)->dGet();
		}
		out << std::endl;
	}
}
   
Elem*
ReadModalForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel)
{
	const Modal *pModal = pDM->ReadElem<const Modal, const Joint, Elem::JOINT>(HP);

	std::vector<unsigned int> modeList;
	if (HP.IsKeyWord("list")) {
		const std::vector<unsigned int>& uModeList = pModal->GetModeList();
		std::vector<bool> gotIt(pModal->uGetNModes());
		for (integer i = 0; i < pModal->uGetNModes(); i++) {
			gotIt[i] = false;
		}

		int iNumModes = HP.GetInt();
		if (iNumModes <= 0 || iNumModes > pModal->uGetNModes()) {
			silent_cerr("ModalForce(" << uLabel << "): "
				"illegal mode number " << iNumModes
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		modeList.resize(iNumModes);

		for (int i = 0; i < iNumModes; i++) {
			int iM = HP.GetInt();
			if (iM <= 0) {
				silent_cerr("ModalForce(" << uLabel << "): "
					"illegal mode " << iM
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			std::vector<unsigned int>::const_iterator
				iv = std::find(uModeList.begin(), uModeList.end(), (unsigned int)iM);
			if (iv == uModeList.end()) {
				silent_cerr("ModalForce(" << uLabel << "): "
					"mode " << iM << " not active "
					"in Modal(" << pModal->GetLabel() << ")"
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			int iModeIdx = iv - uModeList.begin();

			if (gotIt[iModeIdx]) {
				silent_cerr("ModalForce(" << uLabel << "): "
					"mode " << iModeIdx << " already set "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			modeList[i] = iModeIdx + 1;
			gotIt[iModeIdx] = true;
		}

	} else {
		modeList.resize(pModal->uGetNModes());
		for (integer i = 0; i < pModal->uGetNModes(); i++) {
			modeList[i] = i + 1;
		}
	}

	std::vector<DriveCaller *> f(modeList.size());
	Mat3xN *Mt = 0;
	Mat3xN *Mr = 0;
	const StructNode *pNode = pModal->pGetModalNode();
	if (pNode) {
		Mt = new Mat3xN(modeList.size(), 0.);
		Mr = new Mat3xN(modeList.size(), 0.);
	}
	for (unsigned i = 0; i < f.size(); i++) {
		f[i] = HP.GetDriveCaller(false);
		if (f[i] == 0) {
			silent_cerr("ModalForce(" << uLabel << "): "
				"unable to read DriveCaller for mode #" << i + 1
				<< " (mode number " << modeList[i] << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (pNode && HP.IsKeyWord("resultant")) {
			for (unsigned r = 1; r <= 3; r++) {
				(*Mt)(r, i + 1) = HP.GetReal();
			}

			for (unsigned r = 1; r <= 3; r++) {
				(*Mr)(r, i + 1) = HP.GetReal();
			}
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalForce,
		ModalForce(uLabel, pModal, modeList, f, Mt, Mr, fOut));

	return pEl;
}

