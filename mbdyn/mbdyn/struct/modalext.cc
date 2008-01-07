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
#include "modalext.h"

#include <fstream>

/* ModalExt - begin */

/* Costruttore */
ModalExt::ModalExt(unsigned int uL,
	Modal *pmodal,
	bool bOutputAccelerations,
	std::string& fin,
	bool bRemoveIn,
        std::string& fout,
	bool bNoClobberOut,
	int iSleepTime,
	int iCoupling,
	int iPrecision,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, fin, bRemoveIn, fout, bNoClobberOut, iSleepTime, iCoupling, iPrecision, fOut), 
pModal(pmodal),
F(0.),
M(0.),
bOutputAccelerations(bOutputAccelerations)
{
	ASSERT(pModal != 0);
	f.resize(pModal->uGetNModes());
}

ModalExt::~ModalExt(void)
{
	NO_OP;
}

/*
 * Send output to companion software
 */
void
ModalExt::Send(std::ostream& outf)
{
	const StructNode *pNode = pModal->pGetModalNode();
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
}

void
ModalExt::Recv(std::istream& inf)
{
	/* assume unsigned int label */
	unsigned l;

	inf >> l >> F(1) >> F(2) >> F(3) >> M(1) >> M(2) >> M(3);
	if (!inf) {
		silent_cerr("ModalExt(" << GetLabel() << "): "
			"unable to read reference node forces "
			"from file \"" << fin << "\"" << std::endl);
		throw ErrGeneric();
	}

	for (integer i = 0; i < pModal->uGetNModes(); i++) {
		inf >> f[i];
		if (!inf) {
			silent_cerr("ModalExt(" << GetLabel() << "): "
				"unable to read modal force #" << i
				<< " from file \"" << fin << "\"" << std::endl);
			throw ErrGeneric();
		}
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
	std::string fin;
	bool bUnlinkIn;
	std::string fout;
	bool bNoClobberOut;
	int iSleepTime;
	int iCoupling;
	int iPrecision;
	
	ReadExtForce(pDM, HP, uLabel, fin, bUnlinkIn, fout, bNoClobberOut, iSleepTime, iCoupling, iPrecision);

	Modal *pModal = dynamic_cast<Modal *>(pDM->ReadElem(HP, Elem::JOINT));
	if (pModal == 0) {
		silent_cerr("ModalExt(" << uLabel << "): illegal Modal joint "
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	bool bOutputAccelerations(false);
	if (HP.IsKeyWord("accelerations")) {
		bOutputAccelerations = true;
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalExt,
		ModalExt(uLabel, pModal, bOutputAccelerations,
			fin, bUnlinkIn, fout, bNoClobberOut,
			iSleepTime, iCoupling, iPrecision, fOut));

	return pEl;
}

