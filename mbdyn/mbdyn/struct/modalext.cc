/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2011
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

#include "dataman.h"
#include "extedge.h"
#include "modalext.h"
#include "modaledge.h"

#include <fstream>
#include <cerrno>

/* ExtModalForceBase - begin */

ExtModalForceBase::~ExtModalForceBase(void)
{
	NO_OP;
}

/* ExtModalForceBase - end */

/* ExtModalForce - begin */

ExtModalForce::ExtModalForce(void)
{
	NO_OP;
}

ExtModalForce::~ExtModalForce(void)
{
	NO_OP;
}

bool
ExtModalForce::Prepare(ExtFileHandlerBase *pEFH, unsigned uLabel, bool bRigid, unsigned uModes)
{
	bool bResult = true;

	switch (pEFH->NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_NO:
		break;

	case ExtFileHandlerBase::NEGOTIATE_CLIENT: {
		std::ostream *outfp = pEFH->GetOutStream();
		if (outfp) {
			*outfp << MBC_MODAL
				<< ' ' << bRigid
				<< ' ' << uModes
				<< std::endl;

#ifdef USE_SOCKET
		} else {
			char buf[sizeof(uint32_t) + sizeof(uint32_t)];
			uint32_t *uint32_ptr;

			uint32_ptr = (uint32_t *)&buf[0];
			uint32_ptr[0] = MBC_MODAL;
			if (bRigid) {
				uint32_ptr[0] |= MBC_REF_NODE;
			}

			uint32_ptr[1] = uModes;

			ssize_t rc = send(pEFH->GetOutFileDes(),
				(const void *)buf, sizeof(buf),
				pEFH->GetSendFlags());
			if (rc == -1) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("ExtModalForce(" << uLabel << "): "
					"negotiation request send() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("ExtModalForce(" << uLabel << "): "
					"negotiation request send() failed "
					"(sent " << rc << " of " << sizeof(buf) << " bytes)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
#endif // USE_SOCKET
		}
		} break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER: {
		unsigned type;
		bool bR;
		unsigned uM;

		std::istream *infp = pEFH->GetInStream();
		if (infp) {
			*infp >> type >> bR >> uM;
			if (type != MBC_MODAL) {
				return false;
			}

#ifdef USE_SOCKET
		} else {
			char buf[sizeof(uint32_t) + sizeof(uint32_t)];
			uint32_t *uint32_ptr;

			ssize_t rc = recv(pEFH->GetInFileDes(),
				(void *)buf, sizeof(buf),
				pEFH->GetRecvFlags());
			if (rc == -1) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("ExtModalForce(" << uLabel << "): "
					"negotiation response recv() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("ExtModalForce(" << uLabel << "): "
					"negotiation response recv() failed "
					"(got " << rc << " of " << sizeof(buf) << " bytes)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			uint32_ptr = (uint32_t *)&buf[0];
			type = (uint32_ptr[0] & MBC_MODAL_NODAL_MASK);
			bR = (uint32_ptr[0] & MBC_REF_NODE);

			uM = uint32_ptr[1];
#endif // USE_SOCKET
		}

		if (type != MBC_MODAL) {
			silent_cerr("ExtModalForce(" << uLabel << "): "
				"negotiation response failed: expecting MBC_MODAL "
				"(=" << MBC_MODAL << "), got " << type
				<< std::endl);
			bResult = false;
		}

		if ((bRigid && !bR) || (!bRigid && bR)) {
			silent_cerr("ExtModalForce(" << uLabel << "): "
				"negotiation response failed: reference node configuration mismatch "
				"(local=" << (bRigid ? "yes" : "no") << ", remote=" << (bR ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (uM != uModes) {
			silent_cerr("ExtModalForce(" << uLabel << "): "
				"negotiation response failed: mode number mismatch "
				"(local=" << uModes << ", remote=" << uM << ")"
				<< std::endl);
			bResult = false;
		}
		} break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return bResult;
}

unsigned
ExtModalForce::Recv(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned& uLabel,
	Vec3& f, Vec3& m, std::vector<doublereal>& fv)
{
	std::istream *infp = pEFH->GetInStream();
	if (infp) {
		return RecvFromStream(*infp, uFlags, uLabel, f, m, fv);

	} else {
		return RecvFromFileDes(pEFH->GetInFileDes(), pEFH->GetRecvFlags(),
			uFlags, uLabel, f, m, fv);
	}
}

void
ExtModalForce::Send(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned uLabel,
	const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
	const std::vector<doublereal>& q,
	const std::vector<doublereal>& qP)
{
	std::ostream *outfp = pEFH->GetOutStream();
	if (outfp) {
		return SendToStream(*outfp, uFlags, uLabel, x, R, v, w, q, qP);

	} else {
		return SendToFileDes(pEFH->GetOutFileDes(), pEFH->GetSendFlags(),
			uFlags, uLabel, x, R, v, w, q, qP);
	}
}

unsigned
ExtModalForce::RecvFromStream(std::istream& inf, unsigned uFlags, unsigned& uLabel,
	Vec3& f, Vec3& m, std::vector<doublereal>& fv)
{
	if ((uFlags & ExtModalForceBase::EMF_RIGID)) {
		doublereal df[3], dm[3];
		inf
			>> df[0] >> df[1] >> df[2]
			>> dm[0] >> dm[1] >> dm[2];
		f = Vec3(df);
		m = Vec3(dm);
	}

	if ((uFlags & ExtModalForceBase::EMF_MODAL)) {
		for (std::vector<doublereal>::iterator i = fv.begin(); i != fv.end(); ++i) {
			doublereal d;
			inf >> d;
			*i = d;
		}
	}

	return uFlags;
}

unsigned
ExtModalForce::RecvFromFileDes(int infd, int recv_flags,
	unsigned uFlags, unsigned& uLabel,
	Vec3& f, Vec3& m, std::vector<doublereal>& fv)
{
#ifdef USE_SOCKET
	ssize_t rc;
	size_t size;

	if ((uFlags & ExtModalForceBase::EMF_RIGID)) {
		size = 3*sizeof(doublereal);

		rc = recv(infd, (void *)f.pGetVec(), size, recv_flags);
		if (rc != (ssize_t)size) {
			// error
		}
		rc = recv(infd, (void *)m.pGetVec(), size, recv_flags);
		if (rc != (ssize_t)size) {
			// error
		}
	}

	if ((uFlags & ExtModalForceBase::EMF_MODAL)) {
		size = fv.size()*sizeof(doublereal);
		rc = recv(infd, (void *)&fv[0], size, recv_flags);
		if (rc != (ssize_t)size) {
			// error
		}
	}

	return uFlags;
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

void
ExtModalForce::SendToStream(std::ostream& outf, unsigned uFlags, unsigned uLabel,
	const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
	const std::vector<doublereal>& q,
	const std::vector<doublereal>& qP)
{
	if ((uFlags & ExtModalForceBase::EMF_RIGID)) {
		outf
			<< x << std::endl
			<< R << std::endl
			<< v << std::endl
			<< w << std::endl;
	}

	if ((uFlags & ExtModalForceBase::EMF_MODAL)) {
		for (unsigned i = 0; i < q.size(); i++) {
			outf << q[i] << " " << qP[i] << std::endl;
		}
	}
}

void
ExtModalForce::SendToFileDes(int outfd, int send_flags,
	unsigned uFlags, unsigned uLabel,
	const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
	const std::vector<doublereal>& q,
	const std::vector<doublereal>& qP)
{
#ifdef USE_SOCKET
	if ((uFlags & ExtModalForceBase::EMF_RIGID)) {
		send(outfd, (const void *)x.pGetVec(), 3*sizeof(doublereal), send_flags);
		send(outfd, (const void *)R.pGetMat(), 9*sizeof(doublereal), send_flags);
		send(outfd, (const void *)v.pGetVec(), 3*sizeof(doublereal), send_flags);
		send(outfd, (const void *)w.pGetVec(), 3*sizeof(doublereal), send_flags);
	}

	if ((uFlags & ExtModalForceBase::EMF_MODAL)) {
		send(outfd, (const void *)&q[0], q.size()*sizeof(doublereal), send_flags);
		send(outfd, (const void *)&qP[0], qP.size()*sizeof(doublereal), send_flags);
	}
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

/* ExtModalForce - end */

/* ModalExt - begin */

/* Costruttore */
ModalExt::ModalExt(unsigned int uL,
	DataManager *pDM,
	Modal *pmodal,
	const StructNode *pnode,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	ExtModalForceBase* pEMF,
	bool bSendAfterPredict,
	int iCoupling,
	ExtModalForceBase::BitMask bm,
	flag fOut)
: Elem(uL, fOut),
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut),
pModal(pmodal),
pNode(pnode),
bOutputAccelerations(bOutputAccelerations),
pEMF(pEMF),
uFlags(ExtModalForceBase::EMF_NONE),
F(Zero3),
M(Zero3)
{
	ASSERT((pModal != 0) || (pNode != 0));
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
		ASSERT(pNode != 0);
		uFlags |= ExtModalForceBase::EMF_RIGID;
	}
}

ModalExt::~ModalExt(void)
{
	if (pEMF) {
		SAFEDELETE(pEMF);
	}
}

bool
ModalExt::Prepare(ExtFileHandlerBase *pEFH)
{
	return pEMF->Prepare(pEFH, GetLabel(),
		uFlags & ExtModalForceBase::EMF_RIGID,
		pModal ? pModal->uGetNModes() : 0);
}

/*
 * Send output to companion software
 */
void
ModalExt::Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when)
{
	Vec3 x;
	Mat3x3 R;
	Vec3 v;
	Vec3 w;

	if (uFlags & ExtModalForceBase::EMF_RIGID) {
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
	if (uFlags & ExtModalForceBase::EMF_MODAL) {
		const VecN& a = pModal->GetA();
		const VecN& b = pModal->GetB();
		for (unsigned i = 0; i < q.size(); i++) {
			q[i] = a(i + 1);
			qP[i] = b(i + 1);
		}
	}

	pEMF->Send(pEFH, uFlags, GetLabel(), x, R, v, w, q, qP);

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
ModalExt::Recv(ExtFileHandlerBase *pEFH)
{
	unsigned uLabel = 0;
	unsigned uOutFlags = pEMF->Recv(pEFH, uFlags, uLabel, F, M, f);

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

Force::Type
ModalExt::GetForceType(void) const
{
	return Force::EXTERNALMODAL;
}

void
ModalExt::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	if (iCoupling == COUPLING_NONE) {
		*piNumRows = 0;
		*piNumCols = 0;

	} else {
		*piNumRows = (pNode ? 6 : 0) + (pModal ? pModal->uGetNModes() : 0);
		*piNumCols = 1;
	}
}

SubVectorHandler&
ModalExt::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	if (iCoupling == COUPLING_NONE) {
		WorkVec.Resize(0);
		return WorkVec;
	}

	integer iR, iC;
	WorkSpaceDim(&iR, &iC);
	WorkVec.ResizeReset(iR);

	integer iIdx = 1;
	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		integer iFirstIndex = pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Put(1, F);
		WorkVec.Put(4, M);
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

		if (uFlags & ExtModalForceBase::EMF_RIGID) {
			out << GetLabel() << "#" << pNode->GetLabel()
				<< " " << F << " " << M
				<< std::endl;
		}

		if (uFlags & ExtModalForceBase::EMF_MODAL) {
			unsigned cnt = 1;
			for (std::vector<doublereal>::const_iterator i = f.begin(); i != f.end(); ++i) {
				out << GetLabel() << '.' << cnt
					<< " " << *i << std::endl;
				cnt++;
			}
		}
	}
}

void
ModalExt::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	if (pNode) {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	} else {
		connectedNodes.resize(0);
	}
}

Elem*
ReadModalExtForce(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH;
	int iCoupling;

	bool bSendAfterPredict;
	ReadExtForce(pDM, HP, uLabel, pEFH, bSendAfterPredict, iCoupling);

	const StructNode *pNode = 0;
	Modal *pModal = 0;
	if (HP.IsKeyWord("reference" "node")) {
		pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (pNode == 0) {
			silent_cerr("ModalExt(" << uLabel << "): illegal reference node "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		pModal = dynamic_cast<Modal *>(pDM->ReadElem(HP, Elem::JOINT));
		if (pModal == 0) {
			silent_cerr("ModalExt(" << uLabel << "): illegal Modal joint "
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pNode = pModal->pGetModalNode();
	}

	bool bOutputAccelerations(false);
	if (HP.IsKeyWord("accelerations")) {
		bOutputAccelerations = true;
	}

	// safe default: both rigid and modal, if possible
	ExtModalForceBase::BitMask bm = ExtModalForceBase::EMF_ALL;
	if (pNode == 0) {
		bm = ExtModalForceBase::EMF_MODAL;
		
	} else if (pModal == 0) {
		bm = ExtModalForceBase::EMF_RIGID;
	}

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

	if ((bm & ExtModalForceBase::EMF_RIGID) && pNode == 0) {
		silent_cerr("ModalExt(" << uLabel << "): \"rigid\" needs rigid body motion "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if ((bm & ExtModalForceBase::EMF_MODAL) && pModal == 0) {
		silent_cerr("ModalExt(" << uLabel << "): \"modal\" needs modal joint "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	// add more types if specific behavior is needed

	} else {
		// ExtFileHandler
		// ExtSocketHandler

		// wild guess: the default is fine
		SAFENEW(pEMF, ExtModalForce);
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
		ModalExt(uLabel, pDM, pModal, pNode, bOutputAccelerations,
			pEFH, pEMF, bSendAfterPredict, iCoupling, bm, fOut));

	return pEl;
}

