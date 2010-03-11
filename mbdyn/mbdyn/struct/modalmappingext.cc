/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2010
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

// #include <fstream>

#include "dataman.h"
#include "extedge.h"
#include "modaledge.h"
#include "modalmappingext.h"

/* ModalMappingExt - begin */

/* Costruttore */
ModalMappingExt::ModalMappingExt(unsigned int uL,
	DataManager *pDM,
	StructNode *pRefNode,
	std::vector<StructNode *>& n,
	SpMapMatrixHandler *pH,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	ExtModalForceBase* pEMF,
	bool bSendAfterPredict,
	int iCoupling,
	ExtModalForceBase::BitMask bm,
	bool bUseRigidBodyForces,
	bool bRotateRigidBodyForces,
	flag fOut)
: Elem(uL, fOut),
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut),
pEMF(pEMF),
uFlags(ExtModalForceBase::EMF_NONE),
bOutputAccelerations(bOutputAccelerations),
bUseRigidBodyForces(bUseRigidBodyForces),
bRotateRigidBodyForces(bRotateRigidBodyForces),
pRefNode(pRefNode),
pH(pH),
F(0.),
M(0.)
{
	ASSERT(pEMF != 0);

	ASSERT(pH != 0);

	ASSERT((bm & ExtModalForceBase::EMF_ALL) != 0);
	ASSERT((bm & ~ExtModalForceBase::EMF_ALL) == 0);

	ASSERT(n.size() > 0);

	if (bm & ExtModalForceBase::EMF_MODAL) {
		f.resize(6*n.size());
		fill(f.begin(), f.end(), 0.);
		p.resize(pH->iGetNumRows());
		fill(p.begin(), p.end(), 0.);

		x.resize(6*n.size());
		fill(x.begin(), x.end(), 0.);
		xP.resize(6*n.size());
		fill(xP.begin(), xP.end(), 0.);
		q.resize(pH->iGetNumRows());
		fill(q.begin(), q.end(), 0.);
		qP.resize(pH->iGetNumRows());
		fill(qP.begin(), qP.end(), 0.);

		uFlags |= ExtModalForceBase::EMF_MODAL;

	} else {
		f.resize(0);
		p.resize(0);

		x.resize(0);
		xP.resize(0);
		q.resize(0);
		qP.resize(0);
	}

	if (bm & ExtModalForceBase::EMF_RIGID) {
		ASSERT(pRefNode != 0);
		uFlags |= ExtModalForceBase::EMF_RIGID;
	}

	Nodes.resize(n.size());
	for (unsigned i = 0; i < n.size(); i++) {
		Nodes[i].pNode = n[i];
		ASSERT(Nodes[i].pNode != 0);

		if (pRefNode) {
			Nodes[i].X0 = pRefNode->GetRCurr().MulTV(Nodes[i].pNode->GetXCurr() - pRefNode->GetXCurr());
			Nodes[i].R0 = pRefNode->GetRCurr().MulTM(Nodes[i].pNode->GetRCurr());

		} else {
			Nodes[i].X0 = Nodes[i].pNode->GetXCurr();
			Nodes[i].R0 = Nodes[i].pNode->GetRCurr();
		}

		Nodes[i].F = Zero3;
		Nodes[i].M = Zero3;
	}
}

ModalMappingExt::~ModalMappingExt(void)
{
	if (pEMF) {
		SAFEDELETE(pEMF);
	}

	if (pH) {
		SAFEDELETE(pH);
	}
}

bool
ModalMappingExt::Prepare(ExtFileHandlerBase *pEFH)
{
	return pEMF->Prepare(pEFH,
		uFlags & ExtModalForceBase::EMF_RIGID,
		pH->iGetNumRows());
}

/*
 * Send output to companion software
 */
void
ModalMappingExt::Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when)
{
	Vec3 X;
	Mat3x3 R;
	Vec3 V;
	Vec3 W;

	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		X = pRefNode->GetXCurr();
		R = pRefNode->GetRCurr();
		V = pRefNode->GetVCurr();
		W = pRefNode->GetWCurr();

	} else {
		X = Zero3;
		R = Eye3;
		V = Zero3;
		W = Zero3;
	}

	if (pRefNode) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 XRel = Nodes[i].pNode->GetXCurr() - X;

			x.Put(6*i + 1, R.MulTV(XRel) - Nodes[i].X0);
			x.Put(6*i + 4, MatR2LinParam(R.MulTM(Nodes[i].pNode->GetRCurr().MulMT(Nodes[i].R0))));

			xP.Put(6*i + 1, R.MulTV(Nodes[i].pNode->GetVCurr() - V + XRel.Cross(W)));
			xP.Put(6*i + 4, R.MulTV(Nodes[i].pNode->GetWCurr() - W));
		}

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			x.Put(6*i + 1, Nodes[i].pNode->GetXCurr() - Nodes[i].X0);
			x.Put(6*i + 4, MatR2LinParam(Nodes[i].pNode->GetRCurr().MulMT(Nodes[i].R0)));

			xP.Put(6*i + 1, Nodes[i].pNode->GetVCurr());
			xP.Put(6*i + 4, Nodes[i].pNode->GetWCurr());
		}
	}

	pH->MatVecMul(q, x);
	pH->MatVecMul(qP, xP);

	pEMF->Send(pEFH, uFlags, GetLabel(), X, R, V, W, q, qP);
}

void
ModalMappingExt::Recv(ExtFileHandlerBase *pEFH)
{
	unsigned uLabel = 0;
	unsigned uOutFlags = pEMF->Recv(pEFH, uFlags, uLabel, F, M, p);

	if (uOutFlags & ExtModalForceBase::EMF_ERR) {
		silent_cerr("ModalMappingExt(" << GetLabel() << "): "
			"error while reading modal force data" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (uFlags != uOutFlags) {
		silent_cerr("ModalMappingExt(" << GetLabel() << "): "
			"error while reading modal force data "
			"(" << uOutFlags << "!=" << uFlags << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pH->MatTVecMul(f, p);

	if (pRefNode) {
		const Vec3& XRef(pRefNode->GetXCurr());
		const Mat3x3& RRef(pRefNode->GetRCurr());

		if (bUseRigidBodyForces) {
			// use rigid body forces as provided by peer
			if (bRotateRigidBodyForces) {
				F = RRef*F;
				M = RRef*M;
			}

		} else {
			// re-compute rigid body forces using nodal forces
			F = Zero3;
			M = Zero3;
		}

		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 FTmp(&f[6*i]);
			Vec3 MTmp(&f[6*i + 3]);

			Nodes[i].F = RRef*FTmp;
			Nodes[i].M = RRef*MTmp;

			if (!bUseRigidBodyForces && Nodes[i].pNode != pRefNode) {
				F += Nodes[i].F;
				M += Nodes[i].M + (Nodes[i].pNode->GetXCurr() - XRef).Cross(Nodes[i].F);
			}
		}

		

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			Nodes[i].F = Vec3(&f[6*i]);
			Nodes[i].M = Vec3(&f[6*i + 3]);
		}
	}
}

Force::Type
ModalMappingExt::GetForceType(void) const
{
	return Force::EXTERNALMODAL;
}

void
ModalMappingExt::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = (pRefNode ? 6 : 0) + 6*Nodes.size();
	*piNumCols = 1;
}

SubVectorHandler&
ModalMappingExt::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	integer iR, iC;
	WorkSpaceDim(&iR, &iC);
	WorkVec.ResizeReset(iR);

	integer iIdx = 0;
	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		integer iFirstIndex = pRefNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Put(1, F);
		WorkVec.Put(4, M);
	}

	if (pRefNode) {
		iIdx += 6;
	}

	if (uFlags & ExtModalForceBase::EMF_MODAL) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			integer iIndex = Nodes[i].pNode->iGetFirstMomentumIndex();
			for (integer iCnt = 1; iCnt <= 6; iCnt++) {
				WorkVec.PutRowIndex(iIdx + iCnt, iIndex + iCnt);
			}

			WorkVec.Put(iIdx + 1, Nodes[i].F);
			WorkVec.Put(iIdx + 4, Nodes[i].M);

			iIdx += 6;
		}
	}

	return WorkVec;
}

void
ModalMappingExt::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Forces();

		if (uFlags & ExtModalForceBase::EMF_RIGID) {
			out << GetLabel() << "@" << pRefNode->GetLabel()
				<< " " << F << " " << M
				<< std::endl;
		}

		if (uFlags & ExtModalForceBase::EMF_MODAL) {
			for (unsigned i = 0; i < p.size(); i++) {
				out << GetLabel() << '.' << i + 1
					<< " " << p[i]
					<< " " << q[i]
					<< " " << qP[i]
					<< std::endl;
			}

			for (unsigned i = 0; i < Nodes.size(); i++) {
				out << GetLabel() << '#' << Nodes[i].pNode->GetLabel()
					<< " " << Nodes[i].F
					<< " " << Nodes[i].M
					<< std::endl;
			}
		}
	}
}

void
ModalMappingExt::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	unsigned iOff = 0;
	if (pRefNode) {
		iOff = 1;
	}

	connectedNodes.resize(iOff + Nodes.size());

	if (iOff) {
		connectedNodes[0] = pRefNode;
	}

	for (unsigned i = 0; i < Nodes.size(); i++) {
		connectedNodes[iOff + i] = Nodes[i].pNode;
	}
}

Elem*
ReadModalMappingExtForce(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH = 0;
	int iCoupling = -1;
	bool bSendAfterPredict = false;
	ReadExtForce(pDM, HP, uLabel, pEFH, bSendAfterPredict, iCoupling);

	StructNode *pRefNode = 0;
	if (HP.IsKeyWord("reference" "node")) {
		pRefNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (pRefNode == 0) {
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"illegal reference node "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("nodes" "number")) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"\"nodes number\" keyword expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	int nNodes = HP.GetInt();
	if (nNodes <= 0) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"invalid nodes number "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::vector<StructNode *> n(nNodes);
	for (int i = 0; i < nNodes; i++) {
		n[i] = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (n[i] == 0) {
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"illegal node #" << i << "/" << nNodes << " "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("modes" "number")) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"\"modes number\" keyword expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	int nModes = HP.GetInt();
	if (nModes <= 0) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"invalid modes number "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	bool bSparse;
	if (HP.IsKeyWord("full" "mapping" "file")) {
		bSparse = false;

	} else if (HP.IsKeyWord("sparse" "mapping" "file")) {
		bSparse = true;

	} else {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"mapping file expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	const char *sFileName = HP.GetFileName();
	if (sFileName == 0) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"unable to read mapping file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"unable to open mapping file "
			"\"" << sFileName << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	char c = in.get();
	while (c == '#') {
		do {
			c = in.get();
		} while (c != '\n');
		c = in.get();
	}
	in.putback(c);

	SpMapMatrixHandler *pH = 0;
	SAFENEWWITHCONSTRUCTOR(pH, SpMapMatrixHandler,
		SpMapMatrixHandler(nModes, 6*nNodes));

	if (bSparse) {
		integer ir, ic, cnt = 0;
		doublereal d;
		while (in >> ir >> ic >> d) {
			if (ir <= 0 || ir > nModes || ic <= 0 || ic > 6*nNodes) {
				silent_cerr("ModalMappingExt(" << uLabel << "): "
					"invalid row(=" << ir << ")/col(=" << ic <<") number for coefficient #" << cnt << " "
					"from file \"" << sFileName << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (d != 0.) {
				(*pH)(ir, ic) = d;
			}

			cnt++;
		}

		pedantic_cout("ModalMappingExt(" << uLabel << "): "
			"got " << cnt << " nonzeros from file \"" << sFileName << "\"" << std::endl);

	} else {
		for (integer ir = 1; ir <= nModes; ir++) {
			for (integer ic = 1; ic <= 6*nNodes; ic++) {
				doublereal d;
				in >> d;
				if (!in) {
					silent_cerr("ModalMappingExt(" << uLabel << "): "
						"unable to read coefficient(" << ir << ", " << ic << ") "
						"from file \"" << sFileName << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				if (d != 0.) {
					(*pH)(ir, ic) = d;
				}
			}
		}
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
			silent_cerr("ModalMappingExt(" << uLabel << "): unknown ModalMappingExt type "
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	ExtModalForceBase *pEMF = 0;
	if (dynamic_cast<ExtFileHandlerEDGE *>(pEFH) != 0) {
		// EDGE needs two separate ModalMappingExt elements,
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
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"EDGE ExtFileHandler can only be used "
				"when ModalMappingExt is either \"rigid\" "
				"or \"modal\" but not when it is \"all\" "
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
			ASSERT(0);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#ifdef USE_SOCKET
	} else if (dynamic_cast<ExtSocketHandler *>(pEFH) != 0) {
		SAFENEW(pEMF, ExtModalForce);
#endif // USE_SOCKET

	// add more types

	} else {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"unknown external force type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	bool bUseRigidBodyForces(true);
	bool bRotateRigidBodyForces(true);
	if (bm & ExtModalForceBase::EMF_RIGID && HP.IsKeyWord("use" "rigid" "body" "forces")) {
		bUseRigidBodyForces = HP.GetYesNo(bUseRigidBodyForces);

		if (bUseRigidBodyForces && HP.IsKeyWord("rotate" "rigid" "body" "forces")) {
			bRotateRigidBodyForces = HP.GetYesNo(bRotateRigidBodyForces);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);

	if (HP.IsArg()) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalMappingExt,
		ModalMappingExt(uLabel, pDM, pRefNode, n, pH,
			bOutputAccelerations,
			pEFH, pEMF, bSendAfterPredict, iCoupling, bm,
			bUseRigidBodyForces, bRotateRigidBodyForces, fOut));

	return pEl;
}

