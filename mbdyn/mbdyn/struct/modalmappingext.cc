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
#include "modaledge.h"
#include "modalmappingext.h"

/* ModalMappingExt - begin */

/* Costruttore */
ModalMappingExt::ModalMappingExt(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	std::vector<const StructNode *>& n,
	SpMapMatrixHandler *pH,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	ExtModalForceBase* pEMF,
	bool bSendAfterPredict,
	int iCoupling,
	ExtModalForceBase::BitMask bm,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	flag fOut)
: Elem(uL, fOut),
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut),
pEMF(pEMF),
uFlags(ExtModalForceBase::EMF_NONE),
bOutputAccelerations(bOutputAccelerations),
bUseReferenceNodeForces(bUseReferenceNodeForces),
bRotateReferenceNodeForces(bRotateReferenceNodeForces),
pRefNode(pRefNode),
pH(pH),
F0(Zero3), M0(Zero3),
F1(Zero3), M1(Zero3),
F2(Zero3), M2(Zero3)
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
	return pEMF->Prepare(pEFH, GetLabel(),
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
	unsigned uOutFlags = pEMF->Recv(pEFH, uFlags, uLabel, F0, M0, p);

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

		if (bUseReferenceNodeForces) {
			// initialize rigid body forces with values provided by peer
			if (bRotateReferenceNodeForces) {
				F1 = RRef*F0;
				M1 = RRef*M0;

			} else {
				F1 = F0;
				M1 = M0;
			}
		}

		F2 = Zero3;
		M2 = Zero3;
		for (unsigned i = 0; i < Nodes.size(); i++) {
			Nodes[i].F = RRef*Vec3(&f[6*i]);
			Nodes[i].M = RRef*Vec3(&f[6*i + 3]);

			if ((bUseReferenceNodeForces || fToBeOutput()) && Nodes[i].pNode != pRefNode) {
				F2 += Nodes[i].F;
				M2 += Nodes[i].M + (Nodes[i].pNode->GetXCurr() - XRef).Cross(Nodes[i].F);
			}
		}

		if (bUseReferenceNodeForces) {
			F1 -= F2;
			M1 -= M2;
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
	if (iCoupling == COUPLING_NONE) {
		*piNumRows = 0;
		*piNumCols = 0;

	} else {
		*piNumRows = (pRefNode ? 6 : 0) + 6*Nodes.size();
		*piNumCols = 1;
	}
}

SubVectorHandler&
ModalMappingExt::AssRes(SubVectorHandler& WorkVec,
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

	integer iIdx = 0;
	if (uFlags & ExtModalForceBase::EMF_RIGID) {
		integer iFirstIndex = pRefNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Put(1, F1);
		WorkVec.Put(4, M1);
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
				<< " " << F0 << " " << M0
				<< " " << F1 << " " << M1
				<< " " << F2 << " " << M2
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
				const double *px = &x[6*i];
				const double *pxP = &xP[6*i];
				out << GetLabel() << '#' << Nodes[i].pNode->GetLabel()
					<< " " << Nodes[i].F
					<< " " << Nodes[i].M
					<< " " << px[0]
					<< " " << px[1]
					<< " " << px[2]
					<< " " << px[3]
					<< " " << px[4]
					<< " " << px[5]
					<< " " << pxP[0]
					<< " " << pxP[1]
					<< " " << pxP[2]
					<< " " << pxP[3]
					<< " " << pxP[4]
					<< " " << pxP[5]
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

SpMapMatrixHandler *
ReadSparseMappingMatrix(MBDynParser& HP, integer& nRows, integer& nCols)
{
	bool bSparse;
	if (HP.IsKeyWord("full" "mapping" "file")) {
		bSparse = false;

	} else if (HP.IsKeyWord("sparse" "mapping" "file")) {
		bSparse = true;

	} else {
		silent_cerr("mapping file expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dThreshold = 0.;
	if (HP.IsKeyWord("threshold")) {
		dThreshold = HP.GetReal();
		if (dThreshold < 0.) {
			silent_cerr("invalid threshold " << dThreshold
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	const char *sFileName = HP.GetFileName();
	if (sFileName == 0) {
		silent_cerr("unable to read mapping file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("unable to open mapping file "
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

	bool bComputeRows(nRows < 0);
	bool bComputeCols(nCols < 0);

	if (bComputeRows || bComputeCols) {
		std::streampos pos = in.tellg();
		integer nVals;

		if (bSparse) {
			for (nVals = 0; !in.eof(); nVals++) {
				integer ir, ic;
				doublereal d;
				in >> ir >> ic >> d;

				if (bComputeRows) {
					if (ir > nRows) {
						nRows = ir;
					}
	
				} else {
					if (ir > nRows) {
						silent_cerr("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
							"inconsistent row=" << ir << " for coefficient #" << nVals << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

				if (bComputeCols) {
					if (ic > nCols) {
						nCols = ic;
					}
	
				} else {
					if (ic > nCols) {
						silent_cerr("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
							"inconsistent col=" << ic << " for coefficient #" << nVals << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}
			nVals--;

		} else {
			if (bComputeRows && bComputeCols) {
				silent_cerr("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
					"cannot compute both row and col numbers" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			for (nVals = 0; !in.eof(); nVals++) {
				doublereal d;

				in >> d;
			}
			nVals--;

			if (bComputeRows) {
				if ((nVals % nCols) != 0) {
					silent_cerr("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
						"vals=" << nVals << " is not a multiple of cols=" << nCols << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				nRows = nVals/nCols;

			} else if (bComputeCols) {
				if ((nVals % nRows) != 0) {
					silent_cerr("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
						"vals=" << nVals << " is not a multiple of rows=" << nRows << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				nCols = nVals/nRows;
			}
		}

		in.clear();
		in.seekg(pos);

		silent_cout("ReadSparseMappingMatrix(\"" << sFileName << "\"): "
			"rows=" << nRows
			<< " cols=" << nCols
			<< " vals=" << nVals
			<< std::endl);
	}

	SpMapMatrixHandler *pH = 0;
	SAFENEWWITHCONSTRUCTOR(pH, SpMapMatrixHandler,
		SpMapMatrixHandler(nRows, nCols));

	if (bSparse) {
		integer ir, ic, cnt = 0, nzcnt = 0;
		doublereal d;
		while (in >> ir >> ic >> d) {
			if (ir <= 0 || ir > nRows || ic <= 0 || ic > nCols) {
				silent_cerr("invalid row(=" << ir << ")/col(=" << ic <<") number for coefficient #" << cnt << " "
					"from file \"" << sFileName << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (std::abs(d) > dThreshold) {
				(*pH)(ir, ic) = d;
				nzcnt++;
			}

			cnt++;
		}

		pedantic_cout("got " << cnt << " nonzeros (" << nzcnt << " actually stored) from file \"" << sFileName << "\"" << std::endl);

	} else {
		int nzcnt = 0;
		for (integer ir = 1; ir <= nRows; ir++) {
			for (integer ic = 1; ic <= nCols; ic++) {
				doublereal d;
				in >> d;
				if (!in) {
					silent_cerr("unable to read coefficient(" << ir << ", " << ic << ") "
						"from file \"" << sFileName << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				if (std::abs(d) > dThreshold) {
					(*pH)(ir, ic) = d;
					nzcnt++;
				}
			}
		}
		pedantic_cout("got " << nzcnt << " nonzeros from file \"" << sFileName << "\"" << std::endl);
	}

	return pH;
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

	const StructNode *pRefNode = 0;
	if (HP.IsKeyWord("reference" "node")) {
		pRefNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
		if (pRefNode == 0) {
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"illegal reference node "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool bUseReferenceNodeForces(pRefNode != 0 ? true : false);
	bool bRotateReferenceNodeForces(true);
	if (pRefNode != 0 && HP.IsKeyWord("use" "rigid" "body" "forces")) {
		if (!HP.GetYesNo(bUseReferenceNodeForces)) {
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"\"use rigid body forces\" must be either \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bUseReferenceNodeForces && HP.IsKeyWord("rotate" "rigid" "body" "forces")) {
			if (!HP.GetYesNo(bRotateReferenceNodeForces)) {
				silent_cerr("ModalMappingExt(" << uLabel << "): "
					"\"rotate rigid body forces\" must be either \"yes\" or \"no\" "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	std::vector<const StructNode *> n(nNodes);
	for (int i = 0; i < nNodes; i++) {
		n[i] = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
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
	int nModes = -1;
	if (!HP.IsKeyWord("from" "file")) {
		nModes = HP.GetInt();
		if (nModes <= 0) {
			silent_cerr("ModalMappingExt(" << uLabel << "): "
				"invalid modes number "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	integer nCols = 6*nNodes;
	SpMapMatrixHandler *pH = ReadSparseMappingMatrix(HP, nModes, nCols);
	ASSERT(nCols == 6*nNodes);

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);

	if (HP.IsArg()) {
		silent_cerr("ModalMappingExt(" << uLabel << "): "
			"semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		SAFEDELETE(pH);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ModalMappingExt,
		ModalMappingExt(uLabel, pDM, pRefNode, n, pH,
			bOutputAccelerations,
			pEFH, pEMF, bSendAfterPredict, iCoupling, bm,
			bUseReferenceNodeForces, bRotateReferenceNodeForces, fOut));

	return pEl;
}

