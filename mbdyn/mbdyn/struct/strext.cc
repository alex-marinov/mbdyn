/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2009
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
#include "strext.h"

#include <fstream>
#include <cerrno>

/* StructExtForce - begin */

/* Costruttore */
StructExtForce::StructExtForce(unsigned int uL,
	DataManager *pDM,
	std::vector<StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	bool bUnsorted,
	bool bNoLabels,
	bool bOutputAccelerations,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut), 
pRefNode(0),
RefOffset(0.),
bUnsorted(bUnsorted),
bNoLabels(bNoLabels),
bOutputAccelerations(bOutputAccelerations)
{
	ASSERT(nodes.size() == offsets.size());
	Nodes.resize(nodes.size());
	Offsets.resize(nodes.size());
	F.resize(nodes.size());
	M.resize(nodes.size());

	for (unsigned i = 0; i < nodes.size(); i++) {
		Nodes[i] = nodes[i];
		Offsets[i] = offsets[i];
		F[i] = Zero3;
		M[i] = Zero3;
	}

	ASSERT(!(bNoLabels && bUnsorted));
	if (bUnsorted) {
		done.resize(nodes.size());
	}

	if (bOutputAccelerations) {
		for (unsigned i = 0; i < nodes.size(); i++) {
			DynamicStructNode *pDSN = dynamic_cast<DynamicStructNode *>(Nodes[i]);
			if (pDSN == 0) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"StructNode(" << Nodes[i]->GetLabel() << ") "
					"is not dynamic"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			pDSN->ComputeAccelerations(true);
		}
	}
}

StructExtForce::~StructExtForce(void)
{
	NO_OP;
}


/*
 * Send output to companion software
 */
void
StructExtForce::Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when)
{
	std::ostream *outfp = pEFH->GetOutStream();
	if (outfp) {
		SendToStream(*outfp, when);

	} else {
		SendToFileDes(pEFH->GetOutFileDes(), when);
	}
}

void
StructExtForce::SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when)
{
	if (pRefNode) {
#if 0
		// TODO
		Vec3 fRef = pRefNode->GetRCurr()*RefOffset;
		Vec3 xRef = pRefNode->GetXCurr() + fRef;
		Mat3x3 RRefT = pRefNode->GetRCurr().Transpose();

		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 f = Nodes[i]->GetRCurr()*Offsets[i];
			Vec3 x = Nodes[i]->GetXCurr() + f;
			Vec3 v = Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(f);

			// manipulate

			outf << Nodes[i]->GetLabel()
				<< " " << x
				<< " " << Nodes[i]->GetRCurr()
				<< " " << v
				<< " " << Nodes[i]->GetWCurr()
				<< std::endl;
		}
#endif

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			/*
				p = x + f
				R = R
				v = xp + w cross f
				w = w
				a = xpp + wp cross f + w cross w cross f
				wp = wp
			 */

			// Optimization of the above formulas
			const Mat3x3& R = Nodes[i]->GetRCurr();
			Vec3 f = R*Offsets[i];
			Vec3 x = Nodes[i]->GetXCurr() + f;
			const Vec3& w = Nodes[i]->GetWCurr();
			Vec3 wCrossf = w.Cross(f);
			Vec3 v = Nodes[i]->GetVCurr() + wCrossf;

			outf << Nodes[i]->GetLabel()
				<< " " << x
				<< " " << R
				<< " " << v
				<< " " << w;

			if (bOutputAccelerations) {
				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 a = Nodes[i]->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				outf
					<< " " << a
					<< " " << wp;
			}

			outf << std::endl;
		}
	}
}

void
StructExtForce::SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when)
{
	if (pRefNode) {
#if 0
		// TODO
		Vec3 fRef = pRefNode->GetRCurr()*RefOffset;
		Vec3 xRef = pRefNode->GetXCurr() + fRef;
		Mat3x3 RRefT = pRefNode->GetRCurr().Transpose();

		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 f = Nodes[i]->GetRCurr()*Offsets[i];
			Vec3 x = Nodes[i]->GetXCurr() + f;
			Vec3 v = Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(f);

			// manipulate

			outf << Nodes[i]->GetLabel()
				<< " " << x
				<< " " << Nodes[i]->GetRCurr()
				<< " " << v
				<< " " << Nodes[i]->GetWCurr()
				<< std::endl;
		}
#endif

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			/*
				p = x + f
				R = R
				v = xp + w cross f
				w = w
				a = xpp + wp cross f + w cross w cross f
				wp = wp
			 */

			unsigned l = Nodes[i]->GetLabel();
			// Optimization of the above formulas
			const Mat3x3& R = Nodes[i]->GetRCurr();
			Vec3 f = R*Offsets[i];
			Vec3 x = Nodes[i]->GetXCurr() + f;
			const Vec3& w = Nodes[i]->GetWCurr();
			Vec3 wCrossf = w.Cross(f);
			Vec3 v = Nodes[i]->GetVCurr() + wCrossf;

			send(outfd, (void *)&l, sizeof(l), 0);
			send(outfd, (void *)x.pGetVec(), 3*sizeof(doublereal), 0);
			send(outfd, (void *)R.pGetMat(), 9*sizeof(doublereal), 0);
			send(outfd, (void *)v.pGetVec(), 3*sizeof(doublereal), 0);
			send(outfd, (void *)w.pGetVec(), 3*sizeof(doublereal), 0);

			if (bOutputAccelerations) {
				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 a = Nodes[i]->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				send(outfd, (void *)a.pGetVec(), 3*sizeof(doublereal), 0);
				send(outfd, (void *)wp.pGetVec(), 3*sizeof(doublereal), 0);
			}
		}
	}
}

void
StructExtForce::Recv(ExtFileHandlerBase *pEFH)
{
	std::istream *infp = pEFH->GetInStream();
	if (infp) {
		RecvFromStream(*infp);

	} else {
		RecvFromFileDes(pEFH->GetInFileDes());
	}
}

void
StructExtForce::RecvFromStream(std::istream& inf)
{
	if (bUnsorted) {
		done.resize(Nodes.size());

		for (unsigned i = 0; i < Nodes.size(); i++) {
			done[i] = false;
		}

		unsigned cnt;
		for (cnt = 0; inf; cnt++) {
			/* assume unsigned int label */
			unsigned l, i;
			doublereal f[3], m[3];

			inf >> l
				>> f[0] >> f[1] >> f[2]
				>> m[0] >> m[1] >> m[2];

			if (!inf) {
				break;
			}

			for (i = 0; i < Nodes.size(); i++) {
				if (Nodes[i]->GetLabel() == l) {
					break;
				}
			}

			if (i == Nodes.size()) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"unknown label " << l
					<< " as " << cnt << "-th node"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (done[i]) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"label " << l << " already done"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			done[i] = true;

			F[i] = Vec3(f);
			M[i] = Vec3(m);
		}

		if (cnt != Nodes.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"invalid number of nodes " << cnt
				<< std::endl);

			for (unsigned i = 0; i < Nodes.size(); i++) {
				if (!done[i]) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"node " << Nodes[i]->GetLabel()
						<< " not done" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			/* assume unsigned int label */
			unsigned l;
			doublereal f[3], m[3];

			if (!bNoLabels) {
				inf >> l;

				if (Nodes[i]->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			inf
				>> f[0] >> f[1] >> f[2]
				>> m[0] >> m[1] >> m[2];

			if (!inf) {
				break;
			}

			F[i] = Vec3(f);
			M[i] = Vec3(m);
		}
	}
}

void
StructExtForce::RecvFromFileDes(int infd)
{
	if (bUnsorted) {
		done.resize(Nodes.size());

		for (unsigned i = 0; i < Nodes.size(); i++) {
			done[i] = false;
		}

		unsigned cnt;
		for (cnt = 0; cnt < Nodes.size(); cnt++) {
			/* assume unsigned int label */
			unsigned l, i;
			doublereal f[6];
			ssize_t len;

			len = recv(infd, (void *)&l, sizeof(l), 0);
			if (len != sizeof(l)) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"recv() failed (" << save_errno << ": "
					<< err_msg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			len = recv(infd, (void *)f, sizeof(f), MSG_DONTWAIT);
			if (len == -1) {
				int save_errno = errno;
				if (errno == EAGAIN) {
					// would block
					break;
				}

				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"recv() failed (" << save_errno << ": "
					<< err_msg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (len != sizeof(f)) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"recv() failed (" << save_errno << ": "
					<< err_msg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			for (i = 0; i < Nodes.size(); i++) {
				if (Nodes[i]->GetLabel() == l) {
					break;
				}
			}

			if (i == Nodes.size()) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"unknown label " << l
					<< " as " << cnt << "-th node"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (done[i]) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"label " << l << " already done"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			done[i] = true;

			F[i] = Vec3(&f[0]);
			M[i] = Vec3(&f[3]);
		}

		if (cnt != Nodes.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"invalid node number " << cnt
				<< std::endl);

			for (unsigned i = 0; i < Nodes.size(); i++) {
				if (!done[i]) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"node " << Nodes[i]->GetLabel()
						<< " not done" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			/* assume unsigned int label */
			unsigned l;
			doublereal f[6];
			ssize_t len;

			if (!bNoLabels) {
				len = recv(infd, (void *)&l, sizeof(l), 0);
				if (len != sizeof(l)) {
					int save_errno = errno;
					char *err_msg = strerror(save_errno);
					silent_cerr("StructExtForce(" << GetLabel() << "): "
						"recv() failed (" << save_errno << ": "
						<< err_msg << ")" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (Nodes[i]->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			len = recv(infd, (void *)f, sizeof(f), MSG_DONTWAIT);
			if (len == -1) {
				int save_errno = errno;
				if (errno == EAGAIN) {
					// would block
					break;
				}

				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"recv() failed (" << save_errno << ": "
					<< err_msg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (len != sizeof(f)) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"recv() failed (" << save_errno << ": "
					<< err_msg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (Nodes[i]->GetLabel() != l) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"invalid " << i << "-th label " << l
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			F[i] = Vec3(&f[0]);
			M[i] = Vec3(&f[3]);
		}
	}
}

SubVectorHandler&
StructExtForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	WorkVec.ResizeReset(6*Nodes.size());

	if (pRefNode) {
		// manipulate

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			integer iFirstIndex = Nodes[i]->iGetFirstMomentumIndex();
			for (int r = 1; r <= 6; r++) {
				WorkVec.PutRowIndex(i*6 + r, iFirstIndex + r);
			}

			WorkVec.Add(i*6 + 1, F[i]);
			WorkVec.Add(i*6 + 4, M[i] + (Nodes[i]->GetRCurr()*Offsets[i]).Cross(F[i]));
		}
	}

	return WorkVec;
}

void
StructExtForce::Output(OutputHandler& OH) const
{
	std::ostream& out = OH.Forces();

	for (unsigned i = 0; i < Nodes.size(); i++) {
		out << GetLabel() << "@" << Nodes[i]->GetLabel()
			<< " " << F[i]
			<< " " << M[i]
			<< std::endl;
	}
}
 
void
StructExtForce::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(Nodes.size());
	for (unsigned int i = 0; i < Nodes.size(); i++) {
		connectedNodes[i] = Nodes[i];
	}
}

Elem*
ReadStructExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH = 0;
	int iCoupling;

	bool bSendAfterPredict;
	ReadExtForce(pDM, HP, uLabel, pEFH, bSendAfterPredict, iCoupling);

	bool bUnsorted(false);
	bool bNoLabels(false);
	if (HP.IsKeyWord("unsorted")) {
		bUnsorted = true;

	} else if (HP.IsKeyWord("no" "labels")) {
		bNoLabels = true;
	}

	bool bOutputAccelerations(false);
	if (HP.IsKeyWord("accelerations")) {
		bOutputAccelerations = true;
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("StructExtForce(" << uLabel << "): illegal node number " << n <<
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<StructNode *> Nodes(n);
	std::vector<Vec3> Offsets(n);

	for (int i = 0; i < n; i++ ) {
		Nodes[i] = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));
		
		ReferenceFrame RF(Nodes[i]);

		if (HP.IsKeyWord("offset")) {
			Offsets[i] = HP.GetPosRel(RF);
		} else {
			Offsets[i] = Vec3(0.);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, StructExtForce,
		StructExtForce(uLabel, pDM, Nodes, Offsets,
			bUnsorted, bNoLabels, bOutputAccelerations,
			pEFH, bSendAfterPredict, iCoupling, fOut));

	return pEl;
}

