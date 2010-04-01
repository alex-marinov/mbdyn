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

#include "dataman.h"
#include "strext.h"
#include "Rot.hh"

#include <fstream>
#include <cerrno>

/* StructExtForce - begin */

/* Costruttore */
StructExtForce::StructExtForce(unsigned int uL,
	DataManager *pDM,
	StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	bool bSorted,
	bool bLabels,
	bool bOutputAccelerations,
	unsigned uRot,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut), 
pRefNode(pRefNode),
bUseReferenceNodeForces(bUseReferenceNodeForces),
bRotateReferenceNodeForces(bRotateReferenceNodeForces),
F0(0.), M0(0.),
bLabels(bLabels),
bSorted(bSorted),
uRot(uRot),
bOutputAccelerations(bOutputAccelerations),
iobuf_labels(0),
iobuf_x(0),
iobuf_R(0),
iobuf_theta(0),
iobuf_euler_123(0),
iobuf_xp(0),
iobuf_omega(0),
iobuf_xpp(0),
iobuf_omegap(0),
iobuf_f(0),
iobuf_m(0)
{
	ASSERT(nodes.size() == offsets.size());
	Nodes.resize(nodes.size());
	Offsets.resize(nodes.size());
	F.resize(nodes.size());
	M.resize(nodes.size());

	switch (uRot) {
	case MBC_ROT_THETA:
	case MBC_ROT_MAT:
	case MBC_ROT_EULER_123:
		break;

	default:
		silent_cerr("StructExtForce(" << GetLabel() << "): "
			"unknown rotation format " << uRot
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (unsigned i = 0; i < nodes.size(); i++) {
		Nodes[i] = nodes[i];
		Offsets[i] = offsets[i];
		F[i] = Zero3;
		M[i] = Zero3;
	}

	ASSERT(!(!bLabels && !bSorted));
	if (!bSorted) {
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

	// I/O will use filedes
	if (!pEFH->GetOutStream()) {
		node_kinematics_size = 3 + 3 + 3;

		switch (uRot) {
		case MBC_ROT_MAT:
			node_kinematics_size += 9;
			break;

		case MBC_ROT_THETA:
		case MBC_ROT_EULER_123:
			node_kinematics_size += 3;
			break;

		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bOutputAccelerations) {
			node_kinematics_size += 3 + 3;
		}

		node_kinematics_size *= sizeof(doublereal);
		dynamics_size = (3 + 3)*sizeof(doublereal);

		if (bLabels) {
			node_kinematics_size += sizeof(uint32_t);
			dynamics_size += sizeof(uint32_t);
		}

		iobuf.resize(node_kinematics_size*Nodes.size());
		dynamics_size *= Nodes.size();

		char *ptr = &iobuf[0];
		if (bLabels) {
			iobuf_labels = (uint32_t *)ptr;
			ptr += sizeof(uint32_t)*Nodes.size();
		}

		iobuf_x = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*Nodes.size();

		switch (uRot) {
		case MBC_ROT_MAT:
			iobuf_R = (doublereal *)ptr;
			ptr += 9*sizeof(doublereal)*Nodes.size();
			break;

		case MBC_ROT_THETA:
			iobuf_theta = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*Nodes.size();
			break;

		case MBC_ROT_EULER_123:
			iobuf_euler_123 = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*Nodes.size();
			break;
		}

		iobuf_xp = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*Nodes.size();

		iobuf_omega = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*Nodes.size();

		if (bOutputAccelerations) {
			iobuf_xpp = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*Nodes.size();

			iobuf_omegap = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*Nodes.size();
		}

		ptr = (char *)iobuf_x;

		iobuf_f = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*Nodes.size();

		iobuf_m = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*Nodes.size();
	}
}

StructExtForce::~StructExtForce(void)
{
	NO_OP;
}

bool
StructExtForce::Prepare(ExtFileHandlerBase *pEFH)
{
	bool bResult = true;

	switch (pEFH->NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_NO:
		break;

	case ExtFileHandlerBase::NEGOTIATE_CLIENT: {
		std::ostream *outfp = pEFH->GetOutStream();
		if (outfp) {

#ifdef USE_SOCKET
		} else {
			char buf[sizeof(uint32_t) + sizeof(uint32_t)];
			uint32_t *uint32_ptr;

			uint32_ptr = (uint32_t *)&buf[0];
			uint32_ptr[0] = MBC_NODAL | uRot;
			if (pRefNode != 0) {
				uint32_ptr[0] |= MBC_REF_NODE;
			}

			if (bLabels) {
				uint32_ptr[0] |= MBC_LABELS;
			}

			if (bOutputAccelerations) {
				uint32_ptr[0] |= MBC_ACCELS;
			}

			uint32_ptr[1] = Nodes.size();

			ssize_t rc = send(pEFH->GetOutFileDes(),
				(const void *)buf, sizeof(buf),
				pEFH->GetSendFlags());
			if (rc == -1) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"negotiation request send() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"negotiation request send() failed "
					"(sent " << rc << " of " << sizeof(buf) << " bytes)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
#endif // USE_SOCKET
		}
		} break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER: {
		unsigned uN;
		unsigned uNodal;
		bool bRef;
		unsigned uR;
		bool bA;
		bool bL;

		std::istream *infp = pEFH->GetInStream();
		if (infp) {

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
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"negotiation response recv() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"negotiation response recv() failed "
					"(got " << rc << " of " << sizeof(buf) << " bytes)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			uint32_ptr = (uint32_t *)&buf[0];
			uNodal = (uint32_ptr[0] & MBC_MODAL_NODAL_MASK);
			bRef = (uint32_ptr[0] & MBC_REF_NODE);
			uR = (uint32_ptr[0] & MBC_ROT_MASK);
			bL = (uint32_ptr[0] & MBC_LABELS);
			bA = (uint32_ptr[0] & MBC_ACCELS);

			uN = uint32_ptr[1];
#endif // USE_SOCKET
		}

		if (uNodal != MBC_NODAL) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: expecting MBC_NODAL "
				"(=" << MBC_MODAL << "), got " << uNodal
				<< std::endl);
			bResult = false;
		}

		if ((pRefNode != 0 && !bRef) || (pRefNode == 0 && bRef)) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: reference node configuration mismatch "
				"(local=" << (pRefNode != 0 ? "yes" : "no") << ", remote=" << (bRef ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (uR != uRot) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: orientation output mismatch "
				"(local=" << uRot  << ", remote=" << uR << ")"
				<< std::endl);
			bResult = false;
		}

		if (bL != bLabels) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: labels output mismatch "
				"(local=" << (bLabels ? "yes" : "no") << ", remote=" << (bL ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (bA != bOutputAccelerations) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: acceleration output mismatch "
				"(local=" << (bOutputAccelerations ? "yes" : "no") << ", remote=" << (bA ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (uN != Nodes.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: node number mismatch "
				"(local=" << Nodes.size() << ", remote=" << uN << ")"
				<< std::endl);
			bResult = false;
		}
		} break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return bResult;
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
		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();
		const Vec3& xpRef = pRefNode->GetVCurr();
		const Vec3& wRef = pRefNode->GetWCurr();
		const Vec3& xppRef = pRefNode->GetXPPCurr();
		const Vec3& wpRef = pRefNode->GetWPCurr();

		if (bLabels) {
			outf
				<< pRefNode->GetLabel()
				<< " ";
		}
		outf
			<< xRef
			<< " ";
		switch (uRot) {
		case MBC_ROT_MAT:
			outf << RRef;
			break;

		case MBC_ROT_THETA:
			outf << RotManip::VecRot(RRef);
			break;

		case MBC_ROT_EULER_123:
			outf << MatR2EulerAngles123(RRef)*dRaDegr;
			break;
		}
		outf
			<< " " << xpRef
			<< " " << wRef;
		if (bOutputAccelerations) {
			outf
				<< " " << xppRef
				<< " " << wpRef;
		}
		outf << std::endl;

		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 f(Nodes[i]->GetRCurr()*Offsets[i]);
			Vec3 x(Nodes[i]->GetXCurr() + f);
			Vec3 Dx(x - xRef);
			Mat3x3 DR(RRef.MulTM(Nodes[i]->GetRCurr()));
			Vec3 v(Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(f));
			Vec3 Dv(v - xpRef - wRef.Cross(Dx));
			const Vec3& w(Nodes[i]->GetWCurr());

			// manipulate

			if (bLabels) {
				outf
					<< Nodes[i]->GetLabel()
					<< " ";
			}
			outf
				<< RRef.MulTV(Dx)
				<< " ";
			switch (uRot) {
			case MBC_ROT_MAT:
				outf << DR;
				break;

			case MBC_ROT_THETA:
				outf << RotManip::VecRot(DR);
				break;

			case MBC_ROT_EULER_123:
				outf << MatR2EulerAngles123(DR)*dRaDegr;
				break;
			}
			outf
				<< " " << RRef.MulTV(Dv)
				<< " " << RRef.MulTV(w - wRef);
			if (bOutputAccelerations) {
				const Vec3& xpp(Nodes[i]->GetXPPCurr());
				const Vec3& wp(Nodes[i]->GetWPCurr());

				outf
					<< " " << RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2))
					<< " " << RRef.MulTV(wp - wpRef - wRef.Cross(w));
			}
			outf << std::endl;
		}

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

			if (bLabels) {
				outf
					<< Nodes[i]->GetLabel()
					<< " ";
			}
			outf
				<< x
				<< " ";
			switch (uRot) {
			case MBC_ROT_MAT:
				outf << R;
				break;

			case MBC_ROT_THETA:
				outf << RotManip::VecRot(R);
				break;

			case MBC_ROT_EULER_123:
				outf << MatR2EulerAngles123(R)*dRaDegr;
				break;
			}
			outf
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
#ifdef USE_SOCKET
	if (pRefNode) {
		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();
		const Vec3& xpRef = pRefNode->GetVCurr();
		const Vec3& wRef = pRefNode->GetWCurr();
		const Vec3& xppRef = pRefNode->GetXPPCurr();
		const Vec3& wpRef = pRefNode->GetWPCurr();

		if (bLabels) {
			uint32_t l = pRefNode->GetLabel();
			send(outfd, (void *)&l, sizeof(l), 0);
		}

		send(outfd, (void *)xRef.pGetVec(), 3*sizeof(doublereal), 0);
		switch (uRot) {
		case MBC_ROT_MAT:
			send(outfd, (void *)RRef.pGetMat(), 9*sizeof(doublereal), 0);
			break;

		case MBC_ROT_THETA: {
			Vec3 Theta(RotManip::VecRot(RRef));
			send(outfd, (void *)Theta.pGetVec(), 3*sizeof(doublereal), 0);
			} break;

		case MBC_ROT_EULER_123: {
			Vec3 E(MatR2EulerAngles123(RRef)*dRaDegr);
			send(outfd, (void *)E.pGetVec(), 3*sizeof(doublereal), 0);
			} break;
		}
		send(outfd, (void *)xpRef.pGetVec(), 3*sizeof(doublereal), 0);
		send(outfd, (void *)wRef.pGetVec(), 3*sizeof(doublereal), 0);
		if (bOutputAccelerations) {
			send(outfd, (void *)xppRef.pGetVec(), 3*sizeof(doublereal), 0);
			send(outfd, (void *)wpRef.pGetVec(), 3*sizeof(doublereal), 0);
		}

		for (unsigned i = 0; i < Nodes.size(); i++) {
			Vec3 f(Nodes[i]->GetRCurr()*Offsets[i]);
			Vec3 x(Nodes[i]->GetXCurr() + f);
			Vec3 Dx(x - xRef);
			Mat3x3 DR(RRef.MulTM(Nodes[i]->GetRCurr()));
			Vec3 v(Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(f));
			Vec3 Dv(v - xpRef - wRef.Cross(Dx));
			const Vec3& w(Nodes[i]->GetWCurr());

			// manipulate
			if (bLabels) {
				uint32_t l = Nodes[i]->GetLabel();
				iobuf_labels[i] = l;
			}

			Vec3 xTilde(RRef.MulTV(Dx));
			memcpy(&iobuf_x[3*i], xTilde.pGetVec(), 3*sizeof(doublereal));

			switch (uRot) {
			case MBC_ROT_MAT:
				memcpy(&iobuf_R[9*i], DR.pGetMat(), 9*sizeof(doublereal));
				break;

			case MBC_ROT_THETA: {
				Vec3 Theta(RotManip::VecRot(DR));
				memcpy(&iobuf_theta[3*i], Theta.pGetVec(), 3*sizeof(doublereal));
				} break;

			case MBC_ROT_EULER_123: {
				Vec3 E(MatR2EulerAngles123(DR)*dRaDegr);
				memcpy(&iobuf_euler_123[3*i], E.pGetVec(), 3*sizeof(doublereal));
				} break;
			}

			Vec3 vTilde(RRef.MulTV(Dv));
			memcpy(&iobuf_xp[3*i], vTilde.pGetVec(), 3*sizeof(doublereal));

			Vec3 wTilde(RRef.MulTV(w - wRef));
			memcpy(&iobuf_omega[3*i], wTilde.pGetVec(), 3*sizeof(doublereal));

			if (bOutputAccelerations) {
				const Vec3& xpp = Nodes[i]->GetXPPCurr();
				Vec3 xppTilde(RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
					- wRef.Cross(wRef.Cross(Dx) + Dv*2)));
				memcpy(&iobuf_xpp[3*i], xppTilde.pGetVec(), 3*sizeof(doublereal));

				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 wpTilde(RRef.MulTV(wp) - wpRef - wRef.Cross(w));
				memcpy(&iobuf_omegap[3*i], wpTilde.pGetVec(), 3*sizeof(doublereal));
			}
		}

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

			if (bLabels) {
				uint32_t l = Nodes[i]->GetLabel();
				iobuf_labels[i] = l;
			}
			memcpy(&iobuf_x[3*i], x.pGetVec(), 3*sizeof(doublereal));
			switch (uRot) {
			case MBC_ROT_MAT:
				memcpy(&iobuf_R[9*i], R.pGetMat(), 9*sizeof(doublereal));
				break;

			case MBC_ROT_THETA: {
				Vec3 Theta(RotManip::VecRot(R));
				memcpy(&iobuf_theta[3*i], Theta.pGetVec(), 3*sizeof(doublereal));
				} break;

			case MBC_ROT_EULER_123: {
				Vec3 E(MatR2EulerAngles123(R)*dRaDegr);
				memcpy(&iobuf_euler_123[3*i], E.pGetVec(), 3*sizeof(doublereal));
				} break;
			}
			memcpy(&iobuf_xp[3*i], v.pGetVec(), 3*sizeof(doublereal));
			memcpy(&iobuf_omega[3*i], w.pGetVec(), 3*sizeof(doublereal));

			if (bOutputAccelerations) {
				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 xpp = Nodes[i]->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				memcpy(&iobuf_xpp[3*i], xpp.pGetVec(), 3*sizeof(doublereal));
				memcpy(&iobuf_omegap[3*i], wp.pGetVec(), 3*sizeof(doublereal));
			}
		}
	}

	send(outfd, &iobuf[0], iobuf.size(), 0);
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
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
	if (pRefNode) {
		unsigned l;
		doublereal *f = F0.pGetVec(), *m = M0.pGetVec();

		inf >> l
			>> f[0] >> f[1] >> f[2]
			>> m[0] >> m[1] >> m[2];
	}

	if (!bSorted) {
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

			if (bLabels) {
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
#ifdef USE_SOCKET
	if (pRefNode) {
		unsigned l;
		size_t ulen;
		char buf[sizeof(uint32_t) + 6*sizeof(doublereal)];
		doublereal *f;
		ssize_t len;

		if (bLabels) {
			ulen = sizeof(buf);

		} else {
			ulen = 6*sizeof(doublereal);
		}

		len = recv(infd, (void *)buf, ulen, 0);
		if (len == -1) {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"recv() failed (" << save_errno << ": "
				<< err_msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (unsigned(len) != ulen) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"recv() failed " "(got " << len << " of "
				<< ulen << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bLabels) {
			uint32_t *uint32_ptr = (uint32_t *)buf;
			l = uint32_ptr[0];
			f = (doublereal *)&uint32_ptr[1];

		} else {
			f = (doublereal *)buf;
		}

		F0 = Vec3(&f[0]);
		M0 = Vec3(&f[3]);
	}

	ssize_t len = recv(infd, (void *)&iobuf[0], dynamics_size, 0);
	if (len == -1) {
		int save_errno = errno;
		char *err_msg = strerror(save_errno);
		silent_cerr("StructExtForce(" << GetLabel() << "): "
			"recv() failed (" << save_errno << ": "
			<< err_msg << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	} else if (unsigned(len) != dynamics_size) {
		silent_cerr("StructExtForce(" << GetLabel() << "): "
			"recv() failed " "(got " << len << " of "
			<< dynamics_size << " bytes)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!bSorted) {
		if (bLabels) {
			done.resize(Nodes.size());

			fill(done.begin(), done.end(), false);
		}

		unsigned cnt;
		for (cnt = 0; cnt < Nodes.size(); cnt++) {
			unsigned i = cnt;
			if (bLabels) {
				unsigned l = iobuf_labels[cnt];
				std::vector<StructNode *>::const_iterator n;
				for (n = Nodes.begin(); n != Nodes.end(); n++) {
					if ((*n)->GetLabel() == l) {
						break;
					}
				}

				if (n == Nodes.end()) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"unknown label " << l
						<< " as " << cnt << "-th node"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				i = n - Nodes.begin();

				if (done[i]) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"label " << l << " already done"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				done[i] = true;
			}

			F[i] = &iobuf_f[3*i];
			M[i] = &iobuf_m[3*i];
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
			if (bLabels) {
				unsigned l = iobuf_labels[i];
				if (Nodes[i]->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			F[i] = Vec3(&iobuf_f[3*i]);
			M[i] = Vec3(&iobuf_m[3*i]);
		}
	}
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

SubVectorHandler&
StructExtForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	if (pRefNode) {
		integer iSize = Nodes.size();
		if (bUseReferenceNodeForces) {
			iSize++;
		}
		WorkVec.ResizeReset(6*iSize);

		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();

		Vec3 FTmp, MTmp;

		// manipulate
		if (bUseReferenceNodeForces) {
			if (bRotateReferenceNodeForces) {
				FTmp = RRef*F0;
				MTmp = RRef*M0;

			} else {
				FTmp = F0;
				MTmp = M0;
			}
		}

		for (unsigned i = 0; i < Nodes.size(); i++) {
			integer iFirstIndex = Nodes[i]->iGetFirstMomentumIndex();
			for (int r = 1; r <= 6; r++) {
				WorkVec.PutRowIndex(i*6 + r, iFirstIndex + r);
			}

			Vec3 f(RRef*F[i]);
			Vec3 m(RRef*M[i] + (Nodes[i]->GetRCurr()*Offsets[i]).Cross(f));

			WorkVec.Add(i*6 + 1, f);
			WorkVec.Add(i*6 + 4, m);

			if (bUseReferenceNodeForces) {
				FTmp -= f;
				MTmp -= m + (Nodes[i]->GetXCurr() - xRef).Cross(f);
			}
		}

		if (bUseReferenceNodeForces) {
			unsigned i = Nodes.size();
			integer iFirstIndex = pRefNode->iGetFirstMomentumIndex();
			for (int r = 1; r <= 6; r++) {
				WorkVec.PutRowIndex(i*6 + r, iFirstIndex + r);
			}

			WorkVec.Add(i*6 + 1, FTmp);
			WorkVec.Add(i*6 + 4, MTmp);
		}

	} else {
		WorkVec.ResizeReset(6*Nodes.size());

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

	StructNode *pRefNode(0);
	if (HP.IsKeyWord("reference" "node")) {
		pRefNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (pRefNode == 0) {
			silent_cerr("StructExtForce(" << uLabel << "): "
				"illegal reference node "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool bSorted(false);
	bool bLabels(false);
	unsigned uRot = MBC_ROT_MAT;
	bool bOutputAccelerations(false);
	bool bUseReferenceNodeForces(true);
	bool bRotateReferenceNodeForces(true);

	bool bGotSorted(false);
	bool bGotLabels(false);
	bool bGotRot(false);
	bool bGotAccels(false);
	bool bGotUseRefForces(false);

	while (HP.IsArg()) {
		if (HP.IsKeyWord("unsorted")) {
			silent_cerr("StructExtForce(" << uLabel << "): "
				"use of \"unsorted\" deprecated in favor of \"sorted, { yes | no }\" at line "
				<< HP.GetLineData() << std::endl);

			if (bGotSorted) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"unsorted\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bSorted = false;
			bGotSorted = true;

		} else if (HP.IsKeyWord("sorted")) {
			if (bGotSorted) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"sorted\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bSorted = HP.GetYesNo(bSorted);
			bGotSorted = true;

		} else if (HP.IsKeyWord("no" "labels")) {
			silent_cerr("StructExtForce(" << uLabel << "): "
				"use of \"no labels\" deprecated in favor of \"labels, { yes | no }\" at line "
				<< HP.GetLineData() << std::endl);

			if (bGotLabels) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"no labels\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bLabels = false;
			bGotLabels = true;

		} else if (HP.IsKeyWord("labels")) {
			if (bGotLabels) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"labels\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bLabels = HP.GetYesNo(bLabels);
			bGotLabels = true;

		} else if (HP.IsKeyWord("orientation")) {
			if (bGotRot) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"orientation\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (HP.IsKeyWord("orientation" "vector")) {
				uRot = MBC_ROT_THETA;

			} else if (HP.IsKeyWord("orientation" "matrix")) {
				uRot = MBC_ROT_MAT;

			} else if (HP.IsKeyWord("euler" "123")) {
				uRot = MBC_ROT_EULER_123;

			} else {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"unknown \"orientation\" format at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bGotRot = true;	

		} else if (HP.IsKeyWord("accelerations")) {
			if (bGotAccels) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"accelerations\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bOutputAccelerations = HP.GetYesNo(bOutputAccelerations);
			bGotAccels = true;

		} else if (HP.IsKeyWord("use" "reference" "node" "forces")) {
			if (pRefNode == 0) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"use reference node forces\" only meaningful when reference node is used at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (bGotUseRefForces) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"use reference node forces\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			
			bUseReferenceNodeForces = HP.GetYesNo(bUseReferenceNodeForces);
			bGotUseRefForces = true;

			if (bUseReferenceNodeForces && HP.IsKeyWord("rotate" "reference" "node" "forces")) {
				bRotateReferenceNodeForces = HP.GetYesNo(bRotateReferenceNodeForces);
			}

		} else {
			break;
		}
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
		StructExtForce(uLabel, pDM, pRefNode,
			bUseReferenceNodeForces, bRotateReferenceNodeForces,
			Nodes, Offsets,
			bSorted, bLabels, bOutputAccelerations, uRot,
			pEFH, bSendAfterPredict, iCoupling, fOut));

	return pEl;
}

