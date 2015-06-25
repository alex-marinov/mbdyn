/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2015
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
#include "strext.h"
#include "stredge.h"
#include "Rot.hh"

#include <fstream>
#include <cerrno>
#include <algorithm>

/* StructExtForce - begin */

/* Costruttore */
StructExtForce::StructExtForce(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<unsigned >& labels,
	std::vector<const StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	bool bSorted,
	bool bLabels,
	bool bOutputAccelerations,
	unsigned uRot,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	unsigned uOutputFlags,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut), 
pRefNode(pRefNode),
bUseReferenceNodeForces(bUseReferenceNodeForces),
bRotateReferenceNodeForces(bRotateReferenceNodeForces),
F0(Zero3), M0(Zero3),
F1(Zero3), M1(Zero3),
F2(Zero3), M2(Zero3),
uOutputFlags(uOutputFlags),
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
	ASSERT((!bLabels) || (nodes.size() == offsets.size()));

	m_Points.resize(nodes.size());

	switch (uRot) {
	case MBC_ROT_NONE:
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
		m_Points[i].pNode = nodes[i];
		m_Points[i].Offset = offsets[i];
		m_Points[i].F = Zero3;
		m_Points[i].M = Zero3;
		if (bLabels) {
			m_Points[i].uLabel = labels[i];
		} else {
			m_Points[i].uLabel = m_Points[i].pNode->GetLabel();
		}
	}

	ASSERT(!(!bLabels && !bSorted));
	if (!bSorted) {
		done.resize(nodes.size());
	}

	if (bOutputAccelerations) {
		for (unsigned i = 0; i < m_Points.size(); i++) {
			const DynamicStructNode *pDSN = dynamic_cast<const DynamicStructNode *>(m_Points[i].pNode);
			if (pDSN == 0) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"point #" << i << " StructNode(" << m_Points[i].pNode->GetLabel() << ") "
					"is not dynamic"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			const_cast<DynamicStructNode *>(pDSN)->ComputeAccelerations(true);
		}
	}

	unsigned node_kinematics_size = 0;
	unsigned labels_size = 0;

	// I/O will use filedes
	if (!pEFH->GetOutStream()) {
		node_kinematics_size = 3 + 3;

		switch (uRot) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			node_kinematics_size += 9 + 3;
			break;

		case MBC_ROT_THETA:
		case MBC_ROT_EULER_123:
			node_kinematics_size += 3 + 3;
			break;

		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bOutputAccelerations) {
			node_kinematics_size += 3;
			if (uRot != MBC_ROT_NONE) {
				node_kinematics_size += 3;
			}
		}

		node_kinematics_size *= sizeof(doublereal);

		dynamics_size = 3;
		if (uRot != MBC_ROT_NONE) {
			dynamics_size += 3;
		}

		dynamics_size *= sizeof(doublereal);

		if (bLabels) {
			labels_size = sizeof(uint32_t)*(m_Points.size() + m_Points.size()%2);
		}

		iobuf.resize(node_kinematics_size*m_Points.size() + labels_size);
		dynamics_size = dynamics_size*m_Points.size() + labels_size;

		char *ptr = &iobuf[0];
		if (bLabels) {
			iobuf_labels = (uint32_t *)ptr;
			ptr += labels_size;
		}

		iobuf_x = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*m_Points.size();

		switch (uRot) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			iobuf_R = (doublereal *)ptr;
			ptr += 9*sizeof(doublereal)*m_Points.size();
			break;

		case MBC_ROT_THETA:
			iobuf_theta = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*m_Points.size();
			break;

		case MBC_ROT_EULER_123:
			iobuf_euler_123 = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*m_Points.size();
			break;
		}

		iobuf_xp = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*m_Points.size();

		if (uRot != MBC_ROT_NONE) {
			iobuf_omega = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*m_Points.size();
		}

		if (bOutputAccelerations) {
			iobuf_xpp = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*m_Points.size();

			if (uRot != MBC_ROT_NONE) {
				iobuf_omegap = (doublereal *)ptr;
				ptr += 3*sizeof(doublereal)*m_Points.size();
			}
		}

		ptr = (char *)iobuf_x;

		iobuf_f = (doublereal *)ptr;
		ptr += 3*sizeof(doublereal)*m_Points.size();

		if (uRot != MBC_ROT_NONE) {
			iobuf_m = (doublereal *)ptr;
			ptr += 3*sizeof(doublereal)*m_Points.size();
		}
	}
}

StructExtForce::~StructExtForce(void)
{
	NO_OP;
}

void
StructExtForce::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	if (iCoupling == COUPLING_NONE) {
		*piNumRows = 0;
		*piNumCols = 0;

	} else {
		*piNumRows = (pRefNode ? 6 : 0) + 6*m_Points.size();
		*piNumCols = 1;
	}
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

			uint32_ptr[1] = m_Points.size();

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
			// TODO: stream negotiation?

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

		if (uN != m_Points.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"negotiation response failed: node number mismatch "
				"(local=" << m_Points.size() << ", remote=" << uN << ")"
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

		outf << "# reference node" << std::endl;

		if (bLabels) {
			outf
				<< pRefNode->GetLabel()
				<< " ";
		}

		outf
			<< xRef;

		switch (uRot) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			outf
				<< " " << RRef;
			break;

		case MBC_ROT_THETA:
			outf
				<< " " << RotManip::VecRot(RRef);
			break;

		case MBC_ROT_EULER_123:
			outf
				<< " " << MatR2EulerAngles123(RRef)*dRaDegr;
			break;
		}

		outf
			<< " " << xpRef;

		if (uRot != MBC_ROT_NONE) {
			outf
				<< " " << wRef;
		}

		if (bOutputAccelerations) {
			outf
				<< " " << xppRef;

			if (uRot != MBC_ROT_NONE) {
				outf
					<< " " << wpRef;
			}
		}
		outf << std::endl;

		outf << "# regular nodes" << std::endl;

		for (std::vector<PointData>::const_iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
			Vec3 f(point->pNode->GetRCurr()*point->Offset);
			Vec3 x(point->pNode->GetXCurr() + f);
			Vec3 Dx(x - xRef);
			Mat3x3 DR(RRef.MulTM(point->pNode->GetRCurr()));
			Vec3 v(point->pNode->GetVCurr() + point->pNode->GetWCurr().Cross(f));
			Vec3 Dv(v - xpRef - wRef.Cross(Dx));
			const Vec3& w(point->pNode->GetWCurr());

			// manipulate

			if (bLabels) {
				outf
					<< point->pNode->GetLabel()
					<< " ";
			}

			outf
				<< RRef.MulTV(Dx);

			switch (uRot) {
			case MBC_ROT_NONE:
				break;

			case MBC_ROT_MAT:
				outf
					<< " " << DR;
				break;

			case MBC_ROT_THETA:
				outf
					<< " " << RotManip::VecRot(DR);
				break;

			case MBC_ROT_EULER_123:
				outf
					<< " " << MatR2EulerAngles123(DR)*dRaDegr;
				break;
			}

			outf
				<< " " << RRef.MulTV(Dv);
			
			if (uRot != MBC_ROT_NONE) {
				outf
					<< " " << RRef.MulTV(w - wRef);
			}

			if (bOutputAccelerations) {
				const Vec3& xpp(point->pNode->GetXPPCurr());
				const Vec3& wp(point->pNode->GetWPCurr());

				outf
					<< " " << RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2)
							+ wp.Cross(f) + w.Cross(w.Cross(f)));
				if (uRot != MBC_ROT_NONE) {
					outf
						<< " " << RRef.MulTV(wp - wpRef - wRef.Cross(w));
				}
			}
			outf << std::endl;
		}

	} else {
		outf << "# regular nodes" << std::endl;

		for (std::vector<PointData>::const_iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
			/*
				p = x + f
				R = R
				v = xp + w cross f
				w = w
				a = xpp + wp cross f + w cross w cross f
				wp = wp
			 */

			// Optimization of the above formulas
			const Mat3x3& R = point->pNode->GetRCurr();
			Vec3 f = R*point->Offset;
			Vec3 x = point->pNode->GetXCurr() + f;
			const Vec3& w = point->pNode->GetWCurr();
			Vec3 wCrossf = w.Cross(f);
			Vec3 v = point->pNode->GetVCurr() + wCrossf;

			if (bLabels) {
				outf
					<< point->pNode->GetLabel()
					<< " ";
			}

			outf
				<< x;

			switch (uRot) {
			case MBC_ROT_NONE:
				break;

			case MBC_ROT_MAT:
				outf
					<< " " << R;
				break;

			case MBC_ROT_THETA:
				outf
					<< " " << RotManip::VecRot(R);
				break;

			case MBC_ROT_EULER_123:
				outf
					<< " " << MatR2EulerAngles123(R)*dRaDegr;
				break;
			}
			outf
				<< " " << v;

			if (uRot != MBC_ROT_NONE) {
				outf
					<< " " << w;
			}

			if (bOutputAccelerations) {
				const Vec3& wp = point->pNode->GetWPCurr();
				Vec3 a = point->pNode->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				outf
					<< " " << a;

				if (uRot != MBC_ROT_NONE) {
					outf
						<< " " << wp;
				}
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
			uint32_t l[2];
			l[0] = pRefNode->GetLabel();
			l[1] = 0;
         		send(outfd, (void *)&l[0], sizeof(l), 0);
		}

		send(outfd, (void *)xRef.pGetVec(), 3*sizeof(doublereal), 0);
		switch (uRot) {
		case MBC_ROT_NONE:
			break;

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
		if (uRot != MBC_ROT_NONE) {
			send(outfd, (void *)wRef.pGetVec(), 3*sizeof(doublereal), 0);
		}
		if (bOutputAccelerations) {
			send(outfd, (void *)xppRef.pGetVec(), 3*sizeof(doublereal), 0);
			if (uRot != MBC_ROT_NONE) {
				send(outfd, (void *)wpRef.pGetVec(), 3*sizeof(doublereal), 0);
			}
		}

		for (unsigned i = 0; i < m_Points.size(); i++) {
			const PointData& point = m_Points[i];
			Vec3 f(point.pNode->GetRCurr()*point.Offset);
			Vec3 x(point.pNode->GetXCurr() + f);
			Vec3 Dx(x - xRef);
			Mat3x3 DR(RRef.MulTM(point.pNode->GetRCurr()));
			Vec3 v(point.pNode->GetVCurr() + point.pNode->GetWCurr().Cross(f));
			Vec3 Dv(v - xpRef - wRef.Cross(Dx));
			const Vec3& w(point.pNode->GetWCurr());

			// manipulate
			if (bLabels) {
				uint32_t l = point.pNode->GetLabel();
				iobuf_labels[i] = l;
			}

			Vec3 xTilde(RRef.MulTV(Dx));
			memcpy(&iobuf_x[3*i], xTilde.pGetVec(), 3*sizeof(doublereal));

			switch (uRot) {
			case MBC_ROT_NONE:
				break;

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

			if (uRot != MBC_ROT_NONE) {
				Vec3 wTilde(RRef.MulTV(w - wRef));
				memcpy(&iobuf_omega[3*i], wTilde.pGetVec(), 3*sizeof(doublereal));
			}

			if (bOutputAccelerations) {
				const Vec3& xpp = point.pNode->GetXPPCurr();
				const Vec3& wp = point.pNode->GetWPCurr();

				Vec3 xppTilde(RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
					- wRef.Cross(wRef.Cross(Dx) + Dv*2)
					+ wp.Cross(f) + w.Cross(w.Cross(f))));
				memcpy(&iobuf_xpp[3*i], xppTilde.pGetVec(), 3*sizeof(doublereal));

				if (uRot != MBC_ROT_NONE) {
					Vec3 wpTilde(RRef.MulTV(wp) - wpRef - wRef.Cross(w));
					memcpy(&iobuf_omegap[3*i], wpTilde.pGetVec(), 3*sizeof(doublereal));
				}
			}
		}

	} else {
		for (unsigned i = 0; i < m_Points.size(); i++) {
			/*
				p = x + f
				R = R
				v = xp + w cross f
				w = w
				a = xpp + wp cross f + w cross w cross f
				wp = wp
			 */

			// Optimization of the above formulas
			const PointData& point = m_Points[i];
			const Mat3x3& R = point.pNode->GetRCurr();
			Vec3 f = R*point.Offset;
			Vec3 x = point.pNode->GetXCurr() + f;
			const Vec3& w = point.pNode->GetWCurr();
			Vec3 wCrossf = w.Cross(f);
			Vec3 v = point.pNode->GetVCurr() + wCrossf;

			if (bLabels) {
				uint32_t l = point.pNode->GetLabel();
				iobuf_labels[i] = l;
			}
			memcpy(&iobuf_x[3*i], x.pGetVec(), 3*sizeof(doublereal));
			switch (uRot) {
			case MBC_ROT_NONE:
				break;

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
			if (uRot != MBC_ROT_NONE) {
				memcpy(&iobuf_omega[3*i], w.pGetVec(), 3*sizeof(doublereal));
			}

			if (bOutputAccelerations) {
				const Vec3& wp = point.pNode->GetWPCurr();
				Vec3 xpp = point.pNode->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				memcpy(&iobuf_xpp[3*i], xpp.pGetVec(), 3*sizeof(doublereal));
				if (uRot != MBC_ROT_NONE) {
					memcpy(&iobuf_omegap[3*i], wp.pGetVec(), 3*sizeof(doublereal));
				}
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
		doublereal *f = F0.pGetVec(), *m = M0.pGetVec();

		if (bLabels) {
			// TODO: implement "unsorted"
			unsigned l;
			inf >> l;
		}

		inf >> f[0] >> f[1] >> f[2];
		if (uRot != MBC_ROT_NONE) {
			inf >> m[0] >> m[1] >> m[2];
		}
	}

	if (!bSorted) {
		ASSERT(bLabels);

		done.resize(m_Points.size());
		fill(done.begin(), done.end(), false);

		unsigned cnt;
		for (cnt = 0; inf; cnt++) {
			/* assume unsigned int label */
			unsigned l, i;
			doublereal f[3], m[3];

			inf >> l
				>> f[0] >> f[1] >> f[2];
			if (uRot != MBC_ROT_NONE) {
				inf >> m[0] >> m[1] >> m[2];
			}

			if (!inf) {
				break;
			}

			for (i = 0; i < m_Points.size(); i++) {
				if (m_Points[i].uLabel == l) {
					break;
				}
			}

			if (i == m_Points.size()) {
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

			m_Points[i].F = Vec3(f);
			if (uRot != MBC_ROT_NONE) {
				m_Points[i].M = Vec3(m);
			}
		}

		if (cnt != m_Points.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"invalid number of nodes " << cnt
				<< std::endl);

			for (unsigned i = 0; i < m_Points.size(); i++) {
				if (!done[i]) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"label " << m_Points[i].uLabel << " node " << m_Points[i].pNode->GetLabel()
						<< " not done" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		for (unsigned i = 0; i < m_Points.size(); i++) {
			PointData& point = m_Points[i];


			/* assume unsigned int label */
			doublereal f[3], m[3];

			if (bLabels) {
				unsigned l;
				inf >> l;

				if (point.pNode->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			inf
				>> f[0] >> f[1] >> f[2];
			if (uRot != MBC_ROT_NONE) {
				inf >> m[0] >> m[1] >> m[2];
			}

			if (!inf) {
				break;
			}

			point.F = Vec3(f);
			if (uRot != MBC_ROT_NONE) {
				point.M = Vec3(m);
			}
		}
	}
}

void
StructExtForce::RecvFromFileDes(int infd)
{
#ifdef USE_SOCKET
	if (pRefNode) {
		size_t ulen = 0;
		char buf[2*sizeof(uint32_t) + 6*sizeof(doublereal)];
		doublereal *f;
		ssize_t len;

		if (bLabels) {
			// to align with double
			ulen = 2*sizeof(uint32_t);
		}

		if (uRot != MBC_ROT_NONE) {
			ulen += 6*sizeof(doublereal);

		} else {
			ulen += 3*sizeof(doublereal);
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
				"recv() failed (got " << len << " of "
				<< ulen << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bLabels) {
			uint32_t *uint32_ptr = (uint32_t *)buf;
			unsigned l = uint32_ptr[0];
			if (l != pRefNode->GetLabel()) {
				silent_cerr("StructExtForce(" << GetLabel() << "): "
					"invalid reference node label "
					"(wanted " << pRefNode->GetLabel() << ", got " << l << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			f = (doublereal *)&uint32_ptr[2];

		} else {
			f = (doublereal *)buf;
		}

		F0 = Vec3(&f[0]);
		if (uRot != MBC_ROT_NONE) {
			M0 = Vec3(&f[3]);
		}
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
		ASSERT(bLabels);

		done.resize(m_Points.size());
		fill(done.begin(), done.end(), false);

		unsigned cnt;
		for (cnt = 0; cnt < m_Points.size(); cnt++) {
			PointData& point = m_Points[cnt];

			unsigned l = iobuf_labels[cnt];
			std::vector<PointData>::const_iterator p;
			for (p = m_Points.begin(); p != m_Points.end(); ++p) {
				if (p->uLabel == l) {
					break;
				}
			}

			if (p == m_Points.end()) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"unknown label " << l
					<< " as " << cnt << "-th node"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			unsigned i = p - m_Points.begin();

			if (done[i]) {
				silent_cerr("StructExtForce"
					"(" << GetLabel() << "): "
					"label " << l << " already done"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			done[i] = true;

			point.F = &iobuf_f[3*i];
			if (uRot != MBC_ROT_NONE) {
				point.M = &iobuf_m[3*i];
			}
		}

		if (cnt != m_Points.size()) {
			silent_cerr("StructExtForce(" << GetLabel() << "): "
				"invalid node number " << cnt
				<< std::endl);

			for (unsigned i = 0; i < m_Points.size(); i++) {
				if (!done[i]) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"label " << m_Points[i].uLabel << "node " << m_Points[i].pNode->GetLabel()
						<< " not done" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		for (unsigned i = 0; i < m_Points.size(); i++) {
			PointData& point = m_Points[i];

			if (bLabels) {
				unsigned l = iobuf_labels[i];
				if (point.pNode->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			point.F = Vec3(&iobuf_f[3*i]);
			if (uRot != MBC_ROT_NONE) {
				point.M = Vec3(&iobuf_m[3*i]);
			}
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

	if (iCoupling == COUPLING_NONE) {
		WorkVec.Resize(0);
		return WorkVec;
	}

	int iOffset;
	if (uRot != MBC_ROT_NONE) {
		iOffset = 6;

	} else {
		iOffset = 3;
	}

	if (pRefNode) {
		integer iSize = m_Points.size();
		if (bUseReferenceNodeForces) {
			iSize++;
		}

		WorkVec.ResizeReset(iOffset*iSize);

		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();

		// manipulate
		if (bUseReferenceNodeForces) {
			if (bRotateReferenceNodeForces) {
				F1 = RRef*F0;
				if (uRot != MBC_ROT_NONE) {
					M1 = RRef*M0;
				}

			} else {
				F1 = F0;
				if (uRot != MBC_ROT_NONE) {
					M1 = M0;
				}
			}
		}

		F2 = Zero3;
		M2 = Zero3;
		for (unsigned i = 0; i < m_Points.size(); i++) {
			const PointData& point = m_Points[i];

			integer iFirstIndex = point.pNode->iGetFirstMomentumIndex();
			for (int r = 1; r <= iOffset; r++) {
				WorkVec.PutRowIndex(i*iOffset + r, iFirstIndex + r);
			}

			Vec3 f(RRef*point.F);
			WorkVec.Add(i*iOffset + 1, f);

			Vec3 m;
			if (uRot != MBC_ROT_NONE) {
				m = RRef*point.M + (point.pNode->GetRCurr()*point.Offset).Cross(f);
				WorkVec.Add(i*iOffset + 4, m);
			}

			if (bUseReferenceNodeForces || fToBeOutput()) {
				F2 += f;
				if (uRot != MBC_ROT_NONE) {
					M2 += m + (point.pNode->GetXCurr() - xRef).Cross(f);
				}
			}
		}

		if (bUseReferenceNodeForces) {
			unsigned i = m_Points.size();
			integer iFirstIndex = pRefNode->iGetFirstMomentumIndex();
			for (int r = 1; r <= iOffset; r++) {
				WorkVec.PutRowIndex(i*iOffset + r, iFirstIndex + r);
			}

			F1 -= F2;
			M1 -= M2;

			WorkVec.Add(i*iOffset + 1, F1);
			if (uRot != MBC_ROT_NONE) {
				WorkVec.Add(i*iOffset + 4, M1);
			}
		}

	} else {
		WorkVec.ResizeReset(iOffset*m_Points.size());

		for (unsigned i = 0; i < m_Points.size(); i++) {
			const PointData& point = m_Points[i];

			integer iFirstIndex = point.pNode->iGetFirstMomentumIndex();
			for (int r = 1; r <= iOffset; r++) {
				WorkVec.PutRowIndex(i*iOffset + r, iFirstIndex + r);
			}

			WorkVec.Add(i*iOffset + 1, point.F);
			if (uRot != MBC_ROT_NONE) {
				WorkVec.Add(i*iOffset + 4, point.M + (point.pNode->GetRCurr()*point.Offset).Cross(point.F));
			}
		}
	}

	return WorkVec;
}

void
StructExtForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::FORCES)) {
			std::ostream& out = OH.Forces();

			if (pRefNode && (uOutputFlags & ExtFileHandlerBase::OUTPUT_REF_DYN)) {
				out << GetLabel() << "#" << pRefNode->GetLabel()
					<< " " << F0
					<< " " << M0
					<< " " << F1
					<< " " << M1
					<< " " << F2
					<< " " << M2
					<< std::endl;
			}

			if (uOutputFlags & ExtFileHandlerBase::OUTPUT_DYN) {
				for (std::vector<PointData>::const_iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
					out << GetLabel() << "@" << point->uLabel
						<< " " << point->F
						<< " " << point->M
						<< std::endl;
				}
			}

			if (uOutputFlags & ExtFileHandlerBase::OUTPUT_KIN) {
				if (pRefNode) {
					const Vec3& xRef = pRefNode->GetXCurr();
					const Mat3x3& RRef = pRefNode->GetRCurr();
					const Vec3& xpRef = pRefNode->GetVCurr();
					const Vec3& wRef = pRefNode->GetWCurr();
					//const Vec3& xppRef = pRefNode->GetXPPCurr();
					//const Vec3& wpRef = pRefNode->GetWPCurr();

					for (std::vector<PointData>::const_iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
						Vec3 f(point->pNode->GetRCurr()*point->Offset);
						Vec3 x(point->pNode->GetXCurr() + f);
						Vec3 Dx(x - xRef);
						Mat3x3 DR(RRef.MulTM(point->pNode->GetRCurr()));
						Vec3 v(point->pNode->GetVCurr() + point->pNode->GetWCurr().Cross(f));
						Vec3 Dv(v - xpRef - wRef.Cross(Dx));
						const Vec3& w(point->pNode->GetWCurr());

						out << GetLabel() << "." << point->uLabel
							<< " " << RRef.MulTV(Dx);

						switch (uRot) {
						case MBC_ROT_MAT:
							out << " " << DR;
							break;

						case MBC_ROT_THETA:
							out << " " << RotManip::VecRot(DR);
							break;

						case MBC_ROT_NONE:
						case MBC_ROT_EULER_123:
							out << " " << MatR2EulerAngles123(DR)*dRaDegr;
							break;
						}

						out << " " << RRef.MulTV(Dv)
							<< " " << RRef.MulTV(w - wRef)
							<< std::endl;
					}

				} else {
					for (std::vector<PointData>::const_iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
						const Mat3x3& R = point->pNode->GetRCurr();
						Vec3 f = R*point->Offset;
						Vec3 x = point->pNode->GetXCurr() + f;
						const Vec3& w = point->pNode->GetWCurr();
						Vec3 wCrossf = w.Cross(f);
						Vec3 v = point->pNode->GetVCurr() + wCrossf;

						out << GetLabel() << "." << point->uLabel
							<< " " << x;

						switch (uRot) {
						case MBC_ROT_MAT:
							out << " " << R;
							break;

						case MBC_ROT_THETA:
							out << " " << RotManip::VecRot(R);
							break;

						case MBC_ROT_NONE:
						case MBC_ROT_EULER_123:
							out << " " << MatR2EulerAngles123(R)*dRaDegr;
							break;
						}
						out << " " << v
							<< " " << w
							<< std::endl;
					}
				}
			}
		}

		/* TODO: NetCDF */
	}
}
 
void
StructExtForce::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(m_Points.size());
	for (unsigned int i = 0; i < m_Points.size(); i++) {
		connectedNodes[i] = m_Points[i].pNode;
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

	const StructNode *pRefNode(0);
	if (HP.IsKeyWord("reference" "node")) {
		pRefNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	}

	bool bSorted(true);
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

			if (!HP.GetYesNo(bSorted)) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"sorted\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
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

			if (!HP.GetYesNo(bLabels)) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"labels\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bGotLabels = true;

		} else if (HP.IsKeyWord("orientation")) {
			if (bGotRot) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"orientation\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (HP.IsKeyWord("none")) {
				uRot = MBC_ROT_NONE;

			} else if (HP.IsKeyWord("orientation" "vector")) {
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

			if (!HP.GetYesNo(bOutputAccelerations)) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"accelerations\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
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
			
			if (!HP.GetYesNo(bUseReferenceNodeForces)) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"\"use reference node forces\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bGotUseRefForces = true;

			if (bUseReferenceNodeForces && HP.IsKeyWord("rotate" "reference" "node" "forces")) {
				if (!HP.GetYesNo(bRotateReferenceNodeForces)) {
					silent_cerr("StructExtForce(" << uLabel << "): "
						"\"rotate reference node forces\" must be either \"yes\" or \"no\" at line "
						<< HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

		} else {
			break;
		}
	}

	if (!bLabels && !bSorted) {
		silent_cerr("StructExtForce(" << uLabel << "): "
			"\"no labels\" and \"unsorted\" incompatible" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("StructExtForce(" << uLabel << "): illegal node number " << n <<
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<unsigned> Labels;
	std::vector<unsigned>::iterator curr_label;
	if (bLabels) {
		Labels.resize(n);
 		curr_label = Labels.begin();
	}
	std::vector<const StructNode *> Nodes(n);
	std::vector<Vec3> Offsets(n);

	for (int i = 0; i < n; i++ ) {
		Nodes[i] = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
		
		ReferenceFrame RF(Nodes[i]);

		if (bLabels) {
			unsigned ul;
			if (HP.IsKeyWord("label")) {
				int l = HP.GetInt();
				if (l < 0) {
					silent_cerr("StructExtForce(" << uLabel << "): "
						"invalid label for point #" << i << " at line "
						<< HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				ul = unsigned(l);

			} else {
				ul = Nodes[i]->GetLabel();
			}

			if (std::find(Labels.begin(), curr_label, ul) < curr_label) {
				silent_cerr("StructExtForce(" << uLabel << "): "
					"duplicate label for point #" << i << " at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			*curr_label = ul;
			++curr_label;

		}

		if (HP.IsKeyWord("offset")) {
			Offsets[i] = HP.GetPosRel(RF);

		} else {
			Offsets[i] = Vec3(Zero3);
		}

		for (int j = 0; j < i; j++) {
			if (Nodes[j] == Nodes[i]) {
				if (Offsets[j].IsExactlySame(Offsets[i])) {
					silent_cerr("StructExtForce(" << uLabel << "): "
						"warning, point #" << i << " is identical to point #" << j << " (same node, same offset)" << std::endl);
				} else {
					silent_cerr("StructExtForce(" << uLabel << "): "
						"warning, point #" << i << " is based on same node of point #" << j << " (offsets differ)" << std::endl);
				}
			}
		}
	}

	std::ofstream out;
	if (HP.IsKeyWord("echo")) {
		const char *s = HP.GetFileName();
		if (s == NULL) {
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"unable to parse echo file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		out.open(s);
		if (!out) {
			int save_errno = errno;
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"unable to open file \"" << s << "\" (" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#if 0

		std::ofstream out(s);

		out.setf(std::ios::scientific);

		for (unsigned i = 0; i < Nodes.size(); i++) {
			if (bLabels) {
				out << Labels[i];
			} else {
				out << Nodes[i]->GetLabel();
			}

			out << " " << (Nodes[i]->GetXCurr() + Offsets[i]) << " ";
			switch (uRot) {
			case MBC_ROT_MAT:
				out << Nodes[i]->GetRCurr();
				break;

			case MBC_ROT_THETA:
				out << RotManip::VecRot(Nodes[i]->GetRCurr());
				break;

			case MBC_ROT_NONE:
			case MBC_ROT_EULER_123:
				out << MatR2EulerAngles123(Nodes[i]->GetRCurr())*dRaDegr;
				break;
			}
			out << " " << Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(Offsets[i])
				<< " " << Nodes[i]->GetWCurr()
				<< std::endl;
		}
#endif
	}

	unsigned uOutputFlags = ExtFileHandlerBase::OUTPUT_REF_DYN;
	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("kinematics")) {
			uOutputFlags |= ExtFileHandlerBase::OUTPUT_KIN;
			fOut = 1;
		} else {
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"unexpected arg at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	StructExtForce *pEl = 0;

	if (dynamic_cast<ExtFileHandlerEDGE *>(pEFH)) {
		SAFENEWWITHCONSTRUCTOR(pEl, StructExtEDGEForce,
			StructExtEDGEForce(uLabel, pDM, pRefNode,
				bUseReferenceNodeForces, bRotateReferenceNodeForces,
				Labels, Nodes, Offsets,
				bSorted, bLabels, bOutputAccelerations, uRot,
				pEFH, bSendAfterPredict, iCoupling, uOutputFlags, fOut));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, StructExtForce,
			StructExtForce(uLabel, pDM, pRefNode,
				bUseReferenceNodeForces, bRotateReferenceNodeForces,
				Labels, Nodes, Offsets,
				bSorted, bLabels, bOutputAccelerations, uRot,
				pEFH, bSendAfterPredict, iCoupling, uOutputFlags, fOut));
	}

	if (out.is_open()) {
		pEl->SendToStream(out, ExtFileHandlerBase::SEND_FIRST_TIME);
	}

	return pEl;
}

