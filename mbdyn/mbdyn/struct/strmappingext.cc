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
#include "strmappingext.h"
#include "modalmappingext.h"
#include "Rot.hh"

#include <fstream>
#include <cerrno>
#include <cstdlib>
#include <ctime>

/* StructMappingExtForce - begin */

/* Costruttore */
StructMappingExtForce::StructMappingExtForce(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<const StructDispNode *>& nodes,
	std::vector<Vec3>& offsets,
	std::vector<unsigned>& labels,
	SpMapMatrixHandler *pH,
	std::vector<uint32_t>& mappedlabels,
	bool bLabels,
	bool bOutputAccelerations,
	unsigned uRRot,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
ExtForce(uL, pDM, pEFH, bSendAfterPredict, iCoupling, fOut), 
pRefNode(pRefNode),
bUseReferenceNodeForces(bUseReferenceNodeForces),
bRotateReferenceNodeForces(bRotateReferenceNodeForces),
F0(Zero3), M0(Zero3),
F1(Zero3), M1(Zero3),
F2(Zero3), M2(Zero3),
pH(pH),
m_uResSize(0),
uPoints(nodes.size()),
uMappedPoints(pH ? unsigned(pH->iGetNumRows())/3 : 0),
bLabels(bLabels),
bOutputAccelerations(bOutputAccelerations),
uRRot(uRRot),
m_qlabels(pH ? mappedlabels : labels),
m_x(3*uPoints),
m_xP(3*uPoints),
m_xPP(0),
m_q(3*uMappedPoints),
m_qP(3*uMappedPoints),
m_qPP(0),
m_f(3*uPoints),
m_p(3*uMappedPoints)
{
	ASSERT(nodes.size() == offsets.size());
	if (pH) {
		ASSERT(3*uPoints == unsigned(pH->iGetNumCols()));
		ASSERT(3*uMappedPoints == unsigned(pH->iGetNumRows()));
	}

	if (pRefNode) {
		switch (uRRot) {
		case MBC_ROT_THETA:
		case MBC_ROT_MAT:
		case MBC_ROT_EULER_123:
			break;

		default:
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"invalid reference node rotation type " << uRRot << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bUseReferenceNodeForces) {
			m_uResSize += 6;
		}
	}

	const StructDispNode *pNode = 0;
	unsigned uNodes = 0;
	std::vector<const StructDispNode *>::const_iterator p;
	for (p = nodes.begin(); p != nodes.end(); ++p) {
		if (*p != pNode) {
			pNode = *p;
			uNodes++;
		}
	}

	Nodes.resize(uNodes);
	p = nodes.begin();
	std::vector<const StructDispNode *>::const_iterator pPrev = p;
	std::vector<NodeData>::iterator n = Nodes.begin();
	while (true) {
		++p;
		if ((p == nodes.end()) || (*p != *pPrev)) {
			n->pNode = *pPrev;
			n->Offsets.resize(p - pPrev);
			for (std::vector<OffsetData>::iterator i = n->Offsets.begin(); i != n->Offsets.end(); ++i) {
				i->Offset = ::Zero3;
			}
			n->F = Zero3;
			n->M = Zero3;

			if (dynamic_cast<const StructNode *>(n->pNode) == 0) {
				switch (n->Offsets.size()) {
				case 1:		// regular rotationless node
				case 2:		// special rotationless node for "membrane"
					break;

				default:
					silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
						"too many offsets for StructDispNode(" << n->pNode->GetLabel() << ")"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				m_uResSize += 3;

			} else {
				m_uResSize += 6;
			}

			if (p == nodes.end()) {
				break;
			}

			++n;
			pPrev = p;
		}
	}

	unsigned uPts = 0;
	n = Nodes.begin();
	std::vector<Vec3>::const_iterator o = offsets.begin();
	std::vector<uint32_t>::iterator l = labels.begin();
	for (; o != offsets.end(); ++o, uPts++) {
		if (uPts == n->Offsets.size()) {
			++n;
			uPts = 0;

			if (dynamic_cast<const StructNode *>(n->pNode) == 0) {
				for (std::vector<StructMappingExtForce::OffsetData>::const_iterator i = n->Offsets.begin();
					i != n->Offsets.end(); i++)
				{
					if (!i->Offset.IsNull()) {
						silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
							"offset #" << (i - n->Offsets.begin())
							<< " (" << i->Offset << ") "
							"for StructDispNode(" << n->pNode->GetLabel() << ") is not zero"
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}
		}

		// FIXME: pass labels
		n->Offsets[uPts].uLabel = unsigned(-1);
		n->Offsets[uPts].Offset = *o;
		n->Offsets[uPts].F = Zero3;

		if (bLabels) {
			n->Offsets[uPts].uLabel = *l;
			++l;
		}
	}

	if (bOutputAccelerations) {
		for (unsigned i = 0; i < nodes.size(); i++) {
			const DynamicStructDispNode *pDSN = dynamic_cast<const DynamicStructDispNode *>(Nodes[i].pNode);
			if (pDSN == 0) {
				silent_cerr("StructMappingExtForce"
					"(" << GetLabel() << "): "
					"StructNode(" << Nodes[i].pNode->GetLabel() << ") "
					"is not dynamic"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			const_cast<DynamicStructDispNode *>(pDSN)->ComputeAccelerations(true);
		}

		m_xPP.resize(3*uPoints);
		m_qPP.resize(3*uMappedPoints);
		
	}
}

StructMappingExtForce::~StructMappingExtForce(void)
{
	NO_OP;
}

void
StructMappingExtForce::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	if (iCoupling == COUPLING_NONE) {
		*piNumRows = 0;
		*piNumCols = 0;

	} else {
		*piNumRows = (pRefNode ? 6 : 0) + 6*Nodes.size();
		*piNumCols = 1;
	}
}

bool
StructMappingExtForce::Prepare(ExtFileHandlerBase *pEFH)
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
			uint32_ptr[0] = MBC_NODAL | MBC_U_ROT_2_REF_NODE_ROT(uRRot);
			if (pRefNode != 0) {
				uint32_ptr[0] |= MBC_REF_NODE;
			}

			if (bLabels) {
				uint32_ptr[0] |= MBC_LABELS;
			}

			if (bOutputAccelerations) {
				uint32_ptr[0] |= MBC_ACCELS;
			}

			uint32_ptr[1] = uPoints;

			ssize_t rc = send(pEFH->GetOutFileDes(),
				(const void *)buf, sizeof(buf),
				pEFH->GetSendFlags());
			if (rc == -1) {
				int save_errno = errno;
				char *err_msg = strerror(save_errno);
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
					"negotiation request send() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
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
		unsigned uRR;
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
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
					"negotiation response recv() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (rc != sizeof(buf)) {
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
					"negotiation response recv() failed "
					"(got " << rc << " of " << sizeof(buf) << " bytes)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			uint32_ptr = (uint32_t *)&buf[0];
			uNodal = (uint32_ptr[0] & MBC_MODAL_NODAL_MASK);
			bRef = (uint32_ptr[0] & MBC_REF_NODE);
			uRR = (uint32_ptr[0] & MBC_REF_NODE_ROT_MASK);
			uR = (uint32_ptr[0] & MBC_ROT_MASK);
			bL = (uint32_ptr[0] & MBC_LABELS);
			bA = (uint32_ptr[0] & MBC_ACCELS);

			uN = uint32_ptr[1];
#endif // USE_SOCKET
		}

		if (uNodal != MBC_NODAL) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: expecting MBC_NODAL "
				"(=" << MBC_MODAL << "), got " << uNodal
				<< std::endl);
			bResult = false;
		}

		if ((pRefNode != 0 && !bRef) || (pRefNode == 0 && bRef)) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: reference node configuration mismatch "
				"(local=" << (pRefNode != 0 ? "yes" : "no") << ", remote=" << (bRef ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (uR != MBC_ROT_NONE) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: orientation output mismatch "
				"(local=" << MBC_ROT_NONE << ", remote=" << uR << ")"
				<< std::endl);
			bResult = false;
		}

		if (uRR != MBC_U_ROT_2_REF_NODE_ROT(uRRot)) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: reference node orientation output mismatch "
				"(local=" << MBC_U_ROT_2_REF_NODE_ROT(uRRot) << ", remote=" << uR << ")"
				<< std::endl);
			bResult = false;
		}

		if (bL != bLabels) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: labels output mismatch "
				"(local=" << (bLabels ? "yes" : "no") << ", remote=" << (bL ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		if (bA != bOutputAccelerations) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: acceleration output mismatch "
				"(local=" << (bOutputAccelerations ? "yes" : "no") << ", remote=" << (bA ? "yes" : "no") << ")"
				<< std::endl);
			bResult = false;
		}

		unsigned uMP = pH ? uMappedPoints : uPoints;
		if (uN != uMP) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"negotiation response failed: node number mismatch "
				"(local=" << uMP << ", remote=" << uN << ")"
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
StructMappingExtForce::Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when)
{
	std::ostream *outfp = pEFH->GetOutStream();
	if (outfp) {
		SendToStream(*outfp, when);

	} else {
		SendToFileDes(pEFH->GetOutFileDes(), when);
	}
}

void
StructMappingExtForce::SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when)
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
			<< " " << RRef
			<< " " << xpRef
			<< " " << wRef;

		if (bOutputAccelerations) {
			outf
				<< " " << xppRef
				<< " " << wpRef;
		}
		outf << std::endl;

		for (unsigned i = 0; i < Nodes.size(); i++) {
#if 0
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
				const Vec3& xpp(Nodes[i]->GetXPPCurr());

				outf
					<< " " << RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2));
				if (uRot != MBC_ROT_NONE) {
					const Vec3& wp(Nodes[i]->GetWPCurr());

					outf
						<< " " << RRef.MulTV(wp - wpRef - wRef.Cross(w));
				}
			}
			outf << std::endl;
#endif
		}

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
#if 0
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
				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 a = Nodes[i]->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				outf
					<< " " << a;

				if (uRot != MBC_ROT_NONE) {
					outf
						<< " " << wp;
				}
			}

			outf << std::endl;
#endif
		}
	}
}

void
StructMappingExtForce::SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when)
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
		switch (uRRot) {
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

		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Vec3 x(pNode->GetXCurr() + f);
					Vec3 Dx(x - xRef);
					Vec3 v(pNode->GetVCurr() + pNode->GetWCurr().Cross(f));
					Vec3 Dv(v - xpRef - wRef.Cross(Dx));
					const Vec3& w(pNode->GetWCurr());

					Vec3 xTilde(RRef.MulTV(Dx));
					m_x.Put(p3 + 1, xTilde);

					Vec3 vTilde(RRef.MulTV(Dv));
					m_xP.Put(p3 + 1, vTilde);

					if (bOutputAccelerations) {
						const Vec3& xpp = pNode->GetXPPCurr();
						const Vec3& wp = pNode->GetWPCurr();

						Vec3 xppTilde(RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2)
							+ wp.Cross(f) + w.Cross(w.Cross(f))));
						m_xPP.Put(p3 + 1, xppTilde);
					}
				}

			} else {
				Vec3 Dx(Nodes[n].pNode->GetXCurr() - xRef);
				Vec3 Dv(Nodes[n].pNode->GetVCurr() - xpRef - wRef.Cross(Dx));

				Vec3 xTilde(RRef.MulTV(Dx));
				m_x.Put(p3 + 1, xTilde);

				Vec3 vTilde(RRef.MulTV(Dv));
				m_xP.Put(p3 + 1, vTilde);

				if (bOutputAccelerations) {
					Vec3 xppTilde(RRef.MulTV(Nodes[n].pNode->GetXPPCurr()
						- xppRef - wpRef.Cross(Dx)
						- wRef.Cross(wRef.Cross(Dx) + Dv*2)));
					m_xPP.Put(p3 + 1, xppTilde);
				}

				p3 += 3;
			}
		}

	} else {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 +=3 ) {
					/*
						p = x + f
						R = R
						v = xp + w cross f
						w = w
						a = xpp + wp cross f + w cross w cross f
						wp = wp
					 */

					// Optimization of the above formulas
					const Mat3x3& R = pNode->GetRCurr();
					Vec3 f = R*Nodes[n].Offsets[o].Offset;
					Vec3 x = pNode->GetXCurr() + f;
					const Vec3& w = pNode->GetWCurr();
					Vec3 wCrossf = w.Cross(f);
					Vec3 v = pNode->GetVCurr() + wCrossf;

					m_x.Put(p3 + 1, x);
					m_xP.Put(p3 + 1, v);

					if (bOutputAccelerations) {
						const Vec3& wp = pNode->GetWPCurr();
						Vec3 xpp = pNode->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

						m_xPP.Put(p3 + 1, xpp);
					}
				}

			} else {
				m_x.Put(p3 + 1, Nodes[n].pNode->GetXCurr());
				m_xP.Put(p3 + 1, Nodes[n].pNode->GetVCurr());

				if (bOutputAccelerations) {
					m_xPP.Put(p3 + 1, Nodes[n].pNode->GetXPPCurr());
				}

				p3 += 3;
			}
		}
	}

	if (bLabels) {
		send(outfd, &m_qlabels[0], sizeof(uint32_t)*m_qlabels.size(), 0);
	}

	if (pH) {
		pH->MatVecMul(m_q, m_x);
		pH->MatVecMul(m_qP, m_xP);

		send(outfd, &m_q[0], sizeof(double)*m_q.size(), 0);
		send(outfd, &m_qP[0], sizeof(double)*m_qP.size(), 0);

		if (bOutputAccelerations) {
			pH->MatVecMul(m_qPP, m_xPP);
			send(outfd, &m_qPP[0], sizeof(double)*m_qPP.size(), 0);
		}

	} else {
		send(outfd, &m_x[0], sizeof(double)*m_x.size(), 0);
		send(outfd, &m_xP[0], sizeof(double)*m_xP.size(), 0);

		if (bOutputAccelerations) {
			send(outfd, &m_xPP[0], sizeof(double)*m_xPP.size(), 0);
		}
	}

#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

void
StructMappingExtForce::Recv(ExtFileHandlerBase *pEFH)
{
	std::istream *infp = pEFH->GetInStream();
	if (infp) {
		RecvFromStream(*infp);

	} else {
		RecvFromFileDes(pEFH->GetInFileDes());
	}
}

void
StructMappingExtForce::RecvFromStream(std::istream& inf)
{
#if 0
	if (pRefNode) {
		unsigned l;
		doublereal *f = F0.pGetVec(), *m = M0.pGetVec();

		if (bLabels) {
			inf >> l;
		}

		inf >> f[0] >> f[1] >> f[2];
		if (uRRot != MBC_ROT_NONE) {
			inf >> m[0] >> m[1] >> m[2];
		}
	}

	for (unsigned i = 0; i < Nodes.size(); i++) {
		/* assume unsigned int label */
		unsigned l;
		doublereal f[3], m[3];

		if (bLabels) {
			inf >> l;

			if (Nodes[i]->GetLabel() != l) {
				silent_cerr("StructMappingExtForce"
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

		F[i] = Vec3(f);
		if (uRot != MBC_ROT_NONE) {
			M[i] = Vec3(m);
		}
	}
#endif
}

void
StructMappingExtForce::RecvFromFileDes(int infd)
{
#ifdef USE_SOCKET
	if (pRefNode) {
		size_t ulen = 0;
		char buf[sizeof(uint32_t) + 6*sizeof(doublereal)];
		doublereal *f;
		ssize_t len;

		if (bLabels) {
			ulen = sizeof(uint32_t);
		}

		ulen += 6*sizeof(doublereal);

		len = recv(infd, (void *)buf, ulen, pEFH->GetRecvFlags());
		if (len == -1) {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"recv() failed (" << save_errno << ": "
				<< err_msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (unsigned(len) != ulen) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"recv() failed " "(got " << len << " of "
				<< ulen << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bLabels) {
			uint32_t *uint32_ptr = (uint32_t *)buf;
			unsigned l = uint32_ptr[0];
			if (l != pRefNode->GetLabel()) {
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
					"invalid reference node label "
					"(wanted " << pRefNode->GetLabel() << ", got " << l << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			f = (doublereal *)&uint32_ptr[1];

		} else {
			f = (doublereal *)buf;
		}

		F0 = Vec3(&f[0]);
		M0 = Vec3(&f[3]);
	}

	if (bLabels) {
		// Hack!
		ssize_t len = recv(infd, (void *)&m_p[0], sizeof(uint32_t)*m_p.size(),
			pEFH->GetRecvFlags());
		if (len == -1) {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"recv() failed (" << save_errno << ": "
				<< err_msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (unsigned(len) != sizeof(double)*m_p.size()) {
			silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
				"recv() failed " "(got " << len << " of "
				<< sizeof(uint32_t)*m_p.size() << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		uint32_t *labels = (uint32_t *)&m_p[0];
		for (unsigned l = 0; l < m_qlabels.size(); l++) {
			if (labels[l] != m_qlabels[l]) {
				silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
					"label mismatch for point #" << l << "/" << m_qlabels.size()
					<< " local=" << m_qlabels[l] << " remote=" << labels[l] << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	size_t fsize;
	double *fp;

	if (pH) {
		fp = &m_p[0];
		fsize = sizeof(double)*m_p.size();

	} else {
		fp = &m_f[0];
		fsize = sizeof(double)*m_f.size();
	}

	ssize_t len = recv(infd, (void *)fp, fsize, pEFH->GetRecvFlags());
	if (len == -1) {
		int save_errno = errno;
		char *err_msg = strerror(save_errno);
		silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
			"recv() failed (" << save_errno << ": "
			<< err_msg << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	} else if (unsigned(len) != fsize) {
		silent_cerr("StructMappingExtForce(" << GetLabel() << "): "
			"recv() failed " "(got " << len << " of "
			<< fsize << " bytes)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pH) {
		pH->MatTVecMul(m_f, m_p);
	}

	if (pRefNode) {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			Nodes[n].F = Zero3;
			Nodes[n].M = Zero3;
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Nodes[n].Offsets[o].F = pRefNode->GetRCurr()*Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[o].F;
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Nodes[n].M += f.Cross(Nodes[n].Offsets[o].F);
				}

			} else {
				Nodes[n].Offsets[0].F = pRefNode->GetRCurr()*Vec3(&m_f[p3]);
				Nodes[n].F += Nodes[n].Offsets[0].F;
			}
		}

	} else {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			Nodes[n].F = Zero3;
			Nodes[n].M = Zero3;
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Nodes[n].Offsets[o].F = Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[o].F;
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Nodes[n].M += f.Cross(Nodes[n].Offsets[o].F);
				}

			} else {
					Nodes[n].Offsets[0].F = Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[0].F;
			}
		}
	}
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

SubVectorHandler&
StructMappingExtForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	ExtForce::Recv();

	if (iCoupling == COUPLING_NONE) {
		WorkVec.Resize(0);
		return WorkVec;
	}


	if (pRefNode) {
		integer iSize(0);

		WorkVec.ResizeReset(m_uResSize);

		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();

		// manipulate
		if (bUseReferenceNodeForces) {
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
		for (unsigned n = 0; n < Nodes.size(); n++) {
			integer iFirstIndex = Nodes[n].pNode->iGetFirstMomentumIndex();
			integer iDim;

			WorkVec.Add(iSize + 1, Nodes[n].F);
			if (dynamic_cast<const StructNode *>(Nodes[n].pNode)) {
				WorkVec.Add(iSize + 4, Nodes[n].M);
				iDim = 6;

			} else {
				iDim = 3;
			}

			for (int r = 1; r <= iDim; r++) {
				WorkVec.PutRowIndex(iSize + r, iFirstIndex + r);
			}

			if (bUseReferenceNodeForces || fToBeOutput()) {
				// compute Global Reference Node Forces, even if they are not used, for output only :)
				F2 += Nodes[n].F;
				M2 += Nodes[n].M + (Nodes[n].pNode->GetXCurr() - xRef).Cross(Nodes[n].F);
			}

			iSize += iDim;
		}

		if (bUseReferenceNodeForces) {
			integer iFirstIndex = pRefNode->iGetFirstMomentumIndex();
			for (int r = 1; r <= 6; r++) {
				WorkVec.PutRowIndex(iSize + r, iFirstIndex + r);
			}

			F1 -= F2;
			M1 -= M2;
			WorkVec.Add(iSize + 1, F1);
			WorkVec.Add(iSize + 4, M1);

			iSize += 6;
		}

		ASSERT(iSize == m_uResSize);

	} else {
		integer iSize(0);

		WorkVec.ResizeReset(m_uResSize);

		for (unsigned n = 0; n < Nodes.size(); n++) {
			integer iFirstIndex = Nodes[n].pNode->iGetFirstMomentumIndex();
			integer iDim;

			WorkVec.Add(iSize + 1, Nodes[n].F);
			if (dynamic_cast<const StructNode *>(Nodes[n].pNode)) {
				WorkVec.Add(iSize + 4, Nodes[n].M);
				iDim = 6;

			} else {
				iDim = 3;
			}

			for (int r = 1; r <= iDim; r++) {
				WorkVec.PutRowIndex(iSize + r, iFirstIndex + r);
			}

			iSize += iDim;
		}

		ASSERT(iSize == m_uResSize);
	}

	return WorkVec;
}

void
StructMappingExtForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::FORCES)) {
			std::ostream& out = OH.Forces();

			if (pRefNode) {
				out << GetLabel() << "#" << pRefNode->GetLabel()
					<< " " << F0
					<< " " << M0
					<< " " << F1
					<< " " << M1
					<< " " << F2
					<< " " << M2
					<< std::endl;
			}

			for (unsigned n = 0; n < Nodes.size(); n++) {
				out << GetLabel() << "@" << Nodes[n].pNode->GetLabel()
					<< " " << Nodes[n].F
					<< " " << Nodes[n].M
					<< std::endl;
			}
		}

		/* TODO: NetCDF */
	}
}
 
void
StructMappingExtForce::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	unsigned n = Nodes.size();
	if (pRefNode) {
		n++;
	}
	connectedNodes.resize(n);

	for (n = 0; n < Nodes.size(); n++) {
		connectedNodes[n] = Nodes[n].pNode;
	}

	if (pRefNode) {
		connectedNodes[n] = pRefNode;
	}
}

/* StructMappingExtForce - end */


/* StructMembraneMappingExtForce - begin */

/* Costruttore */
StructMembraneMappingExtForce::StructMembraneMappingExtForce(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<const StructDispNode *>& nodes,
	std::vector<Vec3>& offsets,
	std::vector<unsigned>& labels,
	std::vector<NodeConnData>& nodesConn,
	SpMapMatrixHandler *pH,
	std::vector<uint32_t>& mappedlabels,
	bool bLabels,
	bool bOutputAccelerations,
	unsigned uRRot,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
StructMappingExtForce(uL, pDM,
	pRefNode, bUseReferenceNodeForces, bRotateReferenceNodeForces,
	nodes, offsets, labels, pH, mappedlabels,
	bLabels, bOutputAccelerations, uRRot,
	pEFH, bSendAfterPredict, iCoupling,
	fOut)
{
	NodesConn.resize(Nodes.size());
	for (unsigned n = 0; n < NodesConn.size(); n++) {
		NodesConn[n] = nodesConn[n];
	}
}

StructMembraneMappingExtForce::~StructMembraneMappingExtForce(void)
{
	NO_OP;
}

/*
 * Send output to companion software
 */
void
StructMembraneMappingExtForce::SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when)
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
			<< " " << RRef
			<< " " << xpRef
			<< " " << wRef;

		if (bOutputAccelerations) {
			outf
				<< " " << xppRef
				<< " " << wpRef;
		}
		outf << std::endl;

		for (unsigned i = 0; i < Nodes.size(); i++) {
#if 0
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
				const Vec3& xpp(Nodes[i]->GetXPPCurr());

				outf
					<< " " << RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2));
				if (uRot != MBC_ROT_NONE) {
					const Vec3& wp(Nodes[i]->GetWPCurr());

					outf
						<< " " << RRef.MulTV(wp - wpRef - wRef.Cross(w));
				}
			}
			outf << std::endl;
#endif
		}

	} else {
		for (unsigned i = 0; i < Nodes.size(); i++) {
#if 0
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
				const Vec3& wp = Nodes[i]->GetWPCurr();
				Vec3 a = Nodes[i]->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

				outf
					<< " " << a;

				if (uRot != MBC_ROT_NONE) {
					outf
						<< " " << wp;
				}
			}

			outf << std::endl;
#endif
		}
	}
}

void
StructMembraneMappingExtForce::SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when)
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
		switch (uRRot) {
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

		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Vec3 x(pNode->GetXCurr() + f);
					Vec3 Dx(x - xRef);
					Vec3 v(pNode->GetVCurr() + pNode->GetWCurr().Cross(f));
					Vec3 Dv(v - xpRef - wRef.Cross(Dx));
					const Vec3& w(pNode->GetWCurr());

					Vec3 xTilde(RRef.MulTV(Dx));
					m_x.Put(p3 + 1, xTilde);

					Vec3 vTilde(RRef.MulTV(Dv));
					m_xP.Put(p3 + 1, vTilde);

					if (bOutputAccelerations) {
						const Vec3& xpp = pNode->GetXPPCurr();
						const Vec3& wp = pNode->GetWPCurr();

						Vec3 xppTilde(RRef.MulTV(xpp - xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2)
							+ wp.Cross(f) + w.Cross(w.Cross(f))));
						m_xPP.Put(p3 + 1, xppTilde);
					}
				}

			} else {
				if (NodesConn[n].pNode[0] != 0) {
					Vec3 e1(NodesConn[n].pNode[1]->GetXCurr() - NodesConn[n].pNode[0]->GetXCurr());
					Vec3 e2(NodesConn[n].pNode[3]->GetXCurr() - NodesConn[n].pNode[2]->GetXCurr());
					Vec3 e3(e1.Cross(e2));
					e3 /= e3.Norm();

					Nodes[n].Offsets[0].Offset = e3*NodesConn[n].h1;
					Nodes[n].Offsets[1].Offset = e3*NodesConn[n].h2;

					for (unsigned o = 0; o < 2; o++, p3 += 3) {
						Vec3 Dx(Nodes[n].pNode->GetXCurr() + Nodes[n].Offsets[o].Offset - xRef);
						Vec3 Dv(Nodes[n].pNode->GetVCurr() - xpRef - wRef.Cross(Dx));

						Vec3 xTilde(RRef.MulTV(Dx));
						m_x.Put(p3 + 1, xTilde);

						Vec3 vTilde(RRef.MulTV(Dv));
						m_xP.Put(p3 + 1, vTilde);

						if (bOutputAccelerations) {
							Vec3 xppTilde(RRef.MulTV(Nodes[n].pNode->GetXPPCurr()
								- xppRef - wpRef.Cross(Dx)
								- wRef.Cross(wRef.Cross(Dx) + Dv*2)));
							m_xPP.Put(p3 + 1, xppTilde);
						}
					}

				} else {
					Vec3 Dx(Nodes[n].pNode->GetXCurr() - xRef);
					Vec3 Dv(Nodes[n].pNode->GetVCurr() - xpRef - wRef.Cross(Dx));

					Vec3 xTilde(RRef.MulTV(Dx));
					m_x.Put(p3 + 1, xTilde);

					Vec3 vTilde(RRef.MulTV(Dv));
					m_xP.Put(p3 + 1, vTilde);

					if (bOutputAccelerations) {
						Vec3 xppTilde(RRef.MulTV(Nodes[n].pNode->GetXPPCurr()
							- xppRef - wpRef.Cross(Dx)
							- wRef.Cross(wRef.Cross(Dx) + Dv*2)));
						m_xPP.Put(p3 + 1, xppTilde);
					}

					p3 += 3;
				}
			}
		}

	} else {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 +=3 ) {
					/*
						p = x + f
						R = R
						v = xp + w cross f
						w = w
						a = xpp + wp cross f + w cross w cross f
						wp = wp
					 */

					// Optimization of the above formulas
					const Mat3x3& R = pNode->GetRCurr();
					Vec3 f = R*Nodes[n].Offsets[o].Offset;
					Vec3 x = pNode->GetXCurr() + f;
					const Vec3& w = pNode->GetWCurr();
					Vec3 wCrossf = w.Cross(f);
					Vec3 v = pNode->GetVCurr() + wCrossf;

					m_x.Put(p3 + 1, x);
					m_xP.Put(p3 + 1, v);

					if (bOutputAccelerations) {
						const Vec3& wp = pNode->GetWPCurr();
						Vec3 xpp = pNode->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

						m_xPP.Put(p3 + 1, xpp);
					}
				}

			} else {
				if (NodesConn[n].pNode[0] != 0) {
					Vec3 e1(NodesConn[n].pNode[1]->GetXCurr() - NodesConn[n].pNode[0]->GetXCurr());
					Vec3 e2(NodesConn[n].pNode[3]->GetXCurr() - NodesConn[n].pNode[2]->GetXCurr());
					Vec3 e3(e1.Cross(e2));
					e3 /= e3.Norm();

					Nodes[n].Offsets[0].Offset = e3*NodesConn[n].h1;
					Nodes[n].Offsets[1].Offset = e3*NodesConn[n].h2;

					for (unsigned o = 0; o < 2; o++, p3 += 3) {
						m_x.Put(p3 + 1, Nodes[n].pNode->GetXCurr() + Nodes[n].Offsets[o].Offset);
						m_xP.Put(p3 + 1, Nodes[n].pNode->GetVCurr());

						if (bOutputAccelerations) {
							m_xPP.Put(p3 + 1, Nodes[n].pNode->GetXPPCurr());
						}
					}

				} else {
					m_x.Put(p3 + 1, Nodes[n].pNode->GetXCurr());
					m_xP.Put(p3 + 1, Nodes[n].pNode->GetVCurr());

					if (bOutputAccelerations) {
						m_xPP.Put(p3 + 1, Nodes[n].pNode->GetXPPCurr());
					}

					p3 += 3;
				}
			}
		}
	}

	if (bLabels) {
		send(outfd, &m_qlabels[0], sizeof(uint32_t)*m_qlabels.size(), 0);
	}

	if (pH) {
		pH->MatVecMul(m_q, m_x);
		pH->MatVecMul(m_qP, m_xP);

		send(outfd, &m_q[0], sizeof(double)*m_q.size(), 0);
		send(outfd, &m_qP[0], sizeof(double)*m_qP.size(), 0);

		if (bOutputAccelerations) {
			pH->MatVecMul(m_qPP, m_xPP);
			send(outfd, &m_qPP[0], sizeof(double)*m_qPP.size(), 0);
		}

	} else {
		send(outfd, &m_x[0], sizeof(double)*m_x.size(), 0);
		send(outfd, &m_xP[0], sizeof(double)*m_xP.size(), 0);

		if (bOutputAccelerations) {
			send(outfd, &m_xPP[0], sizeof(double)*m_xPP.size(), 0);
		}
	}

#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

void
StructMembraneMappingExtForce::RecvFromStream(std::istream& inf)
{
#if 0
	if (pRefNode) {
		unsigned l;
		doublereal *f = F0.pGetVec(), *m = M0.pGetVec();

		if (bLabels) {
			inf >> l;
		}

		inf >> f[0] >> f[1] >> f[2];
		if (uRRot != MBC_ROT_NONE) {
			inf >> m[0] >> m[1] >> m[2];
		}
	}

	for (unsigned i = 0; i < Nodes.size(); i++) {
		/* assume unsigned int label */
		unsigned l;
		doublereal f[3], m[3];

		if (bLabels) {
			inf >> l;

			if (Nodes[i]->GetLabel() != l) {
				silent_cerr("StructMembraneMappingExtForce"
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

		F[i] = Vec3(f);
		if (uRot != MBC_ROT_NONE) {
			M[i] = Vec3(m);
		}
	}
#endif
}

void
StructMembraneMappingExtForce::RecvFromFileDes(int infd)
{
#ifdef USE_SOCKET
	if (pRefNode) {
		size_t ulen = 0;
		char buf[sizeof(uint32_t) + 6*sizeof(doublereal)];
		doublereal *f;
		ssize_t len;

		if (bLabels) {
			ulen = sizeof(uint32_t);
		}

		ulen += 6*sizeof(doublereal);

		len = recv(infd, (void *)buf, ulen, pEFH->GetRecvFlags());
		if (len == -1) {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);
			silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
				"recv() failed (" << save_errno << ": "
				<< err_msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (unsigned(len) != ulen) {
			silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
				"recv() failed " "(got " << len << " of "
				<< ulen << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bLabels) {
			uint32_t *uint32_ptr = (uint32_t *)buf;
			unsigned l = uint32_ptr[0];
			if (l != pRefNode->GetLabel()) {
				silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
					"invalid reference node label "
					"(wanted " << pRefNode->GetLabel() << ", got " << l << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			f = (doublereal *)&uint32_ptr[1];

		} else {
			f = (doublereal *)buf;
		}

		F0 = Vec3(&f[0]);
		M0 = Vec3(&f[3]);
	}

	if (bLabels) {
		// Hack!
		ssize_t len = recv(infd, (void *)&m_p[0], sizeof(uint32_t)*m_p.size(),
			pEFH->GetRecvFlags());
		if (len == -1) {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);
			silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
				"recv() failed (" << save_errno << ": "
				<< err_msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (unsigned(len) != sizeof(double)*m_p.size()) {
			silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
				"recv() failed " "(got " << len << " of "
				<< sizeof(uint32_t)*m_p.size() << " bytes)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		uint32_t *labels = (uint32_t *)&m_p[0];
		for (unsigned l = 0; l < m_qlabels.size(); l++) {
			if (labels[l] != m_qlabels[l]) {
				silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
					"label mismatch for point #" << l << "/" << m_qlabels.size()
					<< " local=" << m_qlabels[l] << " remote=" << labels[l] << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	size_t fsize;
	double *fp;

	if (pH) {
		fp = &m_p[0];
		fsize = sizeof(double)*m_p.size();

	} else {
		fp = &m_f[0];
		fsize = sizeof(double)*m_f.size();
	}

	ssize_t len = recv(infd, (void *)fp, fsize, pEFH->GetRecvFlags());
	if (len == -1) {
		int save_errno = errno;
		char *err_msg = strerror(save_errno);
		silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
			"recv() failed (" << save_errno << ": "
			<< err_msg << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	} else if (unsigned(len) != fsize) {
		silent_cerr("StructMembraneMappingExtForce(" << GetLabel() << "): "
			"recv() failed " "(got " << len << " of "
			<< fsize << " bytes)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pH) {
		pH->MatTVecMul(m_f, m_p);
	}

	if (pRefNode) {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			Nodes[n].F = Zero3;
			Nodes[n].M = Zero3;
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Nodes[n].Offsets[o].F = pRefNode->GetRCurr()*Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[o].F;
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Nodes[n].M += f.Cross(Nodes[n].Offsets[o].F);
				}

			} else {
				if (NodesConn[n].pNode[0] != 0) {
					for (unsigned o = 0; o < 2; o++, p3 += 3) {
						Nodes[n].Offsets[o].F = pRefNode->GetRCurr()*Vec3(&m_f[p3]);
						Nodes[n].F += Nodes[n].Offsets[o].F;
					}

				} else {
					Nodes[n].Offsets[0].F = pRefNode->GetRCurr()*Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[0].F;
				}
			}
		}

	} else {
		for (unsigned p3 = 0, n = 0; n < Nodes.size(); n++) {
			Nodes[n].F = Zero3;
			Nodes[n].M = Zero3;
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[n].pNode));
			if (pNode != 0) {
				for (unsigned o = 0; o < Nodes[n].Offsets.size(); o++, p3 += 3) {
					Nodes[n].Offsets[o].F = Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[o].F;
					Vec3 f(pNode->GetRCurr()*Nodes[n].Offsets[o].Offset);
					Nodes[n].M += f.Cross(Nodes[n].Offsets[o].F);
				}

			} else {
				if (NodesConn[n].pNode[0] != 0) {
					for (unsigned o = 0; o < 2; o++, p3 += 3) {
						Nodes[n].Offsets[o].F = Vec3(&m_f[p3]);
						Nodes[n].F += Nodes[n].Offsets[o].F;
					}

				} else {
					Nodes[n].Offsets[0].F = Vec3(&m_f[p3]);
					Nodes[n].F += Nodes[n].Offsets[0].F;
				}
			}
		}
	}
#else // ! USE_SOCKET
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}

/* StructMembraneMappingExtForce - end */


Elem*
ReadStructMappingExtForce(DataManager* pDM, 
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

	bool bLabels(false);
	unsigned uRRot = pRefNode ? MBC_ROT_MAT : MBC_ROT_NONE;
	bool bOutputAccelerations(false);
	bool bUseReferenceNodeForces(true);
	bool bRotateReferenceNodeForces(true);

	bool bGotLabels(false);
	bool bGotRot(false);
	bool bGotAccels(false);
	bool bGotUseRefForces(false);

	while (HP.IsArg()) {
		if (HP.IsKeyWord("no" "labels")) {
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"use of \"no labels\" deprecated in favor of \"labels, { yes | no }\" at line "
				<< HP.GetLineData() << std::endl);

			if (bGotLabels) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"no labels\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bLabels = false;
			bGotLabels = true;

		} else if (HP.IsKeyWord("labels")) {
			if (bGotLabels) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"labels\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!HP.GetYesNo(bLabels)) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"labels\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bGotLabels = true;

		} else if (HP.IsKeyWord("orientation")) {
			if (bGotRot) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"orientation\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (HP.IsKeyWord("none")) {
				uRRot = MBC_ROT_NONE;

			} else if (HP.IsKeyWord("orientation" "vector")) {
				uRRot = MBC_ROT_THETA;

			} else if (HP.IsKeyWord("orientation" "matrix")) {
				uRRot = MBC_ROT_MAT;

			} else if (HP.IsKeyWord("euler" "123")) {
				uRRot = MBC_ROT_EULER_123;

			} else {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"unknown \"orientation\" format at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			bGotRot = true;	

		} else if (HP.IsKeyWord("accelerations")) {
			if (bGotAccels) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"accelerations\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (!HP.GetYesNo(bOutputAccelerations)) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"accelerations\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bGotAccels = true;

		} else if (HP.IsKeyWord("use" "reference" "node" "forces")) {
			if (pRefNode == 0) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"use reference node forces\" only meaningful when reference node is used at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (bGotUseRefForces) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"use reference node forces\" already specified at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			
			if (!HP.GetYesNo(bUseReferenceNodeForces)) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"\"use reference node forces\" must be either \"yes\" or \"no\" at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bGotUseRefForces = true;

			if (bUseReferenceNodeForces && HP.IsKeyWord("rotate" "reference" "node" "forces")) {
				if (!HP.GetYesNo(bRotateReferenceNodeForces)) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"rotate reference node forces\" must be either \"yes\" or \"no\" at line "
						<< HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

		} else {
			break;
		}
	}

	if (!HP.IsKeyWord("points" "number")) {
		silent_cerr("StructMappingExtForce(" << uLabel << "): "
			"\"points number\" keyword expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int nPoints = HP.GetInt();
	if (nPoints <= 0) {
		silent_cerr("StructMappingExtForce(" << uLabel << "): illegal points number " << nPoints <<
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<const StructDispNode *> Nodes(nPoints);
	std::vector<Vec3> Offsets(nPoints, ::Zero3);
	std::vector<uint32_t> Labels;
	if (bLabels) {
		Labels.resize(nPoints);
	}
	bool bMembrane(false);
	std::vector<StructMembraneMappingExtForce::NodeConnData> NodesConn;

	std::map<unsigned, bool> Got;

	for (unsigned n = 0, p = 0; p < unsigned(nPoints); n++) {
		const StructDispNode *pNode = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);
		unsigned uL(pNode->GetLabel());
		if (Got[uL]) {
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"StructNode(" << uL << ") out of order "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		Got[uL] = true;

		StructMembraneMappingExtForce::NodeConnData ncd;
		ncd.pNode[0] = 0;

		ReferenceFrame RF(pNode);
		if (HP.IsKeyWord("membrane")) {
			bMembrane = true;

			for (unsigned n = 0; n < 4; n++) {
				const StructDispNode *pNn = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);
				if ((n%2) && pNn == ncd.pNode[n - 1]) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"nodes #" << n << " and #" << n - 1 << " are the same in \"membrane\" mapping for StructNode(" << uL << ") at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				ncd.pNode[n] = pNn;
			}

			// first node
			Nodes[p] = pNode;
			Offsets[p] = ::Zero3;

			doublereal h1 = HP.GetReal();
			if (h1 == 0.) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"invalid h1 in \"membrane\" mapping for StructNode(" << uL << ") at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			ncd.h1 = h1;

			if (bLabels) {
				int l = HP.GetInt();
				if (l < 0) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid (negative) point label " << l
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				Labels[p] = l;
			}

			p++;

			if (p == unsigned(nPoints)) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"second point in \"membrane\" mapping for StructNode(" << uL << ") exceeds points number " << nPoints << " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			// second node
			Nodes[p] = pNode;
			Offsets[p] = ::Zero3;

			doublereal h2 = HP.GetReal();
			if (h2 == 0. || h1*h2 > 0.) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"invalid h2 in \"membrane\" mapping for StructNode(" << uL << ") at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			ncd.h2 = h2;

			if (bLabels) {
				int l = HP.GetInt();
				if (l < 0) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid (negative) point label " << l
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				Labels[p] = l;
			}

			p++;

		} else {
			while (HP.IsKeyWord("offset")) {
				if (p == unsigned(nPoints)) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"point " << p << " offset from StructNode(" << pNode->GetLabel() << ") exceeds expected value " << nPoints
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bLabels) {
					int l = HP.GetInt();
					if (l < 0) {
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"invalid (negative) point label " << l
							<< " at line " << HP.GetLineData() << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					Labels[p] = l;
				}

				Nodes[p] = pNode;
				Offsets[p] = HP.GetPosRel(RF);
				p++;
			}
		}

		NodesConn.push_back(ncd);
	}

	if (HP.IsKeyWord("echo")) {
		const char *s = HP.GetFileName();
		if (s == NULL) {
			silent_cerr("StructMappingExtForce(" << uLabel << "): "
				"unable to parse echo file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::ofstream out(s);

		out.setf(std::ios::scientific);

		bool bGotSurface(false);
		bool bGotOutput(false);
		bool bGotOrder(false);
		bool bGotBaseNode(false);
		bool bGotWeight(false);

		bool bGotPrecision(false);

	/*
	surface: basicgrid.dat
	output: blade1H.dat
	order: 2
	basenode: 12
	weight: 2 

	*/

		std::string surface;
		std::string output;
		int order;
		int basenode;
		int weight;
		bool bWeightInf(false);

		while (true) {
			if (HP.IsKeyWord("precision")) {
				if (bGotPrecision) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"precision\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotPrecision = true;

				int iPrecision = HP.GetInt();
				if (iPrecision <= 0) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid echo precision " << iPrecision
						<< " at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				out.precision(iPrecision);

			} else if (HP.IsKeyWord("surface")) {
				if (bGotSurface) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"surface\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotSurface = true;

				surface = HP.GetFileName();
				if (surface.empty()) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid \"surface\" "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("output")) {
				if (bGotOutput) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"output\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotOutput = true;

				output = HP.GetFileName();
				if (output.empty()) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid \"output\" "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("order")) {
				if (bGotOrder) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"order\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotOrder = true;

				order = HP.GetInt();
				if (order < 1 || order > 3) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid order=" << order << ", must be 1 <= order <= 3 "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("basenode")) {
				if (bGotBaseNode) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"basenode\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotBaseNode = true;

				basenode = HP.GetInt();
				if (basenode <= 0) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"invalid basenode=" << basenode << ", must be > 0 "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("weight")) {
				if (bGotWeight) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"\"weight\" already specified "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				bGotWeight = true;

				if (HP.IsKeyWord("inf")) {
					bWeightInf = true;

				} else {
					weight = HP.GetInt();
					if (weight < -2) {
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"invalid weight=" << weight << ", must be >= -2 or \"inf\" "
							"at line " << HP.GetLineData() << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

			} else {
				break;
			}
		}

		time_t t = std::time(NULL);
		const char *user = std::getenv("USER");
		const char *host = std::getenv("HOSTNAME");
		if (user == 0) user = "nobody";
		if (host == 0) host = "unknown";

		out
			<< "# Generated by MBDyn StructMappingExtForce(" << uLabel << ")" << std::endl
			<< "# " << user << "@" << host << std::endl
			<< "# " << std::ctime(&t)
			<< "# labels: " << (bLabels ? "on" : "off") << std::endl;
		Vec3 xRef(::Zero3);
		Mat3x3 RRef(::Eye3);
		if (pRefNode) {
			xRef = pRefNode->GetXCurr();
			RRef = pRefNode->GetRCurr();

			out 
				<< "# reference: " << pRefNode->GetLabel() << std::endl
			    	<< "# position: " << xRef << std::endl
			    	<< "# orientation: " << MatR2EulerAngles123(RRef)*dRaDegr << std::endl;
		}

		out
			<< "# points: " << nPoints << std::endl;

		if (bGotSurface) {
			out << "# surface: " << surface << std::endl;
		}

		if (bGotOutput) {
			out << "# output: " << output << std::endl;
		}

		if (bGotOrder) {
			out << "# order: " << order << std::endl;
		}

		if (bGotBaseNode) {
			out << "# basenode: " << basenode << std::endl;
		}

		if (bGotWeight) {
			if (bWeightInf) {
				out << "# weight: inf" << std::endl;
			} else {
				out << "# weight: " << weight << std::endl;
			}
		}

		for (unsigned n = 0, p = 0; p < unsigned(nPoints); p++) {
			if (p > 0 && Nodes[p] != Nodes[p - 1]) {
				n++;
			}

			if (bLabels) {
				out << Labels[p] << " ";
			}

			Vec3 x;
			const StructNode *pNode(dynamic_cast<const StructNode *>(Nodes[p]));
			if (pNode != 0) {
				x = RRef.MulTV(pNode->GetXCurr() + pNode->GetRCurr()*Offsets[p] - xRef);

			} else {
				Vec3 h(::Zero3);
				if (NodesConn[n].pNode[0] != 0) {
					Vec3 e1(NodesConn[n].pNode[1]->GetXCurr() - NodesConn[n].pNode[0]->GetXCurr());
					Vec3 e2(NodesConn[n].pNode[3]->GetXCurr() - NodesConn[n].pNode[2]->GetXCurr());
					Vec3 e3(e1.Cross(e2));
					e3 /= e3.Norm();

					if (Nodes[p + 1] == Nodes[p]) {
						h = e3*NodesConn[n].h1;

					} else {
						h = e3*NodesConn[n].h2;
					}
				}

				x = RRef.MulTV(Nodes[p]->GetXCurr() + h - xRef);
			}
			out << x << std::endl;
		}

		if (HP.IsKeyWord("stop")) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	SpMapMatrixHandler *pH = 0;
	std::vector<uint32_t> MappedLabels;
	if (HP.IsKeyWord("mapped" "points" "number")) {
		int nMappedPoints = 0;
		if (HP.IsKeyWord("from" "file")) {
			nMappedPoints = -1;

		} else {
			nMappedPoints = HP.GetInt();
			if (nMappedPoints <= 0) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					"invalid mapped points number " << nMappedPoints
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			nMappedPoints *= 3;
		}

		integer nCols = 3*nPoints;
		pH = ReadSparseMappingMatrix(HP, nMappedPoints, nCols);
		ASSERT((nMappedPoints%3) == 0);
		ASSERT(nCols == 3*nPoints);
		nMappedPoints /= 3;

		if (bLabels) {
			MappedLabels.resize(nMappedPoints);
			if (HP.IsKeyWord("mapped" "labels" "file")) {
				const char *sFileName = HP.GetFileName();
				if (sFileName == 0) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"unable to read mapped labels file name "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				std::ifstream in(sFileName);
				if (!in) {
					silent_cerr("StructMappingExtForce(" << uLabel << "): "
						"unable to open mapped labels file "
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

				for (unsigned l = 0; l < unsigned(nMappedPoints); l++) {
					int i;
					in >> i;
					if (!in) {
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"unable to read mapped label #" << l << "/" << nMappedPoints
							<< " from mapped labels file \"" << sFileName << "\"" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					if (i < 0) {
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"invalid (negative) mapped label #" << l << "/" << nMappedPoints
							<< " from mapped labels file \"" << sFileName << "\"" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					MappedLabels[l] = i;
				}

			} else {
				for (unsigned l = 0; l < unsigned(nMappedPoints); l++) {
					int i = HP.GetInt();
					if (i < 0) {
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"invalid (negative) mapped label #" << l << "/" << nMappedPoints
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					MappedLabels[l] = i;
				}
			}

			int duplicate = 0;
			for (unsigned l = 1; l < unsigned(nMappedPoints); l++) {
				for (unsigned c = 0; c < l; c++) {
					if (MappedLabels[l] == MappedLabels[c]) {
						duplicate++;
						silent_cerr("StructMappingExtForce(" << uLabel << "): "
							"duplicate mapped label " << MappedLabels[l] << ": "
							"#" << l << "==#" << c << std::endl);
					}
				}
			}

			if (duplicate) {
				silent_cerr("StructMappingExtForce(" << uLabel << "): "
					<< duplicate << " duplicate mapped labels"
					<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	if (bMembrane) {
		SAFENEWWITHCONSTRUCTOR(pEl, StructMembraneMappingExtForce,
			StructMembraneMappingExtForce(uLabel, pDM, pRefNode,
				bUseReferenceNodeForces, bRotateReferenceNodeForces,
				Nodes, Offsets, Labels, NodesConn,
				pH,
				MappedLabels,
				bLabels, bOutputAccelerations, uRRot,
				pEFH, bSendAfterPredict, iCoupling, fOut));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, StructMappingExtForce,
			StructMappingExtForce(uLabel, pDM, pRefNode,
				bUseReferenceNodeForces, bRotateReferenceNodeForces,
				Nodes, Offsets, Labels,
				pH,
				MappedLabels,
				bLabels, bOutputAccelerations, uRRot,
				pEFH, bSendAfterPredict, iCoupling, fOut));
	}

	return pEl;
}

