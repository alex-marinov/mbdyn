/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/struct/stredge.cc,v 1.11 2017/01/12 14:46:44 masarati Exp $ */
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

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

//#include <fstream>
//#include <cerrno>
//#include <algorithm>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "dataman.h"
#include "extsharedmem.h"
#include "strextsharedmem.h"
#include "sharedmem.h"
#include "Rot.hh"

using namespace mbdyn;
using namespace boost::interprocess;

/* Constructor */
StructExtSharedMemForce::StructExtSharedMemForce(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<unsigned>& labels,
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
StructExtForce(uL, pDM, pRefNode, bUseReferenceNodeForces, bRotateReferenceNodeForces,
	labels, nodes, offsets, bSorted, bLabels, bOutputAccelerations, uRot,
	pEFH, bSendAfterPredict, iCoupling, uOutputFlags, fOut)
{
    pedantic_cout("In StructExtSharedMemForce::StructExtSharedMemForce" << std::endl);
	ASSERT(dynamic_cast<ExtSharedMemHandler *>(pEFH) != 0);
	ASSERT(bLabels);
	ASSERT(!bSorted);

	pESMH = dynamic_cast<ExtSharedMemHandler *>(pEFH);

	// initialise the external shared memory communicator. If
	// 'create' was true this will allocate the shared memory
	// region, otherwise, it will just open it
	pedantic_cout("In StructExtSharedMemForce::StructExtSharedMemForce, calling pESMH->Initialise" << std::endl);
    pESMH->Initialise (node_kinematics_size, dynamics_size, labels_size, m_Points.size (), uRot);

    buf = pESMH->GetSharedMemBuffer ();

}

StructExtSharedMemForce::~StructExtSharedMemForce(void)
{
	NO_OP;
}

bool
StructExtSharedMemForce::Prepare(ExtFileHandlerBase *pEFH)
{
	bool bResult = true;

    pedantic_cout("StructExtSharedMemForce:in Prepare" << std::endl);

	switch (pEFH->NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_NO:
		break;

	case ExtFileHandlerBase::NEGOTIATE_CLIENT: {

	    scoped_lock<interprocess_mutex> lock(buf->mutex);

        pedantic_cout("[MBDYN] StructExtSharedMemForce::Prepare: checking if peer is ready to receive" << std::endl);

        if (!buf->peer_ready_for_cmd)
        {
            pedantic_cout("[MBDYN] StructExtSharedMemForce::Prepare: waiting for peer be ready to receive" << std::endl);
            buf->cond_peer_ready_for_cmd.wait (lock);
        }

        buf->NodalOrModal = 0;

        buf->RotationType = uRot;

        if (pRefNode != 0) {
            buf->UseReferanceNode = true;
        }

        if (bLabels) {
            buf->UseLabels = true;
        }

        if (bOutputAccelerations) {
            buf->OutputAccelerations = true;
        }

        buf->NumNodes = m_Points.size();

        pedantic_cout("[MBDYN] StructExtSharedMemForce::Prepare: put data in buffer" << std::endl);

        pedantic_cout("[MBDYN] StructExtSharedMemForce::Prepare: notifying new data is available" << std::endl);

        buf->cond_mbdyn_cmd_available.notify_all();
        buf->mbdyn_cmd_available = true;
        // no peer command is available, since we've just put in an mbdyn command
        //buf->peer_cmd_available = false;
        buf->mbdyn_ready_for_cmd = false;

		} break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER: {
		unsigned uN;
		unsigned uNodal;
		bool bRef;
		unsigned uR;
		bool bA;
		bool bL;

		scoped_lock<interprocess_mutex> lock(buf->mutex);

		buf->mbdyn_ready_for_cmd = true;

        pedantic_cout ("[MBDYN] StructExtSharedMemForce::Prepare: notifying that mbdyn is ready to receive" << std::endl);
        buf->cond_mbdyn_ready_for_cmd.notify_all();

        auto time = boost::posix_time::microsec_clock::universal_time();
        pedantic_cout ("[MBDYN] StructExtSharedMemForce::Prepare: time is: " << time << std::endl);

        pedantic_cout ("[MBDYN] StructExtSharedMemForce::Prepare: checking if peer cmd is available, buffer->peer_cmd_available: " << buf->peer_cmd_available << std::endl);

        if (!buf->peer_cmd_available)
        {
            pedantic_cout ("[MBDYN] StructExtSharedMemForce::Prepare: waiting for command from peer" << std::endl);
            buf->cond_peer_cmd_available.wait (lock);
        }

        uNodal = (buf->NodalOrModal & MBC_MODAL_NODAL_MASK);
        bRef = buf->UseReferanceNode;
        uR = (buf->RotationType & MBC_ROT_MASK);
        bL = buf->UseLabels;
        bA = buf->OutputAccelerations;

        uN = buf->NumNodes;

        pedantic_cout ("[MBDYN] StructExtSharedMemForce::Prepare: received data from peer" << std::endl);

        buf->mbdyn_ready_for_cmd = false;
        buf->peer_cmd_available = false; // read command/data, so it is no longer available
        buf->mbdyn_cmd_available = false;

		// check the settings on the remote match local
		bResult = CheckProblemsMatch (uNodal, bRef, uR, bL, bA, uN);

		} break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return bResult;

	// mutex goes out of scope here
}

void
StructExtSharedMemForce::SendToSharedMem(mbdyn::shared_memory_buffer *buf, ExtFileHandlerBase::SendWhen when)
{
    pedantic_cout("StructExtSharedMemForce:in SendToSharedMem" << std::endl);

    // get the data mutex
    scoped_lock<interprocess_mutex> lock(buf->mutex);

    pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: checking if peer is ready to receive" << std::endl);

    if (!buf->peer_ready_for_data)
    {
        pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: waiting for peer be ready to receive" << std::endl);
        buf->cond_peer_ready_for_data.wait (lock);
    }

    pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: beginning to copy data to buffer" << std::endl);

    std::vector<doublereal> buf_t;

	if (pRefNode) {

		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();
		const Vec3& xpRef = pRefNode->GetVCurr();
		const Vec3& wRef = pRefNode->GetWCurr();
		const Vec3& xppRef = pRefNode->GetXPPCurr();
		const Vec3& wpRef = pRefNode->GetWPCurr();

		if (bLabels) {
			pESMH->SetRefLabel(pRefNode->GetLabel());
		}

		switch (uRot) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			//send(outfd, (const char *)RRef.pGetMat(), 9*sizeof(doublereal), 0);
            buf_t = std::vector<doublereal>(RRef.pGetMat(), RRef.pGetMat()+8);
            pESMH->SetRefNodeTheta (buf_t);
			break;

		case MBC_ROT_THETA: {
			Vec3 Theta(RotManip::VecRot(RRef));
			//send(outfd, (const char *)Theta.pGetVec(), 3*sizeof(doublereal), 0);
			buf_t = std::vector<doublereal>(Theta.pGetVec(), Theta.pGetVec()+2);
            pESMH->SetRefNodeTheta (buf_t);
			} break;

		case MBC_ROT_EULER_123: {
			Vec3 E(MatR2EulerAngles123(RRef)*dRaDegr);
			//send(outfd, (const char *)E.pGetVec(), 3*sizeof(doublereal), 0);
			buf_t = std::vector<doublereal>(E.pGetVec(), E.pGetVec()+2);
            pESMH->SetRefNodeTheta (buf_t);
			} break;
		}

		//send(outfd, (const char *)xpRef.pGetVec(), 3*sizeof(doublereal), 0);
		std::vector<doublereal> buf_v (xpRef.pGetVec(), xpRef.pGetVec()+2);
        pESMH->SetRefNodeXP (buf_v);

		if (uRot != MBC_ROT_NONE) {
			//send(outfd, (const char *)wRef.pGetVec(), 3*sizeof(doublereal), 0);
			std::vector<doublereal> buf_w (wRef.pGetVec(), wRef.pGetVec()+2);
            pESMH->SetRefNodeOmega (buf_w);
		}

		if (bOutputAccelerations) {
			//send(outfd, (const char *)xppRef.pGetVec(), 3*sizeof(doublereal), 0);
			std::vector<doublereal> buf_xpp (xppRef.pGetVec(), xppRef.pGetVec()+2);
            pESMH->SetRefNodeXPP (buf_xpp);

			if (uRot != MBC_ROT_NONE) {
				//send(outfd, (const char *)wpRef.pGetVec(), 3*sizeof(doublereal), 0);
				std::vector<doublereal> buf_wp (wpRef.pGetVec(), wpRef.pGetVec()+2);
                pESMH->SetRefNodeOmegaP (buf_wp);
			}
		}

		// output the other node data in the reference frame of the reference node
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
				//uint32_t l = point.pNode->GetLabel();
				//iobuf_labels[i] = l;
				pESMH->SetLabel (i, point.pNode->GetLabel());
			}

			Vec3 xTilde(RRef.MulTV(Dx));
			//memcpy(&iobuf_x[3*i], xTilde.pGetVec(), 3*sizeof(doublereal));
            std::vector<doublereal> buf_x (xTilde.pGetVec(), xTilde.pGetVec()+2);

            pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node x data for node " << i
                  << ": " << buf_x[0] << ", " << buf_x[1] << ", " << buf_x[2] << std::endl);

			pESMH->SetNodeX (i, buf_x);

			switch (uRot) {
			case MBC_ROT_NONE:
				break;

			case MBC_ROT_MAT:
				//memcpy(&iobuf_R[9*i], DR.pGetMat(), 9*sizeof(doublereal));
				buf_t = std::vector<doublereal> (DR.pGetMat(), DR.pGetMat()+8);
				pESMH->SetNodeTheta (i, buf_t);
				break;

			case MBC_ROT_THETA: {
				Vec3 Theta(RotManip::VecRot(DR));
				//memcpy(&iobuf_theta[3*i], Theta.pGetVec(), 3*sizeof(doublereal));
				buf_t = std::vector<doublereal> (Theta.pGetVec(), Theta.pGetVec()+2);
				pESMH->SetNodeTheta (i, buf_t);
				} break;

			case MBC_ROT_EULER_123: {
				Vec3 E(MatR2EulerAngles123(DR)*dRaDegr);
				//memcpy(&iobuf_euler_123[3*i], E.pGetVec(), 3*sizeof(doublereal));
				buf_t = std::vector<doublereal> (E.pGetVec(), E.pGetVec()+2);
				pESMH->SetNodeTheta (i, buf_t);
				} break;
			}

			Vec3 vTilde(RRef.MulTV(Dv));
			//memcpy(&iobuf_xp[3*i], vTilde.pGetVec(), 3*sizeof(doublereal));
			std::vector<doublereal> buf_v (vTilde.pGetVec(), vTilde.pGetVec()+2);
			pESMH->SetNodeXP (i, buf_v);

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

				//memcpy(&iobuf_xpp[3*i], xppTilde.pGetVec(), 3*sizeof(doublereal));
				std::vector<doublereal> buf_xpp (xppTilde.pGetVec(), xppTilde.pGetVec()+2);
				pESMH->SetNodeXPP (i, buf_xpp);

				if (uRot != MBC_ROT_NONE) {
					Vec3 wpTilde(RRef.MulTV(wp) - wpRef - wRef.Cross(w));
					//memcpy(&iobuf_omegap[3*i], wpTilde.pGetVec(), 3*sizeof(doublereal));
					std::vector<doublereal> buf_wp (wpTilde.pGetVec(), wpTilde.pGetVec()+2);
                    pESMH->SetNodeOmegaP (i, buf_wp);
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
				//uint32_t l = point.pNode->GetLabel();
				//iobuf_labels[i] = l;
                pESMH->SetLabel (i, point.pNode->GetLabel());
			}

			//memcpy(&iobuf_x[3*i], x.pGetVec(), 3*sizeof(doublereal));
			std::vector<doublereal> buf_x (x.pGetVec(), x.pGetVec()+3);

			pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node x data for node " << i
                  << ": " << buf_x[0] << ", " << buf_x[1] << ", " << buf_x[2] << std::endl);

			pESMH->SetNodeX (i, buf_x);

			switch (uRot) {
			case MBC_ROT_NONE:
				break;

			case MBC_ROT_MAT:
				//memcpy(&iobuf_R[9*i], R.pGetMat(), 9*sizeof(doublereal));
				buf_t = std::vector<doublereal> (R.pGetMat(), R.pGetMat()+9);
				pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node orientation data for node " << i
                  << ": " << buf_t[0] << ", " << buf_t[1] << ", " << buf_t[2] << ", " << buf_t[3] << ", " << buf_t[4] << ", " << buf_t[5] << ", " << buf_t[6] << ", " << buf_t[7] << ", " << buf_t[8]  << std::endl);
				pESMH->SetNodeTheta (i, buf_t);
				break;

			case MBC_ROT_THETA: {
				Vec3 Theta(RotManip::VecRot(R));
				//memcpy(&iobuf_theta[3*i], Theta.pGetVec(), 3*sizeof(doublereal));
				buf_t = std::vector<doublereal> (Theta.pGetVec(), Theta.pGetVec()+3);
				pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node orientation data for node " << i
                  << ": " << buf_t[0] << ", " << buf_t[1] << ", " << buf_t[2] << std::endl);

				pESMH->SetNodeTheta (i, buf_t);
				} break;

			case MBC_ROT_EULER_123: {
				Vec3 E(MatR2EulerAngles123(R)*dRaDegr);
				//memcpy(&iobuf_euler_123[3*i], E.pGetVec(), 3*sizeof(doublereal));
				buf_t = std::vector<doublereal> (E.pGetVec(), E.pGetVec()+3);
				pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node orientation data for node " << i
                  << ": " << buf_t[0] << ", " << buf_t[1] << ", " << buf_t[2] << std::endl);
				pESMH->SetNodeTheta (i, buf_t);
				} break;
			}

			//memcpy(&iobuf_xp[3*i], v.pGetVec(), 3*sizeof(doublereal));
			std::vector<doublereal> buf_v (v.pGetVec(), v.pGetVec()+3);
			pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node xp data for node " << i
                  << ": " << buf_v[0] << ", " << buf_v[1] << ", " << buf_v[2] << std::endl);
			pESMH->SetNodeXP (i, buf_v);

			if (uRot != MBC_ROT_NONE) {
				//memcpy(&iobuf_omega[3*i], w.pGetVec(), 3*sizeof(doublereal));
				std::vector<doublereal> buf_w (w.pGetVec(), w.pGetVec()+3);
				pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node omega data for node " << i
                  << ": " << buf_w[0] << ", " << buf_w[1] << ", " << buf_w[2] << std::endl);
				pESMH->SetNodeOmega (i, buf_w);
			}

			if (bOutputAccelerations) {
				const Vec3& wp = point.pNode->GetWPCurr();
				Vec3 xpp = point.pNode->GetXPPCurr() + wp.Cross(f) + w.Cross(wCrossf);

                std::vector<doublereal> buf_xpp (xpp.pGetVec(), xpp.pGetVec()+3);
                pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node xpp data for node " << i
                  << ": " << buf_xpp[0] << ", " << buf_xpp[1] << ", " << buf_xpp[2] << std::endl);
				pESMH->SetNodeXPP (i, buf_xpp);

				//memcpy(&iobuf_xpp[3*i], xpp.pGetVec(), 3*sizeof(doublereal));
				if (uRot != MBC_ROT_NONE) {
					//memcpy(&iobuf_omegap[3*i], wp.pGetVec(), 3*sizeof(doublereal));
					std::vector<doublereal> buf_wp (wp.pGetVec(), wp.pGetVec()+3);
					pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Sending node omegaP data for node " << i
                      << ": " << buf_wp[0] << ", " << buf_wp[1] << ", " << buf_wp[2] << std::endl);
                    pESMH->SetNodeOmegaP (i, buf_wp);
				}

			}
		}
	}

    pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: put node data for " << m_Points.size() << " nodes in buffer" << std::endl);

    pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: notifying new command is available" << std::endl);

    buf->mbdyn_data_available = true;
    buf->cond_mbdyn_data_available.notify_all();

    pedantic_cout("[MBDYN] StructExtSharedMemForce::SendToSharedMem: buf->mbdyn_cmd_available is " << buf->mbdyn_cmd_available << std::endl);

	// mutex released here
}

void
StructExtSharedMemForce::RecvFromSharedMem(mbdyn::shared_memory_buffer *buf)
{
    pedantic_cout("StructExtSharedMemForce:in RecvFromSharedMem" << std::endl);

    // get the data mutex
    scoped_lock<interprocess_mutex> lock(buf->mutex);

    buf->mbdyn_ready_for_data = true;

    pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem: notifying that mbdyn is ready to receive" << std::endl);
    buf->cond_mbdyn_ready_for_data.notify_all();

    pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem: checking if peer cmd is available" << std::endl);
    if (!buf->peer_data_available)
    {
        pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem: waiting for data from peer" << std::endl);
        buf->cond_peer_data_available.wait (lock);
    }

    pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem: data from peer available, beginning copy" << std::endl);

	if (pRefNode) {

		if (bLabels) {
			unsigned l = pESMH->GetRefLabel();
			if (l != pRefNode->GetLabel()) {
				silent_cerr("StructExtSharedMemForce(" << GetLabel() << "): "
					"invalid reference node label "
					"(wanted " << pRefNode->GetLabel() << ", got " << l << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		std::vector<doublereal> f = pESMH->GetRefNodeF ();
		F0 = Vec3(f[0], f[1], f[2]);

		if (uRot != MBC_ROT_NONE) {
            std::vector<doublereal> m = pESMH->GetRefNodeM ();
			M0 = Vec3(m[0], m[1], m[2]);
		}

	}

	if (!bSorted) {
		ASSERT(bLabels);

		done.resize(m_Points.size());
		fill(done.begin(), done.end(), false);

		unsigned cnt;
		for (cnt = 0; cnt < m_Points.size(); cnt++) {
			PointData& point = m_Points[cnt];

			unsigned l = pESMH->GetLabel(cnt);

			// check all the points to see if this label exists
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

			std::vector<doublereal> f = pESMH->GetNodeF (i);
			point.F = Vec3(f[0], f[1], f[2]);
			if (uRot != MBC_ROT_NONE) {
                //point.M = &iobuf_m[3*i];
                std::vector<doublereal> m = pESMH->GetNodeM (i);
				point.M = Vec3 (m[0], m[1], m[2]);
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
				unsigned l = pESMH->GetLabel(i);
				if (point.pNode->GetLabel() != l) {
					silent_cerr("StructExtForce"
						"(" << GetLabel() << "): "
						"invalid " << i << "-th label " << l
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			std::vector<doublereal> f = pESMH->GetNodeF (i);
			point.F = Vec3(f[0], f[1], f[2]);

			pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Forces received for node " << i
                  << " were " << f[0] << ", " << f[1] << ", " << f[2] << std::endl);

			if (uRot != MBC_ROT_NONE) {
                std::vector<doublereal> m = pESMH->GetNodeM (i);
				point.M = Vec3 (m[0], m[1], m[2]);

				pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem. Moments received for node " << i
                  << " were " << m[0] << ", " << m[1] << ", " << m[2] << std::endl);
			}
		}
	}

    pedantic_cout ("[MBDYN] StructExtSharedMemForce::RecvFromSharedMem:received data from peer" << std::endl);

    buf->mbdyn_ready_for_data = false;
    // no peer command is available (we've read it)
    buf->peer_data_available = false;
    // no data is available from mbdyn
    buf->mbdyn_data_available = false;

	// mutex released here
}

void
StructExtSharedMemForce::SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when)
{
    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
StructExtSharedMemForce::SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
StructExtSharedMemForce::RecvFromStream(std::istream& inf)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
StructExtSharedMemForce::RecvFromFileDes(int infd)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
