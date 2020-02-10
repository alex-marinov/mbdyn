/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

//#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <cstring>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "sharedmem.h"
#include "mbc.h"
#include "mbcxx.h"
#include "mbcxxshared.h"

using namespace boost::interprocess;
using namespace mbdyn;

static const char *
mbc_cmd2str(uint8_t cmd)
{
	switch (cmd) {
	case ES_REGULAR_DATA: return "REGULAR_DATA";
	case ES_GOTO_NEXT_STEP: return "GOTO_NEXT_STEP";
	case ES_ABORT: return "ABORT";
	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP: return "REGULAR_DATA_AND_GOTO_NEXT_STEP";
	case ES_NEGOTIATION: return "NEGOTIATION";
	case ES_OK: return "OK";
	default:
		break;
	}

	return "UNKNOWN";
}

int
MBCSharedMemBase::CheckCmd(const uint8_t cmd) const
{
	switch (cmd) {
	case ES_REGULAR_DATA:
	case ES_GOTO_NEXT_STEP:
	case ES_ABORT:
	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP:
	case ES_NEGOTIATION:
	case ES_OK:
    case ES_LAST:
		return 0;
	}

	fprintf(stderr, "unknown cmd (%lu) from peer\n", (unsigned long)cmd);
	return -1;
}

MBCType
MBCSharedMemBase::GetRot(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return _Rot;
}

bool
MBCSharedMemBase::bRefNode(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return _bUseReferanceNode;
}

MBCType
MBCSharedMemBase::GetRefNodeRot(void) const {
	assert(GetStatus() >= INITIALIZED);
	return _RefNodeRot;
}

bool
MBCSharedMemBase::bAccelerations(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return _bOutputAccelerations;
}

bool
MBCSharedMemBase::bLabels(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return _bUseLabels;
}

void
MBCSharedMemBase::SetTimeout(int t)
{
	_timeout = t;
}

void
MBCSharedMemBase::SetVerbose(bool bv)
{
	_bVerbose = bv;
}

void
MBCSharedMemBase::SetDataAndNext(bool bd)
{
	_bDataAndNext = bd;
}

bool
MBCSharedMemBase::bVerbose(void) const
{
	return _bVerbose;
}

bool
MBCSharedMemBase::bDataAndNext(void) const
{
	return _bDataAndNext;
}

MBCSharedMemBase::Status
MBCSharedMemBase::GetStatus(void) const
{
	return m_status;
}

void
MBCSharedMemBase::SetStatus(MBCSharedMemBase::Status s)
{
	m_status = s;
}

MBCSharedMemBase::MBCSharedMemBase(void)
: m_status(NOT_READY)
{
}

MBCSharedMemBase::~MBCSharedMemBase(void)
{
}

int
MBCSharedMemBase::Init()
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::Init"  << std::endl; }

    int rc = 0;

    if (GetStatus() != INITIALIZED) return -1;

	SetStatus(SHARED_MEMORY_READY);

	return rc;
}


//uint8_t
//MBCSharedMemBase::GetCmd(void) const
//{
//	return cmd;
//}

uint32_t
MBCSharedMemBase::GetRefNodeKinematicsLabel(void) const
{
	assert(GetStatus() == READY);
	return refnode_label;
}

//uint32_t
//MBCSharedMemBase::KinematicsLabel(void) const
//{
//	assert(GetStatus() == READY);
//	assert(bLabels());
//	return refnode_label;
//}

std::vector<double>
MBCSharedMemBase::GetRefNodeX(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());

    int start = 0;

    std::vector<double> X = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        X.at(i) = refnode_kinematics_data.at(start+i);
    }

    return X;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeR(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_MAT);

    unsigned start = 9;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = refnode_kinematics_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeTheta(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_THETA);

    unsigned start = 9;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = refnode_kinematics_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeEuler123(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_EULER_123);

    unsigned start = 9;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = refnode_kinematics_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeXP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());

    int start = 3;

    std::vector<double> XP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XP.at(i) = refnode_kinematics_data.at(start+i);
    }

    return XP;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeOmega(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);

    int start = 9 + nrotvalspernode;

    std::vector<double> omega = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omega.at(i) = refnode_kinematics_data.at(start+i);
    }

    return omega;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeXPP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(bAccelerations());

    int start = 6;

    std::vector<double> XPP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XPP.at(i) = refnode_kinematics_data.at(start+i);
    }

    return XPP;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeOmegaP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);
	assert(bAccelerations());

    int start = 12 + nrotvalspernode;

    std::vector<double> omegaP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omegaP.at(i) = refnode_kinematics_data.at(start+i);
    }

    return omegaP;
}


double
MBCSharedMemBase::GetRefNodeX(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;

    int start = 0;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeR(uint8_t ir, uint8_t ic) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_MAT);
	if (ir < 1 || ir > 3 || ic < 1 || ic > 3) throw;

    unsigned start = 9;

    return refnode_kinematics_data.at(start+(3*(ic - 1) + ir - 1));
}

double
MBCSharedMemBase::GetRefNodeTheta(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_THETA);
	if (idx < 1 || idx > 3) throw;

    unsigned start = 9;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeEuler123(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MBC_ROT_EULER_123);
	if (idx < 1 || idx > 3) throw;

    unsigned start = 9;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeXP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;

    int start = 3;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeOmega(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);
	if (idx < 1 || idx > 3) throw;

    int start = 9 + nrotvalspernode;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeXPP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(bAccelerations());
	if (idx < 1 || idx > 3) throw;

    int start = 6;

    return refnode_kinematics_data.at(start+(idx-1));
}

double
MBCSharedMemBase::GetRefNodeOmegaP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);
	assert(bAccelerations());
	if (idx < 1 || idx > 3) throw;

    int start = 12 + nrotvalspernode;

    return refnode_kinematics_data.at(start+(idx-1));
}



//const double&
//MBCSharedMemBase::X(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_X(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::R(uint8_t ir, uint8_t ic) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() == MAT);
//	if (ir < 1 || ir > 3 || ic < 1 || ic > 3) throw;
//	return (MBC_R_R(GetRefNodePtr()))[3*(ic - 1) + ir - 1];
//}
//
//const double&
//MBCSharedMemBase::Theta(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() == THETA);
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_THETA(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::Euler123(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() == EULER_123);
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_EULER_123(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::XP(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_XP(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::Omega(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() != NONE);
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_OMEGA(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::XPP(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(bAccelerations());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_XPP(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::OmegaP(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() != NONE);
//	assert(bAccelerations());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_OMEGAP(GetRefNodePtr()))[idx - 1];
//}

uint32_t
MBCSharedMemBase::GetRefNodeDynamicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return refnode_label;
}

//const uint32_t&
//MBCSharedMemBase::DynamicsLabel(void) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	return MBC_R_D_LABEL(GetRefNodePtr());
//}
//
//uint32_t&
//MBCSharedMemBase::DynamicsLabel(void)
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	return MBC_R_D_LABEL(GetRefNodePtr());
//}

std::vector<double>
MBCSharedMemBase::GetRefNodeF(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());

    int start = 0;

    std::vector<double> forces = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        forces.at(i) = refnode_dynamics_data.at(start+i);
    }

    return forces;
}

std::vector<double>
MBCSharedMemBase::GetRefNodeM(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);

    int start = 3;

    std::vector<double> moments = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        moments.at(i) = refnode_dynamics_data.at(start+i);
    }

    return moments;
}


void
MBCSharedMemBase::SetRefNodeF(uint8_t idx, const double &force)
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;

    int start = 0;

    refnode_dynamics_data.at(start+(idx-1)) = force;
}

void
MBCSharedMemBase::SetRefNodeM(uint8_t idx, const double &moment)
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != MBC_ROT_NONE);

    int start = 3;

    refnode_dynamics_data.at(start+(idx-1)) = moment;
}


//
//const double&
//MBCSharedMemBase::F(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_F(GetRefNodePtr()))[idx - 1];
//}
//
//double&
//MBCSharedMemBase::F(uint8_t idx)
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_F(GetRefNodePtr()))[idx - 1];
//}
//
//const double&
//MBCSharedMemBase::M(uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() != NONE);
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_M(GetRefNodePtr()))[idx - 1];
//}
//
//double&
//MBCSharedMemBase::M(uint8_t idx)
//{
//	assert(GetStatus() == READY);
//	assert(bRefNode());
//	assert(GetRefNodeRot() != NONE);
//	if (idx < 1 || idx > 3) throw;
//	return (MBC_R_M(GetRefNodePtr()))[idx - 1];
//}
//
//mbc_t *
//MBCSharedMemNodal::GetBasePtr(void) const
//{
//	return (mbc_t*)&mbc;
//}
//
//mbc_refnode_stub_t *
//MBCSharedMemNodal::GetRefNodePtr(void) const
//{
//	return (mbc_refnode_stub_t *)&mbc;
//}

MBCSharedMemNodal::MBCSharedMemNodal(void)
{
	//memset(&mbc, 0, sizeof(mbc));
}

MBCSharedMemNodal::MBCSharedMemNodal(MBCType refnode_rot, unsigned nodes,
       	bool uselabels, MBCType rot, bool accels, const std::string &shmName)
{
	if (Initialize(refnode_rot, nodes, uselabels, rot, accels, shmName)) {
		throw;
	}
}

MBCSharedMemNodal::~MBCSharedMemNodal(void)
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::~MBCSharedMemNodal"  << std::endl; }
	Close();
}

MBCSharedMemBase::Type
MBCSharedMemNodal::GetType(void) const
{
	return NODAL;
}

int
MBCSharedMemNodal::Initialize( MBCType refnode_rot,
                               unsigned nodes,
                               bool uselabels,
                               MBCType rot,
                               bool accels,
                               const std::string &shmName )
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::~Initialize"  << std::endl; }

	if (GetStatus() != NOT_READY) return -1;

    // Open the managed segment
    if (_bVerbose)
    {
        std::cout << "Opening shared memory with name: '" << shmName << "'" << std::endl;
    }

    shm = managed_shared_memory(open_only, shmName.c_str ());

    // there should only be one buffer in this shared memory region
    std::pair<shared_memory_buffer*, int> temp = shm.find<shared_memory_buffer>(unique_instance);

    buffer = temp.first;

    // get the buffers in shared memory
    S_X_data = shm.find<shmVector_double>("X_data_vector").first;
    S_XP_data = shm.find<shmVector_double>("XP_data_vector").first;
    S_XPP_data = shm.find<shmVector_double>("XPP_data_vector").first;
    S_theta_data = shm.find<shmVector_double>("theta_data_vector").first;
    S_omega_data = shm.find<shmVector_double>("omega_data_vector").first;
    S_omegaP_data = shm.find<shmVector_double>("omegaP_data_vector").first;
    S_force_data = shm.find<shmVector_double>("force_data_vector").first;
    S_moment_data = shm.find<shmVector_double>("moment_data_vector").first;
    S_refnode_kinematics_data = shm.find<shmVector_double>("refnode_kinematics_data_vector").first;
    S_refnode_dynamics_data = shm.find<shmVector_double>("refnode_dynamics_data_vector").first;
    S_labels = shm.find<shmVector_uint32_t>("label_vector").first;

//	int rc = mbc_nodal_init(&mbc, refnode_rot, nodes, labels, unsigned(rot) | MBC_U_ROT_2_REF_NODE_ROT(unsigned(refnode_rot)), accels);

    nnodes = nodes;
    nrotvalspernode = 0;
    nrefnoderotvalspernode = 0;

    // n node labels
    labels.resize (nnodes);

    _bUseReferanceNode = false;
    if (refnode_rot) {
		_bUseReferanceNode = true;

		_RefNodeRot = refnode_rot;

        switch (_RefNodeRot) {
        case MBC_ROT_NONE:
            fprintf(stderr, "rotation must be defined for reference node\n");
            return -1;

        case MBC_ROT_MAT:
            // full orientation matrix for orientations
            nrefnoderotvalspernode = 9;
            break;

        case MBC_ROT_THETA:
        case MBC_ROT_EULER_123:
            // 3 orientation values
            nrefnoderotvalspernode = 3;
            break;
        default:
            fprintf(stderr, "unknown orientation parametrization 0x%lx in flags\n", (unsigned long)rot);
            return -1;
        }
	}

	if ((_bUseReferanceNode == false) && (nnodes == 0)) {
		fprintf(stderr, "need at least 1 node or reference node data\n");
		return -1;
	}

    _bOutputAccelerations = false;
	if (accels) {
		_bOutputAccelerations = true;
	}

	_bUseLabels = false;
	if (uselabels) {
		_bUseLabels = true;
	}

	// refnode motion in 6dof
    refnode_kinematics_data.resize (3 + 3 + 3 + nrefnoderotvalspernode + 3 + 3);

    // refnode motion in 6dof
    refnode_dynamics_data.resize (3 + 3 + 3 + nrefnoderotvalspernode + 3 + 3);

    _Rot = rot;

	if (nnodes > 0) {

		switch (_Rot) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			nrotvalspernode = 9;
			break;

		case MBC_ROT_THETA:
		case MBC_ROT_EULER_123:
			nrotvalspernode = 3;
			break;
		}

		if (_bOutputAccelerations) {
			XPP_data.resize (3 * nnodes);
			if (_Rot != MBC_ROT_NONE) {
				omegaP_data.resize (3 * nnodes);
			}
		}

        // position, velocity and acceleration in 6dof

        // forces and moments
        force_data.resize (3 * nnodes);
		if (_Rot != MBC_ROT_NONE) {
			moment_data.resize (3 * nnodes);
		}

		if (_bUseLabels) {
			labels.resize (3 * nnodes);
		}

		X_data.resize (3 * nnodes);
		XP_data.resize (3 * nnodes);
		if (_Rot != MBC_ROT_NONE) {
            theta_data.resize (nrotvalspernode * nnodes);
            omega_data.resize (3 * nnodes);
		}

	}

    if (_bVerbose){ std::cout << "Completed initialisation" << std::endl; }

	SetStatus(INITIALIZED);

	return 0;
}

int
MBCSharedMemNodal::Negotiate(void)
{

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::Negotiate"  << std::endl; }

    uint8_t u = 0;

	if (!_bUseReferanceNode && nnodes == 0) {
		fprintf(stderr, "need at least 1 node or reference node data\n");
		return -1;
	}

    PutCmd(uint8_t(ES_NEGOTIATION));

    {
        // lock the mutex so we can safely access the data
        scoped_lock<interprocess_mutex> lock(buffer->mutex);

        if (_bVerbose) { std::cout << "MBCSharedMemNodal::Negotiate: checking if mbdyn is ready for command" << std::endl; }

        if (!buffer->mbdyn_ready_for_cmd)
        {
            if (_bVerbose) { std::cout << "MBCSharedMemNodal::Negotiate: waiting for mbdyn to be ready for command" << std::endl; }
            buffer->cond_mbdyn_ready_for_cmd.wait (lock);
        }

        buffer->NumNodes = nnodes;
        buffer->NodalOrModal = MBC_NODAL;
        buffer->UseReferanceNode = _bUseReferanceNode;
        buffer->RotationType = _Rot;
        buffer->RefRotationType = _RefNodeRot;
        buffer->UseLabels = _bUseLabels;
        buffer->OutputAccelerations = _bOutputAccelerations;
        buffer->NumNodes = nnodes;

        if (_bVerbose) { std::cout << "MBCSharedMemNodal::Negotiate: peer put problem data in shared memory" << std::endl; }

        buffer->peer_cmd_available = true;
        buffer->peer_ready_for_cmd = false;
        //buffer->mbdyn_cmd_available = false;

        // notify mbdyn it can check for a command
        buffer->cond_peer_cmd_available.notify_all();

        if (_bVerbose) { std::cout << "MBCSharedMemNodal::Negotiate: peer notified that command is available, buffer->peer_cmd_available: " << buffer->peer_cmd_available << std::endl; }

        if (_bVerbose) { std::cout << "MBCSharedMemNodal::Negotiate: time is: " << boost::posix_time::microsec_clock::universal_time() << std::endl; }

    } // buffer mutex released here

    // get the response when it is ready
    GetCmd(u);

	if (_bVerbose) {
		fprintf(stdout, "cmd from mbdyn: %lu (%s)\n",
			(unsigned long)u, mbc_cmd2str(u));
	}

	switch (u) {
	case ES_ABORT:
		fprintf(stdout, "got ABORT from mbdyn\n");
		return -1;

    case ES_LAST:
        // mbdyn simulation is finished
        SetStatus (FINISHED);
        return -1;

	case ES_OK:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%lu from mbdyn\n", (unsigned long)u);
		return -1;
	}

	SetStatus (READY);

	if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::Negotiate at end"  << std::endl; }

	return 0;
}

int
MBCSharedMemBase::PutCmd(const uint8_t u) const
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutCmd"  << std::endl; }

	if (CheckCmd(u)) {
		return -1;
	}

	if (_bVerbose) {
		fprintf(stdout, "cmd to mbdyn: %lu (%s)\n",
			(unsigned long)u, mbc_cmd2str(u));
	}

	// lock the mutex so we can safely access the data
    scoped_lock<interprocess_mutex> lock(buffer->mutex);

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::PutCmd checking if mbdyn is ready to receive command" << std::endl; }

    if (!buffer->mbdyn_ready_for_cmd)
    {
        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::PutCmd waiting for mbdyn to be ready to receive command" << std::endl; }
        buffer->cond_mbdyn_ready_for_cmd.wait (lock);
    }

	buffer->cmd = u;

	if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::PutCmd peer put command " << unsigned (u) << " in shared memory" << std::endl; }

    buffer->peer_cmd_available = true;
    buffer->mbdyn_cmd_available = false;
    buffer->peer_ready_for_cmd = false;

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::PutCmd peer notifying command " << unsigned(buffer->cmd) << " is available" << std::endl; }
    buffer->cond_peer_cmd_available.notify_all();

	return 0;

	// buffer mutex released here
}

int
MBCSharedMemBase::GetCmd(uint8_t &u) const
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetCmd"  << std::endl; }

	// lock the mutex so we can safely access the data
    scoped_lock<interprocess_mutex> lock(buffer->mutex);

    if (_bVerbose) { std::cout << "MBCSharedMemBase::GetCmd: buffer->peer_cmd_available: " << buffer->peer_cmd_available << std::endl; }

    buffer->peer_ready_for_cmd = true;

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::GetCmd notifying that peer is ready to receive" << std::endl; }

    buffer->cond_peer_ready_for_cmd.notify_all();

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::GetCmd checking if mbdyn command is available" << std::endl; }
    if (!buffer->mbdyn_cmd_available)
    {
        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::GetCmd waiting for mbdyn command to become available" << std::endl; }
        buffer->cond_mbdyn_cmd_available.wait (lock);
    }

	u = buffer->cmd;

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemBase::GetCmd received cmd " << unsigned(u) << " from mbdyn" << std::endl; }

    buffer->peer_ready_for_cmd = false;
    buffer->mbdyn_cmd_available = false;
    buffer->peer_cmd_available = false;

    if (CheckCmd(u)) {
		return -1;
	}

//	if (u == ES_LAST)
//    {
//        // mbdyn simulation is finished
//        SetStatus (FINISHED);
//    }

	return 0;

	// buffer mutex released here
}

int
MBCSharedMemNodal::PutForces(bool bConverged)
{

    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutForces"  << std::endl; }

	if (GetStatus() != READY) return -1;

	cmd = 0;

    if (bConverged) {
		if (_bDataAndNext) {
			cmd = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;

		} else {
			cmd = ES_GOTO_NEXT_STEP;
		}

	} else {
		cmd = ES_REGULAR_DATA;
	}

	if (PutCmd(cmd)) {
		return -1;
	}

	if (cmd != ES_GOTO_NEXT_STEP) {
        // send new forces, copying the data from the local buffers

        // lock the mutex so we can safely access the data
        scoped_lock<interprocess_mutex> lock(buffer->mutex);

        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutForces: checking if mbdyn ready to receive force data: buf->mbdyn_ready_for_cmd is " << buffer->mbdyn_ready_for_cmd << std::endl; }

        if (!buffer->mbdyn_ready_for_data)
        {
            if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutForces: waiting for mbdyn to be ready to receive force data"  << std::endl; }
            buffer->cond_mbdyn_ready_for_data.wait (lock);
        }

        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutForces: writing force data to shared memory buffer"  << std::endl; }

		/* reference node */
		if (_bUseReferanceNode) {
            S_refnode_dynamics_data->assign (refnode_dynamics_data.begin(), refnode_dynamics_data.end());
		}

		/* nodal */
		if (nnodes > 0) {
            S_force_data->assign(force_data.begin(), force_data.end());

            if (_bVerbose)
            {
                int j = 0;
                std::cout << "[PEER] MBCSharedMemNodal::PutForces: S_force_data now contains:"  << std::endl;
                for (int i = 0; i < S_force_data->size (); i++)
                {
                    std::cout << S_force_data->at(i);

                    if (j == 2)
                    {
                        std::cout << std::endl;
                        j = 0;
                    }

                    j++;
                }
            }

            S_moment_data->assign(moment_data.begin(), moment_data.end());

            if (_bVerbose)
            {
                int j = 0;
                std::cout << "[PEER] MBCSharedMemNodal::PutForces: S_moment_data now contains:"  << std::endl;
                for (int i = 0; i < S_moment_data->size (); i++)
                {
                    std::cout << S_moment_data->at(i);

                    if (j == 2)
                    {
                        std::cout << std::endl;
                        j = 0;
                    }

                    j++;
                }
            }
		}

		buffer->peer_data_available = true;
		buffer->mbdyn_data_available = false;
		buffer->peer_ready_for_data = false;
		buffer->cond_peer_data_available.notify_one ();

		if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::PutForces: notifying that peer cmd is available"  << std::endl; }
	}

	return 0;
}

int
MBCSharedMemNodal::GetMotion(void)
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion"  << std::endl; }

	if (GetStatus() != READY) return -1;

	cmd = 0;

	if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion about to call GetCmd"  << std::endl; }
    if (GetCmd(cmd)) {
		return -1;
	}

	if (_bVerbose) {
		fprintf(stdout, "cmd from mbdyn: %lu (%s)\n",
			(unsigned long)cmd, mbc_cmd2str(cmd));
	}

	if (cmd == ES_ABORT) {
		fprintf(stdout, "got ABORT from mbdyn\n");
		return -1;
	}

	if (cmd == ES_LAST) {
		fprintf(stdout, "got ES_LAST from mbdyn (sim finished)\n");
		SetStatus (FINISHED);
		return -1;
	}

	if (cmd != ES_GOTO_NEXT_STEP) {
		// copy over the motion data to local buffers

        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion cmd != ES_GOTO_NEXT_STEP"  << std::endl; }

		// lock the mutex so we can safely access the data
        scoped_lock<interprocess_mutex> lock(buffer->mutex);

        buffer->peer_ready_for_data = true;

        buffer->cond_peer_ready_for_data.notify_one ();

        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion checking if mbdyn data is available: buffer->mbdyn_cmd_available is" << buffer->mbdyn_cmd_available << std::endl; }

        if (!buffer->mbdyn_data_available)
        {
            if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion waiting for mbdyn data to be available"  << std::endl; }
            buffer->cond_mbdyn_data_available.wait (lock);
        }

        if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::GetMotion copying data from shared memory"  << std::endl; }

		if (_bUseReferanceNode) {
            refnode_kinematics_data.assign (S_refnode_dynamics_data->begin(), S_refnode_dynamics_data->end());
            // wait for the server to release the semaphore
            refnode_label = buffer->refnode_label;
		}

		if (nnodes > 0) {
			X_data.assign (S_X_data->begin(), S_X_data->end());
			XP_data.assign (S_XP_data->begin(), S_XP_data->end());
			XPP_data.assign (S_XPP_data->begin(), S_XPP_data->end());
			theta_data.assign (S_theta_data->begin(), S_theta_data->end());
			omega_data.assign (S_omega_data->begin(), S_omega_data->end());
			omegaP_data.assign (S_omegaP_data->begin(), S_omegaP_data->end());
			labels.assign (S_labels->begin(), S_labels->end());
		}

		// we read the command, it is no longer available
		buffer->mbdyn_data_available = false;
		// peer is no longer ready for data (we've just read it)
		buffer->peer_ready_for_data = false;
		// we've just read the motion data, there will not be any peer data available yet
		buffer->peer_data_available = false;

	}

	return 0;
}

int
MBCSharedMemNodal::Close(void) const
{
    if (_bVerbose) { std::cout << "[PEER] MBCSharedMemNodal::Close"  << std::endl; }
	int rc = 0;
	if (GetStatus() == READY) {
		//rc = mbc_nodal_destroy(&mbc);
		PutCmd(ES_ABORT);
		const_cast<MBCSharedMemNodal *>(this)->SetStatus(FINISHED);
	}
	return rc;
}

//uint32_t
//MBCSharedMemNodal::KinematicsLabel(void) const
//{
//	return MBCSharedMemBase::KinematicsLabel();
//}
//
//const double&
//MBCSharedMemNodal::X(uint8_t idx) const
//{
//	return MBCSharedMemBase::X(idx);
//}
//
//const double&
//MBCSharedMemNodal::R(uint8_t ir, uint8_t ic) const
//{
//	return MBCSharedMemBase::R(ir, ic);
//}
//
//const double&
//MBCSharedMemNodal::Theta(uint8_t idx) const
//{
//	return MBCSharedMemBase::Theta(idx);
//}
//
//const double&
//MBCSharedMemNodal::Euler123(uint8_t idx) const
//{
//	return MBCSharedMemBase::Euler123(idx);
//}
//
//const double&
//MBCSharedMemNodal::XP(uint8_t idx) const
//{
//	return MBCSharedMemBase::XP(idx);
//}
//
//const double&
//MBCSharedMemNodal::Omega(uint8_t idx) const
//{
//	return MBCSharedMemBase::Omega(idx);
//}
//
//const double&
//MBCSharedMemNodal::XPP(uint8_t idx) const
//{
//	return MBCSharedMemBase::XPP(idx);
//}
//
//const double&
//MBCSharedMemNodal::OmegaP(uint8_t idx) const
//{
//	return MBCSharedMemBase::OmegaP(idx);
//}
//
//const uint32_t&
//MBCSharedMemNodal::DynamicsLabel(void) const
//{
//	return MBCSharedMemBase::DynamicsLabel();
//}
//
//uint32_t&
//MBCSharedMemNodal::DynamicsLabel(void)
//{
//	return MBCSharedMemBase::DynamicsLabel();
//}
//
//const double&
//MBCSharedMemNodal::F(uint8_t idx) const
//{
//	return MBCSharedMemBase::F(idx);
//}
//
//double&
//MBCSharedMemNodal::F(uint8_t idx)
//{
//	return MBCSharedMemBase::F(idx);
//}
//
//const double&
//MBCSharedMemNodal::M(uint8_t idx) const
//{
//	return MBCSharedMemBase::M(idx);
//}
//
//double&
//MBCSharedMemNodal::M(uint8_t idx)
//{
//	return MBCSharedMemBase::M(idx);
//}

uint32_t
MBCSharedMemNodal::GetNodes(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return nnodes;
}

const std::vector<uint32_t> *
MBCSharedMemNodal::GetKinematicsLabels(void) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	return &labels;
}

const std::vector<double> *
MBCSharedMemNodal::GetX(void) const
{
	assert(GetStatus() == READY);
	return &X_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetR(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_MAT);
	return &theta_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetTheta(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_THETA);
	return &theta_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetEuler123(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_EULER_123);
	return &theta_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetXP(void) const
{
	assert(GetStatus() == READY);
	return &XP_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetOmega(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	return &omega_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetXPP(void) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	return &XPP_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetOmegaP(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	assert(bAccelerations());
	return &omegaP_data;
}
//
//const double&
//MBCSharedMemNodal::X(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_X(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::R(uint32_t n, uint8_t ir, uint8_t ic) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() == MAT);
//	if (ir < 1 || ir > 3 || ic < 1 || ic > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_R(&mbc))[9*(n - 1) + 3*(ic - 1) + ir - 1];
//}
//
//const double&
//MBCSharedMemNodal::Theta(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() == THETA);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_THETA(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::Euler123(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() == EULER_123);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_EULER_123(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::XP(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_XP(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::Omega(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() != NONE);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_OMEGA(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::XPP(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(bAccelerations());
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_XPP(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::OmegaP(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() != NONE);
//	assert(bAccelerations());
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_OMEGAP(&mbc))[3*(n - 1) + idx - 1];
//}

uint32_t
MBCSharedMemNodal::GetKinematicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return labels[n - 1];
}

uint32_t
MBCSharedMemNodal::GetDynamicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return labels[n - 1];
}

//uint32_t
//MBCSharedMemNodal::KinematicsLabel(uint32_t n) const
//{
//	assert(GetStatus() == READY);
//	assert(bLabels());
//	if (n < 1 || n > GetNodes()) throw;
//	return labels[n - 1];
//}

std::vector<double>
MBCSharedMemNodal::GetX(uint32_t n) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> X = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        X.at(i) = X_data.at(start+i);
    }

    return X;
}

std::vector<double>
MBCSharedMemNodal::GetR(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_MAT);
	if (n < 1 || n > GetNodes()) throw;

    unsigned start = (n-1) * nrotvalspernode;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = theta_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemNodal::GetTheta(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_THETA);
	if (n < 1 || n > GetNodes()) throw;

    unsigned start = (n-1) * nrotvalspernode;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = theta_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemNodal::GetEuler123(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_EULER_123);
	if (n < 1 || n > GetNodes()) throw;

    unsigned start = (n-1) * nrotvalspernode;

    std::vector<double> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = theta_data.at(start+i);
    }

    return theta;
}

std::vector<double>
MBCSharedMemNodal::GetXP(uint32_t n) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> XP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XP.at(i) = XP_data.at(start+i);
    }

    return XP;
}

std::vector<double>
MBCSharedMemNodal::GetOmega(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> omega = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omega.at(i) = omega_data.at(start+i);
    }

    return omega;
}

std::vector<double>
MBCSharedMemNodal::GetXPP(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> XPP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XPP.at(i) = XPP_data.at(start+i);
    }

    return XPP;

}

std::vector<double>
MBCSharedMemNodal::GetOmegaP(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	assert(bAccelerations());
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> omegaP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omegaP.at(i) = omegaP_data.at(start+i);
    }

    return omegaP;
}



double
MBCSharedMemNodal::GetX(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    return X_data.at(start+(idx-1));
}

double
MBCSharedMemNodal::GetR(uint32_t n, uint8_t ir, uint8_t ic) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_MAT);
	if (n < 1 || n > GetNodes()) throw;
    if (ir < 1 || ir > 3 || ic < 1 || ic > 3) throw;

    unsigned start = (n-1) * nrotvalspernode;

    return theta_data.at(start+(3*(ic - 1) + ir - 1));
}

double
MBCSharedMemNodal::GetTheta(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_THETA);
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    unsigned start = (n-1) * nrotvalspernode;

    return theta_data.at(start+(idx-1));
}

double
MBCSharedMemNodal::GetEuler123(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MBC_ROT_EULER_123);
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    unsigned start = (n-1) * nrotvalspernode;

    return theta_data.at(start+(idx-1));
}

double
MBCSharedMemNodal::GetXP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    return XP_data.at(start+(idx-1));
}

double
MBCSharedMemNodal::GetOmega(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    return omega_data.at(start+(idx-1));;
}

double
MBCSharedMemNodal::GetXPP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    return XPP_data.at(start+(idx-1));;

}

double
MBCSharedMemNodal::GetOmegaP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	assert(bAccelerations());
	if (n < 1 || n > GetNodes()) throw;
	if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    return omegaP_data.at(start+(idx-1));;
}




const std::vector<uint32_t> *
MBCSharedMemNodal::GetDynamicsLabels(void) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	return &labels;
}

const uint32_t&
MBCSharedMemNodal::DynamicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return labels[n - 1];
}

//uint32_t&
//MBCSharedMemNodal::DynamicsLabel(uint32_t n)
//{
//	assert(GetStatus() == READY);
//	assert(bLabels());
//	if (n < 1 || n > GetNodes()) throw;
//	return labels[n - 1];
//}

const std::vector<double> *
MBCSharedMemNodal::GetF(void) const
{
	assert(GetStatus() == READY);
	return &force_data;
}

const std::vector<double> *
MBCSharedMemNodal::GetM(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	return &moment_data;
}

const std::vector<double>
MBCSharedMemNodal::GetF(uint32_t n) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> forces = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        forces.at(i) = force_data.at(start+i);
    }

    return forces;
}

const std::vector<double>
MBCSharedMemNodal::GetM(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != MBC_ROT_NONE);
	if (n < 1 || n > GetNodes()) throw;

    int start = (n-1) * 3;

    std::vector<double> moments = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        moments.at(i) = moment_data.at(start+i);
    }

    return moments;
}

void
MBCSharedMemNodal::SetF(uint32_t n, std::vector<double> fxyz)
{
    if (n < 1 || n > GetNodes()) throw;
    if (fxyz.size() != 3) throw;

    int start = (n-1) * 3;

    for (int i = 0; i < 3; i++) {
        force_data.at(start+i) = fxyz.at(i);
    }
}

void
MBCSharedMemNodal::SetF(uint32_t n, unsigned idx, double fidx)
{
    if (n < 1 || n > GetNodes()) throw;
    if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    force_data.at(start+(idx-1)) = fidx;

}

void
MBCSharedMemNodal::SetM(uint32_t n, std::vector<double> m123)
{
    if (n < 1 || n > GetNodes()) throw;
    if (m123.size() != 3) throw;

    int start = (n-1) * 3;

    for (int i = 0; i < 3; i++) {
        moment_data.at(start+i) = m123.at(i);
    }

}

void
MBCSharedMemNodal::SetM(uint32_t n, unsigned idx, double midx)
{
    if (n < 1 || n > GetNodes()) throw;
    if (idx < 1 || idx > 3) throw;

    int start = (n-1) * 3;

    moment_data.at(start+(idx-1)) = midx;

}

//const double&
//MBCSharedMemNodal::F(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
//}
//
//double&
//MBCSharedMemNodal::F(uint32_t n, uint8_t idx)
//{
//	assert(GetStatus() == READY);
//	if (idx <   ` 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
//}
//
//const double&
//MBCSharedMemNodal::M(uint32_t n, uint8_t idx) const
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() != NONE);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
//}
//
//double&
//MBCSharedMemNodal::M(uint32_t n, uint8_t idx)
//{
//	assert(GetStatus() == READY);
//	assert(GetRot() != NONE);
//	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
//	return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
//}

//mbc_t *
//MBCModal::GetBasePtr(void) const
//{
//	return (mbc_t*)&mbc;
//}
//
//mbc_refnode_stub_t *
//MBCModal::GetRefNodePtr(void) const
//{
//	return (mbc_refnode_stub_t *)&mbc;
//}
//
//MBCModal::MBCModal(void)
//{
//	memset(&mbc, 0, sizeof(mbc));
//}
//
//MBCModal::MBCModal(MBCType refnode_rot, unsigned modes)
//{
//	if (Initialize(refnode_rot, modes)) {
//		throw;
//	}
//}
//
//MBCModal::~MBCModal(void)
//{
//	Close();
//}
//
//MBCSharedMemBase::Type
//MBCModal::GetType(void) const
//{
//	return MODAL;
//}
//
//int
//MBCModal::Initialize(MBCType refnode_rot, unsigned modes)
//{
//	if (GetStatus() != NOT_READY) return -1;
//	memset(&mbc, 0, sizeof(mbc));
//	int rc = mbc_modal_init(&mbc, refnode_rot, modes);
//	if (rc == 0) const_cast<MBCModal *>(this)->SetStatus(INITIALIZED);
//	return rc;
//}
//
//int
//MBCModal::Negotiate(void) const
//{
//	if (GetStatus() != SHARED_MEMORY_READY) return -1;
//	int rc = mbc_modal_negotiate_request(&mbc);
//	if (rc == 0) const_cast<MBCModal *>(this)->SetStatus(READY);
//	return rc;
//}
//
//int
//MBCModal::PutForces(bool bConverged) const
//{
//	if (GetStatus() != READY) return -1;
//	return mbc_modal_put_forces(&mbc, bConverged);
//}
//
//int
//MBCModal::GetMotion(void) const
//{
//	if (GetStatus() != READY) return -1;
//	return mbc_modal_get_motion(&mbc);
//}
//
//int
//MBCModal::Close(void) const
//{
//	int rc = -1;
//	if (GetStatus() == READY) {
//		rc = mbc_modal_destroy(&mbc);
//		const_cast<MBCModal *>(this)->SetStatus(CLOSED);
//	}
//	return rc;
//}
//
//uint32_t
//MBCModal::KinematicsLabel(void) const
//{
//	return MBCSharedMemBase::KinematicsLabel();
//}
//
//const double&
//MBCModal::X(uint8_t idx) const
//{
//	return MBCSharedMemBase::X(idx);
//}
//
//const double&
//MBCModal::R(uint8_t ir, uint8_t ic) const
//{
//	return MBCSharedMemBase::R(ir, ic);
//}
//
//const double&
//MBCModal::Theta(uint8_t idx) const
//{
//	return MBCSharedMemBase::Theta(idx);
//}
//
//const double&
//MBCModal::Euler123(uint8_t idx) const
//{
//	return MBCSharedMemBase::Euler123(idx);
//}
//
//const double&
//MBCModal::XP(uint8_t idx) const
//{
//	return MBCSharedMemBase::XP(idx);
//}
//
//const double&
//MBCModal::Omega(uint8_t idx) const
//{
//	return MBCSharedMemBase::Omega(idx);
//}
//
//const double&
//MBCModal::XPP(uint8_t idx) const
//{
//	return MBCSharedMemBase::XPP(idx);
//}
//
//const double&
//MBCModal::OmegaP(uint8_t idx) const
//{
//	return MBCSharedMemBase::OmegaP(idx);
//}
//
//const uint32_t&
//MBCModal::DynamicsLabel(void) const
//{
//	return MBCSharedMemBase::DynamicsLabel();
//}
//
//uint32_t&
//MBCModal::DynamicsLabel(void)
//{
//	return MBCSharedMemBase::DynamicsLabel();
//}
//
//const double&
//MBCModal::F(uint8_t idx) const
//{
//	return MBCSharedMemBase::F(idx);
//}
//
//double&
//MBCModal::F(uint8_t idx)
//{
//	return MBCSharedMemBase::F(idx);
//}
//
//const double&
//MBCModal::M(uint8_t idx) const
//{
//	return MBCSharedMemBase::M(idx);
//}
//
//double&
//MBCModal::M(uint8_t idx)
//{
//	return MBCSharedMemBase::M(idx);
//}
//
//uint32_t
//MBCModal::GetModes(void) const
//{
//	assert(GetStatus() == READY);
//	return mbc.modes;
//}
//
//const double *const
//MBCModal::GetQ(void) const
//{
//	assert(GetStatus() == READY);
//	return MBC_Q(&mbc);
//}
//
//const double *const
//MBCModal::GetQP(void) const
//{
//	assert(GetStatus() == READY);
//	return MBC_QP(&mbc);
//}
//
//const double&
//MBCModal::Q(uint32_t m) const
//{
//	assert(GetStatus() == READY);
//	if (m < 1 || m > GetModes()) throw;
//	return (MBC_Q(&mbc))[m - 1];
//}
//
//const double&
//MBCModal::QP(uint32_t m) const
//{
//	assert(GetStatus() == READY);
//	if (m < 1 || m > GetModes()) throw;
//	return (MBC_QP(&mbc))[m - 1];
//}
//
//const double *
//MBCModal::GetP(void) const
//{
//	assert(GetStatus() == READY);
//	return MBC_P(&mbc);
//}
//
//const double&
//MBCModal::P(uint32_t m) const
//{
//	assert(GetStatus() == READY);
//	if (m < 1 || m > GetModes()) throw;
//	return (MBC_P(&mbc))[m - 1];
//}
//
//double&
//MBCModal::P(uint32_t m)
//{
//	assert(GetStatus() == READY);
//	if (m < 1 || m > GetModes()) throw;
//	return (MBC_P(&mbc))[m - 1];
//}

// #endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

