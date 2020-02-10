/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/extforce.cc,v 1.50 2017/01/12 14:46:09 masarati Exp $ */
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
 * MERCHANTABILITY or FITNES[MBDYN] FO[MBDYN] A PARTICULA[MBDYN] PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "extforce.h"
#include "except.h"
#include "solver.h"

#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <sys/stat.h>

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <boost/interprocess/sync/scoped_lock.hpp>

#include "sharedmem.h"
#include "extsharedmem.h"

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

ExtSharedMemHandler::ExtSharedMemHandler (mbsleep_t SleepTime,
                                          std::string shmName,
                                          bool create)
: ExtRemoteHandler(SleepTime, true, false),
shmName(shmName), bCreate(create), bIsInitialised(false)
{
    pedantic_cout("ExtSharedMemHandler in constructor" << std::endl);

}

ExtSharedMemHandler::~ExtSharedMemHandler(void)
{
    pedantic_cout("ExtSharedMemHandler: in destructor" << std::endl);

    SendCmdBySharedMem (ES_LAST);

     if (bCreate) {
        //Remove shared memory on destruction
        bool result = shared_memory_object::remove(shmName.c_str());
        if (result == false)
        {
            pedantic_cout("ExtSharedMemHandler: in destructor, shared_memory_object::remove was false" << std::endl);
        }
     }
}

int ExtSharedMemHandler::roundUp(int numToRound, int multiple)
{
    if (multiple == 0)
        return numToRound;

    int remainder = numToRound % multiple;
    if (remainder == 0)
        return numToRound;

    return numToRound + multiple - remainder;
}

int
ExtSharedMemHandler::Initialise( unsigned int node_kinematics_size,
                                 unsigned int dynamics_size,
                                 unsigned int labels_size,
                                 unsigned int nnodes,
                                 unsigned int rot_type )
{

    pedantic_cout("ExtSharedMemHandler::Initialise, bCreate is " << bCreate << std::endl);

    if (bCreate) {
        // we are the server, note that we use double here, not the doublereal
        // typedef as the other program will not use doublereal
        unsigned node_kinematics_nbytes = sizeof (double) * node_kinematics_size;
        unsigned dynamics_nbytes = sizeof (double) * dynamics_size;
        unsigned labels_nbytes = sizeof (uint32_t) * labels_size;

        // initialise shared memory region, make it always at least one
        // page size (nearly always 4k bytes)
        unsigned int totbytes = mapped_region::get_page_size() + node_kinematics_nbytes + dynamics_nbytes + labels_nbytes;
        // and round up to the nearest page size to the number of bytes needed
        totbytes = roundUp (totbytes, mapped_region::get_page_size());

        // Erase previous shared memory if it exists
        shared_memory_object::remove(shmName.c_str ());

        pedantic_cout("ExtSharedMemHandler: in Initialise,  creating shared memory named '" << shmName << "' of size " << totbytes << " bytes" << std::endl);
        shm = managed_shared_memory(create_only, shmName.c_str (), totbytes);

        // there should only be one buffer in this shared memory region
        buffer = shm.construct<shared_memory_buffer>(unique_instance)(nnodes, rot_type, shmName);

    } else {

        // Open the managed segment
        shm = managed_shared_memory(open_only, shmName.c_str ());

        // there should only be one buffer in this shared memory region
        std::pair<shared_memory_buffer*, int> temp = shm.find<shared_memory_buffer>(unique_instance);

        buffer = temp.first;

        // TODO: check vector sizes?

    }

    // buffers for motion data
    X_data = shm.find<shmVector_double>("X_data_vector").first;
    XP_data = shm.find<shmVector_double>("XP_data_vector").first;
    XPP_data = shm.find<shmVector_double>("XPP_data_vector").first;
    theta_data = shm.find<shmVector_double>("theta_data_vector").first;
    omega_data = shm.find<shmVector_double>("omega_data_vector").first;
    omegaP_data = shm.find<shmVector_double>("omegaP_data_vector").first;

    // buffers for force data
    force_data = shm.find<shmVector_double>("force_data_vector").first;
    moment_data = shm.find<shmVector_double>("moment_data_vector").first;

    // buffer for reference node motion data
    refnode_kinematics_data = shm.find<shmVector_double>("refnode_kinematics_data_vector").first;

    // buffer for reference node force data
    refnode_dynamics_data = shm.find<shmVector_double>("refnode_dynamics_data_vector").first;

    // buffer for label data
    labels = shm.find<shmVector_uint32_t>("label_vector").first;

    nrotvalspernode = 0;

    switch (rot_type) {
    case MBC_ROT_NONE:
        break;

    case MBC_ROT_MAT:
        // full orientation matrix for orientations
        nrotvalspernode = 9;
        break;

    case MBC_ROT_THETA:
    case MBC_ROT_EULER_123:
        // 3 orientation values
        nrotvalspernode = 3;
        break;
    }

    bIsInitialised = true;

}

ExtFileHandlerBase::Negotiate
ExtSharedMemHandler::NegotiateRequest(void) const
{
    pedantic_cout("ExtSharedMemHandler: in NegotiateRequest" << std::endl);

	if (bCreate) {
		return ExtFileHandlerBase::NEGOTIATE_SERVER;

	} else {
		return ExtFileHandlerBase::NEGOTIATE_CLIENT;
	}
}

int
ExtSharedMemHandler::SendCmdBySharedMem(uint8_t cmd)
{
    pedantic_cout("ExtSharedMemHandler: in SendCmdBySharedMem, cmd to be sent is " << unsigned(cmd) << " : " << std::string (mbc_cmd2str(cmd)) << std::endl);

//    // lock the mutex so we can safely access the data
//    scoped_lock<interprocess_mutex> lock(buffer->mutex);
//
//    buffer->cond_peer_ready_for_cmd.wait (lock);
//
//    buffer->cmd = u;
//
//    // Notify to the other process that there is a message from mbdyn available
//    //buffer->cond_mbdyn_cmd_available-.notify_one_one ();
//
//    // Mark mbdyn_cmd_available
////    buffer->mbdyn_cmd_available = true;
//
//    buffer->cond_mbdyn_cmd_available.notify_one ();

    scoped_lock<interprocess_mutex> lock(buffer->mutex);

    pedantic_cout("[MBDYN] ExtSharedMemHandler::SendCmdBySharedMem: checking if peer is ready to receive" << std::endl);
    if (!buffer->peer_ready_for_cmd)
    {
        pedantic_cout("[MBDYN] waiting for peer be ready to receive" << std::endl);
        buffer->cond_peer_ready_for_cmd.wait (lock);
    }

    // send command to client
    buffer->cmd = cmd;

    pedantic_cout("[MBDYN] ExtSharedMemHandler::SendCmdBySharedMem: put cmd " << unsigned(buffer->cmd)  << " : " << std::string (mbc_cmd2str(buffer->cmd)) << " in buffer" << std::endl);

    pedantic_cout("[MBDYN] ExtSharedMemHandler::SendCmdBySharedMem: notifying new command is available" << std::endl);

    buffer->cond_mbdyn_cmd_available.notify_all();
    buffer->mbdyn_cmd_available = true;
    //buffer->peer_ready_for_cmd = false;
    buffer->mbdyn_ready_for_cmd = false;
    buffer->peer_cmd_available = false;

//
//    switch (NegotiateRequest()) {
//
//        case ExtFileHandlerBase::NEGOTIATE_SERVER:
//
//            // the server semaphore is initialised as 0 when constructed.
//            // If this is the server, we must wait for the client to
//            // send the negotiation request (and announce this by
//            // post'ing the server semaphore).
//            pedantic_cout("ExtSharedMemHandler: in SendCmdBySharedMem, NEGOTIATE_SERVER, buffer->server.wait ()" << std::endl);
//            buffer->server.wait ();
//
//            pedantic_cout("ExtSharedMemHandler: in SendCmdBySharedMem, NEGOTIATE_SERVER, writing cmd to shared memory" << std::endl);
//            buffer->cmd = u;
//
//            // allow the client to access the command variable
//            buffer->client.post ();
//
//            break;
//
//        case ExtFileHandlerBase::NEGOTIATE_CLIENT:
//
//            // the client semaphore is initialised as 1 when constructed.
//            // If this is the client, we will wait until the client semaphore
//            // is incremented (with 'post') before writing a command to the
//            // shared memory
//            pedantic_cout("ExtSharedMemHandler: in SendCmdBySharedMem, NEGOTIATE_CLIENT, buffer->client.wait ()" << std::endl);
//            buffer->client.wait ();
//
//            pedantic_cout("ExtSharedMemHandler: in SendCmdBySharedMem, NEGOTIATE_CLIENT, writing cmd to shared memory" << std::endl);
//
//            buffer->cmd = u;
//
//            // allow the server to access the command variable
//            buffer->server.post ();
//
//            break;
//    }



    pedantic_cout("ExtSharedMemHandler::SendCmdBySharedMem, sent cmd" << std::endl);

    return 0;

    // mutex lock goes out of scope here and is released
}

int
ExtSharedMemHandler::GetCmdBySharedMem(uint8_t &cmd)
{

//    // lock the mutex so we can safely access the data
//    scoped_lock<interprocess_mutex> lock(buffer->mutex);
//
//    // announce mbdyn is ready to receive command
//    buffer->cond_mbdyn_ready_for_cmd.notify_one ();
//
//    buffer->cond_peer_cmd_available.wait (lock);
//
//    cmd = buffer->cmd;


    scoped_lock<interprocess_mutex> lock(buffer->mutex);

    buffer->mbdyn_ready_for_cmd = true;

    pedantic_cout ("[MBDYN] ExtSharedMemHandler::GetCmdBySharedMem: notifying that server is ready to receive" << std::endl);
    buffer->cond_mbdyn_ready_for_cmd.notify_all();

    pedantic_cout ("[MBDYN] ExtSharedMemHandler::GetCmdBySharedMem: checking if peer cmd is available" << std::endl);
    if (!buffer->peer_cmd_available)
    {
        pedantic_cout ("[MBDYN] ExtSharedMemHandler::GetCmdBySharedMem waiting for command from peer" << std::endl);
        buffer->cond_peer_cmd_available.wait (lock);
    }

    pedantic_cout ("[MBDYN] ExtSharedMemHandler::GetCmdBySharedMem getting command" << std::endl);
    // first value should be 0
    cmd = buffer->cmd;

    pedantic_cout ("[MBDYN] ExtSharedMemHandler::GetCmdBySharedMem received cmd " << unsigned(cmd) << ": " << std::string(mbc_cmd2str(cmd)) << " from peer" << std::endl);

    buffer->mbdyn_ready_for_cmd = false;
    // no peer command is available (we've read it)
    //buffer->peer_cmd_available = false;
    //buffer->mbdyn_cmd_available = false;

//    peer_cmd_available = false;


//    pedantic_cout("ExtSharedMemHandler: in GetCmdBySharedMem" << std::endl);
//
//    switch (NegotiateRequest()) {
//
//        case ExtFileHandlerBase::NEGOTIATE_SERVER:
//
//            // wait for any other process accessing the shared memory area
//            // to release the semaphore
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_SERVER, about to buffer->server.wait ()" << std::endl);
//            buffer->server.wait ();
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_SERVER, about to buffer->cmd" << std::endl);
//            cmd = buffer->cmd;
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_SERVER, cmd was " << unsigned (cmd) << std::endl);
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_SERVER, about to buffer->client.post ();" << std::endl);
//            buffer->client.post ();
//
//            break;
//
//        case ExtFileHandlerBase::NEGOTIATE_CLIENT:
//
//            // wait for any other process accessing the shared memory area
//            // to release the semaphore
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_CLIENT, about to buffer->client.wait ();" << std::endl);
//            buffer->client.wait ();
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_CLIENT, about to buffer->cmd;" << std::endl);
//            cmd = buffer->cmd;
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_CLIENT, cmd was " << unsigned (cmd) << std::endl);
//            pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem, NEGOTIATE_CLIENT, about to buffer->server.post ()" << std::endl);
//            buffer->server.post ();
//
//            break;
//
//    }
//
//    pedantic_cout("ExtSharedMemHandler::GetCmdBySharedMem received cmd is " << unsigned(cmd) << std::endl);
//
    return 0;

    // mutex released here
}

//int
//ExtSharedMemHandler::GetCmdBySharedMem(uint8_t &cmd, mbsleep_t sleeptime)
//{
//
//    int out = 0;
//    bool success = false;
//
//    doublereal t_dur;
//    mbsleep_sleep2real (&sleeptime, &t_dur);
//
//    boost::posix_time::ptime timeout (boost::posix_time::microsec_clock::local_time() + boost::posix_time::milliseconds(t_dur));
//
//
//    switch (NegotiateRequest()) {
//
//        case ExtFileHandlerBase::NEGOTIATE_SERVER:
//
//            success = buffer->server.timed_wait (timeout);
//
//            if (success)
//            {
//                cmd = buffer->cmd;
//
//                buffer->client.post ();
//
//                out = 0;
//
//            } else {
//                out = -1;
//            }
//
//            break;
//
//        case ExtFileHandlerBase::NEGOTIATE_CLIENT:
//
//            success = buffer->client.timed_wait (timeout);
//
//            if (success)
//            {
//                cmd = buffer->cmd;
//
//                buffer->server.post ();
//
//                out = 0;
//
//            } else {
//                out = -1;
//            }
//
//            break;
//    }
//
//    return out;
//
//}

bool
ExtSharedMemHandler::Prepare_pre(void)
{
	uint8_t u;
	ssize_t rc;

	pedantic_cout("ExtSharedMemHandler::Prepare_pre -- in Prepare_pre" << std::endl);

	switch (NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_CLIENT:
	    /* This the server, send negotiation request */
		u = ES_NEGOTIATION;

		pedantic_cout("ExtSharedMemHandler::Prepare_pre -- NEGOTIATE_CLIENT -- sending negotiation request, cmd " << unsigned(u) << std::endl);

        SendCmdBySharedMem(u);

		pedantic_cout("ExtSharedMemHandler::Prepare_pre -- NEGOTIATE_CLIENT -- negotiation request sent" << std::endl);

		break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER:
        /* This the client, check if negotiation request has been sent */
	    pedantic_cout("ExtSharedMemHandler::Prepare_pre -- NEGOTIATE_SERVER -- trying to receive negotiation request" << std::endl);

		GetCmdBySharedMem(u);

		if (u != ES_NEGOTIATION) {
			silent_cerr("ExtSharedMemHandler::Prepare_pre -- NEGOTIATE_SERVER -- negotiation request failed"
				<< std::endl);
			return (bOK = false);
		}

		break;

	default:
		ASSERT(0);
	}

	pedantic_cout("ExtSharedMemHandler::Prepare_pre -- negotiation started" << std::endl);

	return true;
}

void
ExtSharedMemHandler::Prepare_post(bool ok)
{
    int result = 0;

    pedantic_cout("ExtSharedMemHandler::Prepare_post -- in Prepare_post" << std::endl);

	if (NegotiateRequest()) {

		uint8_t u = ok ? ES_OK : ES_ABORT;

        pedantic_cout("ExtSharedMemHandler::Prepare_post -- sending negotiation response" << std::endl);

		result = SendCmdBySharedMem(u);

		pedantic_cout("ExtSharedMemHandler::Prepare_post -- negotiation response sent" << std::endl);

	} else {
		uint8_t u;
		pedantic_cout("ExtSharedMemHandler::Prepare_post -- receiving negotiation response" << std::endl);

		result = GetCmdBySharedMem(u);

        pedantic_cout("ExtSharedMemHandler::Prepare_post -- negotiation response received" << std::endl);

		if (u != ES_OK) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

void
ExtSharedMemHandler::AfterPredict(void)
{
    pedantic_cout("ExtSharedMemHandler: in AfterPredict" << std::endl);
	bLastReadForce = false;
	bReadForces = true;
}

bool
ExtSharedMemHandler::Send_pre(SendWhen when)
{
    int result = 0;

    pedantic_cout("ExtSharedMemHandler: in Send_pre" << std::endl);
	if (!bReadForces || !bOK) {
		return false;
	}

	uint8_t u;
	if (when == SEND_AFTER_CONVERGENCE) {
		u = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;
	} else {
		u = ES_REGULAR_DATA;
	}
	pedantic_cout("ExtSharedMemHandler: sending when to send data" << std::endl);

	result = SendCmdBySharedMem(u);

//	if (rc == SOCKET_ERROR) {
//		int save_errno = WSAGetLastError();
//		silent_cerr("ExtSharedMemHandler: send() failed "
//			"(" << save_errno << ": " << strerror(save_errno) << ")"
//			<< std::endl);
//		mbdyn_set_stop_at_end_of_iteration();
//		return (bOK = false);
//
//	} else if (rc != sizeof(u)) {
//		silent_cerr("ExtSharedMemHandler: send() failed "
//			"(sent " << rc << " bytes "
//			"instead of " << sizeof(u) << ")"
//			<< std::endl);
//		mbdyn_set_stop_at_end_of_iteration();
//		return (bOK = false);
//	}

	pedantic_cout("ExtSharedMemHandler: sent when to send data" << std::endl);

	return true;
}


bool
ExtSharedMemHandler::Recv_pre(void)
{
    int result = 0;

    pedantic_cout("ExtSharedMemHandler: in Recv_pre" << std::endl);
	if (!bReadForces) {
		return false;
	}

	uint8_t u = 0;

    // wait until the status code is returned
	result = GetCmdBySharedMem(u);

//		if (rc == SOCKET_ERROR) {
//			int save_errno = WSAGetLastError();
//
//			if (WSAGetLastError() != EAGAIN) {
//				silent_cerr("ExtSharedMemHandler: "
//					"recv() failed (" << save_errno << ": "
//					<< sock_err_string(save_errno) << ")"
//					<< std::endl);
//				mbdyn_set_stop_at_end_of_iteration();
//				return (bOK = false);
//			}
//
//		} else if (rc != sizeof(u)) {
//			silent_cerr("ExtSharedMemHandler: "
//				"recv()=" << rc << " (expected " << sizeof(u) << ")"
//				<< std::endl);
//			mbdyn_set_stop_at_end_of_iteration();
//			return (bOK = false);
//		}

	return ActOnCmd(u);
}

mbdyn::shared_memory_buffer*
ExtSharedMemHandler::GetSharedMemBuffer(void)
{
    CheckInit ();

    return buffer;
}

// reference node getters
std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeX (void)
{
    CheckInit ();

    int start = 0;

    std::vector<double> X = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        X.at(i) = refnode_kinematics_data->at(start+i);
    }

    return X;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeXP (void)
{
    CheckInit ();

    int start = 3;

    std::vector<doublereal> XP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XP.at(i) = refnode_kinematics_data->at(start+i);
    }

    return XP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeXPP (void)
{
    CheckInit ();

    int start = 6;

    std::vector<doublereal> XPP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XPP.at(i) = refnode_kinematics_data->at(start+i);
    }

    return XPP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeTheta (void)
{
    CheckInit ();

    unsigned start = 9;

    std::vector<doublereal> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = refnode_kinematics_data->at(start+i);
    }

    return theta;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeOmega (void)
{
    CheckInit ();

    int start = 9 + nrotvalspernode;

    std::vector<doublereal> omega = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omega.at(i) = refnode_kinematics_data->at(start+i);
    }

    return omega;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeOmegaP (void)
{
    CheckInit ();

    int start = 12 + nrotvalspernode;

    std::vector<doublereal> omegaP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omegaP.at(i) = refnode_kinematics_data->at(start+i);
    }

    return omegaP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeF (void)
{
    CheckInit ();

    int start = 0;

    std::vector<doublereal> forces = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        forces.at(i) = refnode_dynamics_data->at(start+i);
    }

    return forces;
}

std::vector<doublereal>
ExtSharedMemHandler::GetRefNodeM (void)
{
    CheckInit ();

    int start = 3;

    std::vector<doublereal> moments = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        moments.at(i) = refnode_dynamics_data->at(start+i);
    }

    return moments;
}

uint32_t
ExtSharedMemHandler::GetRefLabel (void) {

    CheckInit ();

    return buffer->refnode_label;

}

// reference node setters
int
ExtSharedMemHandler::SetRefLabel (uint32_t l) {

    CheckInit ();

    buffer->refnode_label = l;

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeX (const std::vector<doublereal> &x)
{
    CheckInit ();

    int start = 0;

    for (int i = 0; i < 3; i++) {
        refnode_kinematics_data->at(start+i) = double(x.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeXP (const std::vector<doublereal> &xp)
{
    CheckInit ();

    int start = 3;

    for (int i = 0; i < 3; i++) {
        refnode_kinematics_data->at(start+i) = double(xp.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeXPP (const std::vector<doublereal> &xpp)
{
    CheckInit ();

    int start = 6;

    for (int i = 0; i < 3; i++) {
        refnode_kinematics_data->at(start+i) = double(xpp.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeTheta (const std::vector<doublereal> &theta)
{
    CheckInit ();

    int start = 9;

    for (int i = 0; i < nrotvalspernode; i++) {
        refnode_kinematics_data->at(start+i) = double(theta.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeOmega (const std::vector<doublereal> &omega)
{
    CheckInit ();

    int start = 9 + nrotvalspernode;

    for (int i = 0; i < 3; i++) {
        refnode_kinematics_data->at(start+i) = double(omega.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeOmegaP (const std::vector<doublereal> &omegaP)
{

    CheckInit ();

    int start = 12 + nrotvalspernode;

    for (int i = 0; i < 3; i++) {
        refnode_kinematics_data->at(start+i) = double(omegaP.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeF (const std::vector<doublereal> &forces)
{

    CheckInit ();

    int start = 0;

    for (int i = 0; i < 3; i++) {
        refnode_dynamics_data->at(start+i) = double(forces.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetRefNodeM (const std::vector<doublereal> &moments)
{

    CheckInit ();

    int start = 3;

    for (int i = 0; i < 3; i++) {
        refnode_dynamics_data->at(start+i) = double(moments.at(i));
    }

    return 0;
}

// normal node getters
std::vector<doublereal>
ExtSharedMemHandler::GetNodeX (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> X = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        X.at(i) = X_data->at(start+i);
    }

    return X;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeXP (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> XP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XP.at(i) = XP_data->at(start+i);
    }

    return XP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeXPP (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> XPP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        XPP.at(i) = XPP_data->at(start+i);
    }

    return XPP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeTheta (unsigned nodenum)
{

    CheckInit ();

    unsigned start = nodenum * nrotvalspernode;

    std::vector<doublereal> theta;
    theta.resize (nrotvalspernode);

    for (int i = 0; i < nrotvalspernode; i++) {
        theta.at(i) = theta_data->at(start+i);
    }

    return theta;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeOmega (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> omega = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omega.at(i) = omega_data->at(start+i);
    }

    return omega;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeOmegaP (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> omegaP = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        omegaP.at(i) = omegaP_data->at(start+i);
    }

    return omegaP;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeF (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> forces = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        forces.at(i) = force_data->at(start+i);
    }

    return forces;
}

std::vector<doublereal>
ExtSharedMemHandler::GetNodeM (unsigned nodenum)
{

    CheckInit ();

    int start = nodenum * 3;

    std::vector<doublereal> moments = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        moments.at(i) = moment_data->at(start+i);
    }

    return moments;
}

uint32_t
ExtSharedMemHandler::GetLabel (unsigned nodenum)
{

    CheckInit ();

    return labels->at(nodenum);

}

// setters
int
ExtSharedMemHandler::SetNodeX (unsigned nodenum, const std::vector<doublereal> &X) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        X_data->at(start+i) = double(X.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeXP (unsigned nodenum, const std::vector<doublereal> &XP) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        XP_data->at(start+i) = double(XP.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeXPP (unsigned nodenum, const std::vector<doublereal> &XPP) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        XPP_data->at(start+i) = double(XPP.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeTheta (unsigned nodenum, const std::vector<doublereal> &theta) {

    CheckInit ();

    unsigned start = nodenum * nrotvalspernode;

    for (int i = 0; i < nrotvalspernode; i++) {
        theta_data->at(start+i) = double(theta.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeOmega (unsigned nodenum, const std::vector<doublereal> &omega) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        omega_data->at(start+i) = double(omega.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeOmegaP (unsigned nodenum, const std::vector<doublereal> &omegaP) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        omegaP_data->at(start+i) = double(omegaP.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeF (unsigned nodenum, const std::vector<doublereal> &forces) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        force_data->at(start+i) = double(forces.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetNodeM (unsigned nodenum, const std::vector<doublereal> &moments) {

    CheckInit ();

    int start = nodenum * 3;

    for (int i = 0; i < 3; i++) {
        moment_data->at(start+i) = double(moments.at(i));
    }

    return 0;
}

int
ExtSharedMemHandler::SetLabel (unsigned nodenum, uint32_t l) {

    CheckInit ();

    labels->at(nodenum) = l;

    return 0;
}

void
ExtSharedMemHandler::CheckInit (void)
{
    if (!bIsInitialised)
    {
        silent_cerr("ExtSharedMemHandler"
				"A attempt to access shared memory was made before initialisation"
				<< std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
}

#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

ExtFileHandlerBase *
ReadExtSharedMemHandler(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
    pedantic_cout("In ReadExtSharedMemHandler" << std::endl);
#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
	ExtFileHandlerBase *pEFH = 0;

	bool bCreate = true;
	std::string shm_name;

	if (HP.IsKeyWord("create")) {
		if (!HP.GetYesNo(bCreate)) {
			silent_cerr("ExtSharedMemHandler"
				"(" << uLabel << "): "
				"\"create\" must be either \"yes\" or \"no\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("name")) {

		const char *h;

		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("ExtSharedMemHandler"
				"(" << uLabel << "): "
				"unable to read name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		shm_name = h;

	} else {
		silent_cerr("ExtSharedMemHandler"
			"(" << uLabel << "): "
			"shared memory region name undefined "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pedantic_cout("Finished reading shared memory communicator.\n"
               << "Name: " << shm_name << std::endl
               << "Create: " << bCreate << std::endl);

	mbsleep_t SleepTime = mbsleep_init(0);
	std::streamsize Precision = 0;
	ReadExtFileParams(pDM, HP, uLabel, SleepTime, Precision);
	// NOTE: so far, precision is ignored

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtSharedMemHandler,
		ExtSharedMemHandler(SleepTime, shm_name, bCreate));
    pedantic_cout("In ReadExtSocketHandler at end" << std::endl);
	return pEFH;
#else /* HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP */
	silent_cerr("ExtSharedMemHandler not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP */
}



