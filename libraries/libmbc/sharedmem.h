/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/libraries/libmbc/sock.h,v 1.9 2017/01/12 14:43:43 masarati Exp $ */
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

#ifndef SHAREDMEM_H
#define SHAREDMEM_H

//#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <iostream>
#include <string>
#include "mbc.h"

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/vector.hpp>

// see the following references for further information:
//
// http://www.boost.org/doc/libs/1_55_0/doc/html/interprocess/synchronization_mechanisms.html#interprocess.synchronization_mechanisms.conditions
//

namespace bi = boost::interprocess;

namespace mbdyn {

// Define an STL compatible allocator of doubles that allocates from the managed_shared_memory.
// This allocator will allow placing stl containers in the segment
typedef bi::allocator<double, bi::managed_shared_memory::segment_manager>  ShmemAllocator_double;

// Alias a vector that uses the previous STL-like allocator so that it
// allocates its values from the segment
typedef bi::vector<double, ShmemAllocator_double> shmVector_double;

// Define an STL compatible allocator for uint32_t in shared memory
typedef bi::allocator<uint32_t, bi::managed_shared_memory::segment_manager>  ShmemAllocator_uint32_t;

// Shared memory vector for uint32_t type for labels
typedef bi::vector<uint32_t, ShmemAllocator_uint32_t> shmVector_uint32_t;


class shared_memory_buffer
{

private:

    unsigned nrotvalspernode;

    // buffers for motion data
    shmVector_double *X_data;
    shmVector_double *XP_data;
    shmVector_double *XPP_data;
    shmVector_double *theta_data;
    shmVector_double *omega_data;
    shmVector_double *omegaP_data;

    // buffers for force data
    shmVector_double *force_data;
    shmVector_double *moment_data;

    // buffer for reference node motion data
    shmVector_double *refnode_kinematics_data;

    // buffer for reference node force data
    shmVector_double *refnode_dynamics_data;

    // buffer for label data
    shmVector_uint32_t *labels;

    bi::managed_shared_memory segment;


public:
    // semaphore for read/write control
//    bi::interprocess_semaphore server, client, mutex;

    // Mutex to protect access to the queue
    bi::interprocess_mutex      mutex;

    // Condition to wait until mbdyn posts command
    bi::interprocess_condition  cond_mbdyn_ready_for_cmd;

    bi::interprocess_condition  cond_mbdyn_cmd_available;

    // Condition to wait until peer posts command
    bi::interprocess_condition  cond_peer_ready_for_cmd;

    bi::interprocess_condition  cond_peer_cmd_available;

    // Condition to wait until mbdyn posts command
    bool  mbdyn_ready_for_cmd;
    bool  mbdyn_cmd_available;

    // Condition to wait until peer posts command
    bool  peer_ready_for_cmd;
    bool  peer_cmd_available;

    bi::interprocess_condition  cond_peer_ready_for_data;
    bi::interprocess_condition  cond_mbdyn_ready_for_data;
    bi::interprocess_condition  cond_mbdyn_data_available;
    bi::interprocess_condition  cond_peer_data_available;

    bool  peer_data_available;
    bool  mbdyn_data_available;
    bool  peer_ready_for_data;
    bool  mbdyn_ready_for_data;

    // reference node label
    uint32_t refnode_label;

    // variable for command/response values
    uint8_t cmd;

    // information about problem
    unsigned  NodalOrModal;
    bool UseReferanceNode;
    unsigned  RotationType;
    unsigned  RefRotationType;
    bool UseLabels;
    bool OutputAccelerations;
    unsigned NumNodes;

    // Server initialized with one to start, client must wait for it to be released
    shared_memory_buffer ( unsigned int nnodes,
                           unsigned int rot_type,
                           std::string shmName )
    {

        // get the shared memory segment
        segment = bi::managed_shared_memory (bi::open_only, shmName.c_str ());

        // Initialize shared memory STL-compatible allocators
        const ShmemAllocator_double alloc_inst_d (segment.get_segment_manager());
        const ShmemAllocator_uint32_t alloc_inst_u (segment.get_segment_manager());


        // Construct the data vectors
        X_data = segment.construct<shmVector_double>("X_data_vector")(alloc_inst_d);
        XP_data = segment.construct<shmVector_double>("XP_data_vector")(alloc_inst_d);
        XPP_data = segment.construct<shmVector_double>("XPP_data_vector")(alloc_inst_d);
        theta_data = segment.construct<shmVector_double>("theta_data_vector")(alloc_inst_d);
        omega_data = segment.construct<shmVector_double>("omega_data_vector")(alloc_inst_d);
        omegaP_data = segment.construct<shmVector_double>("omegaP_data_vector")(alloc_inst_d);
        force_data = segment.construct<shmVector_double>("force_data_vector")(alloc_inst_d);
        moment_data = segment.construct<shmVector_double>("moment_data_vector")(alloc_inst_d);
        // node labels
        labels = segment.construct<shmVector_uint32_t>("label_vector")(alloc_inst_u);

        // Construct the reference node data vectors
        refnode_kinematics_data = segment.construct<shmVector_double>("refnode_kinematics_data_vector")(alloc_inst_d);
        refnode_dynamics_data = segment.construct<shmVector_double>("refnode_dynamics_data_vector")(alloc_inst_d);

        refnode_label = 0;

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

        // position, velocity and acceleration in 6dof
        X_data->resize (3 * nnodes);
        XP_data->resize (3 * nnodes);
        XPP_data->resize (3 * nnodes);
        theta_data->resize (nrotvalspernode * nnodes);
        omega_data->resize (3 * nnodes);
        omegaP_data->resize (3 * nnodes);

//        std::cout << "In shared_memory_buffer, X_data contains: ";
//        for (unsigned i; i < X_data->size (); i++)
//            std::cout << X_data->at(i) << ' ';
//
//        std::cout << std::endl;

        // force in 6dof
        force_data->resize (3 * nnodes);
        moment_data->resize (3 * nnodes);

        // force in 6dof
        refnode_kinematics_data->resize (3 + 3 + 3 + nrotvalspernode + 3 + 3);

        // n node labels
        labels->resize (nnodes);

        // the following are filled in by the remote process and checked
        // against what mbdyn extracted from the input file
        NodalOrModal = 0;
        UseReferanceNode = false;
        RotationType = 0;
        UseLabels = false;
        OutputAccelerations = false;
        NumNodes = 0;

        cmd = 0;

        peer_ready_for_cmd = false;
        mbdyn_ready_for_cmd = false;
        mbdyn_cmd_available = false;
        peer_cmd_available = false;
        peer_data_available = false;
        mbdyn_data_available = false;
        peer_ready_for_data = false;
        mbdyn_ready_for_data = false;

    }

};

} // namespace mbdyn

//#endif /* HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP */

#endif /* SHAREDMEM_H */
