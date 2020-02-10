/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/extforce.h,v 1.30 2017/01/12 14:46:09 masarati Exp $ */
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

/* Forza */

#ifndef EXTSHAREDMEM_H
#define EXTSHAREDMEM_H

#include <mbconfig.h>

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <vector>
#include <string>

#include "mbsleep.h"
#include "force.h"
#include "converged.h"

#include "mbc.h"

#include "extforce.h"
#include "sharedmem.h"

class ExtSharedMemHandler : public ExtRemoteHandler {
protected:

    boost::interprocess::managed_shared_memory shm;
	mbdyn::shared_memory_buffer * buffer;
	std::string shmName;
    bool bCreate;
    bool bIsInitialised;
    unsigned nrotvalspernode;

    // buffers for motion data
    mbdyn::shmVector_double *X_data;
    mbdyn::shmVector_double *XP_data;
    mbdyn::shmVector_double *XPP_data;
    mbdyn::shmVector_double *theta_data;
    mbdyn::shmVector_double *omega_data;
    mbdyn::shmVector_double *omegaP_data;

    // buffers for force data
    mbdyn::shmVector_double *force_data;
    mbdyn::shmVector_double *moment_data;

    // buffer for reference node motion data
    mbdyn::shmVector_double *refnode_kinematics_data;

    // buffer for reference node force data
    mbdyn::shmVector_double *refnode_dynamics_data;

    // buffer for label data
    mbdyn::shmVector_uint32_t *labels;

	int SendCmdBySharedMem(uint8_t u);
	int GetCmdBySharedMem(uint8_t &cmd);
//	int GetCmdBySharedMem(uint8_t &cmd, mbsleep_t SleepTime);
    int roundUp(int numToRound, int multiple);

    void CheckInit(void);

public:
	ExtSharedMemHandler(mbsleep_t SleepTime,
                        std::string shmName,
                        bool create);
	virtual ~ExtSharedMemHandler(void);

	int Initialise(unsigned int node_kinematics_size,
                   unsigned int dynamics_size,
                   unsigned int labels_size,
                   unsigned int nnodes,
                   unsigned int rot_type );

	virtual bool Prepare_pre(void);
	virtual Negotiate NegotiateRequest(void) const;
	virtual void Prepare_post(bool ok);

	virtual void AfterPredict(void);

	virtual bool Send_pre(SendWhen when);
	virtual void Send_post(SendWhen when) { };

	virtual bool Recv_pre(void);
	//virtual bool Recv_post(void) { return true; };

	virtual ExtFileHandlerBase::Type GetType(void) { return ExtFileHandlerBase::TYPE_SHARED_MEMORY; };

	// returns pointer to mbdyn::shared_memory_buffer *
	virtual mbdyn::shared_memory_buffer* GetSharedMemBuffer(void);

	//reference node getters
	std::vector<doublereal> GetRefNodeX (void);
    std::vector<doublereal> GetRefNodeXP (void);
    std::vector<doublereal> GetRefNodeXPP (void);
    std::vector<doublereal> GetRefNodeTheta (void);
    std::vector<doublereal> GetRefNodeOmega (void);
    std::vector<doublereal> GetRefNodeOmegaP (void);
    std::vector<doublereal> GetRefNodeF (void);
    std::vector<doublereal> GetRefNodeM (void);
    uint32_t GetRefLabel (void);

    // reference node setters
    int SetRefNodeX (const std::vector<doublereal> &x);
    int SetRefNodeXP (const std::vector<doublereal> &xp);
    int SetRefNodeXPP (const std::vector<doublereal> &xpp);
    int SetRefNodeTheta (const std::vector<doublereal> &theta);
    int SetRefNodeOmega (const std::vector<doublereal> &omega);
    int SetRefNodeOmegaP (const std::vector<doublereal> &omegaP);
    int SetRefNodeF (const std::vector<doublereal> &forces);
    int SetRefNodeM (const std::vector<doublereal> &moments);
    int SetRefLabel (uint32_t l);

    // normal node getters
    std::vector<doublereal> GetNodeX (unsigned nodenum);
    std::vector<doublereal> GetNodeXP (unsigned nodenum);
    std::vector<doublereal> GetNodeXPP (unsigned nodenum);
    std::vector<doublereal> GetNodeTheta (unsigned nodenum);
    std::vector<doublereal> GetNodeOmega (unsigned nodenum);
    std::vector<doublereal> GetNodeOmegaP (unsigned nodenum);
    std::vector<doublereal> GetNodeF (unsigned nodenum);
    std::vector<doublereal> GetNodeM (unsigned nodenum);
    uint32_t GetLabel (unsigned nodenum);

    // normal node setters
    int SetNodeX (unsigned nodenum, const std::vector<doublereal> &x);
    int SetNodeXP (unsigned nodenum, const std::vector<doublereal> &xp);
    int SetNodeXPP (unsigned nodenum, const std::vector<doublereal> &xpp);
    int SetNodeTheta (unsigned nodenum, const std::vector<doublereal> &theta);
    int SetNodeOmega (unsigned nodenum, const std::vector<doublereal> &omega);
    int SetNodeOmegaP (unsigned nodenum, const std::vector<doublereal> &omegaP);
    int SetNodeF (unsigned nodenum, const std::vector<doublereal> &forces);
    int SetNodeM (unsigned nodenum, const std::vector<doublereal> &moments);
    int SetLabel (unsigned nodenum, uint32_t l);

};
#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

class DataManager;
class MBDynParser;

extern ExtFileHandlerBase *
ReadExtSharedMemHandler(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel);

#endif // EXTSHAREDMEM_H

