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

#ifndef MBCXXSHARED_H
#define MBCXXSHARED_H

// define  HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP for sharedmem.h
//#ifndef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
//#define HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP 1
//#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <boost/interprocess/managed_shared_memory.hpp>

#include <vector>
#include <sharedmem.h>
#include <mbc.h>

class MBCSharedMemBase {


public:
	enum Type {
		MODAL = MBC_MODAL,
		NODAL = MBC_NODAL
	};

	enum Status {
		NOT_READY,
		INITIALIZED,
		SHARED_MEMORY_READY,
		READY,
		FINISHED
	} m_status;

	virtual MBCSharedMemBase::Type GetType(void) const = 0;

	MBCType GetRot(void) const;
	bool bRefNode(void) const;
	MBCType GetRefNodeRot(void) const;
	bool bAccelerations(void) const;
	bool bLabels(void) const;

	void SetTimeout(int t);
	void SetVerbose(bool bv);
	void SetDataAndNext(bool bd);

	bool bVerbose(void) const;
	bool bDataAndNext(void) const;

protected:


	int _timeout;
	unsigned nnodes;

	bool _bVerbose;
	bool _bDataAndNext;
    bool _bOutputAccelerations;
	bool _bUseLabels;
	bool _bUseReferanceNode;
	bool _bRefNode;

	uint8_t cmd;

	MBCType _Rot;
	MBCType _RefNodeRot;

	boost::interprocess::managed_shared_memory shm;
	mbdyn::shared_memory_buffer* buffer;

    // shared memory buffers for motion data
    mbdyn::shmVector_double *S_X_data;
    mbdyn::shmVector_double *S_XP_data;
    mbdyn::shmVector_double *S_XPP_data;
    mbdyn::shmVector_double *S_theta_data;
    mbdyn::shmVector_double *S_omega_data;
    mbdyn::shmVector_double *S_omegaP_data;

    // shared memory buffers for force data
    mbdyn::shmVector_double *S_force_data;
    mbdyn::shmVector_double *S_moment_data;

    // shared memory buffer for reference node motion data
    mbdyn::shmVector_double *S_refnode_kinematics_data;

    // shared memory buffer for reference node force and moment data
    mbdyn::shmVector_double *S_refnode_dynamics_data;

    // shared memory buffer for label data
    mbdyn::shmVector_uint32_t *S_labels;

    // buffers for motion data
    std::vector<double> X_data;
    std::vector<double> XP_data;
    std::vector<double> XPP_data;
    std::vector<double> theta_data;
    std::vector<double> omega_data;
    std::vector<double> omegaP_data;

    // buffers for force data
    std::vector<double> force_data;
    std::vector<double> moment_data;

    // buffer for reference node motion data
    std::vector<double> refnode_kinematics_data;

    // buffer for reference node force data
    std::vector<double> refnode_dynamics_data;

    // buffer for label data
    std::vector<uint32_t> labels;

    uint32_t refnode_label;

    unsigned nrotvalspernode;
    unsigned nrefnoderotvalspernode;

	void SetStatus(Status s);
    int CheckCmd(const uint8_t cmd) const;

public:
	MBCSharedMemBase(void);
	virtual ~MBCSharedMemBase(void);

	int Init();

    Status GetStatus(void) const;
	virtual int Negotiate(void) = 0;
	virtual int PutForces(bool bConverged) = 0;
	virtual int GetMotion(void) = 0;
	virtual int Close(void) const = 0;

	int PutCmd(const uint8_t cmd) const;
	int GetCmd(uint8_t &cmd) const;

	uint32_t GetRefNodeKinematicsLabel(void) const;

//	uint32_t KinematicsLabel(void) const;

	std::vector<double> GetRefNodeX(void) const;
	std::vector<double> GetRefNodeR(void) const;
	std::vector<double> GetRefNodeTheta(void) const;
	std::vector<double> GetRefNodeEuler123(void) const;
	std::vector<double> GetRefNodeXP(void) const;
	std::vector<double> GetRefNodeOmega(void) const;
	std::vector<double> GetRefNodeXPP(void) const;
	std::vector<double> GetRefNodeOmegaP(void) const;

	double GetRefNodeX(uint8_t idx) const;
	double GetRefNodeR(uint8_t ir, uint8_t ic) const;
	double GetRefNodeTheta(uint8_t idx) const;
	double GetRefNodeEuler123(uint8_t idx) const;
	double GetRefNodeXP(uint8_t idx) const;
	double GetRefNodeOmega(uint8_t idx) const;
	double GetRefNodeXPP(uint8_t idx) const;
	double GetRefNodeOmegaP(uint8_t idx) const;

//	const double& X(uint8_t idx) const;
//	const double& R(uint8_t ir, uint8_t ic) const;
//	const double& Theta(uint8_t idx) const;
//	const double& Euler123(uint8_t idx) const;
//	const double& XP(uint8_t idx) const;
//	const double& Omega(uint8_t idx) const;
//	const double& XPP(uint8_t idx) const;
//	const double& OmegaP(uint8_t idx) const;

	uint32_t GetRefNodeDynamicsLabel(void) const;
//	const uint32_t& DynamicsLabel(void) const;
//	uint32_t& DynamicsLabel(void);

	std::vector<double> GetRefNodeF(void) const;
	std::vector<double> GetRefNodeM(void) const;

	void SetRefNodeF(uint8_t idx, const double &force);
	void SetRefNodeM(uint8_t idx, const double &moment);

//	const double& F(uint8_t idx) const;
//	double& F(uint8_t idx);
//	const double& M(uint8_t idx) const;
//	double& M(uint8_t idx);

    void SetF(uint8_t idx);
    void SetM(uint8_t idx);
};

class MBCSharedMemNodal : public MBCSharedMemBase {
private:

//	virtual mbc_t *GetBasePtr(void) const;
//	virtual mbc_refnode_stub_t *GetRefNodePtr(void) const;

public:
	MBCSharedMemNodal(void);
	MBCSharedMemNodal(MBCType refnode_rot, unsigned nodes,
        	bool labels, MBCType rot, bool accels, const std::string &shmName);
	virtual ~MBCSharedMemNodal(void);

	MBCSharedMemBase::Type GetType(void) const;

	int Initialize(MBCType refnode_rot, unsigned nodes,
        	bool labels, MBCType rot, bool accels, const std::string &shmName);

	virtual int Negotiate(void);
	virtual int PutForces(bool bConverged);
	virtual int GetMotion(void);
	int Close(void) const;

//	uint32_t KinematicsLabel(void) const;
//	const double& X(uint8_t idx) const;
//	const double& R(uint8_t ir, uint8_t ic) const;
//	const double& Theta(uint8_t idx) const;
//	const double& Euler123(uint8_t idx) const;
//	const double& XP(uint8_t idx) const;
//	const double& Omega(uint8_t idx) const;
//	const double& XPP(uint8_t idx) const;
//	const double& OmegaP(uint8_t idx) const;
//	const uint32_t& DynamicsLabel(void) const;
//	uint32_t& DynamicsLabel(void);
//	const double& F(uint8_t idx) const;
//	double& F(uint8_t idx);
//	const double& M(uint8_t idx) const;
//	double& M(uint8_t idx);

	uint32_t GetNodes(void) const;

	const std::vector<uint32_t> *GetKinematicsLabels(void) const;

	const std::vector<double> * GetX(void) const;
	const std::vector<double> * GetR(void) const;
	const std::vector<double> * GetTheta(void) const;
	const std::vector<double> * GetEuler123(void) const;
	const std::vector<double> * GetXP(void) const;
	const std::vector<double> * GetOmega(void) const;
	const std::vector<double> * GetXPP(void) const;
	const std::vector<double> * GetOmegaP(void) const;
//
//	const double& X(uint32_t n, uint8_t idx) const;
//	const double& R(uint32_t n, uint8_t ir, uint8_t ic) const;
//	const double& Theta(uint32_t n, uint8_t idx) const;
//	const double& Euler123(uint32_t n, uint8_t idx) const;
//	const double& XP(uint32_t n, uint8_t idx) const;
//	const double& Omega(uint32_t n, uint8_t idx) const;
//	const double& XPP(uint32_t n, uint8_t idx) const;
//	const double& OmegaP(uint32_t n, uint8_t idx) const;

	uint32_t GetKinematicsLabel(uint32_t n) const;
	uint32_t GetDynamicsLabel(uint32_t n) const;

//	uint32_t KinematicsLabel(uint32_t n) const;

	std::vector<double> GetX(uint32_t n) const;
	std::vector<double> GetR(uint32_t n) const;
	std::vector<double> GetTheta(uint32_t n) const;
	std::vector<double> GetEuler123(uint32_t n) const;
	std::vector<double> GetXP(uint32_t n) const;
	std::vector<double> GetOmega(uint32_t n) const;
	std::vector<double> GetXPP(uint32_t n) const;
	std::vector<double> GetOmegaP(uint32_t n) const;

	double GetX(uint32_t n, uint8_t idx) const;
	double GetR(uint32_t n, uint8_t ir, uint8_t ic) const;
	double GetTheta(uint32_t n, uint8_t idx) const;
	double GetEuler123(uint32_t n, uint8_t idx) const;
	double GetXP(uint32_t n, uint8_t idx) const;
	double GetOmega(uint32_t n, uint8_t idx) const;
	double GetXPP(uint32_t n, uint8_t idx) const;
	double GetOmegaP(uint32_t n, uint8_t idx) const;

	const std::vector<uint32_t> * GetDynamicsLabels(void) const;

	// get all forces, returns pointer to private data member
	const std::vector<double> *GetF(void) const;
	// get all moments, returns pointer to private data member
	const std::vector<double> *GetM(void) const;

	const uint32_t& DynamicsLabel(uint32_t n) const;
//	uint32_t& DynamicsLabel(uint32_t n);

	// get vector of three forces for a node
	const std::vector<double> GetF(uint32_t n) const;
	// get vector of three moments for a node
	const std::vector<double> GetM(uint32_t n) const;

	// get force for a node in the given index
	double GetF(unsigned n, unsigned idx) const;
	// get moments for a node in the given index
	double GetM(unsigned n, unsigned idx) const;

	// set forces for a node
	void SetF(uint32_t n, std::vector<double> fxyz);
	void SetM(uint32_t n, std::vector<double> m123);

	// set forces for a node in given index
	void SetF(uint32_t n, unsigned idx, double fidx);
	// set moments for a node in given index
	void SetM(uint32_t n, unsigned idx, double midx);

//	const double& F(uint32_t n, uint8_t idx) const;
//	double& F(uint32_t n, uint8_t idx);
//	const double& M(uint32_t n, uint8_t idx) const;
//	double& M(uint32_t n, uint8_t idx);
};

//class MBCModal : public MBCSharedMemBase {
//private:
//	mutable mbc_modal_t mbc;
//
//	virtual mbc_t *GetBasePtr(void) const;
//	virtual mbc_refnode_stub_t *GetRefNodePtr(void) const;
//
//public:
//	MBCModal(void);
//	MBCModal(MBCType refnode_rot, unsigned modes);
//	virtual ~MBCModal(void);
//
//	MBCSharedMemBase::Type GetType(void) const;
//
//	int Initialize(MBCType refnode_rot, unsigned modes);
//
//	virtual int Negotiate(void) const;
//	virtual int PutForces(bool bConverged) const;
//	virtual int GetMotion(void) const;
//	int Close(void) const;
//
//	uint32_t KinematicsLabel(void) const;
//	const double& X(uint8_t idx) const;
//	const double& R(uint8_t ir, uint8_t ic) const;
//	const double& Theta(uint8_t idx) const;
//	const double& Euler123(uint8_t idx) const;
//	const double& XP(uint8_t idx) const;
//	const double& Omega(uint8_t idx) const;
//	const double& XPP(uint8_t idx) const;
//	const double& OmegaP(uint8_t idx) const;
//	const uint32_t& DynamicsLabel(void) const;
//	uint32_t& DynamicsLabel(void);
//	const double& F(uint8_t idx) const;
//	double& F(uint8_t idx);
//	const double& M(uint8_t idx) const;
//	double& M(uint8_t idx);
//
//	uint32_t GetModes(void) const;
//
//	std::vector<double> GetQ(void) const;
//	std::vector<double> GetQP(void) const;
//
//	const double& Q(uint32_t m) const;
//	const double& QP(uint32_t m) const;
//
//	const double *GetP(void) const;
//
//	const double& P(uint32_t m) const;
//	double& P(uint32_t m);
//};

#endif // MBCXXSHARED_H

