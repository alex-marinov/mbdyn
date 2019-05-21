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

#ifndef MBCXX_H
#define MBCXX_H

#include <mbc.h>

// hack...
extern "C" {
struct mbc_refnode_stub_t {
	mbc_t mbc;
	mbc_rigid_t mbcr;
};
}

class MBCBase {
protected:
	virtual mbc_t *GetBasePtr(void) const = 0;
	virtual mbc_refnode_stub_t *GetRefNodePtr(void) const = 0;

public:
	enum Type {
		MODAL = MBC_MODAL,
		NODAL = MBC_NODAL
	};

	enum Rot {
		NONE = MBC_ROT_NONE,
		THETA = MBC_ROT_THETA,
		MAT = MBC_ROT_MAT,
		EULER_123 = MBC_ROT_EULER_123
	};

	virtual MBCBase::Type GetType(void) const = 0;

	MBCBase::Rot GetRot(void) const;
	bool bRefNode(void) const;
	MBCBase::Rot GetRefNodeRot(void) const;
	bool bAccelerations(void) const;
	bool bLabels(void) const;

	void SetTimeout(int t);
	void SetVerbose(bool bv);
	void SetDataAndNext(bool bd);

	bool bVerbose(void) const;
	bool bDataAndNext(void) const;

protected:
	enum Status {
		NOT_READY,
		INITIALIZED,
		SOCKET_READY,
		READY,
		CLOSED
	} m_status;


	void SetStatus(Status s);

public:
	MBCBase(void);
	virtual ~MBCBase(void);

	int Init(const char *const path);
	int Init(const char *const host, short unsigned port);

	Status GetStatus(void) const;
	virtual int Negotiate(void) const = 0;
	virtual int PutForces(bool bConverged) const = 0;
	virtual int GetMotion(void) const = 0;
	virtual int Close(void) const = 0;

	int GetCmd(void) const;

	uint32_t GetRefNodeKinematicsLabel(void) const;

	uint32_t KinematicsLabel(void) const;

	const double * GetRefNodeX(void) const;
	const double * GetRefNodeR(void) const;
	const double * GetRefNodeTheta(void) const;
	const double * GetRefNodeEuler123(void) const;
	const double * GetRefNodeXP(void) const;
	const double * GetRefNodeOmega(void) const;
	const double * GetRefNodeXPP(void) const;
	const double * GetRefNodeOmegaP(void) const;

	const double& X(uint8_t idx) const;
	const double& R(uint8_t ir, uint8_t ic) const;
	const double& Theta(uint8_t idx) const;
	const double& Euler123(uint8_t idx) const;
	const double& XP(uint8_t idx) const;
	const double& Omega(uint8_t idx) const;
	const double& XPP(uint8_t idx) const;
	const double& OmegaP(uint8_t idx) const;

	uint32_t GetRefNodeDynamicsLabel(void) const;
	const uint32_t& DynamicsLabel(void) const;
	uint32_t& DynamicsLabel(void);

	const double *GetRefNodeF(void) const;
	const double *GetRefNodeM(void) const;

	const double& F(uint8_t idx) const;
	double& F(uint8_t idx);
	const double& M(uint8_t idx) const;
	double& M(uint8_t idx);
};

class MBCNodal : public MBCBase {
private:
	mutable mbc_nodal_t mbc;

	virtual mbc_t *GetBasePtr(void) const;
	virtual mbc_refnode_stub_t *GetRefNodePtr(void) const;

public:
	MBCNodal(void);
	MBCNodal(MBCBase::Rot refnode_rot, unsigned nodes,
        	bool labels, MBCBase::Rot rot, bool accels);
	virtual ~MBCNodal(void);

	MBCBase::Type GetType(void) const;

	int Initialize(MBCBase::Rot refnode_rot, unsigned nodes,
        	bool labels, MBCBase::Rot rot, bool accels);

	virtual int Negotiate(void) const;
	virtual int PutForces(bool bConverged) const;
	virtual int GetMotion(void) const;
	int Close(void) const;

	uint32_t KinematicsLabel(void) const;
	const double& X(uint8_t idx) const;
	const double& R(uint8_t ir, uint8_t ic) const;
	const double& Theta(uint8_t idx) const;
	const double& Euler123(uint8_t idx) const;
	const double& XP(uint8_t idx) const;
	const double& Omega(uint8_t idx) const;
	const double& XPP(uint8_t idx) const;
	const double& OmegaP(uint8_t idx) const;
	const uint32_t& DynamicsLabel(void) const;
	uint32_t& DynamicsLabel(void);
	const double& F(uint8_t idx) const;
	double& F(uint8_t idx);
	const double& M(uint8_t idx) const;
	double& M(uint8_t idx);

	uint32_t GetNodes(void) const;

	uint32_t *GetKinematicsLabel(void) const;

	const double * GetX(void) const;
	const double * GetR(void) const;
	const double * GetTheta(void) const;
	const double * GetEuler123(void) const;
	const double * GetXP(void) const;
	const double * GetOmega(void) const;
	const double * GetXPP(void) const;
	const double * GetOmegaP(void) const;

	const double& X(uint32_t n, uint8_t idx) const;
	const double& R(uint32_t n, uint8_t ir, uint8_t ic) const;
	const double& Theta(uint32_t n, uint8_t idx) const;
	const double& Euler123(uint32_t n, uint8_t idx) const;
	const double& XP(uint32_t n, uint8_t idx) const;
	const double& Omega(uint32_t n, uint8_t idx) const;
	const double& XPP(uint32_t n, uint8_t idx) const;
	const double& OmegaP(uint32_t n, uint8_t idx) const;

	uint32_t GetKinematicsLabel(uint32_t n) const;
	uint32_t KinematicsLabel(uint32_t n) const;

	const double * GetX(uint32_t n) const;
	const double * GetR(uint32_t n) const;
	const double * GetTheta(uint32_t n) const;
	const double * GetEuler123(uint32_t n) const;
	const double * GetXP(uint32_t n) const;
	const double * GetOmega(uint32_t n) const;
	const double * GetXPP(uint32_t n) const;
	const double * GetOmegaP(uint32_t n) const;

	uint32_t *GetDynamicsLabel(void) const;

	const double *GetF(void) const;
	const double *GetM(void) const;

	const uint32_t& DynamicsLabel(uint32_t n) const;
	uint32_t& DynamicsLabel(uint32_t n);

	const double *GetF(uint32_t n) const;
	const double *GetM(uint32_t n) const;

	const double& F(uint32_t n, uint8_t idx) const;
	double& F(uint32_t n, uint8_t idx);
	const double& M(uint32_t n, uint8_t idx) const;
	double& M(uint32_t n, uint8_t idx);
};

class MBCModal : public MBCBase {
private:
	mutable mbc_modal_t mbc;

	virtual mbc_t *GetBasePtr(void) const;
	virtual mbc_refnode_stub_t *GetRefNodePtr(void) const;

public:
	MBCModal(void);
	MBCModal(MBCBase::Rot refnode_rot, unsigned modes);
	virtual ~MBCModal(void);

	MBCBase::Type GetType(void) const;

	int Initialize(MBCBase::Rot refnode_rot, unsigned modes);

	virtual int Negotiate(void) const;
	virtual int PutForces(bool bConverged) const;
	virtual int GetMotion(void) const;
	int Close(void) const;

	uint32_t KinematicsLabel(void) const;
	const double& X(uint8_t idx) const;
	const double& R(uint8_t ir, uint8_t ic) const;
	const double& Theta(uint8_t idx) const;
	const double& Euler123(uint8_t idx) const;
	const double& XP(uint8_t idx) const;
	const double& Omega(uint8_t idx) const;
	const double& XPP(uint8_t idx) const;
	const double& OmegaP(uint8_t idx) const;
	const uint32_t& DynamicsLabel(void) const;
	uint32_t& DynamicsLabel(void);
	const double& F(uint8_t idx) const;
	double& F(uint8_t idx);
	const double& M(uint8_t idx) const;
	double& M(uint8_t idx);

	uint32_t GetModes(void) const;

	const double * GetQ(void) const;
	const double * GetQP(void) const;

	const double& Q(uint32_t m) const;
	const double& QP(uint32_t m) const;

	const double *GetP(void) const;

	const double& P(uint32_t m) const;
	double& P(uint32_t m);
};

#endif // MBCXX_H

