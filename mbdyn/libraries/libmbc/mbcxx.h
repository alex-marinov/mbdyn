/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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
struct mbc_rigid_stub_t {
	mbc_t mbc;
	mbc_rigid_t mbcr;
};
}

class MBCBase {
protected:
	virtual mbc_t *GetBasePtr(void) const = 0;
	virtual mbc_rigid_stub_t *GetRigidPtr(void) const = 0;

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

	MBCBase::Rot GetRot(void) const { return MBCBase::Rot(MBC_F_ROT(GetRigidPtr())); };
	bool bRefNode(void) const { return MBC_F_REF_NODE(GetRigidPtr()); };
	MBCBase::Rot GetRefNodeRot(void) const { return MBCBase::Rot(MBC_F_REF_NODE_ROT(GetRigidPtr())); };
	bool bAccelerations(void) const { return MBC_F_ACCELS(GetRigidPtr()); };
	bool bLabels(void) const { return MBC_F_LABELS(GetRigidPtr()); };

	void SetTimeout(int t) { GetBasePtr()->timeout = t; };
	void SetVerbose(bool bv) { GetBasePtr()->verbose = bv; };
	void SetDataAndNext(bool bd) { GetBasePtr()->data_and_next = bd; };

	bool bVerbose(void) const { return GetBasePtr()->verbose; };
	bool bDataAndNext(void) const { return GetBasePtr()->data_and_next; };

protected:
	enum Status {
		NOT_READY,
		READY,
		CLOSED
	} m_status;

	Status GetStatus(void) const { return m_status; };
	void SetStatus(Status s) { m_status = s; };

public:
	MBCBase(void) : m_status(NOT_READY) {};
	virtual ~MBCBase(void) {};

	int Init(const std::string& path) {
		if (m_status != NOT_READY) return -1;
		int rc = mbc_unix_init(GetBasePtr(), path.c_str());
		if (rc == 0) m_status = READY;
		return rc;
	};

	int Init(const std::string& host, short unsigned port) {
		if (m_status != NOT_READY) return -1;
		int rc = mbc_inet_init(GetBasePtr(), host.c_str(), port);
		if (rc == 0) m_status = READY;
		return rc;
	};

	virtual int Negotiate(void) const = 0;
	virtual int PutForces(bool bConverged) const = 0;
	virtual int GetMotion(void) const = 0;
	virtual int Close(void) const = 0;

	uint32_t GetRefNodeKinematicsLabel(void) const { return MBC_R_K_LABEL(GetRigidPtr()); };

	uint32_t KinematicsLabel(void) const { return MBC_R_K_LABEL(GetRigidPtr()); };

	const double *const GetRefNodeX(void) const { return MBC_R_X(GetRigidPtr()); };
	const double *const GetRefNodeR(void) const { return MBC_R_R(GetRigidPtr()); };
	const double *const GetRefNodeTheta(void) const { return MBC_R_THETA(GetRigidPtr()); };
	const double *const GetRefNodeEuler123(void) const { return MBC_R_EULER_123(GetRigidPtr()); };
	const double *const GetRefNodeXP(void) const { return MBC_R_XP(GetRigidPtr()); };
	const double *const GetRefNodeOmega(void) const { return MBC_R_OMEGA(GetRigidPtr()); };
	const double *const GetRefNodeXPP(void) const { return MBC_R_XPP(GetRigidPtr()); };
	const double *const GetRefNodeOmegaP(void) const { return MBC_R_OMEGAP(GetRigidPtr()); };

	const double& X(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_X(GetRigidPtr()))[idx - 1];
	};

	const double& R(uint8_t ir, uint8_t ic) const {
		if (ir < 1 || ir > 3 || ic < 1 || ic > 3) throw;
		return (MBC_R_R(GetRigidPtr()))[3*(ic - 1) + ir - 1];
	};

	const double& Theta(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_THETA(GetRigidPtr()))[idx - 1];
	};

	const double& Euler123(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_EULER_123(GetRigidPtr()))[idx - 1];
	};

	const double& XP(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_XP(GetRigidPtr()))[idx - 1];
	};

	const double& Omega(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_OMEGA(GetRigidPtr()))[idx - 1];
	};

	const double& XPP(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_XPP(GetRigidPtr()))[idx - 1];
	};

	const double& OmegaP(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_OMEGAP(GetRigidPtr()))[idx - 1];
	};

	uint32_t GetRefNodeDynamicsLabel(void) const { return MBC_R_D_LABEL(GetRigidPtr()); };
	const uint32_t& DynamicsLabel(void) const { return MBC_R_D_LABEL(GetRigidPtr()); };
	uint32_t& DynamicsLabel(void) { return MBC_R_D_LABEL(GetRigidPtr()); };

	const double *GetRefNodeF(void) const { return MBC_R_F(GetRigidPtr()); };
	const double *GetRefNodeM(void) const { return MBC_R_M(GetRigidPtr()); };

	const double& F(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_F(GetRigidPtr()))[idx - 1];
	};

	double& F(uint8_t idx) {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_F(GetRigidPtr()))[idx - 1];
	};

	const double& M(uint8_t idx) const {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_M(GetRigidPtr()))[idx - 1];
	};

	double& M(uint8_t idx) {
		if (idx < 1 || idx > 3) throw;
		return (MBC_R_M(GetRigidPtr()))[idx - 1];
	};
};

class MBCNodal : public MBCBase {
private:
	mutable mbc_nodal_t mbc;

	virtual mbc_t *GetBasePtr(void) const { return (mbc_t*)&mbc; };
	virtual mbc_rigid_stub_t *GetRigidPtr(void) const { return (mbc_rigid_stub_t *)&mbc; };

public:
	MBCNodal(MBCBase::Rot rigid, unsigned nodes,
        	bool labels, MBCBase::Rot rot, bool accels)
	{
		if (mbc_nodal_init(&mbc, rigid, nodes, labels,
			unsigned(rot) | MBC_U_ROT_REF_NODE(unsigned(rigid)), accels))
		{
			throw;
		}
	};

	virtual ~MBCNodal(void) { Close(); };

	MBCBase::Type GetType(void) const { return NODAL; };

	virtual int Negotiate(void) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_nodal_negotiate_request(&mbc);
	};

	virtual int PutForces(bool bConverged) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_nodal_put_forces(&mbc, bConverged);
	};

	virtual int GetMotion(void) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_nodal_get_motion(&mbc);
	};

	int Close(void) const {
		int rc = -1;
		if (GetStatus() == READY) {
			rc = mbc_nodal_destroy(&mbc);
			const_cast<MBCNodal *>(this)->SetStatus(CLOSED);
		}
		return rc;
	};


	uint32_t KinematicsLabel(void) const { return MBCBase::KinematicsLabel(); };
	const double& X(uint8_t idx) const { return MBCBase::X(idx); };
	const double& R(uint8_t ir, uint8_t ic) const { return MBCBase::R(ir, ic); };
	const double& Theta(uint8_t idx) const { return MBCBase::Theta(idx); };
	const double& Euler123(uint8_t idx) const { return MBCBase::Euler123(idx); };
	const double& XP(uint8_t idx) const { return MBCBase::XP(idx); };
	const double& Omega(uint8_t idx) const { return MBCBase::Omega(idx); };
	const double& XPP(uint8_t idx) const { return MBCBase::XPP(idx); };
	const double& OmegaP(uint8_t idx) const { return MBCBase::OmegaP(idx); };
	const uint32_t& DynamicsLabel(void) const {  return MBCBase::DynamicsLabel(); };
	uint32_t& DynamicsLabel(void) { return MBCBase::DynamicsLabel(); };
	const double& F(uint8_t idx) const { return MBCBase::F(idx); };
	double& F(uint8_t idx) { return MBCBase::F(idx); };
	const double& M(uint8_t idx) const { return MBCBase::M(idx); };
	double& M(uint8_t idx) { return MBCBase::M(idx); };

	uint32_t GetNodes(void) const { return mbc.nodes; };

	uint32_t *GetKinematicsLabel(void) const { return MBC_N_K_LABELS(&mbc); };

	const double *const GetX(void) const { return MBC_N_X(&mbc); };
	const double *const GetR(void) const { return MBC_N_R(&mbc); };
	const double *const GetTheta(void) const { return MBC_N_THETA(&mbc); };
	const double *const GetEuler123(void) const { return MBC_N_EULER_123(&mbc); };
	const double *const GetXP(void) const { return MBC_N_XP(&mbc); };
	const double *const GetOmega(void) const { return MBC_N_OMEGA(&mbc); };
	const double *const GetXPP(void) const { return MBC_N_XPP(&mbc); };
	const double *const GetOmegaP(void) const { return MBC_N_OMEGAP(&mbc); };

	const double& X(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_X(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& R(uint32_t n, uint8_t ir, uint8_t ic) const {
		if (ir < 1 || ir > 3 || ic < 1 || ic > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_R(&mbc))[9*(n - 1) + 3*(ic - 1) + ir - 1];
	};

	const double& Theta(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_THETA(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& Euler123(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_EULER_123(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& XP(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_XP(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& Omega(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_OMEGA(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& XPP(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_XPP(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& OmegaP(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_OMEGAP(&mbc))[3*(n - 1) + idx - 1];
	};

	uint32_t GetKinematicsLabel(uint32_t n) const {
		if (n < 1 || n > GetNodes()) throw;
		return (MBC_N_K_LABELS(&mbc))[n - 1];
	};

	uint32_t KinematicsLabel(uint32_t n) const {
		if (n < 1 || n > GetNodes()) throw;
		return (MBC_N_K_LABELS(&mbc))[n - 1];
	};

	const double *const GetX(uint32_t n) const { return &(MBC_N_X(&mbc)[3*(n - 1)]); };
	const double *const GetR(uint32_t n) const { return &(MBC_N_R(&mbc)[9*(n - 1)]); };
	const double *const GetTheta(uint32_t n) const { return &(MBC_N_THETA(&mbc)[3*(n - 1)]); };
	const double *const GetEuler123(uint32_t n) const { return &(MBC_N_EULER_123(&mbc)[3*(n - 1)]); };
	const double *const GetXP(uint32_t n) const { return &(MBC_N_XP(&mbc)[3*(n - 1)]); };
	const double *const GetOmega(uint32_t n) const { return &(MBC_N_OMEGA(&mbc)[3*(n - 1)]); };
	const double *const GetXPP(uint32_t n) const { return &(MBC_N_XPP(&mbc)[3*(n - 1)]); };
	const double *const GetOmegaP(uint32_t n) const { return &(MBC_N_OMEGAP(&mbc)[3*(n - 1)]); };

	uint32_t *GetDynamicsLabel(void) const { return MBC_N_D_LABELS(&mbc); };

	const uint32_t& DynamicsLabel(uint32_t n) const {
		if (n < 1 || n > GetNodes()) throw;
		return (MBC_N_D_LABELS(&mbc))[n - 1];
	};

	uint32_t& DynamicsLabel(uint32_t n) {
		if (n < 1 || n > GetNodes()) throw;
		return (MBC_N_D_LABELS(&mbc))[n - 1];
	};

	const double *GetF(void) const { return MBC_N_F(&mbc); };
	const double *GetM(void) const { return MBC_N_M(&mbc); };

	const double *GetF(uint32_t n) const {
		if (n < 1 || n > GetNodes()) throw;
		return &(MBC_N_F(&mbc))[3*(n - 1)];
	};

	const double *GetM(uint32_t n) const {
		if (n < 1 || n > GetNodes()) throw;
		return &(MBC_N_M(&mbc))[3*(n - 1)];
	};

	const double& F(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
	};

	double& F(uint32_t n, uint8_t idx) {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
	};

	const double& M(uint32_t n, uint8_t idx) const {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
	};

	double& M(uint32_t n, uint8_t idx) {
		if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
		return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
	};
};

class MBCModal : public MBCBase {
private:
	mutable mbc_modal_t mbc;

	virtual mbc_t *GetBasePtr(void) const { return (mbc_t*)&mbc; };
	virtual mbc_rigid_stub_t *GetRigidPtr(void) const { return (mbc_rigid_stub_t *)&mbc; };

public:
	MBCModal(MBCBase::Rot rigid, unsigned modes)
	{
		if (mbc_modal_init(&mbc, rigid, modes))
		{
			throw;
		}
	};

	virtual ~MBCModal(void) { Close(); };

	MBCBase::Type GetType(void) const { return MODAL; };

	virtual int Negotiate(void) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_modal_negotiate_request(&mbc);
	};

	virtual int PutForces(bool bConverged) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_modal_put_forces(&mbc, bConverged);
	};

	virtual int GetMotion(void) const {
		if (GetStatus() != READY) {
			return -1;
		}
		return mbc_modal_get_motion(&mbc);
	};

	int Close(void) const {
		int rc = -1;
		if (GetStatus() == READY) {
			rc = mbc_modal_destroy(&mbc);
			const_cast<MBCModal *>(this)->SetStatus(CLOSED);
		}
		return rc;
	};

	uint32_t KinematicsLabel(void) const { return MBCBase::KinematicsLabel(); };
	const double& X(uint8_t idx) const { return MBCBase::X(idx); };
	const double& R(uint8_t ir, uint8_t ic) const { return MBCBase::R(ir, ic); };
	const double& Theta(uint8_t idx) const { return MBCBase::Theta(idx); };
	const double& Euler123(uint8_t idx) const { return MBCBase::Euler123(idx); };
	const double& XP(uint8_t idx) const { return MBCBase::XP(idx); };
	const double& Omega(uint8_t idx) const { return MBCBase::Omega(idx); };
	const double& XPP(uint8_t idx) const { return MBCBase::XPP(idx); };
	const double& OmegaP(uint8_t idx) const { return MBCBase::OmegaP(idx); };
	const uint32_t& DynamicsLabel(void) const {  return MBCBase::DynamicsLabel(); };
	uint32_t& DynamicsLabel(void) { return MBCBase::DynamicsLabel(); };
	const double& F(uint8_t idx) const { return MBCBase::F(idx); };
	double& F(uint8_t idx) { return MBCBase::F(idx); };
	const double& M(uint8_t idx) const { return MBCBase::M(idx); };
	double& M(uint8_t idx) { return MBCBase::M(idx); };

	uint32_t GetModes(void) const { return mbc.modes; };

	const double *const GetQ(void) const { return MBC_Q(&mbc); };
	const double *const GetQP(void) const { return MBC_QP(&mbc); };

	const double& Q(uint32_t m) const {
		if (m < 1 || m > GetModes()) throw;
		return (MBC_Q(&mbc))[m - 1];
	};

	const double& QP(uint32_t m) const {
		if (m < 1 || m > GetModes()) throw;
		return (MBC_QP(&mbc))[m - 1];
	};

	const double *GetP(void) const { return MBC_P(&mbc); };

	const double& P(uint32_t m) const {
		if (m < 1 || m > GetModes()) throw;
		return (MBC_P(&mbc))[m - 1];
	};

	double& P(uint32_t m) {
		if (m < 1 || m > GetModes()) throw;
		return (MBC_P(&mbc))[m - 1];
	};
};

#endif // MBCXX_H

