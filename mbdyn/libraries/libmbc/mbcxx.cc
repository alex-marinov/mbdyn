/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifdef USE_SOCKET

#include <stdlib.h>
#include <assert.h>
#include <cstring>
#include "mbcxx.h"

MBCBase::Rot
MBCBase::GetRot(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return MBCBase::Rot(MBC_F_ROT(GetRefNodePtr()));
}

bool
MBCBase::bRefNode(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return MBC_F_REF_NODE(GetRefNodePtr());
}

MBCBase::Rot
MBCBase::GetRefNodeRot(void) const {
	assert(GetStatus() >= INITIALIZED);
	return MBCBase::Rot(MBC_F_ROT_REF_NODE(GetRefNodePtr()));
}

bool
MBCBase::bAccelerations(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return MBC_F_ACCELS(GetRefNodePtr());
}

bool
MBCBase::bLabels(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return MBC_F_LABELS(GetRefNodePtr());
}

void
MBCBase::SetTimeout(int t)
{
	GetBasePtr()->timeout = t;
}

void
MBCBase::SetVerbose(bool bv)
{
	GetBasePtr()->verbose = bv;
}

void
MBCBase::SetDataAndNext(bool bd)
{
	GetBasePtr()->data_and_next = bd;
}

bool
MBCBase::bVerbose(void) const
{
	return GetBasePtr()->verbose;
}

bool
MBCBase::bDataAndNext(void) const
{
	return GetBasePtr()->data_and_next;
}

MBCBase::Status
MBCBase::GetStatus(void) const
{
	return m_status;
}

void
MBCBase::SetStatus(MBCBase::Status s)
{
	m_status = s;
}

MBCBase::MBCBase(void)
: m_status(NOT_READY)
{
}

MBCBase::~MBCBase(void)
{
}

int
MBCBase::Init(const char *const path)
{
	if (GetStatus() != INITIALIZED) return -1;
	int rc = mbc_unix_init(GetBasePtr(), path);
	if (rc == 0) SetStatus(SOCKET_READY);
	return rc;
}

int
MBCBase::Init(const char *const host, short unsigned port)
{
	if (GetStatus() != INITIALIZED) return -1;
	int rc = mbc_inet_init(GetBasePtr(), host, port);
	if (rc == 0) SetStatus(SOCKET_READY);
	return rc;
}

int
MBCBase::GetCmd(void) const
{
	return GetBasePtr()->cmd;
}

uint32_t
MBCBase::GetRefNodeKinematicsLabel(void) const
{
	assert(GetStatus() == READY);
	return MBC_R_K_LABEL(GetRefNodePtr());
}

uint32_t
MBCBase::KinematicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	return MBC_R_K_LABEL(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeX(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_X(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeR(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MAT);
	return MBC_R_R(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeTheta(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == THETA);
	return MBC_R_THETA(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeEuler123(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == EULER_123);
	return MBC_R_EULER_123(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeXP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_XP(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeOmega(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	return MBC_R_OMEGA(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeXPP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(bAccelerations());
	return MBC_R_XPP(GetRefNodePtr());
}

const double *const
MBCBase::GetRefNodeOmegaP(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	assert(bAccelerations());
	return MBC_R_OMEGAP(GetRefNodePtr());
}

const double&
MBCBase::X(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_X(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::R(uint8_t ir, uint8_t ic) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == MAT);
	if (ir < 1 || ir > 3 || ic < 1 || ic > 3) throw;
	return (MBC_R_R(GetRefNodePtr()))[3*(ic - 1) + ir - 1];
}

const double&
MBCBase::Theta(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == THETA);
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_THETA(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::Euler123(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() == EULER_123);
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_EULER_123(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::XP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_XP(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::Omega(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_OMEGA(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::XPP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(bAccelerations());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_XPP(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::OmegaP(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	assert(bAccelerations());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_OMEGAP(GetRefNodePtr()))[idx - 1];
}

uint32_t
MBCBase::GetRefNodeDynamicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_D_LABEL(GetRefNodePtr());
}

const uint32_t&
MBCBase::DynamicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_D_LABEL(GetRefNodePtr());
}

uint32_t&
MBCBase::DynamicsLabel(void)
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_D_LABEL(GetRefNodePtr());
}

const double *
MBCBase::GetRefNodeF(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	return MBC_R_F(GetRefNodePtr());
}

const double *
MBCBase::GetRefNodeM(void) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	return MBC_R_M(GetRefNodePtr());
}

const double&
MBCBase::F(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_F(GetRefNodePtr()))[idx - 1];
}

double&
MBCBase::F(uint8_t idx)
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_F(GetRefNodePtr()))[idx - 1];
}

const double&
MBCBase::M(uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_M(GetRefNodePtr()))[idx - 1];
}

double&
MBCBase::M(uint8_t idx)
{
	assert(GetStatus() == READY);
	assert(bRefNode());
	assert(GetRefNodeRot() != NONE);
	if (idx < 1 || idx > 3) throw;
	return (MBC_R_M(GetRefNodePtr()))[idx - 1];
}

mbc_t *
MBCNodal::GetBasePtr(void) const
{
	return (mbc_t*)&mbc;
}

mbc_refnode_stub_t *
MBCNodal::GetRefNodePtr(void) const
{
	return (mbc_refnode_stub_t *)&mbc;
}

MBCNodal::MBCNodal(void)
{
	memset(&mbc, 0, sizeof(mbc));
}

MBCNodal::MBCNodal(MBCBase::Rot refnode_rot, unsigned nodes,
       	bool labels, MBCBase::Rot rot, bool accels)
{
	if (Initialize(refnode_rot, nodes, labels, rot, accels)) {
		throw;
	}
}

MBCNodal::~MBCNodal(void)
{
	Close();
}

MBCBase::Type
MBCNodal::GetType(void) const
{
	return NODAL;
}

int
MBCNodal::Initialize(MBCBase::Rot refnode_rot, unsigned nodes,
       	bool labels, MBCBase::Rot rot, bool accels)
{
	if (GetStatus() != NOT_READY) return -1;
	memset(&mbc, 0, sizeof(mbc));
	int rc = mbc_nodal_init(&mbc, refnode_rot, nodes, labels,
		unsigned(rot) | MBC_U_ROT_2_REF_NODE_ROT(unsigned(refnode_rot)),
		accels);
	if (rc == 0) const_cast<MBCNodal *>(this)->SetStatus(INITIALIZED);
	return rc;
}

int
MBCNodal::Negotiate(void) const
{
	if (GetStatus() != SOCKET_READY) return -1;
	int rc = mbc_nodal_negotiate_request(&mbc);
	if (rc == 0) const_cast<MBCNodal *>(this)->SetStatus(READY);
	return rc;
}

int
MBCNodal::PutForces(bool bConverged) const
{
	if (GetStatus() != READY) return -1;
	return mbc_nodal_put_forces(&mbc, bConverged);
}

int
MBCNodal::GetMotion(void) const
{
	if (GetStatus() != READY) return -1;
	return mbc_nodal_get_motion(&mbc);
}

int
MBCNodal::Close(void) const
{
	int rc = -1;
	if (GetStatus() == READY) {
		rc = mbc_nodal_destroy(&mbc);
		const_cast<MBCNodal *>(this)->SetStatus(CLOSED);
	}
	return rc;
}

uint32_t
MBCNodal::KinematicsLabel(void) const
{
	return MBCBase::KinematicsLabel();
}

const double&
MBCNodal::X(uint8_t idx) const
{
	return MBCBase::X(idx);
}

const double&
MBCNodal::R(uint8_t ir, uint8_t ic) const
{
	return MBCBase::R(ir, ic);
}

const double&
MBCNodal::Theta(uint8_t idx) const
{
	return MBCBase::Theta(idx);
}

const double&
MBCNodal::Euler123(uint8_t idx) const
{
	return MBCBase::Euler123(idx);
}

const double&
MBCNodal::XP(uint8_t idx) const
{
	return MBCBase::XP(idx);
}

const double&
MBCNodal::Omega(uint8_t idx) const
{
	return MBCBase::Omega(idx);
}

const double&
MBCNodal::XPP(uint8_t idx) const
{
	return MBCBase::XPP(idx);
}

const double&
MBCNodal::OmegaP(uint8_t idx) const
{
	return MBCBase::OmegaP(idx);
}

const uint32_t&
MBCNodal::DynamicsLabel(void) const
{
	return MBCBase::DynamicsLabel();
}

uint32_t&
MBCNodal::DynamicsLabel(void)
{
	return MBCBase::DynamicsLabel();
}

const double&
MBCNodal::F(uint8_t idx) const
{
	return MBCBase::F(idx);
}

double&
MBCNodal::F(uint8_t idx)
{
	return MBCBase::F(idx);
}

const double&
MBCNodal::M(uint8_t idx) const
{
	return MBCBase::M(idx);
}

double&
MBCNodal::M(uint8_t idx)
{
	return MBCBase::M(idx);
}

uint32_t
MBCNodal::GetNodes(void) const
{
	assert(GetStatus() >= INITIALIZED);
	return mbc.nodes;
}

uint32_t *
MBCNodal::GetKinematicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	return MBC_N_K_LABELS(&mbc);
}

const double *const
MBCNodal::GetX(void) const
{
	assert(GetStatus() == READY);
	return MBC_N_X(&mbc);
}

const double *const
MBCNodal::GetR(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MAT);
	return MBC_N_R(&mbc);
}

const double *const
MBCNodal::GetTheta(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == THETA);
	return MBC_N_THETA(&mbc);
}

const double *const
MBCNodal::GetEuler123(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == EULER_123);
	return MBC_N_EULER_123(&mbc);
}

const double *const
MBCNodal::GetXP(void) const
{
	assert(GetStatus() == READY);
	return MBC_N_XP(&mbc);
}

const double *const
MBCNodal::GetOmega(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	return MBC_N_OMEGA(&mbc);
}

const double *const
MBCNodal::GetXPP(void) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	return MBC_N_XPP(&mbc);
}

const double *const
MBCNodal::GetOmegaP(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	assert(bAccelerations());
	return MBC_N_OMEGAP(&mbc);
}

const double&
MBCNodal::X(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_X(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::R(uint32_t n, uint8_t ir, uint8_t ic) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MAT);
	if (ir < 1 || ir > 3 || ic < 1 || ic > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_R(&mbc))[9*(n - 1) + 3*(ic - 1) + ir - 1];
}

const double&
MBCNodal::Theta(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == THETA);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_THETA(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::Euler123(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == EULER_123);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_EULER_123(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::XP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_XP(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::Omega(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_OMEGA(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::XPP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_XPP(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::OmegaP(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	assert(bAccelerations());
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_OMEGAP(&mbc))[3*(n - 1) + idx - 1];
}

uint32_t
MBCNodal::GetKinematicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return (MBC_N_K_LABELS(&mbc))[n - 1];
}

uint32_t
MBCNodal::KinematicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return (MBC_N_K_LABELS(&mbc))[n - 1];
}

const double *const
MBCNodal::GetX(uint32_t n) const
{
	assert(GetStatus() == READY);
	return &(MBC_N_X(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetR(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == MAT);
	return &(MBC_N_R(&mbc)[9*(n - 1)]);
}

const double *const
MBCNodal::GetTheta(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == THETA);
	return &(MBC_N_THETA(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetEuler123(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() == EULER_123);
	return &(MBC_N_EULER_123(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetXP(uint32_t n) const
{
	assert(GetStatus() == READY);
	return &(MBC_N_XP(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetOmega(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	return &(MBC_N_OMEGA(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetXPP(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bAccelerations());
	return &(MBC_N_XPP(&mbc)[3*(n - 1)]);
}

const double *const
MBCNodal::GetOmegaP(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	assert(bAccelerations());
	return &(MBC_N_OMEGAP(&mbc)[3*(n - 1)]);
}

uint32_t *
MBCNodal::GetDynamicsLabel(void) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	return MBC_N_D_LABELS(&mbc);
}

const uint32_t&
MBCNodal::DynamicsLabel(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return (MBC_N_D_LABELS(&mbc))[n - 1];
}

uint32_t&
MBCNodal::DynamicsLabel(uint32_t n)
{
	assert(GetStatus() == READY);
	assert(bLabels());
	if (n < 1 || n > GetNodes()) throw;
	return (MBC_N_D_LABELS(&mbc))[n - 1];
}

const double *
MBCNodal::GetF(void) const
{
	assert(GetStatus() == READY);
	return MBC_N_F(&mbc);
}

const double *
MBCNodal::GetM(void) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	return MBC_N_M(&mbc);
}

const double *
MBCNodal::GetF(uint32_t n) const
{
	assert(GetStatus() == READY);
	if (n < 1 || n > GetNodes()) throw;
	return &(MBC_N_F(&mbc))[3*(n - 1)];
}

const double *
MBCNodal::GetM(uint32_t n) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	if (n < 1 || n > GetNodes()) throw;
	return &(MBC_N_M(&mbc))[3*(n - 1)];
}

const double&
MBCNodal::F(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
}

double&
MBCNodal::F(uint32_t n, uint8_t idx)
{
	assert(GetStatus() == READY);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_F(&mbc))[3*(n - 1) + idx - 1];
}

const double&
MBCNodal::M(uint32_t n, uint8_t idx) const
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
}

double&
MBCNodal::M(uint32_t n, uint8_t idx)
{
	assert(GetStatus() == READY);
	assert(GetRot() != NONE);
	if (idx < 1 || idx > 3 || n < 1 || n > GetNodes()) throw;
	return (MBC_N_M(&mbc))[3*(n - 1) + idx - 1];
}

mbc_t *
MBCModal::GetBasePtr(void) const
{
	return (mbc_t*)&mbc;
}

mbc_refnode_stub_t *
MBCModal::GetRefNodePtr(void) const
{
	return (mbc_refnode_stub_t *)&mbc;
}

MBCModal::MBCModal(void)
{
	memset(&mbc, 0, sizeof(mbc));
}

MBCModal::MBCModal(MBCBase::Rot refnode_rot, unsigned modes)
{
	if (Initialize(refnode_rot, modes)) {
		throw;
	}
}

MBCModal::~MBCModal(void)
{
	Close();
}

MBCBase::Type
MBCModal::GetType(void) const
{
	return MODAL;
}

int
MBCModal::Initialize(MBCBase::Rot refnode_rot, unsigned modes)
{
	if (GetStatus() != NOT_READY) return -1;
	memset(&mbc, 0, sizeof(mbc));
	int rc = mbc_modal_init(&mbc, refnode_rot, modes);
	if (rc == 0) const_cast<MBCModal *>(this)->SetStatus(INITIALIZED);
	return rc;
}

int
MBCModal::Negotiate(void) const
{
	if (GetStatus() != SOCKET_READY) return -1;
	int rc = mbc_modal_negotiate_request(&mbc);
	if (rc == 0) const_cast<MBCModal *>(this)->SetStatus(READY);
	return rc;
}

int
MBCModal::PutForces(bool bConverged) const
{
	if (GetStatus() != READY) return -1;
	return mbc_modal_put_forces(&mbc, bConverged);
}

int
MBCModal::GetMotion(void) const
{
	if (GetStatus() != READY) return -1;
	return mbc_modal_get_motion(&mbc);
}

int
MBCModal::Close(void) const
{
	int rc = -1;
	if (GetStatus() == READY) {
		rc = mbc_modal_destroy(&mbc);
		const_cast<MBCModal *>(this)->SetStatus(CLOSED);
	}
	return rc;
}

uint32_t
MBCModal::KinematicsLabel(void) const
{
	return MBCBase::KinematicsLabel();
}

const double&
MBCModal::X(uint8_t idx) const
{
	return MBCBase::X(idx);
}

const double&
MBCModal::R(uint8_t ir, uint8_t ic) const
{
	return MBCBase::R(ir, ic);
}

const double&
MBCModal::Theta(uint8_t idx) const
{
	return MBCBase::Theta(idx);
}

const double&
MBCModal::Euler123(uint8_t idx) const
{
	return MBCBase::Euler123(idx);
}

const double&
MBCModal::XP(uint8_t idx) const
{
	return MBCBase::XP(idx);
}

const double&
MBCModal::Omega(uint8_t idx) const
{
	return MBCBase::Omega(idx);
}

const double&
MBCModal::XPP(uint8_t idx) const
{
	return MBCBase::XPP(idx);
}

const double&
MBCModal::OmegaP(uint8_t idx) const
{
	return MBCBase::OmegaP(idx);
}

const uint32_t&
MBCModal::DynamicsLabel(void) const
{
	return MBCBase::DynamicsLabel();
}

uint32_t&
MBCModal::DynamicsLabel(void)
{
	return MBCBase::DynamicsLabel();
}

const double&
MBCModal::F(uint8_t idx) const
{
	return MBCBase::F(idx);
}

double&
MBCModal::F(uint8_t idx)
{
	return MBCBase::F(idx);
}

const double&
MBCModal::M(uint8_t idx) const
{
	return MBCBase::M(idx);
}

double&
MBCModal::M(uint8_t idx)
{
	return MBCBase::M(idx);
}

uint32_t
MBCModal::GetModes(void) const
{
	assert(GetStatus() == READY);
	return mbc.modes;
}

const double *const
MBCModal::GetQ(void) const
{
	assert(GetStatus() == READY);
	return MBC_Q(&mbc);
}

const double *const
MBCModal::GetQP(void) const
{
	assert(GetStatus() == READY);
	return MBC_QP(&mbc);
}

const double&
MBCModal::Q(uint32_t m) const
{
	assert(GetStatus() == READY);
	if (m < 1 || m > GetModes()) throw;
	return (MBC_Q(&mbc))[m - 1];
}

const double&
MBCModal::QP(uint32_t m) const
{
	assert(GetStatus() == READY);
	if (m < 1 || m > GetModes()) throw;
	return (MBC_QP(&mbc))[m - 1];
}

const double *
MBCModal::GetP(void) const
{
	assert(GetStatus() == READY);
	return MBC_P(&mbc);
}

const double&
MBCModal::P(uint32_t m) const
{
	assert(GetStatus() == READY);
	if (m < 1 || m > GetModes()) throw;
	return (MBC_P(&mbc))[m - 1];
}

double&
MBCModal::P(uint32_t m)
{
	assert(GetStatus() == READY);
	if (m < 1 || m > GetModes()) throw;
	return (MBC_P(&mbc))[m - 1];
}

#endif // USE_SOCKET
