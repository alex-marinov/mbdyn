/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbc_py.h"

#include <iostream>
#include <cstring>
#include <vector>

/* reference node global data */

uint32_t *mbc_r_k_label;
double *mbc_r_x;
double *mbc_r_theta;
double *mbc_r_r;
double *mbc_r_euler_123;
double *mbc_r_xp;
double *mbc_r_omega;
double *mbc_r_xpp;
double *mbc_r_omegap;
uint32_t *mbc_r_d_label;
double *mbc_r_f;
double *mbc_r_m;

uint32_t mbc_r_k_label_size;
uint32_t mbc_r_x_size;
uint32_t mbc_r_theta_size;
uint32_t mbc_r_r_size;
uint32_t mbc_r_euler_123_size;
uint32_t mbc_r_xp_size;
uint32_t mbc_r_omega_size;
uint32_t mbc_r_xpp_size;
uint32_t mbc_r_omegap_size;
uint32_t mbc_r_d_label_size;
uint32_t mbc_r_f_size;
uint32_t mbc_r_m_size;

/* nodal element global data */

uint32_t *mbc_n_k_labels;
double *mbc_n_x;
double *mbc_n_theta;
double *mbc_n_r;
double *mbc_n_euler_123;
double *mbc_n_xp;
double *mbc_n_omega;
double *mbc_n_xpp;
double *mbc_n_omegap;
uint32_t *mbc_n_d_labels;
double *mbc_n_f;
double *mbc_n_m;

uint32_t mbc_n_k_labels_size;
uint32_t mbc_n_x_size;
uint32_t mbc_n_theta_size;
uint32_t mbc_n_r_size;
uint32_t mbc_n_euler_123_size;
uint32_t mbc_n_xp_size;
uint32_t mbc_n_omega_size;
uint32_t mbc_n_xpp_size;
uint32_t mbc_n_omegap_size;
uint32_t mbc_n_d_labels_size;
uint32_t mbc_n_f_size;
uint32_t mbc_n_m_size;

static std::vector<mbc_nodal_t> n_mbc;

int
mbc_py_nodal_initialize(const char *const path,
	const char *const host, unsigned port,
	int timeout, unsigned verbose, unsigned data_and_next,
	unsigned refnode, unsigned nodes,
	unsigned labels, unsigned rot, unsigned accels)
{
	mbc_nodal_t mbc;
	std::memset(&mbc, 0, sizeof(mbc_nodal_t));

	mbc.mbc.data_and_next = data_and_next;
	mbc.mbc.verbose = verbose;
	mbc.mbc.timeout = timeout;

	if (mbc.mbc.verbose) {
		if (path && path[0]) {
			std::cout << "connecting to path=" << path << std::endl;
		} else {
			std::cout << "connecting to host=" << host << ":" << port << std::endl;
		}

		if (mbc.mbc.timeout < 0) {
			std::cout << "timeout=forever" << std::endl;
		} else {
			std::cout << "timeout=" << mbc.mbc.timeout << std::endl;
		}
	}

	if (path && path[0]) {
		if (mbc_unix_init((mbc_t *)&mbc, path)) {
			return -1;
		}

	} else if (host && host[0]) {
		if (mbc_inet_init((mbc_t *)&mbc, host, port)) {
			return -1;
		}

	} else {
		return -1;
	}

	if (mbc_nodal_init(&mbc, refnode, nodes, labels, rot, accels)) {
		return -1;
	}

	int id = ::n_mbc.size();
	::n_mbc.push_back(mbc);

	return id;
}

int
mbc_py_nodal_negotiate(unsigned id)
{
	if (id >= n_mbc.size()) {
		return -1;
	}

	mbc_nodal_t *mbcp = &::n_mbc[id];

	if (mbc_nodal_negotiate_request(mbcp)) {
		return -1;
	}

	if (MBC_F_REF_NODE(mbcp)) {
		if (MBC_F_LABELS(mbcp)) {
			mbc_r_k_label = &MBC_R_K_LABEL(mbcp);
			mbc_r_k_label_size = 1;

			mbc_r_d_label = &MBC_R_D_LABEL(mbcp);
			mbc_r_d_label_size = 1;
		}

		mbc_r_x = MBC_R_X(mbcp);
		mbc_r_x_size = 3;

		switch (MBC_F_ROT_REF_NODE(mbcp)) {
		case MBC_ROT_THETA:
			mbc_r_theta = MBC_R_THETA(mbcp);
			mbc_r_theta_size = 3;
			break;

		case MBC_ROT_MAT:
			mbc_r_r = MBC_R_R(mbcp);
			mbc_r_r_size = 9;
			break;

		case MBC_ROT_EULER_123:
			mbc_r_euler_123 = MBC_R_EULER_123(mbcp);
			mbc_r_euler_123_size = 3;
			break;

		default:
			std::cerr << "mbc_py_nodal_negotiate: unknown reference node orientation" << std::endl;
			return -1;
		}

		mbc_r_xp = MBC_R_XP(mbcp);
		mbc_r_xp_size = 3;
		mbc_r_omega = MBC_R_OMEGA(mbcp);
		mbc_r_omega_size = 3;

		if (MBC_F_ACCELS(mbcp)) {
			mbc_r_xpp = MBC_R_XPP(mbcp);
			mbc_r_xpp_size = 3;
			mbc_r_omegap = MBC_R_OMEGAP(mbcp);
			mbc_r_omegap_size = 3;
		}

		mbc_r_f = MBC_R_F(mbcp);
		mbc_r_f_size = 3;
		mbc_r_m = MBC_R_M(mbcp);
		mbc_r_m_size = 3;
	}

	if (mbcp->nodes > 0) {
		if (MBC_F_LABELS(mbcp)) {
			mbc_n_k_labels = MBC_N_K_LABELS(mbcp);
			mbc_n_k_labels_size = mbcp->nodes;
			mbc_n_d_labels = MBC_N_D_LABELS(mbcp);
			mbc_n_d_labels_size = mbcp->nodes;
		}

		mbc_n_x = MBC_N_X(mbcp);
		mbc_n_x_size = 3*mbcp->nodes;

		switch (MBC_F_ROT(mbcp)) {
		case MBC_ROT_THETA:
			mbc_n_theta = MBC_N_THETA(mbcp);
			mbc_n_theta_size = 3*mbcp->nodes;
			break;

		case MBC_ROT_MAT:
			mbc_n_r = MBC_N_R(mbcp);
			mbc_n_r_size = 9*mbcp->nodes;
			break;

		case MBC_ROT_EULER_123:
			mbc_n_euler_123 = MBC_N_EULER_123(mbcp);
			mbc_n_euler_123_size = 3*mbcp->nodes;

		case MBC_ROT_NONE:
			break;

		default:
			std::cerr << "mbc_py_nodal_negotiate: unknown node orientation " << MBC_F_ROT(mbcp) << std::endl;
			return -1;
		}

		mbc_n_xp = MBC_N_XP(mbcp);
		mbc_n_xp_size = 3*mbcp->nodes;
		mbc_n_omega = MBC_N_OMEGA(mbcp);
		mbc_n_omega_size = 3*mbcp->nodes;

		if (MBC_F_ACCELS(mbcp)) {
			mbc_n_xpp = MBC_N_XPP(mbcp);
			mbc_n_xpp_size = 3*mbcp->nodes;
			mbc_n_omegap = MBC_N_OMEGAP(mbcp);
			mbc_n_omegap_size = 3*mbcp->nodes;
		}

		mbc_n_f = MBC_N_F(mbcp);
		mbc_n_f_size = 3*mbcp->nodes;
		mbc_n_m = MBC_N_M(mbcp);
		mbc_n_m_size = 3*mbcp->nodes;
	}

	return id;
}

int
mbc_py_nodal_send(unsigned id, int last)
{
	if (id >= n_mbc.size()) {
		return -1;
	}

	return mbc_nodal_put_forces(&::n_mbc[id], last);
}

int
mbc_py_nodal_recv(unsigned id)
{
	if (id >= n_mbc.size()) {
		return -1;
	}

	return mbc_nodal_get_motion(&::n_mbc[id]);
}

int
mbc_py_nodal_destroy(unsigned id)
{
	if (id >= n_mbc.size()) {
		return -1;
	}

	mbc_nodal_destroy(&::n_mbc[id]);

	return 0;
}

/* modal element global data */

double *mbc_m_q;
double *mbc_m_qp;
double *mbc_m_p;

uint32_t mbc_m_q_size;
uint32_t mbc_m_qp_size;
uint32_t mbc_m_p_size;

static std::vector<mbc_modal_t> m_mbc;

int
mbc_py_modal_initialize(const char *const path,
	const char *const host, unsigned port,
	int timeout, unsigned verbose, unsigned data_and_next,
	unsigned refnode, unsigned modes)
{
	mbc_modal_t mbc;
	std::memset(&mbc, 0, sizeof(mbc_modal_t));

	mbc.mbc.data_and_next = data_and_next;
	mbc.mbc.verbose = verbose;
	mbc.mbc.timeout = timeout;

	if (mbc.mbc.verbose) {
		if (path && path[0]) {
			std::cout << "connecting to path=" << path << std::endl;
		} else {
			std::cout << "connecting to host=" << host << ":" << port << std::endl;
		}

		if (mbc.mbc.timeout < 0) {
			std::cout << "timeout=forever" << std::endl;
		} else {
			std::cout << "timeout=" << mbc.mbc.timeout << std::endl;
		}
	}

	if (path && path[0]) {
		if (mbc_unix_init((mbc_t *)&mbc, path)) {
			return -1;
		}

	} else if (host && host[0]) {
		if (mbc_inet_init((mbc_t *)&mbc, host, port)) {
			return -1;
		}

	} else {
		return -1;
	}

	if (mbc_modal_init(&mbc, refnode, modes)) {
		return -1;
	}

	int id = m_mbc.size();
	m_mbc.push_back(mbc);

	return id;
}

int
mbc_py_modal_negotiate(unsigned id)
{
	if (id >= m_mbc.size()) {
		return -1;
	}

	mbc_modal_t *mbcp = &::m_mbc[id];

	if (mbc_modal_negotiate_request(mbcp)) {
		return -1;
	}

	if (MBC_F_REF_NODE(mbcp)) {
		mbc_r_x = MBC_R_X(mbcp);
		mbc_r_x_size = 3;
		mbc_r_r = MBC_R_R(mbcp);
		mbc_r_r_size = 9;

		mbc_r_xp = MBC_R_XP(mbcp);
		mbc_r_xp_size = 3;
		mbc_r_omega = MBC_R_OMEGA(mbcp);
		mbc_r_omega_size = 3;

		mbc_r_f = MBC_R_F(mbcp);
		mbc_r_f_size = 3;
		mbc_r_m = MBC_R_M(mbcp);
		mbc_r_m_size = 3;
	}

	if (mbcp->modes > 0) {
		mbc_m_q = MBC_Q(mbcp);
		mbc_m_q_size = mbcp->modes;
		mbc_m_qp = MBC_QP(mbcp);
		mbc_m_qp_size = mbcp->modes;

		mbc_m_p = MBC_P(mbcp);
		mbc_m_p_size = mbcp->modes;
	}

	return id;
}

int
mbc_py_modal_send(unsigned id, int last)
{
	if (id >= m_mbc.size()) {
		return -1;
	}

	return mbc_modal_put_forces(&m_mbc[id], last);
}

int
mbc_py_modal_recv(unsigned id)
{
	if (id >= m_mbc.size()) {
		return -1;
	}

	return mbc_modal_get_motion(&m_mbc[id]);
}

int
mbc_py_modal_destroy(unsigned id)
{
	if (id >= m_mbc.size()) {
		return -1;
	}

	mbc_modal_destroy(&m_mbc[id]);

	return 0;
}
