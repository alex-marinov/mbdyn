/* $Header$ */
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

#include "mbc_py.h"

#include <iostream>

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

static mbc_nodal_t priv_mbc, *mbc = &priv_mbc;

int
mbc_py_nodal_initialize(const char *const path,
	const char *const host, unsigned port,
	int timeout, unsigned verbose, unsigned data_and_next,
	unsigned rigid, unsigned nodes,
	unsigned labels, unsigned rot, unsigned accels)
{
	::mbc->mbc.data_and_next = data_and_next;
	::mbc->mbc.verbose = verbose;
	::mbc->mbc.timeout = timeout;

	if (::mbc->mbc.verbose) {
		if (path && path[0]) {
			std::cout << "connecting to path=" << path << std::endl;
		} else {
			std::cout << "connecting to host=" << host << ":" << port << std::endl;
		}

		if (::mbc->mbc.timeout < 0) {
			std::cout << "timeout=forever" << std::endl;
		} else {
			std::cout << "timeout=" << ::mbc->mbc.timeout << std::endl;
		}
	}

	if (path && path[0]) {
		if (mbc_unix_init((mbc_t *)::mbc, path)) {
			return -1;
		}

	} else if (host && host[0]) {
		if (mbc_inet_init((mbc_t *)::mbc, host, port)) {
			return -1;
		}

	} else {
		return -1;
	}

	if (mbc_nodal_init(::mbc, rigid, nodes, labels, rot, accels)) {
		return -1;
	}

	if (mbc_nodal_negotiate_request(::mbc)) {
		return -1;
	}

	if (rigid) {
		if (labels) {
			mbc_r_k_label = &MBC_R_K_LABEL(::mbc);
			mbc_r_k_label_size = 1;

			mbc_r_d_label = &MBC_R_D_LABEL(::mbc);
			mbc_r_d_label_size = 1;
		}

		mbc_r_x = MBC_R_X(::mbc);
		mbc_r_x_size = 3;

		switch (rot) {
		case MBC_ROT_THETA:
			mbc_r_theta = MBC_R_THETA(::mbc);
			mbc_r_theta_size = 3;
			break;

		case MBC_ROT_MAT:
			mbc_r_r = MBC_R_R(::mbc);
			mbc_r_r_size = 9;
			break;

		case MBC_ROT_EULER_123:
			mbc_r_euler_123 = MBC_R_EULER_123(::mbc);
			mbc_r_euler_123_size = 3;
			break;
		}

		mbc_r_xp = MBC_R_XP(::mbc);
		mbc_r_xp_size = 3;
		mbc_r_omega = MBC_R_OMEGA(::mbc);
		mbc_r_omega_size = 3;

		if (accels) {
			mbc_r_xpp = MBC_R_XPP(::mbc);
			mbc_r_xpp_size = 3;
			mbc_r_omegap = MBC_R_OMEGAP(::mbc);
			mbc_r_omegap_size = 3;
		}

		mbc_r_f = MBC_R_F(::mbc);
		mbc_r_f_size = 3;
		mbc_r_m = MBC_R_M(::mbc);
		mbc_r_m_size = 3;
	}

	if (nodes > 0) {
		if (labels) {
			mbc_n_k_labels = MBC_N_K_LABELS(::mbc);
			mbc_n_k_labels_size = nodes;
			mbc_n_d_labels = MBC_N_D_LABELS(::mbc);
			mbc_n_d_labels_size = nodes;
		}

		mbc_n_x = MBC_N_X(::mbc);
		mbc_n_x_size = 3*nodes;

		switch (rot) {
		case MBC_ROT_THETA:
			mbc_n_theta = MBC_N_THETA(::mbc);
			mbc_n_theta_size = 3*nodes;
			break;

		case MBC_ROT_MAT:
			mbc_n_r = MBC_N_R(::mbc);
			mbc_n_r_size = 9*nodes;
			break;

		case MBC_ROT_EULER_123:
			mbc_n_euler_123 = MBC_N_EULER_123(::mbc);
			mbc_n_euler_123_size = 3*nodes;
		}

		mbc_n_xp = MBC_N_XP(::mbc);
		mbc_n_xp_size = 3*nodes;
		mbc_n_omega = MBC_N_OMEGA(::mbc);
		mbc_n_omega_size = 3*nodes;

		if (accels) {
			mbc_n_xpp = MBC_N_XPP(::mbc);
			mbc_n_xpp_size = 3*nodes;
			mbc_n_omegap = MBC_N_OMEGAP(::mbc);
			mbc_n_omegap_size = 3*nodes;
		}

		mbc_n_f = MBC_N_F(::mbc);
		mbc_n_f_size = 3*nodes;
		mbc_n_m = MBC_N_M(::mbc);
		mbc_n_m_size = 3*nodes;
	}

	return 0;
}

int
mbc_py_nodal_send(int last)
{
	return mbc_nodal_put_forces(::mbc, last);
}

int
mbc_py_nodal_recv(void)
{
	return mbc_nodal_get_motion(::mbc);
}

void
mbc_py_nodal_destroy(void)
{
	mbc_nodal_destroy(::mbc);
}
