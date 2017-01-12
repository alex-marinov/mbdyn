/* $Header$ */
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

#ifndef MBC_GLOBAL_PY_H
#define MBC_GLOBAL_PY_H

#ifdef SWIG
%import "mbc_py_global.i"
#endif // SWIG

#ifndef extern_t
#define extern_t extern
#endif // extern_t

/* reference node global data */

extern_t unsigned *mbc_r_k_label;
extern_t double *mbc_r_x;
extern_t double *mbc_r_theta;
extern_t double *mbc_r_r;
extern_t double *mbc_r_euler_123;
extern_t double *mbc_r_xp;
extern_t double *mbc_r_omega;
extern_t double *mbc_r_xpp;
extern_t double *mbc_r_omegap;
extern_t unsigned *mbc_r_d_label;
extern_t double *mbc_r_f;
extern_t double *mbc_r_m;

extern_t unsigned mbc_r_k_label_size;
extern_t unsigned mbc_r_x_size;
extern_t unsigned mbc_r_theta_size;
extern_t unsigned mbc_r_r_size;
extern_t unsigned mbc_r_euler_123_size;
extern_t unsigned mbc_r_xp_size;
extern_t unsigned mbc_r_omega_size;
extern_t unsigned mbc_r_xpp_size;
extern_t unsigned mbc_r_omegap_size;
extern_t unsigned mbc_r_d_label_size;
extern_t unsigned mbc_r_f_size;
extern_t unsigned mbc_r_m_size;

/* nodal element global data */

extern_t unsigned *mbc_n_k_labels;
extern_t double *mbc_n_x;
extern_t double *mbc_n_theta;
extern_t double *mbc_n_r;
extern_t double *mbc_n_euler_123;
extern_t double *mbc_n_xp;
extern_t double *mbc_n_omega;
extern_t double *mbc_n_xpp;
extern_t double *mbc_n_omegap;
extern_t unsigned *mbc_n_d_labels;
extern_t double *mbc_n_f;
extern_t double *mbc_n_m;

extern_t unsigned mbc_n_k_labels_size;
extern_t unsigned mbc_n_x_size;
extern_t unsigned mbc_n_theta_size;
extern_t unsigned mbc_n_r_size;
extern_t unsigned mbc_n_euler_123_size;
extern_t unsigned mbc_n_xp_size;
extern_t unsigned mbc_n_omega_size;
extern_t unsigned mbc_n_xpp_size;
extern_t unsigned mbc_n_omegap_size;
extern_t unsigned mbc_n_d_labels_size;
extern_t unsigned mbc_n_f_size;
extern_t unsigned mbc_n_m_size;

/* modal element global data */

extern_t double *mbc_m_q;
extern_t double *mbc_m_qp;
extern_t double *mbc_m_p;

extern_t uint32_t mbc_m_q_size;
extern_t uint32_t mbc_m_qp_size;
extern_t uint32_t mbc_m_p_size;

#endif // MBC_PY_H
