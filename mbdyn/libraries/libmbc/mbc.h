/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

/*
 * NOTE: this is intentionally configuration-independent.
 */

#ifndef MBC_H
#define MBC_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdint.h>

/* legal commands */
enum ESCmd {
	ES_UNKNOWN				= -1,
	ES_REGULAR_DATA				= 2,
	ES_GOTO_NEXT_STEP			= 4,
	ES_ABORT				= 5,
	ES_REGULAR_DATA_AND_GOTO_NEXT_STEP	= 6,
	ES_NEGOTIATION				= 7,
	ES_OK					= 8,

	ES_LAST
};

enum MBCType {
	MBC_MODAL				= 0x0001U,
	MBC_NODAL				= 0x0002U,
	MBC_MODAL_NODAL_MASK			= (MBC_MODAL | MBC_NODAL),

	MBC_REF_NODE				= 0x0004U,

	MBC_ACCELS				= 0x0008U,

	MBC_LABELS				= 0x0010U,

	MBC_ROT_THETA				= 0x0100U,
	MBC_ROT_MAT				= 0x0200U,
	MBC_ROT_EULER_123			= 0x0400U,
	/* add more? */
	MBC_ROT_NONE				= 0x0000U,
	MBC_ROT_MASK				= (MBC_ROT_THETA | MBC_ROT_MAT | MBC_ROT_EULER_123),

	MBC_REF_NODE_ROT_THETA			= (MBC_ROT_THETA << 4),
	MBC_REF_NODE_ROT_MAT			= (MBC_ROT_MAT << 4),
	MBC_REF_NODE_ROT_EULER_123		= (MBC_ROT_EULER_123 << 4),
	MBC_REF_NODE_ROT_MASK			= (MBC_ROT_MASK << 4),

	MBC_LAST
};

/* command structure */
typedef struct {
	int		sock;
	unsigned	sock_flags;
	int		recv_flags;
	int		send_flags;
	uint8_t		cmd;
	char		data_and_next;
	int		verbose;
	int		timeout;
} mbc_t;

/* validate command
 *
 * command needs to be set in mbc->cmd
 */
extern int
mbc_check_cmd(mbc_t *mbc);

/* get command from peer
 *
 * command is stored in mbc->cmd
 */
extern int
mbc_get_cmd(mbc_t *mbc);

/* put command to peer
 *
 * command needs to be set in mbc->cmd
 */
extern int
mbc_put_cmd(mbc_t *mbc);

/* initialize communication using inet socket
 *
 * mbc must be a pointer to a valid mbc_t structure
 * host and port must be defined
 */
extern int
mbc_inet_init(mbc_t *mbc, const char *host, short unsigned port);

/* initialize communication using unix socket
 *
 * mbc must be a pointer to a valid mbc_t structure
 * path must be defined
 */
extern int
mbc_unix_init(mbc_t *mbc, const char *path);

/*
 * reference node (AKA "rigid") stuff
 */
typedef struct {
	uint32_t	flags;
#define MBC_F(mbc)			((mbc)->mbcr.flags)
#define MBC_F_GET(mbc, f)		(MBC_F(mbc) & (f))
#define MBC_F_SET(mbc, f)		(MBC_F(mbc) |= (f))
#define MBC_F_RESET(mbc, f)		(MBC_F(mbc) &= ~(f))

#define MBC_F_REF_NODE(mbc)		MBC_F_GET(mbc, MBC_REF_NODE)
#define MBC_F_LABELS(mbc)		MBC_F_GET(mbc, MBC_LABELS)
#define MBC_F_ACCELS(mbc)		MBC_F_GET(mbc, MBC_ACCELS)
#define MBC_F_ROT(mbc)			MBC_F_GET(mbc, MBC_ROT_MASK)
#define MBC_F_ROT_THETA(mbc)		MBC_F_GET(mbc, MBC_ROT_THETA)
#define MBC_F_ROT_MAT(mbc)		MBC_F_GET(mbc, MBC_ROT_MAT)
#define MBC_F_ROT_EULER_123(mbc)	MBC_F_GET(mbc, MBC_ROT_EULER_123)
#define MBC_U_REF_NODE_ROT_2_ROT(u)	(((u) & MBC_REF_NODE_ROT_MASK) >> 4)
#define MBC_U_ROT_2_REF_NODE_ROT(u)	(((u) & MBC_ROT_MASK) << 4)
#define MBC_F_ROT_REF_NODE(mbc)		MBC_U_REF_NODE_ROT_2_ROT(MBC_F_GET((mbc), MBC_REF_NODE_ROT_MASK))
#define MBC_F_REF_NODE_ROT(mbc)		MBC_F_GET(mbc, MBC_REF_NODE_ROT_MASK)

#define MBC_F_SET_REF_NODE(mbc)		MBC_F_SET(mbc, MBC_REF_NODE)
#define MBC_F_SET_LABELS(mbc)		MBC_F_SET(mbc, MBC_LABELS)
#define MBC_F_SET_ACCELS(mbc)		MBC_F_SET(mbc, MBC_ACCELS)
#define MBC_F_SET_ROT_THETA(mbc)	MBC_F_SET(mbc, MBC_ROT_THETA)
#define MBC_F_SET_ROT_MAT(mbc)		MBC_F_SET(mbc, MBC_ROT_MAT)
#define MBC_F_SET_ROT_EULER_123(mbc)	MBC_F_SET(mbc, MBC_ROT_EULER_123)

	/* reference node data */
	union {
		/* NOTE: wasting an uint32_t for each uint32_t to make sure alignment is on double */
		char		char_r_ptr[(2 + 2) * sizeof(uint32_t) + (3 + 9 + 3 + 3 + 3 + 3 + 3 + 3) * sizeof(double)];
		uint32_t	uint32_t_r_ptr[(2 + 2) + (3 + 9 + 3 + 3 + 3 + 3 + 3 + 3) * sizeof(double) / sizeof(uint32_t)];
		double		double_r_ptr[(2 + 2) * sizeof(uint32_t) / sizeof(double) + (3 + 9 + 3 + 3 + 3 + 3 + 3 + 3) * sizeof(double)];
	} r_ptr;
	uint32_t	k_size;
	int32_t		r_k_label;
	int32_t		r_k_x;
	int32_t		r_k_theta;
	int32_t		r_k_r;
	int32_t		r_k_euler_123;
	int32_t		r_k_xp;
	int32_t		r_k_omega;
	int32_t		r_k_xpp;
	int32_t		r_k_omegap;
	uint32_t	d_size;
	int32_t		r_d_label;
	int32_t		r_d_f;
	int32_t		r_d_m;
#define MBC_R_PTR(mbc, type, off)	((off) < 0 ? NULL : ((type *)&((mbc)->mbcr.r_ptr.type ## _r_ptr[(off)])))
#define	MBC_R_K_LABEL(mbc)		(MBC_R_PTR((mbc), uint32_t, (mbc)->mbcr.r_k_label)[0])
#define	MBC_R_X(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_x))
#define	MBC_R_THETA(mbc)		(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_theta))
#define	MBC_R_R(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_r))
#define	MBC_R_EULER_123(mbc)		(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_euler_123))
#define	MBC_R_XP(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_xp))
#define	MBC_R_OMEGA(mbc)		(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_omega))
#define	MBC_R_XPP(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_xpp))
#define	MBC_R_OMEGAP(mbc)		(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_k_omegap))
#define	MBC_R_D_LABEL(mbc)		(MBC_R_PTR((mbc), uint32_t, (mbc)->mbcr.r_d_label)[0])
#define	MBC_R_F(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_d_f))
#define	MBC_R_M(mbc)			(MBC_R_PTR((mbc), double, (mbc)->mbcr.r_d_m))
#define MBC_R_KINEMATICS_SIZE(mbc)	((mbc)->mbcr.k_size)
#define MBC_R_DYNAMICS_SIZE(mbc)	((mbc)->mbcr.d_size)
#define MBC_R_SIZE(mbc)			(MBC_R_KINEMATICS_SIZE(mbc) + MBC_R_DYNAMICS_SIZE(mbc))
#define MBC_R_KINEMATICS(mbc)		((void *)(mbc)->mbcr.r_ptr.char_r_ptr)
#define MBC_R_DYNAMICS(mbc)		((void *)(MBC_R_KINEMATICS(mbc) + MBC_R_KINEMATICS_SIZE(mbc)))
} mbc_rigid_t;

/*
 * nodal stuff
 */
typedef struct {
	mbc_t		mbc;
	mbc_rigid_t	mbcr;

	/* nodal data */
	uint32_t	nodes;
	uint32_t	k_size;

	char		*n_ptr;
	/* n_k_labels optional */
	uint32_t	*n_k_labels;
	double		*n_k_x;
	/* n_k_theta, n_k_r, n_k_euler_123 mutually exclusive */
	double		*n_k_theta;
	double		*n_k_r;
	double		*n_k_euler_123;
	double		*n_k_xp;
	double		*n_k_omega;
	/* n_k_xpp, n_k_omegap optional */
	double		*n_k_xpp;
	double		*n_k_omegap;

	uint32_t	d_size;
	/* n_d_labels optional */
	uint32_t	*n_d_labels;
	double		*n_d_f;
	double		*n_d_m;

#define MBC_N_K_LABELS(mbc)		((mbc)->n_k_labels)
#define MBC_N_X(mbc)			((mbc)->n_k_x)
#define MBC_N_THETA(mbc)		((mbc)->n_k_theta)
#define MBC_N_R(mbc)			((mbc)->n_k_r)
#define MBC_N_EULER_123(mbc)		((mbc)->n_k_euler_123)
#define MBC_N_XP(mbc)			((mbc)->n_k_xp)
#define MBC_N_OMEGA(mbc)		((mbc)->n_k_omega)
#define MBC_N_XPP(mbc)			((mbc)->n_k_xpp)
#define MBC_N_OMEGAP(mbc)		((mbc)->n_k_omegap)
#define MBC_N_D_LABELS(mbc)		((mbc)->n_d_labels)
#define MBC_N_F(mbc)			((mbc)->n_d_f)
#define MBC_N_M(mbc)			((mbc)->n_d_m)
#define MBC_N_KINEMATICS_SIZE(mbc)	((mbc)->k_size)
#define MBC_N_DYNAMICS_SIZE(mbc)	((mbc)->d_size)
#define MBC_N_KINEMATICS(mbc)		((void *)(mbc)->n_ptr)
#define MBC_N_DYNAMICS(mbc)		((void *)(MBC_N_KINEMATICS(mbc) + MBC_N_KINEMATICS_SIZE(mbc)))
#define MBC_N_SIZE(mbc)			(MBC_N_KINEMATICS_SIZE(mbc) + MBC_N_DYNAMICS_SIZE(mbc))
} mbc_nodal_t;

/* initialize nodal data
 *
 * mbc must be a pointer to a valid mbc_nodal_t structure
 *
 * at least reference node motion must be defined (refnode != 0),
 * or nodes must be > 0
 *
 * if nodes > 0, mallocs memory that needs to be freed calling
 * mbc_nodal_destroy()
 *
 * rot must be one of MBC_ROT_*
 *
 * if accelerations != 0 accelerations are also output
 */
extern int
mbc_nodal_init(mbc_nodal_t *mbc, unsigned refnode, unsigned nodes,
	unsigned labels, unsigned rot, unsigned accels);

/* destroy nodal data
 *
 * does NOT free the mbc structure
 */
extern int
mbc_nodal_destroy(mbc_nodal_t *mbc);

/* negotiate nodal data
 *
 * mbc must be a pointer to a valid mbc_nodal_t structure
 *
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or nodes must be > 0
 *
 * the socket must be initialized and connected
 *
 * this function sends a negotiation request to the master
 *
 * the protocol consists in:
 *
 * - negotiation request:
 *   - the negotiation request tag (ES_NEGOTIATION, uint8_t)
 *   - a bunch of flags (uint32_t), containing:
 *       - the type of interface (MBC_MODAL or MBC_NODAL)
 *       - whether a reference node is defined, MBC_REF_NODE OR-ed to previous
 *       - more...
 *   - the number of nodes (uint32_t)
 *
 * - negotiation response:
 *   - the negotiation response tag (ES_OK or ES_ABORT, uint8_t)
 */
extern int
mbc_nodal_negotiate_request(mbc_nodal_t *mbc);

/* companion of above, provided for completeness; not used
 */
extern int
mbc_nodal_negotiate_response(mbc_nodal_t *mbc);

/* get nodal motion from peer
 *
 * if MBC_F_REF_NODE(mbc), access reference node motion using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->nodes > 0, access nodal motion using macros MBC_N_*
 */
extern int
mbc_nodal_get_motion(mbc_nodal_t *mbc);

/* put forces to peer
 *
 * if MBC_F_REF_NODE(mbc), force and moment must be set in storage pointed to
 *	by macros MBC_F, MBC_M
 * if mbc->nodes > 0, nodal forces must be set in storage pointed to
 *	by macro MBC_N_F, MBC_N_M
 */
extern int
mbc_nodal_put_forces(mbc_nodal_t *mbc, int last);



/*
 * modal stuff
 */

/* modal data structure */
typedef struct {
	mbc_t		mbc;
	mbc_rigid_t	mbcr;

	/* modal data */
	uint32_t	modes;
	double		*m;
#define MBC_Q(mbc)			(&(mbc)->m[0])
#define MBC_QP(mbc)			(&(mbc)->m[(mbc)->modes])
#define MBC_P(mbc)			(&(mbc)->m[2*(mbc)->modes])
#define MBC_M_KINEMATICS(mbc)		MBC_Q((mbc))
#define MBC_M_DYNAMICS(mbc)		MBC_P((mbc))
#define MBC_M_KINEMATICS_SIZE(mbc)	(2*(mbc)->modes*sizeof(double))
#define MBC_M_DYNAMICS_SIZE(mbc)	((mbc)->modes*sizeof(double))
#define MBC_M_SIZE(mbc)			(3*(mbc)->modes*sizeof(double))
} mbc_modal_t;

/* initialize modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or modes must be > 0
 *
 * if modes > 0, mallocs memory that needs to be freed calling
 * mbc_modal_destroy()
 */
extern int
mbc_modal_init(mbc_modal_t *mbc, int refnode, unsigned modes);

/* destroy modal data
 *
 * does NOT free the mbc structure
 */
extern int
mbc_modal_destroy(mbc_modal_t *mbc);

/* negotiate modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or modes must be > 0
 *
 * the socket must be initialized and connected
 * sends a negotiation request to the master
 *
 * the protocol consists in:
 *
 * - negotiation request:
 *   - the negotiation request tag (ES_NEGOTIATION, uint8_t)
 *   - a bunch of flags (uint32_t), containing:
 *       - the type of interface (MBC_MODAL or MBC_NODAL)
 *       - whether a reference node is defined, MBC_REF_NODE OR-ed to previous
 *       - more...
 *   - the number of modes (uint32_t)
 *
 * - negotiation response:
 *   - the negotiation response tag (ES_OK or ES_ABORT, uint8_t)
 */
extern int
mbc_modal_negotiate_request(mbc_modal_t *mbc);

/* companion of above, provided for completeness; not used
 */
extern int
mbc_modal_negotiate_response(mbc_modal_t *mbc);

/* get modal motion from peer
 *
 * if MBC_F_REF_NODE(mbc), access reference node motion using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->modes > 0, access modal motion using macros MBC_Q, MBC_QP
 */
extern int
mbc_modal_get_motion(mbc_modal_t *mbc);

/* put forces to peer
 *
 * if MBC_F_REF_NODE(mbc), force and moment must be set in storage pointed to
 *	by macros MBC_F, MBC_M
 * if mbc->modes > 0, modal forces must be set in storage pointed to
 *	by macro MBC_P
 */
extern int
mbc_modal_put_forces(mbc_modal_t *mbc, int last);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MBC_H */
