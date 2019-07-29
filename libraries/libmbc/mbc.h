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

/*
 * NOTE: this is intentionally configuration-independent.
 */

#ifndef MBC_H
#define MBC_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdint.h>

#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
#else
  #ifndef SOCKET_ERROR
    #define SOCKET_ERROR -1
  #endif /* SOCKET_ERROR */
  /*extern int WSAGetLastError(void);*/
  /* We define the SOCKET type if not defined already */
  #ifndef SOCK_TYPEDEF
    #define SOCK_TYPEDEF
    typedef int SOCKET;
  #endif /* SOCK_TYPEDEF */
#endif /* _WIN32 */


/** \brief Legal commands
 *
 * The values of this enumeration appear as tags at the beginning
 * of each communication.
 */
enum ESCmd {
	ES_UNKNOWN				= -1,
	/* 1 intentionally unused */
	ES_REGULAR_DATA				= 2,
	/* 3 intentionally unused */
	ES_GOTO_NEXT_STEP			= 4,
	ES_ABORT				= 5,
	ES_REGULAR_DATA_AND_GOTO_NEXT_STEP	= 6,
	ES_NEGOTIATION				= 7,
	ES_OK					= 8,

	ES_LAST
};

/** \brief Parameters used to control the communication type and fields. */
enum MBCType {
	MBC_MODAL				= 0x0001U,
	MBC_NODAL				= 0x0002U,
	MBC_MODAL_NODAL_MASK			= (MBC_MODAL | MBC_NODAL),

	MBC_REF_NODE				= 0x0004U,

	MBC_ACCELS				= 0x0008U,

	MBC_LABELS				= 0x0010U,

	/** Regular nodes orientation: orientation vector */
	MBC_ROT_THETA				= 0x0100U,
	/** Regular nodes orientation: orientation matrix */
	MBC_ROT_MAT				= 0x0200U,
	/** Regular nodes orientation: Euler angles (123 sequence) */
	MBC_ROT_EULER_123			= 0x0400U,

	/* add more? */

	/** Regular nodes orientation: suppress orientation
	 * (also suppresses angular velocities and moments) */
	MBC_ROT_NONE				= 0x0000U,

	MBC_ROT_MASK				= (MBC_ROT_THETA | MBC_ROT_MAT | MBC_ROT_EULER_123),

	/** Reference node orientation: orientation vector */
	MBC_REF_NODE_ROT_THETA			= (MBC_ROT_THETA << 4),
	/** Reference node orientation: orientation matrix */
	MBC_REF_NODE_ROT_MAT			= (MBC_ROT_MAT << 4),
	/** Reference node orientation: Euler angles (123 sequence) */
	MBC_REF_NODE_ROT_EULER_123		= (MBC_ROT_EULER_123 << 4),

	MBC_REF_NODE_ROT_MASK			= (MBC_ROT_MASK << 4),

	MBC_LAST
};

/** \brief Connection data structure (partially opaque) */
typedef struct {
	/** Opaque. */
	SOCKET	    sock;

	/** Opaque. */
	unsigned	sock_flags;

	/** Opaque. */
	int		recv_flags;

	/** Opaque. */
	int		send_flags;

	/** Opaque. */
	uint8_t		cmd;

	/** Flag that modifies response behavior:
	 * - when 0, on convergence send(2) passes ES_GOTO_NEXT_STEP (convergence);
	 * - when 1, on convergence send(2) passes ES_REGULAR_DATA_AND_GOTO_NEXT_STEP (convergence along with last set of new data).
	 *
	 * In most cases, should be set to 1.
	 */
	char		data_and_next;

	/** Verbose on stdout */
	int		verbose;

	/** Connect(2) timeout.  When peer is not listening:
	 * - a value of 0 results in immediate return with error;
	 * - a positive value corresponds to a timeout in seconds;
	 * - a negative value corresponds to never timing out.
	 *
	 * During timeout, connect(2) is periodically retried.
	 */
	int		timeout;
} mbc_t;

/** \brief Opaque. */
/* validate command
 *
 * command needs to be set in mbc->cmd
 */
extern int
mbc_check_cmd(mbc_t *mbc);

/** \brief Opaque. */
/* get command from peer
 *
 * command is stored in mbc->cmd
 */
extern int
mbc_get_cmd(mbc_t *mbc);

/** \brief Opaque. */
/* put command to peer
 *
 * command needs to be set in mbc->cmd
 */
extern int
mbc_put_cmd(mbc_t *mbc);

/** \brief Initialize communication using "inet" socket.
 *
 * \param [in,out] mbc a pointer to a valid mbc_t structure
 * \param [in] host hostname
 * \param [in] port port number
 *
 * Connects to peer "inet" socket using hostname and port number;
 * host and port must be defined.  If peer is not listening,
 * the behavior depends on mbc_t::timeout.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_inet_init(mbc_t *mbc, const char *host, short unsigned port);

/** \brief Initialize communication using "unix" socket.
 *
 * \param [in,out] mbc a pointer to a valid mbc_t structure
 * \param [in] path pathname
 *
 * Connects to peer "unix" socket using pathname;
 * path must be defined.  If peer is not listening,
 * the behavior depends on mbc::timeout.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_unix_init(mbc_t *mbc, const char *path);

/**
 * \brief Reference node (AKA "rigid") stuff (partially opaque).
 *
 * Users do not need to access members directly;
 * macros documented in the following should be used instead.
 *
 * Flags:
 * - MBC_F_REF_NODE(mbc): true when reference node is defined
 * - MBC_F_LABELS(mbc): true when labels communication is enabled
 * - MBC_F_ACCELS(mbc): true when accelerations communication is enabled
 * - MBC_F_ROT(mbc): orientation mode (see MBC_ROT_* in ::MBCType)
 * - MBC_F_ROT_THETA(mbc): true when orientation mode is orientation vector
 * - MBC_F_ROT_MAT(mbc): true when orientation mode is orientation matrix
 * - MBC_F_ROT_EULER_123(mbc): true when orientation mode is Euler angles (123 sequence)
 * - MBC_F_ROT_REF_NODE(mbc): orientation mode of reference node (see MBC_ROT_* in ::MBCType; MBC_F_REF_NODE(mbc) must be true)
 *
 * Fields (MBC_F_REF_NODE(mbc) must be true):
 * - MBC_R_K_LABEL(mbc): returns the (uint32_t) label of the reference node (MBC_F_LABELS(mbc) must be true)
 * - MBC_R_X(mbc): returns a (double *) pointer to the three components of reference node position
 * - MBC_R_THETA(mbc): returns a (double *) pointer to the three components of the reference node orientation vector
 *   (MBC_F_ROT_THETA(mbc) must be true)
 * - MBC_R_R(mbc): returns a (double *) pointer to the nine components of the reference node orientation matrix,
 *   column-major (Fortran style; MBC_F_ROT_MAT(mbc) must be true)
 * - MBC_R_EULER_123(mbc): returns a (double *) pointer to the three components of the reference node Euler angles
 *   (123 sequence; MBC_F_ROT_EULER_123(mbc) must be true)
 * - MBC_R_XP(mbc): returns a (double *) pointer to the three components of reference node velocity
 * - MBC_R_OMEGA(mbc): returns a (double *) pointer to the three components of reference node angular velocity
 * - MBC_R_XPP(mbc): returns a (double *) pointer to the three components of reference node acceleration
 *   (MBC_F_ACCELS(mbc) must be true)
 * - MBC_R_OMEGAP(mbc): returns a (double *) pointer to the three components of reference node angular acceleration
 *   (MBC_F_ACCELS(mbc) must be true)
 * - MBC_R_D_LABEL(mbc): the (uint32_t) label of the reference node (in the output buffer; MBC_F_LABELS(mbc) must be true)
 * - MBC_R_F(mbc): returns a (double *) pointer to the three components of reference node force
 * - MBC_R_M(mbc): returns a (double *) pointer to the three components of reference node moment
 * - MBC_R_KINEMATICS_SIZE(mbc): returns the size of the input buffer
 * - MBC_R_DYNAMICS_SIZE(mbc): returns the size of the output buffer
 * - MBC_R_SIZE(mbc): returns the cumulative size of the buffers
 * - MBC_R_KINEMATICS(mbc): returns a (void *) pointer to the input buffer
 * - MBC_R_DYNAMICS(mbc): returns a (void *) pointer to the output buffer
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

/**
 * \brief Nodal stuff (partially opaque).
 *
 * Users do not need to access members directly;
 * macros documented in the following should be used instead.
 *
 * Fields:
 * - MBC_N_K_LABELS(mbc): returns a (uint32_t *) pointer to input node labels (MBC_F_LABELS(mbc) must be true)
 * - MBC_N_X(mbc): returns a (double *) pointer to positions of nodes (first node x,y,z; second node x,y,z; ...)
 * - MBC_N_THETA(mbc): returns a (double *) pointer to orientation vectors of nodes
 *   (first node x,y,z; second node x,y,z; ...; MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_R(mbc): returns a (double *) pointer to orientation vectors of nodes
 *   (first node r11,r21,r31,r12,r22,r32,r13,r23,r33; second node r11,...,r33; ...;
 *   MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_EULER_123(mbc): returns a (double *) pointer to Euler angles of nodes
 *   (first node x,y,z; second node x,y,z; ...; MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_XP(mbc): returns a (double *) pointer to velocities of nodes (first node x,y,z; second node x,y,z; ...)
 * - MBC_N_OMEGA(mbc): returns a (double *) pointer to angular velocities of nodes
 *   (first node x,y,z; second node x,y,z; ...; MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_XPP(mbc): returns a (double *) pointer to accelerations of nodes
 *   (first node x,y,z; second node x,y,z; ...; MBC_F_ACCELS(mbc) must be true)
 * - MBC_N_OMEGAP(mbc): returns a (double *) pointer to angular accelerations of nodes
 *   (first node x,y,z; second node x,y,z; ...; MBC_F_ACCELS(mbc) must be true and MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_D_LABELS(mbc): returns a (uint32_t *) pointer to output node labels (MBC_F_LABELS(mbc) must be true)
 * - MBC_N_F(mbc): returns a (double *) pointer to forces (first node x,y,z; second node x,y,z; ...)
 * - MBC_N_M(mbc): returns a (double *) pointer to moments (first node x,y,z; second node x,y,z; ...;
 *   MBC_F_ROT(mbc) must not be ::MBC_ROT_NONE)
 * - MBC_N_KINEMATICS_SIZE(mbc): returns the size of the input buffer
 * - MBC_N_DYNAMICS_SIZE(mbc): returns the size of the output buffer
 * - MBC_N_SIZE(mbc): returns the cumulative size of the buffers
 * - MBC_N_KINEMATICS(mbc): returns a (void *) pointer to the input buffer
 * - MBC_N_DYNAMICS(mbc): returns a (void *) pointer to the output buffer
 * - MBC_N_TIME(mbc): returns a (double *) pointer to the current simulation time
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

	double		*n_t;

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
#define MBC_N_TIME(mbc)			((mbc)->n_t)
} mbc_nodal_t;

/** \brief Initialize nodal data.
 *
 * \param [in,out] mbc pointer to a valid mbc_nodal_t structure
 * \param [in] refnode non-zero if reference node is defined
 * \param [in] nodes number of nodes
 * \param [in] labels true to enable labels
 * \param [in] rot orientation type
 * \param [in] accels true to enable accelerations
 *
 * Either reference node motion must be defined (refnode != 0),
 * or nodes must be > 0, or both.
 *
 * if nodes > 0, this function calls malloc(3) to alloc memory
 * that needs to be freed by calling mbc_nodal_destroy()
 *
 * if labels != 0, labels are handled as well
 *
 * rot must be one of ::MBCType MBC_ROT_*;
 * if it is set to ::MBC_ROT_NONE, only positions and forces are handled.
 *
 * if accels != 0, accelerations are handled as well.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_nodal_init(mbc_nodal_t *mbc, unsigned refnode, unsigned nodes,
	unsigned labels, unsigned rot, unsigned accels);

/** \brief Destroy nodal data.
 *
 * \param [in,out] mbc pointer to a valid mbc_nodal_t structure
 *
 * NOTE: does NOT free mbc.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_nodal_destroy(mbc_nodal_t *mbc);

/** \brief Negotiate nodal data.
 *
 * \param [in] mbc pointer to a valid mbc_nodal_t structure
 *
 * At least reference node motion must be defined
 * (MBC_F_REF_NODE(mbc) must be true),
 * or nodes > 0
 *
 * The socket must be initialized and connected.
 *
 * This function sends a negotiation request to the master.
 *
 * @return 0 on success, !0 on failure.
 */
/*
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

/** \brief Unused. */
/*
 * companion of above, provided for completeness; not used
 */
extern int
mbc_nodal_negotiate_response(mbc_nodal_t *mbc);

/** \brief Get nodal motion from peer.
 *
 * \param [in,out] pointer to a valid mbc_nodal_t structure
 *
 * After the call to this function succeeds:
 * - if MBC_F_REF_NODE(mbc) is true, reference node motion can be accessed using macros MBC_R_X(mbc) and similar;
 * - if mbc_nodal_t::nodes > 0, nodal motion can be accessed using macros MBC_N_X(mbc) and similar.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_nodal_get_motion(mbc_nodal_t *mbc);

/** \brief Put forces to peer.
 *
 * \param [in,out] pointer to a valid mbc_nodal_t structure
 * \param [in] last true when at convergence
 *
* if last is false, before calling this function:
 * - if MBC_F_REF_NODE(mbc) is true, force and moment must be set in storage pointed to by macros MBC_R_F(mbc), MBC_R_M(mbc)
 * - if mbc_nodal_t::nodes > 0, nodal forces must be set in storage pointed to by macros MBC_N_F(mbc), MBC_N_M(mbc)
 *   (the latter only if MBC_F_ROT(mbc) != ::MBC_ROT_NONE).
 *
 * if last is true and mbc_t::data_and_next is false, the output buffer is not sent;
 * thus, there is no need to set forces and moments;
 * otherwise, if mbc_t::data_and_next is true, the output buffer must be filled as described above.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_nodal_put_forces(mbc_nodal_t *mbc, int last);



/**
 * \brief modal stuff (partially opaque).
 *
 * Users do not need to access members directly;
 * macros documented in the following should be used instead.
 *
 * Fields:
 *
 * - MBC_Q(mbc): returns a (double *) pointer to the array of generalized coordinates
 * - MBC_QP(mbc): returns a (double *) pointer to the array of generalized coordinates derivatives
 * - MBC_P(mbc): returns a (double *) pointer to the array of generalized forces
 * - MBC_M_KINEMATICS_SIZE(mbc): returns the size of the input buffer
 * - MBC_M_DYNAMICS_SIZE(mbc): returns the size of the output buffer
 * - MBC_M_SIZE(mbc): returns the cumulative size of the buffers
 * - MBC_M_KINEMATICS(mbc): returns a (void *) pointer to the input buffer
 * - MBC_M_DYNAMICS(mbc): returns a (void *) pointer to the output buffer
 * - MBC_M_TIME(mbc): returns a (double *) pointer to the current simulation time
 */
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
#define MBC_M_KINEMATICS_SIZE(mbc)	((2*(mbc)->modes+1)*sizeof(double))
#define MBC_M_DYNAMICS_SIZE(mbc)	((mbc)->modes*sizeof(double))
#define MBC_M_SIZE(mbc)			(3*(mbc)->modes*sizeof(double))
#define MBC_M_TIME(mbc)			(&(mbc)->m[2*(mbc)->modes+1])
} mbc_modal_t;

/** \brief Initialize modal data.
 *
 * \param [in,out] mbc pointer to a valid mbc_modal_t structure
 * \param [in] refnode non-zero if reference node is defined
 * \param [in] modes number of modes
 *
 * Either reference node motion must be defined (refnode != 0),
 * or modes must be > 0, or both.
 *
 * if modes > 0, this function calls malloc(3) to alloc memory
 * that needs to be freed by calling mbc_modal_destroy()
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_modal_init(mbc_modal_t *mbc, int refnode, unsigned modes);

/** \brief Destroy modal data.
 *
 * \param [in,out] mbc pointer to a valid mbc_modal_t structure
 *
 * NOTE: does NOT free mbc.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_modal_destroy(mbc_modal_t *mbc);

/** \brief Negotiate modal data.
 *
 * \param [in] mbc pointer to a valid mbc_modal_t structure
 *
 * At least reference node motion must be defined
 * (MBC_F_REF_NODE(mbc) must be true),
 * or modes > 0
 *
 * The socket must be initialized and connected.
 *
 * This function sends a negotiation request to the master.
 *
 * @return 0 on success, !0 on failure.
 */
/*
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

/** \brief Unused.
 *
 * companion of above, provided for completeness; not used
 */
extern int
mbc_modal_negotiate_response(mbc_modal_t *mbc);

/** \brief Get modal motion from peer.
 *
 * \param [in,out] pointer to a valid mbc_modal_t structure
 *
 * After the call to this function succeeds:
 * - if MBC_F_REF_NODE(mbc) is true, reference node motion can be accessed using macros MBC_R_X(mbc) and similar;
 * - if mbc_modal_t::modes > 0, modal motion can be accessed using macros MBC_Q(mbc) and similar.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_modal_get_motion(mbc_modal_t *mbc);

/** \brief Put forces to peer.
 *
 * \param [in,out] mbc pointer to a valid mbc_modal_t structure
 * \param [in] last true when at convergence
 *
 * if last is false, before calling this function:
 * - if MBC_F_REF_NODE(mbc) is true, force and moment must be set in storage pointed to by macros MBC_R_F(mbc), MBC_R_M(mbc)
 * - if mbc_modal_t::modes > 0, nodal forces must be set in storage pointed to by macro MBC_P(mbc)
 *
 * if last is true and mbc_t::data_and_next is false, the output buffer is not sent;
 * thus, there is no need to set generalized forces;
 * otherwise, if mbc_t::data_and_next is true, the output buffer must be filled as described above.
 *
 * @return 0 on success, !0 on failure.
 */
extern int
mbc_modal_put_forces(mbc_modal_t *mbc, int last);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MBC_H */
