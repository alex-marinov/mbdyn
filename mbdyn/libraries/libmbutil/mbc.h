/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* legal commands */
enum ESCmd {
	ES_UNKNOWN				= -1,
	ES_REGULAR_DATA				= 2,
	ES_GOTO_NEXT_STEP			= 4,
	ES_ABORT				= 5,
	ES_REGULAR_DATA_AND_GOTO_NEXT_STEP	= 6,

	ES_LAST
};

/* command structure */
typedef struct {
	int		sock;
	int		recv_flags;
	int		send_flags;
	unsigned	cmd;
	char		data_and_next;
	int		verbose;
} mbc_t;

/* modal data structure */
typedef struct {
	mbc_t		mbc;

	/* rigid body data */
	char		rigid;
	double		r[3 + 9 + 3 + 3 + 3 + 3];
#define	MBC_X(mbc)			(&(mbc)->r[0])
#define	MBC_R(mbc)			(&(mbc)->r[3])
#define	MBC_V(mbc)			(&(mbc)->r[12])
#define	MBC_W(mbc)			(&(mbc)->r[15])
#define	MBC_F(mbc)			(&(mbc)->r[18])
#define	MBC_M(mbc)			(&(mbc)->r[21])
#define MBC_R_KINEMATICS(mbc)		MBC_X((mbc))
#define MBC_R_DYNAMICS(mbc)		MBC_F((mbc))
#define MBC_R_KINEMATICS_SIZE(mbc)	(18*sizeof(double))
#define MBC_R_DYNAMICS_SIZE(mbc)	(6*sizeof(double))
#define MBC_R_SIZE(mbc)			((18 + 6)*sizeof(double))

	/* modal data */
	int		modes;
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

/* get modal motion from peer
 *
 * if mbc->rigid, access rigid motion using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->modes > 0, access modal motion using macros MBC_Q, MBC_QP
 */
extern int
mbc_modal_get_motion(mbc_modal_t *mbc);

/* put command to peer
 *
 * command needs to be set in mbc->cmd
 */
extern int
mbc_put_cmd(mbc_t *mbc);

/* put forces to peer
 *
 * if mbc->rigid, force and moment must be set in storage pointed to
 *	by macros MBC_F, MBC_M
 * if mbc->modes > 0, modal forces must be set in storage pointed to
 *	by macro MBC_P
 */
extern int
mbc_modal_put_forces(mbc_modal_t *mbc, int last);

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

/* destroy communication
 *
 * does NOT free the mbc structure
 */
extern int
mbc_destroy(mbc_t *mbc);

/* initialize modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least rigid body motion must be defined (mbc->rigid != 0),
 * or modes must be > 0
 *
 * if modes > 0, mallocs memory that needs to be freed calling
 * mbc_modal_destroy()
 */
extern int
mbc_modal_init(mbc_modal_t *mbc, int modes);

/* destroy modal data
 *
 * does NOT free the mbc structure
 */
extern int
mbc_modal_destroy(mbc_modal_t *mbc);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MBC_H */
