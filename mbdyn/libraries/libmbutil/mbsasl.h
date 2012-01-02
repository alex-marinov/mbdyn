/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#ifndef mbsasl_h
#define mbsasl_h

#ifdef HAVE_SASL2

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* used for negotiation buffers; should be enough for every known mechanism */
#define MBDYN_SASL_BUFSIZE		1024

#define MBDYN_SASL_SERVICE		"mbdyn"
#define MBDYN_SASL_CONFFILE		"mbdyn"

struct mbdyn_sasl_t {
	int		use_sasl;
#define MBDYN_SASL_NONE	0
#define MBDYN_SASL_SERVER	1
#define MBDYN_SASL_CLIENT	2
	unsigned	sasl_flags;
#define MBDYN_SASL_FLAG_NONE		0x0000
#define MBDYN_SASL_FLAG_CRITICAL	0x0001
#define MBDYN_SASL_FLAG_USERAUTHZ	0x0002
#define MBDYN_SASL_FLAG_INTERACT	0x0004
	unsigned long	sasl_usleep;	/* 0: forever */

	const char	*sasl_mech;	/* preferred; NULL -> all available */
	const char	*sasl_user;	/* if NULL, prompt */
	const char	*sasl_cred;	/* if NULL, prompt */
	const char	*sasl_realm;	/* if NULL, prompt? */
	const char	*sasl_authz;	/* if NULL, prompt? */

	const char	*sasl_hostname;	/* if NULL? */
	const char	*sasl_local_ip;	/* NULL is legal */
	const char	*sasl_remote_ip;/* NULL is legal */
};

#define MBDYN_SASL_INIT \
	{ MBDYN_SASL_NONE, MBDYN_SASL_FLAG_NONE, 0L, \
		NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }

/* init client data; urtnet_sasl is filled by mbdyn_sasl_parse_args() */
extern int
mbdyn_sasl_client_init(struct mbdyn_sasl_t *urtnet_sasl);

/* init server data; urtnet_sasl is filled by mbdyn_sasl_parse_args() */
extern int
mbdyn_sasl_server_init(struct mbdyn_sasl_t *urtnet_sasl);

/* init data; server/client is decided by urtnet_sasl;
 * urtnet_sasl is filled by mbdyn_sasl_parse_args() */
extern int
mbdyn_sasl_init(struct mbdyn_sasl_t *urtnet_sasl);

/* cleanup sasl session */
extern int
mbdyn_sasl_fini(void);

/* perform client auth on sock; urtnet_sasl is filled 
 * by mbdyn_sasl_parse_args() and client must be init'ed
 * by mbdyn_sasl_client_init() */
extern int
mbdyn_sasl_client_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *urtnet_sasl);

/* perform server auth on sock; urtnet_sasl is filled 
 * by mbdyn_sasl_parse_args() and server must be init'ed
 * by mbdyn_sasl_server_init() */
extern int
mbdyn_sasl_server_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *urtnet_sasl);

/* perform auth on sock; urtnet_sasl is filled 
 * by mbdyn_sasl_parse_args() and server/client must be init'ed
 * by mbdyn_sasl_init() */
extern int
mbdyn_sasl_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *urtnet_sasl);

/* validates data; server/client is decided by urtnet_sasl;
 * urtnet_sasl is filled by mbdyn_sasl_parse_args() */
extern int
mbdyn_sasl_validate(struct mbdyn_sasl_t *urtnet_sasl);

#define MBDYN_SASL_OPTIONS	"a:f:h:i:l:m:r:s:u:w:"

/* parses one arg in "opt" based on value in "val";
 * use MBDYN_SASL_OPTIONS in getopt for direct parsing of options,
 * or use "x:" (x arbitrary option) and then feed 
 * mbdyn_sasl_parse_args() with opt = optarg[0] and val = &optarg[2]
 * after checking that optarg[1] == '=' */
extern int
mbdyn_sasl_parse_args(int opt, const char *val,
		struct mbdyn_sasl_t *urtnet_sasl);

/* negotiates the entire set of sockets available and alredy bound,
 * but still in blocking mode ... */
extern int
mbdyn_sasl_negotiate(struct mbdyn_sasl_t *urtnet_sasl);

extern sasl_log_t *log_server_f;
extern sasl_log_t *log_client_f;
extern sasl_getsimple_t *get_user_f;
extern sasl_getsimple_t *get_authname_f;
extern sasl_getsecret_t *get_secret_f;
extern sasl_getrealm_t *get_realm_f;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* HAVE_SASL2 */

#endif /* mbsasl_h */
