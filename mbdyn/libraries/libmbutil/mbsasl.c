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

/*
 * The protocol is very simple.  It is client-initiated (of course),
 * by sending
 *
 * C: M
 * S: M<length><mechanism list>
 * C: S<length><mechanism chosen><length><additional data>
 * 
 * S: 	C<length><additional data>
 * C: 	C<length><additional data>
 *
 * S: O | F
 *
 * where M means "Methods", C means "Continuation" (with data), 
 * O means "OK" and F means "Fail"
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef HAVE_SASL2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <signal.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/poll.h>
#include <sys/socket.h>

#include <sasl/sasl.h>
#include "mbsasl.h"

/* global vars that contain the auth data */
static struct mbdyn_sasl_t	global_mbdyn_sasl = MBDYN_SASL_INIT;

static int
get_secret_default_f(sasl_conn_t *conn,
		void *context,
		int id,
		sasl_secret_t **psecret)
{
	size_t		credlen;
	const char	*cred = global_mbdyn_sasl.sasl_cred;

	if (id != SASL_CB_PASS) return SASL_FAIL;

	if (cred == NULL) {
		cred = getpass("password: ");
	}

	credlen = strlen(cred);
	*psecret = malloc(sizeof(sasl_secret_t) + credlen);
	(*psecret)->len = credlen;
	memcpy(&(*psecret)->data[0], cred, credlen);

	return SASL_OK;
}

static int
get_authname_default_f(void *context,
		int id,
		const char **result,
		unsigned *len)
{
	const char	*authz = global_mbdyn_sasl.sasl_authz;
	char		buf[MBDYN_SASL_BUFSIZE];

	if (id != SASL_CB_AUTHNAME) return SASL_FAIL;

	if (!authz && (global_mbdyn_sasl.sasl_flags & MBDYN_SASL_FLAG_USERAUTHZ)) {
		authz = global_mbdyn_sasl.sasl_user;
	}

	if (!authz) {
		size_t	buflen;
		
		printf("authzname: ");
		if (fgets(buf, sizeof(buf), stdin) == NULL) {
			return SASL_FAIL;
		}
		
		buflen = strlen(buf);
		if (buf[buflen - 1] == '\n') {
			--buflen;
			buf[buflen] = '\0';
		}

		if (buflen > 0) {
			authz = buf;
		}
	}

	if (authz) {
		*result = strdup(authz);
		if (len) {
			*len = strlen(*result);
		}
	} else {
		*result = NULL;
		if (len) {
			*len = 0;
		}
	}
	
	return SASL_OK;
}

static int
get_user_default_f(void *context,
		int id,
		const char **result,
		unsigned *len)
{
	const char	*user = global_mbdyn_sasl.sasl_user;
	char		buf[MBDYN_SASL_BUFSIZE];

	if (id != SASL_CB_USER) return SASL_FAIL;

	if (!user) {
		size_t	buflen;
		
		printf("username: ");
		if (fgets(buf, sizeof(buf), stdin) == NULL) {
			return SASL_FAIL;
		}
		
		buflen = strlen(buf);
		if (buf[buflen - 1] == '\n') {
			--buflen;
			buf[buflen] = '\0';
		}
		user = buf;
	}
		
	*result = strdup(user);
	if (len) {
		*len = strlen(*result);
	}
	
	return SASL_OK;
}

static int
get_realm_default_f(void *context,
		int id,
		const char **availrealms,
		const char **result)
{
	const char	*realm = global_mbdyn_sasl.sasl_realm;
	char		buf[MBDYN_SASL_BUFSIZE];

	if (id != SASL_CB_GETREALM) return SASL_FAIL;

	if (!realm) {
		size_t	buflen;
		
		printf("realm: ");
		if (fgets(buf, sizeof(buf), stdin) == NULL) {
			return SASL_FAIL;
		}
		
		buflen = strlen(buf);
		if (buf[buflen - 1] == '\n') {
			--buflen;
			buf[buflen] = '\0';
		}
		
		if (buflen > 0) {
			realm = buf;
		}
	}

	if (realm) {
		*result = strdup(realm);
	} else {
		*result = NULL;
	}

	return SASL_OK;
}

static int
log_server_default_f(void *context,
		int level,
		const char *message)
{
	fprintf(stderr, "[server %d] %s\n", level, message);
	return 0;
}

static int
log_client_default_f(void *context,
		int level,
		const char *message)
{
	fprintf(stderr, "[client %d] %s\n", level, message);
	return 0;
}

sasl_log_t		*log_server_f = &log_server_default_f;
sasl_log_t		*log_client_f = &log_client_default_f;
sasl_getsimple_t	*get_user_f = &get_user_default_f;
sasl_getsimple_t	*get_authname_f = &get_authname_default_f;
sasl_getsecret_t	*get_secret_f = &get_secret_default_f;
sasl_getrealm_t		*get_realm_f = &get_realm_default_f;

static sasl_callback_t server_callbacks[] = {
	{ SASL_CB_LOG, 		NULL /* log_server_f */ ,	NULL },
	{ SASL_CB_LIST_END,	NULL,				NULL }
};

static sasl_callback_t client_callbacks[] = {
	{ SASL_CB_GETREALM,	NULL /* get_realm_f */ ,	NULL },
	{ SASL_CB_USER, 	NULL /* get_user_f */ ,		NULL },
	{ SASL_CB_AUTHNAME, 	NULL /* get_authname_f */ ,	NULL },
	{ SASL_CB_PASS,		NULL /* get_secret_f */ ,	NULL },
	{ SASL_CB_LOG,		NULL /* log_client_f */ ,	NULL },
	{ SASL_CB_LIST_END,	NULL,				NULL }
};

int
mbdyn_sasl_client_init(struct mbdyn_sasl_t *mbdyn_sasl)
{
	int	rc;

	client_callbacks[0].proc = get_realm_f;
	client_callbacks[1].proc = get_user_f;
	client_callbacks[2].proc = get_authname_f;
	client_callbacks[3].proc = get_secret_f;
	client_callbacks[4].proc = log_client_f;

	if (mbdyn_sasl) {
		global_mbdyn_sasl = *mbdyn_sasl;
	}

	rc = sasl_client_init(client_callbacks); /* Callbacks supported */
	if (rc != SASL_OK) {
		printf("[client] sasl_client_init() failed: %d (%s)\n",
				rc, sasl_errstring(rc, NULL, NULL));
	}

	return rc;
}

int
mbdyn_sasl_server_init(struct mbdyn_sasl_t *mbdyn_sasl)
{
	int	rc;

	server_callbacks[0].proc = log_server_f;

	rc = sasl_server_init(server_callbacks,	/* Callbacks supported */
			MBDYN_SASL_CONFFILE);	/* Name of the application,
						 * used to look for rtai.conf
						 * file in /usr/lib/sasl2/ */ 
	if (rc != SASL_OK) {
		printf("[server] sasl_server_init() failed: %d (%s)\n",
				rc, sasl_errstring(rc, NULL, NULL));
	}

	return rc;
}

int
mbdyn_sasl_init(struct mbdyn_sasl_t *mbdyn_sasl)
{
	switch (mbdyn_sasl->use_sasl) {
	case MBDYN_SASL_SERVER:
		return mbdyn_sasl_server_init(mbdyn_sasl);
		break;
			
	case MBDYN_SASL_CLIENT:
		return mbdyn_sasl_client_init(mbdyn_sasl);
		break;

	default:
		/* if SASL is critical, return SASL_FAIL */
		if (mbdyn_sasl->sasl_flags & MBDYN_SASL_FLAG_CRITICAL) {
			return SASL_FAIL;
		}
		return SASL_OK;
	}
}

int
mbdyn_sasl_fini(void)
{
	sasl_done();

	return SASL_OK;
}

int
mbdyn_sasl_client_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *mbdyn_sasl)
{
	int		rc = 0;
	sasl_conn_t	*conn = NULL;
	sasl_interact_t	*client_interact = NULL;
	const char	*out = NULL, *mechusing = NULL;
	unsigned	outlen = 0;
#if 0	/* will be used later */
	sasl_security_properties_t secprops;
#endif

	int		count = 0;
	char		op;
	char		serverin[MBDYN_SASL_BUFSIZE] = { '\0' };
	const char		*mech = NULL;
	unsigned	serverinlen = 0;

	int		oldflags, currflags;

	oldflags = fcntl (sock, F_GETFL, 0);
	if (oldflags == -1) {
		return SASL_FAIL;
	}
	currflags = oldflags & ~O_NONBLOCK;
	if (fcntl(sock, F_SETFL, currflags) == -1) {
		return SASL_FAIL;
	}
	
	rc = sasl_client_new(MBDYN_SASL_SERVICE,	/* name of service */
			mbdyn_sasl->sasl_hostname,
					/* The fully qualified domain name
					 * of the server we're connecting to */
			mbdyn_sasl->sasl_local_ip,
			mbdyn_sasl->sasl_remote_ip,
					/* Local and remote IP address strings
					 * (NULL disables mechanisms which
					 * require this info)*/
			NULL,		/* connection-specific callbacks */
			0, 		/* security flags */
			&conn);		/* allocated on success */
	if (rc != SASL_OK) {
		printf("[client] sasl_client_new() failed: %d (%s)\n",
				rc, sasl_errstring(rc, NULL, NULL));
		goto cleanup;
	}

retry:;
	count = write(sock, "M", 1);
	if (count != 1) {
		if (count == -1) {
			if (errno == EINVAL && mbdyn_sasl->sasl_usleep) {
				usleep(mbdyn_sasl->sasl_usleep);
				goto retry;
			}
			printf("[client] write() failed errno=%d\n", errno);
		} else {
			printf("[client] write() failed (%d instead of %d)\n", count , 1);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	count = read(sock, &op, 1);
	if (count != 1) {
		if (count == -1) {
			printf("[client] read(M) failed errno=%d\n", errno);
		} else {
			printf("[client] read(M) failed (%d instead of %d)\n", count, 1);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	if (op != 'M') {
		if (op != 'A') {	/* 'A' means abort */
			printf("[client] expecting \"M\" with methods\n");
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	count = read(sock, &serverinlen, sizeof(unsigned));
	if (count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[client] read(Ml) failed errno=%d\n", errno);
		} else {
			printf("[client] read(Ml) failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	if (serverinlen >= sizeof(serverin)) {
		printf("[client] buffer for mechanism list is too small\n");
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = read(sock, serverin, serverinlen);
	if (count < 0 || (unsigned)count != serverinlen) {
		if (count < 0) {
			printf("[client] read(Md) failed errno=%d\n", errno);
		} else {
			printf("[client] read(Md) failed (%d instead of %d)\n",
				count, serverinlen);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	serverin[count] = '\0';

	mech = mbdyn_sasl->sasl_mech;
	if (!mech) {
		mech = serverin;
	}

#if 0	/* will be used later */
	memset(&secprops, 0, sizeof(secprops));
	secprops.min_ssf = 1;
	secprops.max_ssf = 256;
	secprops.maxbufsize = MBDYN_SASL_BUFSIZE;

	secprops.property_names = NULL;
	secprops.property_values = NULL;
	secprops.security_flags = SASL_SEC_NOANONYMOUS | SASL_SEC_NOPLAINTEXT | SASL_SEC_MUTUAL_AUTH; /* as appropriate */

	sasl_setprop(conn, SASL_SEC_PROPS, &secprops);
#endif

	do {
		char buf[MBDYN_SASL_BUFSIZE];

		rc = sasl_client_start(conn,	/* same context from above */
				mech,		/* the list of mechanisms
						 * from the server */
				&client_interact,	/* filled in if an
							 * interaction
							 * is needed */
				&out,		/* filled in on success */
				&outlen,	/* filled in on success */
				&mechusing);	/* selected mech */

		if (rc == SASL_INTERACT) {
			if (!(mbdyn_sasl->sasl_flags & MBDYN_SASL_FLAG_INTERACT)) {
				printf("[client] interaction disabled\n");
				rc = SASL_FAIL;
				goto cleanup;
			}

			/* FIXME: handle interaction ... */
			printf("[client] sasl_client_start() requested "
					"interaction \"%s\"; mech=%s\n",
					client_interact->prompt, mechusing);
		}

		if (outlen >= sizeof(buf)) {
			printf("[client] buf is too small\n");
			rc = SASL_FAIL;
			goto cleanup;
		}
		memcpy(buf, out, outlen);
		buf[outlen] = '\0';
	} while (rc == SASL_INTERACT);

	switch (rc) {
	case SASL_FAIL:
		goto cleanup;
		
	case SASL_OK:
	case SASL_CONTINUE:
		break;

	default:
		printf("[client] sasl start failure\n");
		rc = SASL_FAIL;
		goto cleanup;
	}

	count = write(sock, "S", 1);
	if (count != 1) {
		if (count == -1) {
			printf("[client] write() failed errno=%d\n", errno);
		} else {
			printf("[client] write() failed (%d instead of %d)\n", count, 1);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	serverinlen = strlen(mechusing);
	count = write(sock, &serverinlen, sizeof(unsigned));
	if (count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[client] write() failed errno=%d\n", errno);
		} else {
			printf("[client] write() failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = write(sock, mechusing, serverinlen);
	if (count < 0 || (unsigned)count != serverinlen) {
		if (count == -1) {
			printf("[client] mech=\"%s\" write failed errno=%d\n",
				mechusing, errno);
		} else {
			printf("[client] mech=\"%s\" write failed "
				"(%d instead of %d)\n",
				mechusing, count, serverinlen);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	count = write(sock, &outlen, sizeof(unsigned));
	if (count < 0 || (unsigned)count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[client] write() failed errno=%d\n", errno);
		} else {
			printf("[client] write() failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = write(sock, out ? out : "", outlen);
	if (count < 0 || (unsigned)count != outlen) {
		if (count == -1) {
			printf("[client] mech=\"%s\" write failed errno=%d\n",
				mechusing, errno);
		} else {
			printf("[client] mech=\"%s\" write failed "
				"(%d instead of %d)\n",
				mechusing, count, outlen);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	do {
		serverin[0] = '\0';
		serverinlen = 0;

		count = read(sock, &op, 1);
		if (count != 1) {
			if (count == -1) {
				printf("[client] read(c) failed errno=%d\n", errno);
			} else {
				printf("[client] read(c) failed (%d instead of %d)\n", count, 1);
			}
			rc = SASL_FAIL;
			goto cleanup;
		}

		switch (op) {
		case 'O':	/* ok; but continue, to check the library's response ... */
			if (rc == SASL_OK) {
				goto cleanup;
			}
			break;

		case 'F':	/* fail */
		case 'A':	/* abort */
			rc = SASL_FAIL;
			goto cleanup;

		case 'C':	/* continue ... */
			count = read(sock, serverin, sizeof(unsigned));
			if (count != sizeof(unsigned)) {
				if (count == -1) {
					printf("[client] read(Cd) failed errno=%d\n", errno);
				} else {
					printf("[client] read(Cd) failed "
						"(%d instead of %u)\n",
						count, (unsigned)sizeof(unsigned));
				}
				rc = SASL_FAIL;
				goto cleanup;
			}
			serverinlen = ((unsigned *)(serverin))[0];
			if (serverinlen >= sizeof(serverin)) {
				printf("[client] buffer for "
						"sasl_client_step() "
						"continuation data "
						"is too small\n");
				rc = SASL_FAIL;
				goto cleanup;
			}
			count = read(sock, serverin, serverinlen);
			if (count < 0 || (unsigned)count != serverinlen) {
				if (count == -1) {
					printf("[client] read() failed errno=%d\n", errno);
				} else {
					printf("[client] read() failed (%d instead of %d)\n",
						count, serverinlen);
				}
				rc = SASL_FAIL;
				goto cleanup;
			}
			serverin[serverinlen] = '\0';
			break;

		default:
			break;
		}

		rc = sasl_client_step(conn,	/* our context */
				serverin,	/* the data from the server */
				serverinlen,	/* its length */
				&client_interact,	/* this should be
							 * unallocated 
							 * and NULL */
				&out,		/* filled in on success */
				&outlen);	/* filled in on success */ 

		count = write(sock, "C", 1);
		if (count != 1) {
			if (count == -1) {
				printf("[client] write failed errno=%d\n", errno);
			} else {
				printf("[client] write failed (%d instead of %d)\n", count, 1);
			}
			rc = SASL_FAIL;
			goto cleanup;
		}

		count = write(sock, &outlen, sizeof(unsigned));
		if (count != sizeof(unsigned)) {
			if (count == -1) {
				printf("[client] write failed errno=%d\n", errno);
			} else {
				printf("[client] write failed "
					"(%d instead of %u)\n",
					count, (unsigned)sizeof(unsigned));
			}
			rc = SASL_FAIL;
			goto cleanup;
		}
		count = write(sock, out ? out : "", outlen);
		if (count < 0 || (unsigned)count != outlen) {
			if (count == -1) {
				printf("[client] write failed errno=%d\n", errno);
			} else {
				printf("[client] write failed (%d instead of %d)\n", count, outlen);
			}
			rc = SASL_FAIL;
			goto cleanup;
		}

		switch (rc) {
		case SASL_INTERACT:
			/* FIXME: interaction */
		case SASL_CONTINUE:
			break;
			
		case SASL_OK:
			break;

		case SASL_FAIL:
			goto cleanup;

		default:
			break;
		}
	} while (rc == SASL_INTERACT || rc == SASL_CONTINUE || rc == SASL_OK);

	if (rc != SASL_OK) {
		printf("[client] sasl step failure\n");
		goto cleanup;
	}

cleanup:;

#if 0	/* will use later */
	if (rc == SASL_OK) {
		int i;
		char buf[MBDYN_SASL_BUFSIZE];
		int ssf = 0, outbufsize = MBDYN_SASL_BUFSIZE;
		const int *ssfp = &ssf, *outbufsizep = &outbufsize;

		rc = sasl_getprop(conn, SASL_SSF, (const void **)&ssfp);
		if (rc != SASL_OK) {
			fprintf(stderr, "[client] unable to retrieve ssf\n");
			exit(EXIT_FAILURE);
		}
		if (*ssfp > 0) {
			ssf = *ssfp;
			rc = sasl_getprop(conn, SASL_MAXOUTBUF, (const void **)&outbufsizep);
			if (rc != SASL_OK) {
				fprintf(stderr, "[client] unable to retrieve out buf size\n");
				exit(EXIT_FAILURE);
			}
			outbufsize = *outbufsizep;
		}

		fprintf(stderr, "[client] SASL_OK (ssf=%d, outbufsize=%d)\n", ssf, outbufsize);
	}
#endif

	/* FIXME: don't dispose if using encoding */
	if (conn) {
		sasl_dispose(&conn);
	}

	if (fcntl (sock, F_SETFL, oldflags) == -1) {
		rc = SASL_FAIL;
	}

	return rc;
}

int
mbdyn_sasl_server_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *mbdyn_sasl)
{
	int		rc = SASL_FAIL;
	sasl_conn_t	*conn = NULL;
	const char	*result_string = NULL;
	unsigned	string_length = 0;
	int		number_of_mechanisms = 0;
	const char	*out = NULL;
	unsigned	outlen = 0;
#if 0	/* will be used later */
	sasl_security_properties_t secprops;
#endif
	
	int		count = 0;
	char		op = '\0';
	char		mechanism_client_chose[MBDYN_SASL_BUFSIZE] = { '\0' };
	char		clientin[MBDYN_SASL_BUFSIZE] = { '\0' };
	unsigned	clientinlen = 0;

	int		oldflags, currflags;

	oldflags = fcntl (sock, F_GETFL, 0);
	if (oldflags == -1) {
		return SASL_FAIL;
	}
	currflags = oldflags & ~O_NONBLOCK;
	if (fcntl(sock, F_SETFL, currflags) == -1) {
		return SASL_FAIL;
	}
	
	rc = sasl_server_new(MBDYN_SASL_SERVICE,	/* name of service */
		NULL,			/* my fully qualified domain name; 
					 * NULL says use gethostname() */
		mbdyn_sasl->sasl_realm,/* The user realm used for password
		                         * lookups; NULL means default
					 * to serverFQDN.  Note: This does 
					 * not affect Kerberos */
		mbdyn_sasl->sasl_local_ip,
		mbdyn_sasl->sasl_remote_ip,	/* IP Address information
						 * strings */
		NULL,			/* Callbacks supported only
					 * for this connection */
		0,			/* security flags (security layers 
					   are enabled using security 
					   properties, separately) */
		&conn);			/* allocated on success */
	if (rc != SASL_OK) {
		printf("[server] sasl_server_new() failed: %d (%s)\n",
				rc, sasl_errstring(rc, NULL, NULL));
		goto cleanup;
	}

#if 0	/* will be used later */
	memset(&secprops, 0, sizeof(secprops));
	secprops.min_ssf = 1;
	secprops.max_ssf = 256;
	secprops.maxbufsize = MBDYN_SASL_BUFSIZE;

	secprops.property_names = NULL;
	secprops.property_values = NULL;
	secprops.security_flags = SASL_SEC_NOANONYMOUS | SASL_SEC_NOPLAINTEXT | SASL_SEC_MUTUAL_AUTH; /* as appropriate */

	sasl_setprop(conn, SASL_SEC_PROPS, &secprops);
#endif

retry:;
	count = read(sock, &op, 1);
	if (count != 1) {
		if (count == -1) {
			if (errno == EINVAL && mbdyn_sasl->sasl_usleep) {
				usleep(mbdyn_sasl->sasl_usleep);
				goto retry;
			}
			printf("[server] read failed errno=%d\n", errno);
		} else {
			printf("[server] read failed (%d instead of %d)\n", count, 2);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	if (op != 'M') {
		printf("[server] \"M\" expected\n");
		rc = SASL_FAIL;
		goto cleanup;
	}

	rc = sasl_listmech(conn,	/* The context for this connection */
			NULL,		/* not supported */
			"",		/* What to prepend the string with */
			" ",		/* What to separate mechanisms with */
			"",		/* What to append to the string */
			&result_string,	/* The produced string. */
			&string_length,	/* length of the string */
			&number_of_mechanisms); /* Number of mechanisms 
						 * in the string */
	if (rc != SASL_OK) {
		printf("[server] sasl_listmech() failed: %s\n",
				sasl_errdetail(conn));
		rc = SASL_FAIL;
		goto cleanup;
	}

	count = write(sock, "M", 1);
	if (count != 1) {
		if (count == -1) {
			printf("[server] write() failed errno=%d\n", errno);
		} else {
			printf("[server] write() failed (%d instead of %d)\n", count, 1);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = write(sock, &string_length, sizeof(unsigned));
	if (count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[server] write() failed errno=%d\n", errno);
		} else {
			printf("[server] write() failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = write(sock, result_string, string_length);
	if (count < 0 || (unsigned)count != string_length) {
		if (count == -1) {
			printf("[server] write() failed errno=%d\n", errno);
		} else {
			printf("[server] write() failed (%d instead of %d)\n",
				count, string_length);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}

	/* read client's choice and more */
	count = read(sock, &op, 1);
	if (count != 1) {
		if (count == -1) {
			printf("[server] read() failed errno=%d\n", errno);
		} else {
			printf("[server] read() failed (%d instead of %d)\n",
				count, 1);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	if (op != 'S') {
		/* send abort? */
		if (op != 'A') {
			printf("[server] expecting \"S\" to start auth\n");
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = read(sock, &clientinlen, sizeof(unsigned));
	if (count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[server] read() failed errno=%d\n", errno);
		} else {
			printf("[server] read() failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	if (clientinlen >= sizeof(mechanism_client_chose)) {
		printf("[server] buffer for mechanism too small\n");
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = read(sock, mechanism_client_chose, clientinlen);
	if (count < 0 || (unsigned)count != clientinlen) {
		if (count == -1) {
			printf("[server] read() failed errno=%d\n", errno);
		} else {
			printf("[server] read() failed (%d instead of %u)\n",
				count, clientinlen);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	mechanism_client_chose[clientinlen] = '\0';

	/* read client's additional string ... */
	count = read(sock, &clientinlen, sizeof(unsigned));
	if (count != sizeof(unsigned)) {
		if (count == -1) {
			printf("[server] read() failed errno=%d\n", errno);
		} else {
			printf("[server] read() failed (%d instead of %u)\n",
				count, (unsigned)sizeof(unsigned));
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	if (clientinlen >= sizeof(clientin)) {
		printf("[server] buffer for client optional string "
				"is too small\n");
		rc = SASL_FAIL;
		goto cleanup;
	}
	count = read(sock, clientin, clientinlen);
	if (count < 0 || (unsigned)count != clientinlen) {
		if (count == -1) {
			printf("[server] read() failed errno=%d\n", errno);
		} else {
			printf("[server] read() failed (%d instead of %u)\n",
				count, clientinlen);
		}
		rc = SASL_FAIL;
		goto cleanup;
	}
	clientin[clientinlen] = '\0';

	rc = sasl_server_start(conn,		/* context */
			mechanism_client_chose,	/* selected mechanism */
			clientinlen ? clientin : NULL,	/* the optional string
							 * the client gave us */
			clientinlen,			/* and its length */
			&out,		/* The output of the library.
					 * Might not be NULL terminated */
			&outlen);	/* its lenght */

	do {
		switch (rc) {
		case SASL_OK:
			count = write(sock, "O", 1);
			if (count != 1) {
				if (count == -1) {
					printf("[server] write() failed errno=%d\n", errno);
				} else {
					printf("[server] write() failed (%d instead of %d)\n", count, 1);
				}
				rc = SASL_FAIL;
			}
			goto cleanup;

		case SASL_CONTINUE: {
			char buf[MBDYN_SASL_BUFSIZE];

			if (1 + outlen + sizeof(unsigned) >= sizeof(buf)) {
				printf("[server] buffer for "
						"sasl_server_step() "
						"continuation data "
						"is too small\n");
				rc = SASL_FAIL;
				goto cleanup;
			}

			memcpy(buf, "C", 1);
			memcpy(buf + 1, &outlen, sizeof(unsigned));
			memcpy(buf + 1 + sizeof(unsigned), out, outlen);
			
			count = write(sock, buf, 1 + outlen + sizeof(unsigned));
			if (count < 0 || (unsigned)count != 1 + outlen + sizeof(unsigned)) {
				if (count == -1) {
					printf("[server] write() failed errno=%d\n", count);
				} else {
					printf("[server] write() failed "
						"(%d instead of %u)\n",
						count,
						(unsigned)(1 + outlen + sizeof(unsigned)));
				}
				rc = SASL_FAIL;
				goto cleanup;
			}
			}
			break;

		default:
			count = write(sock, "F", 1);
			if (count != 1) {
				if (count == -1) {
					printf("[server] write() failed errno=%d\n", errno);
				} else {
					printf("[server] write() failed (%d instead of %d)\n", count, 1);
				}
				rc = SASL_FAIL;
			}
			goto cleanup;
		}

		count = read(sock, &op, 1);
		if (count != 1) {
			if (count == -1) {
				printf("[server] read() failed errno=%d\n", errno);
			} else {
				printf("[server] read() failed (%d instead of %d)\n", count, 1);
			}
			rc = SASL_FAIL;
			goto cleanup;
		}
		if (op != 'C') {
			if (op != 'A') {
				printf("[server] expecting \"C\" to continue\n");
			}
			rc = SASL_FAIL;
			goto cleanup;
		}

		clientin[0] = '\0';
		clientinlen = 0;

		count = read(sock, &clientinlen, sizeof(unsigned));
		if (count != sizeof(unsigned)) {
			if (count == -1) {
				printf("[server] read() failed errno=%d\n", errno);
			} else {
				printf("[server] read() failed "
					"(%d instead of %u)\n",
					count, (unsigned)sizeof(unsigned));
			}
			rc = SASL_FAIL;
			goto cleanup;
		}
		if (clientinlen >= sizeof(clientin)) {
			printf("[server] buffer for sasl_server_step() "
					"client data is too small\n");
			rc = SASL_FAIL;
			goto cleanup;
		}
		count = read(sock, clientin, clientinlen);
		if (count < 0 || (unsigned)count != clientinlen) {
			if (count == -1) {
				printf("[server] read() failed errno=%d\n", errno);
			} else {
				printf("[server] read() failed "
					"(%d instead of %u)\n",
					count, clientinlen);
			}
			rc = SASL_FAIL;
			goto cleanup;
		}

		rc = sasl_server_step(conn,	/* context */
				clientin,	/* what the client gave */
				clientinlen,	/* its length */
				&out,		/* allocated by library
						 * on success.  Might not 
						 * be NULL terminated */
				&outlen);	/* its length */
	} while (1);

cleanup:;
	/* FIXME: don't dispose if encoding connection */
	if (conn) {
		sasl_dispose(&conn);
	}

	if (fcntl (sock, F_SETFL, oldflags) == -1) {
		rc = SASL_FAIL;
	}

	return rc;
}

int
mbdyn_sasl_auth(int sock, struct sockaddr *bindaddr,
		struct mbdyn_sasl_t *mbdyn_sasl)
{
	switch (mbdyn_sasl->use_sasl) {
	case MBDYN_SASL_SERVER:
		return mbdyn_sasl_server_auth(sock, bindaddr, mbdyn_sasl);
		break;
			
	case MBDYN_SASL_CLIENT:
		return mbdyn_sasl_client_auth(sock, bindaddr, mbdyn_sasl);
		break;

	default:
		/* if SASL is critical, return SASL_FAIL */
		if (mbdyn_sasl->sasl_flags & MBDYN_SASL_FLAG_CRITICAL) {
			return SASL_FAIL;
		}
		return SASL_OK;
	}
}

int
mbdyn_sasl_validate(struct mbdyn_sasl_t *mbdyn_sasl)
{
	switch (mbdyn_sasl->use_sasl) {
	case MBDYN_SASL_NONE:
		return SASL_FAIL;

	case MBDYN_SASL_SERVER:
		if (mbdyn_sasl->sasl_user != NULL) {
			return SASL_FAIL;
		}

		if (mbdyn_sasl->sasl_cred != NULL) {
			return SASL_FAIL;
		}

		if (mbdyn_sasl->sasl_authz != NULL) {
			return SASL_FAIL;
		}

		return SASL_OK;

	case MBDYN_SASL_CLIENT:
		if (mbdyn_sasl->sasl_hostname == NULL) {
			mbdyn_sasl->sasl_hostname = "localhost";
		}

		return SASL_OK;

	}
	
	return SASL_FAIL;
}

/*
 * a=<authz>	client: authorization identity (optional)
 * f=<flag>[=<value>]
 * i=<remoteip>	remote ip
 * l=<localip>	local ip
 * m=<method>	(list of) acceptable method(s)
 * r=<realm>	client: user realm;
 * 		server: server realm
 * s={server|client}	use SASL to negotiate auth 
 * 			as server or client
 * u=<user>	client: user identity
 * w=<cred>	client: user credential
 */
int
mbdyn_sasl_parse_args(int opt, const char *optarg, 
		struct mbdyn_sasl_t *mbdyn_sasl)
{
	switch (opt) {
	case 'a':
		mbdyn_sasl->sasl_authz = optarg[0] ? optarg : NULL;
		break;

	case 'f':
		if (strcasecmp(optarg, "userauthz") == 0) {
			mbdyn_sasl->sasl_flags |= MBDYN_SASL_FLAG_USERAUTHZ;

		} else if (strcasecmp(optarg, "interact") == 0) {
			mbdyn_sasl->sasl_flags |= MBDYN_SASL_FLAG_INTERACT;

		} else {
			printf("UNKNOWN FLAG '%s'\n", optarg);
		}
		break;

	case 'h':
		mbdyn_sasl->sasl_hostname = optarg[0] ? optarg : NULL;
		break;

	case 'i':
		mbdyn_sasl->sasl_remote_ip = optarg[0] ? optarg : NULL;
		break;

	case 'l':
		mbdyn_sasl->sasl_local_ip = optarg[0] ? optarg : NULL;
		break;

	case 'm':
		mbdyn_sasl->sasl_mech = optarg[0] ? optarg : NULL;
		break;

	case 'r':
		mbdyn_sasl->sasl_realm = optarg[0] ? optarg : NULL;
		break;

	case 's':
		if (strcasecmp(optarg, "server") == 0) {
			mbdyn_sasl->use_sasl = MBDYN_SASL_SERVER;
		} else if (strcasecmp(optarg, "client") == 0) {
			mbdyn_sasl->use_sasl = MBDYN_SASL_CLIENT;
		} else if (strcasecmp(optarg, "none") == 0) {
			mbdyn_sasl->use_sasl = MBDYN_SASL_NONE;
		} else {
			printf("UNKNOWN SASL MODE; SASL DISABLED\n");
			mbdyn_sasl->use_sasl = MBDYN_SASL_NONE;
		}
		break;

	case 't': {
		char		*next = NULL;
		unsigned long	l;

		errno = 0;
		l = strtoul(optarg, &next, 10);
		int save_errno = errno;
		if (next == NULL || next[0] != '\0') {
			printf("ILLEGAL SLEEP TIME '%s'\n", optarg);

		} else if (save_errno == ERANGE) {
			printf("SLEEP TIME '%s' OVERFLOWS\n", optarg);

		} else {
			mbdyn_sasl->sasl_usleep = l;
		}
		}
		break;

	case 'u':
		mbdyn_sasl->sasl_user = optarg[0] ? optarg : NULL;
		break;

	case 'w':
		mbdyn_sasl->sasl_cred = optarg[0] ? optarg : NULL;
		break;

	default:
		return -1;
	}

	return 0;
}

#endif /* HAVE_SASL2 */

