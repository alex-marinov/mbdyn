/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <unistd.h>
#ifdef HAVE_PWD_H
#include <pwd.h>
#endif // HAVE_PWD_H

#include "dataman.h"
#include "auth.h"
#include "crypt.h"

/* NoAuth - begin */

/* from sockdrv.c */
extern int get_auth_token(FILE *fd, char *user, char *cred, char **nextline);

AuthMethod::AuthRes
NoAuth::Auth(const char * /* user */ , const char * /* cred */ ) const
{
	return AuthMethod::AUTH_OK;
}

AuthMethod::AuthRes
NoAuth::Auth(int /* sock */ ) const
{
	return AuthMethod::AUTH_OK;
}

/* NoAuth - end */


/* PasswordAuth - begin */

#ifdef HAVE_CRYPT

PasswordAuth::PasswordAuth(const char *u, const char *c, const char *salt_format)
{
	ASSERT(u != NULL);
	ASSERT(c != NULL);

	strncpy(User, u, sizeof(User));
	User[STRLENOF(User)] = '\0';

	char *tmp = 0;
	if (strncmp(c, "{CRYPT}", STRLENOF("{CRYPT}")) == 0) {
		tmp = &c[STRLENOF("{CRYPT}")];

		if (strlen(tmp) >= sizeof(Cred)) {
			silent_cerr("unable to handle credentials (too long)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		char salt[33];
		tmp = crypt(c, mbdyn_make_salt(salt, sizeof(salt), salt_format));
		if (tmp == NULL) {
			silent_cerr("crypt() failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	strncpy(Cred, tmp, sizeof(Cred));
	Cred[STRLENOF(Cred)] = '\0';
}

AuthMethod::AuthRes
PasswordAuth::Auth(const char *user, const char *cred) const
{
	if (user == NULL || cred == NULL) {
		return AuthMethod::AUTH_ERR;
	}

	char *tmp = crypt(cred, Cred);
	if (tmp == NULL) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (strcmp(User, user) == 0 && strcmp(Cred, tmp) == 0) {
		return AuthMethod::AUTH_OK;
	}

	return AuthMethod::AUTH_FAIL;
}

AuthMethod::AuthRes
PasswordAuth::Auth(int sock) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

#endif /* HAVE_CRYPT */

/* PasswordAuth - end */


/* PAM_Auth - begin */

#ifdef USE_PAM

extern "C" {
#include <security/pam_appl.h>
#ifdef HAVE_PAM_MISC_H
#include <pam_misc.h>
#elif HAVE_SECURITY_PAM_MISC_H
#include <security/pam_misc.h>
#endif /* HAVE_SECURITY_PAM_MISC_H */
#include <string.h>
}

#define INPUTSIZE PAM_MAX_MSG_SIZE           /* maximum length of input+1 */
#define CONV_ECHO_ON  1                            /* types of echo state */
#define CONV_ECHO_OFF 0

static void
pam_misc_conv_delete_binary(void **delete_me)
{
	if (delete_me && *delete_me) {
		unsigned char *packet = *(unsigned char **)delete_me;
		int length;

		length = 4+(packet[0]<<24)+(packet[1]<<16)+(packet[2]<<8)+packet[3];
		memset(packet, 0, length);
		free(packet);
		*delete_me = packet = NULL;
	}
}

int (*mb_pam_bh_fn)(const void *send, void **receive) = NULL;
void (*mb_pam_bh_free)(void **packet_p) = pam_misc_conv_delete_binary;

/*
 * Uso una funzione di conversazione che si limita a restituire il valore
 * di credenziali passato in appdata_ptr (cast a void*)
 */
int
mbdyn_conv(int num_msg, const struct pam_message **msgm,
	   struct pam_response **response, void *appdata_ptr)
{
	int count = 0;
	struct pam_response *reply;

	if (num_msg <= 0) {
		return PAM_CONV_ERR;
	}

	reply = (struct pam_response *) calloc(num_msg, sizeof(struct pam_response));
	if (reply == NULL) {
		return PAM_CONV_ERR;
	}

	for (count = 0; count < num_msg; ++count) {
		char *string = NULL;

		switch (msgm[count]->msg_style) {
		case PAM_PROMPT_ECHO_OFF:
		case PAM_PROMPT_ECHO_ON:
			string = (char *)x_strdup((char *)appdata_ptr);
			if (string == NULL) {
				goto failed_conversation;
			}
			break;

		case PAM_ERROR_MSG:
			if (fprintf(stderr, "%s\n", msgm[count]->msg) < 0) {
				goto failed_conversation;
			}
			break;

		case PAM_TEXT_INFO:
			if (::fSilent < 2) {
				if (fprintf(stdout, "%s\n", msgm[count]->msg) < 0) {
					goto failed_conversation;
				}
			}
			break;

		case PAM_BINARY_PROMPT: {
			void *pack_out = NULL;
			const void *pack_in = msgm[count]->msg;

			if (!mb_pam_bh_fn
					|| mb_pam_bh_fn(pack_in, &pack_out) != PAM_SUCCESS
					|| pack_out == NULL) {
				goto failed_conversation;
			}
			string = (char *) pack_out;
			pack_out = NULL;
			break;
		}
#if 0 /* non-standard message styles */
		case PAM_BINARY_MSG: {
			const void *pack_in = msgm[count]->msg;
			if (!pam_binary_handler_fn
					|| pam_binary_handler_fn(pack_in, NULL) != PAM_SUCCESS) {
				goto failed_conversation;
			}
			break;
		}
#endif /* non-standard message styles */

		default:
			fprintf(stderr, "erroneous conversation (%d)\n",
					msgm[count]->msg_style);
			goto failed_conversation;
		}

		if (string) {	/* must add to reply array */
			/* add string to list of responses */
			reply[count].resp_retcode = 0;
			reply[count].resp = string;
			string = NULL;
		}
	}

	/* New (0.59+) behavior is to always have a reply - this is
	 * compatable with the X/Open (March 1997) spec. */
	*response = reply;
	reply = NULL;

	return PAM_SUCCESS;

failed_conversation:;
	if (reply) {
		for (count = 0; count < num_msg; ++count) {
			if (reply[count].resp == NULL) {
				continue;
			}
			
			switch (msgm[count]->msg_style) {
			case PAM_PROMPT_ECHO_ON:
			case PAM_PROMPT_ECHO_OFF:
				_pam_overwrite(reply[count].resp);
				free(reply[count].resp);
				break;
			
			case PAM_BINARY_PROMPT:
				mb_pam_bh_free((void **) &reply[count].resp);
				break;
			
			case PAM_ERROR_MSG:
			case PAM_TEXT_INFO:
#if 0 /* non-standard message style */
			case PAM_BINARY_MSG:
#endif /* non-standard message style */
				/* should not actually be able to get here ... */
				free(reply[count].resp);
			}
			
			reply[count].resp = NULL;
		}

		/* forget reply too */
		free(reply);
		reply = NULL;
	}

	return PAM_CONV_ERR;
}

PAM_Auth::PAM_Auth(const char *u)
: User(NULL)
{
	if (u == NULL) {
		struct passwd* pw = getpwuid(getuid());

		if (pw == NULL) {
			silent_cerr("PAM_Auth: cannot determine the effective user!" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		u = pw->pw_name;
	}

	SAFESTRDUP(User, u);

	struct pam_conv conv;
	conv.conv = mbdyn_conv;
	conv.appdata_ptr = NULL;
	
	pam_handle_t *pamh = NULL;
	int retval = pam_start("mbdyn", User, &conv, &pamh);

	if (retval != PAM_SUCCESS) {
		silent_cerr("PAM_Auth: user \"" << User 
				<< "\" cannot be authenticated " << std::endl);

		if (pam_end(pamh, retval) != PAM_SUCCESS) {
			silent_cerr("PAM_Auth: unable to release PAM authenticator" << std::endl);
		}
		
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pam_end(pamh, retval) != PAM_SUCCESS) {
		silent_cerr("PAM_Auth: unable to release PAM authenticator" << std::endl);
	}
}

AuthMethod::AuthRes
PAM_Auth::Auth(const char *user, const char *cred) const
{
	pam_handle_t *pamh = NULL;
	int retval;

	AuthMethod::AuthRes r(AuthMethod::AUTH_UNKNOWN);

	if (user == NULL || cred == NULL) {
		return AuthMethod::AUTH_ERR;
	}

	if (strcmp(User, user) != 0) {
		silent_cerr("PAM_Auth::Auth: user \"" << user 
				<< "\" cannot be authenticated " << std::endl);
		return AuthMethod::AUTH_ERR;
	}

	struct pam_conv conv;
	conv.conv = mbdyn_conv;
	conv.appdata_ptr = (void*)cred;
	retval = pam_start("mbdyn", User, &conv, &pamh);
	if (retval == PAM_SUCCESS) {
		retval = pam_authenticate(pamh, 0);
		if (retval == PAM_SUCCESS) {
			r = AuthMethod::AUTH_OK;
		} else {
			r = AuthMethod::AUTH_FAIL;
		}
	} else {
		r = AuthMethod::AUTH_ERR;
	}

	if (pam_end(pamh, retval) != PAM_SUCCESS) {
		silent_cerr("PAM_Auth::Auth: unable to release PAM authenticator" << std::endl);
	}

	return r;
}

AuthMethod::AuthRes
PAM_Auth::Auth(int sock) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

#endif /* USE_PAM */

/* PAM_Auth - end */

#ifdef HAVE_SASL2

#if defined(HAVE_SASL_SASL_H)
#include <sasl/sasl.h>
#elif defined(HAVE_SASL_H)
#include <sasl.h>
#endif /* HAVE_SASL_SASL_H || HAVE_SASL_H */
#include "mbsasl.h"

int
mbdyn_sasl_log(void *context, int level, const char *message)
{
	switch (level) {
	case 0:
	case 1:
		std::cerr << "[mbdyn " << level << "] " << message << std::endl;
		break;

	case 2:
		silent_cerr("[mbdyn " << level << "] " << message << std::endl);
		break;

	case 3:
		std::cout << "[mbdyn " << level << "] " << message << std::endl;
		break;

	default:
		silent_cout("[mbdyn " << level << "] " << message << std::endl);
		break;
	}

	return 0;
}

SASL2_Auth::SASL2_Auth(const mbdyn_sasl_t *ms)
: mbdyn_sasl(*ms)
{
	/* server operations */
	log_server_f = mbdyn_sasl_log;

	if (mbdyn_sasl_init(&mbdyn_sasl) != SASL_OK) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

AuthMethod::AuthRes
SASL2_Auth::Auth(const char *user, const char *cred) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

AuthMethod::AuthRes
SASL2_Auth::Auth(int sock) const
{
	switch (mbdyn_sasl_auth(sock, NULL, &mbdyn_sasl)) {
	case SASL_OK:
		return AuthMethod::AUTH_OK;

	case SASL_FAIL:
		return AuthMethod::AUTH_FAIL;

	default:
		return AuthMethod::AUTH_ERR;
	}
}

#endif /* HAVE_SASL2 */


AuthMethod*
ReadAuthMethod(const DataManager* /* pDM */ , MBDynParser& HP)
{
	AuthMethod* pAuth = NULL;

	const char* sKeyWords[] = {
		"noauth",
		"password",
		"pwdb",
		"pam",
		"sasl",
		NULL
	};

	enum KeyWords {
		UNKNOWN = -1,

		NOAUTH = 0,
		PASSWORD,
		PWDB,
		PAM,
		SASL,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di drive */
	KeyWords CurrKeyWord = KeyWords(HP.GetWord());

	switch (CurrKeyWord) {
	/* auth is always successful */
	case NOAUTH:
		SAFENEW(pAuth, NoAuth);
		break;

	/* auth is based on user id and password */
	case PASSWORD: {
#ifdef HAVE_CRYPT
		if (!HP.IsKeyWord("user")) {
			silent_cerr("ReadAuthMethod: user expected at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const char* tmp = HP.GetStringWithDelims();
		if (strlen(tmp) == 0) {
			silent_cerr("ReadAuthMethod: Need a legal user id at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		char* user = 0;
		SAFESTRDUP(user, tmp);

		if (!HP.IsKeyWord("credentials")) {
			silent_cerr("ReadAuthMethod: credentials expected at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (HP.IsKeyWord("prompt")) {
			tmp = getpass("password: ");
		} else {
			tmp = HP.GetStringWithDelims();
		}

		if (strlen(tmp) == 0) {
			silent_cout("ReadAuthMethod: null credentials at line "
					<< HP.GetLineData() << std::endl);
		}

		char* cred = 0;
		SAFESTRDUP(cred, tmp);
		memset((char *)tmp, '\0', strlen(tmp));

		char *salt_format = 0;
		if (HP.IsKeyWord("salt" "format")) {
			tmp = HP.GetStringWithDelims();
			SAFESTRDUP(salt_format, tmp);
		}

		SAFENEWWITHCONSTRUCTOR(pAuth,
				PasswordAuth,
				PasswordAuth(user, cred, salt_format));
		SAFEDELETEARR(user);
		memset(cred, '\0', strlen(cred));
		SAFEDELETEARR(cred);
		if (salt_format) {
			SAFEDELETEARR(salt_format);
		}

		break;
#else /* !HAVE_CRYPT */
		silent_cerr("ReadAuthMethod: line " << HP.GetLineData()
				<< ": no working crypt(3)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !HAVE_CRYPT */
	}

	case PAM: {
#ifdef USE_PAM
		char* user = NULL;
		if (HP.IsKeyWord("user")) {
			const char *tmp = HP.GetStringWithDelims();
			if (strlen(tmp) == 0) {
				silent_cerr("ReadAuthMethod: Need a legal user id at line "
						<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			SAFESTRDUP(user, tmp);
		}

		SAFENEWWITHCONSTRUCTOR(pAuth, PAM_Auth, PAM_Auth(user));
		break;
#else /* !USE_PAM */
		silent_cerr("ReadAuthMethod: line " << HP.GetLineData()
			<< ": no PAM support" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_PAM */
	}

	case SASL: {
#ifdef HAVE_SASL2
		mbdyn_sasl_t	mbdyn_sasl = MBDYN_SASL_INIT;
		mbdyn_sasl.use_sasl = MBDYN_SASL_SERVER;

		if (HP.IsKeyWord("mechanism") || HP.IsKeyWord("mech")) {
			const char *s = HP.GetStringWithDelims();
			if (s != NULL) {
				SAFESTRDUP(mbdyn_sasl.sasl_mech, s);
			} else {
				silent_cerr("ReadAuthMethod: unable to get SASL mech at line "
						<< HP.GetLineData() << std::endl);
			}
		}
		
		SAFENEWWITHCONSTRUCTOR(pAuth, SASL2_Auth, SASL2_Auth(&mbdyn_sasl));
		break;
#else /* !HAVE_SASL2 */
		silent_cerr("ReadAuthMethod: line " << HP.GetLineData()
			<< ": no SASL2 support" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !HAVE_SASL2 */
	}

	case PWDB:
		silent_cerr("ReadAuthMethod: PWDB not implemented yet" << std::endl);

	default:
		throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
	}

	return pAuth;
}

