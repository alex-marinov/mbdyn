/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include <pwd.h>

#include <dataman.h>
#include <auth.h>

extern "C" char *crypt(const char *, const char *);

/* NoAuth - begin */

AuthMethod::AuthRes
NoAuth::Auth(const char * /* user */ , const char * /* cred */ ) const 
{
   return AuthMethod::AUTH_OK;
}

/* NoAuth - end */


/* PasswordAuth - begin */

static char *
make_salt(void)
{
   static char salt_charset[] = 
     "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789./";      
   static char salt[2];
   
   ASSERT(strlen(salt_charset) == 64);
   
   salt[0] = salt_charset[rand()%64];
   salt[1] = salt_charset[rand()%64];
   
   return salt;
}


PasswordAuth::PasswordAuth(const char *u, const char *c) 
{
   ASSERT(u != NULL);
   ASSERT(c != NULL);
   
   strncpy(User, u, 8);
   User[8] = '\0';
   
   strncpy(Cred, crypt(c, make_salt()), 13);
   Cred[13] = '\0';	
}


AuthMethod::AuthRes
PasswordAuth::Auth(const char *user, const char *cred) const 
{
   if (user == NULL || cred == NULL) {
      return AuthMethod::AUTH_ERR;
   }
   if (strcmp(User, user) == 0 && strcmp(Cred, crypt(cred, Cred)) == 0) {
      return AuthMethod::AUTH_OK;
   }
   return AuthMethod::AUTH_FAIL;
}

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
extern "C" char * strdup(const char *);
}

#define INPUTSIZE PAM_MAX_MSG_SIZE           /* maximum length of input+1 */
#define CONV_ECHO_ON  1                            /* types of echo state */
#define CONV_ECHO_OFF 0

static void pam_misc_conv_delete_binary(void **delete_me)
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

int (*pam_binary_handler_fn)(const void *send, void **receive) = NULL;
void (*pam_binary_handler_free)(void **packet_p) = pam_misc_conv_delete_binary;

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
	 if (!::fSilent) {
	    if (fprintf(stdout, "%s\n", msgm[count]->msg) < 0) {
	       goto failed_conversation;
	    }
	 }
	 break;
       case PAM_BINARY_PROMPT: {
	  void *pack_out = NULL;
	  const void *pack_in = msgm[count]->msg;
	  
	  if (!pam_binary_handler_fn
	      || pam_binary_handler_fn(pack_in, &pack_out) != PAM_SUCCESS
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

      if (string) {                         /* must add to reply array */
	 /* add string to list of responses */
	 
	 reply[count].resp_retcode = 0;
	 reply[count].resp = string;
	 string = NULL;
      }
   }
   
   /* New (0.59+) behavior is to always have a reply - this is
    compatable with the X/Open (March 1997) spec. */
   *response = reply;
   reply = NULL;
   
   return PAM_SUCCESS;
   
failed_conversation:
   
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
	    pam_binary_handler_free((void **) &reply[count].resp);
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
	 cerr << "cannot determine the effective user!" << endl;
	 THROW(ErrGeneric());
      }
      
      u = pw->pw_name;
   }
   
   SAFESTRDUP(User, u, DMmm);
   
   struct pam_conv conv;
   conv.conv = mbdyn_conv;
   conv.appdata_ptr = NULL;

   pam_handle_t *pamh = NULL;
   int retval = pam_start("mbdyn", User, &conv, &pamh);
   
   if (retval != PAM_SUCCESS) {      
      cerr << "user \"" << User << "\" cannot be authenticated " << endl;
      
      if (pam_end(pamh, retval) != PAM_SUCCESS) { 
	 cerr << "unable to release PAM authenticator" << endl;
      }
      
      THROW(ErrGeneric());
   }
   
   if (pam_end(pamh, retval) != PAM_SUCCESS) { 
      cerr << "unable to release PAM authenticator" << endl;
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
      cerr << "user \"" << user << "\" cannot be authenticated " << endl;
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
      cerr << "unable to release PAM authenticator" << endl;
   }
   
   return r;
}

#endif /* USE_PAM */

/* PAM_Auth - end */


AuthMethod* 
ReadAuthMethod(DataManager* /* pDM */ , MBDynParser& HP)
{
   AuthMethod* pAuth = NULL;
   
   const char* sKeyWords[] = {
      "noauth",
	"password",
	"pwdb",
	"pam"
   };
   
   enum KeyWords {
      UNKNOWN = -1,
	
	NOAUTH = 0,
	PASSWORD,
	PWDB,
	PAM,
	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);   
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
   switch (CurrKeyWord) {
      
      /* auth is always successful */
    case NOAUTH: {
       SAFENEW(pAuth, NoAuth, DMmm);
       break;
    }
      
      /* auth is based on user id and password;
       * 
       */
    case PASSWORD: {
       
       if (!HP.IsKeyWord("user")) {
	  cerr << "user expected at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }
       
       const char* tmp = HP.GetStringWithDelims();
       if (strlen(tmp) == 0) {
	  cerr << "Need a legal user id at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }
       
       char* user = NULL;
       SAFESTRDUP(user, tmp, DMmm);
       
       if (!HP.IsKeyWord("credentials")) {
	  cerr << "credentials expected at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }
       if (HP.IsKeyWord("prompt")) {
	  tmp = getpass("password: ");
       } else {
	  tmp = HP.GetStringWithDelims();
       }
       if (strlen(tmp) == 0) {
	  cerr << "Warning: null credentials at line " << HP.GetLineData()  << endl;
       }
       
       char* cred = NULL;
       SAFESTRDUP(cred, tmp, DMmm);
       memset((char *)tmp, '\0', strlen(tmp));
    
       SAFENEWWITHCONSTRUCTOR(pAuth,
			      PasswordAuth,
			      PasswordAuth(user, cred),
			      DMmm);
       SAFEDELETEARR(user, DMmm);
       memset(cred, '\0', strlen(cred));
       SAFEDELETEARR(cred, DMmm);	     
       
       break;
    }
      
    case PAM: {
#ifdef USE_PAM
       char* user = NULL;
       if (HP.IsKeyWord("user")) {
	  const char *tmp = HP.GetStringWithDelims();
	  if (strlen(tmp) == 0) {
	     cerr << "Need a legal user id at line " << HP.GetLineData() << endl;
	     THROW(ErrGeneric());
	  }
       	 
	  SAFESTRDUP(user, tmp, DMmm);
       }
       
       SAFENEWWITHCONSTRUCTOR(pAuth, PAM_Auth, PAM_Auth(user), DMmm);
       break;
#else /* !USE_PAM */       
       cerr << "line " << HP.GetLineData() 
	 << ": sorry, this system does not support PAM" << endl;
       THROW(ErrGeneric());
#endif /* !USE_PAM */       
    }
           
    case PWDB:
      cerr << "not implemented yet" << endl;
    default:
      THROW(ErrNotImplementedYet());
   }      
      
   return pAuth;
}

