/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <termios.h>

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <string.h>

const unsigned short int PORT = 5555;
const char* SERVERHOST = "localhost";

static void
init_sockaddr (struct sockaddr_in *name,
	       const char *hostname,
	       unsigned short int port)
{
   struct hostent *hostinfo;
   
   name->sin_family = AF_INET;
   name->sin_port = htons (port);
   hostinfo = gethostbyname (hostname);
   if (hostinfo == NULL) {
      fprintf (stderr, "Unknown host %s.\n", hostname);
      exit (EXIT_FAILURE);
   }
   name->sin_addr = *(struct in_addr *) hostinfo->h_addr;
}

static void 
keys(FILE * fh)
{
   fprintf(fh, 
	   "\tkeys:\n"
	   "\t\t'i':\tswitches the incremental mode on\n"
	   "\t\t'p':\tincrements the drive\n"
	   "\t\t'm':\tdecrements the drive\n"
	   "\t\t'^D':\tquits\n\n");
}
   

static void
usage(void)
{
   fprintf(stderr,
	   "\n\tusage: autopilot [h:p:D:w:Wx:] label\n\n"
	   "\t\t-h host\t\thost name\n"
	   "\t\t-p port\t\tport number\n"
	   "\t\t-D user\t\tuser name\n"
	   "\t\t-w cred\t\tuser credentials\n"
	   "\t\t-W\t\tprompt for user credentials\n"
	   "\t\t-x value\tincrement\n\n"
	   "\tlabel:\tfile drive index to modify\n\n");
   keys(stderr);
}


/* Use this variable to remember original terminal attributes. */
struct termios saved_attributes;


void
reset_input_mode (void)
{
   tcsetattr(STDIN_FILENO, TCSANOW, &saved_attributes);
}
     

void
set_input_mode (void)
{
   struct termios tattr;
   
   /* Make sure stdin is a terminal. */
   if (!isatty(STDIN_FILENO)) {
      fprintf(stderr, "Not a terminal.\n");
      exit(EXIT_FAILURE);
   }
     
   /* Save the terminal attributes so we can restore them later. */
   tcgetattr(STDIN_FILENO, &saved_attributes);
   atexit(reset_input_mode);
   
   /* Set the funny terminal modes. */
   tcgetattr(STDIN_FILENO, &tattr);
   tattr.c_lflag &= ~(ICANON|ECHO); /* Clear ICANON and ECHO. */
   tattr.c_cc[VMIN] = 0;
   tattr.c_cc[VTIME] = 0;
   if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &tattr) < 0) {
      perror("tcsetattr");
      exit(EXIT_FAILURE);
   }
}


int
send_message(const char *host, unsigned short int port, const char *message)
{
   int sock;
   struct sockaddr_in server_name;
   FILE *fd;   
   
   /* Create the socket. */
   sock = socket (PF_INET, SOCK_STREAM, 0);
   if (sock < 0) {
      return -1;
   }
   
   /* Connect to the server. */
   init_sockaddr (&server_name, host, port);
   if (0 > connect (sock,
		    (struct sockaddr *) &server_name,
		    sizeof (server_name))) {
      return -1;
   }

   fd = fdopen(sock, "w");
   fputs(message, fd);
   fclose (fd);
   
   return 0;
}


int
main (int argc, char *argv[])
{   
   char *host = (char *)SERVERHOST;
   unsigned short int port = PORT;
   
   char *user = NULL;
   char *cred = NULL;
   
   char *increment = "1.";
   char *label = NULL;
   
   char c;
   
   char *auth = NULL;
   char *inc = NULL; 
   char *plus = NULL;
   char *minus = NULL;
   
   int verbose = 0;
     
   while (1) {
      int opt;
      
      opt = getopt (argc, argv, "h:p:D:w:Wx:v");
      
      if (opt == EOF) {
	 break;
      }
      
      switch (opt) {
       case 'h':
	 host = strdup (optarg);
	 break;
	 
       case 'p':
	 port = atoi (optarg);
	 break;
	 
       case 'D':
	 user = strdup (optarg);
	 break;
	 
       case 'w':
	 cred = strdup (optarg);
	 break;
	 
       case 'W': {
	  char *tmp = getpass("password: ");
	  if (tmp) {	     
	     cred = strdup(tmp);
	     memset(tmp, '\0', strlen(tmp));
	  }
	  break;
       }
	 
       case 'x':
	 increment = strdup(optarg);
	 break;
	 
       case 'v':
	 verbose++;
	 break;
      }
   }
   
   if (argc-optind < 1) {
      usage();
      exit(EXIT_SUCCESS);
   }
   
   label = argv[optind];

   
   /* messaggi: */
   if (user) {
      int ul = strlen(user);
      if (cred) {
	 int cl = strlen(cred);
	 auth = (char *)calloc(sizeof(char), ul+cl+7+11+1);
	 sprintf(auth, "user: %s\npassword: %s\n", user, cred);
      } else {
	 auth = (char *)calloc(sizeof(char), ul+7+1);
	 sprintf(auth, "user: %s\n", user);
      }
   }

   if (auth) {
      int l;
      
      l = strlen(auth)+8+strlen(label)+9+2;
      inc = (char *)calloc(sizeof(char), l+1);
      sprintf(inc, "%slabel: %s\ninc: yes\n.\n", auth, label);
      
      l = strlen(auth)+8+strlen(label)+8+strlen(increment)+2;
      plus = (char *)calloc(sizeof(char), l+1);
      minus = (char *)calloc(sizeof(char), l+1+1);
      sprintf(plus, "%slabel: %s\nvalue: %s\n.\n", auth, label, increment);
      sprintf(minus, "%slabel: %s\nvalue: -%s\n.\n", auth, label, increment);
      
   } else {
      int l;
      
      l = 8+strlen(label)+9+2;
      inc = (char *)calloc(sizeof(char), l+1);
      sprintf(inc, "label: %s\ninc: yes\n.\n", label);
      
      l = 8+strlen(label)+8+strlen(increment)+2;
      plus = (char *)calloc(sizeof(char), l+1);
      minus = (char *)calloc(sizeof(char), l+1+1);
      sprintf(plus, "label: %s\nvalue: %s\n.\n", label, increment);
      sprintf(minus, "label: %s\nvalue: -%s\n.\n", label, increment);
   }

   set_input_mode();
   
   if (verbose) {
      fprintf(stdout, "Connecting to host %s:%d\n", host, port);
      if (user) {
	 if (cred) {	    
	    fprintf(stdout, "Accounting as \"%s\" (with creds)\n", user);
	 } else {
	    fprintf(stdout, "Accounting as \"%s\"\n", user);
	 }
      }
      fprintf(stdout, "Incrementing drive %s by %s\n", label, increment);
      keys(stdout);
   }
   
   while (1) {
      size_t i;
      
      i = read(STDIN_FILENO, &c, 1);
      
      if (i > 0) {
	 if (c == '\004') {         /* `C-d' */
	    break;
	 } else {
	    switch (c) {
	       
	     case 'i':
	       if (send_message(host, port, inc) == -1) {
		  fprintf(stderr, "unable to connect to host %s:%d\n", 
			  host, port);
	       }
	       break;
	       
	     case 'p':
	       if (send_message(host, port, plus) == -1) {
		  fprintf(stderr, "unable to connect to host %s:%d\n", 
			  host, port);
	       }
	       break;
	       
	     case 'm':
	       if (send_message(host, port, minus) == -1) {
		  fprintf(stderr, "unable to connect to host %s:%d\n", 
			  host, port);
	       }
	       break;
	    }
	 }
      }
   }
   
   exit (EXIT_SUCCESS);
}
