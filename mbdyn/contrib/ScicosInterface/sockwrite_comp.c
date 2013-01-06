/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 * Author: Tommaso Solcia <tommaso.solcia@mail.polimi.it>
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <machine.h>
#include <math.h>

#include "scilab/scicos_blocks/scicos_block4.h"

/* inputs */
#define in(i)         (GetRealInPortPtrs(block, i+1))[0]

/* outputs */
#define out(i)        (GetRealOutPortPtrs(block, i+1))[0]

/* workspace */
#define workspace     GetWorkPtrs(block)

/* parameters */
#define IparPtrs      GetIparPtrs(block)
#define hostname      Getint8OparPtrs(block, 1)

#define PORT          IparPtrs[0]
#define N_CH          IparPtrs[1]

void
sockwrite(scicos_block *block, int flag)
{
	struct sockaddr_in server_addr;
	struct hostent *host;
	int sockfd, numbytes;
	/* warning: C99 only */
	double num[N_CH];

	if (flag == 4) {
		/* Initializing block */
		int *p_sock;
		p_sock = scicos_malloc(sizeof(int));
		if (p_sock == NULL){
			set_block_error(-16);
		}
		workspace = p_sock;

		/* create a socket: */
		sockfd = socket(PF_INET, SOCK_STREAM, 0);
		if (sockfd == -1) {
			perror("socket in function sockwrite");
			set_block_error(-3);
		}
		*p_sock = sockfd;

		/* get the host by name */
		host = gethostbyname(hostname);

		/* make the address structure */
		server_addr.sin_family = AF_INET;
		server_addr.sin_port = htons(PORT);
		server_addr.sin_addr = *((struct in_addr *)host->h_addr);
		bzero(&(server_addr.sin_zero), 8);

		/* connect it to the server through the port we passed
		 * in serv_addr.sin_port: */
		if (connect(sockfd, (struct sockaddr *)&server_addr,
			sizeof(struct sockaddr)) == -1)
		{
			perror("connect in function sockwrite");
			set_block_error(-3);
		}

	} else if (flag == 1) {
		 /* Updating output */
		int *p_sock = workspace;
		int k;

		sockfd = *p_sock;

		for (k = 0; k < N_CH; k++) {
			num[k] = in(k);
			out(k) = num[k];
		}

		numbytes = send(sockfd, (void *)&num, sizeof(num),
			MSG_NOSIGNAL);
		if (numbytes == -1) {
			perror("send in function sockwrite");
			set_block_error(-3);
		}

	} else if (flag == 5) {
		/* Ending */
		int *p_sock = workspace;

		sockfd = *p_sock;
		close(sockfd);

		scicos_free(workspace);
	}
}
