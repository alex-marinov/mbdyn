/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* socket driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

/* include del programma */

//#include <rtai.h>
//#include <rtai_sched.h>
#include <net_rpc.h>

#include <mbrtai_utils.h>

int
mbdyn_rt_request_port(unsigned long node)
{
	return rt_request_port(node);
}

int
mbdyn_rt_mbx_init(void **__mbx, int size)
{
	MBX	*mbx;
	int	rc;

	mbx = (MBX *)malloc(sizeof(MBX));
	if (mbx == NULL) {
		return -1;
	}
	
	rc = rt_mbx_init(mbx, size);
	if (rc) {
		free(mbx);
		mbx = NULL;
	}

	*__mbx = mbx;

	return rc;
}

int
mbdyn_rt_mbx_delete(void **__mbx)
{
	MBX	*mbx = (MBX *)__mbx;
	int	rc;

	rc = rt_mbx_delete(mbx);

	free(mbx);

	*__mbx = NULL;

	return rc;
}

int
mbdyn_RT_mbx_send_if(unsigned long node, int port, void *__mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)__mbx;
	
	return RT_mbx_send_if(node, port, mbx, msg, msg_size);
}

#endif /* USE_RTAI */

