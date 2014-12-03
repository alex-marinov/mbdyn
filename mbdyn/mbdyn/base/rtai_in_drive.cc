/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <netdb.h>

#include "dataman.h"
#include "rtai_in_drive.h"
#include "mbrtai_utils.h"

RTMBDynInDrive::RTMBDynInDrive(unsigned int uL,
	const DriveHandler* pDH,
	const std::string& sFileName,
	const std::string& host,
	integer nd, const std::vector<doublereal>& v0,
	bool c, unsigned long /*int*/ n,
	bool bNonBlocking)
: StreamDrive(uL, pDH, sFileName, nd, v0, c),
host(host), node(n), port(-1), bNonBlocking(bNonBlocking),
mbx(NULL)
{
	ASSERT(!sFileName.empty());

	if (create) {
		ASSERT(node == 0);

		if (rtmbdyn_rt_mbx_init(sFileName.c_str(), size, &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << sFileName << ") "
				"init failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (node) {
			/* get port ... */
			port = rtmbdyn_rt_request_port(node);
			/* FIXME: what in case of failure? */
		}

		if (rtmbdyn_RT_get_adr(node, port, sFileName.c_str(), &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << sFileName << ") "
				"get_adr failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (bNonBlocking) {
		f_receive = rtmbdyn_RT_mbx_receive_if;

	} else {
		f_receive = rtmbdyn_RT_mbx_receive;
	}
}

RTMBDynInDrive::~RTMBDynInDrive(void) 
{
	/*
	 * destroy mailbox and so on
	 */
	if (mbx) {
		rtmbdyn_rt_mbx_delete(&mbx);
	}
}

void 
RTMBDynInDrive::ServePending(const doublereal& /* t */ )
{
	/*
	 * store in pdVal the values of all the channels
	 * served by the mailbox
	 */
	int rc = f_receive(node, port, mbx, (void *)buf, size);
	if (!rc) {
		doublereal *rbuf = (doublereal *)buf;
		for (int i = 1; i <= iNumDrives; i++) {
			pdVal[i] = rbuf[i-1];
		}	

	} else {
		/* FIXME: error */
	}
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
RTMBDynInDrive::Restart(std::ostream& out) const
{
   	return out << "0. /* RTMBDynInDrive not implemented yet */" << std::endl;
}

Drive *
ReadRTMBDynInDrive(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	unsigned long node = (unsigned long)-1;
	std::string host;
	std::string name;
	bool create = false;

	if (HP.IsKeyWord("stream" "drive" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == NULL) {
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"unable to read mailbox name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (strlen(m) != 6) {
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"illegal mailbox name \"" << m << "\" "
				"(must be exactly 6 chars) "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		name = m;

	} else {
		silent_cerr("RTMBDynInDrive(" << uLabel << "): "
			"missing mailbox name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("create")) {
		if (!HP.GetYesNo(create)) {
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"\"create\" must be \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("non" "blocking")) {
		silent_cout("RTMBDynInDrive(" << uLabel << "): "
			"RTMBDyn mailboxes are always non-blocking"
			<< std::endl);
	}

	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetStringWithDelims();
		
		silent_cout("RTMBDynInDrive(" << uLabel << "): "
			"local path \"" << m << "\" silently ignored"
			<< std::endl);
	}

	if (HP.IsKeyWord("port")){
		int p = HP.GetInt();
		
		silent_cout("RTMBDynInDrive(" << uLabel << "): "
			"port " << p << " silently ignored" << std::endl);
	}
	
	if (HP.IsKeyWord("host")) {
		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == NULL) {
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (create) {
			silent_cout("RTMBDynInDrive(" << uLabel << "): "
				"host name \"" << h << "\" silently ignored"
				<< std::endl);
		} else {
			host = h;

			// resolve host
			// TODO: add support for getnameinfo()
#if defined(HAVE_GETADDRINFO)
			struct addrinfo hints = { 0 }, *res = NULL;
			int rc;

			hints.ai_family = AF_INET;
			hints.ai_socktype = SOCK_STREAM; // FIXME: SOCK_DGRAM?
			rc = getaddrinfo(host.c_str(), NULL, &hints, &res);
			if (rc == 0) {
				node = ((struct sockaddr_in *)res->ai_addr)->sin_addr.s_addr;
				freeaddrinfo(res);
			}
#elif defined(HAVE_GETHOSTBYNAME)
			// FIXME: non-reentrant
			struct hostent *he = gethostbyname(host.c_str());
			if (he != NULL)
			{
				node = ((unsigned long *)he->h_addr_list[0])[0];
			} 
#elif defined(HAVE_INET_ATON)
			struct in_addr addr;
			if (inet_aton(host.c_str(), &addr)) {
				node = addr.s_addr;
			}
#else // ! HAVE_GETADDRINFO && ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"host (RTAI RPC) not supported "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! HAVE_GETADDRINFO && ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON

			if (node == (unsigned long)-1) {
				silent_cerr("RTMBDynInDrive(" << uLabel << "): "
					"unable to convert host \"" << host << "\" to node" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	bool bNonBlocking(true);
	while (HP.IsArg()) {
		if (HP.IsKeyWord("signal")) {
			// ignore

		} else if (HP.IsKeyWord("no" "signal")) {
			// ignore

		} else if (HP.IsKeyWord("blocking")) {
			bNonBlocking = false;

		} else if (HP.IsKeyWord("non" "blocking")) {
			bNonBlocking = true;

		} else {
			break;
		}
	}

	int idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("RTMBDynInDrive(" << uLabel << "): "
			"illegal number of channels "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<doublereal> v0;
	if (HP.IsKeyWord("initial" "values")) {
		v0.resize(idrives);
		for (int i = 0; i < idrives; i++) {
			v0[i] = HP.GetReal();
		}
	}

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr, RTMBDynInDrive,
			RTMBDynInDrive(uLabel, 
			pDM->pGetDrvHdl(),
			name, host, idrives, v0, create, node, bNonBlocking));
	
	return pDr;
}
 
