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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

#include <netdb.h>

#include <dataman.h>
#include <rtai_in_drive.h>
#include <mbrtai_utils.h>

RTAIInDrive::RTAIInDrive(unsigned int uL,
			 const DriveHandler* pDH,
			 const char* const sFileName,
			 const char *h,
			 integer nd ,bool c, int n)
: FileDrive(uL, pDH, sFileName, nd),//host 
host(h), create(c), node(n), port(-1), mbx(NULL)
{
   	ASSERT(nd > 0);
	ASSERT(sFileName != NULL);
	
	/*
	 * initialize mailbox and so on
	 */
	size = sizeof(double)*nd;
	SAFENEWARR(buf, char, size);

	if (create) {
		ASSERT(node == 0);

		if (mbdyn_rt_mbx_init(sFileName, size, &mbx)) {
			std::cerr << "RTAI mailbox(" << sFileName << ") "
				"init failed" << std::endl;
			THROW(ErrGeneric());
		}

	} else {
		if (node) {
			/* get port ... */
			port = mbdyn_rt_request_port(node);
			/* FIXME: what in case of failure? */
		}

		if (mbdyn_RT_get_adr(node, port, sFileName, &mbx)) {
			std::cerr << "RTAI mailbox(" << sFileName << ") "
				"get_adr failed" << std::endl;
			THROW(ErrGeneric());
		}
	}
	
}

RTAIInDrive::~RTAIInDrive(void) 
{
	/*
	 * destroy mailbox and so on
	 */
	if (mbx) {
		mbdyn_rt_mbx_delete(&mbx);
	}
	
	if (host) {
		SAFEDELETEARR(host);
	}
	NO_OP;
}

void 
RTAIInDrive::ServePending(const doublereal& /* t */ )
{
	/*
	 * store in pdVal the values of all the channels
	 * served by the mailbox
	 */
	if(!(mbdyn_RT_mbx_receive_if(node, port, mbx, (void *)buf, size))){
		
		doublereal *rbuf = (doublereal *)buf;
		for (int i = 1; i <= iNumDrives; i++){
			pdVal[i] = rbuf[i-1];
		}	
	}else{
			/*FIXME: error */
	}
	
}

FileDrive::Type 
RTAIInDrive::GetFileDriveType(void) const
{
   	return FileDrive::RTAI_IN;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
RTAIInDrive::Restart(std::ostream& out) const
{
   	return out << "# RTAIInDrive not implemented yet" << std::endl;
}

Drive *
ReadRTAIInDrive(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	unsigned long node = 0;
	const char *host = NULL;
	const char *name = NULL;
	bool create = false;

	if (HP.IsKeyWord("mailbox" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == NULL) {
			std::cerr << "unable to read mailbox name "
				"for RTAIInElem(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl;
			THROW(ErrGeneric());

		} else if (strlen(m) != 6) {
			std::cerr << "illegal mailbox name \"" << m
				<< "\" for RTAIInElem(" << uLabel 
				<< ") (must be 6 char) at line "
				<< HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}

		SAFESTRDUP(name, m);

	} else {
		std::cerr << "missing mailbox name for RTAIInElem(" << uLabel
			<< ") at line " << HP.GetLineData() << std::endl;
		THROW(ErrGeneric());
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			create = true;
		} else if (HP.IsKeyWord("no")) {
			create = false;
		} else {
			std::cerr << "\"create\" must be \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}
	}

	if (HP.IsKeyWord("host")) {
#if defined(HAVE_GETHOSTBYNAME) || defined(HAVE_INET_ATON)
		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == NULL) {
			std::cerr << "unable to read host for "
				<< psElemNames[Elem::RTAI_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}

		if (create) {
			std::cerr << "cannot create mailbox(" << name
				<< ") as remote on host " << h << std::endl;
			THROW(ErrGeneric());
		}

		SAFESTRDUP(host, h);

		/* resolve host
		 * FIXME: non-reentrant ... */
#if defined(HAVE_GETHOSTBYNAME)
		struct hostent *he = gethostbyname(host);
		if (he != NULL)
		{
			node = ((unsigned long *)he->h_addr_list[0])[0];
		} 
#elif defined(HAVE_INET_ATON)
		struct in_addr addr;
		if (inet_aton(host, &addr)) {
			node = addr.s_addr;
		}
#endif /* ! HAVE_GETHOSTBYNAME && HAVE_INET_ATON */
		else {
			std::cerr << "unable to convert host "
				<< host << " at line "
				<< HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}
#else /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
		std::cerr << "host (RTAI RPC) not supported for "
			<< psElemNames[Elem::RTAI_OUTPUT]
			<< "(" << uLabel << ") at line " << HP.GetLineData()
			<< std::endl;
		THROW(ErrGeneric());
#endif /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
	}

	int idrives = HP.GetInt();
	if (idrives <= 0) {
		std::cerr << "illegal number of channels for "
			<< psElemNames[Elem::RTAI_OUTPUT]
			<< "(" << uLabel << ") at line " << HP.GetLineData()
			<< std::endl;
		THROW(ErrGeneric());
	}

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr, RTAIInDrive,
			RTAIInDrive(uLabel, 
			pDM->pGetDrvHdl(),
			name, host, idrives, create, node));
	
	return pDr;
}
   
#endif /* USE_RTAI */

