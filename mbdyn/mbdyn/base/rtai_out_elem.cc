/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/* include del programma */

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <rtai_out_elem.h>
#include <mbrtai_utils.h>
#include <dataman.h>

/* RTMBDynOutElem - begin */

RTMBDynOutElem::RTMBDynOutElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		const char *h, const char *m, unsigned long n, bool c)
: Elem(uL, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(0),
host(h), node(n), name(m), create(c), port(-1), mbx(0)
{
	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*nch;
	SAFENEWARR(buf, char, size);

	/* RATIONALE:
	 *
	 * if host/node is present, the mailbox is "remote";
	 * if it not, we may need to create it
	 */

	if (create) {
		ASSERT(node == 0);

		if (mbdyn_rt_mbx_init(name, size, &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << name << ") "
				"init failed" << std::endl);
			throw ErrGeneric();
		}

	} else {
		if (node) {
			/* get port ... */
			port = mbdyn_rt_request_port(node);
			/* FIXME: what in case of failure? */
		}

		if (mbdyn_RT_get_adr(node, port, name, &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << name << ") "
				"get_adr failed" << std::endl);
			throw ErrGeneric();
		}
	}
}

RTMBDynOutElem::~RTMBDynOutElem(void)
{
	if (mbx) {
		mbdyn_rt_mbx_delete(&mbx);
	}

	if (buf) {
		SAFEDELETEARR(buf);
	}

	if (pNodes) {
		/* FIXME: pNodes leak */
		SAFEDELETEARR(pNodes);
	}

	if (host) {
		SAFEDELETEARR(host);
	}
}

std::ostream&
RTMBDynOutElem::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

Elem::Type
RTMBDynOutElem::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
RTMBDynOutElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
RTMBDynOutElem::AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
RTMBDynOutElem::AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

void
RTMBDynOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	char *curbuf = buf;
	
	for (unsigned int i = 0; i < NumChannels; i++) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = pNodes[i].dGetValue();

		doublereal *dbuf = (doublereal *)curbuf;
		dbuf[0] = v;

		curbuf += sizeof(doublereal);
	}
	if (mbdyn_RT_mbx_send_if(node, -port, mbx, (void *)buf, size) != size) {
		/* error */
	}
	 
}

/* Inverse Dynamics */
void
RTMBDynOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);	
}

Elem *
ReadRTMBDynOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	unsigned long node = 0;
	const char *host = NULL;
	const char *name = NULL;
	bool create = false;

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == NULL) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"unable to read mailbox name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();

		} else if (strlen(m) != 6) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"illegal mailbox name \"" << m << "\" "
				"(must be exactly 6 chars) "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		SAFESTRDUP(name, m);

	} else {
		silent_cerr("RTMBDynOutElem(" << uLabel << "): "
			"missing mailbox name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			create = true;
		} else if (HP.IsKeyWord("no")) {
			create = false;
		} else {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"\"create\" must be \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetStringWithDelims();
		silent_cout("RTMBDynOutElem(" << uLabel << "): "
			"local path \"" << m << "\" silently ignored "
			"at line " << HP.GetLineData() << std::endl);

	}

	if (HP.IsKeyWord("port")){
		int p = HP.GetInt();
		
		silent_cout ("RTMBDynOutElem(" << uLabel << "): "
			"port " << p << " silently ignored "
			"at line " << HP.GetLineData() << std::endl);
	}

	if (HP.IsKeyWord("host")) {
#if defined(HAVE_GETHOSTBYNAME) || defined(HAVE_INET_ATON)
		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == NULL) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		if (create) {
			silent_cout("RTMBDynOutElem(" << uLabel << "): "
				"host name \"" << h << "\" silently ignored "
				"at line " << HP.GetLineData() << std::endl);			
		} else {

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
				silent_cerr("RTMBDynOutElem(" << uLabel << "): "
					"unable to convert host "
					"\"" << host << "\" "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}
#else /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"host (RTAI RPC) not supported "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
#endif /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
		}
	}

	int nch = HP.GetInt();
	if (nch <= 0) {
		silent_cerr("RTMBDynOutElem(" << uLabel << "): "
			"illegal number of channels "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	ScalarDof *pNodes = NULL;
	SAFENEWARRNOFILL(pNodes, ScalarDof, nch);
	for (int i = 0; i < nch; i++) {
		pNodes[i] = ReadScalarDof(pDM, HP, 1);
	}

   	// (void)pDM->fReadOutput(HP, Elem::LOADABLE); 

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("RTMBDynOutElem(" << uLabel << "): "
			"semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}
      
	/* costruzione dell'elemento */
	Elem *pEl = NULL;
	SAFENEWWITHCONSTRUCTOR(pEl, RTMBDynOutElem,
			RTMBDynOutElem(uLabel, nch, pNodes,
				host, name, node, create));
	return pEl;
}

