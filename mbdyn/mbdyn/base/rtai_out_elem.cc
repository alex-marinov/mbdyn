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

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <rtai_out_elem.h>
#include <dataman.h>

/* RTAIOutElem - begin */

RTAIOutElem::RTAIOutElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		const char *h, unsigned long n)
: Elem(uL, Elem::RTAI_OUTPUT, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(NULL),
host(h), node(n), port(-1), mbx(NULL)
{
	/* FIXME: size depends on the type of the output signals */
	size = sizeof(double)*nch;
	SAFENEWARR(buf, char, size);

	SAFENEW(mbx, MBX);
	if (rt_mbx_init(mbx, size)) {
		std::cerr << "RTAI mailbox failed" << std::endl;
		THROW(ErrGeneric());
	}

	if (node) {
		/* get port ... */
	}
}

RTAIOutElem::~RTAIOutElem(void)
{
	if (mbx) {
		rt_mbx_delete(mbx);
		SAFEDELETE(mbx);
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
RTAIOutElem::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

Elem::Type
RTAIOutElem::GetElemType(void) const
{
	return Elem::RTAI_OUTPUT;
}

void
RTAIOutElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
RTAIOutElem::AssRes(SubVectorHandler& WorkVec, double dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
RTAIOutElem::AssJac(VariableSubMatrixHandler& WorkMat, double dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

void
RTAIOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	for (unsigned int i; i < NumChannels; i++) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = pNodes[i].pNode->dGetDofValue(1, pNodes[i].iOrder);
		memcpy(&buf[sizeof(double)*i], &v, sizeof(double));
	}

	if (RT_mbx_send_if(node, port, mbx, buf, size) != size) {
		/* error */
	}
}


Elem *
ReadRTAIOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	unsigned long node = 0;
	const char *host = NULL;

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

	int nch = HP.GetInt();
	if (nch <= 0) {
		std::cerr << "illegal number of channels for "
			<< psElemNames[Elem::RTAI_OUTPUT]
			<< "(" << uLabel << ") at line " << HP.GetLineData()
			<< std::endl;
		THROW(ErrGeneric());
	}

	ScalarDof *pNodes = NULL;
	SAFENEWARR(pNodes, ScalarDof, nch);
	for (int i = 0; i < nch; i++) {
		pNodes[i] = ReadScalarDof(pDM, HP, 0);
	}

	Elem *pEl = NULL;
	SAFENEWWITHCONSTRUCTOR(pEl, RTAIOutElem,
			RTAIOutElem(uLabel, nch, pNodes, host, node));
	return pEl;
}

#endif /* USE_RTAI */

