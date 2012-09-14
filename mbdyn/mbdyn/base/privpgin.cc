/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include <sstream>
#include "privpgin.h"

PrivPlugIn::PrivPlugIn(MathParser& mp, DataManager *pDM)
: MathParser::PlugIn(mp), pSE(0), iIndex(0), pDM(pDM)
{
	ASSERT(pDM != 0);
}

PrivPlugIn::~PrivPlugIn(void) 
{
	NO_OP;
}

int 
PrivPlugIn::Read(int argc, char *argv[])
{
	unsigned int uLabel;

	if (argc < 1 || argv[0] == 0) {
		silent_cerr("PrivPlugIn::Read(): "
			"illegal number of parameters " << argc
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	uLabel = ReadLabel(argv[0]);

	if (argc < 2 || argv[1] == 0) {
		silent_cerr("PrivPlugIn::Read(" << argv[0] << "): "
			<< "illegal number of parameters " << argc
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	ReadSE(uLabel, argv[1]);

	unsigned int iMaxIndex = pSE->iGetNumPrivData();
	switch (iMaxIndex) {
	case 0:
		silent_cerr(*this << "allows no private data" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case 1:
		iIndex = 1;
		if (argc < 3) {
			break;
		}
		/* continue to next case */

	default:
		if (argc > 3) {
			silent_cerr("PrivPlugIn::Read(" << argv[0] << "): "
				<< "illegal number of parameters " << argc
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ReadIndex(iMaxIndex, argv[2]);
		break;
	}

	return 0;
}

TypedValue::Type
PrivPlugIn::GetType(void) const 
{
	return TypedValue::VAR_REAL;
}

TypedValue 
PrivPlugIn::GetVal(void) const 
{
	return TypedValue(pSE->dGetPrivData(iIndex));
}

unsigned int 
PrivPlugIn::ReadLabel(const char* s) 
{
	unsigned int rc;

	/*
	 * deve essere terminato da ';' per essere letto da math parser :(
	 */
	std::istringstream in(std::string(s) + ";");
	InputStream In(in);
	rc = (unsigned int)mp.Get(In);

	return rc;
}

void
PrivPlugIn::ReadIndex(unsigned int iMaxIndex, const char *s) 
{
	bool bRefresh(false);
	if (strncasecmp(s, "string=", STRLENOF("string=")) == 0) {
		sIndexName = &s[STRLENOF("string=")];
		iIndex = pSE->iGetPrivDataIdx(sIndexName.c_str());
		bRefresh = true;

	} else if (strncasecmp(s, "name=", STRLENOF("name=")) == 0) {
		sIndexName = &s[STRLENOF("name=")];
		iIndex = pSE->iGetPrivDataIdx(sIndexName.c_str());
		bRefresh = true;

		silent_cerr("PrivPlugIn: "
			"\"name=" << sIndexName << "\" is deprecated; "
			"use \"string=" << sIndexName << "\" instead" << std::endl);

	} else if (strncasecmp(s, "index=", STRLENOF("index=")) == 0) {
		s = &s[STRLENOF("index=")];
		iIndex = ReadLabel(s);

	} else {
		iIndex = pSE->iGetPrivDataIdx(s);
		if (iIndex > 0) {
			bRefresh = true;
		} else {
			iIndex = ReadLabel(s);
		}
	}

	if (bRefresh) {
		// refresh, as iGetPrivDataIdx could have changed it
		iMaxIndex = pSE->iGetNumPrivData();
	}

	if (iIndex == 0 || iIndex > iMaxIndex) {
		silent_cerr("illegal index " << iIndex << " for "
			<< *this << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

NodePrivPlugIn::NodePrivPlugIn(MathParser& mp, DataManager *pDM)
: PrivPlugIn(mp, pDM)
{
	NO_OP;
}

NodePrivPlugIn::~NodePrivPlugIn(void) 
{
	NO_OP;
}

const char *
NodePrivPlugIn::sName(void) const 
{
	return "node";
}

void
NodePrivPlugIn::ReadSE(unsigned int uLabel, const char *ss) 
{
	unsigned int i;
	char *s = 0;

	/* eat spaces */
	SAFESTRDUP(s, ss);
	for (i = 0; s[i]; i++) {
		if (isspace(s[i])) {
			memmove(&s[i], &s[i + 1], strlen(&s[i]));
		}
	}

	i = str2nodetype(s);
	
	SAFEDELETEARR(s);
	
	if (i == Node::LASTNODETYPE) {
		silent_cerr("unknown node type '" << ss << "'" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pSE = dynamic_cast<SimulationEntity *>(pDM->pFindNode(Node::Type(i), uLabel));
	if (pSE == 0) {
		silent_cerr("NodePrivPlugIn: "
			<< psReadNodesNodes[i] << "(" << uLabel << ") "
			"not defined" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

std::ostream&
NodePrivPlugIn::Err(std::ostream& out) const
{
	Node *pNode = dynamic_cast<Node *>(pSE);
	if (pNode == 0) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return out << psNodeNames[pNode->GetNodeType()]
		<< "(" << pNode->GetLabel() << ")";
}

ElemPrivPlugIn::ElemPrivPlugIn(MathParser& mp, DataManager *pDM)
: PrivPlugIn(mp, pDM)
{
	NO_OP;
}

ElemPrivPlugIn::~ElemPrivPlugIn(void) 
{
	NO_OP;
}

const char *
ElemPrivPlugIn::sName(void) const 
{
	return "element";
}

void
ElemPrivPlugIn::ReadSE(unsigned int uLabel, const char *ss) 
{
	unsigned int i;
	char *s = 0;

	/* eat spaces */
	SAFESTRDUP(s, ss);
	for (i = 0; s[i]; i++) {
		if (isspace(s[i])) {
			memmove(&s[i], &s[i + 1], strlen(&s[i]));
		}
	}

	i = str2elemtype(s);

	SAFEDELETEARR(s);
	
	if (i == Elem::LASTELEMTYPE) {
		silent_cerr("unknown element type \"" << ss << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pSE = dynamic_cast<SimulationEntity *>(pDM->pFindElem(Elem::Type(i), uLabel));
	if (pSE == 0) {
		silent_cerr("ElemPrivPlugIn: "
			<< psReadElemsElems[i] << "(" << uLabel << ") "
			"not defined" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

std::ostream&
ElemPrivPlugIn::Err(std::ostream& out) const
{
	Elem *pElem = dynamic_cast<Elem *>(pSE);
	if (pElem == 0) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return out << psElemNames[pElem->GetElemType()]
		<< "(" << pElem->GetLabel() << ")";
}

std::ostream& 
operator << (std::ostream& out, const PrivPlugIn& p)
{
	return out << p.Err(out);
}

MathParser::PlugIn *
node_priv_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = 0;
	SAFENEWWITHCONSTRUCTOR(p, NodePrivPlugIn,
			NodePrivPlugIn(mp, (DataManager *)arg));
	return p;
}

MathParser::PlugIn *
elem_priv_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = 0;
	SAFENEWWITHCONSTRUCTOR(p, ElemPrivPlugIn,
			ElemPrivPlugIn(mp, (DataManager *)arg));
	return p;
}

