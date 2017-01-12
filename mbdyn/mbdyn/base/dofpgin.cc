/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include <cstring>
#include <sstream>

#include "dofpgin.h"

DofPlugIn::DofPlugIn(MathParser& mp, DataManager *pDM)
: MathParser::PlugIn(mp), pDM(pDM), bPrev(false)
{
	ASSERT(pDM != NULL);
}
	
DofPlugIn::~DofPlugIn(void) 
{
	NO_OP;
}

const char *
DofPlugIn::sName(void) const 
{
	return "dof";
}

int 
DofPlugIn::Read(int argc, char *argv[])
{
	unsigned int uLabel;
	Node *pNode;
	int iOrder = 0;
	int iParams = 2;

	if (argc < 1 || argv[0] == NULL) {
		silent_cerr("DofPlugIn::Read(): "
			"illegal number of parameters #" << argc
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	uLabel = ReadLabel(argv[0]);

	if (argc < 2 || argv[1] == NULL) {
		silent_cerr("DofPlugIn::Read(" << argv[0] << "): "
			"illegal number of parameters " << argc
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pNode = ReadNode(uLabel, argv[1]);

	unsigned int iMaxIndex = pNode->iGetNumDof();
	switch (iMaxIndex) {
	case 0:
		/* parameter? */
		ASSERT(pNode->GetNodeType() == Node::PARAMETER);
		break;

	case 1:
		iParams++;

		if (argc < 3 || argv[2] == NULL) {
			silent_cerr("DofPlugIn::Read(" << argv[0] << ","
				<< argv[1] << "): "
				"illegal number of parameters " << argc 
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		iOrder = ReadDofOrder(pNode, 1, argv[2]);
		break;

	default: {
		iParams += 2;

		if (argc < 3 || argv[2] == NULL) {
			silent_cerr("DofPlugIn::Read(" << argv[0] << ","
				<< argv[1] << "): "
				"illegal number of parameters " << argc 
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		unsigned int iIndex = ReadIndex(pNode, iMaxIndex, argv[2]);

		if (argc < 4 || argv[3] == NULL) {
			silent_cerr("DofPlugIn::Read(" << argv[0] << ","
				<< argv[1] << ","
				<< argv[2] << "): "
				"illegal number of parameters " << argc 
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		iOrder = ReadDofOrder(pNode, iIndex, argv[3]);

		NodeDof nd(iIndex-1, pNode);
		pNode = 0;
		/* Chi dealloca questa memoria? ci vorrebbe l'handle */
		SAFENEWWITHCONSTRUCTOR(pNode, Node2Scalar, Node2Scalar(nd));
		pedantic_cerr(psNodeNames[pNode->GetNodeType()]
			<< "(" << pNode->GetLabel() << "): "
			"possibly allocating a NodeDof "
		 	"that nobody will delete" << std::endl);
		break;
	}
	}
	dof = ScalarDof(dynamic_cast<ScalarNode *>(pNode), iOrder, 0);

	/* legge i parametri (<param>=<value>, separati da '&') */
	if (argc == iParams + 1) {
		char *parm = NULL, *p, *end;

		SAFESTRDUP(parm, argv[iParams]);

		for (p = parm; p != NULL; p = end) {
			end = std::strchr(p, '&');
			while (end != NULL && end[-1] == '\\') {
				end = std::strchr(end + 1, '&');
			}

			if (end != NULL) {
				end[0] = '\0';
				end++;
			}

			char *v = std::strchr(p, '=');
			if (v == NULL) {
				silent_cerr("dof plugin: missing \"=\" "
					"in <param>=<value>" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			v[0] = '\0';
			v++;

			/* prende il valore al passo precedente */
			if (strncasecmp(p, "prev", STRLENOF("prev")) == 0) {
				if (strcasecmp(v, "true") == 0) {
					bPrev = true;
				} else if (strcasecmp(v, "false") == 0) {
					bPrev = false;
				} else {
					silent_cerr("dof plugin: "
						"unknown mode for parameter "
						"\"" << p << "=" << v << "\""
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else {
				silent_cerr("dof plugin: unknown parameter "
					"\"" << p << "=" << v << "\""
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		SAFEDELETEARR(parm);
	}

	return 0;
}

TypedValue::Type
DofPlugIn::GetType(void) const 
{
	return TypedValue::VAR_REAL;
}

TypedValue 
DofPlugIn::GetVal(void) const 
{
	if (bPrev) {
		return TypedValue(dof.dGetValuePrev());
	}
	return TypedValue(dof.dGetValue());
}

unsigned int 
DofPlugIn::ReadLabel(const char* s) 
{
	// must be semicolon-terminated to be read by math parser
	std::istringstream in(std::string(s) + ";");
	InputStream In(in);
	int i = mp.Get(In);
	if (i < 0) {
		silent_cerr("DofPlugIn::ReadLabel(" << s << "): invalid negative label" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return unsigned(i);
}

Node *
DofPlugIn::ReadNode(unsigned int uLabel, const char *s) 
{
	Node *pNode = NULL;
	unsigned int i;
	for (i = 0; i < Node::LASTNODETYPE; i++) {
		if (strcasecmp(s, psReadNodesNodes[i]) == 0) {
			break;
		}
	}
	if (i == Node::LASTNODETYPE) {
		silent_cerr("unknown node type '" << s << "'" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	if ((pNode = pDM->pFindNode(Node::Type(i), uLabel)) == NULL) {
		silent_cerr(psNodeNames[Node::Type(i)] 
			<< "(" << uLabel << ") not defined" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	return pNode;
}

unsigned int 
DofPlugIn::ReadIndex(Node *pNode, unsigned int iMaxIndex, const char *s) 
{
	unsigned int i = ReadLabel(s);
	if (i == 0 || i > iMaxIndex) {
		silent_cerr("illegal index " << i << " for "
			<< psNodeNames[pNode->GetNodeType()]
			<< "(" << pNode->GetLabel() << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	return i;
}

int 
DofPlugIn::ReadDofOrder(Node *pNode, unsigned int iIndex, const char *s) 
{
	if (strcasecmp(s, "differential") == 0) {
		if (pNode->GetDofType(iIndex-1) != DofOrder::DIFFERENTIAL) {
			silent_cerr("cannot take differential value of "
				<< psNodeNames[pNode->GetNodeType()] 
				<< "(" << pNode->GetLabel() << ")"
				"[" << iIndex << "]" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		return 1;
	}
	
	if (strcasecmp(s, "algebraic") == 0) {
		return 0;
	} 

	silent_cerr("unknown dof order '" << s << "'" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

MathParser::PlugIn *
dof_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = NULL;
	SAFENEWWITHCONSTRUCTOR(p, DofPlugIn, DofPlugIn(mp, (DataManager *)arg));
	return p;
}

