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

#include <sstream.h>
#include <dofpgin.h>

DofPlugIn::DofPlugIn(MathParser& mp, DataManager *pDM)
: MathParser::PlugIn(mp), pDM(pDM) 
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
	unsigned int uLabel = ReadLabel(argv[0]);
	Node *pNode = ReadNode(uLabel, argv[1]);
	int iOrder = 0;
		
	unsigned int iMaxIndex = pNode->iGetNumDof();
	if (iMaxIndex > 1) {
		unsigned int iIndex = ReadIndex(pNode, iMaxIndex, argv[2]);
		iOrder = ReadDofOrder(pNode, iIndex, argv[3]);
		NodeDof nd(pNode->GetLabel(), iIndex-1, pNode);
		pNode = NULL;
		/* Chi dealloca questa memoria? ci vorrebbe l'handle */
		SAFENEWWITHCONSTRUCTOR(pNode, 
				       Node2Scalar, 
				       Node2Scalar(nd), 
				       DMmm);
		cerr << "warning, possibly allocating a NodeDof that nobody"
			" will delete until handles will be used" << endl;
	} else {
		iOrder = ReadDofOrder(pNode, 1, argv[2]);
	}
	dof = ScalarDof((ScalarNode*)pNode, iOrder);

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
	return TypedValue(dof.dGetValue());
}

unsigned int 
DofPlugIn::ReadLabel(const char* s) 
{
	unsigned int rc;
	char *stmp = NULL;

	/*
	 * deve essere terminato da ';' per essere letto da math parser :(
	 */
	SAFENEWARR(stmp, char, strlen(s)+2, DMmm);
	strcpy(stmp, s);
	strcat(stmp, ";");
#if defined(HAVE_SSTREAM)
	std::istringstream in(stmp);
#else /* HAVE_STRSTREAM_H */
	istrstream in(stmp);
#endif /* HAVE_STRSTREAM_H */
	InputStream In(in);
	rc = (unsigned int)mp.Get(In);
	SAFEDELETEARR(stmp, DMmm);

	return rc;
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
		cerr << "unknown node type '" << s << "'" << endl;
		THROW(ErrGeneric());
	}
	if ((pNode = pDM->pFindNode(Node::Type(i), uLabel)) == NULL) {
		cerr << "node " << uLabel << " not defined" << endl;
		THROW(ErrGeneric());
	}
	return pNode;
}

unsigned int 
DofPlugIn::ReadIndex(Node *pNode, unsigned int iMaxIndex, const char *s) 
{
	unsigned int i = ReadLabel(s);
	if (i == 0 || i > iMaxIndex) {
		cerr << "illegal index " << i << " for node "
			<< psNodeNames[pNode->GetNodeType()]
			<< "(" << pNode->GetLabel() << ")" << endl;
		THROW(ErrGeneric());
	}
	return i;
}

int 
DofPlugIn::ReadDofOrder(Node *pNode, unsigned int iIndex, const char *s) 
{
	if (strcasecmp(s, "differential") == 0) {
		if (pNode->SetDof(iIndex-1) != DofOrder::DIFFERENTIAL) {
			cerr << "cannot take differential value of "
				<< psNodeNames[pNode->GetNodeType()] 
				<< "(" << pNode->GetLabel() << ")[" 
				<< iIndex << "]" << endl;
			THROW(ErrGeneric());
		}
		return 1;
	} else if (strcasecmp(s, "algebraic") == 0) {
		return 0;
	} else {
		cerr << "unknown dof order '" << s << "'" << endl;
		THROW(ErrGeneric());
	}
	return -1;
}

MathParser::PlugIn *
dof_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = NULL;
	SAFENEWWITHCONSTRUCTOR(p,
			       DofPlugIn,
			       DofPlugIn(mp, (DataManager *)arg),
			       DMmm);
	return p;
}

