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

#include <ac/sstream>
#include <privpgin.h>

PrivPlugIn::PrivPlugIn(MathParser& mp, DataManager *pDM)
: MathParser::PlugIn(mp), pElem(0), iIndex(0), sIndexName(0), pDM(pDM)
{
	ASSERT(pDM != NULL);
}
	
PrivPlugIn::~PrivPlugIn(void) 
{
	if (sIndexName) {
		SAFEDELETEARR(sIndexName);
	}
}

const char *
PrivPlugIn::sName(void) const 
{
	return "element";
}

int 
PrivPlugIn::Read(int argc, char *argv[])
{
	unsigned int uLabel;

	if (argc < 1 || argv[0] == NULL) {
		std::cerr << "PrivPlugIn::Read(): "
			"illegal number of parameters " << argc
			<< std::endl;
		THROW(ErrGeneric());
	}
	uLabel = ReadLabel(argv[0]);

	if (argc < 2 || argv[1] == NULL) {
		std::cerr << "PrivPlugIn::Read(" << argv[0] 
			<< "): illegal number of parameters " << argc
			<< std::endl;
		THROW(ErrGeneric());
	}
	ReadElem(uLabel, argv[1]);

	unsigned int iMaxIndex = pElem->iGetNumPrivData();
	switch (iMaxIndex) {
	case 0:
		silent_cerr(psElemNames[pElem->GetElemType()]
				<< "(" << pElem->GetLabel() << ") "
				"allows no private data" << std::endl);
		THROW(ErrGeneric());

	case 1:
		iIndex = 1;
		if (argc < 3) {
			break;
		}
		/* continue to next case */

	default:
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
	return TypedValue(pElem->dGetPrivData(iIndex));
}

unsigned int 
PrivPlugIn::ReadLabel(const char* s) 
{
	unsigned int rc;
	char *stmp = NULL;

	/*
	 * deve essere terminato da ';' per essere letto da math parser :(
	 */
	int l = strlen(s)+2;
	SAFENEWARR(stmp, char, l);
	strcpy(stmp, s);
	strcat(stmp, ";");
#if defined(HAVE_SSTREAM)
	std::istringstream in(stmp);
#else /* HAVE_STRSTREAM_H */
	istrstream in(stmp);
#endif /* HAVE_STRSTREAM_H */
	InputStream In(in);
	rc = (unsigned int)mp.Get(In);
	SAFEDELETEARR(stmp);

	return rc;
}

void
PrivPlugIn::ReadElem(unsigned int uLabel, const char *ss) 
{
	unsigned int i;
	char *s = NULL;

	/* eat spaces */
	SAFESTRDUP(s, ss);
	for (i = 0; s[i]; i++) {
		if (isspace(s[i])) {
			memmove(&s[i], &s[i + 1], strlen(&s[i]));
		}
	}
	
	for (i = 0; i < Elem::LASTELEMTYPE; i++) {
		if (strcasecmp(s, psReadElemsElems[i]) == 0) {
			break;
		}
	}

	SAFEDELETEARR(s);
	
	if (i == Elem::LASTELEMTYPE) {
		std::cerr << "unknown element type '" << ss << "'" << std::endl;
		THROW(ErrGeneric());
	}

	if ((pElem = (Elem *)pDM->pFindElem(Elem::Type(i), uLabel)) == NULL) {
		std::cerr << psElemNames[Elem::Type(i)] 
			<< "(" << uLabel << ") not defined" << std::endl;
		THROW(ErrGeneric());
	}
}

void
PrivPlugIn::ReadIndex(unsigned int iMaxIndex, const char *s) 
{
	bool bIsName(false);

	if (strncasecmp(s, "name=", sizeof("name=") - 1) == 0) {
		iIndex = pElem->iGetPrivDataIdx(s + sizeof("name=") - 1);
		bIsName = true;

	} else {
		iIndex = ReadLabel(s);
	}

	if (iIndex == 0 || iIndex > iMaxIndex) {
		silent_cerr("illegal index " << iIndex << " for "
			<< psElemNames[pElem->GetElemType()]
			<< "(" << pElem->GetLabel() << ")"
			<< std::endl);
		THROW(ErrGeneric());
	}

	if (bIsName) {
		SAFESTRDUP(sIndexName, s + sizeof("name=") - 1);
	}
}

MathParser::PlugIn *
priv_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = NULL;
	SAFENEWWITHCONSTRUCTOR(p, PrivPlugIn,
			PrivPlugIn(mp, (DataManager *)arg));
	return p;
}

