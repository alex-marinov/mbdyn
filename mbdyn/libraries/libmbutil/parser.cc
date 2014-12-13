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

/* parser */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <stdlib.h>
#include <stack>
#include <map>

#include "mathtyp.h"
#include "parser.h"
#include "filename.h"
#include "Rot.hh"

/* LowParser - begin */

static int
skip_remarks(HighParser& HP, InputStream& In, char &cIn)
{
skip_again:;
	for (cIn = In.get(); isspace(cIn); cIn = In.get()) {
		if (In.eof()) {
			return -1;
		}
	}

	switch (cIn) {
	case -1:
		return -1;
	
	case MathParser::ONE_LINE_REMARK:
		for (cIn = In.get(); cIn != '\n'; cIn = In.get()) {
			if (cIn == '\\') {
				cIn = In.get();
				if (In.eof()) {
					return -1;
				}
				if (cIn == '\r') {
					/* if input file was prepared
					 * under DOS/Windows */
					cIn = In.get();
					if (In.eof()) {
						return -1;
					}
				}
			}
			if (In.eof()) {
				return -1;
			}
		}
		goto skip_again;

	case '/':
		cIn = In.get();
		if (In.eof()) {
			return -1;

		} else if (cIn == '*') {
			for (cIn = In.get(); !In.eof(); cIn = In.get()) {
				if (cIn == '*') {
end_of_comment:;
					cIn = In.get();
					if (In.eof()) {
						return -1;
					}
					if (cIn == '/') {
						goto skip_again;
					}

				} else if (cIn == '/') {
					cIn = In.get();
					if (In.eof()) {
						return -1;
					}
					if (cIn == '*') {
						silent_cerr("warning: '/*' inside a comment "
							"at line " << HP.GetLineData()
							<< std::endl);
						goto end_of_comment;
					}
				}
			}
			if (In.eof()) {
				return -1;
			}

		} else {
			In.putback(cIn);
			return 0;
		}
	}

	return 0;
}

LowParser::LowParser(HighParser& hp)
: HP(hp), sCurrWordBuf(0), iBufSize(iDefaultBufSize)
{
	SAFENEWARR(sCurrWordBuf, char, iBufSize);
}

LowParser::~LowParser(void)
{
	if (sCurrWordBuf) {
		SAFEDELETEARR(sCurrWordBuf);
	}
}

void
LowParser::PackWords(InputStream& In)
{
	unsigned iCur = 0;
	char cIn;

	/* note: no remarks allowed inside words */
	for (cIn = In.get(); !In.eof(); cIn = In.get()) {
		switch (cIn) {
		case COLON:
		case COMMA:
		case SEMICOLON:
			goto end_of_word;

		default:
      			if (!isspace(cIn)) {
				sCurrWordBuf[iCur] = cIn;
				iCur++;
				 if (iCur == iBufSize - 1) {
					 char *s = NULL;
					 unsigned i = 2*iBufSize;

					 /* FIXME: no limit on max size? */

					 SAFENEWARR(s, char, i);
					 memcpy(s, sCurrWordBuf, iBufSize);
					 SAFEDELETEARR(sCurrWordBuf);
					 sCurrWordBuf = s;
					 iBufSize = i;
				 }
			}
		}
	}

	throw EndOfFile(MBDYN_EXCEPT_ARGS);

end_of_word:;

	sCurrWordBuf[iCur] = '\0';
	In.putback(cIn);
}


LowParser::Token
LowParser::GetToken(InputStream& In)
{
	/* toglie gli spazi iniziali e tutti i commenti */
	char cIn;
	if (skip_remarks(HP, In, cIn)) {
		return CurrToken = LowParser::ENDOFFILE;
	}

	if (isalpha(cIn) || cIn == '_') {
		PackWords(In.putback(cIn));
		return CurrToken = LowParser::WORD;
	}

	switch (cIn) {
	case ',':
		return CurrToken = LowParser::COMMA;

	case ':':
		return CurrToken = LowParser::COLON;

	case ';':
		return CurrToken = LowParser::SEMICOLON;

	case '.':
	case '-':
	case '+':
is_digit:;
		In.putback(cIn) >> dCurrNumber;
		return CurrToken = LowParser::NUMBER;

	default:
		if (isdigit(cIn)) {
			goto is_digit;
		}
		In.putback(cIn);
		return CurrToken = LowParser::UNKNOWN;
	}
}


doublereal
LowParser::dGetReal(void) const
{
	return dCurrNumber;
}


integer
LowParser::iGetInt(void) const
{
	return integer(dCurrNumber);
}


char*
LowParser::sGetWord(void)
{
	return sCurrWordBuf;
}

/* LowParser - end */


/* KeyTable - begin */

KeyTable::KeyTable(HighParser& hp, const char* const sTable[])
: sKeyWords(0), oldKey(0), HP(hp) 
{
	sKeyWords = (char* const*)sTable;
	oldKey = HP.PutKeyTable(*this);
}


KeyTable::~KeyTable(void)
{
	if (oldKey) {
		(void)HP.PutKeyTable(*oldKey);
	}
}


int
KeyTable::Find(const char* sToFind) const
{
	for (int iCnt = 0; sKeyWords[iCnt]; iCnt++) {
		if (strcasecmp(sKeyWords[iCnt], sToFind) == 0) {
			return iCnt;
		}
	}

	return -1;
}

/* KeyTable - end */


/* DescRead - begin */

/* bag that contains functions to parse descriptions */

typedef std::map<std::string, DescRead *, ltstrcase> DescFuncMapType;
static DescFuncMapType DescFuncMap;

struct DescWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return DescFuncMap.find(s) != DescFuncMap.end();
	};
};

static DescWordSetType DescWordSet;

bool
SetDescData(const std::string& name, DescRead *rf)
{
	pedantic_cout("registering description \"" << name << "\"" << std::endl);
	return DescFuncMap.insert(DescFuncMapType::value_type(name, rf)).second;
}

/* Reads descriptions */

bool
ReadDescription(HighParser& HP, const std::string& desc)
{
	DEBUGCOUTFNAME("ReadDescription()");

	bool bRC(false);
	DescFuncMapType::iterator func = DescFuncMap.find(desc);
	if (func != DescFuncMap.end()) {
		HP.GotDescription();
		if (!HP.IsArg() && !HP.IsDescription()) {
			silent_cerr("Parser error in ReadDescription(),"
				" colon or semicolon expected after description at line "
				<< HP.GetLineData() << std::endl);
			throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
		}

		bRC = func->second->Read(HP);

		if (HP.IsArg()) {
			silent_cerr("semicolon expected at line " << HP.GetLineData() << std::endl);
			throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
		}
	}

	return bRC;
}

DescRead::~DescRead(void)
{
	NO_OP;
}

struct RemarkDR : public DescRead {
public:
	bool Read(HighParser& HP);
};

bool
RemarkDR::Read(HighParser& HP)
{
	if (!HP.IsArg()) {
		silent_cerr("Parser error in RemarkDR::Read(),"
			" arg(s) expected at line "
			<< HP.GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	/* eat it anyway ;) */
	const char *s = HP.GetStringWithDelims();
	silent_cout("line " << HP.GetLineData() << ": " << s);
	bool bFirst(true);
	while (HP.IsArg()) {
		if (bFirst) {
			silent_cout(": ");
			bFirst = false;

		} else {
			silent_cout(", ");
		}

#if 0
		bool bIgnore(false);
		if (HP.IsKeyWord("ignore")) {
			bIgnore = true;
		}
#endif

		TypedValue v;
		v = HP.GetValue(v);
		silent_cout(v);
	}

	silent_cout(std::endl);

	return true;
}

struct PrintSymbolTableDR : public DescRead {
public:
	bool Read(HighParser& HP);
};

bool
PrintSymbolTableDR::Read(HighParser& HP)
{
	if (!HP.IsArg()) {
		silent_cout( "math parser symbol table at line "
			<< HP.GetLineData() << ":" << std::endl
			<< HP.GetMathParser().GetSymbolTable() << std::endl);
		return true;
	}

	if (HP.IsKeyWord("all")) {
		const MathParser::NameSpaceMap& ns = HP.GetMathParser().GetNameSpaceMap();
		for (MathParser::NameSpaceMap::const_iterator i = ns.begin(); i != ns.end(); ++i) {
			const std::string& sName = i->second->sGetName();
			const Table *pT = i->second->GetTable();
			if (pT != 0) {
				silent_cout( "namespace \"" << sName << "\" symbol table at line "
					<< HP.GetLineData() << ":" << std::endl
					<< *pT << std::endl);
			}
		}

		return true;
	}

	while (HP.IsArg()) {
		const char *sName = HP.GetString();
		MathParser::NameSpace *pN = HP.GetMathParser().GetNameSpace(sName);
		if (pN == 0) {
			silent_cerr("PrintSymbolTableDR::Read(): warning, unable to find namespace \"" << sName << "\" at line " 
				<< HP.GetLineData() << std::endl);

		} else {
			Table *pT = pN->GetTable();
			if (pT == 0) {
				silent_cerr("PrintSymbolTableDR::Read(): warning, namespace \"" << sName << "\" "
					"has no symbol table at line " << HP.GetLineData() << std::endl);

			} else {
				silent_cout( "namespace \"" << sName << "\" symbol table at line "
					<< HP.GetLineData() << ":" << std::endl
					<< *pT << std::endl);
			}
		}
	}

	return true;
}

struct SetDR : public DescRead {
public:
	bool Read(HighParser& HP);
};

bool
SetDR::Read(HighParser& HP)
{
	if (!HP.IsArg()) {
     		silent_cerr("Parser error in SetDR::Read(), "
     			"arg expected at line "
     			<< HP.GetLineData() << std::endl);
     		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	TypedValue v;
	HP.GetValue(v);

	return true;
}

struct SetEnvDR : public DescRead {
public:
	bool Read(HighParser& HP);
};

bool
SetEnvDR::Read(HighParser& HP)
{
#ifdef HAVE_SETENV
	if (!HP.IsArg()) {
     		silent_cerr("Parser error in SetEnvDR::Read(), "
     			"arg(s) expected at line "
     			<< HP.GetLineData() << std::endl);
     		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	int overwrite = 0;
	if (HP.IsKeyWord("overwrite")) {
		bool b = HP.GetYesNoOrBool();
		overwrite = b ? 1 : 0;
	}

	const char *ava = HP.GetStringWithDelims();
	if (ava == NULL) {
		silent_cerr("unable to get AVA for \"setenv\" at line "
				<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	char *avasep = std::strchr(const_cast<char *>(ava), '=');
	if (avasep == NULL) {
#ifdef HAVE_UNSETENV
		unsetenv(ava);
#elif defined(HAVE_PUTENV)
		if (putenv(ava)) {
			silent_cerr("unable to unset the environment variable "
					"\"" << ava << "\" at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#endif	/* !HAVE_UNSETENV && !HAVE_PUTENV */

	} else {
		if (avasep == ava) {
			silent_cerr("illegal AVA \"" << ava
					<< "\" at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		avasep[0] = '\0';
		avasep++;
		bool bPresent(getenv(ava) != NULL);
		int rc = setenv(ava, avasep, overwrite);
		if (rc) {
			silent_cerr("unable to set the environment variable \""
					<< ava << "\" to \"" << avasep 
					<< "\" at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bPresent && overwrite == 0) {
			silent_cout("Environment variable \"" << ava
				<< "\" _not_ overwritten with \"" << avasep
				<< "\" (current value is \"" << getenv(ava)
				<< "\") at line " << HP.GetLineData()
				<< std::endl);

		} else if (!bPresent) {
			silent_cout("Environment variable \"" << ava
				<< "\" set to \"" << avasep
				<< "\" at line " << HP.GetLineData()
				<< std::endl);

		} else {
			silent_cout("Environment variable \"" << ava
				<< "\" overwritten to \"" << avasep
				<< "\" at line " << HP.GetLineData()
				<< std::endl);
		}
	}
#else // ! HAVE_SETENV
	silent_cerr("SetEnvDR::Read(): warning, setenv() not available; ignored at line "
		<< HP.GetLineData() << std::endl);
#endif // !HAVE_SETENV
	return true;
}

struct ExitDR : public DescRead {
public:
	bool Read(HighParser& HP);
};

bool
ExitDR::Read(HighParser& HP)
{
	if (!HP.IsDescription()) {
		silent_cerr("Parser error in ExitDR::Read(),"
			" semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
	}

	/* exits with no error */
	throw NoErr(MBDYN_EXCEPT_ARGS);
}

static unsigned desc_done;

static void
InitDescData(void)
{
	if (::desc_done++ > 0) {
		return;
	}

	SetDescData("remark", new RemarkDR);
	SetDescData("print" "symbol" "table", new PrintSymbolTableDR);
	SetDescData("set", new SetDR);
	SetDescData("setenv", new SetEnvDR);
	SetDescData("exit", new ExitDR);

	/* NOTE: add here initialization of new built-in descriptions;
	 * alternative ways to register new custom descriptions are:
	 * - call SetDescData() from anywhere in the code
	 * - write a module that calls SetDescData() from inside a function
	 *   called module_init(), and run-time load it using "module load"
	 *   in the input file.
	 */
}

static void
DestroyDescData(void)
{
	if (::desc_done == 0) {
		silent_cerr("DestroyDescData() called once too many" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--::desc_done > 0) {
		return;
	}

	/* free stuff */
	for (DescFuncMapType::iterator i = DescFuncMap.begin(); i != DescFuncMap.end(); ++i) {
		delete i->second;
	}
	DescFuncMap.clear();
}

/* DescRead - end */


/* HighParser - begin */

static std::stack<const HighParser *> pHP;

static const HighParser::ErrOut unknownErr = { "(unknown)", "(unknown)", 0 };

HighParser::ErrOut
mbdyn_get_line_data(void)
{
	if (!pHP.empty()) {
		return pHP.top()->GetLineData();
	}

	return unknownErr;
}

std::ostream&
mbdyn_print_line_data(std::ostream& out)
{
	if (!pHP.empty()) {
		out << pHP.top()->GetLineData();
	}

	return out;
}

HighParser::HighParser(MathParser& MP, InputStream& streamIn)
: ESCAPE_CHAR('\\'),
LowP(*this),
pIn(&streamIn),
pf(NULL),
MathP(MP),
KeyT(0)
{
	DEBUGCOUTFNAME("HighParser::HighParser");
	CurrToken = HighParser::DESCRIPTION;

	InitDescData();

	pHP.push(this);
}


HighParser::~HighParser(void)
{
	DEBUGCOUTFNAME("HighParser::~HighParser");
	Close();
	ASSERT(pHP.top() == this);
	pHP.pop();

	DestroyDescData();
}


void
HighParser::Close(void)
{
	NO_OP;
}


const KeyTable*
HighParser::PutKeyTable(const KeyTable& KT)
{
	const KeyTable* oldKey = KeyT;

	KeyT = &KT;

	return oldKey;
}

MathParser&
HighParser::GetMathParser(void)
{
	return MathP;
}

int
HighParser::GetLineNumber(void) const
{
	return ((InputStream*)pIn)->GetLineNumber();
}


HighParser::ErrOut
HighParser::GetLineData(void) const
{
	ErrOut LineData;
	LineData.iLineNumber = GetLineNumber();
	LineData.sFileName = NULL;
	LineData.sPathName = NULL;
	return LineData;
}


bool
HighParser::IsDescription(void) const
{
     	return (CurrToken == HighParser::DESCRIPTION);
}

HighParser::Token
HighParser::GotDescription(void)
{
	return FirstToken();
}

int
HighParser::iGetDescription_int(const char* const s)
{
	int i = -1;
	
	if (KeyT) {
		i = KeyT->Find(s);
	}

     	if (FirstToken() == HighParser::UNKNOWN) {
		silent_cerr("Parser error in HighParser::iGetDescription_int(), "
			"semicolon expected at line "
			<< GetLineData() << std::endl);
	  	throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
     	}

     	return i;
}


void
HighParser::Eof(void)
{
	 throw EndOfFile(MBDYN_EXCEPT_ARGS);
}

int
HighParser::GetDescription(void)
{
	/* Checks if current token is a description */
	if (!IsDescription()) {
		silent_cerr("Parser error in HighParser::GetDescription, "
			"invalid call to GetDescription at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrInvalidCallToGetDescription(MBDYN_EXCEPT_ARGS);
	}

restart_parsing:;

	CurrLowToken = LowP.GetToken(*pIn);
	if (CurrLowToken != LowParser::WORD) {
		if (pIn->eof()) {
			Eof();
			goto restart_parsing;
		}
		
		silent_cerr("Parser error in HighParser::GetDescription, "
			<< "keyword expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
	}

	/* Description corrente */
	char* s = LowP.sGetWord();

 	if (ReadDescription(*this, s)) {
		goto restart_parsing;
	}

	return iGetDescription_int(s);
}


HighParser::Token
HighParser::FirstToken(void)
{
	CurrLowToken = LowP.GetToken(*pIn);

	switch (CurrLowToken) {
	case LowParser::COLON:
		CurrToken = HighParser::ARG;
		break;

	case LowParser::SEMICOLON:
		CurrToken = HighParser::DESCRIPTION;
		break;

	default:
		CurrToken = HighParser::UNKNOWN;
		break;
	}

	return CurrToken;
}

void
HighParser::ExpectDescription(void)
{
	/* forces the next expected token to be a "description"
	 * e.g. a keyword followed by a colon (deprecated) */
	CurrToken = HighParser::DESCRIPTION;
}


void
HighParser::ExpectArg(void)
{
	/* forces the next expected token to be an argument
	 * e.g. a keyword followed by a separator (deprecated) */
	CurrToken = HighParser::ARG;
}


bool
HighParser::IsArg(void)
{
	return (CurrToken == ARG);
}

void
HighParser::PutBackSemicolon(void)
{
	if (CurrLowToken == LowParser::SEMICOLON) {
		pIn->putback(';');
	}
}


void
HighParser::NextToken(const char* sFuncName)
{
	CurrLowToken = LowP.GetToken(*pIn);
	switch (CurrLowToken) {
	case LowParser::COMMA:
		CurrToken = HighParser::ARG;
		break;

	case LowParser::SEMICOLON:
		CurrToken = HighParser::DESCRIPTION;
		break;

	default:
		silent_cerr("Parser error in "
			<< sFuncName << ", missing separator at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrMissingSeparator(MBDYN_EXCEPT_ARGS);
	}
}

int
HighParser::ParseWord(unsigned flags)
{
	char* sBuf = sStringBuf;
	char* sBufWithSpaces = sStringBufWithSpaces;

	char cIn;
	if (skip_remarks(*this, *pIn, cIn)) {
		return CurrToken = HighParser::ENDOFFILE;
	}

	if (!isalpha(cIn) && cIn != '_') {
		pIn->putback(cIn);
		return -1;
	}

	*sBufWithSpaces++ = cIn;

	if (flags & LOWER) {
		*sBuf++ = tolower(cIn);

	} else if (flags & UPPER) {
		*sBuf++ = toupper(cIn);

	} else {
		*sBuf++ = cIn;
	}

	for (cIn = pIn->get(); isalnum(cIn) || cIn == '_' || isspace(cIn); cIn = pIn->get()) {
		*sBufWithSpaces++ = cIn;
		if (sBufWithSpaces >= sStringBufWithSpaces + iDefaultBufSize - 1) {
			break;
		}

		if (isspace(cIn)) {
			continue;
		}
		
		if (flags & LOWER) {
			*sBuf++ = tolower(cIn);

		} else if (flags & UPPER) {
			*sBuf++ = toupper(cIn);

		} else {
			*sBuf++ = cIn;
		}
	}
	pIn->putback(cIn);

	*sBuf = '\0';
	*sBufWithSpaces = '\0';

	return 0;
}

void
HighParser::PutbackWord(void)
{
	char* sBufWithSpaces = sStringBufWithSpaces + strlen(sStringBufWithSpaces);

   
	while (sBufWithSpaces > sStringBufWithSpaces) {
		pIn->putback(*--sBufWithSpaces);
	}
}

bool
HighParser::IsKeyWord(const char* sKeyWord)
{
	const char sFuncName[] = "HighParser::IsKeyWord()";

	if (CurrToken != HighParser::ARG) {
		return false;
	}

	switch (ParseWord()) {
	case 0:
		break;

	case HighParser::ENDOFFILE:
		return true;

	default:
		return false;
	}

	if (!strcasecmp(sStringBuf, sKeyWord)) {
		NextToken(sFuncName);
		return true;
	}

	PutbackWord();

	return false;
}

int
HighParser::IsKeyWord(void)
{
	const char sFuncName[] = "HighParser::IsKeyWord()";

	if (CurrToken != HighParser::ARG) {
		return -1;
	}

	switch (ParseWord()) {
	case 0:
		break;

	case HighParser::ENDOFFILE:
		return HighParser::ENDOFFILE;

	default:
		return -1;
	}

	int iKW = -1;

	if (KeyT) {
		iKW = KeyT->Find(sStringBuf);
	}
   
	if (iKW >= 0) {
		NextToken(sFuncName);
		return iKW;
	}

	PutbackWord();

	return -1;
}

/* 1 se l'argomento successivo e' una parola in un WordSet */
const char *
HighParser::IsWord(const HighParser::WordSet& ws)
{
	const char sFuncName[] = "HighParser::IsWord()";

	if (CurrToken != HighParser::ARG) {
		return 0;
	}

	switch (ParseWord()) {
	case 0:
		break;

	default:
		return 0;
	}

	if (ws.IsWord(std::string(sStringBuf))) {
		NextToken(sFuncName);
		return sStringBuf;
	}

	PutbackWord();

	return 0;
}

TypedValue
HighParser::GetValue(const TypedValue& vDefVal)
{
	const char sFuncName[] = "HighParser::GetValue()";

	if (CurrToken != HighParser::ARG) {
		silent_cerr("Parser error in "
			<< sFuncName << ", arg expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrIntegerExpected(MBDYN_EXCEPT_ARGS);
	}

	TypedValue v(vDefVal);

	try {
		v = MathP.Get(*pIn, v);
	}
	catch (TypedValue::ErrWrongType e) {
		silent_cerr(sFuncName << ": " << e.what() << " at line "
			<< GetLineData() << std::endl);
		throw e;
	}
	catch (MathParser::ErrGeneric e) {
		silent_cerr(sFuncName << ": error return from MathParser at line "
			<< GetLineData() << std::endl);
		throw e;
	}
	catch (...) {
		throw;
	}

	NextToken(sFuncName);
	return v;
}

bool
HighParser::GetBool(bool bDefVal)
{
	TypedValue v(bDefVal);
	v = GetValue(v);
	return v.GetBool();
}

/*
 * Read a boolean as "yes" or "no" and put the result in bRet
 * return true in case of success, false otherwise (bRet undefined)
 */
bool
HighParser::GetYesNo(bool& bRet)
{
	if (IsKeyWord("yes")) {
		bRet = true;

	} else if (IsKeyWord("no")) {
		bRet = false;

	} else {
		return false;
	}

	return true;
}

bool
HighParser::GetYesNoOrBool(bool bDefval)
{
	bool bRet;

	if (!GetYesNo(bRet)) {
		bRet = GetBool(bDefval);
	}

	return bRet;
}

integer
HighParser::GetInt(int iDefVal)
{
	TypedValue v(iDefVal);
	v = GetValue(v);
	return v.GetInt();
}

doublereal
HighParser::GetReal(const doublereal& dDefVal)
{
	TypedValue v(dDefVal);
	v = GetValue(v);
	return v.GetReal();
}

mbsleep_t
HighParser::GetTimeout(const mbsleep_t& DefVal)
{
	doublereal d;
	mbsleep_sleep2real(&DefVal, &d);
	TypedValue v(d);
	v = GetValue(v);
	mbsleep_t newval;
	mbsleep_real2sleep(v.GetReal(), &newval);
	return newval;
}

std::string
HighParser::GetString(const std::string& sDefVal)
{
	TypedValue v(sDefVal);
	v = GetValue(v);
	return v.GetString();
}


int
HighParser::GetWord(void)
{
	const char sFuncName[] = "HighParser::GetWord()";

	if (CurrToken != HighParser::ARG) {
		silent_cerr("Parser error in "
			<< sFuncName << ", keyword arg expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
	}

	CurrLowToken = LowP.GetToken(*pIn);
	if (CurrLowToken != LowParser::WORD) {
		silent_cerr("Parser error in "
			<< sFuncName << ", keyword expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
	}

	int i = -1;
	if (KeyT) {
		i = KeyT->Find(LowP.sGetWord());
	}

	NextToken(sFuncName);

	return i;
}

const char*
HighParser::GetString(unsigned flags)
{
	const char sFuncName[] = "HighParser::GetString()";

	pedantic_cout("use of deprecated method \"GetString\" at line"
		<< GetLineData() << std::endl);

	if (CurrToken != HighParser::ARG) {
		silent_cerr("Parser error in "
			<< sFuncName << ", string arg expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrStringExpected(MBDYN_EXCEPT_ARGS);
	}

	char* s = sStringBuf;
	char* sTmp = s;

	char cIn = '\0';

	while (isspace(cIn = pIn->get())) {
		NO_OP;
	}

	if (pIn->eof()) {
		CurrToken = HighParser::ENDOFFILE;
		return NULL;
	}

	pIn->putback(cIn);
	for (cIn = pIn->get(); cIn != ',' && cIn != ';'; cIn = pIn->get()) {
		/* Attenzione! cosi' la legge tutta,
		 * ma ne tiene solo iBufSize-1 caratteri */
		if (pIn->eof()) {
			CurrToken = HighParser::ENDOFFILE;
			*sTmp = '\0';
			return s;

		} else if (sTmp < s + iDefaultBufSize - 1) {
			if (!(flags & HighParser::EATSPACES) || !isspace(cIn)) {
				if (flags & HighParser::LOWER) {
					cIn = tolower(cIn);

				} else if (flags & HighParser::UPPER) {
					cIn = toupper(cIn);
				}

				*sTmp++ = cIn;
			}
		}
	}

	pIn->putback(cIn);
	*sTmp = '\0';

	NextToken(sFuncName);

	return s;
}

void
HighParser::SetDelims(enum Delims Del, char &cLdelim, char &cRdelim) const
{
	cLdelim = '\0';
	cRdelim = '\0';

	switch (Del) {
	case PLAINBRACKETS:
		cLdelim = '(';
		cRdelim = ')';
		break;

	case SQUAREBRACKETS:
		cLdelim = '[';
		cRdelim = ']';
		break;

	case CURLYBRACKETS:
		cLdelim = '{';
		cRdelim = '}';
		break;

	case SINGLEQUOTE:
		cLdelim = '`';
		cRdelim = '\'';
		break;

	default:
	case UNKNOWNDELIM:
	case DEFAULTDELIM:
	case DOUBLEQUOTE:
		cLdelim = '"';
		cRdelim = '"';
		break;
	}
}

bool
HighParser::IsStringWithDelims(enum Delims Del)
{
	char cLdelim, cRdelim;
	SetDelims(Del, cLdelim, cRdelim);

	char cIn;
	if (skip_remarks(*this, *pIn, cIn)) {
		return false;
	}

	/* put back the first non-remark char */
	pIn->putback(cIn);

	/* if the left delimiter is found, true */
	return (cIn == cLdelim);
}

const char*
HighParser::GetStringWithDelims(enum Delims Del, bool escape)
{
	const char sFuncName[] = "HighParser::GetStringWithDelims()";

	if (CurrToken != HighParser::ARG) {
		silent_cerr("Parser error in "
			<< sFuncName << ", string arg expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrStringExpected(MBDYN_EXCEPT_ARGS);
	}

	char* s = sStringBuf;
	char* sTmp = s;

	char cLdelim, cRdelim;
	SetDelims(Del, cLdelim, cRdelim);

	char cIn;
	if (skip_remarks(*this, *pIn, cIn)) {
		return NULL;
	}

	/* Se trova il delimitatore sinistro, legge la stringa */
	if (cIn == cLdelim) {
		for (cIn = pIn->get(); cIn != cRdelim; cIn = pIn->get()) {
			/* Attenzione! cosi' la legge tutta,
			 * ma ne tiene solo iBufSize-1 caratteri */
			if (pIn->eof()) {
				/* FIXME: this should be an error ... */
				sTmp[0] = '\0';
				return s;

			} else if (sTmp < s + iDefaultBufSize - 1) {
				if (cIn == ESCAPE_CHAR) {
					cIn = pIn->get();
					if (cIn == '\n') {

						/*
						 * eat the newline as well, so that

							"first line\
							second line"

						 * actually results in "first linesecond line"
						 */

						cIn = pIn->get();

					} else if (cIn == '\r') {
						cIn = pIn->get();
						if (cIn != '\n') {
							pIn->putback(cIn);
							goto escaped_generic;
						}
						cIn = pIn->get();

					} else if ((cIn == ESCAPE_CHAR) || (cIn == cRdelim)) {
						if (!escape) {
							sTmp[0] = ESCAPE_CHAR;
							++sTmp;
						}

					} else {
escaped_generic:;
						if (escape) {
							int i, c = 0;
							char hex[3];

							/*
							 * allow non-printable chars in the form "\<hexpair>",
							 * so that "\78" is equivalent to "x";
							 * "\<non-hexpair>" is treated as an error.
							 */

							hex[0] = cIn;
							hex[1] = pIn->get();
							hex[2] = '\0';

							for (i = 0; i < 2; i++) {
								int shift = 4*(1 - i), h = 0;

								/* NOTE: this conversion relies
								 * on 0-9, a-f, A-F being consecutive,
								 * which is true for ASCII, but might
								 * not be for other encodings;
								 * bah, not critical right now */
								if (hex[i] >= '0' && hex[i]  <= '9') {
									h = hex[i] - '0';
								} else if (hex[i] >= 'a' && hex[i] <= 'f') {
									h = hex[i] - 'a';
								} else if (hex[i] >= 'A' && hex[i] <= 'F') {
									h = hex[i] - 'A';
								} else {
									silent_cerr("invalid escape sequence "
										"\"\\" << hex << "\" "
										"at line " << GetLineData()
										<< std::endl);
									throw ErrGeneric(MBDYN_EXCEPT_ARGS);
								}

								c += (h << shift);
							}
							cIn = c;

						} else {
							sTmp[0] = ESCAPE_CHAR;
							++sTmp;
						}
					}
				}
				sTmp[0] = cIn;
				++sTmp;
			}
		}

		/* Se trova una virgola o un punto e virgola, le rimette nello stream
		 * e passa oltre, restituendo un puntatore nullo. Il chiamante deve
		 * occuparsi della gestione del valore di default */
	} else if (cIn == ',' || cIn == ';') {
		pIn->putback(cIn);
		goto nullstring;

		/* Altrimenti c'e' qualcosa senza delimitatore. Adesso da' errore,
		 * forse e' piu' corretto fargli ritornare lo stream intatto */
	} else {
		silent_cerr("Parser error in "
			<< sFuncName << std::endl
			<< "first non-blank char at line "
			<< GetLineData() << " isn't a valid left-delimiter"
			<< std::endl);
		throw HighParser::ErrIllegalDelimiter(MBDYN_EXCEPT_ARGS);
	}

	/* Mette zero al termine della stringa */
	*sTmp = '\0';

nullstring:;
	NextToken(sFuncName);
	return s;
}

/* HighParser - end */

std::ostream&
operator << (std::ostream& out, const HighParser::ErrOut& err)
{
	out << err.iLineNumber;

	if (err.sFileName != 0) {
		out << ", file <";
		if (err.sPathName != 0) {
			out << err.sPathName << DIR_SEP;
		}
		out << err.sFileName << '>';
	}

	return out;
}

/* HighParser - end */

