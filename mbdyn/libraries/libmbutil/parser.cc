/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <parser.h>

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

	if (cIn == REMARK) {
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
				if (cIn != '\n') {
					In.putback(cIn);
				}
			}
			if (In.eof()) {
				return -1;
			}
		}
		goto skip_again;

	} else if (cIn == '/') {
		cIn = In.get();
		if (In.eof()) {
			return -1;

		} else if (cIn == '*') {
			for (; !In.eof(); cIn = In.get()) {
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
						silent_cerr("warning: '/*' inside a comment at line " << HP.GetLineData() << std::endl);
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

	THROW(EndOfFile());

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
LowParser::dGetReal(void)
{
	return dCurrNumber;
}


integer
LowParser::iGetInt(void)
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


/* HighParser - begin */

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
}


HighParser::~HighParser(void)
{
	DEBUGCOUTFNAME("HighParser::~HighParser");
	Close();
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
	return LineData;
}


bool
HighParser::IsDescription(void)
{
     	if (CurrToken != HighParser::DESCRIPTION) {
	  	return false;
     	}
     	return true;
}


int
HighParser::iGetDescription_(const char* const s)
{
	int i = -1;
	
	if (KeyT) {
		i = KeyT->Find(s);
	}

     	if (FirstToken() == HighParser::UNKNOWN) {
		silent_cerr("Parser error in HighParser::iGetDescription_(), "
			"semicolon expected at line "
			<< GetLineData() << std::endl);
	  	THROW(HighParser::ErrSemicolonExpected());
     	}

     	return i;
}


void
HighParser::Set_(void)
{
     	if (FirstToken() == UNKNOWN) {
     		silent_cerr("Parser error in HighParser::Set_(), "
     			"colon expected at line "
     			<< GetLineData() << std::endl);
     		THROW(HighParser::ErrColonExpected());
	}

	GetReal();
}

void
HighParser::SetEnv_(void)
{
     	if (FirstToken() == UNKNOWN) {
     		silent_cerr("Parser error in HighParser::SetEnv_(), "
     			"colon expected at line "
     			<< GetLineData() << std::endl);
     		THROW(HighParser::ErrColonExpected());
	}

	int overwrite = 0;
	if (IsKeyWord("overwrite")) {
		overwrite = 1;
	}

	const char *ava = GetStringWithDelims();
	if (ava == NULL) {
		silent_cerr("unable to get AVA for \"setenv\" at line "
				<< GetLineData() << std::endl);
		THROW(ErrGeneric());
	}

	char *avasep = strchr(ava, '=');
	if (avasep == NULL) {
		int	rc = 0;
#ifdef HAVE_UNSETENV
		rc = unsetenv(ava);
#elif defined(HAVE_PUTENV)
		rc = putenv(ava);
#endif	/* !HAVE_UNSETENV && !HAVE_PUTENV */
		if (rc) {
			silent_cerr("unable to unset the environment variable "
					"\"" << ava << "\" at line "
					<< GetLineData() << std::endl);
			THROW(ErrGeneric());
		}


	} else {
		if (avasep == ava) {
			silent_cerr("illegal AVA \"" << ava
					<< "\" at line "
					<< GetLineData() << std::endl);
			THROW(ErrGeneric());
		}

		unsigned l = avasep - ava;
		char buf[l + 1];
		memcpy(buf, ava, l);
		buf[l] = '\0';
		avasep++;
		if (setenv(buf, avasep, overwrite)) {
			silent_cerr("unable to set the environment variable \""
					<< buf << "\" to \"" << avasep 
					<< "\" at line " << GetLineData()
					<< std::endl);
			THROW(ErrGeneric());
		}
	}
}

void
HighParser::Remark_(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::Remark_(),"
			" colon expected at line "
			<< GetLineData() << std::endl);
		THROW(HighParser::ErrColonExpected());
	}

	/* eat it anyway ;) */
	const char *s = GetStringWithDelims();
	silent_cout("line " << GetLineData() << ": " << s);
	while (IsArg()) {
		/* eat it anyway ;) */
		double d = GetReal();
		silent_cout(", " << d);
	}

	silent_cout(std::endl);
}


void
HighParser::Eof(void)
{
	 THROW(EndOfFile());
}

int
HighParser::GetDescription(void)
{
	/* Checks if current token is a description */
	if (!IsDescription()) {
		silent_cerr("Parser error in HighParser::GetDescription, "
			"invalid call to GetDescription at line "
			<< GetLineData() << std::endl);
		THROW(HighParser::ErrInvalidCallToGetDescription());
	}

restart_parsing:;

	CurrLowToken = LowP.GetToken(*pIn);
	if (CurrLowToken != LowParser::WORD) {
		if (pIn->GetStream().eof()) {
			Eof();
			goto restart_parsing;

		}
		
		silent_cerr("Parser error in HighParser::GetDescription, "
			<< "keyword expected at line "
			<< GetLineData() << std::endl);
		THROW(HighParser::ErrKeyWordExpected());
	}

	/* Description corrente */
	char* s = LowP.sGetWord();

	if (GetDescription_int(s)) {
		goto restart_parsing;
	}

	return iGetDescription_(s);
}


/*
 * special descriptions, that are eaten by the parser
 */
bool
HighParser::GetDescription_int(const char *s)
{
	/* calls the MathParser */
	if (strcmp(s, "set") == 0) {
		Set_();
		return true;

	/* sets environment variable */
	} else if (strcmp(s, "setenv") == 0) {
		SetEnv_();
		return true;

	/* exits with no error */
	} else if (strcmp(s, "exit") == 0) {
		THROW(NoErr());
	}

	return false;
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
	if (CurrToken == ARG) {
		return true;
	}
	return false;
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
		THROW(HighParser::ErrMissingSeparator());
	}
}


bool
HighParser::IsKeyWord(const char* sKeyWord)
{
	const char sFuncName[] = "HighParser::IsKeyWord()";

	if (CurrToken != HighParser::ARG) {
		return false;
	}

	char* sBuf = sStringBuf;
	char* sBufWithSpaces = sStringBufWithSpaces;

	char cIn;
	if (skip_remarks(*this, *pIn, cIn)) {
		CurrToken = HighParser::ENDOFFILE;
		return true;
	}

	if (!isalpha(cIn)) {
		pIn->putback(cIn);
		return false;
	}

	*sBuf++ = cIn;
	*sBufWithSpaces++ = cIn;

	/* Forse e' meglio modificare in modo che digerisca anche gli spazi
	 * tra due parole, magari con due buffer, uno in cui li mangia per fare
	 * il confronto con la keyword, l'altro in cui li tiene per l'eventuale
	 * putback */
	for (cIn = pIn->get(); isalnum(cIn) || isspace(cIn); cIn = pIn->get()) {
		*sBufWithSpaces++ = cIn;
		if (isalnum(cIn)) {
			*sBuf++ = cIn;
		}
		if (sBufWithSpaces >= sStringBufWithSpaces + iDefaultBufSize - 1) {
			break;
		}
	}
	pIn->putback(cIn);

	*sBuf = '\0';
	*sBufWithSpaces = '\0';

	if (!strcasecmp(sStringBuf, sKeyWord)) {
		NextToken(sFuncName);
		return true;
	}

	while (sBufWithSpaces > sStringBufWithSpaces) {
		pIn->putback(*--sBufWithSpaces);
	}

	return false;
}


int
HighParser::IsKeyWord(void)
{
	const char sFuncName[] = "HighParser::IsKeyWord()";

	if (CurrToken != HighParser::ARG) {
		return -1;
	}

	char* sBuf = sStringBuf;
	char* sBufWithSpaces = sStringBufWithSpaces;

	char cIn;
	if (skip_remarks(*this, *pIn, cIn)) {
		return CurrToken = HighParser::ENDOFFILE;
	}

	if (!isalpha(cIn)) {
		pIn->putback(cIn);
		return -1;
	}
	*sBuf++ = cIn;
	*sBufWithSpaces++ = cIn;

	for (cIn = pIn->get(); isalnum(cIn) || isspace(cIn); cIn = pIn->get()) {
		*sBufWithSpaces++ = cIn;
		if (isalnum(cIn)) {
			*sBuf++ = cIn;
		}
		if (sBufWithSpaces >= sStringBufWithSpaces + iDefaultBufSize - 1) {
			break;
		}
	}
	pIn->putback(cIn);

	*sBuf = '\0';
	*sBufWithSpaces = '\0';

	int iKW = -1;

	if (KeyT) {
		iKW = KeyT->Find(sStringBuf);
	}
   
	if (iKW >= 0) {
		NextToken(sFuncName);
		return iKW;
	}

	while (sBufWithSpaces > sStringBufWithSpaces) {
		pIn->putback(*--sBufWithSpaces);
	}

	return -1;
}


integer
HighParser::GetInt(int iDefval)
{
   const char sFuncName[] = "HighParser::GetInt()";

   if (CurrToken != HighParser::ARG) {
	silent_cerr("Parser error in "
	<< sFuncName << ", integer arg expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrIntegerExpected());
   }

   integer iReturnValue;

#ifdef USE_EXCEPTIONS
   try {
#endif

      iReturnValue = int(MathP.Get(*pIn, (doublereal)iDefval));

#ifdef USE_EXCEPTIONS
   }
   catch (MathParser::ErrGeneric e) {
      silent_cerr(sFuncName << ": error return from MathParser at line "
	      << GetLineData() << std::endl);
      throw e;
   }
#endif

   NextToken(sFuncName);
   return iReturnValue;
}


doublereal
HighParser::GetReal(const doublereal& dDefval)
{
   const char sFuncName[] = "HighParser::GetReal()";

   if (CurrToken != HighParser::ARG) {
      silent_cerr("Parser error in "
	<< sFuncName << ", real arg expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrRealExpected());
   }

   doublereal dReturnValue;
#ifdef USE_EXCEPTIONS
   try {
#endif

      dReturnValue = doublereal(MathP.Get(*pIn, dDefval));

#ifdef USE_EXCEPTIONS
   }
   catch (MathParser::ErrGeneric e) {
      silent_cerr(sFuncName << ": error return from MathParser at line "
	      << GetLineData() << std::endl);
      throw e;
   }
#endif

   NextToken(sFuncName);
   return dReturnValue;
}


int
HighParser::GetWord(void)
{
   const char sFuncName[] = "HighParser::GetWord()";

   if (CurrToken != HighParser::ARG) {
      silent_cerr("Parser error in "
	<< sFuncName << ", keyword arg expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrKeyWordExpected());
   }

   if ((CurrLowToken = LowP.GetToken(*pIn)) != LowParser::WORD) {
      silent_cerr("Parser error in "
	<< sFuncName << ", keyword expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrKeyWordExpected());
   }

   int i = -1;
   if (KeyT) {
      i = KeyT->Find(LowP.sGetWord());
   }

   NextToken(sFuncName);
   return i;
}


const char*
HighParser::GetString(void)
{
   const char sFuncName[] = "HighParser::GetString()";

   silent_cout("line " << GetLineData()
     << ": warning, use of deprecated method \"GetString\"" << std::endl);

   if (CurrToken != HighParser::ARG) {
      silent_cerr("Parser error in "
	<< sFuncName << ", string arg expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrStringExpected());
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
	 *sTmp++ = cIn;
      }
   }

   pIn->putback(cIn);
   *sTmp = '\0';

   NextToken(sFuncName);
   return s;
}


const char*
HighParser::GetStringWithDelims(enum Delims Del)
{
   const char sFuncName[] = "HighParser::GetStringWithDelims()";

   if (CurrToken != HighParser::ARG) {
      silent_cerr("Parser error in "
	<< sFuncName << ", string arg expected at line "
	<< GetLineData() << std::endl);
      THROW(HighParser::ErrStringExpected());
   }

   char* s = sStringBuf;
   char* sTmp = s;

   char cLdelim = '\0';
   char cRdelim = '\0';

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
	    sTmp[0] = '\0';
	    return s;
	 } else if (sTmp < s + iDefaultBufSize - 1) {
	    if (cIn == ESCAPE_CHAR) {
#if 0
	       /* FIXME: the escape char must be eaten 
	        * NOTE: the only charthat is worth escaping 
	        * is the right delimiter, I guess */
	       sTmp[0] = cIn;
	       ++sTmp;
#endif
	       cIn = pIn->get();
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
	<< GetLineData() << " isn't a valid left-delimiter;" << std::endl);
      THROW(HighParser::ErrIllegalDelimiter());
   }

   /* Mette zero al termine della stringa */
   *sTmp = '\0';

nullstring:

   NextToken(sFuncName);
   return s;
}


/* Legge un Vec3 */
Vec3
HighParser::GetVec3(void)
{
   Vec3 v(0.);
   return GetVec3(v);
}


/* Legge un Vec3 */
Vec3
HighParser::GetVec3(const Vec3& vDef)
{
   if (IsKeyWord("default")) {
      return vDef;
   }

   if (IsKeyWord("null")) {
      return Zero3;
   }

   doublereal x1 = GetReal(vDef.dGet(1));
   doublereal x2 = GetReal(vDef.dGet(2));
   doublereal x3 = GetReal(vDef.dGet(3));

   Vec3 v(x1, x2, x3);

   if (IsKeyWord("scale")) {
       v *= GetReal(1.);
   }

   return v;
}


/* Legge una matrice R sotto forma di due vettori (oppure eye) */
Mat3x3
HighParser::GetMatR2vec(void)
{
   if (IsKeyWord("eye")) {
      return Eye3;
   }

   if (IsKeyWord("matr")) {
      doublereal r11 = GetReal();
      doublereal r12 = GetReal();
      doublereal r13 = GetReal();
      doublereal r21 = GetReal();
      doublereal r22 = GetReal();
      doublereal r23 = GetReal();
      doublereal r31 = GetReal();
      doublereal r32 = GetReal();
      doublereal r33 = GetReal();

      /* FIXME: check for orthogonality? */
      return Mat3x3(r11, r21, r31, r12, r22, r32, r13, r23, r33);
   }

   if (IsKeyWord("euler" "parameters")) {
#if 0 /* FIXME: this function is TODO */
      doublereal e0 = GetReal();
      doublereal e1 = GetReal();
      doublereal e2 = GetReal();
      doublereal e3 = GetReal();

      return EulerParams2MatR(Vec3(e1, e2, e3));
#else
      silent_cerr("Line " << GetLineData()
	      << ": Euler parameters not allowed yet" << std::endl);
      THROW(ErrGeneric());
#endif
   }

   if (IsKeyWord("euler")) {
      doublereal e1 = GetReal();
      doublereal e2 = GetReal();
      doublereal e3 = GetReal();

      return EulerAngles2MatR(Vec3(e1, e2, e3));
   }

   int i1 = GetInt();
   doublereal x1 = GetReal();
   doublereal x2 = GetReal();
   doublereal x3 = GetReal();
   Vec3 v1(x1, x2, x3);

   int i2 = GetInt();
   x1 = GetReal();
   x2 = GetReal();
   x3 = GetReal();
   Vec3 v2(x1, x2, x3);

   return MatR2vec(i1, v1, i2, v2);
}


/* Legge una matrice 3x3 simmetrica come diagonale o triangolare superiore */
Mat3x3
HighParser::GetMat3x3Sym(void)
{
   if (IsKeyWord("null")) {
      return Mat3x3(0.);
   }
   
   Mat3x3 m;

   if (IsKeyWord("eye")) {
      m = Eye3;

   } else if (IsKeyWord("diag")) {
      doublereal m11 = GetReal();
      doublereal m22 = GetReal();
      doublereal m33 = GetReal();
      m = Mat3x3(m11, 0., 0., 0., m22, 0., 0., 0., m33);

   } else {   
      doublereal m11 = GetReal();
      doublereal m12 = GetReal();
      doublereal m13 = GetReal();
      doublereal m22 = GetReal();
      doublereal m23 = GetReal();
      doublereal m33 = GetReal();
      m = Mat3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);
   }

   if (IsKeyWord("scale")) {
       m *= GetReal(1.);
   }

   return m;
}

/* Legge una matrice 3x3 generica (diagonale o nulla) */
Mat3x3
HighParser::GetMat3x3(void)
{
   Mat3x3 m(0.);
   return GetMat3x3(m);
}


/* Legge una matrice 3x3 generica (diagonale o nulla) */
Mat3x3
HighParser::GetMat3x3(const Mat3x3& mDef)
{
   if (IsKeyWord("default")) {
      return mDef;

   } else if (IsKeyWord("null")) {
      return Zero3x3;
   }

   Mat3x3 m;

   if (IsKeyWord("eye")) {
      m = Eye3;

   } else if (IsKeyWord("diag")) {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      m = Mat3x3(m11, 0., 0., 0., m22, 0., 0., 0., m33);

   } else if (IsKeyWord("sym")) {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      m = Mat3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);

   } else {
      if (IsKeyWord("matr")) {
	      /* eat it; not required */
	      NO_OP;
      }
   
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      m = Mat3x3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
   }

   if (IsKeyWord("scale")) {
       m *= GetReal(1.);
   }

   return m;
}


/* Legge un Vec6 */
Vec6
HighParser::GetVec6(void)
{
   Vec6 v(0.);
   return GetVec6(v);
}


/* Legge un Vec6 */
Vec6
HighParser::GetVec6(const Vec6& vDef)
{
   if (IsKeyWord("null")) {
      return Zero6;
   }

   if (IsKeyWord("default")) {
      return vDef;
   }

   doublereal x1 = GetReal(vDef.dGet(1));
   doublereal x2 = GetReal(vDef.dGet(2));
   doublereal x3 = GetReal(vDef.dGet(3));
   doublereal x4 = GetReal(vDef.dGet(4));
   doublereal x5 = GetReal(vDef.dGet(5));
   doublereal x6 = GetReal(vDef.dGet(6));
   Vec6 v(x1, x2, x3, x4, x5, x6);

   if (IsKeyWord("scale")) {
       v *= GetReal(1.);
   }

   return v;
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6
HighParser::GetMat6x6(void)
{
   Mat6x6 m(0.);
   return GetMat6x6(m);
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6
HighParser::GetMat6x6(const Mat6x6& mDef)
{
   if (IsKeyWord("null")) {
      return Zero6x6;

   } else if (IsKeyWord("default")) {
      return mDef;
   }
   
   Mat6x6 m;

   if (IsKeyWord("eye")) {
      m = Eye6;

   } else if (IsKeyWord("diag")) {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      m = Mat6x6(m11, 0., 0., 0., 0., 0.,
		 0., m22, 0., 0., 0., 0.,
		 0., 0., m33, 0., 0., 0.,
		 0., 0., 0., m44, 0., 0.,
		 0., 0., 0., 0., m55, 0.,
		 0., 0., 0., 0., 0., m66);
 
   } else if (IsKeyWord("sym")) {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      m = Mat6x6(m11, m12, m13, m14, m15, m16,
		 m12, m22, m23, m24, m25, m26,
		 m13, m23, m33, m34, m35, m36,
		 m14, m24, m34, m44, m45, m46,
		 m15, m25, m35, m45, m55, m56,
		 m16, m26, m36, m46, m56, m66);

   } else if (IsKeyWord("anba")) {
      /* Formato ANBA, in cui vale la trasformazione:
       * ex = e2
       * ey = e3
       * ez = e1
       */
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m52 = GetReal(mDef.dGet(5, 2));
      doublereal m53 = GetReal(mDef.dGet(5, 3));
      doublereal m51 = GetReal(mDef.dGet(5, 1));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m54 = GetReal(mDef.dGet(5, 4));
      doublereal m62 = GetReal(mDef.dGet(6, 2));
      doublereal m63 = GetReal(mDef.dGet(6, 3));
      doublereal m61 = GetReal(mDef.dGet(6, 1));
      doublereal m65 = GetReal(mDef.dGet(6, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      doublereal m64 = GetReal(mDef.dGet(6, 4));
      doublereal m42 = GetReal(mDef.dGet(4, 2));
      doublereal m43 = GetReal(mDef.dGet(4, 3));
      doublereal m41 = GetReal(mDef.dGet(4, 1));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m44 = GetReal(mDef.dGet(4, 4));

      m = Mat6x6(m11, m21, m31, m41, m51, m61,
		 m12, m22, m32, m42, m52, m62,
		 m13, m23, m33, m43, m53, m63,
		 m14, m24, m34, m44, m54, m64,
		 m15, m25, m35, m45, m55, m65,
		 m16, m26, m36, m46, m56, m66);

   } else {
      if (IsKeyWord("matr")) {
	      /* eat it; not required */
	      NO_OP;
      }

      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m41 = GetReal(mDef.dGet(4, 1));
      doublereal m42 = GetReal(mDef.dGet(4, 2));
      doublereal m43 = GetReal(mDef.dGet(4, 3));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m51 = GetReal(mDef.dGet(5, 1));
      doublereal m52 = GetReal(mDef.dGet(5, 2));
      doublereal m53 = GetReal(mDef.dGet(5, 3));
      doublereal m54 = GetReal(mDef.dGet(5, 4));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m61 = GetReal(mDef.dGet(6, 1));
      doublereal m62 = GetReal(mDef.dGet(6, 2));
      doublereal m63 = GetReal(mDef.dGet(6, 3));
      doublereal m64 = GetReal(mDef.dGet(6, 4));
      doublereal m65 = GetReal(mDef.dGet(6, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      m = Mat6x6(m11, m21, m31, m41, m51, m61,
	         m12, m22, m32, m42, m52, m62,
	         m13, m23, m33, m43, m53, m63,
	         m14, m24, m34, m44, m54, m64,
	         m15, m25, m35, m45, m55, m65,
	         m16, m26, m36, m46, m56, m66);
   }

   if (IsKeyWord("scale")) {
       m *= GetReal(1.);
   }

   return m;
}

/* provvisoria */
void
HighParser::GetMat6xN(Mat3xN& m1, Mat3xN& m2, integer iNumCols)
{
   ASSERT(iNumCols > 0);
   ASSERT(m1.iGetNumCols() == iNumCols);
   ASSERT(m2.iGetNumCols() == iNumCols);

   if (IsKeyWord("null")) {
      m1.Reset();
      m2.Reset();

   } else {
      int vi[] = { 1, 2, 3 };

      if (IsKeyWord("anba")) {
	 vi[0] = 2;
	 vi[1] = 3;
	 vi[2] = 1;

      } else if (IsKeyWord("matr")) {
	      /* eat it; not required */
	      NO_OP;
      }

      for (int i = 0; i < 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m1.Put(vi[i], j, GetReal());
	 }
      }

      for (int i = 0; i < 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m2.Put(vi[i], j, GetReal());
	 }
      }

      if (IsKeyWord("scale")) {
         doublereal d = GetReal(1.);

	 m1 *= d;
	 m2 *= d;
      }
   }
}


/* HighParser - end */

std::ostream&
operator << (std::ostream& out, const HighParser::ErrOut& err)
{
   out << err.iLineNumber;
   if (err.sFileName != NULL) {
      out << ", file <" << err.sFileName << '>';
   }
   return out;
}
