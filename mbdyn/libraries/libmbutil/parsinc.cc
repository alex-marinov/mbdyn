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

/* parser con include */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#if defined(HAVE_PWD_H)
#include <pwd.h>
#endif /* HAVE_PWD_H */

#include <parsinc.h>

const int PATHBUFSIZE = 1024;

/* IncludeParser - begin */

IncludeParser::IncludeParser(MathParser& MP, 
			     InputStream& streamIn,
			     const char *sInitialFile)
: HighParser(MP, streamIn),
sCurrPath(NULL),
sCurrFile(NULL)
{   
	ASSERT(sInitialFile != NULL);
#ifdef USE_INCLUDE_PARSER
   	char s[PATHBUFSIZE];
   	if (getcwd(s, PATHBUFSIZE) == NULL) {
		std::cerr << "Error in getcwd()" << std::endl;
      		THROW(ErrFileSystem());
   	}
	SAFESTRDUP(sCurrPath, s);
   	DEBUGCOUT("Current directory is <" << sCurrPath << '>' << std::endl);
   
   	SAFESTRDUP(sCurrFile, sInitialFile);
#else /* !USE_INCLUDE_PARSER */
   	NO_OP;
#endif /* !USE_INCLUDE_PARSER */
}


IncludeParser::~IncludeParser(void)
{   
   	IncludeParser::Close();
}
 

void IncludeParser::Close(void)
{
   	MyInput* pmi = NULL;
   	if (MyInStack.Pop(pmi)) {
      		ASSERT(pmi != NULL);
      		/* Nota: deve esserci solo l'ultimo file */
      		ASSERT(pf != NULL);
      		ASSERT(pIn != NULL);

#ifdef USE_INCLUDE_PARSER
      		ASSERT(sCurrPath != NULL);
      		ASSERT(sCurrFile != NULL);
#endif /* USE_INCLUDE_PARSER */

      		if (pf != NULL) {
	 		SAFEDELETE(pf);
      		}
      		if (pIn != NULL) {
	 		SAFEDELETE(pIn);
      		}

#ifdef USE_INCLUDE_PARSER
      		DEBUGCOUT("Leaving directory <" << sCurrPath 
			<< ">, file <" << sCurrFile << '>' << std::endl);
      		if (sCurrPath != NULL) {
	 		SAFEDELETEARR(sCurrPath);
	 		sCurrPath = NULL;
      		}
      		if (sCurrFile != NULL) {
	 		SAFEDELETEARR(sCurrFile);
	 		sCurrFile = NULL;
      		}
#endif /* USE_INCLUDE_PARSER */
      
      		pf = pmi->pfile;
      		pIn = pmi->pis;

#ifdef USE_INCLUDE_PARSER
      		sCurrPath = pmi->sPath;
      		sCurrFile = pmi->sFile;
      		DEBUGCOUT("Entering directory <" << sCurrPath 
			<< ">, file <" << sCurrFile << '>' << std::endl);
      		if (chdir(sCurrPath)) {
			std::cerr << "Error in chdir, path = " 
				<< sCurrPath << std::endl;
	 		THROW(ErrFileSystem());
      		}
#endif /* USE_INCLUDE_PARSER */

      		/* pmi must be non NULL */
      		SAFEDELETE(pmi);     
   	}
   
   	/* sCurrPath can be NULL if Close() has been already called */
#ifdef USE_INCLUDE_PARSER
   	if (sCurrPath != NULL) {
      		SAFEDELETEARR(sCurrPath);
      		sCurrPath = NULL;
   	}
   	if (sCurrFile != NULL) {
      		SAFEDELETEARR(sCurrFile);
      		sCurrFile = NULL;
   	}
#endif /* USE_INCLUDE_PARSER */
}

flag 
IncludeParser::fCheckStack(void)
{
   	MyInput* pmi = NULL;
   	if (MyInStack.Pop(pmi)) {
      		ASSERT(pmi != NULL);
      		/* 
       		 * Nota: se la stack e' piena, allora sia pf che pIn
		 * devono essere diversi da NULL; viceversa, se la stack
		 * e' vuota, pf deve essere NULL.
		 */
      		ASSERT(pf != NULL);
      		ASSERT(pIn != NULL);
#ifdef USE_INCLUDE_PARSER
      		ASSERT(sCurrPath != NULL);
      		ASSERT(sCurrFile != NULL);
#endif /* USE_INCLUDE_PARSER */
      
      		SAFEDELETE(pf); 
      		SAFEDELETE(pIn);
#ifdef USE_INCLUDE_PARSER
      		DEBUGCOUT("Leaving directory <" << sCurrPath 
			<< ">, file <" << sCurrFile << '>' << std::endl);
      		SAFEDELETEARR(sCurrPath);
      		sCurrPath = NULL;
      		SAFEDELETEARR(sCurrFile);
      		sCurrFile = NULL;
#endif /* USE_INCLUDE_PARSER */
      
      		pf = pmi->pfile;
      		pIn = pmi->pis;
#ifdef USE_INCLUDE_PARSER
      		sCurrPath = pmi->sPath;
      		sCurrFile = pmi->sFile;
      		DEBUGCOUT("Entering directory <" << sCurrPath 
			<< ">, file <" << sCurrFile << '>' << std::endl);
      		if (chdir(sCurrPath)) {
			std::cerr << "Error in chdir, path = " 
				<< sCurrPath << std::endl;
	 		THROW(ErrFileSystem());
      		}
#endif /* USE_INCLUDE_PARSER */
      
      		SAFEDELETE(pmi);
      		return flag(1);
   	} else {
      		return flag(0);
   	}
}

void 
IncludeParser::Include_()
{
   	if (FirstToken() == UNKNOWN) {
		std::cerr << "Parser error in IncludeParser::Include_(),"
			" colon expected at line " << GetLineData() 
			<< std::endl;
      		THROW(HighParser::ErrColonExpected());
   	}
   
	std::ifstream *pf_old = pf;
	InputStream *pIn_old = pIn;
	char *sOldPath = sCurrPath;
	char *sOldFile = sCurrFile;
   
   	const char* sfname = GetFileName();

   	pf = NULL;
   	pIn = NULL;

   	SAFENEWWITHCONSTRUCTOR(pf, std::ifstream, std::ifstream(sfname));
   	if (!(*pf)) {
#ifdef DEBUG
		char *buf = getcwd(NULL, 0);
		if (buf != NULL) {
			DEBUGCERR("Current directory <" << buf << ">" 
					<< std::endl);
			free(buf);
		}
#endif /* DEBUG */
		std::cerr << "Invalid file <" << sfname << '>' << std::endl;
      		THROW(ErrFile());
   	}
   
   	SAFENEWWITHCONSTRUCTOR(pIn, InputStream, InputStream(*pf));

   	/* Cambio di directory */
#ifdef USE_INCLUDE_PARSER
   	sCurrPath = NULL;
   	sCurrFile = NULL;
   	char* stmp = NULL;
   	SAFESTRDUP(stmp, sfname);
   	char* s = (char*)stmp+strlen(sfname);
   	while (--s >= stmp) {
      		if (s[0] == '/') {
	 		char c = s[1];
	 		s[1] = '\0';
	 		if (chdir(stmp)) {
				std::cerr << "Error in chdir, path = " 
					<< stmp << std::endl;
	    			THROW(ErrFileSystem());
	 		}
	 		char p[PATHBUFSIZE];
	 		if (getcwd(p, PATHBUFSIZE) == NULL) {
				std::cerr << "Error in getcwd()" << std::endl;
	    			SAFEDELETEARR(stmp);
	    			THROW(ErrFileSystem());
	 		}
			SAFESTRDUP(sCurrPath, p);
	 		DEBUGCOUT("Current directory is <" << sCurrPath 
				<< '>' << std::endl);
	 
	 		s[1] = c;
	 		break;
      		}
   	}
   	s++;
   	SAFESTRDUP(sCurrFile, s);
   	DEBUGCOUT("Opening file <" << sCurrFile << '>' << std::endl);
      
   	SAFEDELETEARR(stmp);
   
   	if (sCurrPath == NULL) {
      		char s[PATHBUFSIZE];
      		sCurrPath = getcwd(s, PATHBUFSIZE);
      		if (sCurrPath == NULL) {
			std::cerr << "Error in getcwd()" << std::endl;
	 		THROW(ErrFileSystem());
      		}
		SAFESTRDUP(sCurrPath, s);
      		DEBUGCOUT("Current directory is <" << sCurrPath 
			<< '>' << std::endl);
   	}
#endif /* USE_INCLUDE_PARSER */

      	MyInput* pmi = NULL;
#ifdef USE_INCLUDE_PARSER
   	SAFENEWWITHCONSTRUCTOR(pmi, 
			       MyInput,
			       MyInput(pf_old, pIn_old, sOldPath, sOldFile));
#else /* !USE_INCLUDE_PARSER */
   	SAFENEWWITHCONSTRUCTOR(pmi, MyInput, MyInput(pf_old, pIn_old));
#endif /* !USE_INCLUDE_PARSER */
 
   	MyInStack.Push(pmi);
 
   	/* FIXME: mettere un test se c'e' il punto e virgola? */
   	CurrToken = HighParser::DESCRIPTION;
}

void
IncludeParser::Eof(void)
{
	if(!fCheckStack()) {
		THROW(EndOfFile());
	}
}

bool
IncludeParser::GetDescription_int(const char *s)
{
   	if (!strcmp(s, "include")) {
      		Include_();
      		return true;

#ifdef USE_INCLUDE_PARSER
	} else if (!strcmp(s, "chdir")) {
  	 	if (FirstToken() == UNKNOWN) {
			std::cerr << "Parser error "
				"in IncludeParser::Include_(), "
				"colon expected at line " << GetLineData() 
				<< std::endl;
      			THROW(HighParser::ErrColonExpected());
   		}
   
   		const char* sfname = GetFileName();

      		if (chdir(sfname)) {
			std::cerr << "Error in chdir, path = " 
				<< sfname << std::endl;
	 		THROW(ErrFileSystem());
      		}
      		return true;

#endif /* USE_INCLUDE_PARSER*/
		
   	}
	
	return HighParser::GetDescription_int(s);
}

char *
resolve_filename(const char *filename)
{
   	if (filename[0] == '~') {
      		filename++;
      		if (filename[0] == '/') {
	 		/* do environment stuff */
	 		char *home;
	 
	 		home = getenv("HOME");
	 		if (home == NULL) {	 	 
	    			return NULL;
	 		}
	 
	 		char *s = NULL;
	 		int l, ll;
	 
	 		l = strlen(home);
			ll = l+strlen(filename)+1;
	 		SAFENEWARR(s, char, ll);
	 
	 		strncpy(s, home, l);
	 		strcpy(s+l, filename);
	 
	 		return s;

#if defined(HAVE_PWD_H)
      		} else {
	 		char *p;
	 
	 		p = strchr(filename, '/');
	 		if (p == NULL) {
	    			return NULL;
	 		} 
	 
	 		char *s = NULL;
	 		int l = p-filename;
	 
	 		SAFENEWARR(s, char, l+1);
	 		strncpy(s, filename, l);
	 		s[l] = '\0';

	 		/* do passwd stuff */
	 		struct passwd *pw;
	 
	 		pw = getpwnam(s);
	 		SAFEDELETEARR(s);
	 
	 		if (pw == NULL ) {
	    			return NULL;
	 		}
	 
	 		l = strlen(pw->pw_dir);
			int ll = l+strlen(p)+1;
	 		s = NULL;
	 		SAFENEWARR(s, char, ll);
	 		strncpy(s, pw->pw_dir, l);
	 		strcpy(s+l, p);
	 
	 		return s;
#endif /* HAVE_PWD_H */
      		}
   	} else {
      		return NULL;
   	}
}

const char* 
IncludeParser::GetFileName(enum Delims Del)
{
   	const char *s = GetStringWithDelims(Del);
   	const char *stmp = resolve_filename(s);
	
   	if (stmp == NULL) {
      		return s;
   	} else {
      		if (strlen(stmp) >= iBufSize) {
	 		/* errore */
	 		strncpy(sStringBuf, stmp, iBufSize-1);
	 		sStringBuf[iBufSize-1] = '\0';
      		} else {
	 		strcpy(sStringBuf, stmp);
      		}
      		SAFEDELETEARR(stmp);
   	}
	
   	return sStringBuf;
}

IncludeParser::ErrOut
IncludeParser::GetLineData(void) const
{      
   	ErrOut LineData;
   	LineData.sFileName = sCurrFile;
   	LineData.iLineNumber = GetLineNumber();
   	return LineData;
}

/* IncludeParser - end */

