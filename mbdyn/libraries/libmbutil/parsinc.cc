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

/* parser con include */

#include <mbconfig.h>

#include <parsinc.h>

#include <unistd.h>
#include <pwd.h>

const int PATHBUFSIZE = 256;

/* IncludeParser - begin */

IncludeParser::IncludeParser(MathParser& MP, 
			     KeyTable& KT, 
			     InputStream& streamIn)
: HighParser(MP, KT, streamIn), sCurrPath(NULL), sCurrFile(NULL)
{   
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
   char* s = NULL;
   SAFENEWARR(s, char, PATHBUFSIZE, MPmm);
   sCurrPath = getcwd(s, PATHBUFSIZE);
   if (sCurrPath == NULL) {
      cerr << "Error in getcwd()" << endl;
      SAFEDELETEARR(s, MPmm);
      THROW(ErrFileSystem());
   }
   DEBUGCOUT("Current directory is <" << sCurrPath << '>' << endl);
   
   const char sInitialFile[] = "initial file";
   SAFENEWARR(sCurrFile, char, strlen(sInitialFile)+1, MPmm);
   strcpy(sCurrFile, sInitialFile); 
#else /* !defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
   NO_OP;
#endif /* !defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
}   


IncludeParser::~IncludeParser(void)
{   
   IncludeParser::Close();
}
 

void IncludeParser::Close(void)
{
   MyInput* pmi = NULL;
   if (MyInStack.iPop(pmi)) {
      ASSERT(pmi != NULL);
      /* Nota: deve esserci solo l'ultimo file */
      ASSERT(pf != NULL);
      ASSERT(pIn != NULL);
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      ASSERT(sCurrPath != NULL);
      ASSERT(sCurrFile != NULL);
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
      if (pf != NULL) {
	 SAFEDELETE(pf, MPmm);
      }
      if(pIn != NULL) {
	 SAFEDELETE(pIn, MPmm);
      }
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      DEBUGCOUT("Leaving directory <" << sCurrPath 
		<< ">, file <" << sCurrFile << '>' << endl);
      if (sCurrPath != NULL) {
	 SAFEDELETEARR(sCurrPath, MPmm);
	 sCurrPath = NULL;
      }
      if (sCurrFile != NULL) {
	 SAFEDELETEARR(sCurrFile, MPmm);
	 sCurrFile = NULL;
      }
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
      
      pf = pmi->pfile;
      pIn = pmi->pis;
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      sCurrPath = pmi->sPath;
      sCurrFile = pmi->sFile;
      DEBUGCOUT("Entering directory <" << sCurrPath 
		<< ">, file <" << sCurrFile << '>' << endl);
      if (chdir(sCurrPath)) {
	 cerr << "Error in chdir, path = " << sCurrPath << endl;
	 THROW(ErrFileSystem());
      };
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */

      /* pmi must be non NULL */
      SAFEDELETE(pmi, MPmm);     
   }
   
   /* sCurrPath can be NULL if Close() has been already called */
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
   if (sCurrPath != NULL) {
      SAFEDELETEARR(sCurrPath, MPmm);
      sCurrPath = NULL;
   }   
   if (sCurrFile != NULL) {
      SAFEDELETEARR(sCurrFile, MPmm);
      sCurrFile = NULL;
   }   
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
}


flag IncludeParser::fCheckStack(void)
{
   MyInput* pmi = NULL;
   if (MyInStack.iPop(pmi)) {
      ASSERT(pmi != NULL);
      /* 
       * Nota: se la stack e' piena, allora sia pf che pIn devono essere
       * diversi da NULL; viceversa, se la stack e' vuota, pf deve essere NULL 
       */
      ASSERT(pf != NULL);
      ASSERT(pIn != NULL);
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      ASSERT(sCurrPath != NULL);
      ASSERT(sCurrFile != NULL);
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
      
      SAFEDELETE(pf, MPmm); 
      SAFEDELETE(pIn, MPmm);
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      DEBUGCOUT("Leaving directory <" << sCurrPath 
		<< ">, file <" << sCurrFile << '>' << endl);
      SAFEDELETEARR(sCurrPath, MPmm);
      sCurrPath = NULL;
      SAFEDELETEARR(sCurrFile, MPmm);
      sCurrFile = NULL;
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
      
      pf = pmi->pfile;
      pIn = pmi->pis;
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
      sCurrPath = pmi->sPath;
      sCurrFile = pmi->sFile;
      DEBUGCOUT("Entering directory <" << sCurrPath 
		<< ">, file <" << sCurrFile << '>' << endl);
      if (chdir(sCurrPath)) {
	 cerr << "Error in chdir, path = " << sCurrPath << endl;
	 THROW(ErrFileSystem());
      };
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
      
      SAFEDELETE(pmi, MPmm);
      return flag(1);
   } else {
      return flag(0);
   }
}


void IncludeParser::Include_()
{
   if (FirstToken() == UNKNOWN) {
      cerr << endl << "Parser error in IncludeParser::Include_(), colon expected at line "
	<< GetLineData() << endl;
      THROW(HighParser::ErrColonExpected());
   }
   
   const char* sfname = this->GetFileName();
   
   MyInput* pmi = NULL;
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
   SAFENEWWITHCONSTRUCTOR(pmi, MyInput, MyInput(pf, pIn, sCurrPath, sCurrFile), MPmm);
#else /* !defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
   SAFENEWWITHCONSTRUCTOR(pmi, MyInput, MyInput(pf, pIn), MPmm);
#endif /* !defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
   MyInStack.Push(pmi);
   
   pf = NULL;
   SAFENEWWITHCONSTRUCTOR(pf, ifstream, ifstream(sfname), MPmm);
   if (!(*pf)) {
      cerr << "Invalid file <" << sfname << '>' << endl;
      THROW(ErrFile());
   }
   
   pIn = NULL;
   SAFENEWWITHCONSTRUCTOR(pIn, InputStream, InputStream(*pf), MPmm);
   
   /* Cambio di directory */
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
   sCurrPath = NULL;
   sCurrFile = NULL;
   int il = strlen(sfname);
   char* stmp = NULL;
   SAFENEWARR(stmp, char, il+1, MPmm);
   strcpy(stmp, sfname);
   char* s = (char*)stmp+il;
   while (--s >= stmp) {
      if (s[0] == '/') {
	 char c = *(s+1);
	 s[1] = '\0';
	 if (chdir(stmp)) {
	    cerr << "Error in chdir, path = " << stmp << endl;
	    THROW(ErrFileSystem());
	 };	 
	 char* p = NULL;
	 SAFENEWARR(p, char, PATHBUFSIZE, MPmm);
	 sCurrPath = getcwd(p, PATHBUFSIZE);
	 if (sCurrPath == NULL) {
	    cerr << "Error in getcwd()" << endl;
	    SAFEDELETEARR(s, MPmm);
	    THROW(ErrFileSystem());
	 }
	 DEBUGCOUT("Current directory is <" << sCurrPath << '>' << endl);
	 
	 s[1] = c;
	 break;
      }	 
   }
   s++;
   il = strlen(s);
   SAFENEWARR(sCurrFile, char, il+1, MPmm);
   strcpy(sCurrFile, s);   
   DEBUGCOUT("Opening file <" << sCurrFile << '>' << endl);
      
   SAFEDELETEARR(stmp, MPmm);
   
   if (sCurrPath == NULL) {
      char* s = NULL;
      
      SAFENEWARR(s, char, PATHBUFSIZE, MPmm);
      sCurrPath = getcwd(s, PATHBUFSIZE);
      if (sCurrPath == NULL) {
	 cerr << "Error in getcwd()" << endl;
	 SAFEDELETEARR(s, MPmm);
	 THROW(ErrFileSystem());
      }
      DEBUGCOUT("Current directory is <" << sCurrPath << '>' << endl);
   }
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */
   
   /* mettere un test se c'e' il punto e virgola? */
   CurrToken = HighParser::DESCRIPTION;
}


int IncludeParser::GetDescription(void)
{
   const char sFuncName[] = "IncludeParser::GetDescription()";

   /* Checks if current token is a description */
   if (!fIsDescription()) {
      THROW(HighParser::ErrInvalidCallToGetDescription());
   }
   
restart:
   
   if ((CurrLowToken = LowP.GetToken(*pIn)) != LowParser::WORD) {
      if (pIn->GetStream().eof()) {
	 if(fCheckStack()) {
	    /* Se la stack e' vuota lancia l'eccezione (errore) */
	    goto restart;
	 } else {
	    THROW(ErrFile());
	 }
      } else {     	 
	 cerr << endl << "Parser error in "
	   << sFuncName << ", keyword expected at line " 
	   << GetLineData() << endl;
	 THROW(HighParser::ErrKeyWordExpected());
      }      
   }
   
   /* Description corrente */
   char* s = LowP.sGetWord();
   
   /* Se trova la direttiva "include", la gestisce direttamente in modo
    * da aprire il nuovo file conservando quello corrente nella stack */
   if (!strcmp(s, "include")) {
      Include_();
      goto restart;      

      /* Se trova un sistema di riferimento, lo gestisce direttamente */
   } else if (!strcmp(s, "set")) {
      Set_();
      goto restart;
   } /* else */   
   return iGetDescription_(s);
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
	 int l;
	 
	 l = strlen(home);
	 SAFENEWARR(s, char, l+strlen(filename)+1, MPmm);
	 
	 strncpy(s, home, l);
	 strcpy(s+l, filename);
	 
	 return s;
	 
      } else {
	 char *p;
	 
	 p = strchr(filename, '/');
	 if (p == NULL) {
	    return NULL;
	 } 
	 
	 char *s = NULL;
	 int l = p-filename;
	 
	 SAFENEWARR(s, char, l+1, MPmm);
	 strncpy(s, filename, l);
	 s[l] = '\0';
	 
	 /* do passwd stuff */
	 struct passwd *pw;
	 
	 pw = getpwnam(s);
	 SAFEDELETEARR(s, MPmm);
	 
	 if (pw == NULL ) {
	    return NULL;
	 }
	 
	 l = strlen(pw->pw_dir);
	 s = NULL;
	 SAFENEWARR(s, char, l+strlen(p)+1, MPmm);
	 strncpy(s, pw->pw_dir, l);
	 strcpy(s+l, p);
	 
	 return s;
      }
   } else {
      return NULL;
   }
}


const char* IncludeParser::GetFileName(enum Delims Del)
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
      SAFEDELETEARR(stmp, MPmm);
   }
   return sStringBuf;
}


IncludeParser::ErrOut IncludeParser::GetLineData(void) const
{      
   ErrOut LineData;
   LineData.sFileName = sCurrFile;
   LineData.iLineNumber = GetLineNumber();
   return LineData;
}

/* IncludeParser - end */

