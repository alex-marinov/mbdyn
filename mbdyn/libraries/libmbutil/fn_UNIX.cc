/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/* Codice relativo alla classe File_name, per l'handling di nomi di files
 * in ambiente DOS e UNIX */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "filename.h"
#include "myassert.h"
#include "mynewmem.h"


FileName::FileName(const char *sFName, int iExtSepNum)
: sName(NULL), sExt(NULL), sRef(NULL)
{ 
   	if (sFName != NULL) { 
      		iInit(sFName, iExtSepNum); 
   	} 
}

FileName::~FileName(void)
{
   	if (sName != NULL) { 
      		/* delete []sName; */
      		SAFEDELETEARR(sName);
   	} 

   	if (sExt != NULL) { 
      		/* delete []sExt; */
      		SAFEDELETEARR(sExt);
   	} 
}

int 
FileName::iInit(const char *sFName, int iExtSepNum)
{
   	ASSERT(sFName != NULL);

   	if (sFName == NULL) { 
      		return 0; 
   	}
   
   	unsigned int iNewSize = strlen(sFName);
   	if ((sName == NULL) || (iNewSize > iMaxSize)) {
      		if (sName != NULL) {      
	 		SAFEDELETEARR(sName);
			sName = NULL;
      		}
      		iMaxSize = iNewSize;
      		SAFENEWARR(sName, char, iMaxSize + 1);
   	}
   
   	strcpy(sName, sFName); 
   	sRef = (char *)sName + strlen(sName);
   	ASSERT(sRef[0] == '\0');
   
   	if (iExtSepNum > 0) {  /* se si parte dall'inizio */     
      		while (--sRef > sName) {  /* cerca il simbolo "/" */	  
	 		if (sRef[0] == DIR_SEP) {	    
	    			break;  
	 		}
      		}
	
      		int iCnt = 0;
      		while (sRef[0] != '\0') {  /* cerca il "." n. iExtSepNum */
	 		if (sRef[0] == EXT_SEP) {
	    			iCnt++; 
	    			if (iCnt == iExtSepNum) {	       
	       				goto label;
	    			}
	 		}
	 		sRef++;
      		}
	
label:
      		iNewSize = strlen(sRef);
      		if (sExt == NULL || strlen(sExt) < iNewSize) {  
			/* se c'e' gia' sExt lo cancella */
	 		if (sExt != NULL) {	    
	    			SAFEDELETEARR(sExt);
	 		}
	 		SAFENEWARR(sExt, char, iNewSize+1);
      		}

      		ASSERT(sRef != NULL);
      		strcpy(sExt, sRef);
      		sRef[0] = '\0';
   	} else if (iExtSepNum < 0) {  /* se si parte dalla fine */   
      		iExtSepNum = -iExtSepNum;
      		while (--sRef > sName) {  /* cerca il simbolo "/" */	 
	 		if (sRef[0] == DIR_SEP) {	    
	    			break;
	 		}
      		}
	
      		int iCnt = 0;
      		char *sLim = sRef;
      		sRef = (char *)sName + strlen(sName);
      		while (--sRef > sLim) {  /* cerca il "." n. iExtSepNum */	  
	 		if (sRef[0] == EXT_SEP) {	    	       
	    			iCnt++; 
	    			if (iCnt == iExtSepNum) {	       
	       				goto label2; 
	    			}
	 		}
      		}
      		sRef = (char *)sName + strlen(sName);
	
label2:
      		iNewSize = strlen(sRef);
      		if (sExt == NULL || strlen(sExt) < iNewSize) { 
	 		if (sExt != NULL) {	    
	    			SAFEDELETEARR(sExt);
	 		}
	 		SAFENEWARR(sExt, char, iNewSize + 1);
      		}
	
      		ASSERT(sRef != NULL);
      		strcpy(sExt, sRef);
      		sRef[0] = '\0';
   	} else { 
		if (sExt == 0) {
      			SAFENEWARR(sExt, char, 1);
		}
      		sExt[0] = '\0'; 
   	}
   
   	iCurSize = iMaxSize - strlen(sExt);
   	return iCurSize;
}

const char *const
FileName::_sPutExt(const char *sEName)
{
   	if (sEName == NULL) {
      		sEName = sExt; 
   	}
 
   	unsigned int uExtLen = strlen(sEName);
   	if (sEName[0] != '\0' && sEName[0] != EXT_SEP) {      
     		uExtLen++;
   	}
   
   	if (iCurSize + uExtLen > iMaxSize) {
      		char *sTmp = NULL;
      		SAFENEWARR(sTmp, char, iCurSize + uExtLen + 1);
      
      		ASSERT(sName != NULL);
      		strcpy(sTmp, sName);
      		sRef = sTmp + (sRef - sName);
      		SAFEDELETEARR(sName);
      		sName = sTmp;
		iMaxSize = iCurSize + uExtLen;
   	}
   
   	ASSERT(sRef != NULL);
   	ASSERT(sEName != NULL);

   	if (sEName[0] != '\0') {
      		char *sTmp = sRef;
      		if (sEName[0] != EXT_SEP) {
	 		sTmp[0] = EXT_SEP;
			sTmp++;
      		} 
      		strcpy(sTmp, sEName);
   	}
   
   	return sName;
}

const char *const
FileName::sGet(void) const
{ 
   	return ((FileName*)this)->_sPutExt(NULL); 
}

int
is_abs_path(const char *const p)
{
	ASSERT(p != 0);

#ifdef _WIN32
	if ((strstr(p, ":\\") != 0) || (strstr(p, ":/") != 0)) {
		return 1;
	}
#else // ! _WIN32
	if (p[0] == '/') {
		return 1;
	}
#endif // ! _WIN32

	return 0;
}

