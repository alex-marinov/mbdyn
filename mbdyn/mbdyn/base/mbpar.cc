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

/* parser */

#include <mbconfig.h>

#include <mbpar.h>

#if defined(USE_HYDRAULIC_NODES)
#include <hfluid.h>
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_AERODYNAMIC_ELEMS)
#include <aerodc81.h>
#include <c81data.h>
#endif /* USE_AERODYNAMIC_ELEMS */

/* MBDynParser - begin */

MBDynParser::MBDynParser(MathParser& MP, KeyTable& KT, InputStream& streamIn)
: IncludeParser(MP, KT, streamIn)
#if defined(USE_STRUCT_NODES)
#ifdef DEBUG_MEMMANAGER
  , RFHD(MPmm)
#else
  , RFHD()
#endif
, RF(RFHD)
#endif /* USE_STRUCT_NODES */
#if defined(USE_HYDRAULIC_NODES)
#ifdef DEBUG_MEMMANAGER
  , HFHD(MPmm)
#else
  , HFHD()
#endif
, HF(HFHD)
#endif /* USE_HYDRAULIC_NODES */
#if defined(USE_AERODYNAMIC_ELEMS)
#ifdef DEBUG_MEMMANAGER
  , ADHD(MPmm)
#else
  , ADHD()
#endif
, AD(ADHD)
#endif /* USE_AERODYNAMIC_ELEMS */
{
   NO_OP;
}   


MBDynParser::~MBDynParser(void)
{   
   NO_OP;
}

#if defined(USE_STRUCT_NODES)
const ReferenceFrame AbsRefFrame(0, 
				 Vec3(0.), 
				 Mat3x3(1., 0., 0., 
					0., 1., 0.,
					0., 0., 1.),
				 Vec3(0.),
				 Vec3(0.));


void 
MBDynParser::Reference_(void)
{
   if (FirstToken() == UNKNOWN) {
      cerr << endl << "Parser error in MBDynParser::Reference_(), colon expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrColonExpected());
   }
   
   unsigned int uLabel(GetInt());
   
   DEBUGLCOUT(MYDEBUG_INPUT, "Reference frame " << uLabel << endl);
   
   Vec3 x(GetPosAbs(AbsRefFrame));
   Mat3x3 R(GetRotAbs(AbsRefFrame));
   Vec3 v(0.);
   Vec3 w(0.);
   if (fIsArg()) {
      v = GetVelAbs(AbsRefFrame, x);
      w = GetOmeAbs(AbsRefFrame);
   }
   
   DEBUGLCOUT(MYDEBUG_INPUT, endl
	      << "\tX = " << x << endl
	      << "\tR = " << R << endl
	      << "\tV = " << v << endl
	      << "\tW = " << w << endl);
   
   ReferenceFrame* pRF = NULL;
   SAFENEWWITHCONSTRUCTOR(pRF, ReferenceFrame,
			  ReferenceFrame(uLabel, x, R, v, w), MPmm);
   if (RF.iAdd(pRF)) {
      cerr << "Reference frame " << uLabel << "already defined at line"
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrReferenceAlreadyDefined());
   }
}
#endif /* USE_STRUCT_NODES */

#if defined(USE_HYDRAULIC_NODES)
void 
MBDynParser::HydraulicFluid_(void)
{
   if (FirstToken() == UNKNOWN) {
      cerr << endl << "Parser error in MBDynParser::HydraulicFluid_(), colon expected at line "
	<< GetLineData() << endl;
      THROW(HighParser::ErrColonExpected());
   }
   
   unsigned int uLabel(GetInt());
   
   KeyTable CurrTable = KeyT;
   
   HydraulicFluid* pHF = ReadHydraulicFluid(*this, uLabel);
   if (pHF == NULL) {
      cerr << "unable to read hydraulic fluid " << uLabel << endl;
      THROW(ErrGeneric());
   }
   PutKeyTable(CurrTable);
   
   if (HF.iAdd(pHF)) {
      cerr << "hydraulic fluid " << uLabel << "already defined at line"
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrGeneric());
   }
}
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_AERODYNAMIC_ELEMS)
void 
MBDynParser::C81Data_(void)
{
   if (FirstToken() == UNKNOWN) {
      cerr << endl << "Parser error in MBDynParser::C81Data_(), colon expected at line "
	<< GetLineData() << endl;
      THROW(HighParser::ErrColonExpected());
   }
   
   unsigned int uLabel(GetInt());

   const char* filename = GetFileName();
   ifstream in(filename);
   if (!in) {
      cerr << "unable to open file '" << filename << "' at line " 
	<< GetLineData() << endl;
      THROW(ErrGeneric());
   }
   
   DEBUGLCOUT(MYDEBUG_INPUT, "reading c81 data " << uLabel 
	      << " from file '" << filename << "'" << endl);
   
   C81Data* data = NULL;
   SAFENEWWITHCONSTRUCTOR(data, C81Data, C81Data(uLabel), DMmm);

   if (read_c81_data(in, data) != 0) {
      cerr << "unable to read c81 data " << uLabel 
	<< " from file '" << filename 
	<< "' at line " << GetLineData() << endl;
      THROW(ErrGeneric());
   }
   
# ifdef DEBUG
   if (DEBUG_LEVEL(MYDEBUG_INPUT)) {
      write_c81_data(cout, data);
   }
# endif

   if (AD.iAdd(data)) {
      cerr << "c81 data " << uLabel << "already defined at line"
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrGeneric());
   }
}
#endif /* USE_AERODYNAMIC_ELEMS */


int 
MBDynParser::GetDescription(void)
{
   const char sFuncName[] = "IncludeParser::GetDescription()";

   /* Checks if current token is a description */
   if (!fIsDescription()) {
      THROW(HighParser::ErrInvalidCallToGetDescription());
   }
   
restart:
   
   if ((CurrLowToken = LowP.GetToken(*pIn)) != LowParser::WORD) {
      if (pIn->GetStream().eof()) {
	 if (fCheckStack()) {
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
      
   } else if (!strcmp(s, "reference")) {
#if defined(USE_STRUCT_NODES)      
      Reference_();
      goto restart;
#else /* USE_STRUCT_NODES */
      THROW(MBDynParser::ErrGeneric());
#endif /* USE_STRUCT_NODES */
   } else if (!strcmp(s, "hydraulic" "fluid")) {
#if defined(USE_HYDRAULIC_NODES)
      HydraulicFluid_();
      goto restart;
#else /* USE_HYDRAULIC_NODES */
      THROW(MBDynParser::ErrGeneric());
#endif /* USE_HYDRAULIC_NODES */
   } else if (!strcmp(s, "c81" "data")) {
#if defined(USE_AERODYNAMIC_ELEMS)
      C81Data_();
      goto restart;
#else /* USE_AERODYNAMIC_ELEMS */
      THROW(MBDynParser::ErrGeneric());
#endif /* USE_AERODYNAMIC_ELEMS */
   } /* else */   
   return iGetDescription_(s);
}


#if defined(USE_STRUCT_NODES)
MBDynParser::Frame 
MBDynParser::GetRef(ReferenceFrame& rf)
{
   if (!IsKeyWord("reference")) {
      return MBDynParser::UNKNOWNFRAME;
   }
      
   if (IsKeyWord("global")) {
      return MBDynParser::GLOBAL;
   }
   
   if (IsKeyWord("node")) {
      return MBDynParser::NODE;
   }
   
   if (IsKeyWord("local")) {
      return MBDynParser::LOCAL;
   }

   unsigned int uLabel((unsigned int)GetInt());
   const ReferenceFrame* pRF = RF.Get(uLabel);
   if (pRF == NULL) {
      cerr << "reference " << uLabel << " is undefined at line " 
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrReferenceUndefined());
   }
   rf = *pRF;
   return MBDynParser::REFERENCE;   
}


Vec3 
MBDynParser::GetPosRel(const ReferenceFrame& rf)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL: {
       return (rf.GetR()).Transpose()*(GetVec3()-rf.GetX());
    }
      
    case NODE:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return R*GetVec3();
    }
      
    case REFERENCE: {	   
       return
	 rf.GetR().Transpose()*((rfOut.GetX()
				 +rfOut.GetR()*GetVec3()
				 )-rf.GetX());
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetPosAbs(const ReferenceFrame& rf)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case NODE: {
       return rf.GetX()+rf.GetR()*GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return rf.GetX()+rf.GetR()*(R*GetVec3());
    }
      
    case REFERENCE: {	   
       return rfOut.GetX()+rfOut.GetR()*GetVec3();
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetVelRel(const ReferenceFrame& rf, const Vec3& x)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL: {
       return (rf.GetR()).Transpose()*(GetVec3()
				       -rf.GetV()
				       -rf.GetW().Cross(x-rf.GetX()));
    }
      
    case NODE:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return R*GetVec3();
    }
      
    case REFERENCE: {	   
       return
	 (rf.GetR()).Transpose()*(
				  rfOut.GetV()
				  +rfOut.GetR()*GetVec3()
				  +rfOut.GetW().Cross(x-rfOut.GetX())
				  -rf.GetV()
				  -rf.GetW().Cross(x-rf.GetX())
				  );
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());       
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetVelAbs(const ReferenceFrame& rf, const Vec3& x)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case NODE: {
       return rf.GetV()+rf.GetR()*GetVec3()+rf.GetW().Cross(x-rf.GetX());
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return rf.GetV()+rf.GetR()*(R*GetVec3())
	 +rf.GetW().Cross(x-rf.GetX());	  
    }
      
    case REFERENCE: {	   
       return rfOut.GetV()+rfOut.GetR()*GetVec3()
	 +rfOut.GetW().Cross(x-rfOut.GetX());
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");      
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetOmeRel(const ReferenceFrame& rf)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL: {
       return (rf.GetR()).Transpose()*(GetVec3()-rf.GetW());
    }
      
    case NODE:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return R*GetVec3();
    }
      
    case REFERENCE: {	   
       return
	 (rf.GetR()).Transpose()*((rfOut.GetW()
				   +rfOut.GetR()*GetVec3()
				   )-rf.GetW());
	}
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");       
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetOmeAbs(const ReferenceFrame& rf)
{
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case NODE: {
       return rf.GetW()+rf.GetR()*GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return rf.GetW()+rf.GetR()*(R*GetVec3());
    }
      
    case REFERENCE: {	   
       return rfOut.GetW()+rfOut.GetR()*GetVec3();
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetVecRel(const ReferenceFrame& rf)
{   
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL: {
       return (rf.GetR()).Transpose()*GetVec3();
    }
      
    case NODE:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return R*GetVec3();
    }
      
    case REFERENCE: {	   
       return (rf.GetR()).Transpose()*(rfOut.GetR()*GetVec3());
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Vec3 
MBDynParser::GetVecAbs(const ReferenceFrame& rf)
{   
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL:
    case UNKNOWNFRAME: {
       return GetVec3();
    }
      
    case NODE: {
       return rf.GetR()*GetVec3();
    }
      
    case LOCAL: {
       Mat3x3 R(GetMatR2vec());
       return rf.GetR()*(R*GetVec3());
    }
      
    case REFERENCE: {	   
       return rfOut.GetR()*GetVec3();
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3;
# endif /* USE_EXCEPTIONS */
}


Mat3x3 
MBDynParser::GetRotRel(const ReferenceFrame& rf)
{   
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL: {
       return rf.GetR().Transpose()*GetMatR2vec();
    }
      
    case NODE:
    case LOCAL:
    case UNKNOWNFRAME: {
       return GetMatR2vec();
    }
      
    case REFERENCE: {
       return rf.GetR().Transpose()*(rfOut.GetR()*GetMatR2vec());
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3x3;
# endif /* USE_EXCEPTIONS */
}


Mat3x3 
MBDynParser::GetRotAbs(const ReferenceFrame& rf)
{   
   ReferenceFrame rfOut;
   switch (GetRef(rfOut)) {
    case GLOBAL:
    case UNKNOWNFRAME: {
       return GetMatR2vec();
    }
      
    case NODE:
    case LOCAL: {
       return rf.GetR()*GetMatR2vec();
    }
      
    case REFERENCE: {
       return rfOut.GetR()*GetMatR2vec();
    }
      
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(MBDynParser::ErrGeneric());
    }
   }
# ifndef USE_EXCEPTIONS
   return Zero3x3;
# endif /* USE_EXCEPTIONS */
}
#endif /* USE_STRUCT_NODES */


#if defined(USE_HYDRAULIC_NODES)
HydraulicFluid* 
MBDynParser::GetHydraulicFluid(void)
{
   /* verifica che sia stato chiamato con "hydraulic" "fluid" */
   if (!IsKeyWord("hydraulic" "fluid") && !IsKeyWord("fluid")) {
      cerr << "hydraulic fluid expected at line " << GetLineData() << endl;
      THROW(ErrGeneric());
   }
   
   /* se non trova "reference", legge direttamente un fluido */
   if (!IsKeyWord("reference")) {
      return ReadHydraulicFluid(*this, 0);
   }
   
   /* altrimenti usa un fluido predefinito, se lo trova */
   unsigned int uLabel = GetInt();
   const HydraulicFluid* pHF = HF.Get(uLabel);
   ASSERT(pHF != NULL);
   if (pHF == NULL) {
      cerr << "hydraulic fluid " << uLabel << " is undefined at line " 
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrGeneric());
   }
   return pHF->pCopy();
}
#endif /* USE_HYDRAULIC_NODES */


#if defined(USE_AERODYNAMIC_ELEMS)
const c81_data* 
MBDynParser::GetC81Data(integer profile)
{
   /* cerca i dati predefiniti, se li trova */
   const c81_data* data = AD.Get(profile);
   ASSERT(data != NULL);
   if (data == NULL) {
      cerr << "c81 data " << profile << " is undefined at line " 
	<< GetLineData() << endl;
      THROW(MBDynParser::ErrGeneric());
   }
   return data;
}
#endif /* USE_AERODYNAMIC_ELEM */

/* MBDynParser - end */
