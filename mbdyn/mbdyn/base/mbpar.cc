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
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "mbpar.h"

#if defined(USE_HYDRAULIC_NODES)
#include "hfluid.h"
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_AERODYNAMIC_ELEMS)
#include "aerodc81.h"
#include "c81data.h"
#endif /* USE_AERODYNAMIC_ELEMS */

#ifdef HAVE_LTDL_H
#include <ltdl.h>
#elif defined(HAVE_DLFCN_H)
#include <dlfcn.h>
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

#include "dataman.h"

/* MBDynParser - begin */

void
mbdyn_license(std::ostream& out)
{
	out << "license not available yet;"
		" see GPL" << std::endl;
}

void
mbdyn_warranty(std::ostream& out)
{
	out << "warranty not available yet;"
		" see warranty coming with GPL"
		<< std::endl;
}

MBDynParser::MBDynParser(MathParser& MP, 
		InputStream& streamIn,
		const char *initial_file)
: IncludeParser(MP, streamIn, initial_file),
#if defined(USE_STRUCT_NODES)
RFHD(),
RF(RFHD),
#endif /* USE_STRUCT_NODES */
#if defined(USE_HYDRAULIC_NODES)
HFHD(),
HF(HFHD),
#endif /* USE_HYDRAULIC_NODES */
#if defined(USE_AERODYNAMIC_ELEMS)
ADHD(),
AD(ADHD),
#endif /* USE_AERODYNAMIC_ELEMS */
C1DHD(),
C1D(C1DHD),
C3DHD(),
C3D(C3DHD),
C6DHD(),
C6D(C6DHD),
DCHD(),
DC(DCHD),
pDM(0)
{
	NO_OP;
}   


MBDynParser::~MBDynParser(void)
{   
	NO_OP;
}

void
MBDynParser::SetDataManager(DataManager *pdm)
{
	ASSERT(pdm != NULL);
	pDM = pdm;
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
MBDynParser::Reference_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::Reference_int(),"
			" colon expected at line " 
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected();
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del reference */
	const char *sName = NULL;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		SAFESTRDUP(sName, sTmp);
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Reference frame " << uLabel << std::endl);
	
	Vec3 x(GetPosAbs(AbsRefFrame));
	Mat3x3 R(GetRotAbs(AbsRefFrame));
	Vec3 v(0.);
	Vec3 w(0.);
	if (IsArg()) {
		v = GetVelAbs(AbsRefFrame, x);
		w = GetOmeAbs(AbsRefFrame);
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, std::endl
		   << "\tX = " << x << std::endl
		   << "\tR = " << R << std::endl
		   << "\tV = " << v << std::endl
		   << "\tW = " << w << std::endl);
	
	ReferenceFrame* pRF = NULL;
	SAFENEWWITHCONSTRUCTOR(pRF,
		ReferenceFrame,
		ReferenceFrame(uLabel, x, R, v, w));
	if (RF.Add(pRF)) {
		silent_cerr("Reference frame " << uLabel
			<< " already defined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrReferenceAlreadyDefined();
	}
	
	if (sName != NULL) {
		pRF->PutName(sName);
		SAFEDELETEARR(sName);
	}
}
#endif /* USE_STRUCT_NODES */

#if defined(USE_HYDRAULIC_NODES)
void 
MBDynParser::HydraulicFluid_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::HydraulicFluid_int(),"
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected();
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del fluido */
	const char *sName = NULL;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		SAFESTRDUP(sName, sTmp);
	}
	
	HydraulicFluid* pHF = ReadHydraulicFluid(*this, uLabel);
	if (pHF == NULL) {
		silent_cerr("unable to read hydraulic fluid " << uLabel
				<< std::endl);
		throw ErrGeneric();
	}
	
	if (HF.Add(pHF)) {
		silent_cerr("hydraulic fluid " << uLabel
			<< " already defined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric();
	}
	
	if (sName != NULL) {
		pHF->PutName(sName);
		SAFEDELETEARR(sName);
	}
}
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_AERODYNAMIC_ELEMS)
void 
MBDynParser::C81Data_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::C81Data_int(),"
			" colon expected at line " << GetLineData()
			<< std::endl);
		throw HighParser::ErrColonExpected();
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del profilo c81 */
	const char *sName = NULL;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		SAFESTRDUP(sName, sTmp);
	}
	
	const char* filename = GetFileName();
	std::ifstream in(filename);
	if (!in) {
		silent_cerr("unable to open file '" << filename << "' at line " 
			<< GetLineData() << std::endl);
		throw ErrGeneric();
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "reading c81 data " << uLabel 
		   << " from file '" << filename << "'" << std::endl);
	
	C81Data* data = NULL;
	SAFENEWWITHCONSTRUCTOR(data, C81Data, C81Data(uLabel));
	
	if (read_c81_data(in, data) != 0) {
		silent_cerr("unable to read c81 data " << uLabel 
			<< " from file '" << filename 
			<< "' at line " << GetLineData() << std::endl);
		throw ErrGeneric();
	}
	
#ifdef DEBUG
	if (DEBUG_LEVEL(MYDEBUG_INPUT)) {
		write_c81_data(std::cout, data);
	}
#endif

	if (AD.Add(data)) {
		silent_cerr("c81 data " << uLabel
			<< " already defined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric();
	}
	
	if (sName != NULL) {
		data->PutName(sName);
		SAFEDELETEARR(sName);
	}
}
#endif /* USE_AERODYNAMIC_ELEMS */

void 
MBDynParser::ConstitutiveLaw_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::ConstitutiveLaw_int(),"
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected();
	}

	if (pDM == 0) {
		silent_cerr("constitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric();
	}
	
	unsigned int uLabel(GetInt());
	
	/* Constitutive law name */
	const char *sName = NULL;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		SAFESTRDUP(sName, sTmp);
	}

	int dim = GetInt();
	ConstLawType::Type clt;
	switch (dim) {
	case 1:
	{
		ConstitutiveLaw1D *pCL = pDM->ReadConstLaw1D(*this, clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 1D " 
					<< uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric();
		}

		pCL->PutLabel(uLabel);
		if (sName) {
			pCL->PutName(sName);
		}
	
		if (C1D.Add(pCL)) {
			silent_cerr("constitutive law 1D " << uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric();
		}
		break;
	}

	case 3:
	{
		ConstitutiveLaw3D *pCL = pDM->ReadConstLaw3D(*this, clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 1D " 
					<< uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric();
		}
	
		pCL->PutLabel(uLabel);
		if (sName) {
			pCL->PutName(sName);
		}
	
		if (C3D.Add(pCL)) {
			silent_cerr("constitutive law 3D " << uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric();
		}
		break;
	}

	case 6:
	{
		ConstitutiveLaw6D *pCL = pDM->ReadConstLaw6D(*this, clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 6D " 
					<< uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric();
		}
	
		pCL->PutLabel(uLabel);
		if (sName) {
			pCL->PutName(sName);
		}
	
		if (C6D.Add(pCL)) {
			silent_cerr("constitutive law 6D " << uLabel);
			if (sName) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric();
		}
		break;
	}

	default:
		silent_cerr("unknown constitutive law dimensionality " 
				<< dim << " at line " << GetLineData()
				<< std::endl);
		throw ErrGeneric();
	}

	if (sName) {
		SAFEDELETEARR(sName);
	}
}

void 
MBDynParser::DriveCaller_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::DriveCaller_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected();
	}

	unsigned int uLabel(GetInt());
	
	/* drive name */
	const char *sName = NULL;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		SAFESTRDUP(sName, sTmp);
	}

	DriveCaller *pDC = ReadDriveData(pDM, *this);
	if (pDC == NULL) {
		silent_cerr("unable to read drive caller " << uLabel);
		if (sName) {
			silent_cerr(" (" << sName << ")");
		}
		silent_cerr(" at line " << GetLineData()
				<< std::endl);
		throw MBDynParser::ErrGeneric();
	}

	pDC->PutLabel(uLabel);
	if (sName) {
		pDC->PutName(sName);
	}
	
	if (DC.Add(pDC)) {
		silent_cerr("drive caller " << uLabel);
		if (sName) {
			silent_cerr(" (" << sName << ")");
		}
		silent_cerr(" already defined at line " 
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	if (sName) {
		SAFEDELETEARR(sName);
	}
}

void 
MBDynParser::ModuleLoad_int(void)
{
#if !defined(HAVE_RUNTIME_LOADING)
	silent_cerr("ModuleLoad_int: dynamic modules not supported"
			<< std::endl);
	throw ErrGeneric();
#else /* HAVE_RUNTIME_LOADING */
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::ModuleLoad_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected();
	}

   	/* nome del modulo */
   	const char* s = GetFileName();
	if (s == NULL) {
		silent_cerr("ModuleLoad_int: unable to get module name"
			<< std::endl);
		throw ErrGeneric();
	}

	char *module_name = 0;
   	SAFESTRDUP(module_name, s);

#ifdef HAVE_LTDL_H
	lt_dlhandle handle = lt_dlopenext(module_name);
#elif defined(HAVE_DLFCN_H)
   	void* handle = dlopen(module_name, RTLD_NOW /* RTLD_LAZY */ );
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

	if (handle == NULL) {
#ifdef HAVE_LTDL_H
      		const char* err = lt_dlerror();
#elif defined(HAVE_DLFCN_H)
      		const char* err = dlerror();
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

      		silent_cerr("ModuleLoad_int: "
			<< "unable to open module <" << module_name 
			<< "> (" << err << ") at line " << GetLineData()
			<< std::endl);
      		throw ErrGeneric();
   	}

	typedef int (*sym_t)(const char *, void *, void *);
#ifdef HAVE_LTDL_H
	sym_t sym = (sym_t)lt_dlsym(handle, "module_init");
#elif defined(HAVE_DLFCN_H)
	sym_t sym = (sym_t)dlsym(handle, "module_init");
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */
	
   	if (sym == NULL) {
#ifdef HAVE_LTDL_H
      		const char* err = lt_dlerror();
#elif defined(HAVE_DLFCN_H)
      		const char* err = dlerror();
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

      		if (err == NULL) {
			silent_cerr("ModuleLoad_int: module_init() "
				"function not available "
				"in module <" << module_name
				<< "> at line " << GetLineData()
				<< std::endl);
      		} else {
	 		silent_cerr("ModuleLoad_int: "
	   			<< "error while binding to symbol "
				"module_init() in module <" << module_name
	   			<< "> (" << err << ") at line " << GetLineData()
				<< std::endl);
      		}
      		throw ErrGeneric();
   	}

	if ((*sym)(module_name, (void *)pDM, (void *)this)) {
		silent_cerr("ModuleLoad_int: module_init() "
				"of module <" << module_name
				<< "> failed at line " << GetLineData()
				<< std::endl);
		throw ErrGeneric();
	}

   	SAFEDELETEARR(module_name);
#endif /* HAVE_RUNTIME_LOADING */
}

bool
MBDynParser::GetDescription_int(const char *s)
{
	/* Se trova un remark, scrive il commento ed eventualmente
	 * quello che segue */
	if (!strcmp(s, "remark")) {
		Remark_int();
		return true;

	/* Se trova un sistema di riferimento, lo gestisce direttamente */
	} else if (!strcmp(s, "reference")) {
#if defined(USE_STRUCT_NODES)      
		Reference_int();
		return true;
#else /* USE_STRUCT_NODES */
		throw MBDynParser::ErrGeneric();
#endif /* USE_STRUCT_NODES */

	/* Se trova un fluido idraulico, lo gestisce direttamente */
	} else if (!strcmp(s, "hydraulic" "fluid")) {
#if defined(USE_HYDRAULIC_NODES)
		HydraulicFluid_int();
		return true;
#else /* USE_HYDRAULIC_NODES */
		throw MBDynParser::ErrGeneric();
#endif /* USE_HYDRAULIC_NODES */

	/* Se trova dati aerodinamici c81, li gestisce direttamente */
	} else if (!strcmp(s, "c81" "data")) {
#if defined(USE_AERODYNAMIC_ELEMS)
		C81Data_int();
		return true;
#else /* USE_AERODYNAMIC_ELEMS */
		throw MBDynParser::ErrGeneric();
#endif /* USE_AERODYNAMIC_ELEMS */

	/* Reads a constitutive law */
	} else if (!strcmp(s, "constitutive" "law")) {
		ConstitutiveLaw_int();
		return true;

	/* Reads a drive caller */
	} else if (!strcmp(s, "drive" "caller")) {
		DriveCaller_int();
		return true;

	/* Loads a dynamic module */
	} else if (!strcmp(s, "module" "load" )) {
		ModuleLoad_int();
		return true;

	/* Scrive la licenza */
	} else if (!strcmp(s, "license")) {
		mbdyn_license(std::cout);
		CurrLowToken = LowP.GetToken(*pIn);
		return true;

	/* Scrive il disclaimer */
	} else if (!strcmp(s, "warranty")) {
		mbdyn_warranty(std::cout);
		CurrLowToken = LowP.GetToken(*pIn);
		return true;
	}

	/* altrimenti e' una description normale */
	return IncludeParser::GetDescription_int(s);
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
		silent_cerr("reference " << uLabel << " is undefined at line " 
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrReferenceUndefined();
	}
	
	rf = *pRF;
	return MBDynParser::REFERENCE;   
}

void 
MBDynParser::OutputFrames(std::ostream& out) const
{
	ReferenceFrame *pRF = NULL;

	if (!RF.GetFirst(pRF)) {
		return;	
	}
	do {
		pRF->Output(out);
	} while (RF.GetNext(pRF));
}

Vec3 
MBDynParser::GetPosRel(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return (rf.GetR()).Transpose()*(GetVec3()-rf.GetX());
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return rf.GetR().Transpose()*((rfOut.GetX()
			+rfOut.GetR()*GetVec3()
			)-rf.GetX());
			
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case UNKNOWNFRAME:
		return GetVec3();
	
	case NODE:
		return rf.GetX()+rf.GetR()*GetVec3();
		
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return rf.GetX()+rf.GetR()*(R*GetVec3());
	}
	
	case REFERENCE:
		return rfOut.GetX()+rfOut.GetR()*GetVec3();
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case GLOBAL:
		return (rf.GetR()).Transpose()*(GetVec3()
			-rf.GetV()
			-rf.GetW().Cross(x-rf.GetX()));
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return (rf.GetR()).Transpose()*(
			rfOut.GetV()
			+rfOut.GetR()*GetVec3()
			+rfOut.GetW().Cross(x-rfOut.GetX())
			-rf.GetV()
			-rf.GetW().Cross(x-rf.GetX())
			);
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case UNKNOWNFRAME:
		return GetVec3();
	
	case NODE:
		return rf.GetV()+rf.GetR()*GetVec3()
			+rf.GetW().Cross(x-rf.GetX());
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return rf.GetV()+rf.GetR()*(R*GetVec3())
			+rf.GetW().Cross(x-rf.GetX());	  
	}
	
	case REFERENCE:
		return rfOut.GetV()+rfOut.GetR()*GetVec3()
			+rfOut.GetW().Cross(x-rfOut.GetX());
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case GLOBAL:
		return (rf.GetR()).Transpose()*(GetVec3()-rf.GetW());
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return (rf.GetR()).Transpose()*((rfOut.GetW()
			+rfOut.GetR()*GetVec3()
			)-rf.GetW());
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case UNKNOWNFRAME:
		return GetVec3();
	
	case NODE:
		return rf.GetW()+rf.GetR()*GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return rf.GetW()+rf.GetR()*(R*GetVec3());
	}
	
	case REFERENCE:
		return rfOut.GetW()+rfOut.GetR()*GetVec3();
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case GLOBAL:
		return (rf.GetR()).Transpose()*GetVec3();
	
	case UNKNOWNFRAME: /* global */
		if (IsKeyWord("fromnode")) {
			/* FIXME */
			silent_cerr("'from node' at line " << GetLineData()
				<< " not implemented yet :)" << std::endl);
			throw MBDynParser::ErrGeneric();
			
			unsigned int uLabel = GetInt();
			StructNode *pNode1 = NULL; /* get node 1 */
			if (IsKeyWord("tonode")) {
				silent_cerr("missing keyword 'to node' at line "
					<< GetLineData() << std::endl);
				throw MBDynParser::ErrGeneric();
			}
			uLabel = GetInt();
			StructNode *pNode2 = NULL; /* get node 2 */

			Vec3 v = pNode2->GetXCurr()-pNode1->GetXCurr();
			return (rf.GetR()).Transpose()*v;
		} /* else local */
	case NODE:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return (rf.GetR()).Transpose()*(rfOut.GetR()*GetVec3());
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case UNKNOWNFRAME: /* global */
		if (IsKeyWord("fromnode")) {
			/* FIXME */
			silent_cerr("'from node' at line " << GetLineData()
				<< " not implemented yet :)" << std::endl);
			throw MBDynParser::ErrGeneric();
			
			unsigned int uLabel = GetInt();
			StructNode *pNode1 = NULL; /* get node 1 */
			if (IsKeyWord("tonode")) {
				silent_cerr("missing keyword 'to node' at line "
					<< GetLineData() << std::endl);
				throw MBDynParser::ErrGeneric();
			}
			uLabel = GetInt();
			StructNode *pNode2 = NULL; /* get node 2 */

			return pNode2->GetXCurr()-pNode1->GetXCurr();
		} /* else global */
	case GLOBAL:
		return GetVec3();
	
	case NODE:
		return rf.GetR()*GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return rf.GetR()*(R*GetVec3());
	}
	
	case REFERENCE:
		return rfOut.GetR()*GetVec3();
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
	}
	
# ifndef USE_EXCEPTIONS
	return Zero3;
# endif /* USE_EXCEPTIONS */
}

Mat3x3
MBDynParser::GetMatRel(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().Transpose()*(GetMat3x3()*rf.GetR());
		
	case NODE:
	case LOCAL:
	case UNKNOWNFRAME:
		return GetMat3x3();
	
	case REFERENCE:
		return rf.GetR().Transpose()
			*(rfOut.GetR()
			*(GetMat3x3()
			*(rfOut.GetR().Transpose()*rf.GetR())));
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
	}
	
# ifndef USE_EXCEPTIONS
	return Zero3x3;
# endif /* USE_EXCEPTIONS */
}

Mat3x3
MBDynParser::GetMatAbs(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
	case UNKNOWNFRAME:
		return GetMat3x3();
	
	case NODE:
	case LOCAL:
		return rf.GetR()*(GetMat3x3()*rf.GetR().Transpose());

	case REFERENCE:
		return rfOut.GetR()*(GetMat3x3()*rfOut.GetR().Transpose());

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
	}
	
# ifndef USE_EXCEPTIONS
	return Zero3x3;
# endif /* USE_EXCEPTIONS */
}

Mat3x3 
MBDynParser::GetRotRel(const ReferenceFrame& rf)
{   
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().Transpose()*GetMatR2vec();
	
	case NODE:
	case LOCAL:
	case UNKNOWNFRAME:
		return GetMatR2vec();
	
	case REFERENCE:
		return rf.GetR().Transpose()*(rfOut.GetR()*GetMatR2vec());
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
	case UNKNOWNFRAME:
		return GetMatR2vec();
	
	case NODE:
	case LOCAL:
		return rf.GetR()*GetMatR2vec();
	
	case REFERENCE:
		return rfOut.GetR()*GetMatR2vec();
	
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric();
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
		silent_cerr("hydraulic fluid expected at line "
			<< GetLineData() << std::endl);
		throw ErrGeneric();
	}
	
	/* se non trova "reference", legge direttamente un fluido */
	if (!IsKeyWord("reference")) {
		return ReadHydraulicFluid(*this, 0);
	}
	
	/* altrimenti usa un fluido predefinito, se lo trova */
	unsigned int uLabel = GetInt();
	const HydraulicFluid* pHF = HF.Get(uLabel);
	if (pHF == NULL) {
		silent_cerr("hydraulic fluid " << uLabel
			<< " is undefined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric();
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
		silent_cerr("c81 data " << profile << " is undefined at line " 
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}
	return data;
}
#endif /* USE_AERODYNAMIC_ELEM */

ConstitutiveLaw1D *
MBDynParser::GetConstLaw1D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	if (!IsKeyWord("reference")) {
		return pDM->ReadConstLaw1D(*this, clt);
	}

	unsigned int uLabel = GetInt();
	const ConstitutiveLaw1D *pCL = C1D.Get(uLabel);
	if (pCL == NULL) {
		silent_cerr("constitutive law 1D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	clt = pCL->GetConstLawType();
	return pCL->pCopy();
}

ConstitutiveLaw3D *
MBDynParser::GetConstLaw3D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	if (!IsKeyWord("reference")) {
		return pDM->ReadConstLaw3D(*this, clt);
	}

	unsigned int uLabel = GetInt();
	const ConstitutiveLaw3D *pCL = C3D.Get(uLabel);
	if (pCL == NULL) {
		silent_cerr("constitutive law 3D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	clt = pCL->GetConstLawType();
	return pCL->pCopy();
}

ConstitutiveLaw6D *
MBDynParser::GetConstLaw6D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	if (!IsKeyWord("reference")) {
		return pDM->ReadConstLaw6D(*this, clt);
	}

	unsigned int uLabel = GetInt();
	const ConstitutiveLaw6D *pCL = C6D.Get(uLabel);
	if (pCL == NULL) {
		silent_cerr("constitutive law 6D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	clt = pCL->GetConstLawType();
	return pCL->pCopy();
}

DriveCaller *
MBDynParser::GetDriveCaller(void)
{
	if (!IsKeyWord("reference")) {
		return ReadDriveData(pDM, *this);
	}

	unsigned int uLabel = GetInt();
	const DriveCaller *pDC = DC.Get(uLabel);
	if (pDC == NULL) {
		silent_cerr("drive caller " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric();
	}

	return pDC->pCopy();
}

/* MBDynParser - end */

