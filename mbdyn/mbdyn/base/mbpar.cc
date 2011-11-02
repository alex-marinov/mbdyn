/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include <limits>
#include <cfloat>
#include <limits>
#include <typeinfo>

#if defined(USE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif // USE_RUNTIME_LOADING && HAVE_LTDL_H

#include "mbpar.h"
#include "Rot.hh"

#include "hfluid.h"

#include "aerodc81.h"
#include "c81data.h"

#include "dataman.h"
#include "modules.h"

/* MBDynParser - begin */

void
mbdyn_license(void)
{
	silent_cout("license not available yet;"
		" see GPL at http://www.gnu.org/" << std::endl);
}

void
mbdyn_warranty(void)
{
	silent_cout("From GPL 2.1:\n"
"\n"
"  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY\n"
"  FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN\n"
"  OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES\n"
"  PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED\n"
"  OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n"
"  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS\n"
"  TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE\n"
"  PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,\n"
"  REPAIR OR CORRECTION.\n"
"\n"
"  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING\n"
"  WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR\n"
"  REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,\n"
"  INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING\n"
"  OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED\n"
"  TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY\n"
"  YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER\n"
"  PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE\n"
"  POSSIBILITY OF SUCH DAMAGES.\n"
<< std::endl);
}

MBDynParser::MBDynParser(MathParser& MP, 
		InputStream& streamIn,
		const char *initial_file)
: IncludeParser(MP, streamIn, initial_file),
moduleInitialized(false),
pDM(0)
{
	/* make sure this is init'ed */
	InitDriveData();
	InitTplDC();
	InitCL();
	InitSF();
}   

MBDynParser::~MBDynParser(void)
{   
	DestroyDriveData();
	DestroyTplDC();
	DestroyCL();
	DestroySF();

	for (SFType::iterator i = SF.begin(); i != SF.end(); ++i) {
		SAFEDELETE(i->second);
	}
	SF.clear();

	for (RFType::iterator i = RF.begin(); i != RF.end(); ++i) {
		SAFEDELETE(i->second);
	}
	RF.clear();

	for (HFType::iterator i = HF.begin(); i != HF.end(); ++i) {
		SAFEDELETE(i->second);
	}
	HF.clear();

	for (ADType::iterator i = AD.begin(); i != AD.end(); ++i) {
		c81_data_destroy(i->second);
		SAFEDELETE(i->second);
	}
	AD.clear();

	for (CL1DType::iterator i = CL1D.begin(); i != CL1D.end(); ++i) {
		SAFEDELETE(i->second);
	}
	CL1D.clear();

	for (CL3DType::iterator i = CL3D.begin(); i != CL3D.end(); ++i) {
		SAFEDELETE(i->second);
	}
	CL3D.clear();

	for (CL6DType::iterator i = CL6D.begin(); i != CL6D.end(); ++i) {
		SAFEDELETE(i->second);
	}
	CL6D.clear();

	for (DCType::iterator i = DC.begin(); i != DC.end(); ++i) {
		SAFEDELETE(i->second);
	}
	DC.clear();

	if (!bEmptyManip()) {
		silent_cerr("MBDynParser::~MBDynParser: "
			"manipulators' stack not empty" << std::endl);
	}

#if defined(USE_RUNTIME_LOADING)
	if (moduleInitialized) {
		module_finalize();
	}
#endif // USE_RUNTIME_LOADING
}

void
MBDynParser::SetDataManager(DataManager *pdm)
{
	ASSERT(pdm != NULL);
	pDM = pdm;
	const DriveHandler *pDH = pDM->pGetDrvHdl();
	if (pDH == 0) {
		silent_cerr("no drive handler is associated to data manager?"
				<< std::endl);

	} else {
		/* add the drive handler to the drive callers... */
		for (DCType::const_iterator i = DC.begin(); i != DC.end(); ++i) {
			const_cast<DriveCaller *>(i->second)->SetDrvHdl(pDH);
		}
	}
}

const ReferenceFrame AbsRefFrame(0, Zero3, Eye3, Zero3, Zero3, EULER_123);


void 
MBDynParser::Reference_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::Reference_int(),"
			" colon expected at line " 
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del reference */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Reference frame " << uLabel << std::endl);
	
	if (!IsKeyWord("position")) {
		pedantic_cerr("ReferenceFrame(" << uLabel
			<< "): missing keyword \"position\" at line "
			<< GetLineData() << std::endl);
	}
	Vec3 x(GetPosAbs(AbsRefFrame));
	if (!IsKeyWord("orientation")) {
		pedantic_cerr("ReferenceFrame(" << uLabel
			<< "): missing keyword \"orientation\" at line "
			<< GetLineData() << std::endl);
	}
	Mat3x3 R(GetRotAbs(AbsRefFrame));
	Vec3 v(Zero3);
	Vec3 w(Zero3);
	if (IsArg()) {
		if (!IsKeyWord("velocity")) {
			pedantic_cerr("ReferenceFrame(" << uLabel
				<< "): missing keyword \"velocity\" at line "
				<< GetLineData() << std::endl);
		}
		v = GetVelAbs(AbsRefFrame, x);
		if (IsArg()) {
			if (!IsKeyWord("angular" "velocity")) {
				pedantic_cerr("ReferenceFrame(" << uLabel
					<< "): missing keyword \"angular velocity\" at line "
					<< GetLineData() << std::endl);
			}
			w = GetOmeAbs(AbsRefFrame);
		}
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after reference " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DEBUGLCOUT(MYDEBUG_INPUT, std::endl
		   << "\tX = " << x << std::endl
		   << "\tR = " << R << std::endl
		   << "\tV = " << v << std::endl
		   << "\tW = " << w << std::endl);
	
	ReferenceFrame* pRF = NULL;
	OrientationDescription od
		= ReadOptionalOrientationDescription(pDM, *this);
	SAFENEWWITHCONSTRUCTOR(pRF,
		ReferenceFrame,
		ReferenceFrame(uLabel, x, R, v, w, od));
	if (!RF.insert(RFType::value_type(uLabel, pRF)).second) {
		silent_cerr("Reference frame " << uLabel
			<< " already defined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrReferenceAlreadyDefined(MBDYN_EXCEPT_ARGS);
	}
	
	if (!sName.empty()) {
		pRF->PutName(sName);
	}
}

void 
MBDynParser::HydraulicFluid_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::HydraulicFluid_int(),"
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del fluido */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}
	
	HydraulicFluid* pHF = ReadHydraulicFluid(*this, uLabel);
	if (pHF == NULL) {
		silent_cerr("unable to read hydraulic fluid " << uLabel
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after hydraulic fluid " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HF.insert(HFType::value_type(uLabel, pHF)).second) {
		silent_cerr("hydraulic fluid " << uLabel
			<< " already defined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	if (!sName.empty()) {
		pHF->PutName(sName);
	}
}

void 
MBDynParser::C81Data_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::C81Data_int(),"
			" colon expected at line " << GetLineData()
			<< std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}
	
	unsigned int uLabel(GetInt());
	
	/* Nome del profilo c81 */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}
	
	std::string filename = GetFileName();
	if (filename.empty()) {
		silent_cerr("C81Data(" << uLabel << "): "
			"invalid file at line " 
			<< GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::ifstream in(filename.c_str());
	if (!in) {
		silent_cerr("C81Data(" << uLabel << "): "
			"unable to open file '" << filename << "' at line " 
			<< GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "reading c81 data " << uLabel 
		   << " from file '" << filename << "'" << std::endl);
	
	C81Data* data = NULL;
	SAFENEWWITHCONSTRUCTOR(data, C81Data, C81Data(uLabel));

	doublereal dcptol = 1e-6;
	if (IsKeyWord("tolerance")) {
		dcptol = GetReal();
		if (dcptol <= 0.) {
			silent_cerr("C81Data(" << uLabel << "): "
				"invalid c81 data tolerance at line "
				<< GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool bFF(false);
	if (IsKeyWord("fc511")) {
		if (c81_data_fc511_read(in, data, dcptol) != 0) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to read c81 data " << uLabel 
				<< " from file '" << filename << "' "
				"in fc511 format at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (IsKeyWord("nrel")) {
		if (c81_data_nrel_read(in, data, dcptol) != 0) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to read c81 data " << uLabel 
				<< " from file '" << filename << "' "
				"in NREL format at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (IsKeyWord("free" "format")) {
		if (c81_data_read_free_format(in, data, dcptol) != 0) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to read c81 data " << uLabel 
				<< " from file '" << filename << "' "
				"in free format at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		bFF = true;

	} else {
		int ff = 0;
		if (c81_data_read(in, data, dcptol, &ff) != 0) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to read c81 data " << uLabel 
				<< " from file '" << filename << "' "
				"at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (ff) {
			bFF = true;
		}
	}

	if (IsKeyWord("flip")) {
		(void)c81_data_flip(data);
	}

	// CL
	if (data->al[0] > -180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CL alpha[0]=" << data->al[0] << " > -180 (error=" << 100*(data->al[0]/180. + 1.) << "%)" << std::endl);
	}
	if (data->al[data->NAL - 1] < 180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CL alpha[" << data->NAL - 1 << "]=" << data->al[data->NAL - 1] << " < 180 (error=" << 100*(data->al[data->NAL - 1]/180. - 1.) << "%)" << std::endl);
	}
	if (data->ml[0] > 0.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CL mach[0]=" << data->ml[0] << " > 0" << std::endl);
	}

	// CD
	if (data->ad[0] > -180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CD alpha[0]=" << data->ad[0] << " > -180 (error=" << 100*(data->ad[0]/180. + 1.) << "%)" << std::endl);
	}
	if (data->ad[data->NAD - 1] < 180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CD alpha[" << data->NAD - 1 << "]=" << data->ad[data->NAD - 1] << " < 180 (error=" << 100*(data->ad[data->NAD - 1]/180. - 1.) << "%)" << std::endl);
	}
	if (data->md[0] > 0.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CD mach[0]=" << data->md[0] << " > 0" << std::endl);
	}

	// CM
	if (data->am[0] > -180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CM alpha[0]=" << data->am[0] << " > -180 (error=" << 100*(data->am[0]/180. + 1.) << "%)" << std::endl);
	}
	if (data->am[data->NAM - 1] < 180.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CM alpha[" << data->NAM - 1 << "]=" << data->am[data->NAM - 1] << " < 180 (error=" << 100*(data->am[data->NAM - 1]/180. - 1.) << "%)" << std::endl);
	}
	if (data->mm[0] > 0.) {
		silent_cerr("C81Data(" << uLabel << "): "
			"warning, CM mach[0]=" << data->mm[0] << " > 0" << std::endl);
	}

	if (IsKeyWord("echo")) {
		const char *sOutName = GetFileName();
		if (sOutName == NULL) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to read output file name "
				"at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::ofstream out(sOutName);
		if (!out) {
			silent_cerr("C81Data(" << uLabel << "): "
				"unable to open output file "
				"\"" << sOutName << "\" "
				"at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (IsKeyWord("free" "format")) {
			c81_data_write_free_format(out, data);

		} else if (!IsArg()) {
			if (bFF) {
				c81_data_write_free_format(out, data);

			} else {
				c81_data_write(out, data);
			}

		} else {
			silent_cerr("C81Data(" << uLabel << "): "
				"unknown output format "
				"at line " << GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after c81 data " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!AD.insert(ADType::value_type(uLabel, data)).second) {
		silent_cerr("C81Data(" << uLabel << "): "
			"redefined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	if (!sName.empty()) {
		data->PutName(sName.c_str());
	}
}

void 
MBDynParser::ConstitutiveLaw_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::ConstitutiveLaw_int(),"
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	if (pDM == 0) {
		silent_cerr("constitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	unsigned int uLabel(GetInt());
	
	/* Constitutive law name */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}

	int dim = GetInt();
	ConstLawType::Type clt;
	switch (dim) {
	case 1:
	{
		/* allow "reference" (copy cached constitutive law) */
		ConstitutiveLaw1D *pCL = GetConstLaw1D(clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 1D " 
					<< uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pCL->PutLabel(uLabel);
		if (!sName.empty()) {
			pCL->PutName(sName);
		}
	
		if (!CL1D.insert(CL1DType::value_type(uLabel, pCL)).second) {
			silent_cerr("constitutive law 1D " << uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		break;
	}

	case 3:
	{
		/* allow "reference" (copy cached constitutive law) */
		ConstitutiveLaw3D *pCL = GetConstLaw3D(clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 3D " 
					<< uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	
		pCL->PutLabel(uLabel);
		if (!sName.empty()) {
			pCL->PutName(sName);
		}
	
		if (!CL3D.insert(CL3DType::value_type(uLabel, pCL)).second) {
			silent_cerr("constitutive law 3D " << uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		break;
	}

	case 6:
	{
		/* allow "reference" (copy cached constitutive law) */
		ConstitutiveLaw6D *pCL = GetConstLaw6D(clt);
		if (pCL == NULL) {
			silent_cerr("unable to read constitutive law 6D " 
					<< uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" at line " << GetLineData()
					<< std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	
		pCL->PutLabel(uLabel);
		if (!sName.empty()) {
			pCL->PutName(sName);
		}
	
		if (!CL6D.insert(CL6DType::value_type(uLabel, pCL)).second) {
			silent_cerr("constitutive law 6D " << uLabel);
			if (!sName.empty()) {
				silent_cerr(" (" << sName << ")");
			}
			silent_cerr(" already defined at line " 
					<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		break;
	}

	default:
		silent_cerr("unknown constitutive law dimensionality " 
				<< dim << " at line " << GetLineData()
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after constitutive law " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void 
MBDynParser::DriveCaller_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::DriveCaller_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	unsigned int uLabel(GetInt());
	
	/* drive name */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}

	bool bDeferred(false);
	if (IsKeyWord("deferred")) {
		bDeferred = true;
	}

	/* allow "reference" (copy cached drive) */
	DriveCaller *pDC = GetDriveCaller(bDeferred);
	if (pDC == NULL) {
		silent_cerr("unable to read drive caller " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after drive caller " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pDC->PutLabel(uLabel);
	if (!sName.empty()) {
		pDC->PutName(sName);
	}
	
	if (!DC.insert(DCType::value_type(uLabel, pDC)).second) {
		silent_cerr("drive caller " << uLabel);
		if (!sName.empty()) {
			silent_cerr(" (" << sName << ")");
		}
		silent_cerr(" already defined at line " 
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void 
MBDynParser::TplDriveCaller_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::TplDriveCaller_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	unsigned int uLabel(GetInt());
	
	/* drive name */
	std::string sName;
	if (IsKeyWord("name")) {
		const char *sTmp = GetStringWithDelims();
		sName = sTmp;
	}

	int dim = GetInt();
	std::string type;
	switch (dim) {
	case 1: {
		type = "doublereal ";
		TplDriveCaller<doublereal> *pDC = GetTplDriveCaller<doublereal>();
		if (pDC == NULL) {
			silent_cerr("unable to read doublereal template drive caller " << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
		pDC->PutLabel(uLabel);
		if (!sName.empty()) {
			pDC->PutName(sName);
		}
#endif
	
		if (!DC1D.insert(DC1DType::value_type(uLabel, pDC)).second) {
			silent_cerr("doublereal template drive caller" << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"already defined at line " 
				<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		} break;

	case 3: {
		type = "Vec3 ";
		TplDriveCaller<Vec3> *pDC = GetTplDriveCaller<Vec3>();
		if (pDC == NULL) {
			silent_cerr("unable to read Vec3 template drive caller " << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
		pDC->PutLabel(uLabel);
		if (!sName.empty()) {
			pDC->PutName(sName);
		}
#endif
	
		if (!DC3D.insert(DC3DType::value_type(uLabel, pDC)).second) {
			silent_cerr("Vec3 template drive caller" << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"already defined at line " 
				<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		} break;

	case 6: {
		type = "Vec6 ";
		TplDriveCaller<Vec6> *pDC = GetTplDriveCaller<Vec6>();
		if (pDC == NULL) {
			silent_cerr("unable to read Vec6 template drive caller " << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
		pDC->PutLabel(uLabel);
		if (!sName.empty()) {
			pDC->PutName(sName);
		}
#endif
	
		if (!DC6D.insert(DC6DType::value_type(uLabel, pDC)).second) {
			silent_cerr("Vec6 template drive caller" << uLabel
				<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
				"already defined at line " 
				<< GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		} break;
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after " << type << "template drive caller " << uLabel
			<< " (" << (sName.empty() ? "unknown" : sName.c_str()) << ") "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void 
MBDynParser::ScalarFunction_int(void)
{
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::ScalarFunction_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

	(void)ParseScalarFunction(*this, pDM);

	if (IsArg()) {
		silent_cerr("semicolon expected after scalar function "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void 
MBDynParser::ModuleLoad_int(void)
{
#ifndef USE_RUNTIME_LOADING
	silent_cerr("ModuleLoad_int: dynamic modules not supported"
			<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);

#else // USE_RUNTIME_LOADING
	if (FirstToken() == UNKNOWN) {
		silent_cerr("Parser error in MBDynParser::ModuleLoad_int(), "
			" colon expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
	}

   	/* nome del modulo */
   	const char* s = GetFileName();
	if (s == NULL) {
		silent_cerr("ModuleLoad_int: unable to get module name"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::string module_name(s);

	if (!moduleInitialized) {
		module_initialize();
		moduleInitialized = true;
	}

	lt_dlhandle handle = lt_dlopenext(module_name.c_str());

	if (handle == NULL) {
      		const char* err = lt_dlerror();
		if (err == 0) {
			err = "";
		}

      		silent_cerr("ModuleLoad_int: "
			<< "unable to open module <" << module_name 
			<< "> (" << err << ") at line " << GetLineData()
			<< std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}

	typedef int (*sym_t)(const char *, void *, void *);
	sym_t sym = (sym_t)lt_dlsym(handle, "module_init");
   	if (sym == NULL) {
      		const char* err = lt_dlerror();
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
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}

	if ((*sym)(module_name.c_str(), (void *)pDM, (void *)this)) {
		silent_cerr("ModuleLoad_int: module_init() "
				"of module <" << module_name
				<< "> failed at line " << GetLineData()
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (IsArg()) {
		silent_cerr("semicolon expected after module "
			"\"" << module_name << "\""
			<< " at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	silent_cout("module \"" << module_name << "\" loaded" << std::endl);
#endif // USE_RUNTIME_LOADING
}

bool
MBDynParser::GetDescription_int(const char *s)
{
	/* Se trova un sistema di riferimento, lo gestisce direttamente */
	if (!strcmp(s, "reference")) {
		Reference_int();
		return true;

	/* Se trova un fluido idraulico, lo gestisce direttamente */
	} else if (!strcmp(s, "hydraulic" "fluid")) {
		HydraulicFluid_int();
		return true;

	/* Se trova dati aerodinamici c81, li gestisce direttamente */
	} else if (!strcmp(s, "c81" "data")) {
		C81Data_int();
		return true;

	/* Reads a constitutive law */
	} else if (!strcmp(s, "constitutive" "law")) {
		ConstitutiveLaw_int();
		return true;

	/* Reads a drive caller */
	} else if (!strcmp(s, "drive" "caller")) {
		DriveCaller_int();
		return true;

	/* Reads a template drive caller */
	} else if (!strcmp(s, "template" "drive" "caller")) {
		TplDriveCaller_int();
		return true;

	/* Reads a scalar function */
	} else if (!strcmp(s, "scalar" "function")) {
		ScalarFunction_int();
		return true;

	/* Loads a dynamic module */
	} else if (!strcmp(s, "module" "load" )) {
		ModuleLoad_int();
		return true;

	/* Scrive la licenza */
	} else if (!strcmp(s, "license")) {
		mbdyn_license();
		CurrLowToken = LowP.GetToken(*pIn);

		if (IsArg()) {
			silent_cerr("semicolon expected after \"license\" "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return true;

	/* Scrive il disclaimer */
	} else if (!strcmp(s, "warranty")) {
		mbdyn_warranty();
		CurrLowToken = LowP.GetToken(*pIn);

		if (IsArg()) {
			silent_cerr("semicolon expected after \"warranty\" "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return true;
	}

	/* altrimenti e' una description normale */
	return IncludeParser::GetDescription_int(s);
}

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

	if (IsKeyWord("other" "position")) {
		return MBDynParser::OTHER_POSITION;
	}
	
	if (IsKeyWord("other" "orientation")) {
		return MBDynParser::OTHER_ORIENTATION;
	}
	
	if (IsKeyWord("other" "node")) {
		return MBDynParser::OTHER_NODE;
	}
	
	unsigned int uLabel((unsigned int)GetInt());
	RFType::const_iterator i = RF.find(uLabel);
	if (i == RF.end()) {
		silent_cerr("reference " << uLabel << " is undefined at line " 
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrReferenceUndefined(MBDYN_EXCEPT_ARGS);
	}

	rf = *(i->second);

	return MBDynParser::REFERENCE;   
}

void 
MBDynParser::OutputFrames(std::ostream& out) const
{
	for (RFType::const_iterator i = RF.begin(); i != RF.end(); ++i) {
		i->second->Output(out);
	}
}

Vec3 
MBDynParser::GetPosRel(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTV(GetVec3() - rf.GetX());
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return rf.GetR().MulTV((rfOut.GetX() + rfOut.GetR()*GetVec3()) - rf.GetX());

	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetPosRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetPosRel(const ReferenceFrame& rf, const ReferenceFrame& other_rf, const Vec3& other_X)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTV(GetVec3() - rf.GetX());
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return rf.GetR().MulTV((rfOut.GetX() + rfOut.GetR()*GetVec3()) - rf.GetX());

	case OTHER_POSITION:
		return rf.GetR().MulTV((other_rf.GetX() + other_rf.GetR()*(other_X + GetVec3())) - rf.GetX());

	case OTHER_NODE:
		return rf.GetR().MulTV((other_rf.GetX() + other_rf.GetR()*GetVec3()) - rf.GetX());

	case OTHER_ORIENTATION:
		silent_cerr("GetPosRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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
		return rf.GetX() + rf.GetR()*GetVec3();
		
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return rf.GetX() + rf.GetR()*(R*GetVec3());
	}
	
	case REFERENCE:
		return rfOut.GetX() + rfOut.GetR()*GetVec3();
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetPosAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetVelRel(const ReferenceFrame& rf, const Vec3& x)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTV(GetVec3()
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
		return rf.GetR().MulTV(
			rfOut.GetV()
			+rfOut.GetR()*GetVec3()
			+rfOut.GetW().Cross(x-rfOut.GetX())
			-rf.GetV()
			-rf.GetW().Cross(x-rf.GetX())
			);
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetVelRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetVelAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetOmeRel(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTV(GetVec3()-rf.GetW());
	
	case NODE:
	case UNKNOWNFRAME:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return rf.GetR().MulTV(rfOut.GetW() + rfOut.GetR()*GetVec3() - rf.GetW());
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetOmeRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetOmeAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetVecRel(const ReferenceFrame& rf)
{   
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTV(GetVec3());
	
	case UNKNOWNFRAME: /* global */
		if (IsKeyWord("from" "node")) {
			/* FIXME */
			silent_cerr("'from node' at line " << GetLineData()
				<< " not implemented yet :)" << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

#if 0
			unsigned int uLabel = GetInt();
			StructNode *pNode1 = NULL; /* get node 1 */


			if (!IsKeyWord("to" "node")) {
				silent_cerr("missing keyword 'to node' at line "
					<< GetLineData() << std::endl);
				throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			uLabel = GetInt();
			StructNode *pNode2 = NULL; /* get node 2 */

			Vec3 v = pNode2->GetXCurr() - pNode1->GetXCurr();
			return rf.GetR().MulTV(v);
#endif // 0
		} /* else local */
	case NODE:
		return GetVec3();
	
	case LOCAL: {
		Mat3x3 R(GetMatR2vec());
		return R*GetVec3();
	}
	
	case REFERENCE:
		return rf.GetR().MulTV(rfOut.GetR()*GetVec3());
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetVecRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetVecAbs(const ReferenceFrame& rf)
{   
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case UNKNOWNFRAME: /* global */
		if (IsKeyWord("from" "node")) {
			/* FIXME */
			silent_cerr("'from node' at line " << GetLineData()
				<< " not implemented yet :)" << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

#if 0
			unsigned int uLabel = GetInt();
			StructNode *pNode1 = NULL; /* get node 1 */
			if (IsKeyWord("to" "node")) {
				silent_cerr("missing keyword 'to node' at line "
					<< GetLineData() << std::endl);
				throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			uLabel = GetInt();
			StructNode *pNode2 = NULL; /* get node 2 */

			return pNode2->GetXCurr() - pNode1->GetXCurr();
#endif // 0
		} /* else global: fallthru */

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
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetVecAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
MBDynParser::GetUnitVecRel(const ReferenceFrame& rf)
{
	Vec3 v = GetVecRel(rf);

	doublereal d = v.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	d = std::sqrt(d);

	if (std::fabs(d - 1.) > std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("warning: non-unit vector (norm=" << d << ") "
			"at line " << GetLineData() << "; "
			"normalized" << std::endl);
	}

	return v /= d;
}

Vec3 
MBDynParser::GetUnitVecAbs(const ReferenceFrame& rf)
{
	Vec3 v = GetVecAbs(rf);

	doublereal d = v.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	d = std::sqrt(d);

	if (std::fabs(d - 1.) > std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("warning: non-unit vector (norm=" << d << ") "
			"at line " << GetLineData() << "; "
			"normalized" << std::endl);
	}

	return v /= d;
}

Mat3x3
MBDynParser::GetMatRel(const ReferenceFrame& rf)
{
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTM(GetMat3x3()*rf.GetR());
		
	case NODE:
	case LOCAL:
	case UNKNOWNFRAME:
		return GetMat3x3();
	
	case REFERENCE:
		return rf.GetR().MulTM(rfOut.GetR()*GetMat3x3()*rfOut.GetR().MulTM(rf.GetR()));
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetMatRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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

	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetMatAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Mat3x3 
MBDynParser::GetRotRel(const ReferenceFrame& rf)
{   
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTM(GetMatR2vec());
	
	case NODE:
	case LOCAL:
	case UNKNOWNFRAME:
		return GetMatR2vec();
	
	case REFERENCE:
		return rf.GetR().MulTM(rfOut.GetR()*GetMatR2vec());
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetRotRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Mat3x3 
MBDynParser::GetRotRel(const ReferenceFrame& rf, const ReferenceFrame& other_rf, const Mat3x3& other_R)
{   
	ReferenceFrame rfOut;
	switch (GetRef(rfOut)) {
	case GLOBAL:
		return rf.GetR().MulTM(GetMatR2vec());
	
	case NODE:
	case LOCAL:
	case UNKNOWNFRAME:
		return GetMatR2vec();
	
	case REFERENCE:
		return rf.GetR().MulTM(rfOut.GetR()*GetMatR2vec());
	
	case OTHER_POSITION:
		silent_cerr("GetRotRel: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	case OTHER_ORIENTATION:
		return rf.GetR().MulTM(other_rf.GetR()*other_R*GetMatR2vec());

	case OTHER_NODE:
		return rf.GetR().MulTM(other_rf.GetR()*GetMatR2vec());

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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
	
	case OTHER_POSITION:
	case OTHER_ORIENTATION:
	case OTHER_NODE:
		silent_cerr("GetRotAbs: \"other\" meaningless in this context "
			"at line " << GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

HydraulicFluid* 
MBDynParser::GetHydraulicFluid(void)
{
	/* verifica che sia stato chiamato con "hydraulic" "fluid" */
	if (!IsKeyWord("hydraulic" "fluid") && !IsKeyWord("fluid")) {
		silent_cerr("hydraulic fluid expected at line "
			<< GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	/* se non trova "reference", legge direttamente un fluido */
	if (!IsKeyWord("reference")) {
		return ReadHydraulicFluid(*this, 0);
	}
	
	/* altrimenti usa un fluido predefinito, se lo trova */
	unsigned int uLabel = GetInt();
	HFType::const_iterator i = HF.find(uLabel);
	if (i == HF.end()) {
		silent_cerr("hydraulic fluid " << uLabel
			<< " is undefined at line " << GetLineData()
			<< std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return i->second->pCopy();
}

const c81_data* 
MBDynParser::GetC81Data(integer profile)
{
	/* cerca i dati predefiniti, se li trova */
	ADType::const_iterator i = AD.find(profile);
	if (i == AD.end()) {
		silent_cerr("c81 data " << profile << " is undefined at line " 
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (i == AD.end()) {
		silent_cerr("c81 data " << profile << " is undefined at line " 
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return i->second;
}

ConstitutiveLaw1D *
MBDynParser::GetConstLaw1D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!IsKeyWord("reference")) {
		return ReadCL1D(pDM, *this, clt);
	}

	unsigned int uLabel = GetInt();
	CL1DType::const_iterator i = CL1D.find(uLabel);
	if (i == CL1D.end()) {
		silent_cerr("constitutive law 1D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	clt = i->second->GetConstLawType();
	return i->second->pCopy();
}

ConstitutiveLaw3D *
MBDynParser::GetConstLaw3D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!IsKeyWord("reference")) {
		return ReadCL3D(pDM, *this, clt);
	}

	unsigned int uLabel = GetInt();
	CL3DType::const_iterator i = CL3D.find(uLabel);
	if (i == CL3D.end()) {
		silent_cerr("constitutive law 3D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	clt = i->second->GetConstLawType();
	return i->second->pCopy();
}

ConstitutiveLaw6D *
MBDynParser::GetConstLaw6D(ConstLawType::Type& clt)
{
	if (pDM == 0) {
		silent_cerr("consitutive law parsing at line "
				<< GetLineData() << " allowed "
				"only after control data block" << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!IsKeyWord("reference")) {
		return ReadCL6D(pDM, *this, clt);
	}

	unsigned int uLabel = GetInt();
	CL6DType::const_iterator i = CL6D.find(uLabel);
	if (i == CL6D.end()) {
		silent_cerr("constitutive law 6D " << uLabel
				<< " is undefined at line "
				<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	clt = i->second->GetConstLawType();
	return i->second->pCopy();
}

const DriveCaller *
MBDynParser::GetDrive(unsigned uLabel) const
{
	DCType::const_iterator i = DC.find(uLabel);
	if (i == DC.end()) {
		return 0;
	}

	return i->second;
}

DriveCaller *
MBDynParser::GetDriveCaller(bool bDeferred)
{
	if (!IsKeyWord("reference")) {
		DriveCaller *pDC = 0;
		try {
			pDC = ReadDriveData(pDM, *this, bDeferred);
		}
		catch (DataManager::ErrNeedDataManager) {
			silent_cerr("the required drive caller must appear "
					"inside or after the \"control data\" "
					"block"
					<< std::endl);
			throw DataManager::ErrNeedDataManager(MBDYN_EXCEPT_ARGS);
		}
		return pDC;
	}

	unsigned int uLabel = GetInt();
	const DriveCaller *pDC = GetDrive(uLabel);
	if (pDC == 0) {
		silent_cerr("drive caller " << uLabel
			<< " is undefined at line "
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pDC->pCopy();
}

template <class T>
const TplDriveCaller<T> *
MBDynParser::GetTplDrive(unsigned uLabel) const
{
	if (typeid(T) == typeid(doublereal)) {
		DC1DType::const_iterator i = DC1D.find(uLabel);
		if (i == DC1D.end()) {
			return 0;
		}

		return dynamic_cast<const TplDriveCaller<T> *>(i->second);

	} else if (typeid(T) == typeid(Vec3)) {
		DC3DType::const_iterator i = DC3D.find(uLabel);
		if (i == DC3D.end()) {
			return 0;
		}

		return dynamic_cast<const TplDriveCaller<T> *>(i->second);

	} else if (typeid(T) == typeid(Vec6)) {
		DC6DType::const_iterator i = DC6D.find(uLabel);
		if (i == DC6D.end()) {
			return 0;
		}

		return dynamic_cast<const TplDriveCaller<T> *>(i->second);
	}

	return 0;
}

template <class T>
TplDriveCaller<T> *
MBDynParser::GetTplDriveCaller(void)
{
	if (!IsKeyWord("reference")) {
		TplDriveCaller<T> *pDC = 0;
		try {
			if (typeid(T) == typeid(doublereal)) {
				pDC = dynamic_cast<TplDriveCaller<T> *>(ReadDC1D(pDM, *this));

			} else if (typeid(T) == typeid(Vec3)) {
				pDC = dynamic_cast<TplDriveCaller<T> *>(ReadDC3D(pDM, *this));

			} else if (typeid(T) == typeid(Vec6)) {
				pDC = dynamic_cast<TplDriveCaller<T> *>(ReadDC6D(pDM, *this));

			} else {
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		catch (DataManager::ErrNeedDataManager) {
			silent_cerr("the required drive caller must appear "
				"inside or after the \"control data\" block"
				<< std::endl);
			throw DataManager::ErrNeedDataManager(MBDYN_EXCEPT_ARGS);
		}
		return pDC;
	}

	unsigned int uLabel = GetInt();
	const TplDriveCaller<T> *pDC = GetTplDrive<T>(uLabel);
	if (pDC == 0) {
		silent_cerr("template drive caller " << uLabel
			<< " is undefined at line "
			<< GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pDC->pCopy();
}

const BasicScalarFunction *
MBDynParser::GetScalarFunction(void)
{
	if (!IsKeyWord("reference")) {
		return ParseScalarFunction(*this, pDM);
	}

	std::string s(GetStringWithDelims());
	return GetScalarFunction(s);
}

const BasicScalarFunction *
MBDynParser::GetScalarFunction(const std::string &s)
{
	SFType::const_iterator i = SF.find(s);
	if (i != SF.end()) {
		return i->second;
	}
	return 0;
}

bool
MBDynParser::SetScalarFunction(const std::string &s, const BasicScalarFunction *sf)
{
	return SF.insert(SFType::value_type(s, sf)).second;
}

void
MBDynParser::PopManip(void)
{
	manip.pop();
}

void
MBDynParser::PushManip(const Manip *m)
{
	manip.push(m);
}

const MBDynParser::Manip *
MBDynParser::GetManip(void) const
{
	if (manip.empty()) {
		return 0;
	}

	return manip.top();
}

bool
MBDynParser::bEmptyManip(void) const
{
	return manip.empty();
}

doublereal
MBDynParser::Get(const doublereal& d)
{
	return GetReal(doublereal(d));
}

/* Legge un Vec3 */
Vec3
MBDynParser::GetVec3(void)
{
	Vec3 v(Zero3);
	return GetVec3(v);
}


/* Legge un Vec3 */
Vec3
MBDynParser::GetVec3(const Vec3& vDef)
{
	if (IsKeyWord("default")) {
		return vDef;
	}

	if (IsKeyWord("null")) {
		return Zero3;
	}

	doublereal x1 = GetReal(vDef(1));
	doublereal x2 = GetReal(vDef(2));
	doublereal x3 = GetReal(vDef(3));

	Vec3 v(x1, x2, x3);

	if (IsKeyWord("scale")) {
		v *= GetReal(1.);
	}

	return v;
}


/* Legge una matrice di orientazione sotto forma di due vettori (oppure eye) */
Mat3x3
MBDynParser::GetMatR2vec(void)
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

		Mat3x3 R(r11, r21, r31, r12, r22, r32, r13, r23, r33);
		doublereal dTol = 1.e-12;
		if (IsKeyWord("tolerance")) {
			dTol = GetReal();
			if (dTol < 0.) {
				silent_cerr("invalid (negative) tolerance while reading orientation matrix in \"matr\" format at line " << GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		if (!Eye3.IsSame(R.MulTM(R), dTol)) {
			silent_cerr("warning: orientation matrix orthogonality not within tolerance at line " << GetLineData() << std::endl);
			// not an error?
			// throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return R;
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
			<< ": Euler parameters not implemented yet"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif
	}

	if (IsKeyWord("euler") || IsKeyWord("euler" "123")) {
		doublereal e1 = GetReal();
		doublereal e2 = GetReal();
		doublereal e3 = GetReal();

		return EulerAngles123_2MatR(Vec3(e1, e2, e3));
	}

	if (IsKeyWord("euler" "313")) {
		doublereal e1 = GetReal();
		doublereal e2 = GetReal();
		doublereal e3 = GetReal();

		return EulerAngles313_2MatR(Vec3(e1, e2, e3));
	}

	if (IsKeyWord("euler" "321")) {
		doublereal e1 = GetReal();
		doublereal e2 = GetReal();
		doublereal e3 = GetReal();

		return EulerAngles321_2MatR(Vec3(e1, e2, e3));
	}

	if (IsKeyWord("vector")) {
		doublereal phi1 = GetReal();
		doublereal phi2 = GetReal();
		doublereal phi3 = GetReal();

		return RotManip::Rot(Vec3(phi1, phi2, phi3));
	}

	int i1 = GetInt();
	Vec3 v1;
	v1(1) = GetReal();
	v1(2) = GetReal();
	v1(3) = GetReal();

	int i2 = GetInt();
	Vec3 v2(Zero3);

	if (IsKeyWord("guess")) {
		int i_max = 1;
		int i_min = 1;

		if (fabs(v1(2)) > fabs(v1(1))) {
			i_max = 2;
		} else {
			i_min = 2;
		}

		if (fabs(v1(3)) > fabs(v1(i_max))) {
			i_max = 3;
		} else if (fabs(v1(3)) < fabs(v1(i_min))) {
			i_min = 3;
		}

		v2(i_min) = 1.;
		v2(i_max) = -v1(i_min)/v1(i_max);

	} else {
		v2(1) = GetReal();
		v2(2) = GetReal();
		v2(3) = GetReal();
	}

	return MatR2vec(i1, v1, i2, v2);
}


/* Legge una matrice 3x3 simmetrica come diagonale o triangolare superiore */
Mat3x3
MBDynParser::GetMat3x3Sym(void)
{
	if (IsKeyWord("null")) {
		return Zero3x3;
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
MBDynParser::GetMat3x3(void)
{
	Mat3x3 m(Zero3x3);
	return GetMat3x3(m);
}


/* Legge una matrice 3x3 generica (diagonale o nulla) */
Mat3x3
MBDynParser::GetMat3x3(const Mat3x3& mDef)
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
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m33 = GetReal(mDef(3, 3));
		m = Mat3x3(m11, 0., 0., 0., m22, 0., 0., 0., m33);

	} else if (IsKeyWord("sym")) {
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m12 = GetReal(mDef(1, 2));
		doublereal m13 = GetReal(mDef(1, 3));
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m23 = GetReal(mDef(2, 3));
		doublereal m33 = GetReal(mDef(3, 3));
		m = Mat3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);

	} else if (IsKeyWord("skew")) {
		doublereal v1 = GetReal(mDef(3, 2));
		doublereal v2 = GetReal(mDef(1, 3));
		doublereal v3 = GetReal(mDef(2, 1));
		m = Mat3x3(0., -v3, v2, v3, 0., -v1, -v2, v1, 0.);

	} else {
		if (IsKeyWord("matr")) {
			/* eat it; not required */
			NO_OP;
		}
   
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m12 = GetReal(mDef(1, 2));
		doublereal m13 = GetReal(mDef(1, 3));
		doublereal m21 = GetReal(mDef(2, 1));
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m23 = GetReal(mDef(2, 3));
		doublereal m31 = GetReal(mDef(3, 1));
		doublereal m32 = GetReal(mDef(3, 2));
		doublereal m33 = GetReal(mDef(3, 3));
		m = Mat3x3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
	}

	if (IsKeyWord("scale")) {
		m *= GetReal(1.);
	}

	return m;
}


/* Legge un Vec6 */
Vec6
MBDynParser::GetVec6(void)
{
	Vec6 v(Zero6);
	return GetVec6(v);
}


/* Legge un Vec6 */
Vec6
MBDynParser::GetVec6(const Vec6& vDef)
{
	if (IsKeyWord("null")) {
		return Zero6;
	}

	if (IsKeyWord("default")) {
		return vDef;
	}

	doublereal x1 = GetReal(vDef(1));
	doublereal x2 = GetReal(vDef(2));
	doublereal x3 = GetReal(vDef(3));
	doublereal x4 = GetReal(vDef(4));
	doublereal x5 = GetReal(vDef(5));
	doublereal x6 = GetReal(vDef(6));
	Vec6 v(x1, x2, x3, x4, x5, x6);

	if (IsKeyWord("scale")) {
		v *= GetReal(1.);
	}

	return v;
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6
MBDynParser::GetMat6x6(void)
{
	Mat6x6 m(Zero6x6);
	return GetMat6x6(m);
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6
MBDynParser::GetMat6x6(const Mat6x6& mDef)
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
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m33 = GetReal(mDef(3, 3));
		doublereal m44 = GetReal(mDef(4, 4));
		doublereal m55 = GetReal(mDef(5, 5));
		doublereal m66 = GetReal(mDef(6, 6));
		m = Mat6x6(
			m11, 0., 0., 0., 0., 0.,
			0., m22, 0., 0., 0., 0.,
			0., 0., m33, 0., 0., 0.,
			0., 0., 0., m44, 0., 0.,
			0., 0., 0., 0., m55, 0.,
			0., 0., 0., 0., 0., m66);
 
	} else if (IsKeyWord("sym")) {
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m12 = GetReal(mDef(1, 2));
		doublereal m13 = GetReal(mDef(1, 3));
		doublereal m14 = GetReal(mDef(1, 4));
		doublereal m15 = GetReal(mDef(1, 5));
		doublereal m16 = GetReal(mDef(1, 6));

		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m23 = GetReal(mDef(2, 3));
		doublereal m24 = GetReal(mDef(2, 4));
		doublereal m25 = GetReal(mDef(2, 5));
		doublereal m26 = GetReal(mDef(2, 6));

		doublereal m33 = GetReal(mDef(3, 3));
		doublereal m34 = GetReal(mDef(3, 4));
		doublereal m35 = GetReal(mDef(3, 5));
		doublereal m36 = GetReal(mDef(3, 6));

		doublereal m44 = GetReal(mDef(4, 4));
		doublereal m45 = GetReal(mDef(4, 5));
		doublereal m46 = GetReal(mDef(4, 6));

		doublereal m55 = GetReal(mDef(5, 5));
		doublereal m56 = GetReal(mDef(5, 6));

		doublereal m66 = GetReal(mDef(6, 6));

		m = Mat6x6(
			m11, m12, m13, m14, m15, m16,
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
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m23 = GetReal(mDef(2, 3));
		doublereal m21 = GetReal(mDef(2, 1));
		doublereal m25 = GetReal(mDef(2, 5));
		doublereal m26 = GetReal(mDef(2, 6));
		doublereal m24 = GetReal(mDef(2, 4));

		doublereal m32 = GetReal(mDef(3, 2));
		doublereal m33 = GetReal(mDef(3, 3));
		doublereal m31 = GetReal(mDef(3, 1));
		doublereal m35 = GetReal(mDef(3, 5));
		doublereal m36 = GetReal(mDef(3, 6));
		doublereal m34 = GetReal(mDef(3, 4));

		doublereal m12 = GetReal(mDef(1, 2));
		doublereal m13 = GetReal(mDef(1, 3));
		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m15 = GetReal(mDef(1, 5));
		doublereal m16 = GetReal(mDef(1, 6));
		doublereal m14 = GetReal(mDef(1, 4));

		doublereal m52 = GetReal(mDef(5, 2));
		doublereal m53 = GetReal(mDef(5, 3));
		doublereal m51 = GetReal(mDef(5, 1));
		doublereal m55 = GetReal(mDef(5, 5));
		doublereal m56 = GetReal(mDef(5, 6));
		doublereal m54 = GetReal(mDef(5, 4));

		doublereal m62 = GetReal(mDef(6, 2));
		doublereal m63 = GetReal(mDef(6, 3));
		doublereal m61 = GetReal(mDef(6, 1));
		doublereal m65 = GetReal(mDef(6, 5));
		doublereal m66 = GetReal(mDef(6, 6));
		doublereal m64 = GetReal(mDef(6, 4));

		doublereal m42 = GetReal(mDef(4, 2));
		doublereal m43 = GetReal(mDef(4, 3));
		doublereal m41 = GetReal(mDef(4, 1));
		doublereal m45 = GetReal(mDef(4, 5));
		doublereal m46 = GetReal(mDef(4, 6));
		doublereal m44 = GetReal(mDef(4, 4));

		m = Mat6x6(
			m11, m21, m31, m41, m51, m61,
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

		doublereal m11 = GetReal(mDef(1, 1));
		doublereal m12 = GetReal(mDef(1, 2));
		doublereal m13 = GetReal(mDef(1, 3));
		doublereal m14 = GetReal(mDef(1, 4));
		doublereal m15 = GetReal(mDef(1, 5));
		doublereal m16 = GetReal(mDef(1, 6));

		doublereal m21 = GetReal(mDef(2, 1));
		doublereal m22 = GetReal(mDef(2, 2));
		doublereal m23 = GetReal(mDef(2, 3));
		doublereal m24 = GetReal(mDef(2, 4));
		doublereal m25 = GetReal(mDef(2, 5));
		doublereal m26 = GetReal(mDef(2, 6));

		doublereal m31 = GetReal(mDef(3, 1));
		doublereal m32 = GetReal(mDef(3, 2));
		doublereal m33 = GetReal(mDef(3, 3));
		doublereal m34 = GetReal(mDef(3, 4));
		doublereal m35 = GetReal(mDef(3, 5));
		doublereal m36 = GetReal(mDef(3, 6));

		doublereal m41 = GetReal(mDef(4, 1));
		doublereal m42 = GetReal(mDef(4, 2));
		doublereal m43 = GetReal(mDef(4, 3));
		doublereal m44 = GetReal(mDef(4, 4));
		doublereal m45 = GetReal(mDef(4, 5));
		doublereal m46 = GetReal(mDef(4, 6));

		doublereal m51 = GetReal(mDef(5, 1));
		doublereal m52 = GetReal(mDef(5, 2));
		doublereal m53 = GetReal(mDef(5, 3));
		doublereal m54 = GetReal(mDef(5, 4));
		doublereal m55 = GetReal(mDef(5, 5));
		doublereal m56 = GetReal(mDef(5, 6));

		doublereal m61 = GetReal(mDef(6, 1));
		doublereal m62 = GetReal(mDef(6, 2));
		doublereal m63 = GetReal(mDef(6, 3));
		doublereal m64 = GetReal(mDef(6, 4));
		doublereal m65 = GetReal(mDef(6, 5));
		doublereal m66 = GetReal(mDef(6, 6));

		m = Mat6x6(
			m11, m21, m31, m41, m51, m61,
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
MBDynParser::GetMat6xN(Mat3xN& m1, Mat3xN& m2, integer iNumCols)
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

Vec3
MBDynParser::Get(const Vec3& v)
{
	bool bManip = false;
	Vec3 vv = v;

	while (!bEmptyManip()) {
		const Manip *m = GetManip();
		const TplManip<Vec3> *v3m = dynamic_cast<const TplManip<Vec3> *>(m);
		if (v3m) {
			bManip = true;
			vv = v3m->Get(vv);
		}

		PopManip();
	}

	if (!bManip) {
		vv = GetVec3(v);
	}

	return vv;
}

Mat3x3
MBDynParser::Get(const Mat3x3& m)
{
	return GetMat3x3(m);
}

Vec6
MBDynParser::Get(const Vec6& v)
{
	return GetVec6(v);
}

Mat6x6
MBDynParser::Get(const Mat6x6& m)
{
	return GetMat6x6(m);
}

template <class T>
MBDynParser::TplManip<T>::~TplManip(void)
{
	NO_OP;
}

template <class T>
MBDynParser::RefTplManip<T>::RefTplManip(
	MBDynParser& HP,
	const ReferenceFrame& rf,
	VecMatOpType type)
: HP(HP),
rf(rf),
type(type)
{
	NO_OP;
}

template <class T>
MBDynParser::RefTplManip<T>::~RefTplManip(void)
{
	NO_OP;
}

MBDynParser::RefVec3Manip::RefVec3Manip(
	MBDynParser& HP,
	const ReferenceFrame& rf,
	VecMatOpType type)
: RefTplManip<Vec3>(HP, rf, type)
{
	NO_OP;
}

MBDynParser::RefVec3Manip::~RefVec3Manip(void)
{
	NO_OP;
}

#if 0
Vec6
MBDynParser::RefVec6Manip::Get(const Vec6& v) const
{
	switch (type) {
	case VM_NONE:
	case VM_POSREL:
	case VM_POSABS:
	case VM_VELREL:
	case VM_VELABS:
	case VM_OMEREL:
	case VM_OMEABS:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case VM_VECREL: {
		ReferenceFrame rfOut;
		switch (HP.GetRef(rfOut)) {
		case GLOBAL: {
			Vec6 v(HP.GetVec6());
			return MultRV(v, rf.GetR());
			}

		case UNKNOWNFRAME: /* global */
		case NODE:
			return HighParser::Get(t);

		case LOCAL: {
			Mat3x3 R(GetMatR2vec());
			Vec6 v(HighParser::Get(t)());
			return Vec6(R*v.GetVec1(), R*v.GetVec2());
			}

		case REFERENCE: {
			Vec6 v(HighParser::Get(t)());
			return Vec6(rf.GetR().MulTV(rfOut.GetR()*v.GetVec1()),
				rf.GetR().MulTV(rfOut.GetR()*v.GetVec2()));
			}

		case OTHER_POSITION:
		case OTHER_ORIENTATION:
		case OTHER_NODE:
			silent_cerr("Get<Vec6>: \"other\" meaningless in this context "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
			// will error out later
			break;
		}
		} break;

	case VM_VECABS: {
		ReferenceFrame rfOut;
		switch (GetRef(rfOut)) {
		case UNKNOWNFRAME: /* global */
		case GLOBAL:
			return HighParser::Get(t)();

		case NODE: {
			Vec6 v(HighParser::Get(t)());
			return Vec6(rf.GetR()*v.GetVec1(), rf.GetR()*v.GetVec2());
			}

		case LOCAL: {
			Mat3x3 R(GetMatR2vec());
			Vec6 v(HighParser::Get(t)());
			return Vec6(rf.GetR()*(R*v.GetVec1()), rf.GetR()*(R*v.GetVec2()));
			}

		case REFERENCE: {
			Vec6 v(HighParser::Get(t)());
			return Vec6(rfOut.GetR()*v.GetVec1(), rfOut.GetR()*v.GetVec2());
		}

		case OTHER_POSITION:
		case OTHER_ORIENTATION:
		case OTHER_NODE:
			silent_cerr("Get<Vec6>: \"other\" meaningless in this context "
				"at line " << GetLineData() << std::endl);
			throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
			// will error out later
			break;
		}
		} break;

	case VM_UNITVECREL:
	case VM_UNITVECABS:
		// fall thru

	case VM_MATREL:
	case VM_MATABS:
	case VM_ROTREL:
	case VM_ROTABS:
		// fall thru

	case VM_LAST:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERTMSG(0, "You shouldn't have reached this point");
	throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
}
#endif

Vec3
MBDynParser::RefVec3Manip::Get(const Vec3& v) const
{
	switch (type) {
	case VM_NONE:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case VM_POSREL:
		return HP.GetPosRel(rf);

	case VM_POSABS:
		return HP.GetPosAbs(rf);

	case VM_VELREL:
	case VM_VELABS:
		// need offset
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case VM_OMEREL:
		return HP.GetOmeRel(rf);

	case VM_OMEABS:
		return HP.GetOmeAbs(rf);

	case VM_VECREL:
		return HP.GetVecRel(rf);

	case VM_VECABS:
		return HP.GetVecAbs(rf);

	case VM_UNITVECREL:
		return HP.GetUnitVecRel(rf);

	case VM_UNITVECABS:
		return HP.GetUnitVecAbs(rf);

	case VM_MATREL:
	case VM_MATABS:
	case VM_ROTREL:
	case VM_ROTABS:
		// fall thru

	default:
	case VM_LAST:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return v;
}

MBDynParser::VecRelManip::VecRelManip(
	MBDynParser& HP,
	const ReferenceFrame& rf)
: RefVec3Manip(HP, rf, VM_VECREL)
{
	NO_OP;
}

MBDynParser::VecRelManip::~VecRelManip(void)
{
	NO_OP;
}

MBDynParser::VecAbsManip::VecAbsManip(
	MBDynParser& HP,
	const ReferenceFrame& rf)
: RefVec3Manip(HP, rf, VM_VECABS)
{
	NO_OP;
}

MBDynParser::VecAbsManip::~VecAbsManip(void)
{
	NO_OP;
}

// icc hack
#if defined(__ICC)
static Table dummy_t(false);
static MathParser dummy_mp(dummy_t);
static InputStream dummy_in(std::cin);
static MBDynParser dummy_hp(dummy_mp, dummy_in, "dummy initialization");

void
MBDynParser_dummy_init(void)
{
	TplDriveCaller<doublereal> *dummy_tpldc1d(dummy_hp.GetTplDriveCaller<doublereal>());
	TplDriveCaller<Vec3> *dummy_tpldc3d(dummy_hp.GetTplDriveCaller<Vec3>());
	TplDriveCaller<Vec6> *dummy_tpldc6d(dummy_hp.GetTplDriveCaller<Vec6>());
}
#endif // __ICC
// end of icc hack

/* MBDynParser - end */

