/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

//#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <fstream>

#include "userelem.h"
#include "socketstream_out_elem.h"
#include "socketstreamdrive.h"
#include "bufmod.h"
#include "drive_.h"
#include "streamdrive.h"

#include "module-flightgear.h"

/*FileDrive content-type-------------------------------------------------------------------------------*/

StreamDrive::Modifier * FlightGearFileDriveReader::Read(std::vector<doublereal> &v0, MBDynParser& HP, int &idrives){
	StreamDrive::Modifier *pMod = 0;

	FieldsDescriptionFG *fieldsDescriptionFGInUse;

	if(HP.IsKeyWord("NetFDM")){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetFDM;
		idrives = fieldsDescriptionFGNetFDM.size();
	}
	else if(HP.IsKeyWord("NetCtrls")){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetCtrls;
		idrives = fieldsDescriptionFGNetCtrls.size();
	}
	else{
		silent_cerr("invalid FG data structure "
		"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if(HP.IsKeyWord("print" "options")){
		printOptionsOnTextFile("FGNetFDM_options.txt", fieldsDescriptionFGNetFDM);
		printOptionsOnTextFile("FGNetCtrls_options.txt", fieldsDescriptionFGNetCtrls);
		silent_cout("Flight Gear options correctly printed on files." << std::endl);
	}

	pMod = this->buildFGStreamDriveModifier(fieldsDescriptionFGInUse);

	return pMod;
}

StreamDrive::Modifier *FlightGearFileDriveReader::buildFGStreamDriveModifier(FieldsDescriptionFG *fieldsDescriptionFGInUse){
	StreamDrive::Modifier *pSDM(0);

	std::vector<BufCast *> data(fieldsDescriptionFGInUse->size());
	buildFGBufCast(data, *fieldsDescriptionFGInUse);
	size_t minsize = data[data.size() - 1]->offset() + data[data.size() - 1]->size();
	size_t size = minsize;

	if(fieldsDescriptionFGInUse == &fieldsDescriptionFGNetFDM){
		pSDM = new FGNetFDMStreamDriveCopyCast(size, data);
	}else if(fieldsDescriptionFGInUse == &fieldsDescriptionFGNetCtrls){
		pSDM = new FGNetCtrlsStreamDriveCopyCast(size, data);
	}else{
		silent_cerr("FlightGearFileDriveReader::buildFGStreamDriveModifier internal error: unknown fieldsDescriptionFGInUse " << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pSDM;
}

/*FileDriveCaller type----------------------------------------------------------------------------------*/
integer FlightGearFileDriveCallerTypeReader::Read(const DataManager* pDM, MBDynParser& HP, FileDrive* pDrv){
	FieldsDescriptionFG *fieldsDescriptionFGInUse;
	HighParser::WordSet *FGWordSetInUse;
  std::string FGStructureName;
	const char *s;

	const StreamDrive *pDrv_StreamDrive = dynamic_cast<StreamDrive*>(pDrv);
	const StreamDrive::Modifier *pMod = pDrv_StreamDrive->pGetModifier();

	if(pMod != 0 && dynamic_cast<const FGNetFDMStreamDriveCopyCast*>(pMod) != 0){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetFDM;
		FGWordSetInUse = &fGNetFDMWordSet;
		FGStructureName = "NetFDM";
	}else if(pMod != 0 && dynamic_cast<const FGNetCtrlsStreamDriveCopyCast*>(pMod) != 0){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetCtrls;
		FGWordSetInUse = &fGNetCtrlsWordSet;
		FGStructureName = "NetCtrls";
	}else{
		silent_cerr("error: specified FileDrive does not receive Flight Gear data structure. Please check it out. "
		"error at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if((s = HP.IsWord(*FGWordSetInUse)) != NULL){
		FieldsDescriptionFG::iterator it = fieldsDescriptionFGInUse->find(std::string(s));
		return it->second.position;
	}else{
		silent_cerr("FlightGearFileDriveCallerTypeReader: unknown FG field of " << FGStructureName
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

FGNetFDMStreamDriveCopyCast::FGNetFDMStreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data):StreamDriveCopyCast(size, data){}
FGNetFDMStreamDriveCopyCast::~FGNetFDMStreamDriveCopyCast(void){}

FGNetCtrlsStreamDriveCopyCast::FGNetCtrlsStreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data):StreamDriveCopyCast(size, data){}
FGNetCtrlsStreamDriveCopyCast::~FGNetCtrlsStreamDriveCopyCast(void){}
