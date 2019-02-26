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
#include "socketstreamdrive.h"
#include "socketstream_out_elem.h"

#include "bufmod.h"
#include "drive_.h"

#include "module-flightgear.h"

HighParser::WordSet *sendToFGWordSetInUse;

bool FGNetFDMWordSet::IsWord(const std::string& s) const{
	return fieldsDescriptionFGNetFDM.find(s) != fieldsDescriptionFGNetFDM.end();
}

bool FGNetCtrlsWordSet::IsWord(const std::string& s) const{
	return fieldsDescriptionFGNetCtrls.find(s) != fieldsDescriptionFGNetCtrls.end();
}

StreamContent* FlightGearStreamOutputReader::Read(DataManager* pDM, MBDynParser& HP){
	StreamContent* pSC;

	FieldsDescriptionFG *fieldsDescriptionFGInUse;
	int totChannels;

	if(HP.IsKeyWord("NetFDM")){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetFDM;
		totChannels = fieldsDescriptionFGNetFDM.size();
		sendToFGWordSetInUse = &fGNetFDMWordSet;
	}
	else if(HP.IsKeyWord("NetCtrls")){
		fieldsDescriptionFGInUse = &fieldsDescriptionFGNetCtrls;
		totChannels = fieldsDescriptionFGNetCtrls.size();
		sendToFGWordSetInUse = &fGNetCtrlsWordSet;
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

	/*map containing descriptions and pointers to ScalarValue specified by user for every field he wants to use*/
	FlightGearUserChannels flightGearUserChannels;
	this->ReadFlightGearScalarValues(pDM, HP, flightGearUserChannels, fieldsDescriptionFGInUse);

	/*bulding of buffer that will be sent to Flight gear*/
	std::vector<ScalarValue*> allValues(totChannels);
	std::vector<bool> assignedByUser(totChannels, false);

	/*modifies user-specified fields in the buffer*/
	for(FlightGearUserChannels::iterator i = flightGearUserChannels.begin() ;	i != flightGearUserChannels.end() ; i++){
		FieldsDescriptionFG::iterator el = fieldsDescriptionFGInUse->find(i->first);
		int pos = (el->second.position) -1;
		allValues.at(pos) = i->second;
		assignedByUser.at(pos) = true;
	}

	/*initializes each remaining field with its default value*/
	this->setDefaultValues(allValues, assignedByUser, *fieldsDescriptionFGInUse);

	StreamContent::Modifier *pMod(0);
	pMod = this->buildFGStreamContentModifier(*fieldsDescriptionFGInUse);

	SAFENEWWITHCONSTRUCTOR(pSC, StreamContentValue, StreamContentValue(allValues, pMod));

	return pSC;
}

void FlightGearStreamOutputReader::ReadFlightGearScalarValues(DataManager *pDM, MBDynParser& HP,
	FlightGearUserChannels &flightGearUserChannels, FieldsDescriptionFG *fieldsDescriptionFGInUse){
	const char *fieldDescription;//field the user has written in MBDyn input file
	ScalarValue *scalarValue;
	std::string *s;

	while((fieldDescription = HP.IsWord(*sendToFGWordSetInUse)) != NULL){
		s = new std::string(fieldDescription);
		if(this->fieldAlreadyUsed(std::string(fieldDescription), flightGearUserChannels)){
			silent_cerr("field of FGNetFDM already used "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		scalarValue = ReadScalarValue(pDM, HP);
		flightGearUserChannels.insert(FlightGearUserChannels::value_type(*s, scalarValue));
	}
}

void FlightGearStreamOutputReader::setDefaultValues(std::vector<ScalarValue *> &allValues, std::vector<bool> &assignedByUser, FieldsDescriptionFG &fieldsDescriptionFGInUse){
	DriveCaller *pDC = 0;

	if(!assignedByUser.at(fieldsDescriptionFGInUse.find("version")->second.position-1)){
		if(sendToFGWordSetInUse == &fGNetFDMWordSet){
			SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(FGNetFDMCurVersion));
			allValues.at(fieldsDescriptionFGInUse.find("version")->second.position-1) = new ScalarDriveValue(pDC); //version
		}else if(sendToFGWordSetInUse == &fGNetCtrlsWordSet){
			SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(FGNetCtrlsCurVersion));
			allValues.at(fieldsDescriptionFGInUse.find("version")->second.position-1) = new ScalarDriveValue(pDC); //version
		}
	}

	for (int i=1 ; i < allValues.size() ; i++){
		if(!assignedByUser.at(i)){
			SAFENEW(pDC, NullDriveCaller);
			allValues.at(i) = new ScalarDriveValue(pDC);
		}
	}
}

bool FlightGearStreamOutputReader::fieldAlreadyUsed(std::string fieldDescription, FlightGearUserChannels &flightGearUserChannels){
	if(flightGearUserChannels.find(fieldDescription) != flightGearUserChannels.end()){
		return true;
	}

	return false;
}

StreamContent::Modifier *	FlightGearStreamOutputReader::buildFGStreamContentModifier(FieldsDescriptionFG &fieldsDescriptionFGInUse){
	StreamContent::Modifier *pSCM(0);

	std::vector<BufCast *> data(fieldsDescriptionFGInUse.size());
	buildFGBufCast(data, fieldsDescriptionFGInUse);

	size_t minsize = data[data.size() - 1]->offset() + data[data.size() - 1]->size();
	size_t size = minsize;

	pSCM = new StreamContentCopyCast(0, 0, size, data);

	return pSCM;
}
